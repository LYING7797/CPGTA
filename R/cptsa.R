#' @title
#' Clinical Proteogenomic Tumor Survival Analysis (cptsa)
#'
#' @description
#' This function performs survival analysis by correlating clinical outcomes with molecular profiles
#' for a specified cancer type, omics data category, and gene of interest.
#' It generates survival curves to visualize the relationship between gene expression and patient prognosis,
#' helping to identify potential biomarkers for patient stratification.
#'
#' @param cancer.type String. The full name or abbreviation of the cancer type (e.g., "BRCA" or "Breast Cancer").
#' @param data.category String. Omics data category; must be either "Transcriptome" or "Proteome".
#' @param gene.name String. The gene symbol or name to be analyzed in the survival analysis.
#'
#' @return A list containing:
#'   \item{km_fit}{A 'survfit' object representing the Kaplan-Meier survival fit.}
#'   \item{logrank_test}{A log-rank test result comparing survival between groups.}
#'   \item{plot}{A 'ggsurvplot' object visualizing the survival curves and statistics.}
#'
#' @examples
#' # Survival analysis for TP53 in BRCA transcriptome data
#' result <- cptsa(cancer.type = "BRCA", data.category = "Transcriptome", gene.name = "TP53")
#'
#' # Survival analysis for EGFR in LUAD proteome data
#' result <- cptsa(cancer.type = "LUAD", data.category = "Proteome", gene.name = "EGFR")


cptsa <- function(cancer.type, data.category, gene.name) {
  # Check the validity of input parameters
  if (missing(cancer.type)) stop("Cancer type is missing")
  if (missing(data.category)) stop("Omics data category is missing")
  if (!data.category %in% c("Transcriptome", "Proteome")) stop("Invalid data category. Only 'Transcriptome' and 'Proteome' are supported")
  if (missing(gene.name)) stop("Gene name is missing")

  # Read cancer-PDC mapping information
  pdc_cancer_info <- read.csv("./data1/cancer_PDC_info.csv",check.names = FALSE)

  valid_cancer_types <- unique(c(pdc_cancer_info$cancer_type, pdc_cancer_info$abbreviation))

  # Validate if cancer.type is within the valid range
  if (!(cancer.type %in% valid_cancer_types)) {
    stop(paste("Error: 'cancer.type' is not in the valid range. Valid values are:", paste(valid_cancer_types, collapse = "\n ")))
  }

  # Validate if data.category is within the valid range
  valid_data_categories <- c("Transcriptome", "Proteome")
  if (!(data.category %in% valid_data_categories)) {
    stop(paste("Error: 'data.category' is not in the valid range. Valid values are:", paste(valid_data_categories, collapse = "\n ")))
  }

  # Dynamically select column names based on data.category
  if (data.category == "Proteome") {
    pdc_ids <- pdc_cancer_info[pdc_cancer_info$cancer_type == cancer.type | pdc_cancer_info$abbreviation == cancer.type, "Proteome"]
  } else if (data.category == "Transcriptome") {
    pdc_ids <- pdc_cancer_info[pdc_cancer_info$cancer_type == cancer.type | pdc_cancer_info$abbreviation == cancer.type, "Transcriptome"]
  } else {
    stop("Invalid data.category. Please specify either 'Proteome' or 'Transcriptome'.")
  }

  # Filter out empty values
  pdc_ids <- pdc_ids[!is.na(pdc_ids) & pdc_ids != ""]

  if(length(pdc_ids) == 0) {
    stop(paste("No PDC IDs found for", cancer.type, "in", data.category))
  }

  # Create a list to store results for each PDC ID
  results_list <- list()

  # Perform analysis separately for each PDC ID
  for (pdc_id in pdc_ids) {
    # Initialize results for current PDC ID
    current_result <- list()

    # Clinical data file path
    clinical_file <- file.path("./data1/Clinical data", pdc_id, paste0(pdc_id, ".csv"))
    if (!file.exists(clinical_file)) {
      warning(sprintf("Clinical data file not found for PDC ID %s. Skipping.", pdc_id))
      next
    }

    # Read clinical data
    clinical <- read.csv(clinical_file, check.names = FALSE)

    # Gene expression data file path
    if (data.category == "Transcriptome") {
      expression_file <- file.path("./data1", "Transcriptome", paste0(pdc_id, "_rna_tumor_tidy.csv"))
      normal_expression_file <- file.path("./data1", "Transcriptome", paste0(pdc_id, "_rna_normal_tidy.csv"))
    } else if (data.category == "Proteome") {
      expression_file <- file.path("./data1", "Proteome", paste0(pdc_id, "_pro_tumor_nor.csv"))
      normal_expression_file <- file.path("./data1", "Proteome", paste0(pdc_id, "_pro_normal_nor.csv"))
    }

    if (!file.exists(expression_file)) {
      warning(sprintf("Gene expression file not found for PDC ID %s. Skipping.", pdc_id))
      next
    }

    # Read tumor gene expression data
    tumor_expression <- read.csv(expression_file, row.names = 1, check.names = FALSE, header = TRUE)

    # Check if normal tissue data is available
    has_normal_data <- file.exists(normal_expression_file)
    normal_expression <- NULL

    if (has_normal_data) {
      normal_expression <- read.csv(normal_expression_file, row.names = 1, check.names = FALSE, header = TRUE)
    }

    # Keep samples with "Vital Status" as Alive or Dead
    clinical1 <- clinical[clinical$`Vital Status` %in% c("Alive", "Dead"), ]
    clinical1$`Vital Status` <- ifelse(clinical1$`Vital Status` == "Dead", 1, 0)

    # Fill missing Days to Death values
    na_and_alive <- is.na(clinical1$`Days to Death`) & (clinical1$`Vital Status` == 0)
    clinical1$`Days to Death`[na_and_alive] <- ifelse(!is.na(clinical1$`Days to Last Known Disease Status`[na_and_alive]),
                                                      clinical1$`Days to Last Known Disease Status`[na_and_alive],
                                                      clinical1$`Days to Last Follow Up`[na_and_alive])

    # Remove cases where death occurred before study start or Days to Death is NA
    clinical2 <- clinical1[clinical1$`Days to Death` >= 0 & !is.na(clinical1$`Days to Death`), ]
    clinical2 <- clinical2[, c("Cases Submitter ID", "Cause of Death", "Days to Death", "Vital Status", "Days to Last Follow Up", "Days to Last Known Disease Status")]

    # Filter target gene and prepare expression data
    if (!(gene.name %in% rownames(tumor_expression))) {
      warning(sprintf("Gene %s not found in expression data for PDC ID %s. Skipping.", gene.name, pdc_id))
      next
    }

    gene_expression <- as.data.frame(t(tumor_expression[rownames(tumor_expression) == gene.name, , drop=FALSE]))
    colnames(gene_expression) <- gene.name

    # Determine up/down based on normal data or quantiles
    if (has_normal_data && gene.name %in% rownames(normal_expression)) {
      # Calculate mean expression for normal samples
      normal_mean <- mean(as.numeric(normal_expression[gene.name, ]))
      gene_expression$group <- ifelse(gene_expression[, 1] > normal_mean, "Up", "Down")
    } else {
      # Use quantiles to define up/down
      quantile_30 <- quantile(gene_expression[, 1], probs = c(0.3, 0.7))
      gene_expression$group <- ifelse(gene_expression[, 1] > quantile_30[2], "Up",
                                      ifelse(gene_expression[, 1] < quantile_30[1], "Down", NA))
      gene_expression <- gene_expression[!is.na(gene_expression$group), ] # Remove middle 40%
    }

    clinical_m <- merge(clinical2, gene_expression, by.x = "Cases Submitter ID", by.y = "row.names")

    # Check if there is sufficient data for survival analysis
    if (nrow(clinical_m) < 5) {
      warning(sprintf("Not enough data for survival analysis for PDC ID %s. Skipping.", pdc_id))
      next
    }

    # Perform survival analysis
    library(survival)
    library(survminer)

    km.by.pro <- survfit(Surv(`Days to Death`, `Vital Status`) ~ group, data = clinical_m)
    logrank_test <- survdiff(Surv(`Days to Death`, `Vital Status`) ~ group, data = clinical_m)


    group_levels <- levels(factor(clinical_m$group))
    color_palette <- ifelse(group_levels == "Down", "skyblue", "lightpink")

    p <- ggsurvplot(
      km.by.pro,
      data = clinical_m,
      pval = TRUE,
      surv.median.line = "hv",
      risk.table = TRUE,
      risk.table.height = .25,
      title = paste("Survival Analysis for", pdc_id,"_",cancer.type,"_",gene.name),
      xlab = "Time in days",
      ylab = "Overall Survival Probability",
      xlim = c(0, NA),
      legend.title = "",
      palette = color_palette,
      ggtheme = theme_bw()+
        theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    )

    # Save results for the current PDC ID
    results_list[[paste(cancer.type,"_",pdc_id)]] <- list(
      km_fit = km.by.pro,
      logrank_test = logrank_test,
      plot = p,
      clinical_data = clinical_m
    )
  }

  # If there are no results, give a warning
  if (length(results_list) == 0) {
    warning("No results were generated for any PDC ID.")
  }

  return(results_list)
}




