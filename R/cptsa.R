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
#'
#' @examples
#' # Survival analysis for TP53 in BRCA transcriptome data
#' result <- cptsa(cancer.type = "BRCA", data.category = "Transcriptome", gene.name = "TP53")
#'
#' # Survival analysis for EGFR in LUAD proteome data
#' result <- cptsa(cancer.type = "LUAD", data.category = "Proteome", gene.name = "EGFR")
#' @importFrom utils read.csv unzip
#' @importFrom survival Surv survfit survdiff
#' @importFrom survminer ggsurvplot
#' @importFrom ggplot2 theme_bw theme element_blank
#'
#' @export
cptsa <- function(cancer.type, data.category, gene.name) {
  # Check the validity of input parameters
  if (missing(cancer.type)) stop("Cancer type is missing")
  if (missing(data.category)) stop("Omics data category is missing")
  if (!data.category %in% c("Transcriptome", "Proteome")) stop("Invalid data category. Only 'Transcriptome' and 'Proteome' are supported")
  if (missing(gene.name)) stop("Gene name is missing")

  # Read cancer-PDC mapping information
  pdc_cancer_info <- read.csv("cancer_PDC_info.csv", check.names = FALSE)

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

  # Helper function to read CSV from zip file
  read_csv_from_zip <- function(zip_path, file_name) {
    # Check if the file exists in the zip
    zip_contents <- unzip(zip_path, list = TRUE)
    if (!(file_name %in% zip_contents$Name)) {
      return(NULL)
    }

    # Create a temporary directory to extract the file
    temp_dir <- tempdir()
    unzip(zip_path, files = file_name, exdir = temp_dir)

    # Read the CSV file
    result <- read.csv(file.path(temp_dir, file_name), check.names = FALSE)

    # Clean up
    unlink(file.path(temp_dir, file_name))

    return(result)
  }

  # Helper function to read CSV with row names from zip file
  read_csv_with_rownames_from_zip <- function(zip_path, file_name) {
    # Check if the file exists in the zip
    zip_contents <- unzip(zip_path, list = TRUE)
    if (!(file_name %in% zip_contents$Name)) {
      return(NULL)
    }

    # Create a temporary directory to extract the file
    temp_dir <- tempdir()
    unzip(zip_path, files = file_name, exdir = temp_dir)

    # Read the CSV file with row names
    result <- read.csv(file.path(temp_dir, file_name), row.names = 1, check.names = FALSE, header = TRUE)

    # Clean up
    unlink(file.path(temp_dir, file_name))

    return(result)
  }

  # Create a list to store results for each PDC ID
  results_list <- list()

  # Perform analysis separately for each PDC ID
  for (pdc_id in pdc_ids) {
    # Initialize results for current PDC ID
    current_result <- list()

    # Clinical data file path (now in zip file with folders)
    clinical_zip_path <- "Clinical data.zip"
    clinical_file_path <- file.path(pdc_id, paste0(pdc_id, ".csv"))

    if (!file.exists(clinical_zip_path)) {
      warning(sprintf("Clinical data zip file not found. Skipping PDC ID %s.", pdc_id))
      next
    }

    # Check if clinical file exists in the zip
    clinical_zip_contents <- unzip(clinical_zip_path, list = TRUE)$Name
    if (!(clinical_file_path %in% clinical_zip_contents)) {
      warning(sprintf("Clinical data file %s not found in %s. Skipping.", clinical_file_path, clinical_zip_path))
      next
    }

    # Read clinical data from zip
    clinical <- read_csv_from_zip(clinical_zip_path, clinical_file_path)
    if (is.null(clinical)) {
      warning(sprintf("Failed to read clinical data for PDC ID %s. Skipping.", pdc_id))
      next
    }

    # Build zip file paths for expression data
    expression_zip_path <- paste0(data.category, ".zip")
    if (!file.exists(expression_zip_path)) {
      warning(sprintf("Zip file %s not found. Skipping.", expression_zip_path))
      next
    }

    # Define file paths within the zip archive
    if (data.category == "Transcriptome") {
      expression_file <- paste0(pdc_id, "_rna_tumor_tidy.csv")
      normal_expression_file <- paste0(pdc_id, "_rna_normal_tidy.csv")
    } else if (data.category == "Proteome") {
      expression_file <- paste0(pdc_id, "_pro_tumor_nor.csv")
      normal_expression_file <- paste0(pdc_id, "_pro_normal_nor.csv")
    }

    # Get the list of zip file contents
    zip_contents <- unzip(expression_zip_path, list = TRUE)$Name

    # Check if tumor expression file exists in the zip
    if (!(expression_file %in% zip_contents)) {
      warning(sprintf("Gene expression file %s not found in %s. Skipping.", expression_file, expression_zip_path))
      next
    }

    # Read tumor gene expression data from zip
    tumor_expression <- read_csv_with_rownames_from_zip(expression_zip_path, expression_file)
    if (is.null(tumor_expression)) {
      warning(sprintf("Failed to read tumor expression data for PDC ID %s. Skipping.", pdc_id))
      next
    }

    # Check if normal tissue data is available in the zip
    has_normal_data <- normal_expression_file %in% zip_contents
    normal_expression <- NULL

    if (has_normal_data) {
      normal_expression <- read_csv_with_rownames_from_zip(expression_zip_path, normal_expression_file)
      if (is.null(normal_expression)) {
        has_normal_data <- FALSE
        warning(sprintf("Failed to read normal expression data for PDC ID %s.", pdc_id))
      }
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


