cptsa <- function(cancer.type, data.category, gene.name) {
  # Check the validity of input parameters
  if (missing(cancer.type)) stop("Cancer type is missing")
  if (missing(data.category)) stop("Omics data category is missing")
  if (!data.category %in% c("Transcriptome", "Proteome")) stop("Invalid data category. Only 'Transcriptome' and 'Proteome' are supported")
  if (missing(gene.name)) stop("Gene name is missing")

  # Read cancer-PDC mapping information
  pdc_cancer_info <- read.csv("./data1/cancer-PDC info.csv")

  valid_cancer_types <- unique(c(cancer_pdc_info$cancer_type, cancer_pdc_info$abbreviation))

  # Validate if cancer.type is within the valid range
  if (!(cancer.type %in% valid_cancer_types)) {
    stop(paste("Error: 'cancer.type' is not in the valid range. Valid values are:", paste(valid_cancer_types, collapse = "\n ")))
  }

  # Validate if data.category is within the valid range
  valid_data_categories <- c("Transcriptome", "Proteome")
  if (!(data.category %in% valid_data_categories)) {
    stop(paste("Error: 'data.category' is not in the valid range. Valid values are:", paste(valid_data_categories, collapse = "\n ")))
  }

  # Get PDC IDs based on cancer type
  pdc_ids <- pdc_cancer_info[pdc_cancer_info$cancer_type == cancer.type | pdc_cancer_info$abbreviation == cancer.type, "PDC_Pro"]

  # Initialize data storage
  combined_clinical <- NULL
  combined_expression <- NULL
  first_iteration <- TRUE

  # Iterate through all PDC IDs, read and combine data
  for (pdc_id in pdc_ids) {
    # Clinical data file path
    clinical_file <- file.path("./data1/Clinical data", pdc_id, paste0(pdc_id, ".csv"))
    if (!file.exists(clinical_file)) {
      warning(sprintf("Clinical data file not found for PDC ID %s. Skipping.", pdc_id))
      next
    }

    # Read clinical data
    clinical <- read.csv(clinical_file, check.names = FALSE)
    combined_clinical <- rbind(combined_clinical, clinical)

    # Gene expression data file path
    if (data.category == "Transcriptome") {
      expression_file <- file.path("./data1", "Transcriptome", paste0(pdc_id, "_rna_tumor_nor.csv"))
    } else if (data.category == "Proteome") {
      expression_file <- file.path("./data1", "Proteome", paste0(pdc_id, "_pro_tumor_nor.csv"))
    }
    if (!file.exists(expression_file)) {
      warning(sprintf("Gene expression file not found for PDC ID %s. Skipping.", pdc_id))
      next
    }

    # Read gene expression data
    expression <- read.csv(expression_file, row.names = 1, check.names = FALSE, header = TRUE)

    # Merge by row names
    if (first_iteration) {
      combined_expression <- expression
      first_iteration <- FALSE
    } else {
      combined_expression <- merge(combined_expression, expression, by = "row.names", all = FALSE)
      rownames(combined_expression) <- combined_expression$Row.names
      combined_expression <- combined_expression[, -1] # Remove merged Row.names column
    }
  }

  clinical <- combined_clinical

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

  # Read gene expression data
  m <- combined_expression

  # Filter target gene and merge clinical data
  gene_expression <- as.data.frame(t(m[rownames(m) == gene.name, ]))
  if (ncol(gene_expression) == 0) stop("The specified gene was not found in the gene expression data")
  median_value <- median(gene_expression[, 1])
  gene_expression$group <- ifelse(gene_expression[, 1] > median_value, "High", "Low")

  clinical_m <- merge(clinical2, gene_expression, by.x = "Cases Submitter ID", by.y = "row.names")

  # Perform survival analysis
  library(survival)
  library(survminer)

  km.by.pro <- survfit(Surv(`Days to Death`, `Vital Status`) ~ group, data = clinical_m)
  logrank_test <- survdiff(Surv(`Days to Death`, `Vital Status`) ~ group, data = clinical_m)
  print(logrank_test)

  # Plot survival curves
  p <- ggsurvplot(km.by.pro,
                  data = clinical_m,
                  pval = TRUE,
                  surv.median.line = "hv",
                  risk.table = TRUE,
                  risk.table.height = .25,
                  xlab = "Time in days",
                  ylab = "Overall Survival Probability",
                  xlim = c(0, NA),
                  legend.title = "")
  return(list(km_fit = km.by.pro, logrank_test = logrank_test, plot = p))
}

cptsa <- function(cancer.type, data.category, gene.name) {
  # Check the validity of input parameters
  if (missing(cancer.type)) stop("Cancer type is missing")
  if (missing(data.category)) stop("Omics data category is missing")
  if (!data.category %in% c("Transcriptome", "Proteome")) stop("Invalid data category. Only 'Transcriptome' and 'Proteome' are supported")
  if (missing(gene.name)) stop("Gene name is missing")

  # Read cancer-PDC mapping information
  pdc_cancer_info <- read.csv("./data1/cancer-PDC info.csv")

  valid_cancer_types <- unique(c(cancer_pdc_info$cancer_type, cancer_pdc_info$abbreviation))

  # Validate if cancer.type is within the valid range
  if (!(cancer.type %in% valid_cancer_types)) {
    stop(paste("Error: 'cancer.type' is not in the valid range. Valid values are:", paste(valid_cancer_types, collapse = "\n ")))
  }

  # Validate if data.category is within the valid range
  valid_data_categories <- c("Transcriptome", "Proteome")
  if (!(data.category %in% valid_data_categories)) {
    stop(paste("Error: 'data.category' is not in the valid range. Valid values are:", paste(valid_data_categories, collapse = "\n ")))
  }

  # Get PDC IDs based on cancer type
  pdc_ids <- pdc_cancer_info[pdc_cancer_info$cancer_type == cancer.type | pdc_cancer_info$abbreviation == cancer.type, "PDC_Pro"]

  # Initialize data storage
  combined_clinical <- NULL
  combined_expression <- NULL
  first_iteration <- TRUE

  # Iterate through all PDC IDs, read and combine data
  for (pdc_id in pdc_ids) {
    # Clinical data file path
    clinical_file <- file.path("./data1/Clinical data", pdc_id, paste0(pdc_id, ".csv"))
    if (!file.exists(clinical_file)) {
      warning(sprintf("Clinical data file not found for PDC ID %s. Skipping.", pdc_id))
      next
    }

    # Read clinical data
    clinical <- read.csv(clinical_file, check.names = FALSE)
    combined_clinical <- rbind(combined_clinical, clinical)

    # Gene expression data file path
    if (data.category == "Transcriptome") {
      expression_file <- file.path("./data1", "Transcriptome", paste0(pdc_id, "_rna_tumor_nor.csv"))
    } else if (data.category == "Proteome") {
      expression_file <- file.path("./data1", "Proteome", paste0(pdc_id, "_pro_tumor_nor.csv"))
    }
    if (!file.exists(expression_file)) {
      warning(sprintf("Gene expression file not found for PDC ID %s. Skipping.", pdc_id))
      next
    }

    # Read gene expression data
    expression <- read.csv(expression_file, row.names = 1, check.names = FALSE, header = TRUE)

    # Merge by row names
    if (first_iteration) {
      combined_expression <- expression
      first_iteration <- FALSE
    } else {
      combined_expression <- merge(combined_expression, expression, by = "row.names", all = FALSE)
      rownames(combined_expression) <- combined_expression$Row.names
      combined_expression <- combined_expression[, -1] # Remove merged Row.names column
    }
  }

  clinical <- combined_clinical

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

  # Read gene expression data
  m <- combined_expression

  # Filter target gene and merge clinical data
  gene_expression <- as.data.frame(t(m[rownames(m) == gene.name, ]))
  if (ncol(gene_expression) == 0) stop("The specified gene was not found in the gene expression data")
  median_value <- median(gene_expression[, 1])
  gene_expression$group <- ifelse(gene_expression[, 1] > median_value, "High", "Low")

  clinical_m <- merge(clinical2, gene_expression, by.x = "Cases Submitter ID", by.y = "row.names")

  # Perform survival analysis
  library(survival)
  library(survminer)

  km.by.pro <- survfit(Surv(`Days to Death`, `Vital Status`) ~ group, data = clinical_m)
  logrank_test <- survdiff(Surv(`Days to Death`, `Vital Status`) ~ group, data = clinical_m)
  print(logrank_test)

  # Plot survival curves
  p <- ggsurvplot(km.by.pro,
                  data = clinical_m,
                  pval = TRUE,
                  surv.median.line = "hv",
                  risk.table = TRUE,
                  risk.table.height = .25,
                  xlab = "Time in days",
                  ylab = "Overall Survival Probability",
                  xlim = c(0, NA),
                  legend.title = "")
  return(list(km_fit = km.by.pro, logrank_test = logrank_test, plot = p))
}





