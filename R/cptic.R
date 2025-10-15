#' @title
#' Clinical Proteogenomic Tumor Immunosubtype Classification (cptic)
#'
#' @description
#' This function classifies tumor samples based on their immunological profiles,
#' facilitating the study and understanding of the tumor microenvironment.
#' It integrates clinical and proteogenomic data to assign immunosubtypes to tumor samples
#' for a specified cancer type and omics data category.
#'
#' @param cancer.type String. The full name or abbreviation of the cancer type (e.g., "BRCA" or "Breast Cancer").
#' @param data.category String. Omics data category; must be either "Transcriptome" or "Proteome".
#'
#' @return A data frame containing classified tumor samples and their immunological subtypes
#'         for the specified cancer type and omics category.

#' @examples
#' # Classify transcriptome samples for LUSC
#' result <- cptic(cancer.type = "LUSC", data.category = "Transcriptome")
#'
#' # Classify proteome samples for BRCA
#' result <- cptic(cancer.type = "BRCA", data.category = "Proteome")
#' @importFrom stats hclust dist kmeans
#' @export

cptic <- function(cancer.type, data.category) {
  # Check if all necessary parameters are provided
  if (missing(cancer.type) || missing(data.category)) {
    stop("Error: Missing required input parameters 'cancer.type' or 'data.category'.")
  }

  # Define file paths
  cancer_pdc_info_file <- "cancer_PDC_info.csv"
  alltumor_file <- "./TME/alltumor.csv"
  tme_folder <- "./TME/"

  # Check if files exist
  if (!file.exists(cancer_pdc_info_file)) {
    stop(paste("Error: File not found:", cancer_pdc_info_file))
  }
  if (!file.exists(alltumor_file)) {
    stop(paste("Error: File not found:", alltumor_file))
  }

  # Read cancer_PDC_info.csv
  cancer_pdc_info <- tryCatch({
    read.csv(cancer_pdc_info_file, stringsAsFactors = FALSE)
  }, error = function(e) {
    stop("Error: Unable to read 'cancer_PDC_info.csv' file.")
  })

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

  # Choose the appropriate file based on data.category
  if (data.category == "Transcriptome") {
    classified_file <- file.path(tme_folder, "classified_samples_rna.tsv")
    pdc_column <- "Transcriptome"
  } else { # Proteome
    classified_file <- file.path(tme_folder, "classified_samples_pro.tsv")
    pdc_column <- "Proteome"
  }

  # Check if classified file exists
  if (!file.exists(classified_file)) {
    stop(paste("Error: File not found:", classified_file))
  }

  # Read classified samples file
  classified_samples <- tryCatch({
    read.delim(classified_file, stringsAsFactors = FALSE)
  }, error = function(e) {
    stop(paste("Error: Unable to read file:", classified_file))
  })

  # Read alltumor.csv
  alltumor <- tryCatch({
    read.csv(alltumor_file, row.names = 1, stringsAsFactors = FALSE)
  }, error = function(e) {
    stop("Error: Unable to read 'alltumor.csv' file.")
  })

  # Get the PDC value corresponding to the column name (based on cancer.type)
  # First, find the relevant row for cancer.type in 'cancer_PDC_info.csv'
  cancer_info_row <- cancer_pdc_info[
    cancer_pdc_info$cancer_type == cancer.type | cancer_pdc_info$abbreviation == cancer.type,]

  # Get the corresponding PDC column name and provide specific error information based on data.category
  if (data.category == "Transcriptome") {
    pdc_column_name <- cancer_info_row$Transcriptome
    if (all(!(pdc_column_name %in% alltumor$Transcriptome))) {
      stop(paste("Error: The cancer type", cancer.type, "does not have available transcriptomic immune subtype data."))
    }
  } else { # Proteome
    pdc_column_name <- cancer_info_row$Proteome
    if (all(!(pdc_column_name %in% alltumor$Proteome))) {
      stop(paste("Error: The cancer type", cancer.type, "does not have available proteomic immune subtype data."))
    }
  }



  # Get corresponding PDC values
  pdc_values <- alltumor[[pdc_column]]

  # Get sample names corresponding to PDC values
  sample_names <- alltumor$sample[pdc_values %in% pdc_column_name]

  # Filter classified samples based on sample names
  result_matrix <- classified_samples[classified_samples$sample %in% sample_names, ]

  # Check if any matching samples are found
  if (nrow(result_matrix) == 0) {
    stop("Error: No matching samples found based on the provided parameters.")
  }

  # Return the result matrix
  return(result_matrix)
}











