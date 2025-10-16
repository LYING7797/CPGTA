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
  # Check if all required parameters are provided
  if (missing(cancer.type) || missing(data.category)) {
    stop("Error: Missing required input parameters 'cancer.type' or 'data.category'.")
  }
  
  # Directly use system.file to get file paths
  cancer_pdc_info_file <- system.file("extdata", "cancer_PDC_info.csv", package = "CPGTA")
  tme_zip_file <- system.file("extdata", "TME.zip", package = "CPGTA")
  
  # Check if files exist
  if (cancer_pdc_info_file == "") {
    stop("Error: File 'cancer_PDC_info.csv' not found in package data directories")
  }
  
  if (tme_zip_file == "") {
    stop("Error: File 'TME.zip' not found in package data directories")
  }
  
  # Helper function to read files from zip
  read_file_from_zip <- function(zip_path, file_name, is_csv = TRUE, has_header = TRUE, row_names = NULL) {
    if (!file.exists(zip_path)) {
      stop(paste("Error: Zip file not found:", zip_path))
    }
    
    # Check if file exists in zip
    zip_contents <- unzip(zip_path, list = TRUE)$Name
    if (!file_name %in% zip_contents) {
      stop(paste("Error: File", file_name, "not found in", zip_path))
    }
    
    # Create temporary directory and extract file
    temp_dir <- tempdir()
    tryCatch({
      unzip(zip_path, files = file_name, exdir = temp_dir, overwrite = TRUE)
      temp_file_path <- file.path(temp_dir, file_name)
      
      if (is_csv) {
        if (!is.null(row_names)) {
          data <- read.csv(temp_file_path, row.names = row_names, stringsAsFactors = FALSE)
        } else {
          data <- read.csv(temp_file_path, stringsAsFactors = FALSE)
        }
      } else {
        data <- read.delim(temp_file_path, stringsAsFactors = FALSE)
      }
      
      unlink(temp_file_path)
      return(data)
    }, error = function(e) {
      stop(paste("Error reading file from zip:", e$message))
    })
  }
  
  # Read cancer_PDC_info.csv
  cancer_pdc_info <- tryCatch({
    read.csv(cancer_pdc_info_file, stringsAsFactors = FALSE)
  }, error = function(e) {
    stop("Error: Unable to read 'cancer_PDC_info.csv' file. Path tried: ", cancer_pdc_info_file)
  })
  
  valid_cancer_types <- unique(c(cancer_pdc_info$cancer_type, cancer_pdc_info$abbreviation))
  
  # Validate if cancer.type is within valid range
  if (!(cancer.type %in% valid_cancer_types)) {
    stop(paste("Error: 'cancer.type' is not in the valid range. Valid values are:", paste(valid_cancer_types, collapse = "\n ")))
  }
  
  # Validate if data.category is within valid range
  valid_data_categories <- c("Transcriptome", "Proteome")
  if (!(data.category %in% valid_data_categories)) {
    stop(paste("Error: 'data.category' is not in the valid range. Valid values are:", paste(valid_data_categories, collapse = "\n ")))
  }
  
  # Choose appropriate file based on data.category
  if (data.category == "Transcriptome") {
    classified_file <- "classified_samples_rna.tsv"
    pdc_column <- "Transcriptome"
  } else { # Proteome
    classified_file <- "classified_samples_pro.tsv"
    pdc_column <- "Proteome"
  }
  
  # Read alltumor.csv from zip
  alltumor <- read_file_from_zip(tme_zip_file, "alltumor.csv", is_csv = TRUE, row_names = 1)
  
  # Read classified samples file from zip
  classified_samples <- read_file_from_zip(tme_zip_file, classified_file, is_csv = FALSE)
  
  # Get PDC values corresponding to column names (based on cancer.type)
  # First, find the row related to cancer.type in 'cancer_PDC_info.csv'
  cancer_info_row <- cancer_pdc_info[
    cancer_pdc_info$cancer_type == cancer.type | cancer_pdc_info$abbreviation == cancer.type,]
  
  # Get corresponding PDC column name and provide specific error message based on data.category
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
  
  # Check if any matching samples were found
  if (nrow(result_matrix) == 0) {
    stop("Error: No matching samples found based on the provided parameters.")
  }
  
  # Return result matrix
  return(result_matrix)
}











