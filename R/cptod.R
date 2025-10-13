#' @title
#' Clinical Proteogenomic Tumor multi-omics Data Download(cptod)
#'
#' @description This function is used to obtain cancer multi omics data from the data directory.
#' You can filter and download relevant data according to parameters such as cancer type,
#' data category, data type, and sample type.
#'
#' @param cancer.type String. The name or abbreviation of the cancer type.
#' @param data.category String. Data category, options include "Biospecimen", "Clinical data",
#'   "Copy Number Variation", "Phosphoproteome", "Proteome",
#'   "Simple Nucleotide Variation", "Transcriptome".
#' @param data.type String. Data type, usually "Normalized"("nor") or "Tidied"("tidy"). For transcriptome data,
#'   it can be "counts". Default is NULL, meaning no restriction on data type.
#' @param sample.type String or string vector. Sample type, such as "tumor" or "normal".
#'   Default is NULL, meaning no restriction on sample type.
#' @param PDC.study.identifier String. Specific PDC study identifier (e.g., "PDC000125").
#'   Default is NULL, meaning all relevant PDC studies will be found according to cancer type.
#'
#' @return Returns a list containing data files filtered by the parameters.
#' @examples
#' # Download all proteome normalized tumor data for BRCA
#' result <- cptod(cancer.type = "BRCA", data.category = "Proteome", data.type = "Normalized", sample.type = "tumor")
#'
#' # Download biospecimen data for a specific PDC study
#' result <- cptod(PDC.study.identifier = "PDC000125", data.category = "Biospecimen")

#' @import googledrive
#' @importFrom utils download.file unzip read.csv
#' @export

cptod <- function(cancer.type = NULL, data.category, data.type = NULL,
                  sample.type = NULL, PDC.study.identifier = NULL) {
  # Base directory path
  base_dir <- "./data1/"

  # Check required parameters
  if (is.null(cancer.type) && is.null(PDC.study.identifier)) {
    stop("Either cancer.type or PDC.study.identifier must be provided.")
  }

  # Define valid data categories
  valid_categories <- c("Biospecimen", "Clinical data", "Copy Number Variation",
                        "Phosphoproteome", "Proteome", "Simple Nucleotide Variation",
                        "Transcriptome")

  # Validate data category parameter
  if (!data.category %in% valid_categories) {
    stop(paste("Invalid data category. Valid options include:", paste(valid_categories, collapse = ", ")))
  }

  # Load cancer-PDC mapping information
  cancer_pdc_info_path <- file.path(base_dir, "cancer_PDC_info.csv")
  if (!file.exists(cancer_pdc_info_path)) {
    stop("The file cancer_PDC_info.csv was not found, unable to retrieve the mapping between cancer types and PDC studies.")
  }
  cancer_pdc_info <- read.csv("./data1/cancer_PDC_info.csv",check.names = FALSE)

  # Validate cancer.type parameter
  # Get unique cancer types and abbreviations
  valid_cancer_types <- unique(data.frame(
    cancer_type = cancer_pdc_info$cancer_type,
    abbreviation = cancer_pdc_info$abbreviation
  ))
  if (!is.null(cancer.type) && !(cancer.type %in% unlist(valid_cancer_types))) {
    # Format output as two columns, ensure alignment
    formatted_types <- apply(valid_cancer_types, 1, function(row) {
      paste(format(row[1], width = 40, justify = "left"),
            format(row[2], width = 10, justify = "left"))
    })
    header <- paste(format("cancer_type", width = 40, justify = "left"),
                    format("abbreviation", width = 10, justify = "left"))
    valid_types_message <- paste(c(header, formatted_types), collapse = "\n")

    stop(paste("Invalid cancer type. Valid options include:\n\n", valid_types_message))
  }

  # Convert abbreviation to full cancer type name
  if (!is.null(cancer.type) && cancer.type %in% cancer_pdc_info$abbreviation) {
    cancer.type <- unique(cancer_pdc_info[cancer_pdc_info$abbreviation == cancer.type,]$cancer_type)
  }

  # Determine PDC study ID
  if (!is.null(PDC.study.identifier)) {
    pdc_column <- data.category
    # Validate whether PDC identifier is valid in the corresponding data.category
    if (!PDC.study.identifier %in% cancer_pdc_info[[pdc_column]]) {
      stop(paste("The provided PDC.study.identifier is invalid in the", data.category))
    }

    # Validate whether PDC.study.identifier matches cancer.type
    if (!is.null(cancer.type)) {
      matching_rows <- cancer_pdc_info[cancer_pdc_info[[pdc_column]] == PDC.study.identifier, ]
      if (!(cancer.type %in% matching_rows$cancer_type)) {
        stop("The provided PDC.study.identifier does not match the specified cancer.type.")
      }
    }

    pdc_ids <- PDC.study.identifier
  } else {
    pdc_column <- data.category
    # Find all matching PDC IDs according to cancer.type
    matching_rows <- cancer_pdc_info[cancer_pdc_info$cancer_type == cancer.type, ]

    if (nrow(matching_rows) == 0) {
      stop("The specified cancer type was not found. Please check the cancer_PDC_info.csv file for valid cancer types.")
    }

    # Get PDC ID corresponding to data.category
    pdc_column <- data.category
    pdc_ids <- matching_rows[[pdc_column]]
    pdc_ids <- pdc_ids[!is.na(pdc_ids)]

    if (length(pdc_ids) == 0) {
      stop(paste("No PDC studies were found for the specified cancer type ", cancer.type, "under", data.category))
    }
  }

  # Handle Biospecimen data
  if (data.category == "Biospecimen") {
    if (!is.null(data.type) || !is.null(sample.type)) {
      warning("For the Biospecimen category, the data.type and sample.type parameters will be ignored.")
    }

    results <- list()
    for (pdc_id in pdc_ids) {
      file_path <- file.path(base_dir, data.category, paste0(pdc_id, "_biospecimen.csv"))
      if (file.exists(file_path)) {
        results[[pdc_id]] <- read.csv(file_path)
      }
    }

    if (length(results) == 0) {
      stop("No biospecimen data was found for the specified parameters.")
    }

    return(results)
  }

  # Handle Clinical data
  if (data.category == "Clinical data") {
    if (!is.null(data.type) || !is.null(sample.type)) {
      warning("For the Clinical data category, the data.type and sample.type parameters will be ignored.")
    }

    results <- list()
    for (pdc_id in pdc_ids) {
      dir_path <- file.path(base_dir, data.category, pdc_id)
      if (dir.exists(dir_path)) {
        files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
        if (length(files) > 0) {
          pdc_clinical_data <- list()
          for (file in files) {
            file_name <- basename(file)
            pdc_clinical_data[[file_name]] <- read.csv(file)
          }
          results[[pdc_id]] <- pdc_clinical_data
        }
      }
    }

    if (length(results) == 0) {
      stop("No clinical data was found for the specified parameters.")
    }

    return(results)
  }

  # Handle Simple Nucleotide Variation data
  if (data.category == "Simple Nucleotide Variation") {
    if (!is.null(data.type) || !is.null(sample.type)) {
      warning("For the Simple Nucleotide Variation category, the data.type and sample.type parameters will be ignored.")
    }

    results <- list()
    for (pdc_id in pdc_ids) {
      file_path <- file.path(base_dir, data.category, paste0(pdc_id, ".rds"))
      if (file.exists(file_path)) {
        results[[pdc_id]] <- readRDS(file_path)
      }
    }

    if (length(results) == 0) {
      stop("No SNV data was found for the specified parameters.")
    }

    return(results)
  }

  # Handle other data categories (Proteome, Phosphoproteome, Copy Number Variation, Transcriptome)
  # Determine data types to obtain
  valid_data_types <- c("Normalized", "Tidied")
  if (data.category == "Transcriptome") {
    valid_data_types <- c(valid_data_types, "counts")
  }

  # If data.type is not specified, use all valid data types
  if (is.null(data.type)) {
    data.type <- valid_data_types
  } else {
    # Validate data.type
    if (!all(data.type %in% valid_data_types)) {
      stop(paste("Invalid data type. For the", data.category, " valid options:", paste(valid_data_types, collapse = ", ")))
    }
  }

  results <- list()

  # Handle counts data for transcriptome
  if (data.category == "Transcriptome" && "counts" %in% data.type) {
    counts_dir <- file.path(base_dir, "Transcriptome", "Counts_data")
    if (!dir.exists(counts_dir)) {
      warning("The Counts_data directory was not found.")
    } else {
      for (pdc_id in pdc_ids) {
        file_pattern <- paste0("^", pdc_id)
        if (!is.null(sample.type)) {
          sample_patterns <- paste0("_", sample.type, "_", collapse = "|")
          file_pattern <- paste0(file_pattern, ".*?(", sample_patterns, ")")
        }
        file_pattern <- paste0(file_pattern, ".*\\.csv$")

        files <- list.files(counts_dir, pattern = file_pattern, full.names = TRUE)

        for (file in files) {
          file_name <- basename(file)
          results[[file_name]] <- read.csv(file)
        }
      }
    }
  }

  # Handle Normalized and tidied data
  for (dt in data.type) {
    if (dt == "counts") next

    # Determine file name suffix
    suffix <- if (dt == "Normalized") "nor.csv" else "tidy.csv"

    # Determine file name prefix
    prefix <- switch(data.category,
                     "Proteome" = "pro",
                     "Phosphoproteome" = "phos",
                     "Copy Number Variation" = "cnv",
                     "Transcriptome" = "rna",
                     stop(paste("Unsupported data category:", data.category)))

    for (pdc_id in pdc_ids) {
      if (!is.null(sample.type)) {
        # With sample.type restriction
        for (st in sample.type) {
          file_pattern <- paste0("^", pdc_id, "_", prefix, "_", st, "_.*", suffix, "$")
          files <- list.files(file.path(base_dir, data.category), pattern = file_pattern, full.names = TRUE)

          for (file in files) {
            file_name <- basename(file)
            results[[file_name]] <- read.csv(file)
          }
        }
      } else {
        # Without sample.type restriction
        file_pattern <- paste0("^", pdc_id, "_", prefix, "_.*", suffix, "$")
        files <- list.files(file.path(base_dir, data.category), pattern = file_pattern, full.names = TRUE)

        for (file in files) {
          file_name <- basename(file)
          results[[file_name]] <- read.csv(file)
        }
      }
    }
  }

  if (length(results) == 0) {
    warning(paste("No", data.category, "data was found for the specified parameters."))
    return(NULL)
  }

  return(results)
}









