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
#'
#' @importFrom utils read.csv unzip
#' @export

cptod <- function(cancer.type = NULL, data.category, data.type = NULL,
                  sample.type = NULL, PDC.study.identifier = NULL) {
  # Base directory path
  base_dir <- "./"

  # Check required parameters
  if (is.null(cancer.type) && is.null(PDC.study.identifier)) {
    stop("Either cancer.type or PDC.study.identifier must be provided.")
  }

  # Define valid data categories
  valid_categories <- c("Biospecimen", "Clinical data", "Copy Number Variation",
                        "Phosphoproteome", "Proteome", "Simple Nucleotide Variation",
                        "Transcriptome", "Counts")

  # Validate data category parameter
  if (!data.category %in% valid_categories) {
    stop(paste("Invalid data category. Valid options include:", paste(valid_categories, collapse = ", ")))
  }

  # Load cancer-PDC mapping information
  cancer_pdc_info_path <- file.path(base_dir, "cancer_PDC_info.csv")
  if (!file.exists(cancer_pdc_info_path)) {
    stop("The file cancer_PDC_info.csv was not found, unable to retrieve the mapping between cancer types and PDC studies.")
  }
  cancer_pdc_info <- read.csv(cancer_pdc_info_path, check.names = FALSE)

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
    if (data.category == "Counts") {
      pdc_column <- "Transcriptome"  # Use Transcriptome column for Counts data
    }

    # Validate whether PDC identifier is valid in the corresponding data.category
    if (!PDC.study.identifier %in% cancer_pdc_info[[pdc_column]]) {
      stop(paste("The provided PDC.study.identifier is invalid in the", pdc_column))
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
    if (data.category == "Counts") {
      pdc_column <- "Transcriptome"  # Use Transcriptome column for Counts data
    }

    # Find all matching PDC IDs according to cancer.type
    matching_rows <- cancer_pdc_info[cancer_pdc_info$cancer_type == cancer.type, ]

    if (nrow(matching_rows) == 0) {
      stop("The specified cancer type was not found. Please check the cancer_PDC_info.csv file for valid cancer types.")
    }

    # Get PDC ID corresponding to data.category
    pdc_ids <- matching_rows[[pdc_column]]
    pdc_ids <- pdc_ids[!is.na(pdc_ids)]

    if (length(pdc_ids) == 0) {
      stop(paste("No PDC studies were found for the specified cancer type ", cancer.type, "under", pdc_column))
    }
  }

  # Helper function:Read CSV from zip file
  read_csv_from_zip <- function(zip_path, file_name) {
    if (!file.exists(zip_path)) {
      warning(paste("Zip file not found:", zip_path))
      return(NULL)
    }

    # Check if the file exists in the zip
    zip_contents <- unzip(zip_path, list = TRUE)$Name
    if (!file_name %in% zip_contents) {
      return(NULL)
    }

    # Create a temporary directory and extract the file
    temp_dir <- tempdir()
    tryCatch({
      unzip(zip_path, files = file_name, exdir = temp_dir, overwrite = TRUE)
      temp_file_path <- file.path(temp_dir, file_name)
      data <- read.csv(temp_file_path, check.names = FALSE)
      unlink(temp_file_path)
      return(data)
    }, error = function(e) {
      warning(paste("Error reading file from zip:", e$message))
      return(NULL)
    })
  }

  # Helper function: Read RDS from zip file
  read_rds_from_zip <- function(zip_path, file_name) {
    if (!file.exists(zip_path)) {
      warning(paste("Zip file not found:", zip_path))
      return(NULL)
    }

    # Check if the file exists in the zip
    zip_contents <- unzip(zip_path, list = TRUE)$Name
    if (!file_name %in% zip_contents) {
      return(NULL)
    }

    # Create a temporary directory and extract the file
    temp_dir <- tempdir()
    tryCatch({
      unzip(zip_path, files = file_name, exdir = temp_dir, overwrite = TRUE)
      temp_file_path <- file.path(temp_dir, file_name)
      data <- readRDS(temp_file_path)
      unlink(temp_file_path)
      return(data)
    }, error = function(e) {
      warning(paste("Error reading file from zip:", e$message))
      return(NULL)
    })
  }

  # Handle Biospecimen data
  if (data.category == "Biospecimen") {
    if (!is.null(data.type) || !is.null(sample.type)) {
      warning("For the Biospecimen category, the data.type and sample.type parameters will be ignored.")
    }

    results <- list()
    zip_path <- file.path(base_dir, "Biospecimen.zip")

    if (!file.exists(zip_path)) {
      stop("Biospecimen.zip file not found.")
    }

    zip_contents <- unzip(zip_path, list = TRUE)$Name

    for (pdc_id in pdc_ids) {
      file_pattern <- paste0("^", pdc_id, "_biospecimen\\.csv$")
      matching_files <- grep(file_pattern, zip_contents, value = TRUE)

      for (file_name in matching_files) {
        data <- read_csv_from_zip(zip_path, file_name)
        if (!is.null(data)) {
          results[[pdc_id]] <- data
        }
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
    zip_path <- file.path(base_dir, "Clinical data.zip")

    if (!file.exists(zip_path)) {
      stop("Clinical data.zip file not found.")
    }

    zip_contents <- unzip(zip_path, list = TRUE)$Name

    for (pdc_id in pdc_ids) {
      file_pattern <- paste0("^", pdc_id, "/")
      matching_files <- grep(file_pattern, zip_contents, value = TRUE)

      if (length(matching_files) > 0) {
        pdc_clinical_data <- list()
        for (file_name in matching_files) {
          if (grepl("\\.csv$", file_name)) {
            data <- read_csv_from_zip(zip_path, file_name)
            if (!is.null(data)) {
              simple_file_name <- basename(file_name)
              pdc_clinical_data[[simple_file_name]] <- data
            }
          }
        }
        if (length(pdc_clinical_data) > 0) {
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
    zip_path <- file.path(base_dir, "Simple Nucleotide Variation.zip")

    if (!file.exists(zip_path)) {
      stop("Simple Nucleotide Variation.zip file not found.")
    }

    for (pdc_id in pdc_ids) {
      file_name <- paste0(pdc_id, ".rds")
      data <- read_rds_from_zip(zip_path, file_name)
      if (!is.null(data)) {
        results[[pdc_id]] <- data
      }
    }

    if (length(results) == 0) {
      stop("No SNV data was found for the specified parameters.")
    }

    return(results)
  }

  # Handle Counts data (special case)
  if (data.category == "Counts") {
    results <- list()

    # Counts now is a folder in current directory
    counts_dir <- file.path(base_dir, "Counts")
    if (!dir.exists(counts_dir)) {
      warning("The Counts directory was not found.")
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
          results[[file_name]] <- read.csv(file, check.names = FALSE)
        }
      }
    }

    if (length(results) == 0) {
      warning("No Counts data was found for the specified parameters.")
      return(NULL)
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
    counts_dir <- file.path(base_dir, "Counts")
    if (!dir.exists(counts_dir)) {
      warning("The Counts directory was not found.")
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
          results[[file_name]] <- read.csv(file, check.names = FALSE)
        }
      }
    }
  }

  # Handle Normalized and tidied data from zip files
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

    # Build the zip file path
    zip_path <- file.path(base_dir, paste0(data.category, ".zip"))
    if (!file.exists(zip_path)) {
      warning(paste(data.category, ".zip file not found."))
      next
    }

    # Get the list of zip file contents
    zip_contents <- unzip(zip_path, list = TRUE)$Name

    for (pdc_id in pdc_ids) {
      if (!is.null(sample.type)) {
        # With sample.type restriction
        for (st in sample.type) {
          file_pattern <- paste0("^", pdc_id, "_", prefix, "_", st, "_.*", suffix, "$")
          matching_files <- grep(file_pattern, zip_contents, value = TRUE)

          for (file_name in matching_files) {
            data <- read_csv_from_zip(zip_path, file_name)
            if (!is.null(data)) {
              results[[file_name]] <- data
            }
          }
        }
      } else {
        # Without sample.type restriction
        file_pattern <- paste0("^", pdc_id, "_", prefix, "_.*", suffix, "$")
        matching_files <- grep(file_pattern, zip_contents, value = TRUE)

        for (file_name in matching_files) {
          data <- read_csv_from_zip(zip_path, file_name)
          if (!is.null(data)) {
            results[[file_name]] <- data
          }
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











