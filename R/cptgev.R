#' @title
#' Clinical Proteogenomic Tumor Gene Expression Violin plot (cptgev)
#'
#' @description
#' This function creates violin plots to compare gene expression levels between tumor and normal tissues
#' across multiple cancer types . It supports visualization
#' of both transcriptomic and proteomic data.
#'
#' @param gene.name Character string specifying the gene symbol to analyze (e.g., "BRCA1", "TP53").
#' @param data.category Character string indicating the data type to analyze, must be either "Proteome" or "Transcriptome".
#'
#' @return A list containing two elements:
#'   \item{plot}{A ggplot object with violin plots showing gene expression in tumor vs. normal tissues across cancer types}
#'   \item{data}{A data frame containing the processed expression data used for plotting}
#'
#' @examples
#' # Visualize BRCA1 expression in transcriptomic data
#' result <- cptgev("BRCA1", "Transcriptome")
#'
#' # Visualize TP53 expression in proteomic data
#' result <- cptgev("TP53", "Proteome")
#'
#' @importFrom ggplot2 ggplot aes geom_violin scale_fill_manual labs theme_bw theme element_text
#' @importFrom ggpubr stat_compare_means
#' @importFrom utils read.csv unzip
#'
#' @export

cptgev <- function(gene.name, data.category) {
  # Parameter check
  if (missing(gene.name) || missing(data.category)) {
    stop("Please provide gene.name and data.category")
  }

  if (!data.category %in% c("Proteome", "Transcriptome")) {
    stop("data.category must be either 'Proteome' or 'Transcriptome'")
  }

  # Set zip file path (在当前目录下)
  zip_path <- paste0(data.category, ".zip")

  # Check if zip file exists
  if (!file.exists(zip_path)) {
    stop(paste0("Zip file not found: ", zip_path))
  }

  # Get all files in the zip archive
  zip_files <- unzip(zip_path, list = TRUE)$Name

  # Filter for tidy.csv files
  tidy_files <- zip_files[grep(".*_tidy\\.csv$", zip_files)]

  # Return error if no files found
  if (length(tidy_files) == 0) {
    stop(paste0("No tidy.csv files found in ", zip_path))
  }

  # Extract PDC codes
  pdc_codes <- unique(substr(gsub("(PDC\\d+)_.*", "\\1", tidy_files), 1, 9))

  # If no directory structure, try different pattern
  if (length(pdc_codes) == 0 || any(is.na(pdc_codes))) {
    pdc_codes <- unique(gsub("(PDC\\d+)_.*", "\\1", tidy_files))
  }

  # Read cancer info file from current directory
  cancer_info_path <- "cancer_PDC_info.csv"
  if (file.exists(cancer_info_path)) {
    cancer_info <- read.csv(cancer_info_path, stringsAsFactors = FALSE)
  } else {
    cancer_info <- NULL
    warning("cancer_PDC_info.csv not found in current directory. Using PDC codes as labels.")
  }

  # Prepare list to store results
  plot_data <- list()

  # Process data for each PDC code
  for (pdc in pdc_codes) {
    # Determine file prefix
    prefix <- if (data.category == "Transcriptome") "rna" else "pro"

    # Build file names for tumor and normal tissue
    tumor_filename <- paste0(pdc, "_", prefix, "_tumor_tidy.csv")
    normal_filename <- paste0(pdc, "_", prefix, "_normal_tidy.csv")

    # Check if files exist in zip
    tumor_exists <- any(grepl(tumor_filename, tidy_files))
    normal_exists <- any(grepl(normal_filename, tidy_files))

    if (!tumor_exists || !normal_exists) {
      next
    }

    # Find full paths in zip
    tumor_file_path <- tidy_files[grepl(tumor_filename, tidy_files)][1]
    normal_file_path <- tidy_files[grepl(normal_filename, tidy_files)][1]

    # Read data from zip file
    tumor_data <- tryCatch({
      read.csv(unz(zip_path, tumor_file_path), row.names = 1)
    }, error = function(e) {
      return(NULL)
    })

    normal_data <- tryCatch({
      read.csv(unz(zip_path, normal_file_path), row.names = 1)
    }, error = function(e) {
      return(NULL)
    })

    # Check if gene exists in both datasets
    if (is.null(tumor_data) || is.null(normal_data) ||
        !gene.name %in% rownames(tumor_data) || !gene.name %in% rownames(normal_data)) {
      next
    }

    # Extract gene expression data
    tumor_expr <- as.numeric(tumor_data[gene.name, ])
    normal_expr <- as.numeric(normal_data[gene.name, ])

    # Get cancer abbreviation if available
    cancer_abbr <- pdc  # Default to PDC code
    if (!is.null(cancer_info)) {
      idx <- which(cancer_info$Proteome == pdc)
      if (length(idx) > 0) {
        cancer_abbr <- cancer_info$abbreviation[idx[1]]
      }
    }

    # Create data frames for plotting
    tumor_df <- data.frame(
      Expression = tumor_expr,
      Group = "Tumor",
      PDC = pdc,
      Cancer = cancer_abbr,
      stringsAsFactors = FALSE
    )

    normal_df <- data.frame(
      Expression = normal_expr,
      Group = "Normal",
      PDC = pdc,
      Cancer = cancer_abbr,
      stringsAsFactors = FALSE
    )

    # Combine data
    plot_data[[pdc]] <- rbind(tumor_df, normal_df)
  }

  # Check if valid data exists
  if (length(plot_data) == 0) {
    stop(paste0("No PDC codes found with both tumor and normal tissue data for gene: ", gene.name))
  }

  # Combine all PDC data
  all_data <- do.call(rbind, plot_data)

  # Load necessary libraries
  library(ggplot2)

  # Create violin plot
  p <- ggplot(all_data, aes(x = Cancer, y = Expression, fill = Group)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.7) +
    scale_fill_manual(values = c("Tumor" = "#FF6666", "Normal" = "#6699CC")) +
    labs(
      title = paste0(gene.name, " Expression in ", data.category),
      x = "Cancer Type",
      y = "Expression Level",
      fill = "Tissue Type"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top"
    )

  # Add statistical significance markers
  library(ggpubr)

  # Perform t-test for each cancer type
  stat_data <- list()
  for (cancer in unique(all_data$Cancer)) {
    cancer_data <- all_data[all_data$Cancer == cancer, ]
    tumor_vals <- cancer_data$Expression[cancer_data$Group == "Tumor"]
    normal_vals <- cancer_data$Expression[cancer_data$Group == "Normal"]

    if (length(tumor_vals) > 1 && length(normal_vals) > 1) {
      t_test <- t.test(tumor_vals, normal_vals)
      p_val <- t_test$p.value

      # Create significance marker
      if (p_val < 0.001) {
        sig <- "***"
      } else if (p_val < 0.01) {
        sig <- "**"
      } else if (p_val < 0.05) {
        sig <- "*"
      } else {
        sig <- "ns"
      }

      # Calculate y position (max value plus some space)
      y_pos <- max(c(tumor_vals, normal_vals), na.rm = TRUE) * 1.1

      stat_data[[length(stat_data) + 1]] <- data.frame(
        Cancer = cancer,
        y = y_pos,
        label = sig,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(stat_data) > 0) {
    stat_df <- do.call(rbind, stat_data)
    p <- p + geom_text(data = stat_df, aes(x = Cancer, y = y, label = label),
                       inherit.aes = FALSE, size = 5)
  }

  # Return plot object and data
  return(list(plot = p, data = all_data))
}










