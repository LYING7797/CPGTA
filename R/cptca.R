#' @title Clinical Proteogenomic Tumor Correlation Analysis(cptca)
#'
#' @description This function calculates the gene-wise correlation between two omics datasets (e.g., Transcriptome and Proteome)
#' for a specified cancer type and gene. It reads, merges, and aligns data from multiple PDC studies if necessary,
#' and visualizes the correlation distribution between tumor and normal samples, highlighting the specified gene.
#'
#' @param gene.name String. The gene name of interest for correlation analysis.
#' @param cancer.type String. Cancer type name or abbreviation.
#' @param data.category Character vector of length 2. Specifies the two omics to compare, e.g., c("Transcriptome", "Proteome").
#'
#' @return Returns a list containing:
#'   \item{correlation_data}{A data frame of correlation values for all genes in tumor and normal samples.}
#'   \item{plot}{A ggplot object visualizing the correlation distributions and highlighting the specified gene.}
#'
#' @examples
#'   result <- cptca(gene.name = "TP53", cancer.type = "BRCA", data.category = c("Transcriptome", "Proteome"))
#' @importFrom utils read.csv unzip
#' @importFrom stats cor t.test
#' @importFrom ggplot2 ggplot aes geom_boxplot scale_fill_manual theme_bw labs theme geom_point geom_text position_jitter annotate
#' @export


cptca <- function(gene.name,
                  cancer.type,
                  data.category) {


  # 1. Parameter checking
  if (missing(gene.name) | missing(cancer.type) | missing(data.category)) {
    stop("Error: Missing required parameters. Please ensure 'gene.name', 'cancer.type', and 'data.category' are provided.")
  }
  if (length(data.category) != 2) {
    stop("Error: 'data.category' must contain exactly 2 omics, e.g., c('Transcriptome','Proteome').")
  }

  valid_categories <- c("Transcriptome", "Proteome", "Phosphoproteome")
  if (!all(data.category %in% valid_categories)) {
    stop("Error: 'data.category' contains invalid types. Valid options: Transcriptome, Proteome, Phosphoproteome.")
  }

  # 2. Read cancer_PDC_info.csv and get PDC codes
  pdc_info_file <- file.path("cancer_PDC_info.csv")
  if (!file.exists(pdc_info_file)) {
    stop("Error: 'cancer_PDC_info.csv' file not found. Please check the path or filename.")
  }

  cancer_info <- read.csv(pdc_info_file, stringsAsFactors = FALSE)

  # Check if cancer.type exists
  matched_rows <- which(cancer_info$cancer_type == cancer.type |
                          cancer_info$abbreviation == cancer.type)
  if (length(matched_rows) == 0) {
    stop("Error: The input 'cancer.type' was not found in 'cancer_PDC_info.csv'.")
  }

  # Return corresponding PDC column name based on omics type
  get_pdc_col_name <- function(omics) {
    if (omics %in% colnames(cancer_info)) {
      return(omics)
    } else {
      stop("Unknown omics type or column not found: ", omics)
    }
  }

  # Get all PDC codes for the omics type in matched_rows
  get_pdc_codes <- function(omics) {
    pdc_col <- get_pdc_col_name(omics)
    codes <- cancer_info[[pdc_col]][matched_rows]
    codes <- unique(codes)
    codes <- codes[!is.na(codes) & codes != ""]
    if (length(codes) == 0) {
      warning(paste("Warning: No corresponding PDC code found for", omics, "in the matched rows."))
      return(character(0))
    }
    return(codes)
  }

  # Get PDC codes for data.category[1] and data.category[2]
  pdc_codes_1 <- get_pdc_codes(data.category[1])
  pdc_codes_2 <- get_pdc_codes(data.category[2])

  if (length(pdc_codes_1) == 0 && length(pdc_codes_2) == 0) {
    stop("Error: No PDC codes found for either omics type, cannot proceed.")
  }

  # 3. Data reading:
  get_omics_file_name <- function(omics, pdc_code, tissue_type = c("tumor","normal")) {
    tissue_type <- match.arg(tissue_type)
    prefix <- switch(omics,
                     "Transcriptome"    = "rna",
                     "Proteome"         = "pro",
                     "Phosphoproteome"  = "phos")

    file_name <- if (tissue_type == "tumor") {
      paste0(pdc_code, "_", prefix, "_tumor_nor.csv")
    } else {
      paste0(pdc_code, "_", prefix, "_normal_nor.csv")
    }
    return(file_name)
  }


  safe_read_from_zip <- function(zip_path, file_name) {
    if (!file.exists(zip_path)) {
      warning("Warning: Zip file does not exist => ", zip_path)
      return(NULL)
    }


    zip_contents <- unzip(zip_path, list = TRUE)$Name


    if (!file_name %in% zip_contents) {
      warning("Warning: File not found in zip => ", file_name)
      return(NULL)
    }


    temp_dir <- tempdir()


    tryCatch({
      unzip(zip_path, files = file_name, exdir = temp_dir, overwrite = TRUE)
      temp_file <- file.path(temp_dir, file_name)
      df <- read.csv(temp_file, stringsAsFactors = FALSE, row.names = 1)
      # 清理临时文件
      unlink(temp_file)
      return(df)
    }, error = function(e) {
      warning("Warning: Error reading file from zip => ", e$message)
      return(NULL)
    })
  }


  read_and_merge_files <- function(omics, tissue_type, pdc_codes) {
    if (length(pdc_codes) == 0) {
      return(NULL)
    }

    zip_path <- file.path(paste0(omics, ".zip"))

    dfs <- list()
    for (code in pdc_codes) {
      fn <- get_omics_file_name(omics, code, tissue_type = tissue_type)
      tmp <- safe_read_from_zip(zip_path, fn)
      if (!is.null(tmp)) {
        dfs[[code]] <- tmp
      }
    }

    if (length(dfs) == 0) {
      return(NULL)
    }

    # Find common genes
    common_genes <- Reduce(intersect, lapply(dfs, rownames))
    if (length(common_genes) == 0) {
      warning("Warning: No common genes found among multiple PDC code data. Returning NULL.")
      return(NULL)
    }

    # Subset common genes
    for (i in seq_along(dfs)) {
      dfs[[i]] <- dfs[[i]][common_genes, , drop = FALSE]
    }

    # Merge columns
    merged_df <- do.call(cbind, dfs)
    return(merged_df)
  }

  # Read omics1
  omics1_tumor   <- read_and_merge_files(data.category[1], "tumor",   pdc_codes_1)
  omics1_normal  <- read_and_merge_files(data.category[1], "normal", pdc_codes_1)

  # Read omics2
  omics2_tumor   <- read_and_merge_files(data.category[2], "tumor",   pdc_codes_2)
  omics2_normal  <- read_and_merge_files(data.category[2], "normal", pdc_codes_2)

  # 4. Correlation calculation function (保持不变)
  calc_correlation_by_gene <- function(df1, df2) {
    if (is.null(df1) || is.null(df2)) {
      return(NULL)
    }

    common_genes <- intersect(rownames(df1), rownames(df2))
    if (length(common_genes) == 0) {
      warning("Warning: No common genes found. Returning empty result.")
      return(NULL)
    }

    common_samples <- intersect(colnames(df1), colnames(df2))
    if (length(common_samples) < 2) {
      warning("Warning: Fewer than 2 common samples. Correlation cannot be computed.")
      return(NULL)
    }

    df1_sub <- df1[common_genes, common_samples, drop = FALSE]
    df2_sub <- df2[common_genes, common_samples, drop = FALSE]

    cor_values <- sapply(common_genes, function(g) {
      x <- as.numeric(df1_sub[g, ])
      y <- as.numeric(df2_sub[g, ])
      cor(x, y, use = "pairwise.complete.obs", method = "pearson")
    })

    res <- data.frame(
      gene = common_genes,
      correlation = cor_values,
      stringsAsFactors = FALSE
    )
    return(res)
  }

  # Calculate gene-wise correlation for tumor and normal samples
  corr_tumor  <- calc_correlation_by_gene(omics1_tumor, omics2_tumor)
  corr_normal <- calc_correlation_by_gene(omics1_normal, omics2_normal)

  # 5. Organize results and visualization (保持不变)
  if (is.null(corr_tumor)) {
    corr_tumor <- data.frame(gene = character(0), correlation = numeric(0), Tissue = character(0))
  } else {
    corr_tumor$Tissue <- "Tumor"
  }

  if (is.null(corr_normal)) {
    corr_normal <- data.frame(gene = character(0), correlation = numeric(0), Tissue = character(0))
  } else {
    corr_normal$Tissue <- "Normal"
  }

  corr_all <- rbind(corr_tumor, corr_normal)

  if (nrow(corr_all) == 0) {
    stop("No available correlation results (files may be missing or no common genes/samples).")
  }

  # Extract the gene of interest
  gene_of_interest <- corr_all[corr_all$gene == gene.name, ]
  if (nrow(gene_of_interest) == 0) {
    warning("Warning: The specified gene was not found in the correlation results: ", gene.name)
  }

  library(ggplot2)

  # Specify colors: Tumor=pink, Normal=light blue
  tissue_colors <- c("Tumor" = "pink", "Normal" = "lightblue")

  p <- ggplot(corr_all, aes(x = Tissue, y = correlation, fill = Tissue)) +
    geom_boxplot(outlier.alpha = 0.3) +
    scale_fill_manual(values = tissue_colors) +
    theme_bw(base_size = 6) +
    labs(
      title = paste0("Correlation between ",
                     paste(data.category, collapse = " & "),
                     " in ", cancer.type,
                     " (Gene: ", gene.name, ")"),
      y = "Correlation (Pearson)"
    ) +
    theme(legend.position = "none")

  # Add points and value labels for the specified gene on the boxplot
  if (nrow(gene_of_interest) > 0) {
    p <- p +
      geom_point(
        data = gene_of_interest,
        aes(x = Tissue, y = correlation),
        color = "black",
        shape = 21,
        stroke = 0.5,
        size = 3,
        position = position_jitter(width = 0.1, height = 0)
      ) +
      geom_text(
        data = gene_of_interest,
        aes(x = Tissue, y = correlation, label = round(correlation, 3)),
        color = "black",
        position = position_jitter(width = 0.1, height = 0),
        vjust = -1,
        size = 3
      )
  }

  # If both Tumor & Normal groups are present, perform significance test and annotate p-value
  tissue_types <- unique(corr_all$Tissue)
  if (length(tissue_types) == 2 && all(c("Tumor", "Normal") %in% tissue_types)) {
    test_res <- t.test(correlation ~ Tissue, data = corr_all)
    p_value <- test_res$p.value

    if (p_value < 1e-5) {
      p_text <- "p < 1e-5"
    } else {
      p_text <- paste0("p = ", format(p_value, digits = 4, scientific = TRUE))
    }

    y_pos <- max(corr_all$correlation, na.rm = TRUE) * 1.05

    p <- p +
      annotate(
        "text",
        x = 1.5,
        y = y_pos,
        label = p_text,
        size = 4,
        color = "black",
        fontface = "bold"
      )
  }

  print(p)

  return(invisible(list(
    correlation_data = corr_all,
    plot = p
  )))
}

