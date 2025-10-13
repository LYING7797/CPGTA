#' @title
#' Clinical Proteogenomic Tumor Differential Analysis(cptda)
#'
#' @description
#' Perform differential analysis between tumor and normal samples for a specific cancer type using RNA-seq data.
#' The function reads count data from tumor and normal samples based on PDC codes, performs differential expression analysis
#' using DESeq2, and generates volcano plots highlighting significant genes based on specified thresholds.
#'
#' @param cancer.type String. The name or abbreviation of the cancer type.
#' @param log2FC Numeric. Threshold for log2 fold change, up to one decimal place.
#' @param top.gene Integer. Number of top differentially expressed genes to label in the volcano plot (maximum 20).
#' @param gene.name Character vector or NULL. Optional. Names of specific genes to highlight in the volcano plot. If provided and the genes exist in the dataset, they will be highlighted with a different color and size regardless of whether they're in the top genes.
#' @return A list of lists, with each sublist containing:
#'   \item{diff_table}{A data frame of differential analysis results, including gene names, log2FC, p-values, adjusted p-values, and group assignment.}
#'   \item{volcano_plot}{A ggplot2 object of the volcano plot with labeled top genes and optional highlighted gene.}
#'
#' @examples
#' # Differential analysis and volcano plot for BRCA transcriptome data, log2FC threshold 1.5, top 15 genes labeled
#' result <- cptda(cancer.type = "UCEC", log2FC = 1.5, top.gene = 15)
#'
#' # Highlight a specific gene of interest
#' result <- cptda(cancer.type = "LUAD", log2FC = 1.0, top.gene = 10, gene.name = c("KRT5", "KRT14"))


cptda <- function(cancer.type,
                  log2FC,     # numeric, up to 1 decimal
                  top.gene,   # integer, <= 20
                  gene.name = NULL  # gene name(s) for special highlighting, can be a vector
) {

  if (missing(cancer.type) || missing(log2FC) || missing(top.gene)) {
    stop("Error: Missing required parameters. Please provide 'cancer.type', 'log2FC', and 'top.gene'.")
  }

  base_path <- "./data1"
  data.category <- "Transcriptome"  # Default parameter is Transcriptome

  info_file <- file.path(base_path, "cancer_PDC_info.csv")
  if (!file.exists(info_file)) {
    stop("Error: 'cancer_PDC_info.csv' does not exist in the specified path.")
  }

  cancer_info <- read.csv(info_file, stringsAsFactors = FALSE)
  matched_rows <- which(cancer_info$cancer_type == cancer.type |
                          cancer_info$abbreviation == cancer.type)
  if (length(matched_rows) == 0) {
    stop("Error: The input 'cancer.type' was not found in 'cancer_PDC_info.csv'.")
  }

  if (!is.numeric(log2FC)) {
    stop("Error: 'log2FC' must be numeric.")
  }
  if (!grepl("^-?\\d+(\\.\\d)?$", as.character(log2FC))) {
    stop("Error: 'log2FC' must have at most one decimal place.")
  }

  if (!is.numeric(top.gene) || top.gene %% 1 != 0) {
    stop("Error: 'top.gene' must be an integer.")
  }
  if (top.gene > 20) {
    stop("Error: 'top.gene' cannot exceed 20.")
  }

  # 1. Find all PDC codes
  get_pdc_col_name <- function(cat) {
    if (cat %in% colnames(cancer_info)) {
      return(cat)
    } else {
      stop("Unknown category or column not found: ", cat)
    }
  }

  pdc_col <- get_pdc_col_name(data.category)

  pdc_codes <- unique(cancer_info[[pdc_col]][matched_rows])
  pdc_codes <- pdc_codes[!is.na(pdc_codes) & pdc_codes != ""]
  if (length(pdc_codes) == 0) {
    stop("Error: No PDC codes found for this 'cancer.type' in info file.")
  }

  # 2. Read count data for each PDC code and perform differential analysis using DESeq2
  library(ggplot2)
  library(DESeq2)  # Load DESeq2 package
  plots_list <- list()

  for (pdc_code in pdc_codes) {
    # Read count data
    counts_folder <- file.path(base_path, data.category, "Counts_data")

    # Read tumor sample count data
    tumor_counts_file <- file.path(counts_folder, paste0(pdc_code, "_rna_tumor_counts.csv"))
    if (!file.exists(tumor_counts_file)) {
      warning(paste("No tumor counts file found for PDC code:", pdc_code))
      next
    }
    tumor_counts <- read.csv(tumor_counts_file, row.names = 1)
    colnames(tumor_counts) <- paste0("T_",colnames(tumor_counts))

    # Read normal sample count data
    normal_counts_file <- file.path(counts_folder, paste0(pdc_code, "_rna_normal_counts.csv"))
    if (!file.exists(normal_counts_file)) {
      warning(paste("No normal counts file found for PDC code:", pdc_code))
      next
    }
    normal_counts <- read.csv(normal_counts_file, row.names = 1)
    colnames(normal_counts) <- paste0("N_",colnames(normal_counts))


    # Ensure both datasets have common genes
    common_genes <- intersect(rownames(tumor_counts), rownames(normal_counts))
    if (length(common_genes) == 0) {
      warning(paste("No common genes found between tumor and normal counts for PDC code:", pdc_code))
      next
    }

    # Merge count data
    counts_data <- cbind(tumor_counts[common_genes, ], normal_counts[common_genes, ])

    # Create sample information
    condition <- factor(c(rep("tumor", ncol(tumor_counts)), rep("normal", ncol(normal_counts))))
    coldata <- data.frame(row.names = colnames(counts_data), condition = condition)

    # Create DESeq2 dataset object
    dds <- DESeqDataSetFromMatrix(
      countData = round(counts_data),  # DESeq2 requires integer counts
      colData = coldata,
      design = ~ condition
    )

    # Filter low expression genes
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]

    # Run DESeq2 analysis
    dds <- DESeq(dds)

    # Extract results
    res <- results(dds, contrast = c("condition", "tumor", "normal"))
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)

    # Organize results
    res_df <- res_df[, c("gene", "log2FoldChange", "pvalue", "padj")]
    colnames(res_df) <- c("gene", "log2FC", "p.value", "padj")
    res_df$negLog10P <- -log10(res_df$p.value)

    # Filter NA values
    res_df <- res_df[complete.cases(res_df), ]

    # Add group column
    threshold <- log2FC
    res_df$group <- "none"
    res_df$group[res_df$log2FC < -threshold & res_df$padj < 0.05] <- "down"
    res_df$group[res_df$log2FC > threshold & res_df$padj < 0.05] <- "up"

    # Color mapping
    res_df$color_group <- ifelse(res_df$group == "up", "pink",
                                 ifelse(res_df$group == "down", "lightblue", "grey"))

    # Handle multiple genes for highlighting
    res_df$highlight <- FALSE
    if (!is.null(gene.name)) {
      res_df$highlight[res_df$gene %in% gene.name] <- TRUE
    }

    # Select top genes with smallest adjusted p-values for labeling
    filtered_res_df <- res_df[res_df$padj < 0.05, ]
    if (nrow(filtered_res_df) == 0) {
      warning(paste("No genes meet the criteria of padj < 0.05 for PDC code:", pdc_code))
      label_df <- data.frame()
    } else {
      filtered_res_df <- filtered_res_df[order(abs(filtered_res_df$log2FC), decreasing = TRUE), ]
      label_df <- head(filtered_res_df, top.gene)

      # Add any highlighted genes that aren't already in the top genes
      if (!is.null(gene.name)) {
        for (gene in gene.name) {
          if (gene %in% res_df$gene && !(gene %in% label_df$gene)) {
            gene_row <- res_df[res_df$gene == gene, ]
            label_df <- rbind(label_df, gene_row)
          }
        }
      }
    }

    # Create volcano plot - using adjusted p-values
    p <- ggplot(res_df, aes(x = log2FC, y = -log10(padj))) +
      geom_point(aes(color = color_group, size = highlight)) +
      scale_color_identity() +
      scale_size_manual(values = c("FALSE" = 1, "TRUE" = 3)) +
      theme_bw(base_size = 14) +
      labs(
        title = paste0("Volcano Plot for ", cancer.type, " (", pdc_code, ")"),
        x = "log2(FC)",
        y = "-log10(adjusted p-value)"
      ) +
      theme(legend.position = "none") +
      geom_vline(xintercept =  threshold, linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = -threshold, linetype = "dashed", color = "grey40") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40")

    # Add labels for top genes that are not highlighted
    if (nrow(label_df) > 0) {
      p <- p +
        geom_text(data = label_df[!label_df$highlight, ],
                  aes(label = gene, y = -log10(padj)),
                  color = "black",
                  size = 3,
                  vjust = -0.5)
    }

    # Add special highlighting for requested genes
    if (!is.null(gene.name)) {
      highlight_df <- res_df[res_df$highlight, ]
      if (nrow(highlight_df) > 0) {
        p <- p +
          geom_point(data = highlight_df,
                     aes(x = log2FC, y = -log10(padj)),
                     color = "darkred",
                     size = 4) +
          geom_text(data = highlight_df,
                    aes(label = gene, y = -log10(padj)),
                    color = "darkred",
                    size = 4,
                    fontface = "bold",
                    vjust = -0.8)
      }
    }

    # Store results in list
    plots_list[[pdc_code]] <- list(
      diff_table = res_df,
      volcano_plot = p
    )
  }

  # 4. Return results
  return(invisible(plots_list))
}




