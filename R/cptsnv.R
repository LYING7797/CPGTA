#' @title
#' Clinical Proteogenomic Tumor SNV Analysis (cptsnv)
#'
#' @description
#' This function enables the analysis of somatic single-nucleotide variations (SNVs) in cancer.
#' It supports the identification of key mutations by integrating data from multiple PDC identifiers or cancer types.
#' The function generates comprehensive visualizations, such as mutation frequency plots, oncoplots, and heatmaps,
#' to help uncover potential mutation-based biomarkers for patient stratification and therapeutic insights.
#'
#' @param cancer.type String (optional). The full name or abbreviation of the cancer type (e.g., "BRCA" or "Breast Cancer").
#' @param PDC.study.identifier String or vector (optional). The PDC study identifier(s) to be analyzed.
#' @param top_n_genes Integer. Number of top mutated genes to display in summary plots and co-occurrence analysis. Default is 20.
#' @param min_mut_freq Numeric. Minimum mutation frequency threshold for gene selection. Default is 0.05.
#'
#' @return A list containing analysis results, summary statistics, and generated visualization objects.
#'
#' @examples
#' # SNV analysis for BRCA cancer type
#' result <- cptsnv(cancer.type = "BRCA")
#'
#' # SNV analysis for specific PDC identifiers
#' result <- cptsnv(PDC.study.identifier = c("PDC001", "PDC002"))
#' @importFrom maftools read.maf oncoplot lollipopPlot
#' @importFrom ggplot2 ggplot aes geom_bar theme_bw labs
#' @export


cptsnv <- function(cancer.type = NULL,
                        PDC.study.identifier = NULL,
                        top_n_genes = 20,
                        min_mut_freq = 0.05) {
  cancer_pdc_info_path = "./data1/cancer_PDC_info.csv"
  # Create output directories
  output_dir <- file.path(getwd(), "snv_analysis_results")
  plots_dir <- file.path(output_dir, "plots")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  # Parameter validation
  if (is.null(cancer.type) && is.null(PDC.study.identifier)) {
    stop("At least one of the parameters, cancer.type or PDC.study.identifier, must be provided.")
  }

  # Load required packages
  required_packages <- c("dplyr", "ggplot2", "maftools", "ComplexHeatmap",
                         "BSgenome.Hsapiens.UCSC.hg19", "RColorBrewer",
                         "tibble", "tidyr")

  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing package:", pkg))
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }

  cancer_pdc_info <- read.csv(cancer_pdc_info_path, check.names = FALSE)

  # Determine PDC identifiers for analysis
  if (!is.null(cancer.type)) {
    # Find corresponding PDC identifiers by cancer.type
    matched_pdcs <- cancer_pdc_info %>%
      filter(cancer_type == cancer.type | abbreviation == cancer.type) %>%
      pull(`Simple Nucleotide Variation`) %>%
      strsplit(", ") %>%
      unlist()

    if (length(matched_pdcs) == 0) {
      stop(paste("No matching PDC identifiers found for cancer type:", cancer.type))
    }

    # Get full name and abbreviation of cancer type
    cancer_info <- cancer_pdc_info %>%
      filter(cancer_type == cancer.type | abbreviation == cancer.type)
    cancer_full_name <- cancer_info$cancer_type[1]
    cancer_abbr <- cancer_info$abbreviation[1]

    analysis_title <- paste0(cancer_full_name, " (", cancer_abbr, ")")
  } else {
    # Use provided PDC study identifier directly
    matched_pdcs <- PDC.study.identifier
    analysis_title <- paste0("PDC分析: ", paste(matched_pdcs, collapse=", "))
  }

  # If multiple PDCs are provided, merge data for the same cancer type
  if (length(matched_pdcs) > 1) {
    message(paste("Provided", length(matched_pdcs), "PDC identifiers, merging for analysis:",
                  paste(matched_pdcs, collapse=", ")))

    # Merge multiple MAF files
    combined_maf <- NULL
    successful_pdcs <- c()

    for (pdc in matched_pdcs) {
      maf_file <- file.path("./data1/Simple Nucleotide Variation", paste0(pdc, ".rds"))
      if (file.exists(maf_file)) {
        message(paste("Loading MAF data:", maf_file))
        current_maf <- readRDS(maf_file)

        # Convert to MAF object if not already
        if (!inherits(current_maf, "MAF")) {
          current_maf <- read.maf(maf = current_maf)
        }

        if (is.null(combined_maf)) {
          combined_maf <- current_maf
        } else {
          combined_maf <- maftools::merge_mafs(mafs = list(combined_maf, current_maf))
        }

        successful_pdcs <- c(successful_pdcs, pdc)
      } else {
        warning(paste("MAF file not found:", maf_file))
      }
    }

    # Convert to MAF object if not already
    if (!is.null(combined_maf)) {
      maf_data <- combined_maf
      pdc_id <- paste(successful_pdcs, collapse="_")
      message(paste("Successfully merged", length(successful_pdcs), "PDC identifiers' data:",
                    paste(successful_pdcs, collapse=", ")))
    } else {
      stop("Unable to load MAF data for any of the provided PDC identifiers.")
    }

  } else {
    pdc_id <- matched_pdcs
    # Load single MAF file
    maf_file <- file.path("./data1/Simple Nucleotide Variation", paste0(pdc_id, ".rds"))
    message(paste("Loading MAF data:", maf_file))

    if (!file.exists(maf_file)) {
      stop(paste("File does not exist:", maf_file))
    }

    maf_data <- readRDS(maf_file)

    if (!inherits(maf_data, "MAF")) {
      maf_data <- read.maf(maf = maf_data)
    }
  }

  # Initialize result list
  results <- list()

  # Save MAF object and analysis metadata
  results$maf_data <- maf_data
  results$analysis_info <- list(
    pdc_id = pdc_id,
    cancer_type = cancer.type,
    analysis_title = analysis_title,
    date = Sys.Date(),
    is_merged = length(matched_pdcs) > 1
  )

  # 1. Basic mutation statistics - summary analysis
  message("Performing basic mutation statistics (summary analysis)...")
  results$summary <- getSampleSummary(maf_data)
  results$gene_summary <- getGeneSummary(maf_data)

  # Calculate overall statistics
  total_samples <- length(unique(results$summary$Tumor_Sample_Barcode))
  total_mutations <- sum(results$summary$total)
  total_genes <- length(unique(results$gene_summary$Hugo_Symbol))

  # Calculate gene mutation frequency
  message("Calculating mutation frequency for genes...")
  mut_freq <- results$gene_summary %>%
    mutate(mut_freq = AlteredSamples / total_samples) %>%
    arrange(desc(mut_freq))

  # Create an empty data frame if no mutation data
  if (nrow(mut_freq) == 0) {
    mut_freq <- data.frame(Hugo_Symbol = character(0), mut_freq = numeric(0))
  }

  results$overall_stats <- list(
    total_samples = total_samples,
    total_mutations = total_mutations,
    total_genes = total_genes,
    mutations_per_sample = total_mutations / total_samples,
    median_mutations = median(results$summary$total)
  )

  message(paste("Total number of samples:", total_samples))
  message(paste("Total number of mutations:", total_mutations))
  message(paste("Average mutations per sample:", round(total_mutations / total_samples, 2)))

  # 2. Tumor mutation burden (TMB) analysis - summary view
  message("Analyzing tumor mutation burden (summary analysis)...")
  tmb <- tibble(
    Tumor_Sample_Barcode = results$summary$Tumor_Sample_Barcode,
    TMB = results$summary$total / 38
  )

  # Add TMB summary statistics
  tmb_stats <- list(
    median_tmb = median(tmb$TMB),
    mean_tmb = mean(tmb$TMB),
    min_tmb = min(tmb$TMB),
    max_tmb = max(tmb$TMB),
    q1_tmb = quantile(tmb$TMB, 0.25),
    q3_tmb = quantile(tmb$TMB, 0.75)
  )

  # Plot TMB distribution - summary view with median line
  tmb_plot <- ggplot(tmb, aes(x = reorder(Tumor_Sample_Barcode, -TMB), y = TMB)) +
    geom_bar(stat = "identity", fill = "lightblue") +
    geom_hline(yintercept = tmb_stats$median_tmb, linetype = "dashed", color = "red") +
    annotate("text", x = total_samples/2, y = tmb_stats$median_tmb * 1.1,
             label = paste0("Median:  ", round(tmb_stats$median_tmb, 2)), color = "red") +
    theme_minimal() +
    theme(axis.text.x = element_blank()) +  # 隐藏横坐标文本
    labs(title = paste0("Tumor Mutation Burden Distribution - ", analysis_title),
         subtitle = paste0("Number of samples:  ", total_samples, ", median TMB: ", round(tmb_stats$median_tmb, 2)),
         x = "sample", y = "Mutations per megabase (TMB)")


  # Create TMB density distribution plot
  tmb_density_plot <- ggplot(tmb, aes(x = TMB)) +
    geom_histogram(aes(y = ..density..),
                   bins = 30,
                   fill = "lightblue",
                   color = "black") +
    geom_density(alpha = 0.2, fill = "lightpink",color = "grey40") +
    geom_vline(xintercept = tmb_stats$median_tmb,
               linetype = "dashed",
               color = "red") +
    annotate("text",
             x = tmb_stats$median_tmb * 1.2,
             y = 0.8 * max(density(tmb$TMB)$y),
             label = paste0("Median: ", round(tmb_stats$median_tmb, 2)),
             color = "red") +
    theme_minimal() +
    labs(title = paste0("Tumor Mutation Burden Density Distribution -  ", analysis_title),
         x = "Mutations per megabase (TMB)",
         y = "Density")

  results$tmb_plot <- tmb_plot
  results$tmb_density_plot <- tmb_density_plot
  results$tmb_data <- tmb
  results$tmb_stats <- tmb_stats

  ggsave(file.path(plots_dir, "tmb_distribution.pdf"), tmb_plot, width = 12, height = 6)
  ggsave(file.path(plots_dir, "tmb_density.pdf"), tmb_density_plot, width = 10, height = 6)

  # 3. Add mutation type data summary ----
  if (!is.null(maf_data@data)) {
    var_type_summary <- as.data.frame(table(maf_data@data$Variant_Type))
    colnames(var_type_summary) <- c("Variant_Type", "Count")
    var_type_summary$Percentage <- var_type_summary$Count / sum(var_type_summary$Count) * 100

    var_class_summary <- as.data.frame(table(maf_data@data$Variant_Classification))
    colnames(var_class_summary) <- c("Variant_Classification", "Count")
    var_class_summary$Percentage <- var_class_summary$Count / sum(var_class_summary$Count) * 100

    results$variant_type_summary <- var_type_summary
    results$variant_classification_summary <- var_class_summary

    # Create mutation type pie chart
    vt_pie <- ggplot(var_type_summary, aes(x = "", y = Count, fill = Variant_Type)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) +
      theme_minimal() +
      theme(axis.text = element_blank(), axis.title = element_blank()) +
      scale_fill_brewer(palette = "Set2") +
      labs(title = paste0("Mutation Type Distribution - ", analysis_title),
           fill = "Mutation Type") +
      geom_text(aes(label = paste0(round(Percentage, 1), "%")),
                position = position_stack(vjust = 0.5))

    ggsave(file.path(plots_dir, "variant_type_pie.pdf"), vt_pie, width = 8, height = 6)
    results$variant_type_pie <- vt_pie

    # Create mutation classification bar chart
    vc_bar <- ggplot(var_class_summary, aes(x = reorder(Variant_Classification, -Count), y = Count)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste0("Mutation Classification Distribution - ", analysis_title),
           x = "Mutation Classification",
           y = "Number")

    ggsave(file.path(plots_dir, "variant_classification_bar.pdf"), vc_bar, width = 10, height = 6)

    results$variant_classification_bar <- vc_bar
  }

  # 4. Oncoplot for high-frequency mutated genes - summary analysis ----
  message("Generating oncoplot for high-frequency mutated genes (summary analysis)...")
  tryCatch({
    pdf(file.path(plots_dir, "oncoplot.pdf"), width = 12, height = 8)
    oncoplot(
      maf = maf_data,
      top = 30,
      fontSize = 0.8,
      titleFontSize = 1.2,
      legendFontSize = 1,
      annotationFontSize = 1,
      showTitle = TRUE,
      titleText = paste0("Oncoplot - ", analysis_title),
      showTumorSampleBarcodes = FALSE,
      barcodeSrt = 45,
      drawRowBar = TRUE,
      drawColBar = TRUE,
      removeNonMutated = FALSE,
      colors = NULL,
      bgCol = "#FFFFFF",
      borderCol = "gray90",
      sepwd_genes = 0.5,
      sepwd_samples = 0.25,
      annotationColor = NULL,
      sortByMutation = TRUE,
      showPct = TRUE,
      draw_titv = TRUE,
      gene_mar = 5,
      barcode_mar = 6
    )
    dev.off()

    results$oncoplot_generated <- TRUE
  }, error = function(e) {
    warning(paste("Failed to generate oncoplot:", e$message))
    results$oncoplot_generated <- FALSE
    results$oncoplot_error <- e$message
  })

  # 5. Exclusivity and co-occurrence analysis - summary analysis ----
  message("Analyzing mutation exclusivity and co-occurrence (summary analysis)...")
  tryCatch({
    pdf(file.path(plots_dir, "somatic_interactions.pdf"), width = 12, height = 12)

    # 使用参数调整
    si_result <- somaticInteractions(
      maf_data,
      top = top_n_genes,
      pvalue = 0.05,
      leftMar = 6,
      topMar = 6,
      fontSize = 0.8,
      showSigSymbols = TRUE,
      sigSymbolsSize = 1.5
    )

    dev.off()

    si_result <- somaticInteractions(maf_data, top = top_n_genes, pvalue = 0.05, returnAll = TRUE)
    results$somatic_interactions <- si_result

    if (!is.null(si_result)) {
      results$somatic_interactions_generated <- TRUE
      message("Successfully generated exclusivity/co-occurrence plot.")
    } else {
      results$somatic_interactions_generated <- FALSE
      message("No significant exclusivity/co-occurrence relationships found.")
    }
  }, error = function(e) {
    warning(paste("Failed to analyze mutation exclusivity and co-occurrence:", e$message))
    results$somatic_interactions_generated <- FALSE
    results$somatic_interactions_error <- e$message
  })



  # 8. Chromosome mutation distribution (summary analysis) ----
  message("Generating chromosome mutation distribution data (summary analysis)...")
  chrom_data <- tryCatch({
    chromosome_data <- maf_data@data %>%
      mutate(
        Chromosome = gsub("chr", "", Chromosome)
      ) %>%
      select(Chromosome, Start_Position, End_Position, Variant_Type) %>%
      as.data.frame()

    # Summarize mutation counts for each chromosome
    chrom_summary <- chromosome_data %>%
      group_by(Chromosome) %>%
      summarize(
        Mutations = n(),
        SNP = sum(Variant_Type == "SNP"),
        DEL = sum(Variant_Type == "DEL"),
        INS = sum(Variant_Type == "INS"),
        .groups = 'drop'
      ) %>%
      mutate(Proportion = Mutations / sum(Mutations) * 100)

    # Set chromosome order
    chrom_summary$Chromosome <- factor(chrom_summary$Chromosome,
                                       levels = c(1:22, "X", "Y") %>% as.character())

    chrom_summary
  }, error = function(e) {
    warning(paste("Failed to generate chromosome mutation distribution data:", e$message))
    NULL
  })
  results$chromosome_data <- chrom_data

  # Generate chromosome mutation distribution plot
  if(!is.null(chrom_data)) {
    # Sort by chromosome number (handle chromosomes X, Y)
    chrom_data$Chromosome <- factor(chrom_data$Chromosome,
                                    levels = c(1:22, "X", "Y", "MT") %>% as.character())

    # Chromosome mutation total plot
    chrom_plot <- ggplot(chrom_data, aes(x = Chromosome, y = Mutations)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      labs(title = paste0("Chromosome Mutation Distribution - ", analysis_title),
           subtitle = paste0("Total mutations: ", total_mutations),
           x = "Chromosome", y = "Number of mutations") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(plots_dir, "chromosome_mutations.pdf"), chrom_plot, width = 10, height = 6)

    results$chromosome_plot <- chrom_plot

    # Save chromosome mutation distribution plot
    chrom_type_plot <- ggplot(chrom_data %>%
                                select(Chromosome, SNP, DEL, INS) %>%
                                tidyr::pivot_longer(cols = c(SNP, DEL, INS),
                                                    names_to = "Type",
                                                    values_to = "Count")) +
      geom_bar(aes(x = Chromosome, y = Count, fill = Type), stat = "identity") +
      scale_fill_brewer(palette = "Set1") +
      theme_minimal() +
      labs(title = paste0("Chromosome Mutation Type Distribution - ", analysis_title),
           x = "Chromosome", y = "Number of mutations") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(plots_dir, "chromosome_mutation_types.pdf"), chrom_type_plot, width = 10, height = 6)

    results$chromosome_type_plot <- chrom_type_plot
  }

  # Generate summary plot and save as PDF
  pdf(file.path(plots_dir, "maf_summary.pdf"), width = 12, height = 12)
  plotmafSummary(maf_data,
                 rmOutlier = TRUE,
                 addStat = 'median',
                 dashboard = TRUE,
                 titvRaw = TRUE,
                 top = 20,
                 textSize = 0.8)
  dev.off()


  # Create summary table for results
  summary_df <- data.frame(
    "Analysis Title" = analysis_title,
    "PDC Identifier" = pdc_id,
    "Number of Samples" = total_samples,
    "Total Mutations" = total_mutations,
    "Average Mutations per Sample" = round(total_mutations / total_samples, 2),
    "Median TMB" = round(tmb_stats$median_tmb, 2),
    "Gene with Highest Mutation Frequency" = ifelse(nrow(mut_freq) > 0, mut_freq$Hugo_Symbol[1], "No Data"),
    "Highest Mutation Frequency" = ifelse(nrow(mut_freq) > 0,
                    paste0(round(mut_freq$mut_freq[1] * 100, 2), "%"), "No Data"),
    "Analysis Date" = as.character(Sys.Date())
  )

  # Add mutation type data to summary table if available
  if (exists("var_type_summary") && !is.null(var_type_summary)) {
    snp_row <- which(var_type_summary$Variant_Type == "SNP")
    if (length(snp_row) > 0) {
      summary_df$SNP比例 <- paste0(round(var_type_summary$Percentage[snp_row], 2), "%")
    } else {
      summary_df$SNP比例 <- "No Data"
    }
  } else {
    summary_df$SNP比例 <- "No Data"
  }

  results$summary_table <- summary_df

  saveRDS(results, file = file.path(output_dir, "analysis_results.rds"))
  write.csv(summary_df,
            file = file.path(output_dir, "analysis_summary.csv"),
            row.names = FALSE)

  # Create analysis report text file
  report_text <- paste0(
    "SNV Analysis Report\n",
    "Analysis Time: ", Sys.time(), "\n\n",
    "Basic Information:\n",
    "- Analysis Title: ", analysis_title, "\n",
    "- PDC Identifier: ", pdc_id, "\n",
    "- Number of Samples: ", total_samples, "\n",
    "- Total Mutations: ", total_mutations, "\n",
    "- Average Mutations per Sample:  ", round(total_mutations/total_samples, 2), "\n",
    "- Median TMB:  ", round(tmb_stats$median_tmb, 2), "\n\n",
    "- Output Files:\n",
    "- Analysis Results: analysis_results.rds\n",
    "- Analysis Summary: analysis_summary.csv\n",
    "- - Plot files are located in the plots directory.\n\n",
    "Generated Plot Files:\n"
  )

  # List all generated plot files
  plot_files <- list.files(plots_dir, pattern = "\\.pdf$")
  for (file in plot_files) {
    report_text <- paste0(report_text, "- ", file, "\n")
  }

  writeLines(report_text, file.path(output_dir, "analysis_report.txt"))

  message(paste0("Analysis completed! Results have been saved to: ", output_dir))

  # Return results list
  return(results)
}









