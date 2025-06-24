cptsnv <- function(cancer.type = NULL,
                        PDC.study.identifier = NULL,
                        top_n_genes = 20,
                        min_mut_freq = 0.05) {
  cancer_pdc_info_path = "./data1/cancer_PDC_info.csv"
  # 创建输出目录
  output_dir <- file.path(getwd(), "snv_analysis_results")
  plots_dir <- file.path(output_dir, "plots")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  # 参数验证
  if (is.null(cancer.type) && is.null(PDC.study.identifier)) {
    stop("At least one of the parameters, cancer.type or PDC.study.identifier, must be provided.")
  }

  # 加载必要的包
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

  # 读取cancer_PDC_info.csv文件
  cancer_pdc_info <- read.csv(cancer_pdc_info_path, check.names = FALSE)

  # 确定需要分析的PDC标识符
  if (!is.null(cancer.type)) {
    # 根据cancer.type找到对应的PDC标识符
    matched_pdcs <- cancer_pdc_info %>%
      filter(cancer_type == cancer.type | abbreviation == cancer.type) %>%
      pull(`Simple Nucleotide Variation`) %>%
      strsplit(", ") %>%
      unlist()

    if (length(matched_pdcs) == 0) {
      stop(paste("No matching PDC identifiers found for cancer type:", cancer.type))
    }

    # 获取癌症类型的完整名称和缩写
    cancer_info <- cancer_pdc_info %>%
      filter(cancer_type == cancer.type | abbreviation == cancer.type)
    cancer_full_name <- cancer_info$cancer_type[1]
    cancer_abbr <- cancer_info$abbreviation[1]

    analysis_title <- paste0(cancer_full_name, " (", cancer_abbr, ")")
  } else {
    # 直接使用提供的PDC标识符
    matched_pdcs <- PDC.study.identifier
    analysis_title <- paste0("PDC分析: ", paste(matched_pdcs, collapse=", "))
  }

  # 如果提供了多个PDC，将属于同一癌种的标识符数据合并
  if (length(matched_pdcs) > 1) {
    message(paste("Provided", length(matched_pdcs), "PDC identifiers, merging for analysis:",
                  paste(matched_pdcs, collapse=", ")))

    # 合并多个MAF文件
    combined_maf <- NULL
    successful_pdcs <- c()

    for (pdc in matched_pdcs) {
      maf_file <- file.path("./data1/Simple Nucleotide Variation", paste0(pdc, ".rds"))
      if (file.exists(maf_file)) {
        message(paste("Loading MAF data:", maf_file))
        current_maf <- readRDS(maf_file)

        # 如果数据不是MAF对象，转换为MAF对象
        if (!inherits(current_maf, "MAF")) {
          current_maf <- read.maf(maf = current_maf)
        }

        if (is.null(combined_maf)) {
          combined_maf <- current_maf
        } else {
          # 使用maftools的merge_mafs函数合并MAF对象
          combined_maf <- maftools::merge_mafs(mafs = list(combined_maf, current_maf))
        }

        successful_pdcs <- c(successful_pdcs, pdc)
      } else {
        warning(paste("MAF file not found:", maf_file))
      }
    }

    # 使用合并后的MAF对象
    if (!is.null(combined_maf)) {
      maf_data <- combined_maf
      pdc_id <- paste(successful_pdcs, collapse="_")  # 创建组合ID
      message(paste("Successfully merged", length(successful_pdcs), "PDC identifiers' data:",
                    paste(successful_pdcs, collapse=", ")))
    } else {
      stop("Unable to load MAF data for any of the provided PDC identifiers.")
    }

  } else {
    pdc_id <- matched_pdcs
    # 加载单个MAF文件
    maf_file <- file.path("./data1/Simple Nucleotide Variation", paste0(pdc_id, ".rds"))
    message(paste("Loading MAF data:", maf_file))

    if (!file.exists(maf_file)) {
      stop(paste("File does not exist:", maf_file))
    }

    maf_data <- readRDS(maf_file)

    # 如果数据不是MAF对象，转换为MAF对象
    if (!inherits(maf_data, "MAF")) {
      maf_data <- read.maf(maf = maf_data)
    }
  }

  # 初始化结果列表
  results <- list()

  # 保存MAF对象和分析元数据
  results$maf_data <- maf_data
  results$analysis_info <- list(
    pdc_id = pdc_id,
    cancer_type = cancer.type,
    analysis_title = analysis_title,
    date = Sys.Date(),
    is_merged = length(matched_pdcs) > 1
  )

  # 1. 基本突变统计 - 汇总分析
  message("Performing basic mutation statistics (summary analysis)...")
  results$summary <- getSampleSummary(maf_data)
  results$gene_summary <- getGeneSummary(maf_data)

  # 计算总体统计信息
  total_samples <- length(unique(results$summary$Tumor_Sample_Barcode))
  total_mutations <- sum(results$summary$total)
  total_genes <- length(unique(results$gene_summary$Hugo_Symbol))

  # 计算基因突变频率
  message("Calculating mutation frequency for genes...")
  mut_freq <- results$gene_summary %>%
    mutate(mut_freq = AlteredSamples / total_samples) %>%
    arrange(desc(mut_freq))

  # 如果没有基因突变数据，创建一个空的数据框
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

  # 2. 突变负荷分析 (TMB) - 汇总视图
  message("Analyzing tumor mutation burden (summary analysis)...")
  tmb <- tibble(
    Tumor_Sample_Barcode = results$summary$Tumor_Sample_Barcode,
    TMB = results$summary$total / 38  # 假设使用38MB作为人类外显子组大小
  )

  # 添加TMB汇总统计
  tmb_stats <- list(
    median_tmb = median(tmb$TMB),
    mean_tmb = mean(tmb$TMB),
    min_tmb = min(tmb$TMB),
    max_tmb = max(tmb$TMB),
    q1_tmb = quantile(tmb$TMB, 0.25),
    q3_tmb = quantile(tmb$TMB, 0.75)
  )

  # 绘制TMB分布图 - 汇总视图，添加中位数线
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


  # 创建TMB密度分布图
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

  # 保存TMB柱状图
  ggsave(file.path(plots_dir, "tmb_distribution.pdf"), tmb_plot, width = 12, height = 6)
  # 保存TMB密度图
  ggsave(file.path(plots_dir, "tmb_density.pdf"), tmb_density_plot, width = 10, height = 6)

  # 3 添加突变类型的数据汇总----
  if (!is.null(maf_data@data)) {
    var_type_summary <- as.data.frame(table(maf_data@data$Variant_Type))
    colnames(var_type_summary) <- c("Variant_Type", "Count")
    var_type_summary$Percentage <- var_type_summary$Count / sum(var_type_summary$Count) * 100

    var_class_summary <- as.data.frame(table(maf_data@data$Variant_Classification))
    colnames(var_class_summary) <- c("Variant_Classification", "Count")
    var_class_summary$Percentage <- var_class_summary$Count / sum(var_class_summary$Count) * 100

    results$variant_type_summary <- var_type_summary
    results$variant_classification_summary <- var_class_summary

    # 创建突变类型饼图
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

    # 保存突变类型饼图
    ggsave(file.path(plots_dir, "variant_type_pie.pdf"), vt_pie, width = 8, height = 6)
    results$variant_type_pie <- vt_pie

    # 创建突变分类柱状图
    vc_bar <- ggplot(var_class_summary, aes(x = reorder(Variant_Classification, -Count), y = Count)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste0("Mutation Classification Distribution - ", analysis_title),
           x = "Mutation Classification",
           y = "数量")

    # 保存突变分类柱状图
    ggsave(file.path(plots_dir, "variant_classification_bar.pdf"), vc_bar, width = 10, height = 6)

    results$variant_classification_bar <- vc_bar
  }

  # 4. 突变瀑布图(Oncoplot) - 汇总分析----
  message("Generating oncoplot for high-frequency mutated genes (summary analysis)...")
  tryCatch({
    # 保存为PDF
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

  # 5. 互斥性和共现性分析 - 汇总分析----
  # 正确的方法 - 直接在PDF设备中绘图
  message("Analyzing mutation exclusivity and co-occurrence (summary analysis)...")
  tryCatch({
    # 直接在PDF设备中绘制
    pdf(file.path(plots_dir, "somatic_interactions.pdf"), width = 12, height = 12)

    # 使用参数调整
    si_result <- somaticInteractions(
      maf_data,
      top = top_n_genes,
      pvalue = 0.05,
      leftMar = 6,       # 增加左边距，为左侧的基因名留出更多空间（原默认值为4）
      topMar = 6,        # 增加顶部边距，为顶部的基因名留出更多空间（原默认值为4）
      fontSize = 0.8,     # 略微减小字体大小（原默认值为0.8）
      showSigSymbols = TRUE,
      sigSymbolsSize = 1.5  # 调整显著性符号大小以适应更紧凑的排布
    )

    dev.off()

    # 分别保存结果对象（不绘图）
    si_result <- somaticInteractions(maf_data, top = top_n_genes, pvalue = 0.05, returnAll = TRUE)
    results$somatic_interactions <- si_result

    # 检查结果
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



  # 8. 染色体突变分布（汇总分析）----
  message("Generating chromosome mutation distribution data (summary analysis)...")
  chrom_data <- tryCatch({
    # 提取染色体位置信息
    chromosome_data <- maf_data@data %>%
      mutate(
        # 移除"chr"前缀
        Chromosome = gsub("chr", "", Chromosome)
      ) %>%
      select(Chromosome, Start_Position, End_Position, Variant_Type) %>%
      as.data.frame()

    # 汇总每条染色体上的突变数量
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

    # 设置染色体顺序
    chrom_summary$Chromosome <- factor(chrom_summary$Chromosome,
                                       levels = c(1:22, "X", "Y") %>% as.character())

    chrom_summary
  }, error = function(e) {
    warning(paste("Failed to generate chromosome mutation distribution data:", e$message))
    NULL
  })
  results$chromosome_data <- chrom_data

  # 生成染色体突变分布图
  if(!is.null(chrom_data)) {
    # 按染色体编号排序（处理染色体X、Y）
    chrom_data$Chromosome <- factor(chrom_data$Chromosome,
                                    levels = c(1:22, "X", "Y", "MT") %>% as.character())

    # 染色体突变总量图
    chrom_plot <- ggplot(chrom_data, aes(x = Chromosome, y = Mutations)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      labs(title = paste0("Chromosome Mutation Distribution - ", analysis_title),
           subtitle = paste0("Total mutations: ", total_mutations),
           x = "Chromosome", y = "Number of mutations") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # 保存染色体突变分布图
    ggsave(file.path(plots_dir, "chromosome_mutations.pdf"), chrom_plot, width = 10, height = 6)

    results$chromosome_plot <- chrom_plot

    # 添加按突变类型堆叠的柱状图
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

    # 保存染色体突变类型分布图
    ggsave(file.path(plots_dir, "chromosome_mutation_types.pdf"), chrom_type_plot, width = 10, height = 6)

    results$chromosome_type_plot <- chrom_type_plot
  }

  # 生成汇总图并保存为PDF
  pdf(file.path(plots_dir, "maf_summary.pdf"), width = 12, height = 12)
  plotmafSummary(maf_data,
                 rmOutlier = TRUE,  # 移除异常值
                 addStat = 'median', # 添加中位数统计
                 dashboard = TRUE,   # 使用仪表盘式布局
                 titvRaw = TRUE,     # 显示原始Ti/Tv数据
                 top = 20,           # 显示top 20突变基因
                 textSize = 0.8)     # 文本大小
  dev.off()


  # 创建结果的摘要表格
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

  # 如果存在突变类型数据，添加到摘要表
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

  # 保存分析结果为RDS文件
  saveRDS(results, file = file.path(output_dir, "analysis_results.rds"))

  # 将summary_df保存为CSV文件
  write.csv(summary_df,
            file = file.path(output_dir, "analysis_summary.csv"),
            row.names = FALSE)

  # 创建分析报告文本文件
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

  # 列出所有生成的图形文件
  plot_files <- list.files(plots_dir, pattern = "\\.pdf$")
  for (file in plot_files) {
    report_text <- paste0(report_text, "- ", file, "\n")
  }

  writeLines(report_text, file.path(output_dir, "analysis_report.txt"))

  message(paste0("Analysis completed! Results have been saved to: ", output_dir))

  # 返回结果列表
  return(results)
}




res <- cptsnv(
  cancer.type = "LUSC",
  PDC.study.identifier = NULL,
  top_n_genes = 20,
  min_mut_freq = 0.05
)




