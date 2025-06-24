cptca <- function(gene.name,
                  cancer.type,
                  data.category) {
  # -----------------------
  # 0. 设置基础路径
  # -----------------------
  base_path <- "./data1"  # 假设所有数据都在 ./data1 目录下

  # -----------------------
  # 1. 参数检查
  # -----------------------
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

  # -----------------------
  # 2. 读取 cancer_PDC_info.csv 并获取 PDC 编码
  # -----------------------
  pdc_info_file <- file.path(base_path, "cancer_PDC_info.csv")
  if (!file.exists(pdc_info_file)) {
    stop("Error: 'cancer_PDC_info.csv' file not found. Please check the path or filename.")
  }

  cancer_info <- read.csv(pdc_info_file, stringsAsFactors = FALSE)

  # 检查 cancer.type 是否存在
  matched_rows <- which(cancer_info$cancer_type == cancer.type |
                          cancer_info$abbreviation == cancer.type)
  if (length(matched_rows) == 0) {
    stop("Error: The input 'cancer.type' was not found in 'cancer_PDC_info.csv'.")
  }

  # 函数：根据组学类型返回对应的 PDC 列名 - 修改后的函数
  get_pdc_col_name <- function(omics) {
    # 现在列名与组学类型相同
    if (omics %in% colnames(cancer_info)) {
      return(omics)
    } else {
      stop("Unknown omics type or column not found: ", omics)
    }
  }

  # 函数：根据 omics 在 matched_rows 中获取所有 PDC code（可能有多个）
  get_pdc_codes <- function(omics) {
    pdc_col <- get_pdc_col_name(omics)
    codes <- cancer_info[[pdc_col]][matched_rows]
    codes <- unique(codes)  # 去重
    codes <- codes[!is.na(codes) & codes != ""]  # 去掉 NA 或空字符串
    if (length(codes) == 0) {
      warning(paste("Warning: No corresponding PDC code found for", omics, "in the matched rows."))
      return(character(0))
    }
    return(codes)
  }

  # 分别取出 data.category[1], data.category[2] 的 PDC codes
  pdc_codes_1 <- get_pdc_codes(data.category[1])
  pdc_codes_2 <- get_pdc_codes(data.category[2])

  if (length(pdc_codes_1) == 0 && length(pdc_codes_2) == 0) {
    stop("Error: No PDC codes found for either omics type, cannot proceed.")
  }

  # -----------------------
  # 3. 数据读取：可能需要合并多个 PDC code 的文件
  # -----------------------
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

  # 安全读取
  safe_read <- function(fp) {
    if (!file.exists(fp)) {
      warning("Warning: File does not exist => ", fp)
      return(NULL)
    }
    df <- read.csv(fp, stringsAsFactors = FALSE, row.names = 1)
    return(df)
  }

  # 同一 Omics + Tissue 类型下，可能有多个 PDC code，需要合并（按行名对齐，列合并）
  read_and_merge_files <- function(omics, tissue_type, pdc_codes) {
    if (length(pdc_codes) == 0) {
      return(NULL)
    }

    dfs <- list()
    for (code in pdc_codes) {
      folder <- file.path(base_path, omics)
      fn <- get_omics_file_name(omics, code, tissue_type = tissue_type)
      fp <- file.path(folder, fn)
      tmp <- safe_read(fp)
      if (!is.null(tmp)) {
        dfs[[code]] <- tmp
      }
    }
    if (length(dfs) == 0) {
      return(NULL)
    }

    # 找共同基因
    common_genes <- Reduce(intersect, lapply(dfs, rownames))
    if (length(common_genes) == 0) {
      warning("Warning: No common genes found among multiple PDC code data. Returning NULL.")
      return(NULL)
    }
    # 截取共同基因
    for (i in seq_along(dfs)) {
      dfs[[i]] <- dfs[[i]][common_genes, , drop = FALSE]
    }
    # 合并列
    merged_df <- do.call(cbind, dfs)
    return(merged_df)
  }

  # 读取 omics1
  omics1_tumor   <- read_and_merge_files(data.category[1], "tumor",   pdc_codes_1)
  omics1_normal  <- read_and_merge_files(data.category[1], "normal", pdc_codes_1)

  # 读取 omics2
  omics2_tumor   <- read_and_merge_files(data.category[2], "tumor",   pdc_codes_2)
  omics2_normal  <- read_and_merge_files(data.category[2], "normal", pdc_codes_2)

  # -----------------------
  # 4. 相关性计算函数
  # -----------------------
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

  # 分别计算肿瘤和正常的基因相关性
  corr_tumor  <- calc_correlation_by_gene(omics1_tumor, omics2_tumor)
  corr_normal <- calc_correlation_by_gene(omics1_normal, omics2_normal)

  # -----------------------
  # 5. 整理结果并可视化
  # -----------------------
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

  # 提取用户关注的基因
  gene_of_interest <- corr_all[corr_all$gene == gene.name, ]
  if (nrow(gene_of_interest) == 0) {
    warning("Warning: The specified gene was not found in the correlation results: ", gene.name)
  }

  # 简单可视化
  library(ggplot2)

  # 指定颜色：Tumor=粉色，Normal=浅蓝
  tissue_colors <- c("Tumor" = "pink", "Normal" = "lightblue")

  p <- ggplot(corr_all, aes(x = Tissue, y = correlation, fill = Tissue)) +
    geom_boxplot(outlier.alpha = 0.3) +
    scale_fill_manual(values = tissue_colors) +
    theme_bw(base_size = 14) +
    labs(
      title = paste0("Correlation between ",
                     paste(data.category, collapse = " & "),
                     " in ", cancer.type,
                     " (Gene: ", gene.name, ")"),
      y = "Correlation (Pearson)"
    ) +
    theme(legend.position = "none")

  # 对指定基因在箱线图上加点及数值标注
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

  # 如果同时包含 Tumor & Normal 两组，则进行显著性差异检验并标注 p 值
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
        x = 1.5,  # 两组中间位置
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



