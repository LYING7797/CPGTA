cptda <- function(cancer.type,
                  data.category,
                  log2FC,     # numeric, up to 1 decimal
                  top.gene    # integer, <= 20
) {
  # -----------------------
  # 0. 参数检查
  # -----------------------

  if (missing(cancer.type) || missing(data.category) ||
      missing(log2FC) || missing(top.gene)) {
    stop("Error: Missing required parameters. Please provide 'cancer.type', 'data.category', 'log2FC', and 'top.gene'.")
  }

  base_path <- "./data1"

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

  valid_categories <- c("Transcriptome","Proteome")
  if (! data.category %in% valid_categories) {
    stop("Error: 'data.category' must be either 'Transcriptome' or 'Proteome'.")
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

  # -----------------------
  # 1. 找到所有 PDC code，并全部合并
  # -----------------------
  # 修改后的函数 - 直接返回组学类型作为列名
  get_pdc_col_name <- function(cat) {
    # 检查是否是有效的组学类型且列在数据框中存在
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
    stop("Error: No PDC codes found for this 'cancer.type' & 'data.category' in info file.")
  }

  # -----------------------
  # 2. 读取并合并多个 PDC code 的 tumor 和 normal
  # -----------------------
  get_file_name <- function(cat, pdc, tissue = c("tumor","normal")) {
    tissue <- match.arg(tissue)
    prefix <- if (cat == "Transcriptome") "rna" else "pro"
    if (tissue == "tumor") {
      paste0(pdc, "_", prefix, "_tumor_nor.csv")
    } else {
      paste0(pdc, "_", prefix, "_normal_nor.csv")
    }
  }

  safe_read <- function(fp) {
    if (!file.exists(fp)) {
      warning("Warning: File does not exist => ", fp)
      return(NULL)
    }
    df <- read.csv(fp, stringsAsFactors = FALSE, row.names = 1)
    return(df)
  }

  read_and_merge_files <- function(cat, pdc_list, tissue_type) {
    dfs <- list()
    for (code in pdc_list) {
      folder <- file.path(base_path, cat)
      fn <- get_file_name(cat, code, tissue_type)
      fp <- file.path(folder, fn)
      tmp <- safe_read(fp)
      if (!is.null(tmp)) {
        dfs[[code]] <- tmp
      }
    }
    if (length(dfs) == 0) {
      return(NULL)
    }
    # 找交集行并 cbind
    common_genes <- Reduce(intersect, lapply(dfs, rownames))
    if (length(common_genes) == 0) {
      warning("Warning: No common genes found among all PDC code data. Returning NULL.")
      return(NULL)
    }
    for (i in seq_along(dfs)) {
      dfs[[i]] <- dfs[[i]][common_genes, , drop = FALSE]
    }
    merged_df <- do.call(cbind, dfs)
    return(merged_df)
  }

  tumor_df <- read_and_merge_files(data.category, pdc_codes, "tumor")
  normal_df <- read_and_merge_files(data.category, pdc_codes, "normal")

  if (is.null(tumor_df)) {
    stop("Error: No valid tumor data frames were found or no common genes in tumor data.")
  }
  if (is.null(normal_df)) {
    stop("Error: No valid normal data frames were found or no common genes in normal data.")
  }

  # -----------------------
  # 3. 简易差异分析 (t.test + log2FC)
  # -----------------------
  common_genes <- intersect(rownames(tumor_df), rownames(normal_df))
  if (length(common_genes) == 0) {
    stop("Error: No common genes found between tumor and normal data after merging all PDC codes.")
  }

  res_list <- lapply(common_genes, function(g) {
    tumor_vec <- as.numeric(tumor_df[g, ])
    normal_vec <- as.numeric(normal_df[g, ])

    if (all(is.na(tumor_vec)) || all(is.na(normal_vec))) {
      return(c(g, NA, NA, NA))
    }

    tt <- tryCatch(t.test(tumor_vec, normal_vec),
                   error = function(e) NULL)
    if (is.null(tt)) {
      return(c(g, NA, NA, NA))
    }
    pval <- tt$p.value

    mean_tum <- mean(tumor_vec, na.rm = TRUE)
    mean_norm <- mean(normal_vec, na.rm = TRUE)
    if (mean_norm == 0 || is.na(mean_norm)) {
      logfc <- NA
    } else {
      logfc <- log2(mean_tum / mean_norm)
    }

    return(c(g, logfc, pval, -log10(pval)))
  })

  res_df <- data.frame(do.call(rbind, res_list), stringsAsFactors = FALSE)
  colnames(res_df) <- c("gene","log2FC","p.value","negLog10P")
  res_df$log2FC <- as.numeric(res_df$log2FC)
  res_df$p.value <- as.numeric(res_df$p.value)
  res_df$negLog10P <- as.numeric(res_df$negLog10P)

  # 去除 NA
  res_df <- res_df[!(is.na(res_df$log2FC) | is.na(res_df$p.value)), ]

  # -----------------------
  # (1) 在返回表中新增一列 group
  # -----------------------
  # 如果 log2FC < 阈值 => down
  # 如果 log2FC > 阈值 => up
  # 否则 => none
  threshold <- log2FC  # 给定阈值
  res_df$group <- "none"
  res_df$group[res_df$log2FC < -threshold & res_df$p.value < 0.05] <- "down"
  res_df$group[res_df$log2FC > threshold & res_df$p.value < 0.05] <- "up"

  # -----------------------
  # 4. 绘制火山图 (Volcano plot)
  # -----------------------
  # 参考 group 列:
  # up => pink, down => lightblue, none => grey

  # 为了可视化方便，先映射一下颜色
  # 也可以用 scale_color_manual(...)
  res_df$color_group <- ifelse(res_df$group == "up", "pink",
                               ifelse(res_df$group == "down", "lightblue", "grey"))

  # 选 p 值最小的 top.gene 个基因做标注
  # 筛选 p.value < 0.05 且按 log2FC 的绝对值排序
  filtered_res_df <- res_df[res_df$p.value < 0.05, ]
  if (nrow(filtered_res_df) == 0) {
    warning("No genes meet the criteria of p.value < 0.05.")
    label_df <- data.frame()  # 如果没有满足条件的基因，返回空数据框
  } else {
    # 按 log2FC 的绝对值降序排序
    filtered_res_df <- filtered_res_df[order(abs(filtered_res_df$log2FC), decreasing = TRUE), ]
    # 选取前 top.gene 个基因
    label_df <- head(filtered_res_df, top.gene)
  }


  library(ggplot2)
  p <- ggplot(res_df, aes(x = log2FC, y = negLog10P)) +
    geom_point(aes(color = color_group)) +
    # 使用 scale_color_identity 直接使用已有的颜色名称
    scale_color_identity() +
    theme_bw(base_size = 14) +
    labs(
      title = paste0("Volcano Plot for ", cancer.type),
      x = "log2(FC)",
      y = "-log10(p-value)"
    ) +
    # 若想区分 legend，可自行改成 fill=..., or color=... + scale_color_manual()
    theme(legend.position = "none") +
    # 在图中添加两条垂直参考线 (一条在 +threshold, 一条在 -threshold)
    geom_vline(xintercept =  threshold, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = -threshold, linetype = "dashed", color = "grey40")

  # 标注前 top.gene
  p <- p +
    geom_text(data = label_df,
              aes(label = gene),
              color = "black",
              size = 3,
              vjust = -0.5)

  # -----------------------
  # 5. 返回结果
  # -----------------------
  return(invisible(list(
    diff_table = res_df,  # 返回的table里已包含 group 列
    volcano_plot = p
  )))
}

# rna_diff <- cptda(
#   cancer.type = "LUSC",
#   data.category = "Transcriptome",
#   log2FC = 1.5,
#   top.gene = 20
# )
#
# pro_diff <- cptda(
#   cancer.type = "LUSC",
#   data.category = "Proteome",
#   log2FC = 1.5,
#   top.gene = 20
# )
