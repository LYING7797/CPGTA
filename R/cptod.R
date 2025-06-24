
#' @description 该函数用于从数据目录中获取癌症蛋白质组学数据，可以根据癌症类型、
#' 数据类别、数据类型和样本类型等参数筛选并下载相应的数据。
#'
#' @param cancer.type 字符串。癌症类型名称或缩写。
#' @param data.category 字符串。数据类别，可选值包括"Biospecimen", "Clinical data",
#'   "Copy Number Variation", "Phosphoproteome", "Proteome",
#'   "Simple Nucleotide Variation", "Transcriptome"。
#' @param data.type 字符串。数据类型，一般为"Normalized"("nor")或"Tidied"("tidy")。对于转录组数据，
#'   可以是"TPM"或"Counts"。默认为NULL，表示不限制数据类型。
#' @param sample.type 字符串或字符串向量。样本类型，如"tumor"或"normal"。
#'   默认为NULL，表示不限制样本类型。
#' @param PDC.study.identifier 字符串。特定的PDC研究标识符(如"PDC000125")。
#'   默认为NULL，表示根据癌症类型查找相关的所有PDC研究。
#'
#' @return 返回一个列表，包含根据参数筛选的数据文件。
#'

cptod <- function(cancer.type = NULL, data.category, data.type = NULL,
                  sample.type = NULL, PDC.study.identifier = NULL) {
  # 基础目录路径
  base_dir <- "./data1/"

  # 验证必填参数
  if (is.null(cancer.type) && is.null(PDC.study.identifier)) {
    stop("Either cancer.type or PDC.study.identifier must be provided.")
  }

  # 定义有效的数据类别
  valid_categories <- c("Biospecimen", "Clinical data", "Copy Number Variation",
                        "Phosphoproteome", "Proteome", "Simple Nucleotide Variation",
                        "Transcriptome")

  # 验证数据类别参数
  if (!data.category %in% valid_categories) {
    stop(paste("Invalid data category. Valid options include:", paste(valid_categories, collapse = ", ")))
  }

  # 加载癌症-PDC映射信息
  cancer_pdc_info_path <- file.path(base_dir, "cancer_PDC_info.csv")
  if (!file.exists(cancer_pdc_info_path)) {
    stop("The file cancer_PDC_info.csv was not found, unable to retrieve the mapping between cancer types and PDC studies.")
  }
  cancer_pdc_info <- read.csv("./data1/cancer_PDC_info.csv",check.names = FALSE)

  # 验证cancer.type参数
  # 获取唯一的癌症类型和缩写
  valid_cancer_types <- unique(data.frame(
    cancer_type = cancer_pdc_info$cancer_type,
    abbreviation = cancer_pdc_info$abbreviation
  ))
  if (!is.null(cancer.type) && !(cancer.type %in% unlist(valid_cancer_types))) {
    # 格式化输出为两列，确保对齐
    formatted_types <- apply(valid_cancer_types, 1, function(row) {
      paste(format(row[1], width = 40, justify = "left"),  # cancer_type 左对齐，宽度40
            format(row[2], width = 10, justify = "left"))  # abbreviation 左对齐，宽度10
    })

    # 添加标题行
    header <- paste(format("cancer_type", width = 40, justify = "left"),
                    format("abbreviation", width = 10, justify = "left"))

    # 拼接标题和内容
    valid_types_message <- paste(c(header, formatted_types), collapse = "\n")

    # 输出错误信息
    stop(paste("Invalid cancer type. Valid options include:\n\n", valid_types_message))
  }

  # 将缩写转换为完整癌症类型名称
  if (!is.null(cancer.type) && cancer.type %in% cancer_pdc_info$abbreviation) {
    cancer.type <- unique(cancer_pdc_info[cancer_pdc_info$abbreviation == cancer.type,]$cancer_type)
  }

  # 确定PDC研究ID
  if (!is.null(PDC.study.identifier)) {
    # 验证PDC标识符是否在对应data.category中有效
    pdc_column <- data.category
    if (!PDC.study.identifier %in% cancer_pdc_info[[pdc_column]]) {
      stop(paste("The provided PDC.study.identifier is invalid in the", data.category))
    }

    # 验证PDC.study.identifier是否与cancer.type匹配
    if (!is.null(cancer.type)) {
      matching_rows <- cancer_pdc_info[cancer_pdc_info[[pdc_column]] == PDC.study.identifier, ]
      if (!(cancer.type %in% matching_rows$cancer_type)) {
        stop("The provided PDC.study.identifier does not match the specified cancer.type.")
      }
    }

    pdc_ids <- PDC.study.identifier
  } else {
    # 根据cancer.type查找所有匹配的PDC ID
    matching_rows <- cancer_pdc_info[cancer_pdc_info$cancer_type == cancer.type, ]

    if (nrow(matching_rows) == 0) {
      stop("The specified cancer type was not found. Please check the cancer_PDC_info.csv file for valid cancer types.")
    }

    # 获取对应data.category的PDC ID
    pdc_column <- data.category
    pdc_ids <- matching_rows[[pdc_column]]
    pdc_ids <- pdc_ids[!is.na(pdc_ids)]

    if (length(pdc_ids) == 0) {
      stop(paste("No PDC studies were found for the specified cancer type ", cancer.type, "under", data.category))
    }
  }

  # 处理Biospecimen数据
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

  # 处理Clinical data
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

  # 处理Simple Nucleotide Variation数据
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

  # 处理其他数据类别(Proteome, Phosphoproteome, Copy Number Variation, Transcriptome)

  # 确定要获取的数据类型
  valid_data_types <- c("Normalized", "Tidied")
  if (data.category == "Transcriptome") {
    valid_data_types <- c(valid_data_types, "counts")
  }

  # 如果data.type未指定，使用所有有效的数据类型
  if (is.null(data.type)) {
    data.type <- valid_data_types
  } else {
    # 验证data.type
    if (!all(data.type %in% valid_data_types)) {
      stop(paste("Invalid data type. For the", data.category, " valid options:", paste(valid_data_types, collapse = ", ")))
    }
  }

  results <- list()

  # 处理转录组的counts数据
  if (data.category == "Transcriptome" && "counts" %in% data.type) {
    counts_dir <- file.path(base_dir, "Transcriptome", "Counts_data")
    if (!dir.exists(counts_dir)) {
      warning("The Counts_data directory was not found.")
    } else {
      for (pdc_id in pdc_ids) {
        # 构建文件模式
        file_pattern <- paste0("^", pdc_id)
        if (!is.null(sample.type)) {
          sample_patterns <- paste0("_", sample.type, "_", collapse = "|")
          file_pattern <- paste0(file_pattern, ".*?(", sample_patterns, ")")
        }
        file_pattern <- paste0(file_pattern, ".*\\.csv$")

        # 查找匹配的文件
        files <- list.files(counts_dir, pattern = file_pattern, full.names = TRUE)

        for (file in files) {
          file_name <- basename(file)
          results[[file_name]] <- read.csv(file)
        }
      }
    }
  }

  # 处理Normalized和tidied数据
  for (dt in data.type) {
    if (dt == "counts") next  # 已经处理过counts数据

    # 确定文件名后缀
    suffix <- if (dt == "Normalized") "nor.csv" else "tidy.csv"

    # 确定文件名前缀
    prefix <- switch(data.category,
                     "Proteome" = "pro",
                     "Phosphoproteome" = "phos",
                     "Copy Number Variation" = "cnv",
                     "Transcriptome" = "rna",
                     stop(paste("Unsupported data category:", data.category)))

    for (pdc_id in pdc_ids) {
      # 构建文件模式
      if (!is.null(sample.type)) {
        # 有sample.type限制
        for (st in sample.type) {
          file_pattern <- paste0("^", pdc_id, "_", prefix, "_", st, "_.*", suffix, "$")
          files <- list.files(file.path(base_dir, data.category), pattern = file_pattern, full.names = TRUE)

          for (file in files) {
            file_name <- basename(file)
            results[[file_name]] <- read.csv(file)
          }
        }
      } else {
        # 无sample.type限制
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








