cptod <- function(cancer.type,
                       data.category,
                       PDC.study.identifier = NULL,
                       data.type = NULL,
                       sample.type = NULL) {
  # 首先加载dplyr库
  library(dplyr)

  # 定义可接受的癌症类型范围
  cancer_pdc_info <- read.csv("./data1/cancer-PDC info.csv")
  valid_cancer_types <- unique(c(cancer_pdc_info$cancer_type, cancer_pdc_info$abbreviation))
  valid_data_categories <- c("Biospecimen", "Clinical data", "Copy Number Variation","Transcriptome", "Proteome", "Phosphoproteome")
  valid_study_identifiers <- unique(c(cancer_pdc_info$PDC_Pro, cancer_pdc_info$PDC_Phospro))
  # 定义可接受的数据类型和样本类型范围
  valid_data_types <- c('Harmonized',"Normalized")
  valid_sample_types <- c("Primary tumor", "Normal tissue")


  # Validate if cancer.type is within the valid range
  if (is.null(cancer.type) ||!(cancer.type %in% valid_cancer_types)) {
    stop(paste("Error: 'cancer.type' must be provided and 'cancer.type' is not in the valid range. Valid values are:", paste(valid_cancer_types, collapse = "\n ")))
  }

  if(cancer.type %in% cancer_pdc_info$abbreviation){
    cancer.type = unique(cancer_pdc_info[cancer_pdc_info$abbreviation == cancer.type,]$cancer_type)
  }


  if (!(data.category %in% valid_data_categories)) {
    stop(paste("Error: 'data.category' must be provided and 'data.category' is not in the valid range. Valid values are:", paste(valid_data_categories, collapse = "\n ")))
  }

  base_path <- file.path("./data/clinical tumor sample", cancer.type, data.category)

  # "Biospecimen" 和 "Clinical data" 的特殊逻辑
  if (data.category %in% c("Biospecimen", "Clinical data")) {
    # 检查 PDC.study.identifier 并禁止 data.type 和 sample.type
    if (!is.null(data.type) || !is.null(sample.type)) {
      stop("Error: 'data.type' and 'sample.type' should not be provided when data.category is 'Biospecimen' or 'Clinical data'.")
    }

    # 验证 PDC.study.identifier 是否匹配 cancer.type
    if (!is.null(PDC.study.identifier)) {

      matching_rows <- cancer_pdc_info[cancer_pdc_info$PDC_Pro == PDC.study.identifier, ]

      # 验证 PDC.study.identifier 是否与 cancer.type 对应
      if (!(cancer.type %in% matching_rows$cancer_type)) {
        stop("Error: Provided PDC.study.identifier does not match the specified cancer.type.")
      }

      # 构造 PDC.study.identifier 对应的目录路径
      directory_path <- file.path(base_path, PDC.study.identifier)

      # 检查目录是否存在
      if (!dir.exists(directory_path)) {
        stop("Error: Directory does not exist for the specified PDC.study.identifier.")
      }

      # 读取该 PDC 编号目录下的所有文件
      files <- list.files(directory_path, recursive = TRUE, full.names = TRUE)
      if (length(files) == 0) {
        stop("No files found in the specified directory.")
      }

      # 将文件路径的上层文件夹名称和文件名设置为元素名
      data_list <- setNames(
        lapply(files, function(f) {
          read.delim(f, sep = ",", header = TRUE, check.names = FALSE)
        }),
        sapply(files, function(file) {
          paste0(basename(dirname(file)), "/", basename(file))
        })
      )
    } else {
      # 遍历所有 PDC 编号的文件夹
      study_identifier_paths <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)
      data_list <- list()
      for (study_path in study_identifier_paths) {
        files <- list.files(study_path, recursive = TRUE, full.names = TRUE)
        if (length(files) > 0) {
          study_data <- setNames(
            lapply(files, function(f) {
              read.delim(f, sep = ",", header = TRUE, check.names = FALSE)
            }),
            sapply(files, function(file) {
              paste0(basename(dirname(file)), "/", basename(file))
            })
          )
          data_list <- c(data_list, study_data)
        }
      }
    }

    return(list(data_list = data_list))
  }

  # 其他 data.category 的处理逻辑
  if (!(data.category %in% c("Biospecimen", "Clinical data"))) {
    if (is.null(data.type) || !all(data.type %in% valid_data_types)) {
      stop("Error: data.type must be provided and must be valid for this data.category. Valid values are: Harmonized, Normalized")
    }
    if (is.null(sample.type) || !all(sample.type %in% valid_sample_types)) {
      stop("Error: sample.type must be provided and must be valid for this data.category. Valid values are: Primary tumor, Normal tissue")
    }

    data_list <- list()  # 初始化结果列表

    if (!is.null(PDC.study.identifier)) {
      # 确定有效的 study_identifier
      if (!is.null(data.category) & !is.null(cancer.type)) {
        # 检查 PDC.study.identifier 是否在有效范围内
        if (!(PDC.study.identifier %in% valid_study_identifiers)) {
          stop("Error: Provided PDC.study.identifier is not valid for the specified data.category and cancer.type.")
        }
      }

      # 如果 PDC.study.identifier 合理，继续处理文件
      matching_rows <- cancer_pdc_info[cancer_pdc_info$PDC_Pro == PDC.study.identifier, ]

      if (!(cancer.type %in% matching_rows$cancer_type)) {
        stop("Error: Provided PDC.study.identifier does not match the specified cancer.type.")
      }
      for (dtype in data.type) {
        for (stype in sample.type) {
          directory_path <- file.path(base_path, PDC.study.identifier, dtype, stype)
          if (dir.exists(directory_path)) {
            files <- list.files(directory_path, recursive = TRUE, full.names = TRUE)
            if (length(files) > 0) {
              study_data <- setNames(
                lapply(files, function(f) {
                  read.delim(f, sep = ",", header = TRUE, check.names = FALSE)
                }),
                sapply(files, function(file) {
                  paste0(basename(dirname(file)), "/", basename(file))
                })
              )
              data_list <- c(data_list, study_data)
            }
          }
        }
      }
    } else {
      study_identifier_paths <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)
      for (study_path in study_identifier_paths) {
        for (dtype in data.type) {
          for (stype in sample.type) {
            type_path <- file.path(study_path, dtype, stype)
            if (dir.exists(type_path)) {
              files <- list.files(type_path, recursive = TRUE, full.names = TRUE)
              if (length(files) > 0) {
                study_data <- setNames(
                  lapply(files, function(f) {
                    read.delim(f, sep = ",", header = TRUE, check.names = FALSE)
                  }),
                  sapply(files, function(file) {
                    paste0(basename(dirname(file)), "/", basename(file))
                  })
                )
                data_list <- c(data_list, study_data)
              }
            }
          }
        }
      }
    }

    if (length(data_list) == 0) {
      message(paste("No files found: This cancer.type has no", data.category, "PDC.study.identifier with",
                    paste(data.type, collapse = " and "), "and", paste(sample.type, collapse = " and "), "data files."))
    }

    return(data_list)
  }
}




