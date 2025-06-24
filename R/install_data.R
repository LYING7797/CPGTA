install_data <- function() {
  # 加载googledrive包
  if (!requireNamespace("googledrive", quietly = TRUE)) {
    install.packages("googledrive")
  }
  library(googledrive)

  # Google Drive文件夹的ID
  folder_id <- "1Yr_msBFVwu26F4MwGBgbAu0255vb02Eq"

  # 获取文件夹中的所有文件
  files <- drive_find(n_max = 10, path = as_id(folder_id))

  # 下载文件到本地目录（确保目录存在）
  for (file in files$name) {
    # 获取文件的路径
    file_info <- drive_get(file)
    file_path <- file_info$path

    # 确保本地目录存在
    local_dir <- dirname(file_path)
    if (!dir.exists(local_dir)) {
      dir.create(local_dir, recursive = TRUE)
    }

    # 下载文件到本地目录
    drive_download(file, path = file_path, overwrite = TRUE)
  }
}
