install_data <- function() {
  # 安装并加载 googledrive 包
  if (!requireNamespace("googledrive", quietly = TRUE)) {
    install.packages("googledrive")
  }
  library(googledrive)

  # Google Drive文件夹的ID
  folder_ids <- list(
    data = "1Vl5YwTx1q7tTx9RHdNMqnQptH8SNUqm4",  # data文件夹ID
    data1 = "1sExi5ZjhfP8QwD9tVvBPLfdaI4HmGCMF"  # data1文件夹ID
  )
  
  # 遍历每个文件夹并下载文件
  for (folder in names(folder_ids)) {
    folder_id <- folder_ids[[folder]]
    
    # 获取文件夹中的所有文件
    files <- drive_find(path = as_id(folder_id))
    
    # 确保数据目录存在
    data_dir <- file.path(getwd(), folder)  # 目标目录为data1和data
    if (!dir.exists(data_dir)) {
      dir.create(data_dir)
    }
    
    # 下载文件到本地目录
    for (file in files$name) {
      file_path <- file.path(data_dir, file)  # 下载路径
      drive_download(file, path = file_path, overwrite = TRUE)
    }
  }
}

