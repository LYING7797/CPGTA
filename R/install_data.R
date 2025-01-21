install_data <- function() {
# 安装并加载 googledrive 包
  if (!requireNamespace("googledrive", quietly = TRUE)) {
    install.packages("googledrive")
  }
  library(googledrive)

  # Google Drive文件夹的ID
  folder_id <- "1Yr_msBFVwu26F4MwGBgbAu0255vb02Eq"
  
  # 获取文件夹中的所有文件
  files <- drive_find(path = as_id(folder_id))  # 不限制最大文件数，获取文件夹中的所有文件
  
  # 确保安装目录存在
  project_dir <- getwd()  # 获取当前工作目录，即包的根目录
  
  # 下载文件到包的根目录
  for (file in files$name) {
    file_path <- file.path(project_dir, file)  # 文件的本地存储路径，放在包的根目录下
    drive_download(file, path = file_path, overwrite = TRUE)
  }
}
