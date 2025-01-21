.onAttach <- function(libname, pkgname) {
  # 检查数据是否已经存在，如果没有则下载
  project_dir <- getwd()  # 获取当前工作目录
  data_files <- list.files(project_dir, pattern = "\\.csv$")  # 查找是否有CSV文件
  
  if (length(data_files) == 0) {  # 如果没有CSV文件，下载数据
    install_data()
  }
}
