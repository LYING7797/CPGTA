.onAttach <- function(libname, pkgname) {
    # 自动下载数据文件
  install_data()  # 调用 install_data 函数下载数据
}
