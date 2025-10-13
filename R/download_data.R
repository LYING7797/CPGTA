#' Recursively download all files from a Google Drive folder and keep the original folder structure
#'
#' @param folder_id Google Drive folder ID to download (e.g., "1Yr_msBFVwu26F4MwGBgbAu0255vb02Eq")
#' @param local_dir Local directory to save the files. Default is current working directory.
#' @examples
#' download_data("1Yr_msBFVwu26F4MwGBgbAu0255vb02Eq")
#' @import googledrive
#' @export


download_data <- function(folder_id) {
  # Ensure googledrive is installed and loaded
  if (!requireNamespace("googledrive", quietly = TRUE)) {
    install.packages("googledrive")
  }
  library(googledrive)
  
  # Get CPGTA package directory and set download location to data1 directory
  pkg_dir <- find.package("CPGTA")
  local_dir <- file.path(pkg_dir)
  
  # Create data1 directory if it doesn't exist
  if (!dir.exists(local_dir)) {
    dir.create(local_dir, recursive = TRUE)
  }
  
  # Inform the user where data will be downloaded
  message("Downloading data to: ", local_dir)

  # Helper function: recursively download
  download_recursive <- function(current_folder_id, current_local_dir) {
    # List all items in the current folder
    items <- drive_ls(as_id(current_folder_id))
    for (i in seq_len(nrow(items))) {
      item <- items[i, ]
      item_name <- item$name
      item_id <- item$id
      item_type <- item$drive_resource[[1]]$mimeType

      # If item is a folder
      if (item_type == "application/vnd.google-apps.folder") {
        # Create local folder if not exists
        sub_dir <- file.path(current_local_dir, item_name)
        if (!dir.exists(sub_dir)) dir.create(sub_dir, recursive = TRUE)
        # Recursive call
        download_recursive(item_id, sub_dir)
      } else {
        # It's a file, download it
        dest_path <- file.path(current_local_dir, item_name)
        drive_download(as_id(item_id), path = dest_path, overwrite = TRUE)
        message("Downloaded: ", item_name)
      }
    }
  }

  # Start recursive download from the root folder
  tryCatch({
    download_recursive(folder_id, local_dir)
    message("Download completed successfully!")
    return(local_dir) # Return the path where data was downloaded
  }, error = function(e) {
    message("Error during download: ", e$message)
    return(NULL)
  })
}
