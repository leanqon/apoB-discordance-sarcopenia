# ===================================================================
# NHANES Data Download Script
# ===================================================================
#
# This script helps download the required NHANES XPT files
# Data source: https://www.cdc.gov/nchs/nhanes/
#
# Usage:
#   1. Set your desired download directory below
#   2. Run this script: source("download_nhanes_data.R")
#   3. Files will be downloaded to the specified directory
# ===================================================================

# Set download directory
# NOTE: Update this path to your desired location
download_dir <- "data"  # or use: "path/to/your/data/folder"

# Create directory if it doesn't exist
if(!dir.exists(download_dir)) {
  dir.create(download_dir, recursive = TRUE)
  cat("Created directory:", download_dir, "\n")
}

# Base URL for NHANES data
base_url <- "https://wwwn.cdc.gov/Nchs/Nhanes"

# Define required data files for each cycle
data_files <- list(
  # 2011-2012 (G)
  "2011-2012" = c(
    "DEMO_G", "BMX_G", "DXX_G", "APOB_G", "TCHOL_G", "HDL_G",
    "TRIGLY_G", "GLU_G", "GHB_G", "BIOPRO_G", "DIQ_G", "BPQ_G",
    "MCQ_G", "SMQ_G", "ALQ_G", "PAQ_G"
  ),

  # 2013-2014 (H)
  "2013-2014" = c(
    "DEMO_H", "BMX_H", "DXX_H", "APOB_H", "TCHOL_H", "HDL_H",
    "TRIGLY_H", "GLU_H", "INS_H", "GHB_H", "BIOPRO_H", "DIQ_H",
    "BPQ_H", "MCQ_H", "SMQ_H", "ALQ_H", "PAQ_H"
  ),

  # 2015-2016 (I)
  "2015-2016" = c(
    "DEMO_I", "BMX_I", "DXX_I", "APOB_I", "TCHOL_I", "HDL_I",
    "TRIGLY_I", "GLU_I", "INS_I", "GHB_I", "BIOPRO_I", "DIQ_I",
    "BPQ_I", "MCQ_I", "SMQ_I", "ALQ_I", "PAQ_I"
  )
)

#' Download a single NHANES file
#' @param file_name Name of the XPT file (without extension)
#' @param cycle NHANES cycle (e.g., "2011-2012")
#' @param dest_dir Destination directory
download_nhanes_file <- function(file_name, cycle, dest_dir) {

  # Construct URL
  url <- paste0(base_url, "/", cycle, "/", file_name, ".XPT")

  # Construct destination path
  dest_file <- file.path(dest_dir, paste0(file_name, ".xpt"))

  # Skip if file already exists
  if(file.exists(dest_file)) {
    cat("  [SKIP] File already exists:", file_name, "\n")
    return(TRUE)
  }

  # Download file
  cat("  [DOWNLOAD]", file_name, "... ")

  tryCatch({
    download.file(url, dest_file, mode = "wb", quiet = TRUE)
    cat("SUCCESS\n")
    return(TRUE)
  }, error = function(e) {
    cat("FAILED -", e$message, "\n")
    return(FALSE)
  })
}

# ===================================================================
# Download all files
# ===================================================================

cat("\n=== NHANES Data Download ===\n")
cat("Download directory:", download_dir, "\n\n")

total_files <- 0
success_files <- 0
failed_files <- c()

for(cycle in names(data_files)) {
  cat("\nDownloading", cycle, "data:\n")
  cat(rep("-", 50), "\n", sep="")

  files <- data_files[[cycle]]

  for(file_name in files) {
    total_files <- total_files + 1

    result <- download_nhanes_file(file_name, cycle, download_dir)

    if(result) {
      success_files <- success_files + 1
    } else {
      failed_files <- c(failed_files, paste0(file_name, " (", cycle, ")"))
    }

    # Small delay to avoid overwhelming the server
    Sys.sleep(0.5)
  }
}

# ===================================================================
# Summary
# ===================================================================

cat("\n", rep("=", 60), "\n", sep="")
cat("=== Download Summary ===\n")
cat(rep("=", 60), "\n", sep="")
cat("\nTotal files:", total_files)
cat("\nSuccessful:", success_files)
cat("\nFailed:", length(failed_files), "\n")

if(length(failed_files) > 0) {
  cat("\nFailed files:\n")
  for(f in failed_files) {
    cat("  -", f, "\n")
  }
  cat("\nNote: Failed files may need to be downloaded manually from:")
  cat("\nhttps://wwwn.cdc.gov/nchs/nhanes/Default.aspx\n")
} else {
  cat("\n✓ All files downloaded successfully!\n")
  cat("\nYou can now run the analysis script:\n")
  cat("  source('analysis_code_github.R')\n")
}

cat("\n", rep("=", 60), "\n", sep="")

# ===================================================================
# Verify downloaded files
# ===================================================================

cat("\n=== Verifying Files ===\n")

all_files <- unlist(data_files)
existing_files <- list.files(download_dir, pattern = "\\.xpt$", ignore.case = TRUE)

cat("Expected files:", length(all_files), "\n")
cat("Found files:", length(existing_files), "\n")

missing_files <- setdiff(
  tolower(paste0(all_files, ".xpt")),
  tolower(existing_files)
)

if(length(missing_files) > 0) {
  cat("\n⚠ Missing files:\n")
  for(f in missing_files) {
    cat("  -", f, "\n")
  }
} else {
  cat("\n✓ All required files are present!\n")
}

cat("\nDownload complete!\n")
