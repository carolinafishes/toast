#' Assess completeness scores from a multisample busco run
#'
#' Get scores from BUSCO summary files. Function assumes a parent directory, and each run as a directory within that directory
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny fasta ortholog busco
#' @param folder_path the location of your busco "short summary" file
#' @return Returns a dataframe of the summary scores, using each directory as an input ID
#' @export
#' @examples
#' compile_busco_summaries()

# Main function to process all directories
compile_busco_summaries <- function(folder_path) {
  # List all directories in the folder
  directories <- list.dirs(folder_path, full.names = TRUE, recursive = FALSE)
  
  # Initialize an empty list to store dataframes
  results_list <- list()
  
  # Loop over directories
  for (dir in directories) {
    # Extract the directory name
    dir_name <- basename(dir)
    
    # Find the short_summary*.txt file in the directory
    summary_file <- list.files(dir, pattern = "^short_summary.*\\.txt$", full.names = TRUE)
    
    # Ensure there's exactly one file
    if (length(summary_file) == 1) {
      # Parse the file and append the result
      results_list[[dir_name]] <- ParseScores(summary_file, dir_name)
    } else {
      warning(paste("No unique short_summary*.txt file found in:", dir))
    }
  }
  
  # Combine all dataframes into one
  final_df <- bind_rows(results_list)
  return(final_df)
}

