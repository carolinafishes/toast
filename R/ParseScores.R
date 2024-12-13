#' Scrutinize missing data patterns in your sequence data
#'
#' Get scores from BUSCO summary file
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny fasta ortholog busco
#' @param input_id the name of your sample, can be anything “Genus_species”
#' @param file_path the location of your busco "short summary" file
#' @return Returns a dataframe of the summary scores
#' @export
#' @examples
#' ParseScores()

ParseScores<- function(file_path, input_id) {
  # Read the file
  lines <- readLines(file_path)
  
  # Find the line containing the Results summary
  results_line <- grep("Results:", lines, value = TRUE)
  
  # Extract the percentages from the Results line
  if (length(results_line) > 0) {
    percent_line <- lines[grep("\\%\\[", lines)]
    percentages <- regmatches(percent_line, gregexpr("[0-9]+\\.[0-9]+", percent_line))[[1]]
    
    # Create a dataframe with the user input ID and percentages
    df <- data.frame(
      ID = input_id,
      Complete = as.numeric(percentages[1]),
      Complete_Single_Copy = as.numeric(percentages[2]),
      Complete_Duplicated = as.numeric(percentages[3]),
      Fragmented = as.numeric(percentages[4]),
      Missing = as.numeric(percentages[5]),
      stringsAsFactors = FALSE
    )
    return(df)
  } else {
    stop("Results line not found in the file.")
  }
}