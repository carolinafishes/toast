HMMParse <- function(hmmer_output, output_csv="HMMSummary", evalue_threshold=1e-5,type = "AA",extract=TRUE) {
  # Ensure the HMMER output file exists
  if (!file.exists(hmmer_output)) {
    stop(paste("Error: HMMER output file", hmmer_output, "not found!"))
  }
  
  # Read the HMMER output file
  lines <- readLines(hmmer_output)
  
  # Extract the target file name
  target_line <- grep("target sequence database:", lines, value = TRUE)
  if (length(target_line) == 0) stop("Error: 'target sequence database' not found in the HMMER output.")
  target_file <- trimws(sub(".*target sequence database:\\s+", "", target_line))
  
  # Extract accession (everything before the first '.')
  accession <- sub("\\..*", "", basename(target_file))
  
  # Identify lines that likely match the hit format (more specific matching)
  hit_lines <- grep("^\\s*[0-9.eE+-]+\\s+\\S+.*\\s+\\S+\\s+\\S+\\s+\\S+$", lines, value = TRUE)
  
  # Initialize an empty list for results
  results <- list()
  
  # Parse valid hit lines into a data frame
  for (line in hit_lines) {
    fields <- unlist(strsplit(trimws(line), "\\s+"))
    
    # Validate and parse the line
    if (length(fields) >= 9) {
      evalue <- suppressWarnings(as.numeric(fields[1]))
      sequence <- fields[9]
      
      # Only include hits with E-value below the threshold
      if (!is.na(evalue) && evalue <= evalue_threshold) {
        results[[length(results) + 1]] <- data.frame(
          Filename = target_file,
          Accession = accession,
          hits = sequence,
          Evalue = evalue,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # Combine all valid results into a single data frame
  results_df <- do.call(rbind, results)
  
  # Write results to CSV if any hits exist
  if (!is.null(results_df) && nrow(results_df) > 0) {
    write.csv(results_df, file = output_csv, row.names = FALSE,quote=FALSE)
    message(paste("CSV file created:", output_csv))
  } else {
    warning("No hits found below the specified E-value threshold.")
  }
  if (extract==TRUE){
  extractSequences(results_df,type=type, external=FALSE)
  }
}