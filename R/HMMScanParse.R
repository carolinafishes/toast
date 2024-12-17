HMMScanParse<-function(input_file="hmmscan.out", outputfile=TRUE, user_threshold=1e-5, fasta_file="extracted_sequences.fasta", clean=TRUE){
#Parse HMMScan output
# Define output file
output_file <- sub("\\.out$", "_queries.csv", input_file)

# Read the input file
lines <- readLines(input_file)

# Extract indices of lines starting with "Query:"
query_indices <- grep("^Query:", lines)

# Extract the query names using a regular expression
query_names <- sub("^Query:\\s+([^\\s]+)\\s+\\[L=.*$", "\\1", lines[query_indices])

# Initialize results data frame
results <- data.frame(Query = character(), Domain = character(), stringsAsFactors = FALSE)

# Loop through each query and extract associated domains
for (i in seq_along(query_indices)) {
  # Extract the current query name
  query_name <- query_names[i]
  query_line_index <- query_indices[i]
  
  # Find the "Scores for complete sequence" section
  scores_index <- grep("Scores for complete sequence", lines[(query_line_index + 1):length(lines)])
  
  if (length(scores_index) > 0) {
    # Adjust the scores_index relative to the full file
    scores_index <- scores_index[1] + query_line_index
    
    # Extract table rows starting 3 lines below "Scores for complete sequence"
    table_start <- scores_index + 3
    table_lines <- lines[table_start:length(lines)]
    
    # Stop at the "------ inclusion threshold ------" line or an empty line
    stop_index <- which(grepl("^\\s*-+\\s*inclusion threshold", table_lines) | table_lines == "")
    if (length(stop_index) > 0) {
      table_lines <- table_lines[1:(stop_index[1] - 1)]
    }
    
    # Process each valid table row
    for (line in table_lines) {
      # Only process lines that resemble table rows
      if (grepl("^\\s*[0-9eE.+-]+\\s", line)) {  # Starts with numeric E-value
        # Split by 2 or more spaces to get columns
        parts <- unlist(strsplit(line, "\\s{2,}"))
        
        if (length(parts) >= 2) {
          # Extract the description (last column)
          domain_desc <- parts[length(parts)]
          evalue<-parts[2]
          model<-parts[length(parts)-1]
          results <- rbind(results, data.frame(Query = query_name, evalue= evalue, Model=model, Domain = domain_desc, stringsAsFactors = FALSE))
        }
      }
    }
  } else {
    cat("No scores table found for query:", query_name, "\n")  # Debugging message
  }
}
#cleanup
results <- results[!grepl("-{2,}", results$Domain), ]
results <- results[as.numeric(results$evalue) <= user_threshold, ]

  # Add columns for alifrom and ali to
  results$alifrom <- NA
  results$alito <- NA
  
  # Loop through each query in the results
  for (i in seq_len(nrow(results))) {
    query <- results$Query[i]
    model <- results$Model[i]
    
    # Find the section for the query in the hmmscan file
    query_line_index <- grep(paste0("^Query:\\s+", query), lines)
    
    if (length(query_line_index) > 0) {
      # Extract the relevant block for the query
      query_lines <- lines[(query_line_index + 1):length(lines)]
      
      # Find the model within this query's block
      model_index <- grep(paste0("^>>\\s+", model, "\\s"), query_lines)
      
      if (length(model_index) > 0) {
        # Extract the relevant block for the model
        model_lines <- query_lines[(model_index[1] + 1):length(query_lines)]
        domain_lines <- model_lines[grep("^\\s*\\d+\\s+.+\\s+.+\\s+.+\\s+.+\\s+.+\\s+\\d+\\s+\\d+\\s+", model_lines)]
        
        if (length(domain_lines) > 0) {
          # Extract the alifrom and ali to from the first matching line
          line_parts <- unlist(strsplit(domain_lines[1], "\\s+"))
          alifrom <- as.numeric(line_parts[11])
          alito <- as.numeric(line_parts[12])
          
          # Update the results dataframe
          results$alifrom[i] <- alifrom
          results$alito[i] <- alito
        }
      }
    }
  }

# Trim whitespace from the Query column 
results$Query <- trimws(results$Query)
results<-add_sequence_column_protein(results, fasta_file)
#remove duplicate domains with minor wobble in start and stop
if (clean==TRUE){
results<-filter_overlapping_domains(results)
}
if (outputfile ==TRUE){
write.csv(results, output_file, row.names = FALSE, quote=FALSE)
}
#return result
return(results)
}


