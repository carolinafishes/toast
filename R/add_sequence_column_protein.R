add_sequence_column_protein <- function(results, fasta_file) {
  # Load the target FASTA file as protein sequences
  fasta_sequences <- tryCatch(
    readAAStringSet(fasta_file, format = "fasta", use.names = TRUE),
    error = function(e) stop("Error reading FASTA file: ", e$message)
  )
  
  # Print loaded sequences for debugging
  cat("Loaded sequences:\n")
  print(names(fasta_sequences))
  
  # Create a new column for the sequence
  results$sequence <- NA
  
  # Loop through each row in the results
  for (i in seq_len(nrow(results))) {
    query <- results$Query[i]
    alifrom <- results$alifrom[i]
    alito <- results$alito[i]
    
    # Check if alifrom and alito are valid
    if (!is.na(alifrom) && !is.na(alito)) {
      # Find the corresponding sequence in the FASTA file
      fasta_query <- fasta_sequences[names(fasta_sequences) == query]
      
      if (length(fasta_query) > 0) {
        # Extract the subsequence using alifrom and alito
        subseq <- subseq(fasta_query, start = alifrom, end = alito)
        results$sequence[i] <- as.character(subseq)
      } else {
        cat("Query", query, "not found in the FASTA file.\n")
      }
    } else {
      cat("Invalid alifrom or alito for Query", query, "\n")
    }
  }
  return(results)
}