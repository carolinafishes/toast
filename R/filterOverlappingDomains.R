filter_overlapping_domains <- function(results) {
  # Ensure numeric values
  results$alifrom <- as.numeric(results$alifrom)
  results$alito <- as.numeric(results$alito)
  results$evalue <- as.numeric(results$evalue)
  
  # Split the results by Query
  results_split <- split(results, results$Query)
  
  # Process each Query independently
  results_with_keep <- lapply(results_split, function(query_df) {
    # Sort by alifrom and then by evalue (ascending)
    query_df <- query_df[order(query_df$alifrom, query_df$evalue), ]
    
    # Initialize a logical column to mark rows to keep
    query_df$keep <- rep(TRUE, nrow(query_df))
    
    # Compare adjacent rows for overlap
    for (i in seq_len(nrow(query_df) - 1)) {
      if (query_df$alito[i] >= query_df$alifrom[i + 1]) {  # Overlapping ranges
        query_df$keep[i + 1] <- FALSE  # Mark the next row with higher evalue as FALSE
      }
    }
    
    return(query_df)
  })
  
  # Combine the results with `keep` column into a single dataframe
  results_with_keep <-do.call(rbind, results_with_keep)
  final_results <- results_with_keep[results_with_keep$keep, ]
  final_results <- subset(final_results, select = -keep)
}