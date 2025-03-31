HMMScanParse3.1 <- function(input_file, user_threshold = 1e-5, outputfile = TRUE, fasta_file = "extracted_sequences.fasta", clean = TRUE) {
  
  lines <- readLines(input_file)
  
  results <- data.frame(
    Query = character(),
    Model = character(),
    Domain = character(),
    evalue = numeric(),
    alifrom = integer(),
    alito = integer(),
    sequence = character(),
    stringsAsFactors = FALSE
  )
  
  current_query <- NA
  current_model <- NA
  current_domain <- NA
  
  for (i in seq_along(lines)) {
    
    # ✅ Extract the full query name
    if (grepl("^Query:", lines[i])) {
      current_query <- sub("^Query:\\s+(.+)", "\\1", lines[i])
      current_query <- trimws(current_query)
      next
    }
    
    # ✅ Match model and domain after ">>"
    if (grepl("^>>", lines[i]) && !is.na(current_query)) {
      current_model <- sub("^>>\\s+([^\\s]+).*", "\\1", lines[i])
      current_domain <- sub("^>>\\s+[^\\s]+\\s+(.+)", "\\1", lines[i])
      current_model <- trimws(current_model)
      # 🚀 Remove ">>" from Domain name and clean up
      current_domain <- gsub("^>>\\s*", "", current_domain)
      current_domain <- trimws(current_domain)
      # 🚫 Remove commas from Domain (to avoid CSV issues)
      current_domain <- gsub(",", "", current_domain)
      next
    }
    
    # ✅ Extract numeric data from valid result lines using REGEX
    match <- regexec(
      "^\\s*\\d+\\s+!\\s+[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?\\s+[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?\\s+([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\\s+[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?\\s+\\d+\\s+\\d+\\s+\\.{0,2}\\s+(\\d+)\\s+(\\d+)\\s+",
      lines[i]
    )
    
    result <- regmatches(lines[i], match)
    
    if (length(result) > 0 && length(result[[1]]) == 4) {
      evalue <- as.numeric(result[[1]][2])   # Evalue from the regex group
      alifrom <- as.integer(result[[1]][3])  # alifrom from regex group
      alito <- as.integer(result[[1]][4])    # alito from regex group
      
      # ✅ Save only if threshold is met
      if (!is.na(evalue) && !is.na(alifrom) && !is.na(alito) && evalue <= user_threshold) {
        results <- rbind(results, data.frame(
          Query = current_query,
          Model = current_model,
          Domain = current_domain,
          evalue = evalue,
          alifrom = alifrom,
          alito = alito,
          stringsAsFactors = FALSE
        ))
        #cat("✅ SAVED:", current_query, current_model, current_domain, evalue, alifrom, alito, "\n")
      }
    }
  }
  
  # ✅ Remove overlapping domains if requested
  if (clean && nrow(results) > 0) {
    results <- results[!duplicated(results[c("Query", "Model", "alifrom", "alito")]), ]
  }
  
  # ✅ Trim whitespace from the Query column 
  results$Query <- sapply(strsplit(results$Query, split = " "), function(x) x[1])
  
  # ✅ Add sequence data from FASTA file
  results <- add_sequence_column_protein(results, fasta_file)
  
  # ✅ Save to file
  if (outputfile && nrow(results) > 0) {
    output_file <- sub("\\.out$", "_parsed.csv", input_file)
    write.csv(results, output_file, row.names = FALSE, quote = FALSE)
    cat("✅ Results saved to:", output_file, "\n")
  } else {
    cat("❌ No results to save — check parsing!\n")
  }
  
  return(results)
}