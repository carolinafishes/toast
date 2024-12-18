  extractSequences <- function(input_csv, type = "AA", external=TRUE) {
  # Check if the CSV file exists

if (external==TRUE){
  if (!file.exists(input_csv)) {
    stop(paste("Error: CSV file", input_csv, "not found!"))
  }
  
  # Read the CSV file
  csv_data <- read.csv(input_csv, stringsAsFactors = FALSE)
}
if (external==FALSE){
csv_data<-input_csv
}  
  # Rename columns to match expected names
  names(csv_data) <- tolower(names(csv_data))  # Convert column names to lowercase
  if (!all(c("filename", "accession", "hits") %in% names(csv_data))) {
    stop("Error: CSV file must contain columns: 'Filename', 'Accession', 'hits'")
  }
  
  # Loop through each row of the CSV file
  for (i in seq_len(nrow(csv_data))) {
    # Extract and clean values
    filename <- trimws(csv_data$filename[i])
    accession <- trimws(csv_data$accession[i])
    transcript_id <- trimws(csv_data$hits[i])
    
    # Check if the FASTA file exists
    if (!file.exists(filename)) {
      warning(paste("Warning: FASTA file", filename, "not found, skipping!"))
      next
    }
    
    # Output file name
    output_file <- paste0(accession, "_extracted_sequences.fasta")
    
    # Read the FASTA file
    fasta_data <- if (type == "AA") {
      readAAStringSet(filename)  # Protein sequences
    } else if (type == "DNA") {
      readDNAStringSet(filename)  # Nucleotide sequences
    } else {
      stop("Error: Invalid 'type'. Use 'AA' or 'DNA'.")
    }
    
    # Extract the sequence with matching transcript_id
    matching_seq <- fasta_data[grepl(transcript_id, names(fasta_data))]
    
    # Check if matching sequence exists
    if (length(matching_seq) > 0) {
      # Clean the header to include only the transcript_id
      names(matching_seq) <- transcript_id
      
      # Write to the output file (append mode)
      writeXStringSet(matching_seq, filepath = output_file, append = TRUE)
      
      message(paste("Extracted and cleaned sequence for", transcript_id, 
                    "from", filename, "to", output_file))
    } else {
      warning(paste("No matching sequence found for transcript ID:", transcript_id, "in", filename))
    }
  }
  
  message("Processing complete.")
}