library(tidyverse)

# Define the function to convert CSV to FASTA format
csv2Fasta <- function(input_csv, output_fasta="mysequences.fasta", external=TRUE,seqID="Genus",seq="sequence") {
#if a file path is given
  if (external==TRUE){
  # Read the CSV file
  data <- read.csv(input_csv, header = TRUE, stringsAsFactors = FALSE)
  }
# if a dataframe already in R memory
  if (external==FALSE){
  	data<-input_csv
  }
  
  # Create the FASTA format
  fasta_data <- data %>%
    select(seqID, seq) %>%
    mutate(Fasta_Format = paste0(">", seqID, "\n", seq))
  
  # Write the FASTA format to the output file
  writeLines(fasta_data$Fasta_Format, output_fasta)
}