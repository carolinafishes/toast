#' Utility for extracting duplicated BUSCO sequences into individual fasta files
#'
#'
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @author Phillip Souza, \email{psouza1@@uncc.edu}
#' @keywords toast busco fasta extraction utility
#' @param tsvtable the resulting tsv table from a busco analysis
#' @param fasta original fasta sequence file
#' @param ed extracted directory where extracted sequences will be written
#' @param seqID name of sequence, defaults to Genus_species
#' @export
#' @return This function uses the output of a busco analysis and specified fasta 
#'         file to extract duplicated busco sequences and write these into fasta files. 
#'         Fasta files are written into the directory specified by the parameter ed.
#'         Each sequence will be named based on the seqID specified
#' @examples
#' extractBuscoSeqsDup(tsvtable=tsvtable, fasta=fasta,ed="path/to/extracted/, seqID="Genus_species")

extractBuscoSeqsDup<-function(tsvtable, fasta,ed, seqID="Genus_species"){
  #gets position of objects labeled as Duplicated 
  dupseqs<-tsvtable[which(tsvtable$Status == "Duplicated"),] #Get all the duplicated sequences
  #dplyr way to do this
  dupseqsTrimmed<-dupseqs%>%dplyr::group_by(X..Busco.id)%>%dplyr::slice(which.max(Length))
  #Initialize progress bar for the impatient :-)
  progbar <- txtProgressBar(min = 0, max = length(dupseqsTrimmed[,1]), style = 3)
  #target sequences  
  splitnames <- stringr::str_split(dupseqsTrimmed$Sequence,":", n=2,simplify = TRUE)[,1]
  #Extract each sequence name
  for (i in 1:length(splitnames)) {
    if (splitnames[i] != "NA") {
      #trycatch will skip errors and report if there is a problem
      tryCatch({
        
        #Here we are just subsampling the fasta, if the file exists already we can append to it
        seq_index <- grep(paste0(splitnames[i], "\\b"), names(fasta))
        gene_start <- as.numeric(completeseqs$Gene.Start[i])
        gene_end <- as.numeric(completeseqs$Gene.End[i])
        strand     <- completeseqs$Strand[i]
        full_seq <- fasta[[seq_index]]
        
        # Extract subsequence based on strand
        if (strand == "+") {
          subseq_obj <- Biostrings::subseq(full_seq, start = gene_start, end = gene_end)
        } else if (strand == "-") {
          # Ensure start < end when slicing
          subseq_obj <- Biostrings::subseq(full_seq, start = gene_end, end = gene_start)
          subseq_obj <- Biostrings::reverseComplement(subseq_obj)
        } else {
          stop("Unknown strand direction")
        }
        
        seq_to_write <- as.character(subseq_obj)
        
        cat(">", seqID, "\n", seq_to_write, 
            "\n", file = paste0(ed, "/", dupseqsTrimmed[i, 
                                                        1], ".fasta"), append = TRUE, sep = "")}, error=function(e){cat("Warning there is an error here:", splitnames[i])})
      setTxtProgressBar(progbar, i) 
    }
  }
}
