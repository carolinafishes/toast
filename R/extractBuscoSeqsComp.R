#' Utility for extracting complete BUSCO sequences into individual fasta files
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
#'         file to extract complete busco sequences and write these into fasta files. 
#'         Fasta files are written into the directory specified by the parameter ed.
#'         Each sequence will be named based on the seqID specified
#' @examples
#' extractBuscoSeqsComp(tsvtable=tsvtable, fasta=fasta,ed="path/to/extracted/, seqID="Genus_species")

extractBuscoSeqsComp<-function(tsvtable, fasta,ed, seqID="Genus_species"){
  #gets position of objects labeled as complete, note I changed this to now reflect the header in 1 line
  completeseqs<-tsvtable[which(tsvtable$Status == "Complete"),]#only complete sequences
  #Initialize progress bar for the impatient :-)
  progbar <- txtProgressBar(min = 0, max = length(completeseqs[,1]), style = 3)
  #target sequences  
  splitnames <- str_split(completeseqs$Sequence,":", n=2,simplify = TRUE)[,1]
  #Extract each sequence name
  for (i in 1:length(splitnames)) {
    if (splitnames[i] != "NA") {
      #trycatch will skip errors and report if there is a problem
      tryCatch({
        #Here we are just subsampling the fasta, if the file exists already we can append to it
        seq_to_write <-as.character(fasta[[grep(paste0(splitnames[i],"\\b"), names(fasta))]])
        cat(">", seqID, "\n", seq_to_write, 
            "\n", file = paste0(ed, "/", completeseqs[i, 
                                                      1], ".fasta"), append = TRUE, sep = "")}, error=function(e){cat("Warning there is an error here:", splitnames[i])})
      setTxtProgressBar(progbar, i) 
    }
  }
}