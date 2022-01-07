#' Utility for extracting fragmented BUSCO sequences into individual fasta files
#'
#'
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @author Phillip Souza, \email{psouza1@@uncc.edu}
#' @keywords toast busco fasta extraction utility
#' @param tsvtable the resulting tsv table from a busco analysis
#' @param fasta original fasta sequence file
#' @param ed extracted directory where extracted sequences will be written
#' @param seqID name of sequence, defaults to Genus_species
#' @param  threshold minimum number of base pairs required for a fragmented sequence to be extracted. Defaults to 300
#' @export
#' @return This function uses the output of a busco analysis and specified fasta 
#'         file to extract fragmented busco sequences and write these into fasta files. 
#'         Fasta files are written into the directory specified by the parameter ed.
#'         Each sequence will be named based on the seqID specified
#' @examples
#' extractBuscoSeqsFrag(tsvtable=tsvtable, fasta=fasta,ed="path/to/extracted/", seqID="Genus_species",threshold=300)

extractBuscoSeqsFrag<-function(tsvtable, fasta,ed, seqID="Genus_species", threshold=300){
  #gets position of objects labeled as Fragmented
  fragseqs<-tsvtable[which(tsvtable$Status == "Fragmented" & tsvtable$Length > threshold),]#all fragemented sequences greater than 300 bp
  #Initialize progress bar for the impatient :-)
  progbar <- txtProgressBar(min = 0, max = length(fragseqs[,1]), style = 3)
  #target sequences  
  splitnames <- stringr::str_split(fragseqs$Sequence,":", n=2,simplify = TRUE)[,1]
  #Extract each sequence name
  for (i in 1:length(splitnames)) {
    if (splitnames[i] != "NA") {
      #trycatch will skip errors and report if there is a problem
      tryCatch({
        #Here we are just subsampling the fasta, if the file exists already we can append to it
        seq_to_write <-as.character(fasta[[grep(paste0(splitnames[i],"\\b"), names(fasta))]])
        cat(">", seqID, "\n", seq_to_write, 
            "\n", file = paste0(ed, "/", fragseqs[i, 
                                                  1], ".fasta"), append = TRUE, sep = "")}, error=function(e){cat("Warning there is an error here:", splitnames[i])})
      setTxtProgressBar(progbar, i) 
    }
  }
}


