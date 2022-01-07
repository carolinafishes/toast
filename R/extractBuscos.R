#' Utility for extracting BUSCO sequences into individual fasta files from multiple directories.
#'
#'
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @author Phillip Souza, \email{psouza1@@uncc.edu}
#' @keywords toast busco fasta extraction utility
#' @param tsvLocations  the paths that contain the tsv files from a busco analysis.
#' @param fastaLocations the paths to the original fasta sequence files
#' @param ed extracted directory where extracted sequences will be written
#' @param sampleIDs name of sequences
#' @param complete whether or not to include complete sequences
#' @param fragmented whether or not to include fragmented sequences
#' @param duplicated whether or not to include duplicated sequences
#' @param threshold minimum number of base pairs required for a fragmented sequence to be extracted.Not used when fragmented = false
#' @export
#' @return This function uses the output of a busco analysis and specified fasta 
#'         file to extract  busco sequences and write these into fasta files. 
#'         Fasta files are written into the directory specified by the parameter ed.
#'         Each sequence will be named based on the seqID specified
#' @examples
#' extractBuscos(tsvLocations=c("path/to/first/tsvfile","path/to/second/tsvfile",...), fasta=c("path/to/first/fastafile","path/to/second/fastafile",...),ed="path/to/extracted/", sampleIDs=c("Genus_species","Othergenus_otherspecies",...),complete=TRUE, fragmented=TRUE, threshold=300, duplicated=TRUE)


extractBuscos<-function(tsvLocations, fastaLocations, ed, SampleIDs,complete=TRUE, fragmented=TRUE, threshold=300, duplicated=TRUE){
  for (i in 1:length(fastaLocations)) {
    #read in the fasta
    fasta<-readDNAStringSet(paste0(fastaLocations[i]))
    #read in the table, note that we need to use skip=2 to remove the two info lines and use the read.delim function
    tsvtable <- read.delim(file= tsvLocations[i], sep = '\t', header = TRUE, skip=2)
    marker<-NULL
    if (complete==TRUE)	{
      extractBuscoSeqsComp(tsvtable, fasta, ed, seqID=SampleIDs[i])
    }
    if (fragmented==TRUE)	{
      extractBuscoSeqsFrag(tsvtable, fasta, ed, seqID=SampleIDs[i],threshold=threshold)
    }
    if (duplicated==TRUE)	{
      extractBuscoSeqsDup(tsvtable, fasta, ed, seqID=SampleIDs[i])
    }
    cat("\n","\n","Extracted BUSCO orthologs from", SampleIDs[i],"\n","\n")  #I'm adding this it will be nice to wrap this line into an across species function next so we can track progress
  }
}
