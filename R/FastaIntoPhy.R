#' Utility for turning fasta file into phylip file
#'
#'
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast DNA phylogeny fasta phylip conversion utility
#' @param filename defaults to file_name
#' @param fastafile defaults to myfastadata
#' @param externalfile defaults to TRUE
#' @export
#' @return This function uses an aligned fasta file that is either in your working directory 
#'         externalfile=TRUE or that has been read into memory, externalfile=false. 
#'         This then converts it to a phylip format for iqtree using the user specifed filename
#' @examples
#' FastaIntoPhy(fastafile=myfastadata, filename="file_name",externalfile=TRUE)

FastaIntoPhy<-function(fastafile="myfastadata", filename="file_name", externalfile=TRUE){
	if (externalfile==TRUE){
	fasta<-read.fasta(fastafile) }
		else if (externalfile==FALSE){
			fasta<-fastafile
	} 
    taxa<-length(names((fasta)))
    x<-summary(fasta[[1]])
    characters<-x$length
    initialize_blankfile<-paste("touch ",filename)
    system(initialize_blankfile)
    #filecon<- file(filename)
    sink(filename)
    Line1 <- paste(taxa, characters, sep="    " )
    cat(Line1,"\n")
	for (i in 1:taxa){
		sequence<-str_c(as.character(fasta[[i]]))
		cat(names(fasta)[i], sequence,sep="    ")
		cat("\n")
	}
    sink()
   }