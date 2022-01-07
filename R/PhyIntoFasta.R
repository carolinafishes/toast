#' Utility for turning phylip file into fasta file
#'
#'
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast DNA phylogeny fasta phylip conversion utility
#' @param filename defaults to file_name
#' @param phylipfile defaults to myphylipdata
#' @param externalfile defaults to TRUE
#' @export
#' @return This function uses an aligned phylip file that is either in your working directory 
#'         externalfile=TRUE or that has been read into memory, externalfile=false. 
#'         This then converts it to a fasta format writing to file with user specified filename
#' @examples
#' PhyIntoFasta(fastafile= myphylipdata, filename="file_name",externalfile=TRUE)
PhyIntoFasta<-function(phylipfile="myphylipdata", filename="file_name", externalfile=TRUE){
	if (externalfile==TRUE){
	phylip<-read.table(phylipfile, header=TRUE) }
		else if (externalfile==FALSE){
			phylip <-phylipfile
	} 
    initialize_blankfile<-paste("touch ",filename)
    system(initialize_blankfile)
    #filecon<- file(filename)
    sink(filename)
	for (i in 1:length(phylip[,1])){
		cat(">",phylip[i,1])
		cat("\n")
		cat(phylip[i,2])	
		cat("\n")
	}
    sink()
   }