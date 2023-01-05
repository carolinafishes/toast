#' Scrutinize missing data patterns in your sequence data
#'
#' Assess which taxa have no missing data
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny fasta ortholog busco
#' @param tsv A data frame of missing data coverage such as the generated “parsed_busco_results.tsv”
#' @return Returns the names of taxa which have no missing data
#' @export
#' @examples
#' CompleteCoverage(tsv)

CompleteCoverage <- function(tsv){
   cet_data <- tsv
    data_names <- names(cet_data)
    actual_names <- data_names[2:length(data_names)]
    actual_data <- cet_data[2:length(cet_data[,1]),]
    
    actual_data[is.na(actual_data)] <- 0
    actual_data[actual_data != 0]<-1
		actual_data <- actual_data[, 2:length(actual_data[1,])]
    #Plotting with a threshold
    dimensions <- dim(actual_data)
    loci <- dimensions[1]
	actual_data <-data.frame(apply(actual_data, 2, function(x) as.numeric(as.character(x))))
    new.mm <- actual_data[,colSums(actual_data)==loci]
    if (dim(new.mm)[2]==0)
    {
    	return("No taxa have complete coverage")
    }
    if (dim(new.mm)[2]>0)
	{
    return(new.mm)
    }
}
