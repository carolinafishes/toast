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
    actual_data <- cet_data[2:length(data_names)]
    mm <- matrix(nrow = 1, ncol = length(actual_names))

    for  (i in 1:length(actual_names)){
    	subset <- actual_data[,i]
		list_of_nas <- is.na(subset)
        missing <- table(list_of_nas)
        mm[,i] <- missing[2]
    }
    mm[is.na(mm)] <- 0
	colnames(mm) <- actual_names
    mm2 <- as.data.frame(mm)
    mm3 <- cbind(colnames(mm2),t(mm2))
    complete <- mm3[which(mm3[,2]==0)]
    return(complete)
}
