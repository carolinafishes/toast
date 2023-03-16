#' Summarize Missing Data Patterns By Taxa
#'
#' This function allows you summarize missing data by taxon
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @param tsv a data frame of missing data coverage such as the generated “parsed_busco_results.tsv”
#' @keywords toast barplot missing data sequence DNA phylogeny
#' @import ggplot2 viridisLite circlepackeR packcircles viridis gridExtra data.tree foreach doParallel streamgraph reshape
#' @export
#' @examples
#' missingDataSummary(tsv)

missingDataSummary <- function(tsv){
    data_names <- names(tsv)
    #if reading in the external saved version
    if ("X"%in%names(tsv)){
    tsv<-tsv[,-1]
    actual_names <- data_names[2:length(data_names)]
    actual_data <- tsv[2:length(tsv[,1]),]
    }
    if (!"X"%in%names(tsv)) {
    actual_names <-data_names
    actual_data<-tsv[2:length(tsv[,1]),]
    }
    actual_data[is.na(actual_data)] <- 0
    actual_data[actual_data != 0]<-1
    #Plotting summarize
    dimensions <- dim(actual_data)
    loci <- dimensions[1]
	actual_data <-data.frame(apply(actual_data, 2, function(x) as.numeric(as.character(x))))
    new.mm <- cbind(colSums(actual_data), loci-colSums(actual_data))
    colnames(new.mm)<-c("present", "absent")
    return(new.mm)
}

