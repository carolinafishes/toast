#' Explore Missing Data Interactively Based on a Hierarchical Structure
#'
#' This function allows you generate circlepack plots to assess patterns of missing data at different hierarchical levels interactively
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @param tsv a data frame of missing data coverage such as the generated “parsed_busco_results.tsv”
#' @param taxonomy user made data.frame of taxonomy levels following the format of cetacean_taxonomy
#' @param threshold cutoff value for minimum number of sequences allowed per taxon
#' @keywords toast barplot missing data sequence DNA phylogeny
#' @import ggplot2 viridisLite circlepackeR packcircles viridis gridExtra data.tree foreach doParallel streamgraph reshape
#' @export
#' @examples
#' VisualizeTaxonomyInteractive(tsv, taxonomy, threshold)

VisualizeTaxonomyInteractive <- function(tsv, taxonomy, threshold){

    #tsv is your data_file
    #taxonomy is a csv of taxonomic groups
    #threshold is the maximum number of missing loci for a species
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
    #get the number of columns for later
    tt <- dim(taxonomy)
    columns <- tt[2]

    #match missing data to taxonomy by leaves
    mm[is.na(mm) ]< -0.00001
    mm_names <- cbind(actual_names,mm[1,])
    colnames(mm_names) <- c("actual_names","missing_data")
    mm_names <- as.data.frame(mm_names)
    taxonomy$missing <- mm_names$missing_data[match(taxonomy$leaves, mm_names$actual_names)]

    #work in threshold
    actual_data <- cet_data[2:length(data_names)]
    dimensions <- dim(actual_data)
    loci <- dimensions[1]
    limit <- loci-threshold

    if (threshold > 0){
        taxonomy <- taxonomy[which(as.numeric(as.character(taxonomy$missing))<limit),]
    } else { taxonomy <- taxonomy}

    #much of this is modified from the circlepackR tutorial
    pathString <- paste("world", taxonomy$level1, sep = "/" )

    for(i in 2:columns){
    	pathString <- paste(pathString, taxonomy[,i], sep = "/" )
    }

    taxonomy$pathString <- pathString
    population <- as.Node(taxonomy)

    #You can custom the minimum and maximum value of the color range
    circlepackeR(population, size = "missing", color_min = "hsl(206,56%,34%)", color_max = "hsl(367,62%,79%)")
}
