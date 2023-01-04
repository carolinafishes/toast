#' Determine and Appropriate Missing Data Threshold
#'
#' Graphically assess the impact of a missing data threshold on changes in the distribution of missing data
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny
#' @param tsv a data frame of missing data coverage such as the generated “parsed_busco_results.tsv”
#' @param threshold cutoff value for minimum number of sequences allowed per taxon
#' @import ggplot2 viridisLite circlepackeR packcircles viridis gridExtra data.tree foreach doParallel streamgraph reshape
#' @export
#' @return Returns graphic of missing data based on defined threshold
#' @examples
#' VisualizeThreshold("parsed_busco_results.tsv", 1000)

VisualizeThreshold <- function(tsv, threshold){
 	cet_data <- tsv
    data_names <- names(cet_data)
    actual_names <- data_names[2:length(data_names)]
    actual_data <- cet_data[2:length(cet_data[,1]),]
    
    actual_data[is.na(actual_data)] <- 0
    actual_data[actual_data != 0]<-1
	actual_data <- actual_data[, 2:length(actual_data[1,])]
	actual_data <-data.frame(apply(actual_data, 2, function(x) as.numeric(as.character(x))))

    ###Getting post-Threshold Plot
    dimensions <- dim(actual_data)
    loci <- dimensions[1]
    limit <- loci-threshold


    ###Getting pre-Threshold data
    mm2 <- as.data.frame(actual_data)

    ###Getting post-Threshold data
    filtered_rows <- actual_data[,colSums(actual_data)<limit]
    kept_rows <- actual_data[,colSums(actual_data)>limit]

    ###Getting ready to plot
    data.pre.filter <- data.frame(group = colnames(mm2), value = loci-colSums(mm2), treatment=rep("pre", length(colnames(mm2))))
    data.post.filter <- data.frame(group = colnames(kept_rows), value = loci-colSums(kept_rows), treatment=rep("post", length(colnames(kept_rows))))
    data.filter <- data.frame(group = colnames(filtered_rows), value = -1*(loci-colSums(filtered_rows)), treatment=rep("filt", length(colnames(filtered_rows))))

    stacks <- rbind(data.pre.filter, data.filter, data.post.filter)
    colnames(stacks) <- c("Species", "TotalMissing", "Treatment")
    
    n <- length(unique(stacks[,2]))
    viridis_qualitative_pal7 <- viridis_pal()(n)

    scale_fill_discrete <- function(...) {
       scale_fill_manual(..., values = viridis_qualitative_pal7)
    }
#plot
       ggplot(stacks, aes(x=Species, y= TotalMissing)) +
      geom_bar(aes(fill = Species), position = "dodge", stat="identity")+
      facet_grid(~factor(Treatment, levels=c("pre", "filt", "post"))) +
      theme_classic() +
      theme(axis.text.x=element_blank()) +
      ylab("Loci")
}
