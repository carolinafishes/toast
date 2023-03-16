#' Determine and Appropriate Missing Data Threshold for Loci
#'
#' Graphically assess the impact of a % coverage missing data threshold for loci on changes in the distribution of missing data
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny
#' @param tsv a data frame of missing data coverage such as the generated “parsed_busco_results.tsv”
#' @param coverage cutoff value for desired % coverage per locus
#' @import ggplot2 viridisLite circlepackeR packcircles viridis gridExtra data.tree foreach doParallel streamgraph reshape
#' @export
#' @return Returns graphic of missing data based on defined threshold
#' @examples
#' VisualizeCoverageThreshold("parsed_busco_results.tsv", coverage=0.5)

VisualizeCoverageThreshold <- function(tsv, coverage=0.5){
    data_names <- names(tsv)
    #if reading in the external saved version
    if ("X"%in%names(missing)){
        actual_names <- data_names[2:length(data_names)]
        actual_data <- tsv[2:length(tsv[,1]),]
    } else {
        actual_names <-data_names
        actual_data<-tsv
    }

    actual_data[is.na(actual_data)] <- 0
    actual_data[actual_data != 0]<-1
	actual_data <- actual_data[, 2:length(actual_data[1,])]
	actual_data <-data.frame(apply(actual_data, 2, function(x) as.numeric(as.character(x))))

    ###Getting post-Threshold Plot
    dimensions <- dim(actual_data)
    taxa <- dimensions[2]
    loci <- dimensions[1]

    ###Getting pre-Threshold data
    mm2 <- as.data.frame(actual_data)

    ###Getting post-Threshold data
    filtered_rows <- actual_data[rowSums(actual_data)/taxa<coverage,]
    kept_rows <- actual_data[rowSums(actual_data)/taxa> coverage,]

    ###Getting ready to plot
    data.pre.filter <- data.frame(group = colnames(mm2), value = (loci-colSums(mm2))/loci, treatment=rep("pre", length(colnames(mm2))))
    data.post.filter <- data.frame(group = colnames(kept_rows), value = (loci-colSums(kept_rows)-dim(filtered_rows)[1])/(loci-dim(filtered_rows)[1]), treatment=rep("post", length(colnames(kept_rows))))
    data.filter <- data.frame(group = colnames(filtered_rows), value = (-1*dim(filtered_rows)[1])/loci, treatment=rep("filt", length(colnames(filtered_rows))))

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
      ylab("Percent Missing Data")
}
