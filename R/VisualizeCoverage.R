#' Explore Global Patterns of Missing Data
#'
#' Assess coarse missing data patterns using circlepack plots and barplots
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @param tsv A data frame of missing data coverage such as the generated “parsed_busco_results.tsv”
#' @param threshold Cutoff value for minimum number of sequences allowed per taxon
#' @keywords toast missing decisiveness sequence DNA phylogenomics
#' @export
#' @import ggplot2 viridisLite circlepackeR packcircles viridis gridExtra data.tree foreach doParallel streamgraph reshape
#' @examples
#' VisualizeCoverage(tsv, threshold)

VisualizeCoverage <- function(tsv, threshold = 100){
    #threshold is the minimum number of missing loci
    #tsv is your data_file
    cet_data <- tsv
    data_names <- names(cet_data)
    actual_names <- data_names[2:length(data_names)]
    actual_data <- cet_data[2:length(data_names)]
    mm <- matrix(nrow=1, ncol=length(actual_names))

    for  (i in 1:length(actual_names)){
    	subset <- actual_data[,i]
    	list_of_nas <- is.na(subset)
    	missing <- table(list_of_nas)
    	mm[,i] <- missing[2]
    }
    mm[is.na(mm)] <- 0

    #Plotting with a threshold
    dimensions <- dim(actual_data)
    loci <- dimensions[1]
    limit <- loci-threshold
    colnames(mm) <- actual_names
    mm2 <- as.data.frame(mm)
    new.mm <- mm2[which(mm2<limit)]

    #Getting ready to plot
    data <- data.frame(group = colnames(new.mm), value = as.numeric(new.mm))
    packing <- circleProgressiveLayout(data$value, sizetype='area')
    packing$radius=0.95*packing$radius
    data = cbind(data, packing)
    dat.gg <- circleLayoutVertices(packing, npoints=50)
    number <- nlevels(data$group)

    #Plotting Graph A, code modified from R graph gallery
    p1<-ggplot() +

    #Make the bubbles
    geom_polygon(data = dat.gg, aes(x, y, group = id, fill=as.factor(id)), colour = "black", alpha = 0.3) +
    scale_fill_manual(values = viridis(n = number, option = "D")) +

    #Add text in the center of each bubble + control its size
    geom_text(data = data, aes(x, y, size=log(value), label = paste(group,value))) +
    scale_size_continuous(range = c(2,4)) +

    #General theme:
    theme_void() +
    theme(legend.position="none") +
    coord_equal()

    #Preparing Data for plot B
    missing_taxa <- t(mm2[which(mm2>limit)])
    test <- as.data.frame(cbind(row.names(missing_taxa), as.numeric(missing_taxa)))
    V2 <- as.numeric(loci-missing_taxa)
    test2 <- as.data.frame(cbind(row.names(missing_taxa), V2))
    group1 <- rep("missing",length(test[,1]))
    group2 <- rep("present",length(test[,1]))
    missing_data <- cbind(group1,test)
    present_data <- cbind(group2,test2)
    colnames(missing_data) <- c("group","species","loci")
    colnames(present_data) <- c("group","species","loci")
    pa <- rbind(missing_data, present_data)
    pa[,3] <- as.numeric(as.character(pa[,3]))

    #Creating Plot B
    plot2 <- ggplot(data = pa, aes(x = species,y = as.numeric(loci),fill = group)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = viridis(n = 2, option = "D"))

    #Plotting Both Graphics
    grid.arrange(p1, plot2, nrow = 2)
}
