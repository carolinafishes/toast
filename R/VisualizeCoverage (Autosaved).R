#' Explore Global Patterns of Missing Data
#'
#' Assess coarse missing data patterns using circlepack plots and barplots
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @param tsv A data frame of missing data coverage such as the generated “parsed_busco_results.tsv”
#' @param threshold Cutoff value for the maximum number of missing loci allowed
#' @keywords toast missing decisiveness sequence DNA phylogenomics
#' @export
#' @import ggplot2 viridisLite circlepackeR packcircles viridis gridExtra data.tree foreach doParallel streamgraph reshape
#' @examples
#' VisualizeCoverage(tsv, threshold)

VisualizeCoverage <- function(tsv, threshold = 300){
    #threshold is the maximum number of missing loci
    #tsv is your data_file
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
    limit <- loci-threshold
    actual_data <-data.frame(apply(actual_data, 2, function(x) as.numeric(as.character(x))))
    new.mm <- actual_data[,colSums(actual_data)>limit]
	
	if (dim(new.mm)[[2]]==0){
		return("None of your samples have enough coverage. Adjust your threshold (the maximum number of missing loci) to allow for more missing loci")
	}

    #Getting ready to plot
    data <- data.frame(group = colnames(new.mm), missing = loci-colSums(new.mm))
    packing <- circleProgressiveLayout(data$value, sizetype='area')
    packing$radius=0.95*packing$radius
    data = cbind(data, packing)
    dat.gg <- circleLayoutVertices(packing, npoints=50)
    number <- length(data$group)

    #Plotting Graph A, code modified from R graph gallery
    p1<-ggplot() +

    #Make the bubbles
    geom_polygon(data = dat.gg, aes(x, y, group = id, fill=as.factor(id)), colour = "black", alpha = 0.3) +
    scale_fill_manual(values = viridis(n = number, option = "B")) +

    #Add text in the center of each bubble + control its size
    geom_text(data = data, aes(x, y, size=log(value), label = paste(group,value))) +
    scale_size_continuous(range = c(2,3.5)) +

    #General theme:
    theme_void() +
    theme(legend.position="none") +
    coord_equal()

    #Preparing Data for plot B
    missing_taxa <- actual_data[,colSums(actual_data)<limit]
    	if (dim(missing_taxa)[[2]]==0){
		print("All of your samples meet or exceed your threshold (the maximum number of missing loci allowed)")
	}
    md<-colSums(missing_taxa)
    md2<-rbind(names(md),md)
    md3<-rbind(rep("present",dim(md2)[2]) , md2)
    rownames(md3) <- c("group","species","loci")
    md3<-t(md3)
    altmd<-as.data.frame(md3)
    altmd$loci<-loci-as.numeric(altmd$loci)
    altmd$group<-"missing"
    pa <- rbind(md3, altmd)
    pa$loci <- as.numeric(as.character(pa$loci))

    #Creating Plot B
    plot2 <- ggplot(data = pa, aes(x = species,y = as.numeric(loci),fill = group)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = rev(viridis(n = 2, option = "D"))) +
    ylab("Number of Loci")

    #Plotting Both Graphics
    grid.arrange(p1, plot2, nrow = 2)
}
