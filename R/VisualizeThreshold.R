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
    actual_data <- cet_data[2:length(data_names)]
    mm <- matrix(nrow = 1, ncol = length(actual_names))

    for (i in 1:length(actual_names)){
        subset <- actual_data[,i]
        list_of_nas <- is.na(subset)
        missing <- table(list_of_nas)
        mm[,i] <- missing[2]
    }
    mm[is.na(mm)] <- 0

    ###Getting post-Threshold Plot
    dimensions <- dim(actual_data)
    loci <- dimensions[1]
    limit <- loci-threshold
    colnames(mm) <- actual_names

    ###Getting pre-Threshold data
    mm2 <- as.data.frame(mm)

    ###Getting post-Threshold data
    filtered_rows <- which(mm2>limit)

    ###Getting ready to plot
    data.pre.filter = data.frame(group = colnames(mm2), value = as.numeric(mm2))
    data.post.filter <- data.pre.filter
    data.post.filter[filtered_rows,2] <- 0
    data.filter <- data.pre.filter
    data.filter$value <- -1*(data.pre.filter$value-data.post.filter$value)

    stacks <- cbind(data.pre.filter, data.filter$value, data.post.filter$value)
    colnames(stacks) <- c("species", "Total Missing", "Removed", "Remaining")
    stacks2 <- t(stacks)
    stacks2 <- as.data.frame(stacks2)
    identity <- c("species","Total Missing", "Removed", "Remaining")
    stacks2.1 <- cbind(identity, stacks2)
    params <- dim(stacks2.1)
    edge <- params[1]
    edge2 <- params[2]
    st2 <- stacks2.1[2:edge,]
    species <- stacks[,1]
    species <- c("identity", as.character(species))
    colnames(st2) <- species
    stacks.m <- melt(st2, id.vars='identity')
    stacks.m$value <- as.numeric(as.character(stacks.m$value))
    n <- length(unique(stacks.m[,2]))
    viridis_qualitative_pal7 <- viridis_pal()(n)

    scale_fill_discrete <- function(...) {
       scale_fill_manual(..., values = viridis_qualitative_pal7)
    }
    ggplot(stacks.m, aes(x=identity, y=value)) +
      geom_bar(aes(fill = variable), position = "dodge", stat="identity")+
      theme_bw()
}
