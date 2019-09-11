#' Assess Data Coverage Across All Loci and Taxa
#'
#' Compare missing data between loci and taxa at different hierarchical levels
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @param tsv A data frame of missing data coverage such as the generated “parsed_busco_results.tsv”
#' @param taxonomy User made data.frame of taxonomy levels following the format of cetacean_taxonomy
#' @param taxonomy_level Target column name to define level
#' @param threshold Cutoff value for minimum number of sequences allowed per taxon
#' @param type Type of plot to generate, either “bar” or “stream”
#' @keywords toast missing data transcriptome sequence DNA streambar experimental design
#' @return Generate either stacked barplots or stream plots to look at missing data patterns
#'         between loci and taxa at user defined hierarchical levels
#' @import ggplot2 viridisLite circlepackeR packcircles viridis gridExtra data.tree foreach doParallel streamgraph reshape
#' @export
#' @examples
#' MissingStreambar(tsv, taxonomy, taxonomy_level, threshold, type="bar")

MissingStreambar <- function(tsv, taxonomy, taxonomy_level, threshold, type){
#uses similar start to interactive circlepack plot but makes a 0,1 matrix
cet_data <- tsv
data_names <- names(cet_data)
actual_names <- data_names[2:length(data_names)]
actual_data <- cet_data[2:length(data_names)]
#work in threshold
dimensions <- dim(actual_data)
loci <- dimensions[1]
limit <- loci-threshold
#get table of missing for threshold
actual_data[, ] <- lapply(actual_data[, ], as.character)

mm <- matrix(nrow = 1, ncol = length(actual_names))
for (i in 1:length(actual_names)){
    subset <- actual_data[,i]
    subset[is.na(subset)] = 1
    list_of_nas <- subset[subset==1]
    missing <- table(list_of_nas)
    mm[,i] <- as.numeric(missing[1])
}

#match missing data to taxonomy by leaves
actual_data <- t(actual_data)
pool <- names(taxonomy)
ii <- which(pool == taxonomy_level)
level <- as.character(taxonomy[,ii][match(taxonomy$leaves,row.names(actual_data))])

#subsample by threshold
if (threshold > 0){
    get <- which(mm<limit)
    actual_data <- actual_data[get,]
    level <- level[get]
}

else { actual_data <- actual_data
	level <- level
}

#convert to 0 and 1 matrix
actual_data[!is.na(actual_data)] = 0
actual_data[is.na(actual_data)] = 1

#assemble for streamgraph
newdf <- foreach(ind = 1:loci, .combine=rbind) %dopar%
{
    cbind(rep(ind,length(level)),level, actual_data[, ind])
}
data.in <- as.data.frame(newdf)
data.in[,1] <- as.numeric(as.character(data.in[,1]))
data.in[,3] <- as.numeric(as.character(data.in[,3]))
data.in[,2] <- as.character(data.in[,2])
colnames(data.in) <- c("locus_name","level","value")
#data.in[which(data.in$value>0),]->data.in.missing

aggregatedf <- aggregate(.~locus_name+level, data.in,sum)
aggregatedf <- aggregatedf[order(aggregatedf$locus_name),]
if (type=="bar"){
p = ggplot(aggregatedf, aes(x=locus_name, y=value, fill=level))	+ geom_bar(stat='identity')
p	} else {
    aggregatedf%>%
	streamgraph("level", "value", "locus_name",scale = "continuous",offset = "zero", interpolate = "step")%>%
    sg_fill_brewer("PuOr")
    }
}
