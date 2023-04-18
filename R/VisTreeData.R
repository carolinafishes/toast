#' Summarize Missing Data on a Phylogeny
#'
#' This function allows you summarize missing data by taxon on a supplied phylogenetic tree
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @param tsv a data frame of missing data coverage such as the generated “parsed_busco_results.tsv”
#' @param tree a tree object of class phylo
#' @keywords toast barplot missing data sequence DNA phylogeny
#' @import ggplot2 viridisLite circlepackeR packcircles viridis gridExtra data.tree foreach doParallel streamgraph reshape
#' @export
#' @examples
#' missingDataSummary(tsv)


VisTreeData<-function(tsv, tree, option="both",fsize=0.7,lwd=1,color="#336666",...){
	#summarize missing data
	dataSummary<-missingDataSummary(tsv)
	#get the tree to match the data and vice versa
	matchedSet<-treedata(tree, dataSummary)
	if (option=="both"){
	plotTree.barplot(matchedSet$phy, matchedSet$data, args.plotTree=list(fsize=0.7,lwd=1,color="#336666"), args.barplot=list(
	xlab="Loci",
	col=viridis(2)))}
	if (option=="present"){
	plotTree.barplot(matchedSet$phy, matchedSet$data[,1], args.plotTree=list(fsize=0.7,lwd=1,color="#336666"), args.barplot=list(
	xlab="Loci",
	col=viridis(1)))}
	if (option=="absent"){
	plotTree.barplot(matchedSet$phy, matchedSet$data[,2], args.plotTree=list(fsize=0.7,lwd=1,color="#336666"), args.barplot=list(
	xlab="Loci",
	col="#fff033"))}
}
