#' Utility for Assessing Missing Data in a Set of FASTA Files
#'
#' Parse through a directory of FASTA files to assess missing data patterns between taxa and loci
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny fasta ortholog busco
#' @param missing missing data object from TOAST
#' @param tree Tree in newick format
#' @export
#' @return Plots a missmap plot of your data sorted by the tree tips
#' @examples
#' missing_data_df <- phyloMissmap(missing, tree)

PhyloMissmap <- function(missing, tree,...) {
    is_tip <- tree$edge[,2] <= length(tree$tip.label)
	ordered_tips <- tree$edge[is_tip, 2]
	tips<-tree$tip.label[ordered_tips]
	missing2<-missing[,tips]
	#make missmap
	
	missmap(missing2, rank.order=FALSE,...)
}
