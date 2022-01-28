#' Utility for filtering genetrees based on ingroup monophyly or coverage
#'
#'
#' @author Derrek Pope, \email{dpope8@@uncc.edu}
#' @keywords toast phylogeny genetree filtration missing data utility
#' @param ingroup input vector of ingroups
#' @param outgroup input vector of vector of outgroups
#' @param threshold the minimum number of ingroup sequences, defaults to 1
#' @param tree_file_dir string of the directory path for tree files
#' @param tree_file_ext string extension of tree files defaults to ".contree" 
#' @export
#' @return This function uses a user supplied vector of a target ingroup and outgroup taxa along with a path 
#'         to a directory containing newick format genetrees and their extension (e.g. ".tre", ".contree", etc) to assess sufficient taxon coverage 
#'         within a specific ingroup using a user defined threshold, and whether the ingroup of interest is resolved as monophyletic.
#'         The function returns a dataframe that gives each genetree name, the number of ingroup taxa present, and whether
#'         these pass inspection. Trees that do not pass due to low coverage, lack of outgroups, or non-monophyly are indicated by the term "omit".
#' @examples
#' ingroupTaxa<-c("TaxonA", "TaxonB", "TaxonC", "TaxonD", "TaxonE")
#' outgroupTaxa<-c("TaxonZ", "TaxonY", "TaxonX")
#' geneTreeFilter(ingroup=ingroupTaxa, outgroup=outgroupTaxa, threshold=1, tree_file_dir="path/to/files", tree_file_ext=".tre")

geneTreeFilter <- function(ingroup, outgroup, threshold=1, tree_file_dir, tree_file_ext){ 
  trees<-list.files(tree_file_dir, pattern=tree_file_ext)
  tree_info<-matrix(nrow=length(trees), ncol=3)
  for (i in 1:length(trees)){
    file<-trees[i]
    tree_info[i,1]<-file
    tree<-read.tree(file)
    taxa<-tree$tip.label
    coverage<-length(which(taxa %in% ingroups)==TRUE)
    out_coverage<-length(which(taxa %in% outgroups)==TRUE)
    if (coverage<threshold){
      tree_info[i,2]<-coverage
      tree_info[i,3]<-"Omit"
    }
    if (out_coverage==0){
      tree_info[i,2]<-"No_outs"
      tree_info[i,3]<-"Omit"
    } 
    if (coverage=>threshold) {
      outs<-taxa[taxa %in% outgroups]
      rooted<-root(tree, outs[1])
      ins<-outgroups[outgroups %in% taxa]	
      verdict<-is.monophyletic(rooted, ins)	
      if (verdict == TRUE){
        tree_info[i,2]<-coverage
        tree_info[i,3]<-"Pass"
      }
      if (verdict == FALSE){
        tree_info[i,2]<-coverage
        tree_info[i,3]<-"Omit"
      }
    }
  }
  return(tree_info)
}