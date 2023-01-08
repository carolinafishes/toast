#' Prune a phylip alignment by partial text string matches
#'
#' This function allows you to prune a phylip formatted alignment by a vector of full or partial string matches
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @param filepath The full path to your phylip formatted alignment
#' @param targets a vector of partial or full name matches you want to keep 
#' @param filename a newfile name for your pruned alignment
#' @keywords toast barplot missing data sequence DNA phylogeny
#' @import ggplot2 viridisLite circlepackeR packcircles viridis gridExtra data.tree foreach doParallel streamgraph reshape stringr
#' @export
#' @examples
#' PruneSuperAlign(filepath="yourpath/superalign.phy", targets=c("Danio",  "Canis", "Balistes", "Equu"))

PruneSuperAlign <- function(filepath='superalign.txt', targets, fileName="NewFileName"){
   con <- file(filepath, open = 'r')
    header <- readLines(con, n = 1)
    write(header, file = fileName, append = TRUE)
    seqlength<-str_split(header, "\t")[[1]][2]
    write(header, file = "NewAlign.txt", append = TRUE)
    i1<-0
   while(TRUE) {
    line <- readLines(con, n = 1)
    if(length(line) == 0) {break}
    test<-sum(as.integer(str_detect(line, targets)))
     if(test>0){
        write(line, file = fileName, append = TRUE)
        	i1<-i1+1
    } 
  }
  print(c("Replace your first line with", paste(i1, seqlength, sep=" ")))
}

