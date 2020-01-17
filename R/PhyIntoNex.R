#' Utility for turning TOAST superalign.txt into Nexus file
#'
#'
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny fasta ortholog busco
#' @param filename defaults to superalign.txt
#' @export
#' @return Reads through a dataframe of missing and extracts sequences from the aligned_dir
#'         and writes to threshold_fasta_folder. Fasta files must be realigned and this can be
#'         accomplished with MafftOrientAlign()
#' @examples
#' PhyIntoNex(filename = "superalign.txt")

PhyIntoNex<-function(filename="superalign.txt"){
    #get the first line
    first_line <- readLines(filename,n=1)
    #get the sites
    sites<-sub(".*\t","", first_line)
    #get the taxa
    taxa<-sub("\t.*","", first_line)
    fileroot<-sub("txt.*","", filename)
    remove_first_line<-paste("sed '1d'", filename,">", sep=" ")
    newfile<-paste(fileroot,"nex",sep="")
    remove_first_lineb<-paste(remove_first_line, "temp", sep=" ")
    system(remove_first_lineb)
    nexus_tax<-paste("NTAX=",taxa, sep="")
    nexus_sites<-paste("NCHAR=",sites, sep="")
    Line3 = paste("DIMENSIONS",nexus_tax,nexus_sites,";", sep=" " )
    cat("#NEXUS","BEGIN DATA;", Line3,"FORMAT DATATYPE = DNA GAP = - MISSING = ?;","MATRIX", file= newfile, sep="\n", append="TRUE")
    makenewfile<-paste("cat", "temp",">>", newfile, sep=" ")
    system(makenewfile)
    system("rm temp")
    cat(";","end;",file= newfile,sep="\n", append=TRUE)
}
