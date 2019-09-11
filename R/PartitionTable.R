#' Utility for Creating a Partition Table from Fasta Alignments
#'
#' Parse through a directory of FASTA files and writes "table.partition" to working directory
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny fasta ortholog busco
#' @param aligned_dir Directory where fasta alignments can be found
#' @export
#' @return Reads through the full table found in folders ./busco_results/run_Busco* and
#'         returns a dataframe of the results
#' @examples
#' PartitionTable(aligned_dir, species_list, threshold)


PartitionTable <- function(aligned_dir){
    old_max <- 1
    cat("#nexus\nbegin sets;\n", file = "table.partition", append = TRUE)
    for (i in 1:length(row.names(missing_df))){
        if (max(missing_df[i,], na.rm = TRUE) > 0) {
            new_max <- old_max + max(missing_df[i,], na.rm = TRUE) #ignores any lines were all are NA
            new_start <- new_max - 1
            cat("charset ", row.names(missing_df)[i], " = ", old_max, "-", new_start, ";\n", sep = "", file = "table.partition", append = TRUE)
            old_max <- new_max
        }
    }
    cat("end;\n", file = "table.partition", append = TRUE)
}
