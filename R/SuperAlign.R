#' Utility for Creating a RAxML-style Partition Table from Fasta Alignments
#'
#' Parse through a directory of FASTA files to assess
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny fasta ortholog busco
#' @param aligned_dir Directory where fasta alignments can be found
#' @param missing_df Missing data object made by the missing function
#' @param filename string for what you want to call your output
#' @export
#' @return Reads through the full table found in folders ./busco_results/run_Busco* and
#'         returns a dataframe of the results
#' @import seqinr
#' @examples
#' partition_df <- SuperAlign(aligned_dir)

#read in "MissingDataTable.R" results

SuperAlign <- function(aligned_dir, missing_df, filename="superalign.phy"){
    num_species <- length(colnames(missing_df))
    old_max <- 0
        for (i in 1:length(row.names(missing_df))){
        if (max(missing_df[i,], na.rm = TRUE) > 0) {
            new_max <- old_max + max(missing_df[i,], na.rm = TRUE) #ignores any lines were all are NA
            old_max <- new_max
        }
    }

    max_align_length <- old_max
    header <- paste0(num_species, "\t", max_align_length)
    cat(file = filename, header, "\n", append = TRUE, sep = '')

    for (i in 1:length(colnames(missing_df))){
        species_to_find <- colnames(missing_df)[i] #change to i later
        thing_to_append <- paste0(species_to_find, "\t")
        for (j in 1:nrow(missing_df)){
            if (is.na(missing_df[j,i]) == FALSE){ #grab the sequence
                fasta_file <- read.fasta(paste0(aligned_dir, "/", row.names(missing_df[j,]), ".fasta"), as.string = TRUE)
                target_species <- match(species_to_find, names(fasta_file))
                seq_to_append <- fasta_file[[target_species]][[1]]
                thing_to_append <- paste0(thing_to_append, seq_to_append)
            }
            if (is.na(missing_df[j,i]) == TRUE){
                dashes <- max(missing_df[j,], na.rm = TRUE) #add in this number of dashes
                dashes_to_append <- strrep("-", dashes)
                thing_to_append <- paste0(thing_to_append, dashes_to_append)
            }
        }
        cat(file = filename, thing_to_append, "\n", append = TRUE, sep = '')
    }
}
