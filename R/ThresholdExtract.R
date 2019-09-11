#' Utility for retrieving species sequences givin a dataframe of missing data
#'
#' User defines a threshold and is returned a dataframe containing
#'        only the samples/species with at least that many sequences.
#'        This dataframe can be fed into ExtractBuscoSeqs() as the busco_table
#'        parameter
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny fasta ortholog busco
#' @param aligned_dir Directory where fasta alignments can be found
#' @param missing_df Dataframe created by ParseBuscoResults() or MissingDataTable()
#' @param threshold_fasta_folder Folder to store extracted fasta files
#' @export
#' @return Reads through a dataframe of missing and extracts sequences from the aligned_dir
#'         and writes to threshold_fasta_folder. Fasta files must be realigned and this can be
#'         accomplished with MafftOrientAlign()
#' @import seqinr
#' @examples
#' ThresholdExtract(aligned_dir = ad, missing_df = threshold_df, threshold_fasta_folder = "/path/to/store/fastas/threshold100")

ThresholdExtract <- function(aligned_dir, missing_df, threshold_fasta_folder){
    for(j in 1:nrow(missing_df)){ #parse through each row, by column
        for(i in 1:ncol(missing_df)){
            cell_value <- toString(missing_df[j,i])
            if(cell_value != "NA" ){
                partial_fasta_name <- rownames(missing_df)[j]
                full_fasta_name <- paste0(partial_fasta_name, ".fasta")
                fasta_file <- read.fasta(file = paste0(aligned_dir, "/", full_fasta_name), as.string = TRUE)
                species_to_get <- colnames(missing_df)[i]
                seq_to_write <-  gsub("(.{60})", "\\1\n", fasta_file[[species_to_get]][[1]])
                cat(">", full_fasta_name, "\n", seq_to_write, "\n", file = paste0(threshold_fasta_folder, "/", full_fasta_name), append = TRUE, sep = "")
            }
        }
    }
}
