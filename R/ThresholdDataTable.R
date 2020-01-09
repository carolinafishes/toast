#' Utility for retrieving species sequences givin a missing data threshold
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
#' @param missing_df Directory where fasta alignments can be found
#' @param threshold Minimum number of sequences a sample/species must contain
#' @export
#' @return Reads through the full table found in folders ./busco_results/run_Busco* and
#'         returns a dataframe of the results
#' @import seqinr
#' @examples
#' threshold_df <- ThresholdDataTable(missing_df = mdf, threshold = 50)

ThresholdDataTable <- function(missing_df, threshold) {

    appended <- NULL
    max <- nrow(missing_df)
    summary <- colSums(is.na(missing_df)) #make a list of how many null variables are found
    for (i in 1:length(summary)){
        if (summary[i] >= (max - threshold)){
            appended <- append(appended, summary[i])
        }
    }

    wanted <- names(appended)

    for (j in 1:length(wanted)){
            if (exists("threshold_df") == FALSE){ #if it doesn't exist, create it
                threshold_df <- missing_df[appended[j]]
            } else {
                threshold_df <- cbind(threshold_df, missing_df[wanted[j]])
            }
    }
    return(threshold_df)
}




