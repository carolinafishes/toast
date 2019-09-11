#' Utility for Assessing Missing Data in a Set of FASTA Files
#'
#' Parse through a directory of FASTA files to assess missing data patterns between taxa and loci
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny fasta ortholog busco
#' @param aligned_dir Directory where fasta alignments can be found
#' @export
#' @return Reads through the full table found in folders ./busco_results/run_Busco* and
#'         returns a dataframe of the results
#' @import seqinr
#' @examples
#' missing_data_df <- MissingDataTable(aligned_dir)

MissingDataTable <- function(aligned_dir) {
    busco_list <- list.files(aligned_dir) #list of busco file names, path directory is required
    busco_ids <- gsub(".fasta", "", busco_list) #list of busco ids
    headers <- NULL

    for (n in 1:length(busco_list)){
        each_line <- readLines(paste0(aligned_dir, "/", busco_list[n]))
        headers <- append(headers, each_line[grep(">", each_line)])
    }
    uni_names <- unique(headers)
    uni_names <- gsub(">", "", uni_names)

    #initialize the dataframe with col.names = species, and rows = busco_ids
    empty_matrix <- matrix(NA, length(busco_ids), length(uni_names))
    final_df <- as.data.frame(empty_matrix, row.names = busco_ids)
    colnames(final_df) <- uni_names

    #read in fasta files and begin populating the dataframe with length of alignments
    for (i in 1:length(busco_list)){
    fasta_file <- read.fasta(paste0(aligned_dir, "/", busco_list[i]), as.string = TRUE)

        for (j in 1:length(names(fasta_file))){
            matched <- match(names(fasta_file)[j], colnames(final_df))
            final_df[busco_ids[i], matched] <- nchar(fasta_file[[j]][[1]])
        }
    }
    return(final_df)
}
