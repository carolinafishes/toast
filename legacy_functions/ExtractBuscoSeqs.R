#' Harvest FASTA Sequences of BUSCO Orthologs and Write to File
#'
#' Extracts BUSCO IDs found in parsed_busco_results.tsv and writes them to
#'         fasta files. Fasta files are named after the BUSCO ID and each
#'         sequence is named after the Genus_species
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny
#' @param busco_table Resulting dataframe from ParseBuscoResults
#' @param fasta_dir Directory containing matching fasta files for each busco run
#' @param extract_dir Directory to write fasta results to
#' @seealso \code{\link{?ParseBuscoResults()}}
#' @export
#' @import seqinr
#' @examples
#' ExtractBuscoSeqs(busco_table = "tsv_file", fasta_dir = "path/to/fasta/directory")

ExtractBuscoSeqs <- function(busco_table, fasta_dir, extract_dir){

pb4 <- txtProgressBar(min = 0, max = length(names(busco_table)), style = 3)
for (i in 2:length(busco_table[1,])) { #first column is busco.id so start the count at 2
    fasta_file <- read.fasta(file = paste0(fasta_dir, "/", colnames(busco_table)[i], ".fasta"), as.string = TRUE)
    setTxtProgressBar(pb4, i)
    for (j in 1:length(busco_table[,1])) {
        cell_value <- toString(busco_table[j,i])

        if (cell_value != "NA") {
            seq_to_write <- gsub("(.{60})", "\\1\n", fasta_file[[cell_value]][[1]])
            cat(">", colnames(busco_table)[i], "\n", seq_to_write, "\n", file = paste0(extract_dir, "/", busco_table[j,1], ".fasta"), append = TRUE, sep = "")
        }
    }
}
close(pb4)
}
