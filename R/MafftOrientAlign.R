#' Align BUSCO Orthologs
#'
#' This function allows you align each ortholog between all species
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny phy alignment align orient orientate mafft fasta ortholog busco
#' @param extract_dir Directory containing unaligned/unorientated fasta files
#' @param mafft_dir Directory to write orientated and aligned fasta files
#' @param threads Number of cpu threads to use
#' @export
#' @examples
#' MafftOrientAlign(extract_dir = "path/to/extracted/", mafft_dir = "/path/to/mafft/", threads = 1)


MafftOrientAlign <- function(extract_dir, mafft_dir, threads = 1) {
    busco_id_matched_list <- list.files(path = extract_dir) #what happens if mafft tries to align a single sequence?
    pb5 <- txtProgressBar(min = 0, max = length(busco_id_matched_list), style = 3)
    for (i in 1:length(busco_id_matched_list)) {
        setTxtProgressBar(pb5, i)
        file_to_use <- busco_id_matched_list[i]
        system(paste0("mafft --adjustdirection --thread ", threads, " --quiet ", extract_dir, "/", file_to_use,
                      " > ", mafft_dir, "/", file_to_use))
        x <- readLines(paste0(mafft_dir, "/", file_to_use))
        y <- gsub("_R_", "", x)
        cat(y, file = paste0(mafft_dir, "/", file_to_use), sep = "\n")
    }
    close(pb5)
}
