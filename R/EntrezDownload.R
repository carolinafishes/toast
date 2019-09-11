#' Download Public Sequence Data
#'
#' Utility for downloading sequence data from NCBI for a user specified group of organisms
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny ncbi entrez fasta taxa taxonomy
#' @param txid NCBI taxon ID. See \url{https://ncbi.nlm.nih.gov} and manual
#' @param fasta_dir Place to store fasta sequences. Must end with */fasta/
#' @param minimumSeq Minimum number of sequences a species must have to bother downloading
#'                   defaults to 350
#' @param maximumSeq Maximum number of sequences to download per species.
#'                   Useful for testing scripts and troubleshooting TOAST/BUSCO
#' @export
#' @import rentrez
#' @return Returns graphic of missing data based on defined threshold
#' @examples
#' SequenceDownload(txid = 9721, fasta_dir = "path/to/fasta/", minimumSeq = 350)

EntrezDownload <- function(txid , fasta_dir, minimumSeq = 350, maximumSeq = NULL){

    pull1_again <- function(){
        Sys.sleep(0.5)
        tryCatch(pull1(), error = function(e) {pull1_again()})
    }

    pull1 <- function(){
        write_me <- entrez_fetch(db = "nucleotide", web_history = cookie, rettype = "fasta", retstart = j, query_key = qk, retmax = 250)
    }

    pull2_again <- function(){
        Sys.sleep(0.5)
        tryCatch(pull2(), error = function(e) {pull2_again()})
    }

    pull2 <- function(){
        write_me <- entrez_fetch(db = "nucleotide", web_history = cookie, rettype = "fasta", retstart = j, query_key = qk, retmax = dif)
    }

    cat("STEP 1 - GATHERING DATA FROM NCBI\n")
    txid_int <- paste0("txid", txid, "[Organism:exp]") #txid comes from sourcing config.R
    taxon_search <- entrez_search(db = "taxonomy", term = txid_int, retmax = 31000) #look into retmax
    cat("The search term", txid_int, "returned", taxon_search$count, "species hits
        DOWNLOADING SEQUENCES FOR ANY SPECIES WITH >", minimumSeq, "SEQUENCES\n")

    for (i in 1:length(taxon_search$ids)) {
        building_search_term <- paste("txid", taxon_search$ids[i], sep = "")
        search_term <- paste(building_search_term, " AND biomol mrna[prop]", sep = "")
        nuc_search <- entrez_search(db = "nucleotide", term = search_term, use_history = TRUE) #look into the retmax
        cookie <- nuc_search$web_history
        qk <- cookie$QueryKey

        if (is.null(maximumSeq) == FALSE && (nuc_search$count > minimumSeq)) {
            dif <- maximumSeq %% 250
            taxize_summ <- entrez_summary(db = "taxonomy", id = taxon_search$ids[i])
            genus <- taxize_summ$genus
            species <- taxize_summ$species
            file_name <- paste(genus, species, sep = "_")
            fasta_file_name <- paste(file_name, "fasta", sep = ".")
            cat("\nDownloading", maximumSeq,"sequences for", file_name,
                "species number", i, "of", taxon_search$count, "\n")
            pb1b <- txtProgressBar(min = 0, max = maximumSeq, style = 3)

            for (j in seq(from=1, to=maximumSeq, by=250)) { #retmax 500 downloads faster but times out occaisionally
                j <<- j
                tryCatch(write_me <- entrez_fetch(db = "nucleotide", web_history = cookie, rettype = "fasta", retstart = j, query_key = qk, retmax = 250), error = function(e) {pull1_again()})
                write(write_me, file = paste0(fasta_dir, "/", fasta_file_name), append=TRUE)
                setTxtProgressBar(pb1b, j)
            }
            close(pb1b)

        } else {

            if (nuc_search$count > minimumSeq){ #probably change this number to make sure only things with a large number of sequences are included
                dif <- nuc_search$count %% 250
                taxize_summ <- entrez_summary(db = "taxonomy", id = taxon_search$ids[i])
                genus <- taxize_summ$genus
                species <- taxize_summ$species
                file_name <- paste(genus, species, sep = "_")
                fasta_file_name <- paste(file_name, "fasta", sep = ".")
                cat("\nDownloading", nuc_search$count,"sequences for", file_name,
                    "species number", i, "of", taxon_search$count, "\n")
                pb1b <- txtProgressBar(min = 0, max = nuc_search$count, style = 3)

                for (j in seq(from=1, to=nuc_search$count-250, by=250)) { #retmax 500 downloads faster but times out occaisionally
                    j <<- j
                    tryCatch(write_me <- entrez_fetch(db = "nucleotide", web_history = cookie, rettype = "fasta", retstart = j, query_key = qk, retmax = 250), error = function(e) {pull1_again()})
                    write(write_me, file = paste0(fasta_dir, "/", fasta_file_name), append=TRUE)
                    setTxtProgressBar(pb1b, j)
                }

                if (nuc_search$count %% 250 != 0) {
                    j <<- nuc_search$count - dif + 1
                    tryCatch(write_me <- entrez_fetch(db = "nucleotide", web_history = cookie, rettype = "fasta", retstart = j, query_key = qk, retmax = dif), error = function(e) {pull2_again()})
                    write(write_me, file = paste0(fasta_dir, "/", fasta_file_name), append=TRUE)
                    setTxtProgressBar(pb1b, j)
                }
                close(pb1b)
            }
        }

    }
    rm(pull1, pull1_again, pull2, pull2_again, j)
}
