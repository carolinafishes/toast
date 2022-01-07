#' Run BUSCO to Assign Ortholog Ids
#'
#' This function allows you to run BUSCO and harvest orthologs.
#'        Note that this function is limited to linux users as BUSCO will only run on linux
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny busco fasta odb9 ortholog
#' @param fasta_dir directory of fastas; Genus_species.fasta
#' @param path_to_run_busco.py full path where the python script can be found
#' @param toast_dir working directory where busco runs will be placed
#' @param path_to_orthDB directory where orthoDB can be found
#' @param threads number of threads to use, defaults to 1
#' @export
#' @return Runs BUSCO on all fasta sequences in fasta_dir
#' @examples
#' RunBusco(fasta_dir = "path/to/fasta/", toast_dir = "/path/to/toast/dir/",
#'           path_to_run_busco.py = "path/to/run/busco.py/", path_to_orthoDB = "path/to/orthoDB")

#fasta_dir <- paste0(toast_directory, "/fasta/")
RunBusco <- function(fasta_dir, toast_dir, path_to_run_busco.py, path_to_orthoDB, threads = 1){
    reset_wd <- getwd()
    if (length(fasta_dir) == 0) {
        stop("We could not find any fasta files, exiting now\n")
    }

    species_fasta <- dir(fasta_dir)[grep(".fasta", dir(fasta_dir))]
    species_list <- gsub(".fasta", "", species_fasta)

    #check if any species have already been analysed with busco
    busco_dir <- paste0(toast_dir, "/busco_results/")
    busco_folders <- dir(busco_dir)[grep("run_", dir(busco_dir))]
    busco_species_prior <- gsub("run_", "", busco_folders)

    if (length(setdiff(species_list, busco_species_prior)) > 0) {
        differences <- setdiff(species_list, busco_species_prior)

        cat("\nPerforming BUSCO on", length(differences), "fasta file(s)\n\nSTEP 2 - RUNNING BUSCO\n")
        setwd(paste0(toast_dir, "/busco_results/"))
        pb2 <- txtProgressBar(min = 0, max = length(differences), style = 3)
        for (i in 1:length(differences)) {
            busco_fasta <- paste0(differences[i], ".fasta")
            system(paste0("python3 ", path_to_run_busco.py, " -c ", threads, " -i ", toast_dir,
                          "/fasta/", busco_fasta, " -o ", differences[i], " -l ", path_to_orthoDB,
                          " -m tran", " >> ", toast_dir, "busco_run_notes.txt"))
            setTxtProgressBar(pb2, i)
        }
        close(pb2)
    } else { cat("\nIt appears BUSCO has already been run on all given fasta files\n") }
    setwd(reset_wd)
}
