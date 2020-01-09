#' Assemble Single Ortholog IDs from BUSCO
#'
#' Reads through the full table found in folders ./busco_results/run_Busco* and
#'         returns a .tsv of the results
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny fasta busco ortholog
#' @param busco_dir Directory where busco results were written. FIles must be in the format
#'                  run_Busco_Genus_species
#' @param ortho_dir = Directory where orthoDB resides
#' @param  length_threshold length_threshold should be a value between 0 and 1.
#'                          A value of 0 will not set a required length for each busco id.
#'                          A value of 1 will not include sequences that are less than the
#'                          length of the busco id.
#' @export
#' @examples
#' parsed_busco_results <- ParseBuscoResults(busco_dir = "path/to/busco/results" ortho_dir = "path_to_od", length_threshold = 0)

ParseBuscoResults <- function(busco_dir, ortho_dir, length_threshold = 0){

    if(length_threshold > 1) {
        stop("length_threshold should be a value between 0 and 1. A value of 0 will not set a required length for each busco id. A value of 1 will not include sequences that are less than the length of the busco id.")
    }

    busco_folders <- dir(busco_dir)[grep("run_", dir(busco_dir))] #refresh this here as new runs may have been created at STEP 2
    final_columns <- 'busco_id' #initialize the final column names

    busco_lengths <- read.table(paste0(ortho_dir, "/lengths_cutoff")) #lengths of each busco_id

    pb3 <- txtProgressBar(min = 0, max = length(busco_folders), style = 3)
    for (n in 1:length(busco_folders)){

        sub_busco <- gsub("run_", "", busco_folders[n])
        table_to_read <- paste(busco_dir, "/", busco_folders[n], "/full_table_", sub_busco, ".tsv", sep = "") #need to get species list in this function
        MyData <- read.delim(table_to_read, skip = 4, header = TRUE, sep = "\t", fill = TRUE, na.strings = "NA")

        if (exists("row_names") == FALSE){ #if it doesn't exist, create it
            row_names <- unique(MyData[,1]) #get the row names only the first pass through
        }

        if (exists("parsed_busco_df") == FALSE){ #if it doesn't exist, create it
            MyInitial <- data.frame(row_names, NA) #add the row_names to the final data
            matrix_holder <- MyInitial[,1]
            parsed_busco_df <- data.frame(matrix_holder)
        }

        final_columns <- append(final_columns, sub_busco)
        MyCleanData <- data.frame(NULL)
        targets <- unique(MyData[,1])

        for (i in 1:length(targets)){
            specific_species <- MyData[which(MyData[,1]==targets[i]),] #Shows sequences for each busco hit
            target_value <- max(specific_species$Score) #Select sequence with the better score
            row_needed <- specific_species[which(specific_species$Score==target_value),]

            if(!is.na(row_needed[1,1])) {
                current_length_threshold <- busco_lengths[which(busco_lengths[,1] == targets[i]), 4] * length_threshold #Create minimum length sequence should be to include

                if(row_needed$Length[1] < current_length_threshold) {
                    row_needed[1,] <- NA
                }
            }

            MyCleanData <- rbind(MyCleanData,row_needed[1,]) #if scores above were identical, this grabs just the first one

        }
        setTxtProgressBar(pb3, n)
        parsed_busco_df <- cbind(parsed_busco_df, MyCleanData[,3]) #tosses out warnings when all are missing but this is fine to ignore
    }

    close(pb3)
    names(parsed_busco_df) <- final_columns #rename the columns based on the order they were parsed
    return(parsed_busco_df)
}

