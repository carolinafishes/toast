#' Explore Global Patterns of Missing Data
#'
#' Returns a dataframe that can be used to generated graphics of missing data patterns
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @param tsv A data frame of missing data coverage such as the generated “parsed_busco_results.tsv”
#' @param threshold Cutoff value for minimum number of sequences allowed per taxon
#' @keywords toast missing decisiveness sequence DNA phylogenomics
#' @export
#' @examples
#' VPBR <- isualizePBR(busco_dir = bd)

VisualizePBR <- function(busco_dir){

    busco_folders <- dir(busco_dir)[grep("run_", dir(busco_dir))] #refresh this here as new runs may have been created at STEP 2
    final_columns <- 'busco_id' #initialize the final column names
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
        targets<-unique(MyData[,1])

        for (i in 1:length(targets)){
            specific_species <- MyData[which(MyData[,1]==targets[i]),]
            target_value <- max(specific_species$Score)
            row_needed <- specific_species[which(specific_species$Score==target_value),]
            MyCleanData <- rbind(MyCleanData,row_needed[1,]) #if scores above were identical, this grabs just the first one
        }
        setTxtProgressBar(pb3, n)
        parsed_busco_df <- cbind(parsed_busco_df,MyCleanData[,2]) #[,2] should grab the status (missing, fragmented, complete, duplicated)
    }

    close(pb3)
    names(parsed_busco_df) <- final_columns #rename the columns based on the order they were parsed

    parsed_busco_df[] <- lapply(parsed_busco_df, as.character)
    #change the values in the dataframe to numeric
#        which(parsed_busco_df == "Fragmented") <- 0.5
#        which(parsed_busco_df == "Complete") <- 1
#        which(parsed_busco_df == 'Duplicated') <- 1
#        which(parsed_busco_df == "NA") <- 0
    return(parsed_busco_df)
}
