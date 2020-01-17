#' Utility for returning TOAST missing values as a dataframe
#'
#'
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny fasta ortholog busco
#' @param filename defaults to superalign.txt
#' @export
#' @return converts toast mdf into a dataframe containing only the missing values
#' @examples
#' MissingValue(missing_data = mdf)

MissingValue <- function(missing_data) {
    data_names <- names(missing_data)
    actual_names <- data_names[2:length(data_names)]
    actual_data <- missing_data[2:length(data_names)]
    mm <- matrix(nrow = 1, ncol = length(actual_names))
    for (i in 1:length(actual_names)) {
        subset <- actual_data[, i]
        list_of_nas <- is.na(subset)
        missing <- table(list_of_nas)
        mm[, i] <- missing[2]
    }
    mm[is.na(mm)] <- 0
    dimensions <- dim(actual_data)
    loci <- dimensions[1]
    limit <- loci - threshold
    colnames(mm) <- actual_names
    mm2 <- as.data.frame(mm)
    return(mm2)
}
