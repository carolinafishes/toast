#' Assign Go_Slim IDs to BUSCO IDs
#'
#' this function takes the Gene Ontology information contained in the orthoDB and generalizes it
#'         into as as few GO terms as possible given three perspectives: Molecular function,
#'         Biological Process and Cellular Component, each returned as a dataframe
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny fasta busco ortholog gene ontology go slim
#' @param orthogroup_info assing the uncompressed ontology file found in orthoDB/info/*orthogroup_info.txt.gz
#' @param obo download desired ontology information (.obo file) from http://geneontology.org/docs/download-ontology/#go_obo_and_owl
#'            these are consistently updated so make sure to grab the newest
#'            we recommend the GO slim AGR subset (goslim_agr.obo) which can be obtained using the command
#'            wget http://current.geneontology.org/ontology/subsets/goslim_agr.obo
#' @param perspective options are BP(biological process), MF(Molecular function), or CC(cellular component)
#' @import GOstats GSEABase BiocManager
#' @export
#' @examples
#' BP <- ToastGoSlim(orthogroup_info = "path/to/orthoDB/info/*orthogroup_info.txt.gz", obo = "path/to/goslim_agr.obo", perspective = "BP")

ToastGoSlim <- function(orthogroup_info, obo, perspective = "BP"){
    orthogroup_info <- read.delim(orthogroup_info, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
    pbcounter <- 0 #will help keep track of progression in progress bar through loops
    pb7 <- txtProgressBar(min = 0, max = nrow(orthogroup_info), style = 3)
    for (i in 1:nrow(orthogroup_info)){
        go_ids <- unlist(strsplit(as.character(orthogroup_info[i,"BiologicalProcesses"]), ";"))#gene ontology IDs imported from orthoDB info file
        if (length(go_ids) == 1 && is.na(go_ids) == TRUE){ #NA is still length 1; throws a bunch of warnings if it has to evaluate chr of length > 1
            if (exists("empty") == FALSE){ #start the variable if it doesn't exist
                empty <- unlist(as.character(orthogroup_info[i, "OrthoGroupID"]))
            } else { #make a list of empty orthogroups to append to the end
                empty <- append(empty, unlist(as.character(orthogroup_info[i, "OrthoGroupID"])))
            }
            pbcounter <- pbcounter + 1 #keep track of progress through the three perspectives loops
            next()
            } else {
            myCollection <- GOCollection(go_ids)
            slim <- getOBOCollection(obo)
            slimmed <- goSlim(myCollection, slim, perspective)
            row.names(slimmed) <- slimmed$Term
            if (exists("headers") == FALSE) { #start the variable if it doesn't exist
                headers <- unlist(as.character(orthogroup_info[i, "OrthoGroupID"]))
            } else {
            headers <- append(headers, unlist(as.character(orthogroup_info[i, "OrthoGroupID"])))
            }
            if (exists("appended") ==  FALSE) { #start the variable if it doesn't exist
                appended <- slimmed[,1]
            } else {
               appended <- cbind(appended, slimmed[,1])
            }
            pbcounter <- pbcounter + 1 #keep track of progress through the three perspectives loops
        }
        setTxtProgressBar(pb7, pbcounter)
    }
    colnames(appended) <- headers
    #need to merge back in the orthogroup_info[i,"BiologicalProcesses"] that were empty
    zeros <- nrow(appended) * length(empty) #will help populate the emptry matrix
    empty_matrix <- matrix(rep(0, zeros), nrow = nrow(appended)) #create an matrix full of zeros with
    colnames(empty_matrix) <- empty                              #colnames to merge with appended
    merger <- cbind(appended, empty_matrix) #this is the desired matrix
    row.names(merger) <- slimmed$Term
    return(merger)
    close(pb7)
}

