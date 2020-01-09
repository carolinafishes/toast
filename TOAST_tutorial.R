rm(list=ls()) ##clear out any previous data from other R sessions

install.packages("devtools")
library("devtools")
devtools::install_github("carolinafishes/toast")
library("toast")

#set up some variables now to make things easier down the road
#these directories must exist - TOAST will not make them for you
td <- "/home/dustin/projects/toast_revisions/camel_toast" #toast_directory
fd <- "/home/dustin/projects/toast_revisions/camel_toast/fasta" #fasta_dir
bs <- "/home/dustin/software/busco/scripts/run_BUSCO.py" #path to busco_script
bd <- "/home/dustin/projects/toast_revisions/camel_toast/busco_results" #path to busco results directory
ed <- "/home/dustin/projects/toast_revisions/camel_toast/extracted" #extracted_dir
md <- "/home/dustin/projects/toast_revisions/camel_toast/mafft_aligned" #mafft_dir
od <- "/home/dustin/projects/toast_revisions/camel_toast/laurasiatheria_odb9" #path to orthoDB directory
ad <- "/home/dustin/projects/toast_revisions/camel_toast/mafft_aligned" #mafft_dir, which is a directory of aligned fastas
go_oro <- "/home/dustin/projects/go_slim/goslim_agr.obo" #this must be downloaded online
ortho <- "/home/dustin/projects/toast_revisions/camel_toast/laurasiatheria_odb9/info/laurasiatheria_314145_OrthoDB9_orthogroup_info.txt" #make sure you unzip this file
cpu <- 12 #number of threads to use at various steps

setwd(td)

#step 1 download the sequences
EntrezDownload(txid = 9721, fasta_dir = fd, minimumSeq = 350, maximumSeq = NULL)

#step 2 run BUSCO
RunBusco(fasta_dir = fd, toast_dir = td, path_to_run_busco.py = bs, path_to_orthoDB = od, threads = cpu)

#step 3 parse full busco tables
parsed_busco_results <- ParseBuscoResults(busco_dir = bd)

write.table(parsed_busco_results, file = paste0(td, "/parsed_busco_results.tsv"), sep = "\t", row.names = FALSE) #write it to a tsv inside the working directory

#step 4 extract busco sequences
ExtractBuscoSeqs(busco_table = parsed_busco_results, fasta_dir = fd, extract_dir = ed) #parsed_busco_results from previous step

#step 5 align using mafft
MafftOrientAlign(extract_dir = ed, mafft_dir = md, threads = cpu) #important as some sequences may be 3'->5' direction

#missing data check
mdf <- MissingDataTable(aligned_dir = ad)

#have a look at how much data you are missing

colSums(is.na(mdf))

#remove species/samples that are missing too much data

threshold_df <- ThresholdDataTable(missing_df = missing_df, threshold = 100)

#threshold_df can now serve as "parsed_busco_table" from previous steps to generate alignment

ThresholdExtract(aligned_dir = ad, missing_df = threshold_df, threshold_fasta_folder = "home/dustin/temp/toast/trial1/threshold100")

#table.partition

PartitionTable(aligned_dir = md, missing_df = mdf) #writes table.partition to working directory

#multi-gene alignment

SuperAlign(aligned_dir, missing_df = mdf)

#############################
#############################
#### VISUALIZE YOUR DATA ####
#############################
#############################

### have a look at the biological processes assigned to the orthoDB BUSCO Ids via Gene Ontology
BP <- ToastGoSlim(orthogroup_info = ortho, obo = go_oro, perspective = "BP") #converts BUSCO_IDs into a dataframe of Biological Process GoSlim Terms
#MF and CC are currently not assigned in orthoDB and therefore not supported in TOAST
#Visualize the overlapping GoSlim terms using the nifty package UpSetR
# check out their project here: https://github.com/hms-dbmi/UpSetR
if ("UpSetR" %in% rownames(installed.packages()) == FALSE) {
    install.packages("UpSetR")
}

library(UpSetR)
#Create intersection plot
#format BP data set to plot using the upset() function
rownames(BP)[21] <- "carb. derivative metabolic"
to_plot <- t(BP)
to_plot[to_plot > 0] <- 1
to_plot <- data.frame(to_plot)
to_plot <- cbind(rownames(to_plot), to_plot) #set rownames to be the first row, remove rownames later
colnames(to_plot)[1] <- "Identifier"
rownames(to_plot) <- NULL #remove rownames

#format BP data set to plot
to_plot[,2:ncol(to_plot)] <- sapply(to_plot[2:ncol(to_plot)], as.integer)

#create intersection plot
upset(to_plot, order.by = "freq", sets = colnames(to_plot)[2:ncol(to_plot)])

### Create an occumpany matrix plot to display which Busco Ids were found in each sample

VPBR <- VisualizePBR(busco_dir = bd) #creates a data.frame, similar to the function ParseBuscoResults()

dim(VPBR) #check out the dimensions

vpbr <- as.matrix(VPBR[,2:length(colnames(VPBR))]) #convert VPBR data.frame to vpbr matrix
vpbr <- gsub("Fragmented", 0.5, vpbr)              # format the data frame into a matrix and replace Fragmented,
vpbr <- gsub("Complete", 1.0, vpbr)                # Complete, Duplicated and NA (ie missing) with numeric values
vpbr <- gsub("Duplicated", 1.0, vpbr)
vpbr[is.na(vpbr)] <- 0

class(vpbr) <- "numeric" #image function below needs numeric or logical values to plot

vpbr <- vpbr[,c("busco_Cambac_brain.fasta", "busco_Camdro_brain.fasta",
                        "busco_Cambac_hypothalamus.fasta", "busco_Camdro_hypothalamus.fasta",
                        "busco_Cambac_kidney.fasta", "busco_Camdro_kidney.fasta",
                        "busco_Cambac_liver.fasta", "busco_Camdro_liver.fasta",
                        "busco_Cambac_lung.fasta", "busco_Camdro_lung.fasta",
                        "busco_Cambac_muscle.fasta", "busco_Camdro_muscle.fasta",
                        "busco_Cambac_skin.fasta", "busco_Camdro_skin.fasta",
                        "busco_Cambac_testis.fasta", "busco_Camdro_testis.fasta")]

# Create y axis labels
yLabels <- gsub("busco_", "", colnames(vpbr))#remove the leading busco_ and trailing .fasta
yLabels <- gsub(".fasta", "", yLabels)

ColorRamp <- gray.colors(3, start = 1, end = 0, gamma = 2.2, rev = FALSE) #white = absent, gray = fragmented, black = duplicated/complete
Colors <- c("#FFFFFF", "#BABABA", "#57A0D3")
# Set layout
par(mar = c(5,15,2.5,1), font = 2) #mar=C(bot, left, top, righ 45-t)

# Make the plot
image(1:nrow(vpbr), 1:ncol(vpbr), vpbr,
      col=Colors, xlab="", ylab="",
      axes=FALSE, main= NA)

# Annotate the y axis
box()
axis(side = 2, at=seq(1,length(yLabels),1), labels=yLabels, las= 1,
     cex.axis=1)

### Create Stacked Bar Grapsh showing how many orthologs were found belonging to the 21 GoTerms

# create color palette:
library(RColorBrewer)
coul <- brewer.pal(8, "Pastel1")
coul <- c(coul, brewer.pal(8, "Pastel2"))

data_percentage <- as.matrix(sorted_pie_df/rowSums(sorted_pie_df)*100)

class(data_percentage)

# Make a stacked barplot--> it will be in %!
end_point = 0.5 + nrow(sorted_pie_df) + nrow(sorted_pie_df)-1
par(mar = c(15,5,1,15))
barplot(t(data_percentage), col=coul , border="white", space = 1.0, xlab = "", ylab = "Percent",
        main = "Percent Sample within GoTerm Category",
        axisnames = F)
text(seq(1.5,end_point,by=2), par("usr")[3]-0.25,
     srt = 60, adj= 1, xpd = TRUE,
     labels = paste(rownames(sorted_pie_df)), cex=0.9)

legend(x = 45, y = 100, legend = rev(colnames(sorted_pie_df)), fill = rev(coul), xpd = 1)
