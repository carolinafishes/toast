rm(list=ls()) ##clear out any previous data from other R sessions

install.packages("devtools")
library("devtools")
devtools::install_github("carolinafishes/toast")
library("toast")

#set up some variables now to make things easier down the road
#these directories must exist - TOAST will not make them for you
td <- "/home/dustin/temp/toast/trial1" #toast_directory
fd <- "/home/dustin/temp/toast/trial1/fasta" #fasta_dir
bs <- "/home/dustin/software/busco/scripts/run_BUSCO.py" #path to busco_script
bd <- "/home/dustin/temp/toast/trial1/busco_results" #path to busco results directory
ed <- "/home/dustin/temp/toast/trial1/extracted" #extracted_dir
md <- "/home/dustin/temp/toast/trial1/mafft_aligned" #mafft_dir
od <- "/home/dustin/temp/toast/trial1/350_laurasiatheria_odb9" #path to orthoDB directory
ad <- "/home/dustin/temp/toast/trial1/mafft_aligned" #mafft_dir, which is a directory of aligned fastas
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
missing_df <- MissingDataTable(aligned_dir = ad)

#have a look at how much data you are missing

colSums(is.na(missing_df))

#remove species/samples that are missing too much data

threshold_df <- ThresholdDataTable(missing_df = missing_df, threshold = 100)

#threshold_df can now serve as "parsed_busco_table" from previous steps to generate alignment

ThresholdExtract(aligned_dir = ad, missing_df = threshold_df, threshold_fasta_folder = "home/dustin/temp/toast/trial1/threshold100")

#table.partition

PartitionTable(aligned_dir = md) #writes table.partition to working directory

#multi-gene alignment

SuperAlign(aligned_dir)
