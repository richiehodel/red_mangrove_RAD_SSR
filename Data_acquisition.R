####################################################################################
### Script by Richie Hodel #########################################################
####################################################################################

###  This script loads in a large SNP data set -- over 25,000 loci for 96 individuals  ###
###  This script also trims the large data set by applying cutoffs for amount of  ###
###  missing data allowed  ###

setwd("/Users/richiehodel/Documents/mangroves/RAD-Seq_vs_msat/red_mangrove_RAD_SSR/")

### import raw data in generic .gdv format
### these data were produced by STACKS in genepop format
### the data were converted to .gdv for convenience

raw <- read.table(file="batch_17d_March2017.gdv", header=TRUE, skip=14, sep="", na.strings="00")

library("stringr")

working <- raw

working$Population <- NULL
working$Individual <- NULL

working_loci <- colnames(working)
working_loci <- data.frame(str_split_fixed(working_loci, "_", 2))

list_loci <- as.vector(working_loci$X1)
colnames(working) <- list_loci


### read in large fasta file (~650MB), then reduce it so that there is 
### only one sequence per locus

fasta <- read.table(file="batch_17_March2017.fa", sep=" ", fill=TRUE)
fasta$V2 <- NULL

odds <- seq(1,nrow(fasta),2)
evens <- seq(2,nrow(fasta),2)

fasta_new <- data.frame(fasta[odds,])
fasta_new <- data.frame(str_split_fixed(fasta_new$fasta.odds..., "_Sample_", 2))
fasta_new$X2 <- NULL
fasta_new$X2 <- fasta[evens,]

fasta_newer <- subset(fasta_new, !duplicated(fasta_new$X1))

### print out fasta file

write.table(fasta_newer, file="all_loci.fasta", quote=FALSE, sep="\n", row.names=FALSE,
            col.names=FALSE)

###  this file is used to search against databases for contamination  ###

### remove Homo sapiens hits

human <- read.table(file="human.tabular")
human_hits <- data.frame(str_split_fixed(human$V1, "_", 2))
human_hits$X1 <- "X"
human_hits$combined <- paste(human_hits$X1, human_hits$X2, sep="")
human_loci <- unique(human_hits$combined)

human_removed <- working[, !colnames(working) %in% human_loci]

### remove kraken microbial hits

kraken <- read.table(file="classified.txt", sep="_", fill=TRUE)

deletions <- seq(2,nrow(kraken),2)
kraken[deletions,] <- 'NA'
kraken <- na.omit(kraken)
kraken$V1 <- c("X")
kraken$combined <- paste(kraken$V1, kraken$V2, sep="")

microbial <- unique(kraken$combined)

microbial_human_removed <- human_removed[, !colnames(human_removed) %in% microbial]

###  now remove hits of fungal/microbial origin

contaminants <- read.table(file="contaminants.tabular", sep="\t")

contaminant_hits <- data.frame(str_split_fixed(contaminants$V1, "_", 2))
contaminant_hits$X1 <- "X"
contaminant_hits$combined <- paste(contaminant_hits$X1, contaminant_hits$X2, sep="")
contaminant_loci <- unique(contaminant_hits$combined)

contaminant_microbial_human_removed <- microbial_human_removed[, 
            !colnames(microbial_human_removed) %in% contaminant_loci]

working <- contaminant_microbial_human_removed



###  identify how many individuals are missing data for each locus  ###
na_count <-sapply(working, function(y) sum(length(which(is.na(y)))))

missing <- data.frame(na_count)

###  set how many individuals can have no data for a given locus  ###
###  this cutoff = 30, meaning 30 (of 96) individuals must have  
###  allele data to keep a locus  ###
na_cutoff <- 66
cutoff <- 96-na_cutoff

###  select loci to keep that have data presence above the threshold  ###
keep <- subset(missing, missing$na_count<=na_cutoff)
keepers <- row.names(keep)


###  creating the new data sets with certain cutoffs  ###
### no cutoff
RAD_25198 <- working                 ### cutoff=0,  loci

### cutoffs shown directly above (cutoff = 30)  ###

RAD_6255 <- working[keepers]

###
na_cutoff <- 46
cutoff <- 96-na_cutoff
keep <- subset(missing, missing$na_count<=na_cutoff)
keepers <- row.names(keep)

RAD_3831 <- working[keepers]
###

###
na_cutoff <- 31
cutoff <- 96-na_cutoff
keep <- subset(missing, missing$na_count<=na_cutoff)
keepers <- row.names(keep)

RAD_2317 <- working[keepers]
###

###
na_cutoff <- 21
cutoff <- 96-na_cutoff
keep <- subset(missing, missing$na_count<=na_cutoff)
keepers <- row.names(keep)

RAD_1180 <- working[keepers]
###

###
na_cutoff <- 13
cutoff <- 96-na_cutoff
keep <- subset(missing, missing$na_count<=na_cutoff)
keepers <- row.names(keep)

RAD_239 <- working[keepers]
###



#######################################################################
#######################################################################
###    Calculating missing data
#######################################################################
#######################################################################

size_239 <- ncol(RAD_239)*nrow(RAD_239)
na_239 <- table(is.na(RAD_239)) 
na_239 <- na_239[2]
percent_239 <- 100*(na_239/size_239)
percent_239

size_1180 <- ncol(RAD_1180)*nrow(RAD_1180)
na_1180 <- table(is.na(RAD_1180)) 
na_1180 <- na_1180[2]
percent_1180 <- 100*(na_1180/size_1180)
percent_1180

size_2317 <- ncol(RAD_2317)*nrow(RAD_2317)
na_2317 <- table(is.na(RAD_2317)) 
na_2317 <- na_2317[2]
percent_2317 <- 100*(na_2317/size_2317)
percent_2317

size_3831 <- ncol(RAD_3831)*nrow(RAD_3831)
na_3831 <- table(is.na(RAD_3831)) 
na_3831 <- na_3831[2]
percent_3831 <- 100*(na_3831/size_3831)
percent_3831

size_6255 <- ncol(RAD_6255)*nrow(RAD_6255)
na_6255 <- table(is.na(RAD_6255)) 
na_6255 <- na_6255[2]
percent_6255 <- 100*(na_6255/size_6255)
percent_6255

size_25198 <- ncol(RAD_25198)*nrow(RAD_25198)
na_25198 <- table(is.na(RAD_25198)) 
na_25198 <- na_25198[2]
percent_25198 <- 100*(na_25198/size_25198)
percent_25198


