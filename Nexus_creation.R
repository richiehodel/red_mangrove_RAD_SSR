####################################################################################
### Script by Richie Hodel #########################################################
####################################################################################


## This script reads in a haplotype file that was outputted by the STACKS
## program "populations" and converts it to a nexus file that is used in 
## subsequent analyses

raw_haplo <- read.table(file="batch_17d.haplotypes.tsv", header=TRUE, sep="\t", as.is=TRUE, 
                        row.names=1, na.strings = "consensus") 


working_haplo <- na.omit(raw_haplo)
working_haplo$Cnt <- NULL

# convert snps to ambiguities

working_haplo[working_haplo=="A/T"]<-"W"
working_haplo[working_haplo=="A/C"]<-"M"
working_haplo[working_haplo=="A/G"]<-"R"
working_haplo[working_haplo=="C/T"]<-"Y"
working_haplo[working_haplo=="C/A"]<-"M"
working_haplo[working_haplo=="C/G"]<-"S"
working_haplo[working_haplo=="G/T"]<-"K"
working_haplo[working_haplo=="G/C"]<-"S"
working_haplo[working_haplo=="G/A"]<-"R"
working_haplo[working_haplo=="T/A"]<-"W"
working_haplo[working_haplo=="T/C"]<-"Y"
working_haplo[working_haplo=="T/G"]<-"K"
working_haplo[working_haplo=="-"]<-"?"

transposed <- t(working_haplo)
transposed <- data.frame(transposed)
library(stringr)


###  trimming down of the data sets  ###

###  first data set, 239 loci  ###

locus_ID <- colnames(RAD_239)
loci_239 <- transposed[locus_ID]
loci_239$space <- " "
loci_239 <- loci_239[,c(ncol(loci_239),1:(ncol(loci_239)-1))]

nexus <- paste('#NEXUS
               BEGIN DATA;
               DIMENSIONS NTAX=96 NCHAR=',ncol(loci_239)-1,';
               FORMAT DATATYPE=DNA MISSING=? GAP=- 	INTERLEAVE;
               MATRIX')



nexus_next <- ";
END;"


#################################

write.table(nexus, file="snps_239.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)

write.table(loci_239, file="snps_239.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=TRUE, append=TRUE, quote=FALSE)

write.table(nexus_next, file="snps_239.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE)

##################


###  third data set, 1180 loci  ###

locus_ID <- colnames(RAD_1180)
loci_1180 <- transposed[locus_ID]
loci_1180$space <- " "
loci_1180 <- loci_1180[,c(ncol(loci_1180),1:(ncol(loci_1180)-1))]

nexus <- paste('#NEXUS
               BEGIN DATA;
               DIMENSIONS NTAX=96 NCHAR=',ncol(loci_1180)-1,';
               FORMAT DATATYPE=DNA MISSING=? GAP=- 	INTERLEAVE;
               MATRIX')



nexus_next <- ";
END;"


#################################

write.table(nexus, file="snps_1180.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)

write.table(loci_1180, file="snps_1180.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=TRUE, append=TRUE, quote=FALSE)

write.table(nexus_next, file="snps_1180.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE)

##################

###  fourth data set, 2317 loci  ###

locus_ID <- colnames(RAD_2317)
loci_2317 <- transposed[locus_ID]
loci_2317$space <- " "
loci_2317 <- loci_2317[,c(ncol(loci_2317),1:(ncol(loci_2317)-1))]

nexus <- paste('#NEXUS
               BEGIN DATA;
               DIMENSIONS NTAX=96 NCHAR=',ncol(loci_2317)-1,';
               FORMAT DATATYPE=DNA MISSING=? GAP=- 	INTERLEAVE;
               MATRIX')



nexus_next <- ";
END;"


#################################

write.table(nexus, file="snps_2317.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)

write.table(loci_2317, file="snps_2317.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=TRUE, append=TRUE, quote=FALSE)

write.table(nexus_next, file="snps_2317.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE)

##################

###  fifth data set, 3831 loci  ###

locus_ID <- colnames(RAD_3831)
loci_3831 <- transposed[locus_ID]
loci_3831$space <- " "
loci_3831 <- loci_3831[,c(ncol(loci_3831),1:(ncol(loci_3831)-1))]

nexus <- paste('#NEXUS
               BEGIN DATA;
               DIMENSIONS NTAX=96 NCHAR=',ncol(loci_3831)-1,';
               FORMAT DATATYPE=DNA MISSING=? GAP=- 	INTERLEAVE;
               MATRIX')



nexus_next <- ";
END;"


#################################

write.table(nexus, file="snps_3831.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)

write.table(loci_3831, file="snps_3831.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=TRUE, append=TRUE, quote=FALSE)

write.table(nexus_next, file="snps_3831.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE)

##################

###  sixth data set, 6255 loci  ###

locus_ID <- colnames(RAD_6255)
loci_6255 <- transposed[locus_ID]
loci_6255$space <- " "
loci_6255 <- loci_6255[,c(ncol(loci_6255),1:(ncol(loci_6255)-1))]

nexus <- paste('#NEXUS
               BEGIN DATA;
               DIMENSIONS NTAX=96 NCHAR=',ncol(loci_6255)-1,';
               FORMAT DATATYPE=DNA MISSING=? GAP=- 	INTERLEAVE;
               MATRIX')



nexus_next <- ";
END;"


#################################

write.table(nexus, file="snps_6255.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)

write.table(loci_6255, file="snps_6255.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=TRUE, append=TRUE, quote=FALSE)

write.table(nexus_next, file="snps_6255.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE)

##################

###  seventh data set, 25198 loci  ###

locus_ID <- colnames(RAD_25198)
loci_25198 <- transposed[locus_ID]
loci_25198$space <- " "
loci_25198 <- loci_25198[,c(ncol(loci_25198),1:(ncol(loci_25198)-1))]

nexus <- paste('#NEXUS
               BEGIN DATA;
               DIMENSIONS NTAX=96 NCHAR=',ncol(loci_25198)-1,';
               FORMAT DATATYPE=DNA MISSING=? GAP=- 	INTERLEAVE;
               MATRIX')



nexus_next <- ";
END;"


#################################

write.table(nexus, file="snps_25198.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)

write.table(loci_25198, file="snps_25198.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=TRUE, append=TRUE, quote=FALSE)

write.table(nexus_next, file="snps_25198.nexus", sep = "", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE)

##################





