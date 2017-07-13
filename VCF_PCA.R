####################################################################################
### Script by Richie Hodel #########################################################
####################################################################################

## This script uses a VCF file to conduct a PCA for each dataset.  The PCAs for 
## the RAD data were completed in SNPRelate, and the PCA for the SSR data was
## run in Genodive.


source("http://bioconductor.org/biocLite.R")
biocLite("gdsfmt")
biocLite("SNPRelate")

library("gdsfmt")
library("SNPRelate")

library(vegan)
library(ade4)
library(MASS)

getwd()
setwd("/Users/richiehodel/Documents/mangroves/RAD-Seq_vs_msat/red_mangrove_RAD_SSR/")

vcf_header1 <- paste('##fileformat=VCFv4.0 
##fileDate=20160915
##source="Stacks v1.35"
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Allele Depth">
##FORMAT=<ID=GL,Number=.,Type=Float,Description="Genotype Likelihood">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	MID34_BHKFl_R2	MID46_BHKFl_R5	MID67_BHKFl_R4	MID71_BHKFl_R7	MID83_BHKFl_R8	MID88_BHKFl_R6	MID91_BHKFl_R1	MID95_BHKFl_R9	ILL_8nt_10_CpCFl_R1	ILL_8nt_111_CpCFl_R2	ILL_8nt_14_CpCFl_R3	ILL_8nt_27_CpCFl_R5	ILL_8nt_39_CpCFl_R6	ILL_8nt_5_CpCFl_R7	ILL_8nt_62_CpCFl_R8	ILL_8nt_74_CpCFl_R9	ILL_8nt_85_CvPFl_R1	ILL_8nt_97_CvPFl_R2	ILL_9nt_109_CvPFl_R3	ILL_9nt_12_CvPFl_R4	ILL_9nt_24_CvPFl_R5	ILL_9nt_35_CvPFl_R6	ILL_9nt_46_CvPFl_R7	ILL_9nt_60_CvPFl_R8	ILL_8nt_102_HwdFl_R1	ILL_8nt_114_HwdFl_R2	ILL_8nt_17_HwdFl_R3	ILL_8nt_3_HwdFl_R4	ILL_8nt_41_HwdFl_R5	ILL_8nt_66_HwdFl_R7	ILL_8nt_77_HwdFl_R10	MID50_HwdFl_R6	ILL_8nt_87_IsmFl_R1	ILL_9nt_10_IsmFl_R2	ILL_9nt_110_IsmFl_R3	ILL_9nt_13_IsmFl_R4	ILL_9nt_26_IsmFl_R5	ILL_9nt_48_IsmFl_R9	MID14_IsmFl_R10	MID43_IsmFl_R6	ILL_8nt_118_KyLFl_R4	ILL_8nt_98_KyLFl_R2	ILL_9nt_11_KyLFl_R3	ILL_9nt_25_KyLFl_R6	ILL_9nt_36_KyLFl_R7	ILL_9nt_47_KyLFl_R9	MID1_KyLFl_R10	MID28_KyLFl_R5	ILL_8nt_100_MlBFl_R1	ILL_8nt_112_MlBFl_R2	ILL_8nt_15_MlBFl_R5	ILL_8nt_28_MlBFl_R6	ILL_8nt_4_MlBFl_R7	ILL_8nt_50_MlBFl_R8	ILL_8nt_63_MlBFl_R9	ILL_8nt_75_MlBFl_R10	MID11_NPRFl_R2	MID24_NPRFl_R3	MID36_NPRFl_R4	MID48_NPRFl_R5	MID61_NPRFl_R6	MID73_NPRFl_R7	MID76_NPRFl_R8	MID85_NPRFl_R9	MID13_ShKFl_R1	MID25_ShKFl_R2	MID37_ShKFl_R3	MID49_ShKFl_R4	MID62_ShKFl_R5	MID74_ShKFl_R7	MID86_ShKFl_R8	MID98_ShKFl_R9	MID10_TCBFl_R1	MID23_TCBFl_R2	MID35_TCBFl_R3	MID47_TCBFl_R5	MID60_TCBFl_R6	MID72_TCBFl_R7	MID84_TCBFl_R8	MID96_TCBFl_R9	ILL_8nt_88_VKyFl_R1	ILL_9nt_100_VKyFl_R2	ILL_9nt_111_VKyFl_R3	ILL_9nt_38_VKyFl_R6	ILL_9nt_49_VKyFl_R7	ILL_9nt_63_VKyFl_R8	MID81_VKyFl_R9	MID93_VKyFl_R5	ILL_8nt_16_WPBFl_R3	ILL_8nt_29_WPBFl_R4	ILL_8nt_40_WPBFl_R5	ILL_8nt_51_WPBFl_R6	ILL_8nt_65_WPBFl_R7	ILL_8nt_76_WPBFl_R10	MID78_WPBFl_R2	MID90_WPBFl_R1')

vcf_header1

vcf_raw <- read.table(file="batch_17d_April2017.vcf", skip=9)

####   vcf is loaded--now need to remove loci that were contaminant blasts

vcf_working <- vcf_raw

vcf_working$X <- "X"
vcf_working$keep <- vcf_working$V3
vcf_working$index <- paste(vcf_working$X, vcf_working$keep, sep='')

vcf_working <- data.frame(t(vcf_working))
use_loci <- colnames(RAD_25198)

vcf_loci <- vcf_working[, use_loci]

vcf_working <- data.frame(t(vcf_working))

vcf_working$X <- NULL
vcf_working$index <- NULL
vcf_working$keep <- NULL

vcf_25198 <- vcf_working

library(stringr)

split_1 <- str_split_fixed(vcf_working$V8, ";", 2)
split_1 <- data.frame(split_1)
split_1$X2 <- NULL
split_2 <- str_split_fixed(split_1$X1, "=", 2)
split_3 <- as.numeric(split_2[,2])
vcf_working$X2 <- split_3


vcf_239 <- subset(vcf_working, vcf_working$X2>=83)
vcf_1180 <- subset(vcf_working, vcf_working$X2>=75)
vcf_2317 <- subset(vcf_working, vcf_working$X2>=65)
vcf_3831 <- subset(vcf_working, vcf_working$X2>=50)
vcf_6255 <- subset(vcf_working, vcf_working$X2>=30)

vcf_239$X2 <- NULL
write.table(vcf_header1, file="239.vcf", eol= "\n", quote=FALSE,
            col.names=FALSE, row.names=FALSE, append=FALSE)
write.table(vcf_239, file="239.vcf", sep = "\t", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE)

vcf_1180$X2 <- NULL
write.table(vcf_header1, file="1180.vcf", eol= "\n", quote=FALSE,
            col.names=FALSE, row.names=FALSE, append=FALSE)
write.table(vcf_1180, file="1180.vcf", sep = "\t", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE)

vcf_2317$X2 <- NULL
write.table(vcf_header1, file="2317.vcf", eol= "\n", quote=FALSE,
            col.names=FALSE, row.names=FALSE, append=FALSE)
write.table(vcf_2317, file="2317.vcf", sep = "\t", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE)

vcf_3831$X2 <- NULL
write.table(vcf_header1, file="3831.vcf", eol= "\n", quote=FALSE,
            col.names=FALSE, row.names=FALSE, append=FALSE)
write.table(vcf_3831, file="3831.vcf", sep = "\t", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE)

vcf_6255$X2 <- NULL
write.table(vcf_header1, file="6255.vcf", eol= "\n", quote=FALSE,
            col.names=FALSE, row.names=FALSE, append=FALSE)
write.table(vcf_6255, file="6255.vcf", sep = "\t", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE)

vcf_25198$X2 <- NULL
write.table(vcf_header1, file="25198.vcf", eol= "\n", quote=FALSE,
            col.names=FALSE, row.names=FALSE, append=FALSE)
write.table(vcf_25198, file="25198.vcf", sep = "\t", eol= "\n",
            col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE)


####### SNPRelate  #############################################
################################################################

###  239 loci  ###

vcf.fn239 <- "/Users/richiehodel/Documents/mangroves/RAD-Seq_vs_msat/red_mangrove_RAD_SSR/239.vcf"
snpgdsVCF2GDS(vcf.fn239, "239.gds", method="biallelic.only")
snpgdsSummary("239.gds")
genofile <- snpgdsOpen("239.gds")
print(genofile)

###  1180 loci  ###

vcf.fn1180 <- "/Users/richiehodel/Documents/mangroves/RAD-Seq_vs_msat/red_mangrove_RAD_SSR/1180.vcf"
snpgdsVCF2GDS(vcf.fn1180, "1180a.gds", method="biallelic.only")
snpgdsSummary("1180a.gds")
genofile <- snpgdsOpen("1180a.gds")
print(genofile)

###  2317 loci  ###

vcf.fn2317 <- "/Users/richiehodel/Documents/mangroves/RAD-Seq_vs_msat/red_mangrove_RAD_SSR/2317.vcf"
snpgdsVCF2GDS(vcf.fn2317, "2317e.gds", method="biallelic.only")
snpgdsSummary("2317e.gds")
genofile <- snpgdsOpen("2317e.gds")
print(genofile)

###  3831 loci  ###

vcf.fn3831 <- "/Users/richiehodel/Documents/mangroves/RAD-Seq_vs_msat/red_mangrove_RAD_SSR/3831.vcf"
snpgdsVCF2GDS(vcf.fn3831, "3831c.gds", method="biallelic.only")
snpgdsSummary("3831c.gds")
genofile <- snpgdsOpen("3831c.gds")
print(genofile)

###  6255 loci  ###

vcf.fn6255 <- "/Users/richiehodel/Documents/mangroves/RAD-Seq_vs_msat/red_mangrove_RAD_SSR/6255.vcf"
snpgdsVCF2GDS(vcf.fn6255, "6255b.gds", method="biallelic.only")
snpgdsSummary("6255b.gds")
genofile <- snpgdsOpen("6255b.gds")
print(genofile)

###  25198 loci  ###

vcf.fn25198 <- "/Users/richiehodel/Documents/mangroves/RAD-Seq_vs_msat/red_mangrove_RAD_SSR/25198.vcf"
snpgdsVCF2GDS(vcf.fn25198, "25198d.gds", method="biallelic.only")
snpgdsSummary("25198d.gds")
genofile <- snpgdsOpen("25198d.gds")
print(genofile)


#########################
############
######

###   PCAs   ###

pca <- snpgdsPCA(genofile, autosome.only = FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

###  need to use appropriate popmap; for publication use popmap2.txt  ###

color_list <- c("blue","black","darkorange1")

pop_code <- scan("popmap2.txt", what=character())

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop <- cbind(sample.id, pop_code)

tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

pop_code <- as.factor(pop_code)

tab$pop <- pop_code

###  need to change title (main =) and file name  ###

setEPS()
postscript('pca25198.eps')
#postscript('pca6255.eps')
#postscript('pca3831.eps')
#postscript('pca2317.eps')
#postscript('pca1180.eps')
#postscript('pca239.eps')

plot(tab$EV2, tab$EV1, col=color_list[as.integer(tab$pop)], main="25198 RAD Loci",
     xlab="Eigenvector 2", ylab="Eigenvector 1",
     xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))
legend("topright", legend=levels(tab$pop), pch="o", col=c("blue","black","darkorange1"), cex=1)

dev.off()


#### msat

msat_tab <- tab

msat_tab$pop <- pop_code

msat_input <- read.table(file="PCA_msat.gdv", skip=1, sep='\t')

msat_tab$EV1 <- msat_input$V2

msat_tab$EV2 <- msat_input$V3

setEPS()
postscript('pca_msat8.eps')
plot(msat_tab$EV2, msat_tab$EV1, col= color_list[as.integer(msat_tab$pop)], 
     main="8 SSR Loci", 
     xlab="Eigenvector 2", ylab="Eigenvector 1",
     xlim=c(-1,1), ylim=c(-1,1))
legend("topright", legend=levels(msat_tab$pop), pch="o",
       col=c("blue","black","darkorange1"), cex=1)
dev.off()






