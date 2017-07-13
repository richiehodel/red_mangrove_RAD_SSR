####################################################################################
### Script by Richie Hodel #########################################################
####################################################################################

####  This script samples loci from the largest dataset 100 times for each 
####  marker type (RAD_25198 or SSR_8). The script outputs histograms showing
####  the distributions of the 100 sampled loci for each dataset, with a 
####  95% confidence interval.  This script generates figures for Fst, Fis, 
####  and Ho, although only Fst figures are shown in the manuscript.

library("hierfstat")

setwd("/Users/richiehodel/Documents/mangroves/RAD-Seq_vs_msat/red_mangrove_RAD_SSR/")


######   msat supbsampling   ####
###
SSR_8 <- read.table(file="Rhiman96_8loci.txt", skip=14, header = TRUE, sep="", na.strings="00")
SSR_8$Population <- NULL
SSR_8$Individual <- NULL
SSR_test <- cbind(heading, SSR_8)

SSR_stats <- basic.stats(SSR_test)
SSR_stats$overall

### sampling 6 loci

j = 1
for (j in 1:100) {
  
heading <- read.csv(file="head.csv", header=TRUE)
samp <- sample(SSR_8, 6, replace=FALSE)
SSR_set <- cbind(heading, samp)

f_stats_SSR6 <- basic.stats(SSR_set)
summary_SSR6 <- f_stats_SSR6$overall
write.table(summary_SSR6[7], file="output_fst_SSR6.csv", sep=",", append=TRUE, col.names = FALSE)
write.table(summary_SSR6[1], file="output_ho_SSR6.csv", sep=",", append=TRUE, col.names = FALSE)
write.table(summary_SSR6[2], file="output_hs_SSR6.csv", sep=",", append=TRUE, col.names = FALSE)
rm(samp, SSR_set)

}

###

### sampling 7 loci

j = 1
for (j in 1:100) {
  
  heading <- read.csv(file="head.csv", header=TRUE)
  samp <- sample(SSR_8, 7, replace=FALSE)
  SSR_set <- cbind(heading, samp)
  f_stats_SSR7 <- basic.stats(SSR_set)
  summary_SSR7 <- f_stats_SSR7$overall
  write.table(summary_SSR7[7], file="output_fst_SSR7.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(summary_SSR7[1], file="output_ho_SSR7.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(summary_SSR7[2], file="output_hs_SSR7.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(samp, SSR_set)
  
}

#######

###  bootstrapping the 8 SSR loci to get confidence interval

j = 1
for (j in 1:100) {
  
  heading <- read.csv(file="head.csv", header=TRUE)
  samp <- sample(SSR_8, 8, replace=TRUE)
  SSR_set <- cbind(heading, samp)
  
  f_stats_SSR8 <- basic.stats(SSR_set)
  summary_SSR8 <- f_stats_SSR8$overall
  write.table(summary_SSR8[7], file="output_fst_SSR8.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(summary_SSR8[1], file="output_ho_SSR8.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(summary_SSR8[2], file="output_hs_SSR8.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(samp, SSR_set)
  
}



######
###   making the histograms for SSRs
######

SSR_fst <- read.csv(file="output_fst_SSR8.csv", row.names=NULL,header=FALSE)
SSR_ho <- read.csv(file="output_ho_SSR8.csv", row.names=NULL,header=FALSE)
SSR_hs <- read.csv(file="output_hs_SSR8.csv", row.names=NULL,header=FALSE)

SSR_fst_7 <- read.csv(file="output_fst_SSR7.csv", row.names=NULL,header=FALSE)
SSR_ho_7 <- read.csv(file="output_ho_SSR7.csv", row.names=NULL,header=FALSE)
SSR_hs_7 <- read.csv(file="output_hs_SSR7.csv", row.names=NULL,header=FALSE)

SSR_fst_6 <- read.csv(file="output_fst_SSR6.csv", row.names=NULL,header=FALSE)
SSR_ho_6 <- read.csv(file="output_ho_SSR6.csv", row.names=NULL,header=FALSE)
SSR_hs_6 <- read.csv(file="output_hs_SSR6.csv", row.names=NULL,header=FALSE)


#### output histograms for SSRs with 7 sampled loci

setEPS()
postscript("Sampled_hs_SSR7.eps")
hist(SSR_hs_7$V2, 35, xlim=c(0, 0.7), xlab=('Hs'),
     main="7 SSR Loci Sampled")
#abline(v=summary_SSR7[2], col="darkorange1", lty=1, lwd=1.5)
abline(v=SSR_stats$overall[2], col="blue", lty=1, lwd=1.5)
#abline(v=max(SSR_hs_7$V2), col="darkorange1", lty=2, lwd=1)
#abline(v=min(SSR_hs_7$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(SSR_hs$V2), col="blue", lty=2, lwd=1)
abline(v=min(SSR_hs$V2), col="blue", lty=2, lwd=1)
dev.off()

setEPS()
postscript("Sampled_ho_SSR7.eps")
hist(SSR_ho_7$V2, 35, xlim=c(0, 0.7), xlab=('Ho'),
     main="7 SSR Loci Sampled")
#abline(v=summary_SSR7[1], col="darkorange1", lty=1, lwd=1.5)
abline(v=SSR_stats$overall[1], col="blue", lty=1, lwd=1.5)
#abline(v=max(SSR_ho_7$V2), col="darkorange1", lty=2, lwd=1)
#abline(v=min(SSR_ho_7$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(SSR_ho$V2), col="blue", lty=2, lwd=1)
abline(v=min(SSR_ho$V2), col="blue", lty=2, lwd=1)
dev.off()

setEPS()
postscript("Sampled_fst_SSR7.eps")
hist(SSR_fst_7$V2, 35, xlim=c(0, 0.3), xlab=('Fst'),
     main="7 SSR Loci Sampled")
#abline(v=summary_SSR7[7], col="darkorange1", lty=1, lwd=1.5)
abline(v=SSR_stats$overall[7], col="blue", lty=1, lwd=1.5)
#abline(v=max(SSR_fst_7$V2), col="darkorange1", lty=2, lwd=1)
#abline(v=min(SSR_fst_7$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(SSR_fst$V2), col="blue", lty=2, lwd=1)
abline(v=min(SSR_fst$V2), col="blue", lty=2, lwd=1)
dev.off()


#### output histograms for SSRs with 6 sampled loci

setEPS()
postscript("Sampled_hs_SSR6.eps")
hist(SSR_hs_6$V2, 35, xlim=c(0, 0.7), xlab=('Hs'),
     main="6 SSR Loci Sampled")
#abline(v=summary_SSR6[2], col="darkorange1", lty=1, lwd=1.5)
abline(v=SSR_stats$overall[2], col="blue", lty=1, lwd=1.5)
#abline(v=max(SSR_hs_6$V2), col="darkorange1", lty=2, lwd=1)
#abline(v=min(SSR_hs_6$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(SSR_hs$V2), col="blue", lty=2, lwd=1)
abline(v=min(SSR_hs$V2), col="blue", lty=2, lwd=1)
dev.off()

setEPS()
postscript("Sampled_ho_SSR6.eps")
hist(SSR_ho_6$V2, 35, xlim=c(0, 0.7), xlab=('Ho'),
     main="6 SSR Loci Sampled")
#abline(v=summary_SSR6[1], col="darkorange1", lty=1, lwd=1.5)
abline(v=SSR_stats$overall[1], col="blue", lty=1, lwd=1.5)
#abline(v=max(SSR_ho_6$V2), col="darkorange1", lty=2, lwd=1)
#abline(v=min(SSR_ho_6$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(SSR_ho$V2), col="blue", lty=2, lwd=1)
abline(v=min(SSR_ho$V2), col="blue", lty=2, lwd=1)
dev.off()

setEPS()
postscript("Sampled_fst_SSR6.eps")
hist(SSR_fst_6$V2, 35, xlim=c(0, 0.3), xlab=('Fst'),
     main="6 SSR Loci Sampled")
#abline(v=summary_SSR6[7], col="darkorange1", lty=1, lwd=1.5)
abline(v=SSR_stats$overall[7], col="blue", lty=1, lwd=1.5)
#abline(v=max(SSR_fst_6$V2), col="darkorange1", lty=2, lwd=1)
#abline(v=min(SSR_fst_6$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(SSR_fst$V2), col="blue", lty=2, lwd=1)
abline(v=min(SSR_fst$V2), col="blue", lty=2, lwd=1)
dev.off()

######



### import raw RAD-Seq data in generic .gdv format
### these data were produced by STACKS in genepop format
### the data were converted to .gdv for convenience

sub_raw  <- contaminant_microbial_human_removed
all_data <- sub_raw

#######    sampling   #########

j = 1
for (j in 1:100) {
  
  heading <- read.csv(file="head.csv", header=TRUE)
  samp <- sample(all_data, 239, replace=FALSE)
  set <- cbind(heading, samp)
  
  f_stats_239 <- basic.stats(set)
  summary_239 <- f_stats_239$overall
  write.table(summary_239[7], file="output_fst_239.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(summary_239[1], file="output_ho_239.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(summary_239[2], file="output_hs_239.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(samp, set)

}


#######

j <- 1
for (j in 1:100) {
  
  heading <- read.csv(file="head.csv", header=TRUE)
  samp <- sample(all_data, 1180, replace=FALSE)
  set <- cbind(heading, samp)
  
  f_stats_1180 <- basic.stats(set)
  summary_1180 <- f_stats_1180$overall
  write.table(summary_1180[7], file="output_fst_1180.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(summary_1180[1], file="output_ho_1180.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(summary_1180[2], file="output_hs_1180.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(samp, set)
  
}

#######


j <- 1
for (j in 1:100) {
  
  heading <- read.csv(file="head.csv", header=TRUE)
  samp <- sample(all_data, 2317, replace=FALSE)
  set <- cbind(heading, samp)
  
  f_stats_2317 <- basic.stats(set)
  summary_2317 <- f_stats_2317$overall
  write.table(summary_2317[7], file="output_fst_2317.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(summary_2317[1], file="output_ho_2317.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(summary_2317[2], file="output_hs_2317.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(samp, set)
  
}

#######

j <- 1
for (j in 1:100) {
  
  heading <- read.csv(file="head.csv", header=TRUE)
  samp <- sample(all_data, 3831, replace=FALSE)
  set <- cbind(heading, samp)
  
  f_stats_3831 <- basic.stats(set)
  summary_3831 <- f_stats_3831$overall
  write.table(summary_3831[7], file="output_fst_3831.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(summary_3831[1], file="output_ho_3831.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(summary_3831[2], file="output_hs_3831.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(samp, set)
  
}

#######

j <- 1
for (j in 1:100) {
  
  heading <- read.csv(file="head.csv", header=TRUE)
  samp <- sample(all_data, 6255, replace=FALSE)
  set <- cbind(heading, samp)
  
  f_stats_6255 <- basic.stats(set)
  summary_6255 <- f_stats_6255$overall
  write.table(summary_6255[7], file="output_fst_6255.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(summary_6255[1], file="output_ho_6255.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(summary_6255[2], file="output_hs_6255.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(samp, set)
  
}

#######
##########    bootstrapping for confidence intervals    239    ###########

data <- RAD_239

j = 1
for (j in 1:100) {
  
  heading <- read.csv(file="head.csv", header=TRUE)
  boot <- sample(data, length(data), replace=TRUE)
  bootstrapped <- cbind(heading, boot)
  
  f_stats_boot_239 <- basic.stats(bootstrapped)
  boot_summary_239 <- f_stats_boot_239$overall
  write.table(boot_summary_239[7], file="output_boot_fst_239.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(boot_summary_239[1], file="output_boot_ho_239.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(boot_summary_239[2], file="output_boot_hs_239.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(boot, bootstrapped)
  
}


##########    bootstrapping for confidence intervals    1180      ###########

data <- RAD_1180

j = 1
for (j in 1:100) {
  
  heading <- read.csv(file="head.csv", header=TRUE)
  boot <- sample(data, length(data), replace=TRUE)
  bootstrapped <- cbind(heading, boot)
  
  f_stats_boot_1180 <- basic.stats(bootstrapped)
  boot_summary_1180 <- f_stats_boot_1180$overall
  write.table(boot_summary_1180[7], file="output_boot_fst_1180.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(boot_summary_1180[1], file="output_boot_ho_1180.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(boot_summary_1180[2], file="output_boot_hs_1180.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(boot, bootstrapped)
  
}

##########    bootstrapping for confidence intervals     2317       ###########

data <- RAD_2317

j = 1
for (j in 1:100) {
  
  heading <- read.csv(file="head.csv", header=TRUE)
  boot <- sample(data, length(data), replace=TRUE)
  bootstrapped <- cbind(heading, boot)
  
  f_stats_boot_2317 <- basic.stats(bootstrapped)
  boot_summary_2317 <- f_stats_boot_2317$overall
  write.table(boot_summary_2317[7], file="output_boot_fst_2317.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(boot_summary_2317[1], file="output_boot_ho_2317.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(boot_summary_2317[2], file="output_boot_hs_2317.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(boot, bootstrapped)
  
}

##########    bootstrapping for confidence intervals      3831             ###########

data <- RAD_3831

j = 1
for (j in 1:100) {
  
  heading <- read.csv(file="head.csv", header=TRUE)
  boot <- sample(data, length(data), replace=TRUE)
  bootstrapped <- cbind(heading, boot)
  
  f_stats_boot_3831 <- basic.stats(bootstrapped)
  boot_summary_3831 <- f_stats_boot_3831$overall
  write.table(boot_summary_3831[7], file="output_boot_fst_3831.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(boot_summary_3831[1], file="output_boot_ho_3831.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(boot_summary_3831[2], file="output_boot_hs_3831.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(boot, bootstrapped)
  
}

##########    bootstrapping for confidence intervals         6255          ###########

data <- RAD_6255

j = 1
for (j in 1:100) {
  
  heading <- read.csv(file="head.csv", header=TRUE)
  boot <- sample(data, length(data), replace=TRUE)
  bootstrapped <- cbind(heading, boot)
  
  f_stats_boot_6255 <- basic.stats(bootstrapped)
  boot_summary_6255 <- f_stats_boot_6255$overall
  write.table(boot_summary_6255[7], file="output_boot_fst_6255.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(boot_summary_6255[1], file="output_boot_ho_6255.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(boot_summary_6255[2], file="output_boot_hs_6255.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(boot, bootstrapped)
  
}

##########    bootstrapping for confidence intervals           25198            ###########

data <- RAD_25198

j = 1
for (j in 1:100) {
  
  heading <- read.csv(file="head.csv", header=TRUE)
  boot <- sample(data, length(data), replace=TRUE)
  bootstrapped <- cbind(heading, boot)
  
  f_stats_boot_25198 <- basic.stats(bootstrapped)
  boot_summary_25198 <- f_stats_boot_25198$overall
  write.table(boot_summary_25198[7], file="output_boot_fst_25198.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(boot_summary_25198[1], file="output_boot_ho_25198.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(boot_summary_25198[2], file="output_boot_hs_25198.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(boot, bootstrapped)
  
}


#######
#######
#######























#######


boot_25198_fst <- read.csv(file="output_boot_fst_25198.csv",row.names=NULL,header=FALSE)
boot_25198_ho <- read.csv(file="output_boot_ho_25198.csv",row.names=NULL,header=FALSE)
boot_25198_hs <- read.csv(file="output_boot_hs_25198.csv",row.names=NULL,header=FALSE)


#######

loci_239_fst <- read.csv(file="output_fst_239.csv", row.names=NULL,header=FALSE)
loci_239_ho <- read.csv(file="output_ho_239.csv", row.names=NULL,header=FALSE)
loci_239_hs <- read.csv(file="output_hs_239.csv", row.names=NULL,header=FALSE)
boot_239_fst <- read.csv(file="output_boot_fst_239.csv",row.names=NULL,header=FALSE)
boot_239_ho <- read.csv(file="output_boot_ho_239.csv",row.names=NULL,header=FALSE)
boot_239_hs <- read.csv(file="output_boot_hs_239.csv",row.names=NULL,header=FALSE)

setEPS()
postscript("Sampled_fst_239.eps")
hist(loci_239_fst$V2, 35, xlim=c(0, 0.3), xlab=('Fst'),
     main="239 RAD Loci Sampled")
abline(v=stats_239$overall[7], col="darkorange1", lty=1, lwd=1.5)
abline(v=stats_25198$overall[7], col="blue", lty=1, lwd=1.5)
abline(v=max(boot_239_fst$V2), col="darkorange1", lty=2, lwd=1)
abline(v=min(boot_239_fst$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(boot_25198_fst$V2), col="blue", lty=2, lwd=1)
abline(v=min(boot_25198_fst$V2), col="blue", lty=2, lwd=1)
dev.off()

setEPS()
postscript("Sampled_ho_239.eps")
hist(loci_239_ho$V2, 35, xlim=c(0, 0.7), xlab=('Ho'),
     main="239 RAD Loci Sampled")
abline(v=stats_239$overall[1], col="darkorange1", lty=1, lwd=1.5)
abline(v=stats_25198$overall[1], col="blue", lty=1, lwd=1.5)
abline(v=max(boot_239_ho$V2), col="darkorange1", lty=2, lwd=1)
abline(v=min(boot_239_ho$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(boot_25198_ho$V2), col="blue", lty=2, lwd=1)
abline(v=min(boot_25198_ho$V2), col="blue", lty=2, lwd=1)
dev.off()

setEPS()
postscript("Sampled_hs_239.eps")
hist(loci_239_hs$V2, 35, xlim=c(0, 0.7), xlab=('Hs'),
     main="239 RAD Loci Sampled")
abline(v=stats_239$overall[2], col="darkorange1", lty=1, lwd=1.5)
abline(v=stats_25198$overall[2], col="blue", lty=1, lwd=1.5)
abline(v=max(boot_239_hs$V2), col="darkorange1", lty=2, lwd=1)
abline(v=min(boot_239_hs$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(boot_25198_hs$V2), col="blue", lty=2, lwd=1)
abline(v=min(boot_25198_hs$V2), col="blue", lty=2, lwd=1)
dev.off()


#######
# 1180 loci
#######

loci_1180_fst <- read.csv(file="output_fst_1180.csv", row.names=NULL,header=FALSE)
loci_1180_ho <- read.csv(file="output_ho_1180.csv", row.names=NULL,header=FALSE)
loci_1180_hs <- read.csv(file="output_hs_1180.csv", row.names=NULL,header=FALSE)
boot_1180_fst <- read.csv(file="output_boot_fst_1180.csv",row.names=NULL,header=FALSE)
boot_1180_ho <- read.csv(file="output_boot_ho_1180.csv",row.names=NULL,header=FALSE)
boot_1180_hs <- read.csv(file="output_boot_hs_1180.csv",row.names=NULL,header=FALSE)

setEPS()
postscript("Sampled_fst_1180.eps")
hist(loci_1180_fst$V2, 35, xlim=c(0, 0.3), xlab=('Fst'),
     main="1180 RAD Loci Sampled")
abline(v=stats_1180$overall[7], col="darkorange1", lty=1, lwd=1.5)
abline(v=stats_25198$overall[7], col="blue", lty=1, lwd=1.5)
abline(v=max(boot_1180_fst$V2), col="darkorange1", lty=2, lwd=1)
abline(v=min(boot_1180_fst$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(boot_25198_fst$V2), col="blue", lty=2, lwd=1)
abline(v=min(boot_25198_fst$V2), col="blue", lty=2, lwd=1)
dev.off()

setEPS()
postscript("Sampled_ho_1180.eps")
hist(loci_1180_ho$V2, 35, xlim=c(0, 0.7), xlab=('Ho'),
     main="1180 RAD Loci Sampled")
abline(v=stats_1180$overall[1], col="darkorange1", lty=1, lwd=1.5)
abline(v=stats_25198$overall[1], col="blue", lty=1, lwd=1.5)
abline(v=max(boot_1180_ho$V2), col="darkorange1", lty=2, lwd=1)
abline(v=min(boot_1180_ho$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(boot_25198_ho$V2), col="blue", lty=2, lwd=1)
abline(v=min(boot_25198_ho$V2), col="blue", lty=2, lwd=1)
dev.off()

setEPS()
postscript("Sampled_hs_1180.eps")
hist(loci_1180_hs$V2, 35, xlim=c(0, 0.7), xlab=('Hs'),
     main="1180 RAD Loci Sampled")
abline(v=stats_1180$overall[2], col="darkorange1", lty=1, lwd=1.5)
abline(v=stats_25198$overall[2], col="blue", lty=1, lwd=1.5)
abline(v=max(boot_1180_hs$V2), col="darkorange1", lty=2, lwd=1)
abline(v=min(boot_1180_hs$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(boot_25198_hs$V2), col="blue", lty=2, lwd=1)
abline(v=min(boot_25198_hs$V2), col="blue", lty=2, lwd=1)
dev.off()


#######
# 2317 loci
#######

loci_2317_fst <- read.csv(file="output_fst_2317.csv", row.names=NULL,header=FALSE)
loci_2317_ho <- read.csv(file="output_ho_2317.csv", row.names=NULL,header=FALSE)
loci_2317_hs <- read.csv(file="output_hs_2317.csv", row.names=NULL,header=FALSE)
boot_2317_fst <- read.csv(file="output_boot_fst_2317.csv",row.names=NULL,header=FALSE)
boot_2317_ho <- read.csv(file="output_boot_ho_2317.csv",row.names=NULL,header=FALSE)
boot_2317_hs <- read.csv(file="output_boot_hs_2317.csv",row.names=NULL,header=FALSE)

setEPS()
postscript("Sampled_fst_2317.eps")
hist(loci_2317_fst$V2, 20, xlim=c(0, 0.3), xlab=('Fst'),
     main="2317 RAD Loci Sampled")
abline(v=stats_2317$overall[7], col="darkorange1", lty=1, lwd=1.5)
abline(v=stats_25198$overall[7], col="blue", lty=1, lwd=1.5)
abline(v=max(boot_2317_fst$V2), col="darkorange1", lty=2, lwd=1)
abline(v=min(boot_2317_fst$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(boot_25198_fst$V2), col="blue", lty=2, lwd=1)
abline(v=min(boot_25198_fst$V2), col="blue", lty=2, lwd=1)
dev.off()

setEPS()
postscript("Sampled_ho_2317.eps")
hist(loci_2317_ho$V2, 20, xlim=c(0, 0.7), xlab=('Ho'),
     main="2317 RAD Loci Sampled")
abline(v=stats_2317$overall[1], col="darkorange1", lty=1, lwd=1.5)
abline(v=stats_25198$overall[1], col="blue", lty=1, lwd=1.5)
abline(v=max(boot_2317_ho$V2), col="darkorange1", lty=2, lwd=1)
abline(v=min(boot_2317_ho$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(boot_25198_ho$V2), col="blue", lty=2, lwd=1)
abline(v=min(boot_25198_ho$V2), col="blue", lty=2, lwd=1)
dev.off()

setEPS()
postscript("Sampled_hs_2317.eps")
hist(loci_2317_hs$V2, 20, xlim=c(0, 0.7), xlab=('Hs'),
     main="2317 RAD Loci Sampled")
abline(v=stats_2317$overall[2], col="darkorange1", lty=1, lwd=1.5)
abline(v=stats_25198$overall[2], col="blue", lty=1, lwd=1.5)
abline(v=max(boot_2317_hs$V2), col="darkorange1", lty=2, lwd=1)
abline(v=min(boot_2317_hs$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(boot_25198_hs$V2), col="blue", lty=2, lwd=1)
abline(v=min(boot_25198_hs$V2), col="blue", lty=2, lwd=1)
dev.off()


#######
# 3831 loci
#######

loci_3831_fst <- read.csv(file="output_fst_3831.csv", row.names=NULL,header=FALSE)
loci_3831_ho <- read.csv(file="output_ho_3831.csv", row.names=NULL,header=FALSE)
loci_3831_hs <- read.csv(file="output_hs_3831.csv", row.names=NULL,header=FALSE)
boot_3831_fst <- read.csv(file="output_boot_fst_3831.csv",row.names=NULL,header=FALSE)
boot_3831_ho <- read.csv(file="output_boot_ho_3831.csv",row.names=NULL,header=FALSE)
boot_3831_hs <- read.csv(file="output_boot_hs_3831.csv",row.names=NULL,header=FALSE)

setEPS()
postscript("Sampled_fst_3831.eps")
hist(loci_3831_fst$V2, 12, xlim=c(0, 0.3), xlab=('Fst'),
     main="3831 RAD Loci Sampled")
abline(v=stats_3831$overall[7], col="darkorange1", lty=1, lwd=1.5)
abline(v=stats_25198$overall[7], col="blue", lty=1, lwd=1.5)
abline(v=max(boot_3831_fst$V2), col="darkorange1", lty=2, lwd=1)
abline(v=min(boot_3831_fst$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(boot_25198_fst$V2), col="blue", lty=2, lwd=1)
abline(v=min(boot_25198_fst$V2), col="blue", lty=2, lwd=1)
dev.off()

setEPS()
postscript("Sampled_ho_3831.eps")
hist(loci_3831_ho$V2, 12, xlim=c(0, 0.7), xlab=('Ho'),
     main="3831 RAD Loci Sampled")
abline(v=stats_3831$overall[1], col="darkorange1", lty=1, lwd=1.5)
abline(v=stats_25198$overall[1], col="blue", lty=1, lwd=1.5)
abline(v=max(boot_3831_ho$V2), col="darkorange1", lty=2, lwd=1)
abline(v=min(boot_3831_ho$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(boot_25198_ho$V2), col="blue", lty=2, lwd=1)
abline(v=min(boot_25198_ho$V2), col="blue", lty=2, lwd=1)
dev.off()

setEPS()
postscript("Sampled_hs_3831.eps")
hist(loci_3831_hs$V2, 12, xlim=c(0, 0.7), xlab=('Hs'),
     main="3831 RAD Loci Sampled")
abline(v=stats_3831$overall[2], col="darkorange1", lty=1, lwd=1.5)
abline(v=stats_25198$overall[2], col="blue", lty=1, lwd=1.5)
abline(v=max(boot_3831_hs$V2), col="darkorange1", lty=2, lwd=1)
abline(v=min(boot_3831_hs$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(boot_25198_hs$V2), col="blue", lty=2, lwd=1)
abline(v=min(boot_25198_hs$V2), col="blue", lty=2, lwd=1)
dev.off()


#######
# 6255 loci
#######

loci_6255_fst <- read.csv(file="output_fst_6255.csv", row.names=NULL,header=FALSE)
loci_6255_ho <- read.csv(file="output_ho_6255.csv", row.names=NULL,header=FALSE)
loci_6255_hs <- read.csv(file="output_hs_6255.csv", row.names=NULL,header=FALSE)
boot_6255_fst <- read.csv(file="output_boot_fst_6255.csv",row.names=NULL,header=FALSE)
boot_6255_ho <- read.csv(file="output_boot_ho_6255.csv",row.names=NULL,header=FALSE)
boot_6255_hs <- read.csv(file="output_boot_hs_6255.csv",row.names=NULL,header=FALSE)

setEPS()
postscript("Sampled_fst_6255.eps")
hist(loci_6255_fst$V2, 10, xlim=c(0, 0.3), xlab=('Fst'),
     main="6255 RAD Loci Sampled")
abline(v=stats_6255$overall[7], col="darkorange1", lty=1, lwd=1.5)
abline(v=stats_25198$overall[7], col="blue", lty=1, lwd=1.5)
abline(v=max(boot_6255_fst$V2), col="darkorange1", lty=2, lwd=1)
abline(v=min(boot_6255_fst$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(boot_25198_fst$V2), col="blue", lty=2, lwd=1)
abline(v=min(boot_25198_fst$V2), col="blue", lty=2, lwd=1)
dev.off()

setEPS()
postscript("Sampled_ho_6255.eps")
hist(loci_6255_ho$V2, 10, xlim=c(0, 0.7), xlab=('Ho'),
     main="6255 RAD Loci Sampled")
abline(v=stats_6255$overall[1], col="darkorange1", lty=1, lwd=1.5)
abline(v=stats_25198$overall[1], col="blue", lty=1, lwd=1.5)
abline(v=max(boot_6255_ho$V2), col="darkorange1", lty=2, lwd=1)
abline(v=min(boot_6255_ho$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(boot_25198_ho$V2), col="blue", lty=2, lwd=1)
abline(v=min(boot_25198_ho$V2), col="blue", lty=2, lwd=1)
dev.off()

setEPS()
postscript("Sampled_hs_6255.eps")
hist(loci_6255_hs$V2, 10, xlim=c(0, 0.7), xlab=('Hs'),
     main="6255 RAD Loci Sampled")
abline(v=stats_6255$overall[2], col="darkorange1", lty=1, lwd=1.5)
abline(v=stats_25198$overall[2], col="blue", lty=1, lwd=1.5)
abline(v=max(boot_6255_hs$V2), col="darkorange1", lty=2, lwd=1)
abline(v=min(boot_6255_hs$V2), col="darkorange1", lty=2, lwd=1)
abline(v=max(boot_25198_hs$V2), col="blue", lty=2, lwd=1)
abline(v=min(boot_25198_hs$V2), col="blue", lty=2, lwd=1)
dev.off()





