####################################################################################
### Script by Richie Hodel #########################################################
####################################################################################

###  This script computes some basic population genetic stats for all the SNP  ###
###  data sets, and prints the stats to files  ###


###  install necessary packages  ###
install.packages("hierfstat")
library("hierfstat")

####################################

### Run the script Data_acquisition.R first so that all
### of the variables will exist

genodive <- data.frame(raw$Population, raw$Individual)

genodive$Population <- genodive$raw.Population
genodive$Individual <- genodive$raw.Individual
genodive$raw.Individual <- NULL
genodive$raw.Population <- NULL

######
# you will need to change directories as needed

setwd("/Users/richiehodel/Documents/mangroves/RAD-Seq_vs_msat/red_mangrove_RAD_SSR/")

head <- read.csv(file="head.csv", header=TRUE)

head_239 <- cbind(head, RAD_239)
stats_239 <- basic.stats(head_239)
stats_239$overall

######

geno_239 <- head_239
geno_239$Population <- NULL
geno_239 <- cbind(genodive, geno_239)
write.table(geno_239, file="239.gdv", quote=FALSE, row.names=FALSE, na="00", sep="\t")

######

head_1180 <- cbind(head, RAD_1180)
stats_1180 <- basic.stats(head_1180)
stats_1180$overall

geno_1180 <- head_1180
geno_1180$Population <- NULL
geno_1180 <- cbind(genodive, geno_1180)
write.table(geno_1180, file="1180.gdv", quote=FALSE, row.names=FALSE, na="00", sep="\t")

head_2317 <- cbind(head, RAD_2317)
stats_2317 <- basic.stats(head_2317)
stats_2317$overall

geno_2317 <- head_2317
geno_2317$Population <- NULL
geno_2317 <- cbind(genodive, geno_2317)
write.table(geno_2317, file="2317.gdv", quote=FALSE, row.names=FALSE, na="00", sep="\t")

head_3831 <- cbind(head, RAD_3831)
stats_3831 <- basic.stats(head_3831)
stats_3831$overall

geno_3831 <- head_3831
geno_3831$Population <- NULL
geno_3831 <- cbind(genodive, geno_3831)
write.table(geno_3831, file="3831.gdv", quote=FALSE, row.names=FALSE, na="00", sep="\t")

head_6255 <- cbind(head, RAD_6255)
stats_6255 <- basic.stats(head_6255)
stats_6255$overall

geno_6255 <- head_6255
geno_6255$Population <- NULL
geno_6255 <- cbind(genodive, geno_6255)
write.table(geno_6255, file="6255.gdv", quote=FALSE, row.names=FALSE, na="00", sep="\t")

head_25198 <- cbind(head, RAD_25198)
stats_25198 <- basic.stats(head_25198)
stats_25198$overall

geno_25198 <- head_25198
geno_25198$Population <- NULL
geno_25198 <- cbind(genodive, geno_25198)
write.table(geno_25198, file="25198.gdv", quote=FALSE, row.names=FALSE, na="00", sep="\t")



SSR_8 <- read.table(file="Rhiman96_8loci.txt", skip=14, header = TRUE, sep="", na.strings="00")
SSR_8$Population <- NULL
SSR_8$Individual <- NULL
SSR_8_stats <- cbind(head, SSR_8)

SSR_8_stats <- basic.stats(SSR_8_stats)
SSR_8_stats$overall

mean(SSR_8_stats$Ho[1:8])



stats <- data.frame(stats_239$overall)
stats <- cbind(stats, stats_1180$overall)
stats <- cbind(stats, stats_2317$overall)
stats <- cbind(stats, stats_3831$overall)
stats <- cbind(stats, stats_6255$overall)
stats <- cbind(stats, stats_25198$overall)
stats <- cbind(stats, SSR_8_stats$overall)

stats <- t(stats)
stats <- data.frame(stats)

row.names(stats) <- c("RAD_239", "RAD_1180", "RAD_2317", "RAD_3831", 
                      "RAD_6255","RAD_25198", "SSR_8")
stats_table <- stats

stats_table$Ht <- NULL
stats_table$Dst <- NULL
stats_table$Htp <- NULL
stats_table$Dstp <- NULL
stats_table$Dest <- NULL


write.csv(stats_table, file="table2_March2017.csv")



##### getting confidence intervals for estimates of basic stats

##### CI for SSR  ######

j = 1
for (j in 1:100) {
  
  head_SSR <- cbind(head, SSR_8)
  
  confidence <- sample(head_SSR, 8, replace=TRUE)
  interval <- cbind(head, confidence)
  
  stats_SSR <- basic.stats(interval)
  write.table(stats_SSR$overall[9], file="output_CI_fis_SSR.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(head_SSR, confidence, interval)
  
}

CI_fst_SSR_input <- read.table("output_CI_fst_SSR.csv", sep=",")
CI_fst_SSR_upper <- max(CI_fst_SSR_input$V2)
CI_fst_SSR_lower <- min(CI_fst_SSR_input$V2)
CI_fst_SSR <- CI_fst_SSR_upper - CI_fst_SSR_lower
CI_fst95_SSR <- (CI_fst_SSR*0.95)/2

CI_ho_SSR_input <- read.table("output_CI_ho_SSR.csv", sep=",")
CI_ho_SSR_upper <- max(CI_ho_SSR_input$V2)
CI_ho_SSR_lower <- min(CI_ho_SSR_input$V2)
CI_ho_SSR <- CI_ho_SSR_upper - CI_ho_SSR_lower
CI_ho95_SSR <- (CI_ho_SSR*0.95)/2

CI_hs_SSR_input <- read.table("output_CI_hs_SSR.csv", sep=",")
CI_hs_SSR_upper <- max(CI_hs_SSR_input$V2)
CI_hs_SSR_lower <- min(CI_hs_SSR_input$V2)
CI_hs_SSR <- CI_hs_SSR_upper - CI_hs_SSR_lower
CI_hs95_SSR <- (CI_hs_SSR*0.95)/2

CI_fis_SSR_input <- read.table("output_CI_fis_SSR.csv", sep=",")
CI_fis_SSR_upper <- max(CI_fis_SSR_input$V2)
CI_fis_SSR_lower <- min(CI_fis_SSR_input$V2)
CI_fis_SSR <- CI_fis_SSR_upper - CI_fis_SSR_lower
CI_fis95_SSR <- (CI_fis_SSR*0.95)/2



##### CI for 239  ######

j = 1
for (j in 1:100) {
  
  head_239 <- cbind(head, RAD_239)
  
  confidence <- sample(head_239, 239, replace=TRUE)
  interval <- cbind(head, confidence)
  
  stats_239 <- basic.stats(interval)
  write.table(stats_239$overall[9], file="output_CI_fis_239.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(head_239, confidence, interval)
  
}

CI_fst_239_input <- read.table("output_CI_fst_239.csv", sep=",")
CI_fst_239_upper <- max(CI_fst_239_input$V2)
CI_fst_239_lower <- min(CI_fst_239_input$V2)
CI_fst_239 <- CI_fst_239_upper - CI_fst_239_lower
CI_fst95_239 <- (CI_fst_239*0.95)/2

CI_ho_239_input <- read.table("output_CI_ho_239.csv", sep=",")
CI_ho_239_upper <- max(CI_ho_239_input$V2)
CI_ho_239_lower <- min(CI_ho_239_input$V2)
CI_ho_239 <- CI_ho_239_upper - CI_ho_239_lower
CI_ho95_239 <- (CI_ho_239*0.95)/2

CI_hs_239_input <- read.table("output_CI_hs_239.csv", sep=",")
CI_hs_239_upper <- max(CI_hs_239_input$V2)
CI_hs_239_lower <- min(CI_hs_239_input$V2)
CI_hs_239 <- CI_hs_239_upper - CI_hs_239_lower
CI_hs95_239 <- (CI_hs_239*0.95)/2

CI_fis_239_input <- read.table("output_CI_fis_239.csv", sep=",")
CI_fis_239_upper <- max(CI_fis_239_input$V2)
CI_fis_239_lower <- min(CI_fis_239_input$V2)
CI_fis_239 <- CI_fis_239_upper - CI_fis_239_lower
CI_fis95_239 <- (CI_fis_239*0.95)/2



##### CI for 1180  ######

j = 1
for (j in 1:100) {
  
  head_1180 <- cbind(head, RAD_1180)
  
  confidence <- sample(head_1180, 1180, replace=TRUE)
  interval <- cbind(head, confidence)
  
  stats_1180 <- basic.stats(interval)
  write.table(stats_1180$overall[9], file="output_CI_fis_1180.csv", sep=",", append=TRUE, col.names = FALSE)
  
  rm(head_1180, confidence, interval)
  
}

CI_fst_1180_input <- read.table("output_CI_fst_1180.csv", sep=",")
CI_fst_1180_upper <- max(CI_fst_1180_input$V2)
CI_fst_1180_lower <- min(CI_fst_1180_input$V2)
CI_fst_1180 <- CI_fst_1180_upper - CI_fst_1180_lower
CI_fst95_1180 <- (CI_fst_1180*0.95)/2

CI_ho_1180_input <- read.table("output_CI_ho_1180.csv", sep=",")
CI_ho_1180_upper <- max(CI_ho_1180_input$V2)
CI_ho_1180_lower <- min(CI_ho_1180_input$V2)
CI_ho_1180 <- CI_ho_1180_upper - CI_ho_1180_lower
CI_ho95_1180 <- (CI_ho_1180*0.95)/2

CI_hs_1180_input <- read.table("output_CI_hs_1180.csv", sep=",")
CI_hs_1180_upper <- max(CI_hs_1180_input$V2)
CI_hs_1180_lower <- min(CI_hs_1180_input$V2)
CI_hs_1180 <- CI_hs_1180_upper - CI_hs_1180_lower
CI_hs95_1180 <- (CI_hs_1180*0.95)/2

CI_fis_1180_input <- read.table("output_CI_fis_1180.csv", sep=",")
CI_fis_1180_upper <- max(CI_fis_1180_input$V2)
CI_fis_1180_lower <- min(CI_fis_1180_input$V2)
CI_fis_1180 <- CI_fis_1180_upper - CI_fis_1180_lower
CI_fis95_1180 <- (CI_fis_1180*0.95)/2



##### CI for 2317  ######

j = 1
for (j in 1:100) {
  
  head_2317 <- cbind(head, RAD_2317)
  
  confidence <- sample(head_2317, 2317, replace=TRUE)
  interval <- cbind(head, confidence)
  
  stats_2317 <- basic.stats(interval)
   write.table(stats_2317$overall[9], file="output_CI_fis_2317.csv", sep=",", append=TRUE, col.names = FALSE)
  
  rm(head_2317, confidence, interval)
  
}

CI_fst_2317_input <- read.table("output_CI_fst_2317.csv", sep=",")
CI_fst_2317_upper <- max(CI_fst_2317_input$V2)
CI_fst_2317_lower <- min(CI_fst_2317_input$V2)
CI_fst_2317 <- CI_fst_2317_upper - CI_fst_2317_lower
CI_fst95_2317 <- (CI_fst_2317*0.95)/2

CI_ho_2317_input <- read.table("output_CI_ho_2317.csv", sep=",")
CI_ho_2317_upper <- max(CI_ho_2317_input$V2)
CI_ho_2317_lower <- min(CI_ho_2317_input$V2)
CI_ho_2317 <- CI_ho_2317_upper - CI_ho_2317_lower
CI_ho95_2317 <- (CI_ho_2317*0.95)/2

CI_hs_2317_input <- read.table("output_CI_hs_2317.csv", sep=",")
CI_hs_2317_upper <- max(CI_hs_2317_input$V2)
CI_hs_2317_lower <- min(CI_hs_2317_input$V2)
CI_hs_2317 <- CI_hs_2317_upper - CI_hs_2317_lower
CI_hs95_2317 <- (CI_hs_2317*0.95)/2

CI_fis_2317_input <- read.table("output_CI_fis_2317.csv", sep=",")
CI_fis_2317_upper <- max(CI_fis_2317_input$V2)
CI_fis_2317_lower <- min(CI_fis_2317_input$V2)
CI_fis_2317 <- CI_fis_2317_upper - CI_fis_2317_lower
CI_fis95_2317 <- (CI_fis_2317*0.95)/2




##

##### CI for 3831  ######

j = 1
for (j in 1:100) {
  
  head_3831 <- cbind(head, RAD_3831)
  
  confidence <- sample(head_3831, 3831, replace=TRUE)
  interval <- cbind(head, confidence)
  
  stats_3831 <- basic.stats(interval)
  write.table(stats_3831$overall[9], file="output_CI_fis_3831.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(head_3831, confidence, interval)
  
}

CI_fst_3831_input <- read.table("output_CI_fst_3831.csv", sep=",")
CI_fst_3831_upper <- max(CI_fst_3831_input$V2)
CI_fst_3831_lower <- min(CI_fst_3831_input$V2)
CI_fst_3831 <- CI_fst_3831_upper - CI_fst_3831_lower
CI_fst95_3831 <- (CI_fst_3831*0.95)/2

CI_ho_3831_input <- read.table("output_CI_ho_3831.csv", sep=",")
CI_ho_3831_upper <- max(CI_ho_3831_input$V2)
CI_ho_3831_lower <- min(CI_ho_3831_input$V2)
CI_ho_3831 <- CI_ho_3831_upper - CI_ho_3831_lower
CI_ho95_3831 <- (CI_ho_3831*0.95)/2

CI_hs_3831_input <- read.table("output_CI_hs_3831.csv", sep=",")
CI_hs_3831_upper <- max(CI_hs_3831_input$V2)
CI_hs_3831_lower <- min(CI_hs_3831_input$V2)
CI_hs_3831 <- CI_hs_3831_upper - CI_hs_3831_lower
CI_hs95_3831 <- (CI_hs_3831*0.95)/2

CI_fis_3831_input <- read.table("output_CI_fis_3831.csv", sep=",")
CI_fis_3831_upper <- max(CI_fis_3831_input$V2)
CI_fis_3831_lower <- min(CI_fis_3831_input$V2)
CI_fis_3831 <- CI_fis_3831_upper - CI_fis_3831_lower
CI_fis95_3831 <- (CI_fis_3831*0.95)/2

##

##### CI for 6255  ######

j = 1
for (j in 1:100) {
  
  head_6255 <- cbind(head, RAD_6255)
  
  confidence <- sample(head_6255, 6255, replace=TRUE)
  interval <- cbind(head, confidence)
  
  stats_6255 <- basic.stats(interval)
  write.table(stats_6255$overall[9], file="output_CI_fis_6255.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(head_6255, confidence, interval)
  
}

CI_fst_6255_input <- read.table("output_CI_fst_6255.csv", sep=",")
CI_fst_6255_upper <- max(CI_fst_6255_input$V2)
CI_fst_6255_lower <- min(CI_fst_6255_input$V2)
CI_fst_6255 <- CI_fst_6255_upper - CI_fst_6255_lower
CI_fst95_6255 <- (CI_fst_6255*0.95)/2

CI_ho_6255_input <- read.table("output_CI_ho_6255.csv", sep=",")
CI_ho_6255_upper <- max(CI_ho_6255_input$V2)
CI_ho_6255_lower <- min(CI_ho_6255_input$V2)
CI_ho_6255 <- CI_ho_6255_upper - CI_ho_6255_lower
CI_ho95_6255 <- (CI_ho_6255*0.95)/2

CI_hs_6255_input <- read.table("output_CI_hs_6255.csv", sep=",")
CI_hs_6255_upper <- max(CI_hs_6255_input$V2)
CI_hs_6255_lower <- min(CI_hs_6255_input$V2)
CI_hs_6255 <- CI_hs_6255_upper - CI_hs_6255_lower
CI_hs95_6255 <- (CI_hs_6255*0.95)/2

CI_fis_6255_input <- read.table("output_CI_fis_6255.csv", sep=",")
CI_fis_6255_upper <- max(CI_fis_6255_input$V2)
CI_fis_6255_lower <- min(CI_fis_6255_input$V2)
CI_fis_6255 <- CI_fis_6255_upper - CI_fis_6255_lower
CI_fis95_6255 <- (CI_fis_6255*0.95)/2


##


##### CI for 25198  ######

j = 1
for (j in 1:100) {
  
  head_25198 <- cbind(head, RAD_25198)
  
  confidence <- sample(head_25198, 25198, replace=TRUE)
  interval <- cbind(head, confidence)
  
  stats_25198 <- basic.stats(interval)
  #write.table(stats_25198$overall[7], file="output_CI_fst_25198.csv", sep=",", append=TRUE, col.names = FALSE)
  #write.table(stats_25198$overall[1], file="output_CI_ho_25198.csv", sep=",", append=TRUE, col.names = FALSE)
  #write.table(stats_25198$overall[2], file="output_CI_hs_25198.csv", sep=",", append=TRUE, col.names = FALSE)
  write.table(stats_25198$overall[9], file="output_CI_fis_25198.csv", sep=",", append=TRUE, col.names = FALSE)
  rm(head_25198, confidence, interval)
  
}

CI_fst_25198_input <- read.table("output_CI_fst_25198.csv", sep=",")
CI_fst_25198_upper <- max(CI_fst_25198_input$V2)
CI_fst_25198_lower <- min(CI_fst_25198_input$V2)
CI_fst_25198 <- CI_fst_25198_upper - CI_fst_25198_lower
CI_fst95_25198 <- (CI_fst_25198*0.95)/2

CI_ho_25198_input <- read.table("output_CI_ho_25198.csv", sep=",")
CI_ho_25198_upper <- max(CI_ho_25198_input$V2)
CI_ho_25198_lower <- min(CI_ho_25198_input$V2)
CI_ho_25198 <- CI_ho_25198_upper - CI_ho_25198_lower
CI_ho95_25198 <- (CI_ho_25198*0.95)/2

CI_hs_25198_input <- read.table("output_CI_hs_25198.csv", sep=",")
CI_hs_25198_upper <- max(CI_hs_25198_input$V2)
CI_hs_25198_lower <- min(CI_hs_25198_input$V2)
CI_hs_25198 <- CI_hs_25198_upper - CI_hs_25198_lower
CI_hs95_25198 <- (CI_hs_25198*0.95)/2

CI_fis_25198_input <- read.table("output_CI_fis_25198.csv", sep=",")
CI_fis_25198_upper <- max(CI_fis_25198_input$V2)
CI_fis_25198_lower <- min(CI_fis_25198_input$V2)
CI_fis_25198 <- CI_fis_25198_upper - CI_fis_25198_lower
CI_fis95_25198 <- (CI_fis_25198*0.95)/2




 


##### Saving the CIs to file

CIs <- data.frame(CI_fst95_SSR, CI_ho95_SSR, CI_hs95_SSR,
                  CI_fst95_239, CI_ho95_239, CI_hs95_239,
                  CI_fst95_1180, CI_ho95_1180, CI_hs95_1180,
                  CI_fst95_2317, CI_ho95_2317, CI_hs95_2317,
                  CI_fst95_3831, CI_ho95_3831, CI_hs95_3831,
                  CI_fst95_6255, CI_ho95_6255, CI_hs95_6255,
                  CI_fst95_25198, CI_ho95_25198, CI_hs95_25198)

CIs <- data.frame(t(CIs))

write.table(CIs, "CI_April2017.txt", quote = FALSE)

add_CIs <- data.frame(CI_fis_SSR_lower, CI_fis_SSR_upper,
                      CI_fis_239_lower, CI_fis_239_upper,
                      CI_fis_1180_lower, CI_fis_1180_upper,
                      CI_fis_2317_lower, CI_fis_2317_upper,
                      CI_fis_3831_lower, CI_fis_3831_upper,
                      CI_fis_6255_lower, CI_fis_6255_upper,
                      CI_fis_25198_lower, CI_fis_25198_upper)

add_CIs <- data.frame(t(add_CIs))
                      
write.csv(add_CIs, "CI_fis_April2017.csv", quote=FALSE)

#########################################################################
#########################################################################
###  Below are Fis and Ho calculations, by population
#########################################################################
#########################################################################

### Fis

SSR_8_stats$Fis

mean(SSR_8_stats$Fis, na.rm=TRUE)

Fis_SSR <- c(mean(SSR_8_stats$Fis[1:8,1], na.rm=TRUE), mean(SSR_8_stats$Fis[1:8,2], na.rm=TRUE), 
             mean(SSR_8_stats$Fis[1:8,3], na.rm=TRUE), mean(SSR_8_stats$Fis[1:8,4], na.rm=TRUE), 
             mean(SSR_8_stats$Fis[1:8,5], na.rm=TRUE), mean(SSR_8_stats$Fis[1:8,6], na.rm=TRUE),
             mean(SSR_8_stats$Fis[1:8,7], na.rm=TRUE), mean(SSR_8_stats$Fis[1:8,8], na.rm=TRUE), 
             mean(SSR_8_stats$Fis[1:8,9], na.rm=TRUE), mean(SSR_8_stats$Fis[1:8,10], na.rm=TRUE), 
             mean(SSR_8_stats$Fis[1:8,11], na.rm=TRUE), mean(SSR_8_stats$Fis[1:8,12], na.rm=TRUE))

Fis_239 <- c(mean(stats_239$Fis[1:239,1], na.rm=TRUE), mean(stats_239$Fis[1:239,2], na.rm=TRUE), 
             mean(stats_239$Fis[1:239,3], na.rm=TRUE), mean(stats_239$Fis[1:239,4], na.rm=TRUE), 
             mean(stats_239$Fis[1:239,5], na.rm=TRUE), mean(stats_239$Fis[1:239,6], na.rm=TRUE),
             mean(stats_239$Fis[1:239,7], na.rm=TRUE), mean(stats_239$Fis[1:239,8], na.rm=TRUE), 
             mean(stats_239$Fis[1:239,9], na.rm=TRUE), mean(stats_239$Fis[1:239,10], na.rm=TRUE), 
             mean(stats_239$Fis[1:239,11], na.rm=TRUE), mean(stats_239$Fis[1:239,12], na.rm=TRUE))

Fis_1180 <- c(mean(stats_1180$Fis[1:1180,1], na.rm=TRUE), mean(stats_1180$Fis[1:1180,2], na.rm=TRUE), 
             mean(stats_1180$Fis[1:1180,3], na.rm=TRUE), mean(stats_1180$Fis[1:1180,4], na.rm=TRUE), 
             mean(stats_1180$Fis[1:1180,5], na.rm=TRUE), mean(stats_1180$Fis[1:1180,6], na.rm=TRUE),
             mean(stats_1180$Fis[1:1180,7], na.rm=TRUE), mean(stats_1180$Fis[1:1180,8], na.rm=TRUE), 
             mean(stats_1180$Fis[1:1180,9], na.rm=TRUE), mean(stats_1180$Fis[1:1180,10], na.rm=TRUE), 
             mean(stats_1180$Fis[1:1180,11], na.rm=TRUE), mean(stats_1180$Fis[1:1180,12], na.rm=TRUE))

Fis_2317 <- c(mean(stats_2317$Fis[1:2317,1], na.rm=TRUE), mean(stats_2317$Fis[1:2317,2], na.rm=TRUE), 
             mean(stats_2317$Fis[1:2317,3], na.rm=TRUE), mean(stats_2317$Fis[1:2317,4], na.rm=TRUE), 
             mean(stats_2317$Fis[1:2317,5], na.rm=TRUE), mean(stats_2317$Fis[1:2317,6], na.rm=TRUE),
             mean(stats_2317$Fis[1:2317,7], na.rm=TRUE), mean(stats_2317$Fis[1:2317,8], na.rm=TRUE), 
             mean(stats_2317$Fis[1:2317,9], na.rm=TRUE), mean(stats_2317$Fis[1:2317,10], na.rm=TRUE), 
             mean(stats_2317$Fis[1:2317,11], na.rm=TRUE), mean(stats_2317$Fis[1:2317,12], na.rm=TRUE))

Fis_3831 <- c(mean(stats_3831$Fis[1:3831,1], na.rm=TRUE), mean(stats_3831$Fis[1:3831,2], na.rm=TRUE), 
             mean(stats_3831$Fis[1:3831,3], na.rm=TRUE), mean(stats_3831$Fis[1:3831,4], na.rm=TRUE), 
             mean(stats_3831$Fis[1:3831,5], na.rm=TRUE), mean(stats_3831$Fis[1:3831,6], na.rm=TRUE),
             mean(stats_3831$Fis[1:3831,7], na.rm=TRUE), mean(stats_3831$Fis[1:3831,8], na.rm=TRUE), 
             mean(stats_3831$Fis[1:3831,9], na.rm=TRUE), mean(stats_3831$Fis[1:3831,10], na.rm=TRUE), 
             mean(stats_3831$Fis[1:3831,11], na.rm=TRUE), mean(stats_3831$Fis[1:3831,12], na.rm=TRUE))

Fis_6255 <- c(mean(stats_6255$Fis[1:6255,1], na.rm=TRUE), mean(stats_6255$Fis[1:6255,2], na.rm=TRUE), 
              mean(stats_6255$Fis[1:6255,3], na.rm=TRUE), mean(stats_6255$Fis[1:6255,4], na.rm=TRUE), 
              mean(stats_6255$Fis[1:6255,5], na.rm=TRUE), mean(stats_6255$Fis[1:6255,6], na.rm=TRUE),
              mean(stats_6255$Fis[1:6255,7], na.rm=TRUE), mean(stats_6255$Fis[1:6255,8], na.rm=TRUE), 
              mean(stats_6255$Fis[1:6255,9], na.rm=TRUE), mean(stats_6255$Fis[1:6255,10], na.rm=TRUE), 
              mean(stats_6255$Fis[1:6255,11], na.rm=TRUE), mean(stats_6255$Fis[1:6255,12], na.rm=TRUE))

Fis_25198 <- c(mean(stats_25198$Fis[1:25198,1], na.rm=TRUE), mean(stats_25198$Fis[1:25198,2], na.rm=TRUE), 
             mean(stats_25198$Fis[1:25198,3], na.rm=TRUE), mean(stats_25198$Fis[1:25198,4], na.rm=TRUE), 
             mean(stats_25198$Fis[1:25198,5], na.rm=TRUE), mean(stats_25198$Fis[1:25198,6], na.rm=TRUE),
             mean(stats_25198$Fis[1:25198,7], na.rm=TRUE), mean(stats_25198$Fis[1:25198,8], na.rm=TRUE), 
             mean(stats_25198$Fis[1:25198,9], na.rm=TRUE), mean(stats_25198$Fis[1:25198,10], na.rm=TRUE), 
             mean(stats_25198$Fis[1:25198,11], na.rm=TRUE), mean(stats_25198$Fis[1:25198,12], na.rm=TRUE))

Fis_all <- data.frame(Fis_SSR)
Fis_all <- cbind(Fis_all, Fis_239, Fis_1180, Fis_2317, Fis_3831, Fis_6255, Fis_25198)
Fis_all <- t(Fis_all)
Fis_all <- data.frame(Fis_all)


Fis_mean <- c(mean(Fis_all$X1), mean(Fis_all$X2), mean(Fis_all$X3), mean(Fis_all$X4),
              mean(Fis_all$X5), mean(Fis_all$X6), mean(Fis_all$X7), mean(Fis_all$X8),
              mean(Fis_all$X9), mean(Fis_all$X10), mean(Fis_all$X11), mean(Fis_all$X12))

        
Fis_all <- rbind(Fis_all,Fis_mean)
Fis_all <- t(Fis_all)


write.csv(x=Fis_all, file="Fis_April2017.csv")


Ho_SSR <- c(mean(SSR_8_stats$Ho[1:8,1], na.rm=TRUE), mean(SSR_8_stats$Ho[1:8,2], na.rm=TRUE), 
             mean(SSR_8_stats$Ho[1:8,3], na.rm=TRUE), mean(SSR_8_stats$Ho[1:8,4], na.rm=TRUE), 
             mean(SSR_8_stats$Ho[1:8,5], na.rm=TRUE), mean(SSR_8_stats$Ho[1:8,6], na.rm=TRUE),
             mean(SSR_8_stats$Ho[1:8,7], na.rm=TRUE), mean(SSR_8_stats$Ho[1:8,8], na.rm=TRUE), 
             mean(SSR_8_stats$Ho[1:8,9], na.rm=TRUE), mean(SSR_8_stats$Ho[1:8,10], na.rm=TRUE), 
             mean(SSR_8_stats$Ho[1:8,11], na.rm=TRUE), mean(SSR_8_stats$Ho[1:8,12], na.rm=TRUE))

Ho_239 <- c(mean(stats_239$Ho[1:239,1], na.rm=TRUE), mean(stats_239$Ho[1:239,2], na.rm=TRUE), 
             mean(stats_239$Ho[1:239,3], na.rm=TRUE), mean(stats_239$Ho[1:239,4], na.rm=TRUE), 
             mean(stats_239$Ho[1:239,5], na.rm=TRUE), mean(stats_239$Ho[1:239,6], na.rm=TRUE),
             mean(stats_239$Ho[1:239,7], na.rm=TRUE), mean(stats_239$Ho[1:239,8], na.rm=TRUE), 
             mean(stats_239$Ho[1:239,9], na.rm=TRUE), mean(stats_239$Ho[1:239,10], na.rm=TRUE), 
             mean(stats_239$Ho[1:239,11], na.rm=TRUE), mean(stats_239$Ho[1:239,12], na.rm=TRUE))

Ho_1180 <- c(mean(stats_1180$Ho[1:1180,1], na.rm=TRUE), mean(stats_1180$Ho[1:1180,2], na.rm=TRUE), 
              mean(stats_1180$Ho[1:1180,3], na.rm=TRUE), mean(stats_1180$Ho[1:1180,4], na.rm=TRUE), 
              mean(stats_1180$Ho[1:1180,5], na.rm=TRUE), mean(stats_1180$Ho[1:1180,6], na.rm=TRUE),
              mean(stats_1180$Ho[1:1180,7], na.rm=TRUE), mean(stats_1180$Ho[1:1180,8], na.rm=TRUE), 
              mean(stats_1180$Ho[1:1180,9], na.rm=TRUE), mean(stats_1180$Ho[1:1180,10], na.rm=TRUE), 
              mean(stats_1180$Ho[1:1180,11], na.rm=TRUE), mean(stats_1180$Ho[1:1180,12], na.rm=TRUE))

Ho_2317 <- c(mean(stats_2317$Ho[1:2317,1], na.rm=TRUE), mean(stats_2317$Ho[1:2317,2], na.rm=TRUE), 
              mean(stats_2317$Ho[1:2317,3], na.rm=TRUE), mean(stats_2317$Ho[1:2317,4], na.rm=TRUE), 
              mean(stats_2317$Ho[1:2317,5], na.rm=TRUE), mean(stats_2317$Ho[1:2317,6], na.rm=TRUE),
              mean(stats_2317$Ho[1:2317,7], na.rm=TRUE), mean(stats_2317$Ho[1:2317,8], na.rm=TRUE), 
              mean(stats_2317$Ho[1:2317,9], na.rm=TRUE), mean(stats_2317$Ho[1:2317,10], na.rm=TRUE), 
              mean(stats_2317$Ho[1:2317,11], na.rm=TRUE), mean(stats_2317$Ho[1:2317,12], na.rm=TRUE))

Ho_3831 <- c(mean(stats_3831$Ho[1:3831,1], na.rm=TRUE), mean(stats_3831$Ho[1:3831,2], na.rm=TRUE), 
               mean(stats_3831$Ho[1:3831,3], na.rm=TRUE), mean(stats_3831$Ho[1:3831,4], na.rm=TRUE), 
               mean(stats_3831$Ho[1:3831,5], na.rm=TRUE), mean(stats_3831$Ho[1:3831,6], na.rm=TRUE),
               mean(stats_3831$Ho[1:3831,7], na.rm=TRUE), mean(stats_3831$Ho[1:3831,8], na.rm=TRUE), 
               mean(stats_3831$Ho[1:3831,9], na.rm=TRUE), mean(stats_3831$Ho[1:3831,10], na.rm=TRUE), 
               mean(stats_3831$Ho[1:3831,11], na.rm=TRUE), mean(stats_3831$Ho[1:3831,12], na.rm=TRUE))

Ho_6255 <- c(mean(stats_6255$Ho[1:6255,1], na.rm=TRUE), mean(stats_6255$Ho[1:6255,2], na.rm=TRUE), 
             mean(stats_6255$Ho[1:6255,3], na.rm=TRUE), mean(stats_6255$Ho[1:6255,4], na.rm=TRUE), 
             mean(stats_6255$Ho[1:6255,5], na.rm=TRUE), mean(stats_6255$Ho[1:6255,6], na.rm=TRUE),
             mean(stats_6255$Ho[1:6255,7], na.rm=TRUE), mean(stats_6255$Ho[1:6255,8], na.rm=TRUE), 
             mean(stats_6255$Ho[1:6255,9], na.rm=TRUE), mean(stats_6255$Ho[1:6255,10], na.rm=TRUE), 
             mean(stats_6255$Ho[1:6255,11], na.rm=TRUE), mean(stats_6255$Ho[1:6255,12], na.rm=TRUE))

Ho_25198 <- c(mean(stats_25198$Ho[1:25198,1], na.rm=TRUE), mean(stats_25198$Ho[1:25198,2], na.rm=TRUE), 
               mean(stats_25198$Ho[1:25198,3], na.rm=TRUE), mean(stats_25198$Ho[1:25198,4], na.rm=TRUE), 
               mean(stats_25198$Ho[1:25198,5], na.rm=TRUE), mean(stats_25198$Ho[1:25198,6], na.rm=TRUE),
               mean(stats_25198$Ho[1:25198,7], na.rm=TRUE), mean(stats_25198$Ho[1:25198,8], na.rm=TRUE), 
               mean(stats_25198$Ho[1:25198,9], na.rm=TRUE), mean(stats_25198$Ho[1:25198,10], na.rm=TRUE), 
               mean(stats_25198$Ho[1:25198,11], na.rm=TRUE), mean(stats_25198$Ho[1:25198,12], na.rm=TRUE))

Ho_all <- data.frame(Ho_SSR)
Ho_all <- cbind(Ho_all, Ho_239, Ho_1180, Ho_2317, Ho_3831, Ho_6255, Ho_25198)
Ho_all <- t(Ho_all)
Ho_all <- data.frame(Ho_all)


Ho_mean <- c(mean(Ho_all$X1), mean(Ho_all$X2), mean(Ho_all$X3), mean(Ho_all$X4),
              mean(Ho_all$X5), mean(Ho_all$X6), mean(Ho_all$X7), mean(Ho_all$X8),
              mean(Ho_all$X9), mean(Ho_all$X10), mean(Ho_all$X11), mean(Ho_all$X12))


Ho_all <- rbind(Ho_all,Ho_mean)
Ho_all <- t(Ho_all)


write.csv(x=Ho_all, file="Ho_April2017.csv")



