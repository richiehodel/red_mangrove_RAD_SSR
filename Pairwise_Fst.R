####################################################################################
### Script by Richie Hodel #########################################################
####################################################################################

## This script reads in pairwise FST matrices for each dataset from Genodive and 
## computes average pairwise FST for each dataset and writes a summary table to a file


Pairwise_239 <- read.table(file="GenoDive_pairwise_March2017.gdv", skip=10, nrows = 11,
                           header=TRUE, na.strings = "--")

Pairwise_1180 <- read.table(file="GenoDive_pairwise_March2017.gdv", skip=96, nrows = 11,
                           header=TRUE, na.strings = "--")

Pairwise_2317 <- read.table(file="GenoDive_pairwise_March2017.gdv", skip=139, nrows = 11,
                           header=TRUE, na.strings = "--")

Pairwise_3831 <- read.table(file="GenoDive_pairwise_March2017.gdv", skip=182, nrows = 11,
                           header=TRUE, na.strings = "--")

Pairwise_6255 <- read.table(file="GenoDive_pairwise_March2017.gdv", skip=225, nrows = 11,
                           header=TRUE, na.strings = "--")

Pairwise_25198 <- read.table(file="GenoDive_pairwise_March2017.gdv", skip=268, nrows = 11,
                           header=TRUE, na.strings = "--")

Pairwise_SSR <- read.table(file="GenoDive_SSR_pairwise.gdv", skip=10, nrows = 11,
                           header=TRUE, na.strings = "--")


Avg_239 <- c(mean(na.omit(Pairwise_239$BHKFl)), mean(na.omit(Pairwise_239$CpCFl)), 
              mean(na.omit(Pairwise_239$CvPFl)), mean(na.omit(Pairwise_239$HwdFl)),
              mean(na.omit(Pairwise_239$IsmFl)), mean(na.omit(Pairwise_239$KyLFl)),
              mean(na.omit(Pairwise_239$MlBFl)), mean(na.omit(Pairwise_239$NPRFl)),
              mean(na.omit(Pairwise_239$ShKFl)), mean(na.omit(Pairwise_239$TCBFl)),
              mean(na.omit(Pairwise_239$VKyFl)), mean(na.omit(Pairwise_239$WPBFl)))

Avg_1180 <- c(mean(na.omit(Pairwise_1180$BHKFl)), mean(na.omit(Pairwise_1180$CpCFl)), 
              mean(na.omit(Pairwise_1180$CvPFl)), mean(na.omit(Pairwise_1180$HwdFl)),
              mean(na.omit(Pairwise_1180$IsmFl)), mean(na.omit(Pairwise_1180$KyLFl)),
              mean(na.omit(Pairwise_1180$MlBFl)), mean(na.omit(Pairwise_1180$NPRFl)),
              mean(na.omit(Pairwise_1180$ShKFl)), mean(na.omit(Pairwise_1180$TCBFl)),
              mean(na.omit(Pairwise_1180$VKyFl)), mean(na.omit(Pairwise_1180$WPBFl)))

Avg_2317 <- c(mean(na.omit(Pairwise_2317$BHKFl)), mean(na.omit(Pairwise_2317$CpCFl)), 
              mean(na.omit(Pairwise_2317$CvPFl)), mean(na.omit(Pairwise_2317$HwdFl)),
              mean(na.omit(Pairwise_2317$IsmFl)), mean(na.omit(Pairwise_2317$KyLFl)),
              mean(na.omit(Pairwise_2317$MlBFl)), mean(na.omit(Pairwise_2317$NPRFl)),
              mean(na.omit(Pairwise_2317$ShKFl)), mean(na.omit(Pairwise_2317$TCBFl)),
              mean(na.omit(Pairwise_2317$VKyFl)), mean(na.omit(Pairwise_2317$WPBFl)))

Avg_3831 <- c(mean(na.omit(Pairwise_3831$BHKFl)), mean(na.omit(Pairwise_3831$CpCFl)), 
              mean(na.omit(Pairwise_3831$CvPFl)), mean(na.omit(Pairwise_3831$HwdFl)),
              mean(na.omit(Pairwise_3831$IsmFl)), mean(na.omit(Pairwise_3831$KyLFl)),
              mean(na.omit(Pairwise_3831$MlBFl)), mean(na.omit(Pairwise_3831$NPRFl)),
              mean(na.omit(Pairwise_3831$ShKFl)), mean(na.omit(Pairwise_3831$TCBFl)),
              mean(na.omit(Pairwise_3831$VKyFl)), mean(na.omit(Pairwise_3831$WPBFl)))

Avg_6255 <- c(mean(na.omit(Pairwise_6255$BHKFl)), mean(na.omit(Pairwise_6255$CpCFl)), 
              mean(na.omit(Pairwise_6255$CvPFl)), mean(na.omit(Pairwise_6255$HwdFl)),
              mean(na.omit(Pairwise_6255$IsmFl)), mean(na.omit(Pairwise_6255$KyLFl)),
              mean(na.omit(Pairwise_6255$MlBFl)), mean(na.omit(Pairwise_6255$NPRFl)),
              mean(na.omit(Pairwise_6255$ShKFl)), mean(na.omit(Pairwise_6255$TCBFl)),
              mean(na.omit(Pairwise_6255$VKyFl)), mean(na.omit(Pairwise_6255$WPBFl)))

Avg_25198 <- c(mean(na.omit(Pairwise_25198$BHKFl)), mean(na.omit(Pairwise_25198$CpCFl)), 
              mean(na.omit(Pairwise_25198$CvPFl)), mean(na.omit(Pairwise_25198$HwdFl)),
              mean(na.omit(Pairwise_25198$IsmFl)), mean(na.omit(Pairwise_25198$KyLFl)),
              mean(na.omit(Pairwise_25198$MlBFl)), mean(na.omit(Pairwise_25198$NPRFl)),
              mean(na.omit(Pairwise_25198$ShKFl)), mean(na.omit(Pairwise_25198$TCBFl)),
              mean(na.omit(Pairwise_25198$VKyFl)), mean(na.omit(Pairwise_25198$WPBFl)))
 
Avg_SSR <- c(mean(na.omit(Pairwise_SSR$BHKFl)), mean(na.omit(Pairwise_SSR$CpCFl)), 
                mean(na.omit(Pairwise_SSR$CvPFl)), mean(na.omit(Pairwise_SSR$HwdFl)),
                mean(na.omit(Pairwise_SSR$IsmFl)), mean(na.omit(Pairwise_SSR$KyLFl)),
                mean(na.omit(Pairwise_SSR$MlbFl)), mean(na.omit(Pairwise_SSR$NPRFl)),
                mean(na.omit(Pairwise_SSR$ShKFl)), mean(na.omit(Pairwise_SSR$TCBFl)),
                mean(na.omit(Pairwise_SSR$VKyFl)), mean(na.omit(Pairwise_SSR$WPBFl)))
 

Pairwise_Avg <- data.frame(Avg_239, Avg_1180, Avg_2317, Avg_3831,
                           Avg_6255, Avg_25198, Avg_SSR) 


write.table(Pairwise_Avg, file="Pairwise.txt", row.names=FALSE, quote=FALSE)
 