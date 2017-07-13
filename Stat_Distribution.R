####################################################################################
### Script by Richie Hodel #########################################################
####################################################################################

### This script uses Fst, Fis and Ho for each locus to create histograms that 
### show the relative frequency of observed values of each statistic across loci

### This script creates plots to visualize this frequency distribution


library("ggplot2")

### Fst

largest <- na.omit(stats_25198$perloc$Fst)
largest <- abs(largest)
large <- subset(largest, largest<=1)

bucket_25198 <- data.frame(fst=large, data_set="RAD_25198")

largest <- na.omit(stats_6255$perloc$Fst)
largest <- abs(largest)
large <- subset(largest, largest<=1)

bucket_6255 <- data.frame(fst=large, data_set="RAD_6255")

largest <- na.omit(stats_3831$perloc$Fst)
largest <- abs(largest)
large <- subset(largest, largest<=1)

bucket_3831 <- data.frame(fst=large, data_set="RAD_3831")

largest <- na.omit(stats_2317$perloc$Fst)
largest <- abs(largest)
large <- subset(largest, largest<=1)

bucket_2317 <- data.frame(fst=large, data_set="RAD_2317")

largest <- na.omit(stats_1180$perloc$Fst)
largest <- abs(largest)
large <- subset(largest, largest<=1)

bucket_1180 <- data.frame(fst=large, data_set="RAD_1180")

largest <- na.omit(stats_239$perloc$Fst)
largest <- abs(largest)
large <- subset(largest, largest<=1)

bucket_239 <- data.frame(fst=large, data_set="RAD_239")
rm(bucket)
bucket <- rbind(bucket_25198, bucket_6255, bucket_3831, bucket_2317,
                bucket_1180 , bucket_239)

ggplot(bucket, aes(x=fst, fill=data_set)) + geom_histogram(binwidth = 0.02, colour="black") + 
  labs(x="", y="") + scale_fill_brewer(palette="Set1") + theme(legend.title=element_blank(), 
   panel.border = element_rect(colour="black", fill=NA), panel.grid.major = element_line(
     colour="black", size=0.2))


ggsave("fst.eps", width=12, height=12, dpi = 300)

ggsave("fst.tiff", width=6, height=6, dpi = 300)


###   Fis


largest <- na.omit(stats_25198$perloc$Fis)
bucket_fis_25198 <- data.frame(fis=largest, data_set="RAD_25198")

largest <- na.omit(stats_6255$perloc$Fis)
bucket_fis_6255 <- data.frame(fis=largest, data_set="RAD_6255")

largest <- na.omit(stats_3831$perloc$Fis)
bucket_fis_3831 <- data.frame(fis=largest, data_set="RAD_3831")

largest <- na.omit(stats_2317$perloc$Fis)
bucket_fis_2317 <- data.frame(fis=largest, data_set="RAD_2317")

largest <- na.omit(stats_1180$perloc$Fis)
bucket_fis_1180 <- data.frame(fis=largest, data_set="RAD_1180")

largest <- na.omit(stats_239$perloc$Fis)
bucket_fis_239 <- data.frame(fis=largest, data_set="RAD_239")


bucket <- rbind(bucket_fis_25198, bucket_fis_6255, bucket_fis_3831, bucket_fis_2317,
                bucket_fis_1180,bucket_fis_239)


ggplot(bucket, aes(x=fis, fill=data_set)) + geom_histogram(binwidth = 0.04, colour="black") + 
  labs(x="", y="") + scale_fill_brewer(palette="Set1") + theme(legend.title=element_blank(), 
 panel.border = element_rect(colour="black", fill=NA), panel.grid.major = element_line(
     colour="black", size=0.2))


ggplot(bucket, aes(x=fis, fill=data_set)) + geom_histogram(binwidth = 0.04, colour="black") +
  scale_fill_brewer(palette="Set1") + labs(title="Histogram for Fis", hjust=0.5) +
  labs(x="Fis", y="Count") #+ xlim(c(18,52)) + ylim(c(0,30))

ggsave("fis.eps", width=12, height=12, dpi = 300)

ggsave("fis.tiff", width=6, height=6, dpi = 300)


#### Observed heterozygosity  

largest <- na.omit(stats_25198$perloc$Ho)
bucket_Ho_25198 <- data.frame(Ho=largest, data_set="RAD_25198")

largest <- na.omit(stats_6255$perloc$Ho)
bucket_Ho_6255 <- data.frame(Ho=largest, data_set="RAD_6255")

largest <- na.omit(stats_3831$perloc$Ho)
bucket_Ho_3831 <- data.frame(Ho=largest, data_set="RAD_3831")

largest <- na.omit(stats_2317$perloc$Ho)
bucket_Ho_2317 <- data.frame(Ho=largest, data_set="RAD_2317")

largest <- na.omit(stats_1180$perloc$Ho)
bucket_Ho_1180 <- data.frame(Ho=largest, data_set="RAD_1180")

largest <- na.omit(stats_239$perloc$Ho)
bucket_Ho_239 <- data.frame(Ho=largest, data_set="RAD_239")


bucket <- rbind(bucket_Ho_25198, bucket_Ho_6255, bucket_Ho_3831, bucket_Ho_2317,
                bucket_Ho_1180,bucket_Ho_239)


ggplot(bucket, aes(x=Ho, fill=data_set)) + geom_histogram(binwidth = 0.02, colour="black") + 
  labs(x="", y="") + scale_fill_brewer(palette="Set1") + theme(legend.title=element_blank(), 
   panel.border = element_rect(colour="black", fill=NA), panel.grid.major = element_line(
  colour="black", size=0.2))



ggplot(bucket, aes(x=Ho, fill=data_set)) + geom_histogram(binwidth = 0.02, colour="black") +
  scale_fill_brewer(palette="Set1") + labs(title="Histogram for Ho", hjust=0.5) +
  labs(x="Ho", y="Count") #+ xlim(c(18,52)) + ylim(c(0,30))

ggsave("ho.eps", width=12, height=12, dpi = 300)

ggsave("ho.tiff", width=6, height=6, dpi = 300)


