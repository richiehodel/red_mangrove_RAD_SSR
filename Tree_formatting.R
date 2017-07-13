####################################################################################
### Script by Richie Hodel #########################################################
####################################################################################

# This script takes in trees generate by SVDQuartets, annontates and 
# reformats them, and then outputs a tree for each dataset

library("ape")
source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
library("Biostrings")
library("ggplot2")
library("ggtree")

setwd("/Users/richiehodel/Documents/mangroves/RAD-Seq_vs_msat/red_mangrove_RAD_SSR/SVD_trees_April2017/")


####


MyTree <- read.nexus(file="25198.tre")

gulf <- c(MyTree$tip.label[57:80])
atlantic <- c(MyTree$tip.label[1:56], MyTree$tip.label[81:96])
groups <- list(Gulf=c(gulf[1:24]), Atlantic=c(atlantic[1:72]))

MyTree <- groupOTU(MyTree, groups)

ggtree(MyTree, aes(color=group), branch.length = "none") + geom_tiplab(size=2.5, color="black") +
  scale_color_manual(values=c("blue","darkorange1")) + ggtitle("25198 Loci") + 
  theme(plot.title = element_text(hjust=0.5)) +
  geom_label2(aes(label=branch.length, subset=branch.length>20 & !isTip) , color="black", size=2)

ggsave("25198.eps", width=12, height=12, dpi = 300)

####

MyTree <- read.nexus(file="2317.tre")
MyTree <- groupOTU(MyTree, groups)

ggtree(MyTree, aes(color=group), branch.length = "none") + geom_tiplab(size=2.5, color="black") +
 scale_color_manual(values=c("blue","darkorange1")) + ggtitle("2317 Loci") + 
  theme(plot.title = element_text(hjust=0.5)) +
  geom_label2(aes(label=branch.length, subset=branch.length>20 & !isTip) , color="black", size=2)

ggsave("2317.eps", width=12, height=12, dpi = 300)

####

MyTree <- read.nexus(file="3831.tre")
MyTree <- groupOTU(MyTree, groups)

ggtree(MyTree, aes(color=group), branch.length = "none") + geom_tiplab(size=2.5, color="black") +
  scale_color_manual(values=c("blue","darkorange1")) + ggtitle("3831 Loci") + 
  theme(plot.title = element_text(hjust=0.5)) +
  geom_label2(aes(label=branch.length, subset=branch.length>20 & !isTip) , color="black", size=2)

ggsave("3831.eps", width=12, height=12, dpi = 300)


####

MyTree <- read.nexus(file="6255.tre")
MyTree <- groupOTU(MyTree, groups)

ggtree(MyTree, aes(color=group), branch.length = "none") + geom_tiplab(size=2.5, color="black") +
  scale_color_manual(values=c("blue","darkorange1")) + ggtitle("6255 Loci") + 
  theme(plot.title = element_text(hjust=0.5)) +
  geom_label2(aes(label=branch.length, subset=branch.length>70 & !isTip) , color="black", size=2)

ggsave("6255.eps", width=12, height=12, dpi = 300)


####

MyTree <- read.nexus(file="1180.tre")
MyTree <- groupOTU(MyTree, groups)

ggtree(MyTree, aes(color=group), branch.length = "none") + geom_tiplab(size=2.5, color="black") +
  scale_color_manual(values=c("blue","darkorange1")) + ggtitle("1180 Loci") + 
  theme(plot.title = element_text(hjust=0.5)) +
  geom_label2(aes(label=branch.length, subset=branch.length>70 & !isTip) , color="black", size=2)

ggsave("1180.eps", width=12, height=12, dpi = 300)


####

MyTree <- read.nexus(file="548.tre")
MyTree <- groupOTU(MyTree, groups)

ggtree(MyTree, aes(color=group), branch.length = "none") + geom_tiplab(size=2.5, color="black") +
  scale_color_manual(values=c("blue","darkorange1")) + ggtitle("548 Loci") + 
  theme(plot.title = element_text(hjust=0.5)) +
  geom_label2(aes(label=branch.length, subset=branch.length>70 & !isTip) , color="black", size=2)

ggsave("548.eps", width=12, height=12, dpi = 300)


####

MyTree <- read.nexus(file="239.tre")
MyTree <- groupOTU(MyTree, groups)

ggtree(MyTree, aes(color=group), branch.length = "none") + geom_tiplab(size=2.5, color="black") +
  scale_color_manual(values=c("blue","darkorange1")) + ggtitle("239 Loci") + 
  theme(plot.title = element_text(hjust=0.5)) +
  geom_label2(aes(label=branch.length, subset=branch.length>70 & !isTip) , color="black", size=2)

ggsave("239.eps", width=12, height=12, dpi = 300)













