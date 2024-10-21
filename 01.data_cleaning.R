#mental health microbiome paper 
library(dplyr)
library(stringr)
library(ggplot2)
library(ggstatsplot)
library(vegan)
library(ape)
library(phyloseq)

setwd("~/Library/CloudStorage/Box-Box/Tibshirani Lab/Microbiome/")

meatdata <- read.csv("./Microbiome/metadata_full.csv", header = T)
depression_data <- read.table("./Microbiome/depression_dada2.txt", header = T, sep = "\t")
psychometric_data <- read.csv("./Microbiome/Psychometric Data.csv", header = T)

#clean microbiome data 
TAX <- data.frame(t(depression_data[1:6,]))
dim(TAX)

colnames(TAX) <- c("Phylum", "Class","Order","Family","Genus","Species")
TAX2 <- TAX[-1,]

HMP_Taxa <- TAX2

HMP_Taxa$Species[is.na(HMP_Taxa$Species)&!is.na(HMP_Taxa$Genus)] <- paste("Unclassified", HMP_Taxa$Genus[is.na(HMP_Taxa$Species)&!is.na(HMP_Taxa$Genus)], sep="_")

HMP_Taxa$Species[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&!is.na(HMP_Taxa$Family)] <- paste("Unclassified", HMP_Taxa$Family[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&!is.na(HMP_Taxa$Family)], sep="_")
HMP_Taxa$Genus[is.na(HMP_Taxa$Genus)&!is.na(HMP_Taxa$Family)] <- paste("Unclassified", HMP_Taxa$Family[is.na(HMP_Taxa$Genus)&!is.na(HMP_Taxa$Family)], sep="_")

HMP_Taxa$Species[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&!is.na(HMP_Taxa$Order)] <- paste("Unclassified", HMP_Taxa$Order[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&!is.na(HMP_Taxa$Order)], sep="_")
HMP_Taxa$Genus[is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&!is.na(HMP_Taxa$Order)] <- paste("Unclassified", HMP_Taxa$Order[is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&!is.na(HMP_Taxa$Order)], sep="_")
HMP_Taxa$Family[is.na(HMP_Taxa$Family)&!is.na(HMP_Taxa$Order)] <- paste("Unclassified", HMP_Taxa$Order[is.na(HMP_Taxa$Family)&!is.na(HMP_Taxa$Order)], sep="_")

HMP_Taxa$Species[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&!is.na(HMP_Taxa$Class)] <- paste("Unclassified", HMP_Taxa$Class[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&!is.na(HMP_Taxa$Class)], sep="_")
HMP_Taxa$Genus[is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&!is.na(HMP_Taxa$Class)] <- paste("Unclassified", HMP_Taxa$Class[is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&!is.na(HMP_Taxa$Class)], sep="_")
HMP_Taxa$Family[is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&!is.na(HMP_Taxa$Class)] <- paste("Unclassified", HMP_Taxa$Class[is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&!is.na(HMP_Taxa$Class)], sep="_")
HMP_Taxa$Order[is.na(HMP_Taxa$Order)&!is.na(HMP_Taxa$Class)] <- paste("Unclassified", HMP_Taxa$Class[is.na(HMP_Taxa$Order)&!is.na(HMP_Taxa$Class)], sep="_")

HMP_Taxa$Species[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)] <- paste("Unclassified", HMP_Taxa$Phylum[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)], sep="_")
HMP_Taxa$Genus[is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)] <- paste("Unclassified", HMP_Taxa$Phylum[is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)], sep="_")
HMP_Taxa$Family[is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)] <- paste("Unclassified", HMP_Taxa$Phylum[is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)], sep="_")
HMP_Taxa$Order[is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)] <- paste("Unclassified", HMP_Taxa$Phylum[is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)], sep="_")
HMP_Taxa$Class[is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)]<- paste("Unclassified", HMP_Taxa$Phylum[is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)], sep="_")

HMP_Taxa$Species[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)] <- paste("Unclassified", HMP_Taxa$Kingdom[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)], sep="_")
HMP_Taxa$Genus[is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)] <- paste("Unclassified", HMP_Taxa$Kingdom[is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)], sep="_")
HMP_Taxa$Family[is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)] <- paste("Unclassified", HMP_Taxa$Kingdom[is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)], sep="_")
HMP_Taxa$Order[is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)] <- paste("Unclassified", HMP_Taxa$Kingdom[is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)], sep="_")
HMP_Taxa$Class[is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)] <- paste("Unclassified", HMP_Taxa$Kingdom[is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)], sep="_")
HMP_Taxa$Phylum[is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)] <- paste("Unclassified", HMP_Taxa$Kingdom[is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)], sep="_")

TAX2 <- HMP_Taxa
rownames(TAX2) <- paste("ASV", 1:5169, sep="_")

TAX2$Genus[which(str_detect(TAX2$Genus, "Prevotella"))] <- "Prevotella"

ASV <- data.frame(depression_data[7:162,])
dim(ASV)
ASV[1:8, 1:5]

#identical rownames for replacement
colnames(ASV) <- paste("ASV", 0:5169, sep="_")
identical(row.names(TAX2[5120:5169,]), colnames(ASV)[-1][5120:5169])

#table(duplicated(SAMPLE$kitid))
meatdata[1:5,1:5]

ASV$kitid <- str_extract(ASV$ASV_0, "[:digit:]*")

(ASV %>% filter(kitid == "559299082"))[1:5]

#####need to change this once the metadata is sorted out
rownames(ASV) <- ASV$ASV_0
ASV <- select(ASV, -ASV_0,-kitid)
dim(ASV)

ASV_table <- sapply(ASV, as.numeric, simplify = T)
row.names(ASV_table) <- rownames(ASV) 
ASV_table <- as.data.frame(ASV_table)

seq.depth <- data.frame(rowSums(ASV_table))
colnames(seq.depth) <- "depth"

arrange(seq.depth,depth)

gghistostats(
  data       = seq.depth,
  x          = depth,
  title      = "Sequencing Depth",
  binwidth   = 3000
)

######
ASV_BCdist <- vegdist(ASV_table, method = "bray")
ASV_BCdist_PCOA <- pcoa(ASV_BCdist)
ASV_PCOA_AXIS <- ASV_BCdist_PCOA$vectors %>% data.frame()

ggplot(ASV_PCOA_AXIS, aes(x=Axis.1, y=Axis.2)) + geom_point()

####
colnames(TAX2) <- c("phylum", "class", "order", "family", "genus", "species")
TAX2$ASV <- rownames(TAX2)

###this is a duplicate
(ASV_table[row.names(ASV_table)=="526299004_NA0021494121_5_gut",1:10])
(ASV_table[row.names(ASV_table)=="526299004_NA0021476928_5_gut",1:10])

#remove 4 samples that are duplicate or do not have metadata
ASV_table_clean <-filter(ASV_table, !row.names(ASV_table) %in% c("297298359_NA0021476929_6_spare","480323674_NA0021458739_5_gut","526299004_NA0021494121_5_gut", "559299082_NA0021502516_6_spare"))
dim(ASV_table)

write.csv(ASV_table_clean, file = "../../XinZhouFiles/Projects/ZXE12_Mental_Microbiome/Data/ASV_Table.csv")
write.csv(TAX2, file = "../../XinZhouFiles/Projects/ZXE12_Mental_Microbiome/Data/TAX_Table.csv")

