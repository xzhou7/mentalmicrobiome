#02.data_inte

library(phyloseq)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggstatsplot)
library(vegan)
library(ape)
library(cluster)
library(clusterSim)
library(ade4)

setwd("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXE12_Mental_Microbiome/")

ASV_table <- read.csv("./Data/ASV_Table.csv", header = T, row.names = 1)
TAX2 <- read.csv("./Data/TAX_Table.csv", header = T, row.names = 1)
Sample <- read.csv("./Data/SAMPLE.csv", header = T, row.names = 1)
Sample$id <- as.character(Sample$id)

OTU <- otu_table(ASV_table, taxa_are_rows = F)
TAX <- tax_table(as.matrix(TAX2))
physeq = phyloseq(OTU, TAX,sample_data(Sample))
physeq

physeq_genus <- tax_glom(physeq, taxrank = "genus")
GP.ord <- ordinate(physeq_genus, "PCoA", "bray")
p1 = plot_ordination(physeq_genus, GP.ord, type="sample", color="id", title="taxa")
print(p1)
p1 <- p1 + ggtitle("subjectID")+ coord_equal()
p1


physeq_phylum <- tax_glom(physeq, taxrank = "phylum")
physeq_phylum_freq <-  transform_sample_counts(physeq_phylum, function(x) x / sum(x))
data_phylum <- data.frame(t(otu_table(physeq_phylum_freq)))
data_phylum$phylum <- TAX2$phylum[match(rownames(data_phylum),TAX2$ASV)] 
rownames(data_phylum) <- data_phylum$phylum
data_phylum <- dplyr::select(data_phylum, -phylum)

write.csv(file = "./Data/Phylum.table.csv",data_phylum)
###########################
#calculate enterotype

#01.prepare the data
physeq_genus_freq <-  transform_sample_counts(physeq_genus, function(x) x / sum(x))
data_genus <- data.frame(t(otu_table(physeq_genus_freq)))
data_genus$genus <- TAX2$genus[match(rownames(data_genus),TAX2$ASV)] 
which(data_genus$genus == "Unclassified_Actinobacteria")
data_genus$genus[303] <- "Unclassified_Actinobacteria.2"
rownames(data_genus) <- data_genus$genus
data_genus <- dplyr::select(data_genus, -genus)

write.csv(file = "./Data/Genus.table.csv",data_genus)

#02.Define the JSD/KSD
JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
KLD <- function(x,y) sum(x * log(x/y))

dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD <- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}

#03.make the calculation
data.dist=dist.JSD(data_genus)
data.cluster=pam.clustering(data.dist, k=3)

nclusters = index.G1(t(data_genus), data.cluster, d = data.dist, centrotypes = "medoids")
nclusters=NULL

for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(data_genus),data.cluster_temp,d = data.dist,
                          centrotypes = "medoids")
  }
}

plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
data.cluster=pam.clustering(data.dist, k=3)

obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
data.denoized=noise.removal(data_genus, percent=0.01)

obs.pca=dudi.pca(data.frame(t(data_genus)), scannf=F, nf=10)
obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1) 

s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F)
s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F, col=c(4,2,3))
s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(4,2,3))

enterotype.result <- data.frame(obs.bet$ls)
enterotype.result$Cluster <- data.cluster

#pcoa, do not use 
obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F)
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(3,2,4))

row.names(Sample) <- paste("X",row.names(Sample), sep="")
Sample2 <- merge(Sample, enterotype.result, by="row.names")
Sample2
rownames(Sample2) <- Sample2$Row.names

write.csv(Sample2, file = "./Data/SAMPLE_with_entero.csv")





