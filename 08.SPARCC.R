#08.SPARCC

library(SpiecEasi)
library(igraph)
library(tibble)
library(mediation)
library(dplyr)
library(tidyverse)
library(compositions)
library(ggsankey)
library(dplyr)
library(ggplot2)
library(future)

#windows
setwd(dir = "C:/Users/zhoux/Box/XinZhouFiles/Projects/ZXE12_Mental_Microbiome/")

getwd()

ASV_table <- read.csv("./Data/ASV_Table.csv", header = T, row.names = 1)
TAX_table <- read.csv("./Data/TAX_Table.csv", header = T)
Sample <- read.csv("./Data/SAMPLE_with_entero.csv", header = T)

Sample <- Sample %>% mutate(Row.names = substr(Row.names, 2, nchar(Row.names)))

tASV <- t(ASV_table) %>% data.frame() %>%  rownames_to_column("ASV")

data.df <- left_join(select(TAX_table, ASV, genus), tASV, by = c("ASV" = "ASV")) %>% select(-ASV)

# Perform the aggregation
aggregated_data <- data.df %>%
  group_by(genus) %>%
  summarise(across(everything(), sum), .groups = "drop")

# Transpose the numeric columns and keep the genus names as a separate vector
genus_names <- aggregated_data$genus
numeric_data <- as.matrix(aggregated_data[-1])
transposed_data <- t(numeric_data)

# Convert the transposed data to a data frame and set the column names
aggregated_data <- as.data.frame(transposed_data)
colnames(aggregated_data) <- genus_names

# Convert to a matrix if desired
aggregated_data <- as.matrix(aggregated_data)

colSums(aggregated_data)

# Remove the last column
aggregated_data <- aggregated_data[, -ncol(aggregated_data), drop = FALSE]

#make sure all column are matched
rownames(aggregated_data) == Sample$X

#
aggregated_data1 <- aggregated_data[rownames(aggregated_data) %in% filter(Sample, Cluster == "1")$X,]
aggregated_data2 <- aggregated_data[rownames(aggregated_data) %in% filter(Sample, Cluster == "2")$X,]
aggregated_data3 <- aggregated_data[rownames(aggregated_data) %in% filter(Sample, Cluster == "3")$X,]

#se.gl.amgut1 <- spiec.easi(aggregated_data1, method='glasso', lambda.min.ratio=0.02,nlambda=20, pulsar.params=list(rep.num=50))
#save(se.gl.amgut1,file = "./Students results/SPARCC.RData")
load("./Students results/SPARCC.RData")

adj.mat1 <- getRefit(se.gl.amgut1)

adj.mat1


















