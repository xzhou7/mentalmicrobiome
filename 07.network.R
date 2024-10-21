#07 Building microbiome network
#https://biovcnet.github.io/_pages/NetworkScience_igraphviz

library(SpiecEasi)
library(igraph)
library(dplyr)
library(tibble)

#mac
setwd("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXE12_Mental_Microbiome/")

#windows
setwd(dir = "C:/Users/zhoux/Box/XinZhouFiles/Projects/ZXE12_Mental_Microbiome/")

getwd()

customvermillion<-rgb(213/255,94/255,0/255)
custombluegreen<-rgb(0/255,158/255,115/255)
customblue<-rgb(0/255,114/255,178/255)
customskyblue<-rgb(86/255,180/255,233/255)
customreddishpurple<-rgb(204/255,121/255,167/255)


ASV_table <- read.csv("./Data/ASV_Table.csv", header = T, row.names = 1)
TAX_table <- read.csv("./Data/TAX_Table.csv", header = T)

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

se.gl.amgut <- spiec.easi(aggregated_data, method='glasso', lambda.min.ratio=0.02,
                          nlambda=20, pulsar.params=list(rep.num=100))


adj.mat <- getRefit(se.gl.amgut)
table(as.numeric(adj.mat))

se.cor  <- cov2cor(as.matrix(getOptCov(se.gl.amgut)))
weighted.adj.mat <- se.cor*getRefit(se.gl.amgut)

heatmap(as.matrix(weighted.adj.mat))

grph.unweighted <- adj2igraph(adj.mat)
grph <- adj2igraph(weighted.adj.mat)

plot(grph.unweighted,vertex.size=1,vertex.label=NA)

plot(grph,vertex.size=1, vertex.label=NA)
plot(grph,vertex.size=1,
     vertex.label=NA,
     edge.width=1,
     layout=layout.circle(grph))

V(grph)
E(grph)

V(grph)$name <- colnames(aggregated_data)
V(grph)

V(grph)$size <- (degree(grph) + 0.01) # the +1 is to avoid size zero vertices
V(grph)$color <- "black"
plot(grph,
     vertex.label=NA,
     layout=layout.circle(grph))

E(grph)$color <- custombluegreen
E(grph)$color[E(grph)$weight<0] <- customreddishpurple
E(grph)$width <- abs(E(grph)$weight)*10
plot(grph,
     vertex.label=NA,
     layout=layout.circle(grph))
plot(density((E(grph)$weight)),xlab="Edge Weight",main="")

boxplot(abs(E(grph)$weight)~(E(grph)$weight>0),
        xlab="Positive Interaction?",
        ylab="Strength of Interaction")

E(grph)$width[E(grph)$weight<0] <- E(grph)$width[E(grph)$weight<0]*10
plot(grph,
     vertex.label=NA,
     layout=layout.circle(grph))

#Remove edges with very low weight 
weight_threshold <- 0.01
grph <- delete.edges(grph,which(abs(E(grph)$weight)<weight_threshold))
grph.pos <- delete.edges(grph,which(E(grph)$weight<0))
plot(grph.pos,
     vertex.label=NA)

#Remove unconnected vertices
grph.pos <- delete.vertices(grph.pos,which(degree(grph.pos)<1))
plot(grph.pos,
     vertex.label=NA)

#Cleanup a little
V(grph.pos)$size <- V(grph.pos)$size/3
E(grph.pos)$color <- "gray"
plot(grph.pos,
     vertex.label=NA,
     edge.curved=0.5)

plot(grph.pos,
    vertex.label=NA,
    layout=layout_with_fr(grph.pos))
plot(grph.pos,
     vertex.label=NA,
     layout=layout_with_kk(grph.pos))

graph_components <- components(grph.pos)
graph_components

grph.largest.component <- 
  induced.subgraph(grph.pos,V(grph.pos)[which(graph_components$membership == which.max(graph_components$csize))])
plot(grph.largest.component,vertex.label=NA)

