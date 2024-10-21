# if (!require("igraph"))
#   install.packages("igraph")
# if (!require("BiocManager"))
#   install.packages("BiocManager")
# if (!require("phyloseq"))
#   BiocManager::install("phyloseq")
# if (!require("devtools"))
#   install.packages("devtools")
# if (!require("ggClusterNet"))
#   devtools::install_github("taowenmicro/ggClusterNet")
# if (!require("ggplot2"))
#   install.packages("ggplot2")
# if (!require("sna"))
#   install.packages("sna")
# if (!require("tidyfst"))
#   install.packages("tidyfst")

library(igraph)
library(ggplot2)
library(phyloseq)
library(sna)
library(ggClusterNet)

setwd("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXE12_Mental_Microbiome/")
base_theme =
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.minor = element_blank())

phynotype.df <- read.csv("./Data/Psychometric Data.csv", header = T)
sample.df <- read.csv("./Data/SAMPLE_with_entero.csv", header=T)
Genus.df <- read.csv("./Data/Genus.table.csv", header = T)
Phylum.df <- read.csv("./Data/Phylum.table.csv", header = T)

phynotype.df$id <- as.character(phynotype.df$id)
phynotype.df$sampleID <- paste(phynotype.df$Time, phynotype.df$id, sep="_")
sample.df$id <- as.character(sample.df$id)
sample.df$sampleID <- paste(sample.df$Time, sample.df$id, sep = "_")
meta.data <- merge(sample.df, phynotype.df, by="sampleID")
meta.data$count <- "1"

subject.freq <- data.frame(table(meta.data$id.x))

genus.table <- data.frame(t(Genus.df))
colnames(genus.table) <- genus.table[1,]
genus.table <- filter(genus.table, rownames(genus.table) != "X")
genus.meta <- merge(meta.data, genus.table,by.x="Row.names", by.y="row.names")

ASV_table <- read.csv("./Data/ASV_Table.csv", header = T, row.names = 1)
row.names(ASV_table) <- paste("X", rownames(ASV_table), sep="")
TAX2 <- read.csv("./Data/TAX_Table.csv", header = T, row.names = 1)
Sample <- genus.meta
row.names(Sample) <- genus.meta$Row.names

#write.csv(file = "./Result/Result.table/Sample.table.csv",Sample)
Sample <- read.csv(file = "./Result/Result.table/Sample.table.csv", header = T,row.names = 1)

OTU <- otu_table(ASV_table, taxa_are_rows = F)
TAX <- tax_table(as.matrix(TAX2))
physeq = phyloseq(OTU, TAX,sample_data(Sample))
physeq

result = corMicro(ps = physeq, 
                  N = 30, # 根据相关系数选取top100进行可视化
                  method.scale = "TMM", # TMM标准化
                  r.threshold = 0.2, # 相关系数阀值
                  p.threshold = 0.05, # p value阀值
                  method = "pearson")  #这个method需要https://github.com/zdk123/SpiecEasi

# 提取相关矩阵
cor = result[[1]]

igraph <- graph_from_adjacency_matrix(cor, diag = F, mode="undirected",weighted=TRUE)

# 网络中包含的OTU的phyloseq文件提取
ps_net = result[[3]]

# 导出otu表格
otu_table = ps_net %>% 
  vegan_otu() %>%
  t() %>%
  as.data.frame()

clusters(igraph)

#构建分组，可以根据图的最大连接分组，通过clusters(igraph)得到分组信息；也可以自定义分组，这里随机地将100个过滤后的otu分成三组
gp = data.frame(ID = rownames(otu_table), group = sample(1:1, 30, replace = T))

layout = PolygonClusterG(cor = cor, nodeGroup = gp) # 生成网络图布局，'PolygonClusterG'是该论文中的布局
node = layout[[1]] # 提取节点
tax_table = ps_net %>%
  vegan_tax() %>%
  as.data.frame()

# node节点注释
nodes = nodeadd(plotcord  = node, otu_table = otu_table, tax_table = tax_table)
edge = edgeBuild(cor = cor, node = node)  # 构建边

p1 <- ggplot() + geom_segment(data = edge, aes(x = X1, y = Y1, xend = X2, yend = Y2), size = 0.4, color = 'red') +
  geom_point(data = nodes, aes(X1, X2, color = genus), size = 5) + 
  geom_text(data = nodes, aes(X1, X2, label=genus)) +
  #scale_colour_brewer(palette = "Set1") +
  scale_size_continuous(range = c(2, 5)) + 
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p1

ggsave("./plot.png", p1, width = 10, height = 10) # 保存图片


