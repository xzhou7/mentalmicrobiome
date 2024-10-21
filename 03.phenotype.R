#03.phenotype analysis

library(phyloseq)
#library(ggstatsplot)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(vegan)
library(ape)
library(cluster)
library(clusterSim)
library(ade4)
library(patchwork)
library(networkD3)
library(tidyverse)
library(vegan)
library(reshape2)
library(ggpubr)
library(HMP)
library(ComplexHeatmap)
library(factoextra)
library(GLMMadaptive)
library(ggpubr)
library(dunn.test)
library(phylosmith)
library(ggraph)
library(sjPlot)

#stop sign
library(x)

#Mac
#setwd("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXE12_Mental_Microbiome/")

#Windows
setwd("C:/Users/zhoux/Box/XinZhouFiles/Projects/ZXE12_Mental_Microbiome/")

source("./Analysis Code/00.tools.R")

phynotype.df <- read.csv("./Data/Psychometric Data.csv", header = T)
sample.df <- read.csv("./Data/SAMPLE_with_entero.csv", header=T)
Genus.df <- read.csv("./Data/Genus.table.csv", header = T)
Phylum.df <- read.csv("./Data/Phylum.table.csv", header = T)
Cytokine <- read.csv("./Data/corrected_cytokines_data.csv", header = T, row.names = 1)

phynotype.df$id <- as.character(phynotype.df$id)
phynotype.df$sampleID <- paste(phynotype.df$Time, phynotype.df$id, sep="_")
sample.df$id <- as.character(sample.df$id)
sample.df$sampleID <- paste(sample.df$Time, sample.df$id, sep = "_")

####
length(intersect(phynotype.df$kitid, sample.df$kitid))
meta.data <- merge(sample.df, phynotype.df, by="sampleID")
meta.data$count <- "1"

table(meta.data$Time.x)

table(meta.data$id.x, meta.data$depressed)
table(meta.data$depressed)

subject.freq <- data.frame(table(meta.data$id.x))

# filter(subject.freq, Freq == 4)$Var1
# 
# p1 <- ggplot(filter(meta.data, Time.x == "T1") %>% filter(id.x %in% filter(subject.freq, Freq == 4)$Var1), aes(x=depressed,y=count, fill=as.character(Cluster))) + geom_col(position = 'fill') 
# p2 <- ggplot(filter(meta.data, Time.x == "T3") %>% filter(id.x %in% filter(subject.freq, Freq == 4)$Var1), aes(x=depressed,y=count, fill=as.character(Cluster))) + geom_col(position = 'fill')
# p3 <- ggplot(filter(meta.data, Time.x == "T4") %>% filter(id.x %in% filter(subject.freq, Freq == 4)$Var1), aes(x=depressed,y=count, fill=as.character(Cluster))) + geom_col(position = 'fill')
# p4 <- ggplot(filter(meta.data, Time.x == "T5") %>% filter(id.x %in% filter(subject.freq, Freq == 4)$Var1), aes(x=depressed,y=count, fill=as.character(Cluster))) + geom_col(position = 'fill')
# 
# p1 + p2 +p3 + p4 + plot_layout(guides = "collect")
# 

p1.1 <- ggplot(filter(meta.data, Time.x == "T1"), aes(x=depressed,y=count, fill=as.character(Cluster))) + geom_col(position = 'fill') 
p2.1 <- ggplot(filter(meta.data, Time.x == "T3"), aes(x=depressed,y=count, fill=as.character(Cluster))) + geom_col(position = 'fill')
p3.1 <- ggplot(filter(meta.data, Time.x == "T4"), aes(x=depressed,y=count, fill=as.character(Cluster))) + geom_col(position = 'fill')
p4.1 <- ggplot(filter(meta.data, Time.x == "T5"), aes(x=depressed,y=count, fill=as.character(Cluster))) + geom_col(position = 'fill')

p1.1 + p2.1 +p3.1 + p4.1 + plot_layout(guides = "collect")

p1.2 <- ggplot(filter(meta.data, depressed == "Depressed"), aes(x=Time.x,y=count, fill=as.character(Cluster))) + geom_col(position = 'fill') + scale_fill_manual(values = entero_color_1) + ggtitle("Depressed")
p2.2 <- ggplot(filter(meta.data, depressed != "Depressed"), aes(x=Time.x,y=count, fill=as.character(Cluster))) + geom_col(position = 'fill') + scale_fill_manual(values = entero_color_1) + ggtitle("Non-Depressed")
p.entero <- p1.2 + p2.2 + plot_layout(guides = "collect")
p.entero
#ggsave("./Result/Enter.by.T.pdf", p.entero, height = 5, width = 5, dpi = 300)


#revision
#find delta BDI score
library(dplyr)
library(tidyr)

# Assuming 'meta.data' is your data frame
meta.data.delta <- meta.data %>%
  dplyr::select(id.x, Time.x, bdi_total, depressed) %>%
  pivot_wider(names_from = Time.x, values_from = bdi_total) %>%
  mutate(delta_bdi_total = coalesce(T3, T4) - T1)

# View the resulting data frame
print(meta.data.delta)
table(meta.data.delta$depressed)

p.r1 <- ggplot(meta.data.delta, aes(x= depressed, y=-delta_bdi_total, fill= depressed)) + geom_jitter(aes(color=depressed)) + 
  geom_boxplot(alpha=0.3, outliers = F) + theme_minimal() + scale_color_manual(values = depress_color) + scale_fill_manual(values = depress_color)
p.r1

#ggsave("./Result/revision/Reviewer1.2.pdf",p.r1, width = 5, height = 5, dpi = 300)

#get effect size
meta.data.delta_clean <- meta.data.delta %>%
  filter(!is.na(delta_bdi_total), !is.na(depressed))
meta.data.delta_clean$depressed <- factor(meta.data.delta_clean$depressed, levels = c("Not Depressed", "Depressed"))
t_test_result <- t.test(delta_bdi_total ~ depressed, data = meta.data.delta_clean)
print(t_test_result)

shapiro_test_dep <- shapiro.test(meta.data.delta_clean$delta_bdi_total[meta.data.delta_clean$depressed == "Depressed"])
shapiro_test_not_dep <- shapiro.test(meta.data.delta_clean$delta_bdi_total[meta.data.delta_clean$depressed == "Not Depressed"])

print(shapiro_test_dep)
print(shapiro_test_not_dep)

library(rcompanion)
wilcox_test_result <- wilcox.test(delta_bdi_total ~ depressed, data = meta.data.delta_clean, exact = FALSE)
print(wilcox_test_result)
rank_biserial_result <- wilcoxonR(x = meta.data.delta_clean$delta_bdi_total, g = meta.data.delta_clean$depressed)
print(rank_biserial_result)

library(effsize)
cliff_delta_result <- cliff.delta(delta_bdi_total ~ depressed, data = meta.data.delta_clean)
print(cliff_delta_result)

###########################ChiSqure Test
#select those who only have T1 and T3
meta.data.T1T3 <- filter(filter(meta.data, Time.x %in% c("T1","T3")))

meta.data.D <- filter(meta.data, depressed == "Depressed")
freq_De <- table( meta.data.D$Cluster,meta.data.D$Time.x) %>% data.frame() %>% dcast(Var1~Var2, value.var = "Freq") %>%  data.frame()
freq_De

meta.data.N <- filter(meta.data, depressed != "Depressed")
freq_ND <- table(meta.data.N$Cluster,meta.data.N$Time.x) %>% data.frame() %>% dcast(Var1~Var2, value.var = "Freq") %>%  data.frame()
freq_ND

chisq.test(cbind(freq_De$T1, freq_De$T3))
chisq.test(cbind(freq_De$T1, freq_De$T4))
chisq.test(cbind(freq_De$T1, freq_De$T5))

chisq.test(cbind(freq_De$T1, freq_ND$T1))
chisq.test(cbind(freq_De$T3, freq_ND$T3))
chisq.test(cbind(freq_De$T4, freq_ND$T4))
chisq.test(cbind(freq_De$T5, freq_ND$T5))

#no significant difference between the distribution of enterotypes, nor the treatment cause any change

###########################
genus.table <- data.frame(t(Genus.df))
colnames(genus.table) <- genus.table[1,]
genus.table <- filter(genus.table, rownames(genus.table) != "X")
genus.meta <- merge(meta.data, genus.table,by.x="Row.names", by.y="row.names")

phylum.table <- data.frame(t(Phylum.df))
colnames(phylum.table) <- phylum.table[1,]
phylum.table <- filter(phylum.table, rownames(phylum.table) != "X")
phylum.meta <- merge(meta.data, phylum.table,by.x="Row.names", by.y="row.names")


###revision 1 Sample Overlapping
phynotype.df
genus.meta
Cytokine.df

colnames(phynotype.df)
colnames(genus.meta)
colnames(Cytokine.df)

phynotype.df2 <- phynotype.df %>%
  rename(ID = id, Time = Time) %>%  
  dplyr::select(ID, Time) %>%
  mutate(Pheno = TRUE) %>%
  distinct()

genus.meta2 <- genus.meta %>%
  rename(ID = id.x, Time = Time.x) %>%
  dplyr::select(ID, Time) %>%
  mutate(Microbiome = TRUE) %>%
  distinct()

Cytokine.df2 <- Cytokine.df %>%
  rename(ID = id, Time = Time) %>%
  dplyr::select(ID, Time) %>%
  mutate(Cytokine = TRUE) %>%
  distinct()

# Convert ID and Time columns to character in phynotype.df
phynotype.df2 <- phynotype.df2 %>%
  mutate(
    ID = as.character(ID),
    Time = as.character(Time)
  )

# Convert ID and Time columns to character in genus.meta2
genus.meta2 <- genus.meta2 %>%
  mutate(
    ID = as.character(ID),
    Time = as.character(Time)
  )

# Convert ID and Time columns to character in Cytokine.df2
Cytokine.df2 <- Cytokine.df2 %>%
  mutate(
    ID = as.character(ID),
    Time = as.character(Time)
  )

# Merge data frames and specify suffixes to handle overlapping columns
combined_data_1C <- full_join(phynotype.df2, genus.meta2, by = c("ID", "Time"), suffix = c(".pheno", ".micro"))
combined_data_1C <- full_join(combined_data_1C, Cytokine.df2, by = c("ID", "Time"), suffix = c("", ".cytokine"))

# Replace NA with FALSE in indicator columns
combined_data_1C <- combined_data_1C %>%
  mutate(across(c(Pheno, Microbiome, Cytokine), ~ ifelse(is.na(.), FALSE, .)))

upset_data <- combined_data_1C %>%
  dplyr::select(ID,Time, Pheno, Microbiome, Cytokine) %>%
  distinct() %>%
  mutate(across(c(Pheno, Microbiome, Cytokine), as.integer))

upset_matrix <- as.data.frame(upset_data)
rownames(upset_matrix) <- upset_matrix$Sample
upset_matrix$Sample <- NULL

library(UpSetR)
upset(upset_matrix, nsets = 3, nintersects = NA, order.by = "freq")

upset_matrix <- upset_matrix %>% mutate(ID_Time = paste(ID, Time, sep = "_"))

long_data <- upset_matrix %>%
  pivot_longer(cols = c(Pheno, Microbiome, Cytokine),
               names_to = "Data_Type",
               values_to = "Available")

upset_data <- upset_matrix %>%
  mutate(Sample = ID_Time) %>%
  dplyr::select(Sample, Pheno, Microbiome, Cytokine) %>%
  distinct()

upset_data_df <- as.data.frame(upset_data)
rownames(upset_data_df) <- upset_data_df$Sample
upset_data_df$Sample <- NULL

upset(upset_data_df, nsets = 3, order.by = "freq", 
      main.bar.color = "steelblue", sets.bar.color = "steelblue")

upset_matrix <- upset_matrix %>% mutate(Total = Pheno + Microbiome + Cytokine)

upset_matrix <- filter(upset_matrix,Time %in% c("T1", "T3", "T4", "T5"))
upset_matrix_renamed <- upset_matrix %>%
  mutate(Time = recode(Time,
                       "T1" = "T1",  # Keep T1 unchanged
                       "T3" = "T2",  # Rename T3 to T2
                       "T4" = "T3",  # Rename T4 to T3
                       "T5" = "T4"))  # Rename T5 to T4
listsingle <- c(17, 45, 47,49, 60, 63, 68, 70, 71,72)
upset_matrix_renamed <- filter(upset_matrix_renamed, ! ID %in% listsingle)

length(unique(upset_matrix_renamed$ID))

#write.csv(file = "./Result/revision/SampleOverlap.csv",upset_matrix_renamed)

ptimepoint <- ggplot(upset_matrix_renamed, aes(x = ID, y = Time, fill = factor(Total))) +
  geom_tile(color = "white") +  # Creates the tiles with white borders
  scale_fill_manual(
    values = c("1" = "grey", "2" = "blue", "3" = "darkblue"),
    name = "Total Data Types",
    labels = c("1" = "1", "2" = "2", "3" = "3")
  ) +
  theme_minimal() +
  labs(
    title = "Data Availability by Subject and Time Point",
    x = "Subject ID",
    y = "Sample Time Point"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),  # Rotate x-axis labels for readability
    axis.text.y = element_text(size = 10),
    legend.position = "bottom"
  )
ptimepoint
#ggsave(filename = "./Result/revision/SampleOverlapping_A.pdf", ptimepoint, width = 10, height = 4, dpi = 300)

ptimepoint2 <- upset(upset_matrix_renamed, nsets = 3, order.by = "freq", 
      main.bar.color = "steelblue", sets.bar.color = "steelblue") 

# 
# pdf(file = "./Result/revision/SampleOverlapping_B.pdf",width = 4, height = 3)
# print(ptimepoint2)
# dev.off()


#####
genus.table.nm <- sapply(genus.table, as.numeric, simplify=T)
(colSums(genus.table.nm) %>% sort(decreasing=T))[1:20]

p.bac <- ggplot(genus.meta, aes(x=Cluster, y=as.numeric(Bacteroides))) + ylab("Bacteroides") +
  geom_jitter() + geom_boxplot(aes(group=Cluster),alpha=0.3, outlier.alpha = 0) + base_theme
p.bac

p.prevo <- ggplot(genus.meta, aes(x=Cluster, y=as.numeric(Prevotella))) + ylab("Prevotella") +
  geom_jitter() + geom_boxplot(aes(group=Cluster),alpha=0.3,outlier.alpha = 0) + base_theme
p.prevo

p.Ru <- ggplot(genus.meta, aes(x=Cluster, y=(as.numeric(Unclassified_Ruminococcaceae) + as.numeric(Ruminococcus_2)+as.numeric(Ruminococcus_1)))) + ylab("Ruminococcus\nUnclassified_Ruminococcaceae") +
  geom_jitter() + geom_boxplot(aes(group=Cluster),alpha=0.3,outlier.alpha = 0) + base_theme + ylim(0,0.1)
p.Ru

p.fir.tx <- ggplot(phylum.meta, aes(x=as.character(Cluster),y=as.numeric(Firmicutes))) + 
  geom_jitter() + geom_boxplot(aes(group=Cluster),alpha=0.3, outlier.alpha = 0) + base_theme
p.fir.tx

#ggsave(filename = "./Result/Bac.Enter.pdf", p.bac, width = 5, height = 4, dpi = 300)
#ggsave(filename = "./Result/Pre.Enter.pdf", p.prevo, width = 5, height = 4, dpi = 300)
#ggsave(filename = "./Result/Ru.Enter.pdf", p.Ru, width = 5, height = 4, dpi = 300)

p_depre <- genus.meta %>% ggplot(aes(x = Time.x, y= bdi_total, group=id.x, color=depressed)) + geom_point() + geom_line(aes(color=depressed), linetype= "dashed")
p_depre <- p_depre + scale_color_manual(values = depress_color) + base_theme + scale_x_discrete(labels = c("T1", "T2", "T3", "T4")) + ylab("Beck Depression Inventory (BDI)") +xlab(NULL) 
p_depre
#ggsave(filename = "./Result/BDI_Total_H.pdf", p_depre, width = 6, height = 3, dpi=300)

p.proteo <- ggplot(phylum.meta, aes(x=Time.y,y=as.numeric(Proteobacteria))) + geom_boxplot()
p.proteo <- p.proteo + facet_grid(depressed~.) + stat_compare_means(ref.group = "T1") + base_theme
p.proteo
#ggsave(filename = "./Result/Proteobacteria.bytime.pdf", p.proteo, width = 4, height = 3, dpi=300)

proteo_data <- phylum.meta %>%
  mutate(Proteobacteria = as.numeric(Proteobacteria)) %>% 
  group_by(Time.y, depressed) %>%
  summarise(
    mean_proteo = mean(Proteobacteria, na.rm = TRUE),
    sem_proteo = sd(Proteobacteria, na.rm = TRUE) / sqrt(n())) 

p.proteo2 <- proteo_data %>% 
  ggplot(aes(x = Time.y, y = mean_proteo, group = depressed)) +
  geom_line(aes(color = depressed), linewidth = 1.2) +
  geom_point(aes(color = depressed)) +
  geom_errorbar(aes(ymin = mean_proteo - sem_proteo, ymax = mean_proteo + sem_proteo,color = depressed), width = 0.2) +
  labs(y = "Proteobacteria", x = "Time") +
  theme_minimal()
p.proteo2
#ggsave(filename = "./Result/Proteobacteria.bytime_2.pdf", p.proteo2, width = 5, height = 2, dpi=300)

p.bac <- ggplot(phylum.meta, aes(x=Time.y,y=as.numeric(Bacteroidetes))) + geom_boxplot()
p.bac <- p.bac + facet_grid(depressed~.)  + stat_compare_means(ref.group = "T1") + base_theme
p.bac
#ggsave(filename = "./Result/Bacteroidetes.bytime.pdf", p.bac, width = 4, height = 3, dpi=300)

bac_data <- phylum.meta %>%
  mutate(Bacteroidetes = as.numeric(Bacteroidetes)) %>% 
  group_by(Time.y, depressed) %>%
  summarise(
    mean_bac = mean(Bacteroidetes, na.rm = TRUE),
    sem_bac = sd(Bacteroidetes, na.rm = TRUE) / sqrt(n())) 

p.bac2 <- bac_data %>%
  ggplot(aes(x = Time.y, y = mean_bac, group = depressed)) +
  geom_line(aes(color = depressed), linewidth = 1.2) +
  geom_point(aes(color = depressed)) +
  geom_errorbar(aes(ymin = mean_bac - sem_bac, ymax = mean_bac + sem_bac, color = depressed), width = 0.2) +
  labs(y = "Bacteroidetes", x = "Time") +
  theme_minimal()
p.bac2
#ggsave(filename = "./Result/Bacteroidetes.bytime_2.pdf", p.bac2, width = 5, height = 2, dpi=300)

p.act <- ggplot(phylum.meta, aes(x=Time.y,y=as.numeric(Actinobacteria))) + geom_boxplot()
p.act <- p.act + facet_grid(depressed~.) + stat_compare_means(ref.group = "T1") + base_theme
p.act
#ggsave(filename = "./Result/Actinobacteria.bytime.pdf", p.act, width = 4, height = 3, dpi=300)

actino_data <- phylum.meta %>%
  mutate(Actinobacteria = as.numeric(Actinobacteria)) %>% 
  group_by(Time.y, depressed) %>%
  summarise(
    mean_actino = mean(Actinobacteria, na.rm = TRUE),
    sem_actino = sd(Actinobacteria, na.rm = TRUE) / sqrt(n()))

# Creating the plot p.actino
p.actino2 <- actino_data %>%
  ggplot(aes(x = Time.y, y = mean_actino, group = depressed)) +
  geom_line(aes(color = depressed), linewidth = 1.2) +
  geom_point(aes(color = depressed)) +
  geom_errorbar(aes(ymin = mean_actino - sem_actino, ymax = mean_actino + sem_actino, color = depressed), width = 0.2) +
  labs(y = "Actinobacteria", x = "Time") +
  theme_minimal()
p.actino2

# Saving the plot as PDF
#ggsave(filename = "./Result/Actinobacteria.bytime_2.pdf", p.actino2, width = 5, height = 2, dpi = 300)

p.fir <- ggplot(phylum.meta, aes(x=Time.y,y=as.numeric(Firmicutes))) + geom_boxplot()
p.fir <- p.fir + facet_grid(depressed~.)  + stat_compare_means(ref.group = "T1") + base_theme
p.fir
#ggsave(filename = "./Result/Firmicutes.bytime.pdf", p.fir, width = 4, height = 3, dpi=300)

firmicutes_data <- phylum.meta %>%
  mutate(Firmicutes = as.numeric(Firmicutes)) %>% 
  group_by(Time.y, depressed) %>%
  summarise(
    mean_firmicutes = mean(Firmicutes, na.rm = TRUE),
    sem_firmicutes = sd(Firmicutes, na.rm = TRUE) / sqrt(n()))

# Creating the plot p.firmicutes
p.firmicutes <- firmicutes_data %>%
  ggplot(aes(x = Time.y, y = mean_firmicutes, group = depressed)) +
  geom_line(aes(color = depressed), linewidth = 1.2) +
  geom_point(aes(color = depressed)) +
  geom_errorbar(aes(ymin = mean_firmicutes - sem_firmicutes, ymax = mean_firmicutes + sem_firmicutes, color = depressed), width = 0.2) +
  labs(y = "Firmicutes", x = "Time") +
  theme_minimal()

# Displaying the plot
p.firmicutes

# Saving the plot as PDF
#ggsave(filename = "./Result/Firmicutes.bytime_2.pdf", p.firmicutes, width = 5, height = 2, dpi = 300)


phylum.meta <- phylum.meta %>% mutate(pther_phylum = 1-as.numeric(Firmicutes)-as.numeric(Actinobacteria)-as.numeric(Bacteroidetes)-as.numeric(Proteobacteria))

p.other <- ggplot(phylum.meta, aes(x=Time.y,y=as.numeric(pther_phylum))) + geom_boxplot()
p.other <- p.other + facet_grid(depressed~.)  + stat_compare_means(ref.group = "T1") + base_theme
p.other

other_phylum_data <- phylum.meta %>%
  mutate(pther_phylum = as.numeric(pther_phylum)) %>% 
  group_by(Time.y, depressed) %>%
  summarise(
    mean_other_phylum = mean(pther_phylum, na.rm = TRUE),
    sem_other_phylum = sd(pther_phylum, na.rm = TRUE) / sqrt(n()))

# Creating the plot p.other_phylum
p.other_phylum <- other_phylum_data %>%
  ggplot(aes(x = Time.y, y = mean_other_phylum, group = depressed)) +
  geom_line(aes(color = depressed), linewidth = 1.2) +
  geom_point(aes(color = depressed)) +
  geom_errorbar(aes(ymin = mean_other_phylum - sem_other_phylum, ymax = mean_other_phylum + sem_other_phylum, color = depressed), width = 0.2) +
  labs(y = "Other Phylum", x = "Time") +
  theme_minimal()

p.other_phylum

#ggsave(filename = "./Result/Other_Phylum.bytime_2.pdf", p.other_phylum, width = 5, height = 2, dpi = 300)

phylum.meta <- phylum.meta %>% mutate(BF_Ratio = as.numeric(Bacteroidetes)/as.numeric(Firmicutes))

p.bf <- ggplot(phylum.meta, aes(x=Time.y,y=as.numeric(BF_Ratio))) + geom_boxplot()
p.bf <- p.bf + facet_grid(depressed~.)  + stat_compare_means(ref.group = "T1") + base_theme
p.bf

BF_Ratio_data <- phylum.meta %>%
  group_by(Time.y, depressed) %>%
  summarise(
    mean_BF_Ratio = mean(BF_Ratio, na.rm = TRUE),
    sem_BF_Ratio = sd(BF_Ratio, na.rm = TRUE) / sqrt(n()))

# Creating the plot p.BF_Ratio
p.BF_Ratio <- BF_Ratio_data %>%
  ggplot(aes(x = Time.y, y = mean_BF_Ratio, group = depressed)) +
  geom_line(aes(color = depressed), linewidth = 1.2) +
  geom_point(aes(color = depressed)) +
  geom_errorbar(aes(ymin = mean_BF_Ratio - sem_BF_Ratio, ymax = mean_BF_Ratio + sem_BF_Ratio, color = depressed), width = 0.2) +
  labs(y = "BF Ratio (Bacteroidetes/Firmicutes)", x = "Time") +
  theme_minimal()

p.BF_Ratio
#ggsave(filename = "./Result/BF_Ratio.bytime_2.pdf", p.BF_Ratio, width = 5, height = 2, dpi = 300)

phylum.meta$Time.y <- as.factor(phylum.meta$Time.y)
phylum.meta$depressed <- as.factor(phylum.meta$depressed)

bac_anova <- aov(Bacteroidetes ~ Time.y * depressed, data = phylum.meta)
summary(bac_anova)

firmicutes_anova <- aov(Firmicutes ~ Time.y * depressed, data = phylum.meta)
summary(firmicutes_anova)

actinobacteria_anova <- aov(Actinobacteria ~ Time.y * depressed, data = phylum.meta)
summary(actinobacteria_anova)

proteobacteria_anova <- aov(Proteobacteria ~ Time.y * depressed, data = phylum.meta)
summary(proteobacteria_anova)

other_anova <- aov(pther_phylum ~ Time.y * depressed, data = phylum.meta)
summary(other_anova)

BF_Ratio_anova <- aov(BF_Ratio ~ Time.y * depressed, data = phylum.meta)
summary(proteobacteria_anova)

#Plot genus with Time 
#ggplot(genus.meta, aes(x=Time.y, y=as.numeric(Prevotella))) + geom_boxplot() +stat_compare_means(ref.group = "T1")+ facet_grid(depressed~.) + base_theme 

###################################
#revision
# Assuming 'meta.data' is your data frame
genus.meta.collection <- genus.meta %>%
  dplyr::select(id.x, Time.x, bdi_total, depressed, Cluster) %>%
  pivot_wider(names_from = Time.x, values_from = bdi_total) %>%
  mutate(delta_bdi_total = T3 - T1)

genus.clean.table <- dplyr::select(genus.meta,id.x, Bacteroides:Acetobacter) %>% mutate(across(-id.x, as.numeric))
genus.clean.table.mean <- aggregate(genus.clean.table, FUN = mean, by = list(genus.clean.table$id.x)) %>% dplyr::select(-id.x)
meta.data.delta

#keep the top 20 bacteria and combined with delta BDI
correla.1.df <- genus.clean.table.mean %>% dplyr::select(Group.1, 
                                  names(colSums(genus.clean.table.mean[-1]) %>% sort(decreasing = T))[1:20]) %>% 
  left_join(meta.data.delta, by = c( "Group.1"="id.x" ))

#Find correlations between delta BDI and bacteria genus
bacteria_cols <- names(correla.1.df)[which(names(correla.1.df) == "Bacteroides"):
    which(names(correla.1.df) == "Lachnoclostridium")]
correla.1.df[bacteria_cols] <- lapply(correla.1.df[bacteria_cols], as.numeric)
correla.1.df$delta_bdi_total <- as.numeric(correla.1.df$delta_bdi_total)
correla.1.df$depressed <- factor(correla.1.df$depressed, levels = c("Not Depressed", "Depressed"))

correlation_results <- data.frame(
  Bacteria = character(),
  Spearman_rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (bacteria in bacteria_cols) {
  cor_test <- cor.test(
    correla.1.df[[bacteria]],
    correla.1.df$delta_bdi_total,
    method = "spearman",
    exact = FALSE)  
  correlation_results <- rbind(
      correlation_results,
      data.frame(
        Bacteria = bacteria,
        Spearman_rho = cor_test$estimate,
        p_value = cor_test$p.value
      )
    )
}

correlation_results$adj_p_value <- p.adjust(correlation_results$p_value, method = "BH")
print(correlation_results)

pcore.Subdoligranulum <- correla.1.df %>% ggplot(aes(x= delta_bdi_total, y= Subdoligranulum)) + geom_point(aes(color=depressed)) + geom_smooth(method = "lm")
pcore.Subdoligranulum <- pcore.Subdoligranulum + theme_minimal() + scale_color_manual(values = depress_color)
pcore.Subdoligranulum

#ggsave(filename = "./Result/revision/Reviewer1.3.1.pdf", pcore.Subdoligranulum, width = 6, height = 5, dpi = 300)

pcore.Subdoligranulum.depre <- correla.1.df %>% ggplot(aes(x= depressed, y= Subdoligranulum)) + geom_boxplot()
pcore.Subdoligranulum.depre <- pcore.Subdoligranulum.depre + stat_compare_means(method = "t.test")
pcore.Subdoligranulum.depre

p.adjust(0.003, n=20)

#######
ASV_table <- read.csv("./Data/ASV_Table.csv", header = T, row.names = 1)
row.names(ASV_table) <- paste("X", rownames(ASV_table), sep="")
TAX2 <- read.csv("./Data/TAX_Table.csv", header = T, row.names = 1)
Sample <- genus.meta
row.names(Sample) <- genus.meta$Row.names

#write.csv(file = "./Result/Result.table/Sample.table.csv",Sample)

OTU <- otu_table(ASV_table, taxa_are_rows = F)
TAX <- tax_table(as.matrix(TAX2))
physeq = phyloseq(OTU, TAX,sample_data(Sample))
physeq

physeq_ord <- ordinate(physeq, "PCoA", "bray")
p1 = plot_ordination(physeq, physeq_ord, type="sample", color="depressed", title="taxa")
print(p1)
p1 <- p1 + stat_ellipse(aes(color=depressed, group=depressed), level=0.95)
p1 <- p1 + facet_wrap(.~Time.y) +ggtitle("Depression by time")+ coord_equal() + base_theme
p1

#clean taxa 1.do not have phylum info, 2. do not show three counts in at least 20% samples
taxa_to_remove <- which(is.na(tax_table(physeq)[,"phylum"]))
physeq_clean <-  prune_taxa(taxa_names(physeq)[-taxa_to_remove], physeq)
physeq_clean2 = filter_taxa(physeq_clean, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

physeq_ord_clean2 <- ordinate(physeq_clean2, "PCoA", "bray")
p2 = plot_ordination(physeq_clean2, physeq_ord_clean2, type="sample", color="depressed", title="taxa")
print(p2)
p2 <- p2 + stat_ellipse(aes(color=depressed, group=depressed), level=0.95)
p2 <- p2 + facet_wrap(.~Time.y) +ggtitle("Depression by time (major taxa)")+ coord_equal() + base_theme
p2

#get diversity estimate
set.seed(777)
physeq_rarefy <- rarefy_even_depth(physeq, sample.size = 10000)
richness_df <- estimate_richness(physeq_rarefy)
richness_meta <- merge(richness_df, meta.data, by.x="row.names", by.y = "Row.names")

############total 
# https://github.com/schuyler-smith/phylosmith/
filtered_obj <- conglomerate_taxa(physeq_clean2, "family")

#network_layout_ps(filtered_obj, treatment = NULL, subset = NULL, co_occurrence_table = table_co, algorithm = 'kk')

#calculate table of correlation
table_co <- co_occurrence(filtered_obj, treatment =  "depressed", rho = 0, p = 0.05, cores = 0)

#filter table
filtered_table_co <- table_co %>% filter((rho > 0.49) | (rho < -0.42))
filtered_table_co$p.adjust <- p.adjust(filtered_table_co$p, method = "BH", n=length(table_co$p))
max(filtered_table_co$p.adjust)
#plot depressed net
p.co <- co_occurrence_network(filtered_obj, co_occurrence_table = filter(filtered_table_co, Treatment == "Depressed"),
                      classification = 'family') + geom_node_text(aes(label = name)) + ggtitle("Depressed")
p.co

#plot Un-depressed net
p.co2 <- co_occurrence_network(filtered_obj, co_occurrence_table = filter(filtered_table_co, Treatment != "Depressed"),
                              classification = 'family') + geom_node_text(aes(label = name)) + ggtitle("Not Depressed")
p.co2

#combine plot
p.cor.network <- p.co + p.co2 + plot_layout(guides = "collect")
p.cor.network

#ggsave(filename = "./Result/Network_dep_none.pdf", p.cor.network, width = 14, height = 7, dpi = 300)

#create rho null distribution
permute_rho_data <- permute_rho(filtered_obj, treatment = NULL, replicate_samples = 'id.x', permutations = 1000, cores = 0)

#10000 permutation gave same value, ran on Aug.18
#permute_rho_data <- permute_rho(filtered_obj, treatment = NULL, replicate_samples = 'id.x', permutations = 10000, cores = 0)
p.hist_roh <- histogram_permuted_rhos(permute_rho_data, p = 0.0001, x_breaks = 0.05) 
p.hist_roh <- p.hist_roh + scale_fill_manual(values = c("lightgrey"))
p.hist_roh
#create a rho table
permute_rho_data <- permute_rho_data[order(permute_rho_data$rho), ]
permute_rho_data$cumulative_count <- cumsum(permute_rho_data$Count)

threshold <- 0.0001 * max(permute_rho_data$cumulative_count)
positive_rho_threshold <- min(permute_rho_data$rho[permute_rho_data$cumulative_count > (max(permute_rho_data$cumulative_count) - threshold)])
negative_rho_threshold <- max(permute_rho_data$rho[permute_rho_data$cumulative_count < threshold])
positive_rho_threshold
negative_rho_threshold
# Print the thresholds
print(paste("Positive rho threshold:", positive_rho_threshold))
print(paste("Negative rho threshold:", negative_rho_threshold))

#diversity_compare <- list(c("T1","T3"),c("T1","T4"),c("T1","T5"))
#ggplot(filter(richness_meta), aes(x=Time.x, y=Chao1)) + geom_boxplot() + stat_compare_means(ref.group = "T1") + facet_wrap(depressed~.)
p.ace <- ggplot(filter(richness_meta), aes(x=Time.x, y=ACE)) + geom_jitter()  + facet_wrap(depressed~.) + geom_smooth()
# ggplot(filter(richness_meta), aes(x=Time.x, y=Observed)) + geom_boxplot() + stat_compare_means(ref.group = "T5") + facet_wrap(depressed~.)
# ggplot(filter(richness_meta), aes(x=Time.x, y=InvSimpson)) + geom_boxplot() + stat_compare_means(ref.group = "T5") + facet_wrap(depressed~.)
# ggplot(filter(richness_meta), aes(x=Time.x, y=Shannon)) + geom_boxplot() + stat_compare_means(ref.group = "T5") + facet_wrap(depressed~.)

p.ace <- p.ace + base_theme + ylab("ACE")
p.ace

#ggsave(filename = "./Result/Diversity_Time.pdf", p.ace, width = 5, height = 4, dpi = 300)

p.ace <- ggplot(filter(richness_meta), aes(x=depressed, y=Shannon)) + geom_boxplot() + stat_compare_means(label = "p.format",ref.group = "Depressed", method = "t.test") + facet_grid(Time.x~Cluster)
p.ace <- p.ace + base_theme + ylab("Shannon")
p.ace

richness_meta

#Observed
Observed_data <- richness_meta %>%
  group_by(Time.y, depressed) %>%
  summarise(
    mean_Observed = mean(Observed, na.rm = TRUE),
    sem_Observed = sd(Observed, na.rm = TRUE) / sqrt(n()))

p.Observed <- Observed_data %>%
  ggplot(aes(x = Time.y, y = mean_Observed, group = depressed)) +
  geom_line(aes(color = depressed), linewidth = 1.2) +
  geom_point(aes(color = depressed)) +
  geom_errorbar(aes(ymin = mean_Observed - sem_Observed, ymax = mean_Observed + sem_Observed, color = depressed), width = 0.2) +
  labs(y = "Observed Richness", x = "Time") +
  theme_minimal()
p.Observed

#ACE
ACE_data <- richness_meta %>%
  group_by(Time.y, depressed) %>%
  summarise(
    mean_ACE = mean(ACE, na.rm = TRUE),
    sem_ACE = sd(ACE, na.rm = TRUE) / sqrt(n()))

p.ACE <- ACE_data %>%
  ggplot(aes(x = Time.y, y = mean_ACE, group = depressed)) +
  geom_line(aes(color = depressed), linewidth = 1.2) +
  geom_point(aes(color = depressed)) +
  geom_errorbar(aes(ymin = mean_ACE - sem_ACE, ymax = mean_ACE + sem_ACE, color = depressed), width = 0.2) +
  labs(y = "ACE", x = "Time") +
  theme_minimal()
p.ACE

##Chao1
Chao1_data <- richness_meta %>%
  group_by(Time.y, depressed) %>%
  summarise(
    mean_Chao1 = mean(Chao1, na.rm = TRUE),
    sem_Chao1 = sd(Chao1, na.rm = TRUE) / sqrt(n()))

p.Chao1 <- Chao1_data %>%
  ggplot(aes(x = Time.y, y = mean_Chao1, group = depressed)) +
  geom_line(aes(color = depressed), linewidth = 1.2) +
  geom_point(aes(color = depressed)) +
  geom_errorbar(aes(ymin = mean_Chao1 - sem_Chao1, ymax = mean_Chao1 + sem_Chao1, color = depressed), width = 0.2) +
  labs(y = "Chao1", x = "Time") +
  theme_minimal()

p.Chao1

#ggsave(filename = "./Result/CHAO1_by_depression.pdf", p.Chao1, width = 6, height = 5, dpi = 300)

Chao1_anova1 <- aov(Chao1 ~ Time.y * depressed, data = richness_meta)
summary(Chao1_anova1)

#shannon
Shannon_data <- richness_meta %>%
  group_by(Time.y, depressed) %>%
  summarise(
    mean_Shannon = mean(Shannon, na.rm = TRUE),
    sem_Shannon = sd(Shannon, na.rm = TRUE) / sqrt(n()))

p.Shannon <- Shannon_data %>%
  ggplot(aes(x = Time.y, y = mean_Shannon, group = depressed)) +
  geom_line(aes(color = depressed), linewidth = 1.2) +
  geom_point(aes(color = depressed)) +
  geom_errorbar(aes(ymin = mean_Shannon - sem_Shannon, ymax = mean_Shannon + sem_Shannon, color = depressed), width = 0.2) +
  labs(y = "Shannon Diversity", x = "Time") +
  theme_minimal()
p.Shannon

#sim
InvSimpson_data <- richness_meta %>%
  group_by(Time.y, depressed) %>%
  summarise(
    mean_InvSimpson = mean(InvSimpson, na.rm = TRUE),
    sem_InvSimpson = sd(InvSimpson, na.rm = TRUE) / sqrt(n()))

p.InvSimpson <- InvSimpson_data %>%
  ggplot(aes(x = Time.y, y = mean_InvSimpson, group = depressed)) +
  geom_line(aes(color = depressed), linewidth = 1.2) +
  geom_point(aes(color = depressed)) +
  geom_errorbar(aes(ymin = mean_InvSimpson - sem_InvSimpson, ymax = mean_InvSimpson + sem_InvSimpson, color = depressed), width = 0.2) +
  labs(y = "Inverse Simpson Index", x = "Time") +
  theme_minimal()
p.InvSimpson

richness_meta$Time.x <- factor(richness_meta$Time.x, levels = c("T1", "T3", "T4", "T5"))
richness_meta$Cluster <- factor(richness_meta$Cluster, levels =c("3","2","1"))
richness_meta$days <- 0
richness_meta$days[richness_meta$Time.x == "T2"] <- 9
richness_meta$days[richness_meta$Time.x == "T3"] <- 30
richness_meta$days[richness_meta$Time.x == "T4"] <- 90

#analysis on ACE richness and its relationship to other parameters
lm.rt <- lmerTest::lmer(ACE ~ 1 + Time.y + (1|id.x), data=filter(richness_meta, depressed == "Depressed"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(ACE ~ 1 + Time.y  + (1|id.x), data=filter(richness_meta, depressed != "Depressed"))
summary(lm.rt)

lm.rt <- lmerTest::lmer(ACE ~ 1 + Time.y * depressed + (1|id.x), data=filter(richness_meta, Cluster == "1"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(ACE ~ 1 + Time.y * depressed + (1|id.x), data=filter(richness_meta, Cluster == "2"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(ACE ~ 1 + Time.y * depressed + (1|id.x), data=filter(richness_meta, Cluster == "3"))
summary(lm.rt)

lm.rt <- lmerTest::lmer(ACE ~ 1 + Time.y * depressed  + Cluster +  (1|id.x), data=filter(richness_meta))
summary(lm.rt)

lm.rt <- lmerTest::lmer(ACE ~ 1 + days * depressed  : Cluster + (1|id.x), data=filter(richness_meta))
summary(lm.rt)

p.ace_d <- richness_meta %>% filter(depressed == "Depressed") %>%  
  ggline(x = "Time.x", y = "ACE", color = "Cluster", add = c("mean_se", "jitter"),size = 0.7,point.size = 0.8) + 
  scale_color_manual(values = entero_color_1) + ylim(0, 500)

p.ace_n <- richness_meta %>% filter(depressed != "Depressed") %>%  
  ggline(x = "Time.x", y = "ACE", color = "Cluster",linetype = 2, add = c("mean_se", "jitter"),size = 0.7,point.size = 0.8) + 
  scale_color_manual(values = entero_color_1)+ ylim(0, 500)

p.ace_d+ p.ace_n


#analysis on Chao1 richness and its relationship to other parameters
# Linear mixed models for Depressed group
lm.rt <- lmerTest::lmer(Chao1 ~ 1 + Time.y + (1|id.x), data=filter(richness_meta, depressed == "Depressed"))
summary(lm.rt)

# Linear mixed models for Not Depressed group
lm.rt <- lmerTest::lmer(Chao1 ~ 1 + Time.y  + (1|id.x), data=filter(richness_meta, depressed != "Depressed"))
summary(lm.rt)

# Linear mixed models for each Cluster
lm.rt <- lmerTest::lmer(Chao1 ~ 1 + Time.y * depressed + (1|id.x), data=filter(richness_meta, Cluster == "1"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(Chao1 ~ 1 + Time.y * depressed + (1|id.x), data=filter(richness_meta, Cluster == "2"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(Chao1 ~ 1 + Time.y * depressed + (1|id.x), data=filter(richness_meta, Cluster == "3"))
summary(lm.rt)

# Linear mixed model considering Cluster as a fixed effect
lm.rt <- lmerTest::lmer(Chao1 ~ 1 + Time.y * depressed  + Cluster +  (1|id.x), data=filter(richness_meta))
summary(lm.rt)

# Linear mixed model considering days and interaction with Cluster
lm.rt <- lmerTest::lmer(Chao1 ~ 1 + Time.y * depressed  : Cluster + (1|id.x), data=filter(richness_meta))
summary(lm.rt)

# Plots
p.chao1_d <- richness_meta %>% filter(depressed == "Depressed") %>%  
  ggline(x = "Time.x", y = "Chao1", color = "Cluster", add = c("mean_se", "jitter"),size = 0.7,point.size = 0.8) + 
  scale_color_manual(values = entero_color_1) + ylim(0, 500)

p.chao1_n <- richness_meta %>% filter(depressed != "Depressed") %>%  
  ggline(x = "Time.x", y = "Chao1", color = "Cluster",linetype = 2, add = c("mean_se", "jitter"),size = 0.7,point.size = 0.8) + 
  scale_color_manual(values = entero_color_1)+ ylim(0, 500)

p.chao1_d + p.chao1_n

#analysis on Observed richness and its relationship to other parameters

lm.rt1 <- lmerTest::lmer(Observed ~ 1 + days * Cluster  + (1|id.x), data=filter(richness_meta, depressed == "Depressed"))
summary(lm.rt1)
lm.rt2 <- lmerTest::lmer(Observed ~ 1 + days * Cluster + (1|id.x), data=filter(richness_meta, depressed != "Depressed"))
summary(lm.rt2)

lm.rt <- lmerTest::lmer(Observed ~ 1 +  depressed + (1|id.x), data=filter(richness_meta, Cluster == "1"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(Observed ~ 1 +  depressed + (1|id.x), data=filter(richness_meta, Cluster == "2"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(Observed ~ 1 +  depressed + (1|id.x) , data=filter(richness_meta, Cluster == "3"))
summary(lm.rt)

lm.rt <- lmerTest::lmer(Observed ~ 1 + days* depressed : Cluster  + (1|id.x), data=filter(richness_meta))
summary(lm.rt)

p.obs_d <- richness_meta %>% filter(depressed == "Depressed") %>%  
  ggline(x = "Time.x", y = "Observed", color = "Cluster", add = c("mean_se", "jitter"),size = 0.7,point.size = 0.8) + 
  scale_color_manual(values = entero_color_1) + ylim(0, 500)

p.obs_n <- richness_meta %>% filter(depressed != "Depressed") %>%  
  ggline(x = "Time.x", y = "Observed", color = "Cluster",linetype = 2, add = c("mean_se", "jitter"),size = 0.7,point.size = 0.8) + 
  scale_color_manual(values = entero_color_1)+ ylim(0, 500)

p.obs_d+ p.obs_n

#analysis on Simpson Diversity and its relationship to other parameters

lm.rt <- lmerTest::lmer(Simpson ~ 1 + days * Cluster  + (1|id.x), data=filter(richness_meta, depressed == "Depressed"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(Simpson ~ 1 + days * Cluster + (1|id.x), data=filter(richness_meta, depressed != "Depressed"))
summary(lm.rt)

lm.rt <- lmerTest::lmer(Simpson ~ 1 + depressed + (1|id.x), data=filter(richness_meta, Cluster == "1"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(Simpson ~ 1 + depressed + (1|id.x), data=filter(richness_meta, Cluster == "2"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(Simpson ~ 1 +  depressed + (1|id.x) , data=filter(richness_meta, Cluster == "3"))
summary(lm.rt)

lm.rt <- lmerTest::lmer(Simpson ~ 1 + depressed * Cluster  + (1|id.x), data=filter(richness_meta, Cluster != "0"))
summary(lm.rt)

p.smp_d <- richness_meta %>% filter(depressed == "Depressed") %>%  
  ggline(x = "Time.x", y = "Simpson", color = "Cluster", add = c("mean_se", "jitter"),size = 0.7,point.size = 0.8) + 
  scale_color_manual(values = entero_color_1) + ylim(0.9, 1) + ggtitle("Depressed")

p.smp_n <- richness_meta %>% filter(depressed != "Depressed") %>%  
  ggline(x = "Time.x", y = "Simpson", color = "Cluster",linetype = 2, add = c("mean_se", "jitter"),size = 0.7,point.size = 0.8) + 
  scale_color_manual(values = entero_color_1) + ylim(0.9, 1) + ggtitle("None-Depressed")

p.smp_d + p.smp_n

p.smp_d <- richness_meta %>% filter(depressed == "Depressed") %>%  
  ggline(x = "Time.x", y = "Simpson", color = "Cluster", add = c("mean_se"),size = 0.7,point.size = 0.8) + 
  scale_color_manual(values = entero_color_1) + ylim(0.9, 1) + ggtitle("Depressed")

p.smp_n <- richness_meta %>% filter(depressed != "Depressed") %>%  
  ggline(x = "Time.x", y = "Simpson", color = "Cluster",linetype = 2, add = c("mean_se"),size = 0.7,point.size = 0.8) + 
  scale_color_manual(values = entero_color_1) + ylim(0.9, 1) + ggtitle("None-Depressed")

p.smp_d + p.smp_n

#analysis on InvSimpson Diversity and its relationship to other parameters
lm.rt <- lmerTest::lmer(InvSimpson ~ 1 + days * Cluster  + (1|id.x), data=filter(richness_meta, depressed == "Depressed"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(InvSimpson ~ 1 + days * Cluster + (1|id.x), data=filter(richness_meta, depressed != "Depressed"))
summary(lm.rt)

lm.rt <- lmerTest::lmer(InvSimpson ~ 1 + days * depressed + (1|id.x), data=filter(richness_meta, Cluster == "1"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(InvSimpson ~ 1 + days * depressed + (1|id.x), data=filter(richness_meta, Cluster == "2"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(InvSimpson ~ 1 + days * depressed + (1|id.x), data=filter(richness_meta, Cluster == "3"))
summary(lm.rt)

lm.rt <- lmerTest::lmer(InvSimpson ~ 1 + days* depressed + Cluster  + (1|id.x), data=filter(richness_meta))
summary(lm.rt)

p.ismp_d <- richness_meta %>% filter(depressed == "Depressed") %>%  
  ggline(x = "Time.x", y = "InvSimpson", color = "Cluster", add = c("mean_se", "jitter"),size = 0.7,point.size = 0.8) + 
  scale_color_manual(values = entero_color_1) 

p.ismp_n <- richness_meta %>% filter(depressed != "Depressed") %>%  
  ggline(x = "Time.x", y = "InvSimpson", color = "Cluster",linetype = 2, add = c("mean_se", "jitter"),size = 0.7,point.size = 0.8) + 
  scale_color_manual(values = entero_color_1) 

p.ismp_d + p.ismp_n

#analysis on Shannon Diversity and its relationship to other parameters
lm.rt <- lmerTest::lmer(Shannon ~ 1 + days * Cluster  + (1|id.x), data=filter(richness_meta, depressed == "Depressed"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(Shannon ~ 1 + days * Cluster + (1|id.x), data=filter(richness_meta, depressed != "Depressed"))
summary(lm.rt)

lm.rt <- lmerTest::lmer(Shannon ~ 1 + days * depressed + (1|id.x), data=filter(richness_meta, Cluster == "1"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(Shannon ~ 1 + days * depressed + (1|id.x), data=filter(richness_meta, Cluster == "2"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(Shannon ~ 1 + days * depressed + (1|id.x), data=filter(richness_meta, Cluster == "3"))
summary(lm.rt)

lm.rt <- lmerTest::lmer(Shannon ~ 1 + days* depressed + Cluster  + (1|id.x), data=filter(richness_meta))
summary(lm.rt)

p.sha_d <- richness_meta %>% filter(depressed == "Depressed") %>%  
  ggline(x = "Time.x", y = "Shannon", color = "Cluster", add = c("mean_se", "jitter"),size = 0.7,point.size = 0.8) + 
  scale_color_manual(values = entero_color_1) 

p.sha_n <- richness_meta %>% filter(depressed != "Depressed") %>%  
  ggline(x = "Time.x", y = "Shannon", color = "Cluster",linetype = 2, add = c("mean_se", "jitter"),size = 0.7,point.size = 0.8) + 
  scale_color_manual(values = entero_color_1) 

p.sha_d + p.sha_n

#analysis on Fisher Diversity and its relationship to other parameters
lm.rt <- lmerTest::lmer(Fisher ~ 1 + days * Cluster  + (1|id.x), data=filter(richness_meta, depressed == "Depressed"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(Fisher ~ 1 + days * Cluster + (1|id.x), data=filter(richness_meta, depressed != "Depressed"))
summary(lm.rt)

lm.rt <- lmerTest::lmer(Fisher ~ 1 + days * depressed + (1|id.x), data=filter(richness_meta, Cluster == "1"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(Fisher ~ 1 + days * depressed + (1|id.x), data=filter(richness_meta, Cluster == "2"))
summary(lm.rt)
lm.rt <- lmerTest::lmer(Fisher ~ 1 + days * depressed + (1|id.x), data=filter(richness_meta, Cluster == "3"))
summary(lm.rt)

lm.rt <- lmerTest::lmer(Fisher ~ 1 + days* depressed + Cluster  + (1|id.x), data=filter(richness_meta))
summary(lm.rt)

p.fis_d <- richness_meta %>% filter(depressed == "Depressed") %>%  
  ggline(x = "Time.x", y = "Fisher", color = "Cluster", add = c("mean_se", "jitter"),size = 0.7,point.size = 0.8) + 
  scale_color_manual(values = entero_color_1) + ylim(10, 80)

p.fis_n <- richness_meta %>% filter(depressed != "Depressed") %>%  
  ggline(x = "Time.x", y = "Fisher", color = "Cluster",linetype = 2, add = c("mean_se", "jitter"),size = 0.7,point.size = 0.8) + 
  scale_color_manual(values = entero_color_1) + ylim(10, 80)

p.fis_d + p.fis_n

############wil test
Observed_wide <- richness_meta %>% dplyr::select(Observed,Time.x, id.x,depressed) %>% dcast(id.x + depressed ~ Time.x, value.var = "Observed")
wilcox.test(dplyr::filter(Observed_wide, depressed == "Depressed")$T1,dplyr::filter(Observed_wide, depressed == "Depressed")$T4, paired = T)

Chao_wide <- richness_meta %>% dplyr::select(Chao1,Time.x, id.x,depressed) %>% dcast(id.x + depressed ~ Time.x, value.var = "Chao1")
wilcox.test(dplyr::filter(Chao_wide, depressed != "Depressed")$T1,dplyr::filter(Chao_wide, depressed != "Depressed")$T4, paired = T)

Sim_wide <- richness_meta %>% dplyr::select(Simpson,Time.x, id.x,depressed) %>% dcast(id.x + depressed ~ Time.x, value.var = "Simpson")

#overall difference
wilcox.test(dplyr::filter(Sim_wide, depressed != "Depressed")$T1,dplyr::filter(Sim_wide, depressed != "Depressed")$T3, paired = F)
wilcox.test(dplyr::filter(Sim_wide, depressed != "Depressed")$T1,dplyr::filter(Sim_wide, depressed != "Depressed")$T4, paired = F)
wilcox.test(dplyr::filter(Sim_wide, depressed != "Depressed")$T1,dplyr::filter(Sim_wide, depressed != "Depressed")$T5, paired = F)

#By Cluster difference
Sim_wide$Cluster <- richness_meta$Cluster[match(Sim_wide$id.x,richness_meta$id.x)]
Cluster_num <- 3
wilcox.test(dplyr::filter(Sim_wide, depressed != "Depressed"& Cluster == Cluster_num)$T1,dplyr::filter(Sim_wide, depressed == "Depressed"& Cluster == Cluster_num)$T1, paired = F)
wilcox.test(dplyr::filter(Sim_wide, depressed != "Depressed"& Cluster == Cluster_num)$T3,dplyr::filter(Sim_wide, depressed == "Depressed"& Cluster == Cluster_num)$T3, paired = F)
wilcox.test(dplyr::filter(Sim_wide, depressed != "Depressed"& Cluster == Cluster_num)$T1,dplyr::filter(Sim_wide, depressed != "Depressed"& Cluster == Cluster_num)$T5, paired = F)

ACE_wide <- richness_meta %>% dplyr::select(ACE,Time.x, id.x,depressed) %>% dcast(id.x + depressed ~ Time.x, value.var = "ACE")
wilcox.test(dplyr::filter(ACE_wide, depressed != "Depressed")$T1,dplyr::filter(ACE_wide, depressed != "Depressed")$T3, paired = F)
wilcox.test(dplyr::filter(ACE_wide, depressed == "Depressed")$T1,dplyr::filter(ACE_wide, depressed == "Depressed")$T3, paired = F)

ACE_wide

ACE_wide$Cluster <- richness_meta$Cluster[match(ACE_wide$id.x,richness_meta$id.x)]

Cluster_num = 3
t.test(dplyr::filter(ACE_wide, Cluster == Cluster_num)$T1,dplyr::filter(ACE_wide, Cluster ==Cluster_num)$T3, paired = F)

Shannon_wide <- richness_meta %>% dplyr::select(Shannon,Time.x, id.x,depressed) %>% dcast(id.x + depressed ~ Time.x, value.var = "Shannon")
wilcox.test(dplyr::filter(Shannon_wide, depressed != "Depressed")$T1,dplyr::filter(Shannon_wide, depressed != "Depressed")$T4, paired = T)

wilcox.test(dplyr::filter(Shannon_wide, depressed != "Depressed")$T1,dplyr::filter(Shannon_wide, depressed != "Depressed")$T3, paired = F)
wilcox.test(dplyr::filter(Shannon_wide, depressed == "Depressed")$T1,dplyr::filter(Shannon_wide, depressed == "Depressed")$T3, paired = F)

Shannon_wide$Cluster <- richness_meta$Cluster[match(Shannon_wide$id.x,richness_meta$id.x)]

Cluster_num = 3
wilcox.test(dplyr::filter(Shannon_wide, depressed != "Depressed" & Cluster == Cluster_num)$T1,dplyr::filter(Shannon_wide, depressed != "Depressed" & Cluster == Cluster_num)$T2, paired = F)
wilcox.test(dplyr::filter(Shannon_wide, depressed == "Depressed" & Cluster == Cluster_num)$T1,dplyr::filter(Shannon_wide, depressed == "Depressed" & Cluster == Cluster_num)$T4, paired = F)
wilcox.test(dplyr::filter(Shannon_wide, Cluster == Cluster_num)$T1, dplyr::filter(Shannon_wide, Cluster == Cluster_num)$T5)

summary(aov(Shannon ~ Time.x, data = filter(richness_meta, Cluster== 3)))

#Diversity group comparsion by enterotype

ggplot(richness_meta,aes(x = as.character(Time.x), y=Observed, color=depressed)) + geom_boxplot() + facet_grid(Cluster~depressed)

filter(richness_meta, depressed != "Depressed") %>% ggplot(aes(x = Time.x, y=Observed)) + geom_boxplot() + facet_grid(Cluster~.) + stat_compare_means(method = "anova")
filter(richness_meta, depressed == "Depressed") %>% ggplot(aes(x = Time.x, y=Observed)) + geom_boxplot() + facet_grid(Cluster~.) + stat_compare_means(method = "anova")

ggline(richness_meta, x = "Time.x", y = "ACE", color = "Cluster", add = c("mean_se", "jitter"),size = 0.7,point.size = 0.8) + scale_color_manual(values = entero_color_1) + facet_grid(.~depressed)
ggline(richness_meta, x = "Time.x", y = "Shannon", color = "Cluster", add = c("mean_se", "jitter"),size = 0.7,point.size = 0.8) + scale_color_manual(values = entero_color_1)  + facet_grid(.~depressed)

res.aov2 <- aov(ACE ~ 1 + Cluster * depressed, data = filter(richness_meta, Time.x == "T1"&Cluster != 2))
res.aov2 <- aov(Chao1 ~ 1 + Cluster * depressed, data = filter(richness_meta, Time.x == "T3"&Cluster != 2))
res.aov2 <- aov(Chao1 ~ 1 + Cluster * depressed, data = filter(richness_meta, Time.x == "T4"&Cluster != 2))
res.aov2 <- aov(Chao1 ~ 1 + Cluster * depressed, data = filter(richness_meta, Time.x == "T5"&Cluster != 2))

res.aov2 <- aov(Shannon ~ 1 +  Time.x, data = filter(richness_meta,  Cluster == 3))
summary(res.aov2)

#Anova(res.aov2, type = "III")

###########################################
#Baseline level age, sex, BMI
# 
# Age.data <- read.csv("./Data/Metab_psych_data.csv", header = T) %>% dplyr::select(id:safe)
# BMI.data <- read.csv("./Data/SchoolForTheWork-BaselineDemographics_DATA_LABELS.csv", header = T)
# 
# BMI.data$Height_m <- sapply(BMI.data$Height, function(x) {
#   # Split height into feet and inches
#   parts <- unlist(strsplit(x, "'"))
#   feet <- as.numeric(parts[1])
#   inches <- as.numeric(parts[2])
#   # Convert to meters: (feet * 12 + inches) * 0.0254
#   height_m <- (feet * 12 + inches) * 0.0254
#   return(height_m)
# })
# BMI.data$Weight <- as.numeric(BMI.data$Weight)
# BMI.data$Weight_kg <- BMI.data$Weight * 0.453592
# BMI.data$BMI <- BMI.data$Weight_kg / (BMI.data$Height_m ^ 2)
# BMI.data$BMI
# 
# extendedmeta <- BMI.data %>% dplyr::select(ID.number, Sex, BMI) %>% left_join(Age.data, by=c("ID.number"= "id")) %>% 
#   dplyr::select(ID.number,Sex,BMI,age,sex) %>%  distinct(ID.number, Sex, BMI, age, sex, .keep_all = TRUE) 
# extendedmeta <- extendedmeta %>% filter(!is.na(ID.number))
# extendedmeta <- extendedmeta %>% filter(!duplicated(extendedmeta$ID.number))
#
#write.csv(file = "./Data/Age_SEX_BMI.csv",extendedmeta)

extendedmeta <- read.csv(file = "./Data/Age_SEX_BMI.csv", header = T, row.names = 1)
extendedmeta$ID.number <- as.character(extendedmeta$ID.number)
meta.data$id.x <- as.character(meta.data$id.x)
new_meta <- dplyr::left_join(extendedmeta, meta.data, by=c("ID.number" = "id.x"))


sum(is.na(extendedmeta$age))
sum(is.na(extendedmeta$BMI))

extendedmeta
unique(upset_matrix_renamed$ID)

sort(extendedmeta$ID.number)
sort(unique(upset_matrix_renamed$ID))

difference.id <- c("45","47","49", "60", "63","66", "67")

filter(extendedmeta, ID.number %in% difference.id)
# Summary statistics for age, BMI, and sex distribution across clusters
new_meta %>%
  group_by(Cluster) %>%
  summarize(
    avg_age = mean(age, na.rm = TRUE),
    avg_BMI = mean(BMI, na.rm = TRUE),
    count_male = sum(Sex == "Male", na.rm = TRUE),
    count_female = sum(Sex == "Female", na.rm = TRUE)
  )

# Visualization: Age distribution across clusters
p.age <- ggplot(new_meta, aes(x = as.factor(Cluster), y = age)) +
  geom_boxplot() + theme_minimal() +
  labs(title = "Age Distribution Across Clusters", x = "Cluster", y = "Age")
p.age

# Visualization: BMI distribution across clusters
p.BMI <- ggplot(new_meta, aes(x = as.factor(Cluster), y = BMI)) +
  geom_boxplot() +theme_minimal() +
  labs(title = "BMI Distribution Across Clusters", x = "Cluster", y = "BMI")
p.BMI

# Visualization: Sex distribution across clusters
p.sex <- ggplot(new_meta, aes(x = as.factor(Cluster), fill = Sex)) +
  geom_bar(position = "fill") +theme_minimal() +
  labs(title = "Sex Distribution Across Clusters", x = "Cluster", y = "Proportion")
p.sex

p.reviewr2.1 <- p.age + p.BMI + p.sex
p.reviewr2.1

min(new_meta$age, na.rm = T)
max(new_meta$age, na.rm = T)
median(new_meta$age, na.rm = T)
mean(new_meta$age, na.rm = T)
sd(new_meta$age, na.rm = T)

min(new_meta$BMI, na.rm = T)
max(new_meta$BMI, na.rm = T)
median(new_meta$BMI, na.rm = T)
mean(new_meta$BMI, na.rm = T)
sd(new_meta$BMI, na.rm = T)

table(new_meta$Sex, new_meta$X.x)

#ggsave(filename = "./Result/revision/reviewer2.1.pdf",p.reviewr2.1, width = 12, height = 5, dpi = 300)

# Create a contingency table for Sex and Cluster
sex_cluster_table <- table(new_meta$Sex, new_meta$Cluster)

# Perform Chi-square test
chisq.test(sex_cluster_table)

# ANOVA test for age across clusters
summary(aov(age ~ as.factor(Cluster), data = new_meta))

# If data is not normally distributed, use Kruskal-Wallis test
kruskal.test(age ~ as.factor(Cluster), data = new_meta)

# ANOVA test for BMI across clusters
summary(aov(BMI ~ as.factor(Cluster), data = new_meta))

# Kruskal-Wallis test for BMI if not normally distributed
kruskal.test(BMI ~ as.factor(Cluster), data = new_meta)

# Summary statistics for age, BMI, and sex distribution for depression status
new_meta %>%
  group_by(depressed) %>%
  summarize(
    avg_age = mean(age, na.rm = TRUE),
    avg_BMI = mean(BMI, na.rm = TRUE),
    count_male = sum(Sex == "Male", na.rm = TRUE),
    count_female = sum(Sex == "Female", na.rm = TRUE)
  )

# Chi-square test for sex and depression status
sex_depressed_table <- table(new_meta$Sex, new_meta$depressed)
chisq.test(sex_depressed_table)

# ANOVA/Kruskal-Wallis for age and BMI by depression status
ggplot(new_meta, aes(x=depressed, y=age)) + geom_boxplot() + facet_wrap(.~Cluster)

summary(aov(age ~ depressed + Cluster, data = new_meta))
kruskal.test(age ~ depressed, data = new_meta)

summary(aov(BMI ~ depressed, data = new_meta))
kruskal.test(BMI ~ depressed, data = new_meta)

#######between and within distance
dist.ASV.bc <- vegan::vegdist(ASV_table, "bray", upper = T) %>% as.matrix() %>% melt()

dist.ASV.bc$from_id <- Sample$id.x[match(dist.ASV.bc$Var1, Sample$Row.names)]
dist.ASV.bc$to_id <- Sample$id.x[match(dist.ASV.bc$Var2, Sample$Row.names)]

within_distance <- filter(dist.ASV.bc, to_id==from_id) %>% filter(Var1 != Var2)
bet_distance <- filter(dist.ASV.bc, to_id != from_id)

within_distance$from_time <- meta.data$Time.x[match(within_distance$Var1, meta.data$Row.names)]
within_distance$to_time <- meta.data$Time.x[match(within_distance$Var2, meta.data$Row.names)]

within_distance$dist_raw <- paste(within_distance$from_time,within_distance$to_time, sep="_")

dplyr::filter(within_distance, dist_raw == "T1_T4")

T1T_list <- c(dplyr::filter(within_distance, dist_raw == "T1_T4")$Var1,dplyr::filter(within_distance, dist_raw == "T1_T4")$Var2)


upper.list <- c("T1_T3", "T1_T4", "T1_T5","T3_T4", "T3_T5","T4_T5")

within_distance$dist <- "lower"
within_distance$dist[within_distance$dist_raw %in% upper.list] <- "upper"

within_distance$entero <- meta.data$Cluster[match(within_distance$from_id, meta.data$id.x)]
within_distance$depre <- meta.data$depressed[match(within_distance$from_id, meta.data$id.x)]

p_within_all <- ggplot(within_distance, aes(x=as.factor(depre), y=value, fill=depre)) + 
  geom_boxplot()  + stat_compare_means() + base_theme +  ylab("Within Individual Distance") + facet_wrap(dist_raw~.)
p_within_all
#ggsave(filename = "./Result/within_dis_all.pdf", p_within_all, width = 4, height = 4, dpi = 300)

p_within <- ggplot(within_distance, aes(x=as.factor(depre), y=value, fill=depre)) + 
  geom_boxplot() + facet_wrap(.~entero) + stat_compare_means(label = "p.signif") + base_theme + scale_color_manual(values = depress_color) + ylab("Within Individual Distance")
p_within <- p_within + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p_within

#ggsave(filename = "./Result/within_dis.pdf", p_within, width = 4, height = 4, dpi = 300)

within_dist_data <- within_distance %>%
  filter(dist_raw %in% c("T1_T3", "T1_T4", "T1_T5")) %>% 
  group_by(dist_raw, depre) %>%
  summarise(
    mean_BCdistance = mean(value, na.rm = TRUE),
    sem_BCdistance = sd(value, na.rm = TRUE) / sqrt(n()))

p.within_dist <- within_dist_data %>%
  ggplot(aes(x = dist_raw, y = mean_BCdistance, group = depre)) +
  geom_line(aes(color = depre), linewidth = 1.2) +
  geom_point(aes(color = depre)) +
  geom_errorbar(aes(ymin = mean_BCdistance - sem_BCdistance, ymax = mean_BCdistance + sem_BCdistance, color = depre), width = 0.2) +
  labs(y = "within_dist", x = "Distance") +
  theme_minimal()+ scale_color_manual(values = depress_color)

p.within_dist

p.Chao1 / p.within_dist

#ggsave(filename = "./Result/within_dist_Time.pdf", p.within_dist,width = 6, height = 3, dpi = 300)

withindist_anova <- aov(value  ~ dist_raw * depre, data =  filter(within_distance,dist_raw %in% c("T1_T3", "T1_T4", "T1_T5")))
summary(withindist_anova)

# my_comparisons <- list( c("2", "1"), c("1", "3"), c("3", "2") )
 # ggplot(filter(within_distance, dist == "upper"), aes(x=as.factor(entero), y=value)) + geom_boxplot() + facet_wrap(.~depre) +
 #   stat_compare_means() + base_theme

########################################################PERMANOVA
physeqGenus <- tax_glom(physeq, taxrank="genus")

physeq_ord <- ordinate(physeqGenus,"PCoA", "bray")
p3 = plot_ordination(physeqGenus, physeq_ord, type="sample", color="depressed", title="taxa")
p3 <- p3 + facet_wrap(.~Time.y) +ggtitle("Depression by time")+ coord_equal() + base_theme + 
  stat_ellipse(aes(color=depressed, group=depressed), level=0.95)
p3

physeq_freq <- transform_sample_counts(physeq, function(x) x / sum(x))

 # physeq_freq_T1 <- subset_samples(physeq_freq, Time.x == "T1")
 # physeq_freq_T3 <- subset_samples(physeq_freq, Time.x == "T3")
 # physeq_freq_T4 <- subset_samples(physeq_freq, Time.x == "T4")
 # physeq_freq_T5 <- subset_samples(physeq_freq, Time.x == "T5")
#When separating, none of the single time point show significance

physeq_PERMANOVA <- physeq_freq
df = as(sample_data(physeq_PERMANOVA), "data.frame")
d = phyloseq::distance(physeq_PERMANOVA,"bray", "sample")
physeq.ST.adonis = adonis2(d ~  Time.x + depressed, df, permutations= 9999) 
physeq.ST.adonis

#######Prameters by enterotype
Sample$Cluster <- as.factor(Sample$Cluster)
phy.para <- colnames(Sample)[14:62]

############
my_comparisons <- list( c("1", "2"), c("1", "3"), c("2", "3"))
# for (i in 1:49){
#   print(phy.para[i])
#   marker <- phy.para[i]
#   p.compare.De <- Sample %>%
#     filter(depressed == "Depressed") %>%
#     ggplot(aes_string(x="Cluster", y=marker)) + geom_jitter(size=0.5) +
#     geom_boxplot(outlier.alpha = 0, alpha=0.5) + base_theme + xlab("Depressed") +
#     stat_compare_means(comparisons = my_comparisons,color="red")
#   p.compare.De
# 
#   p.compare.ND <- Sample %>%
#     filter(depressed != "Depressed") %>%
#     ggplot(aes_string(x="Cluster", y=marker)) + geom_jitter(size=0.5) +
#     geom_boxplot(outlier.alpha = 0, alpha=0.5) +  base_theme + xlab("Not Depressed") +
#     stat_compare_means(comparisons = my_comparisons, color="red")
#   p.compare.ND
# 
#   p.para <- p.compare.De + p.compare.ND
#   p.para <- p.para & ggtitle(marker)
#   path <- paste0("./Result/para.by.entero/", marker, ".pdf")
#   ggsave(filename = path, p.para, width = 4, height = 4, dpi = 300)
# }


###########Only T1
# my_comparisons <- list( c("1", "2"), c("1", "3"), c("2", "3"))
# for (i in 1:49){
#   print(phy.para[i])
#   marker <- phy.para[i]
#   p.compare.De2 <- Sample %>%
#     filter(depressed == "Depressed") %>% filter(Time.x == "T1") %>%
#     ggplot(aes_string(x="Cluster", y=marker)) + geom_jitter(size=0.5) +
#     geom_boxplot(outlier.alpha = 0, alpha=0.5) + base_theme + xlab("Depressed") +
#     stat_compare_means(comparisons = my_comparisons,color="red") 
#   p.compare.De2
# 
#   p.compare.ND2 <- Sample %>%
#     filter(depressed != "Depressed") %>% filter(Time.x == "T1") %>%
#     ggplot(aes_string(x="Cluster", y=marker)) + geom_jitter(size=0.5) +
#     geom_boxplot(outlier.alpha = 0, alpha=0.5) +  base_theme + xlab("Not Depressed") +
#     stat_compare_means(comparisons = my_comparisons, color="red") 
#   p.compare.ND2
# 
#   p.para2 <- p.compare.De2 + p.compare.ND2
#   p.para2 <- p.para2 & ggtitle(marker)
#   path <- paste0("./Result/para.by.entero/T1Only", marker, ".pdf")
#   ggsave(filename = path, p.para2, width = 4, height = 4, dpi = 300)
# }
# 
# my_comparisons <- list( c("1", "2"), c("1", "3"), c("2", "3"))
# p.compare.De <- Sample %>%
#   filter(depressed == "Depressed") %>% 
#   ggplot(aes(x=as.factor(Cluster), y=bdi_total)) + geom_jitter() +
#   geom_boxplot(outlier.alpha = 0, alpha=0.5) + base_theme + xlab("Depressed") +
#   stat_compare_means(comparisons = my_comparisons,color="red") + ylim(0,50)
# p.compare.De
# 
# p.compare.ND <- Sample %>%
#   filter(depressed != "Depressed") %>% 
#   ggplot(aes(x=as.factor(Cluster), y=bdi_total)) + geom_jitter() +
#   geom_boxplot(outlier.alpha = 0, alpha=0.5) +  base_theme + xlab("Not Depressed") +
#   stat_compare_means(comparisons = my_comparisons, color="red")+ ylim(0,50)
# p.compare.ND
# 
# p.compare.De + p.compare.ND
# 
# p.compare.De1 <- Sample %>%
#   filter(depressed == "Depressed") %>% filter(Time.x == "T1") %>% 
#   ggplot(aes(x=as.factor(Cluster), y=bdi_total)) + geom_jitter() +
#   geom_boxplot(outlier.alpha = 0, alpha=0.5) + base_theme + xlab("Depressed") +
#   stat_compare_means(comparisons = my_comparisons,color="red") + ylim(0,50)
# p.compare.De1
# 
# p.compare.ND1 <- Sample %>%
#   filter(depressed != "Depressed") %>% filter(Time.x == "T1") %>% 
#   ggplot(aes(x=as.factor(Cluster), y=bdi_total)) + geom_jitter() +
#   geom_boxplot(outlier.alpha = 0, alpha=0.5) +  base_theme + xlab("Not Depressed") +
#   stat_compare_means(comparisons = my_comparisons, color="red")+ ylim(0,50)
# p.compare.ND1
# 
# p.compare.De1 + p.compare.ND1
# 
# p.compare.combined <- Sample %>%
#   filter(Time.x == "T1") %>% 
#   ggplot(aes(x=as.factor(Cluster), y=bdi_total)) + geom_jitter() +
#   geom_boxplot(outlier.alpha = 0, alpha=0.5) +  base_theme  +
#   stat_compare_means(comparisons = my_comparisons, color="red")+ ylim(0,60) + ggtitle("BDI score Only T1")
# p.compare.combined
# 
# result <- Sample %>% filter(Time.x == "T1") %>% 
#   group_by(Cluster) %>%
#   summarize(
#     mean_bdi_total = mean(bdi_total, na.rm = TRUE),  # Calculate mean, excluding NA values
#     sd_bdi_total = sd(bdi_total, na.rm = TRUE)       # Calculate standard deviation, excluding NA values
#   )
# 
# p.compare.combined.all <- Sample %>%
#   ggplot(aes(x=as.factor(Cluster), y=bdi_total)) + geom_jitter() +
#   geom_boxplot(outlier.alpha = 0, alpha=0.5) +  base_theme  +
#   stat_compare_means(comparisons = my_comparisons, color="red")+ ylim(0,60)+ ggtitle("BDI score all timepoint")
# p.compare.combined.all
# 
# p.compare.combined.No_T1 <- Sample %>%
#   filter(Time.x != "T1") %>% 
#   ggplot(aes(x=as.factor(Cluster), y=bdi_total)) + geom_jitter() +
#   geom_boxplot(outlier.alpha = 0, alpha=0.5) +  base_theme  +
#   stat_compare_means(comparisons = my_comparisons, color="red")+ ylim(0,60)+ ggtitle("BDI score without T1")
# p.compare.combined.No_T1
# p.compare.combined.all + p.compare.combined.No_T1
# 


# Print the result
print(result)

list(Sample %>% filter(depressed == "Depressed") %>% filter(Cluster==3) %>% dplyr::select(bdi_total),
Sample %>% filter(depressed != "Depressed") %>% filter(Cluster==3) %>% dplyr::select(bdi_total))

list1 <- c(0, 1, 0, 0, 0, 7, 15, 1, 21, 0, 0, 0, 1, 0, 15, 1, 5, 26, 0, 2)
list2 <- c(0, 4, 4, 4, 3, 3, 11, 0, 6, 0, 2, 5, 2, 1, 10, 0)

t.test(list1, list2)

#per Ariel's request
p.compare.test <- Sample %>% filter(Row.names %in% T1T_list) %>%
  ggplot(aes(color=as.factor(Cluster), y=bdi_total, x=Time.x)) + geom_jitter() +
  geom_boxplot(outlier.alpha = 0, alpha=0.5) + base_theme + xlab("Time") + scale_color_manual(values = entero_color_1)
p.compare.test <- p.compare.test + theme(legend.position = "none") + facet_wrap(.~depressed)
p.compare.test
#ggsave(filename = "Result/bdi.bytime.entero.pdf", p.compare.test, width = 6, height = 5, dpi = 300)

p.compare.test1 <- Sample %>% filter(Row.names %in% T1T_list) %>%
  ggplot(aes(color=as.factor(Cluster), y=positive_emotion, x=Time.x)) + geom_jitter() +
  geom_boxplot(outlier.alpha = 0, alpha=0.5) + base_theme + xlab("Time") + scale_color_manual(values = entero_color_1)
p.compare.test1 <- p.compare.test1 + theme(legend.position = "none") + facet_wrap(.~depressed)
p.compare.test1
#ggsave(filename = "Result/poemotion.bytime.entero.pdf", p.compare.test1, width = 6, height = 5, dpi = 300)

p.compare.test2 <- Sample %>% 
  ggplot(aes(color=as.factor(Cluster), y=negative_emotion, x=Time.x)) + geom_jitter() +
  geom_boxplot(outlier.alpha = 0, alpha=0.5) + base_theme + xlab("Time") + scale_color_manual(values = entero_color_1)
p.compare.test2 <- p.compare.test2 + theme(legend.position = "none") + facet_wrap(.~depressed)
p.compare.test2
#ggsave(filename = "Result/neemotion.bytime.entero.pdf", p.compare.test2, width = 6, height = 5, dpi = 300)
# 
# #change to same grid
# my_comparisons <- list(c("1", "2"), c("1", "3"), c("2", "3"))
# for (i in 1:49){
#   if (i==19){
#     i = i + 1
#     } 
#   else {
#   print(phy.para[i])
#   marker <- phy.para[i]
#   p.compare.grid <- Sample %>%
#   ggplot(aes_string(x="Cluster", y=marker)) + geom_jitter(size=0.3) +
#   geom_boxplot(outlier.alpha = 0, alpha=0.5) + base_theme + xlab("Depressed") +
#   stat_compare_means(comparisons = my_comparisons,color="red") + facet_grid(.~depressed)
#   p.compare.grid
#   path <- paste0("./Result/para.by.entero/", marker, "_grid.pdf")
#   ggsave(filename = path, p.compare.grid, width = 3, height = 3, dpi = 300)
#   }
# }

list1 <- c("alive","enticing","loneliness","safe")
list2 <- c("das17_total", "dependency", "engagement", "negative_emotion","pathway","perfectionism","positive_emotion","relationships")
list <- c(list1, list2)
# 
# for (i in 1:12){
#   marker <- list[i]
#   print(i)
#   p.compare.grid <- Sample %>% mutate(depressed = ifelse(depressed == "Not Depressed", "Not_Depre", depressed)) %>% 
#     ggplot(aes_string(x="Cluster", y=marker)) + geom_jitter(size=0.3) +
#     geom_boxplot(aes(color=Cluster), outlier.alpha = 0, alpha=0.5) + base_theme + xlab("Depressed") +
#     stat_compare_means(comparisons = my_comparisons,color="red") + facet_grid(.~depressed)+ scale_color_manual(values = entero_color_1)+
#     scale_x_discrete(labels=c("1" = "Ba", "2" = "Ru", "3" = "Pr"))
# 
#   p.compare.grid
# 
#   path <- paste0("./Result/para.by.entero/Final_", marker, "_grid.pdf")
#   path2 <- paste0("./Result/para.by.entero/Final_", marker, "_grid.png")
#   ggsave(filename = path, p.compare.grid, width = 3, height = 3, dpi = 300)
#   ggsave(filename = path2, p.compare.grid, width = 3, height = 3, dpi = 300)
#   }



#############PCA
set.seed(777)
library(factoextra)
pca.df <- dplyr::filter(Sample, Cluster == 1) %>% dplyr::select(all_of(phy.para), -total_score) 
pca.df <- pca.df[complete.cases(pca.df),]
res.pca <- prcomp(pca.df, scale = TRUE)
fviz_eig(res.pca)
p1.b <- fviz_pca_ind(res.pca,
             col.ind = Sample$depressed[match(rownames(pca.df), Sample$Row.names)],
             geom = c("point"),
             palette = c("#4DBBD5FF", "#DC0000FF"),
             addEllipses = TRUE,
             ellipse.type = "confidence",
             legend.title = "Depression",
             title ="Bacteroides",
             repel = TRUE)

pca.df <- filter(Sample, Cluster == 2) %>% dplyr::select(all_of(phy.para), -total_score) 
pca.df <- pca.df[complete.cases(pca.df),]
res.pca <- prcomp(pca.df, scale = TRUE)
fviz_eig(res.pca)
p2.f <- fviz_pca_ind(res.pca,
                     col.ind = Sample$depressed[match(rownames(pca.df), Sample$Row.names)],
                     geom = c("point"),
                     palette = c("#4DBBD5FF", "#DC0000FF"),
                     addEllipses = TRUE,
                     ellipse.type = "confidence",
                     legend.title = "Depression",
                     title ="Firmicutes",
                     repel = TRUE)


pca.df <- filter(Sample, Cluster == 3) %>% dplyr::select(all_of(phy.para), -total_score) 
pca.df <- pca.df[complete.cases(pca.df),]
res.pca <- prcomp(pca.df, scale = TRUE)
fviz_eig(res.pca)
p3.p <- fviz_pca_ind(res.pca,
                     col.ind = Sample$depressed[match(rownames(pca.df), Sample$Row.names)],
                     geom = c("point"),
                     palette = c("#4DBBD5FF", "#DC0000FF"),
                     addEllipses = TRUE,
                     ellipse.type = "confidence",
                     legend.title = "Depression",
                     title ="Prevotella",
                     repel = TRUE)

p.ent.de <- p1.b + p2.f + p3.p  + plot_layout(guides = "collect")
p.ent.de
#ggsave("./Result/Enterotype.depre.pdf",p.ent.de, height = 3.5, width = 10, dpi = 300)

#####################Per George's request, annotate the same plot via the severty of depression
#https://en.wikipedia.org/wiki/Beck_Depression_Inventory
Sample <- Sample %>% mutate(depression_severity = case_when(
    bdi_total >= 0 & bdi_total <= 13 ~ "minimal_depression",
    bdi_total >= 14 & bdi_total <= 19 ~ "mild_depression",
    bdi_total >= 20 & bdi_total <= 28 ~ "moderate_depression",
    bdi_total >= 29 ~ "severe_depression",
    TRUE ~ NA_character_ ))

Sample$depression_severity <- factor(Sample$depression_severity, levels = c("minimal_depression","mild_depression",
                                                                            "moderate_depression","severe_depression"))

pca.df1 <- dplyr::filter(Sample, Cluster == 1) %>% dplyr::select(all_of(phy.para), -total_score) 
pca.df1 <- pca.df1[complete.cases(pca.df1),]
res.pca1 <- prcomp(pca.df1, scale = TRUE)
fviz_eig(res.pca1)
p1.bc <- fviz_pca_ind(res.pca1,
                     col.ind = Sample$depression_severity[match(rownames(pca.df1), Sample$Row.names)],
                     geom = c("point"),
                     palette = c("#4DBBD5FF",  "#FF9999FF", "#FF6666FF", "#CC0000FF"),
                    # addEllipses = TRUE,
                    # ellipse.type = "confidence",
                    legend.title = "Depression",
                     title ="Bacteroides by severty",
                     repel = TRUE)
p1.bc

pca.df2 <- dplyr::filter(Sample, Cluster == 2) %>% dplyr::select(all_of(phy.para), -total_score) 
pca.df2 <- pca.df2[complete.cases(pca.df2),]
res.pca2 <- prcomp(pca.df2, scale = TRUE)
fviz_eig(res.pca2)
p2.bc <- fviz_pca_ind(res.pca2,
                     col.ind = Sample$depression_severity[match(rownames(pca.df2), Sample$Row.names)],
                     geom = c("point"),
                     palette = c("#4DBBD5FF",  "#FF9999FF", "#FF6666FF", "#CC0000FF"),
                     # addEllipses = TRUE,
                     # ellipse.type = "confidence",
                     legend.title = "Depression", 
                     title ="Firmicutes by severty",
                     repel = TRUE)
p2.bc

pca.df3 <- dplyr::filter(Sample, Cluster == 3) %>% dplyr::select(all_of(phy.para), -total_score) 
pca.df3 <- pca.df3[complete.cases(pca.df3),]
res.pca3 <- prcomp(pca.df3, scale = TRUE)
fviz_eig(res.pca3)
p3.bc <- fviz_pca_ind(res.pca3,
                      col.ind = Sample$depression_severity[match(rownames(pca.df2), Sample$Row.names)],
                      geom = c("point"),
                      palette = c("#4DBBD5FF",  "#FF9999FF", "#FF6666FF", "#CC0000FF"),
                      # addEllipses = TRUE,
                      # ellipse.type = "confidence",
                      legend.title = "Depression",
                      title ="Prevotella by severty",
                      repel = TRUE)
p3.bc

p.ent.de.c <- p1.bc + p2.bc + p3.bc + plot_layout(guides = "collect")
p.ent.de.c
#ggsave("./Result/Enterotype.depre.byserve.pdf",p.ent.de.c, height = 3.5, width = 10, dpi = 300)

Sample_T1 <- Sample %>% filter(Time.x == "T1") %>% 
  mutate(depression_severity_ini = case_when(
    bdi_total >= 0 & bdi_total <= 13 ~ "minimal_depression",
    bdi_total >= 14 & bdi_total <= 19 ~ "mild_depression",
    bdi_total >= 20 & bdi_total <= 28 ~ "moderate_depression",
    bdi_total >= 29 ~ "severe_depression",
    TRUE ~ NA_character_ )) %>% dplyr::select(id.x,depression_severity_ini)

Sample_T1

Sample$depression_severity_ini <- Sample_T1$depression_severity_ini[match(Sample$id.x, Sample_T1$id.x)]
Sample$depression_severity_ini <- factor(Sample$depression_severity_ini, levels = c("minimal_depression","mild_depression",
                                                                            "moderate_depression","severe_depression"))
p1.b.ic <- fviz_pca_ind(res.pca1,
                      col.ind = Sample$depression_severity_ini[match(rownames(pca.df1), Sample$Row.names)],
                      geom = c("point"),
                      palette = c("#4DBBD5FF",  "#FF9999FF", "#FF6666FF", "#CC0000FF"),
                      # addEllipses = TRUE,
                      # ellipse.type = "confidence",
                      legend.title = "Depression",
                      title ="Bacteroides by severty",
                      repel = TRUE)
p1.b.ic

p2.b.ic <- fviz_pca_ind(res.pca2,
                        col.ind = Sample$depression_severity_ini[match(rownames(pca.df2), Sample$Row.names)],
                        geom = c("point"),
                        palette = c("#4DBBD5FF",  "#FF9999FF", "#FF6666FF", "#CC0000FF"),
                        # addEllipses = TRUE,
                        # ellipse.type = "confidence",
                        legend.title = "Depression",
                        title ="Firmicutes by severty",
                        repel = TRUE)
p2.b.ic

p3.b.ic <- fviz_pca_ind(res.pca3,
                        col.ind = Sample$depression_severity_ini[match(rownames(pca.df3), Sample$Row.names)],
                        geom = c("point"),
                        palette = c("#4DBBD5FF",  "#FF9999FF", "#FF6666FF", "#CC0000FF"),
                        # addEllipses = TRUE,
                        # ellipse.type = "confidence",
                        legend.title = "Depression",
                        title ="Prevotella by severty",
                        repel = TRUE)
p3.b.ic

p.b.ic <- p1.b.ic+p2.b.ic+p3.b.ic + plot_layout(guides = "collect")
p.b.ic
#ggsave("./Result/Enterotype.depre.byserve.ini.pdf",p.b.ic, height = 3.5, width = 10, dpi = 300)

#######################
#Per George's request, severty and enterotype plot

p.severty.entero <- Sample %>%
  filter(!is.na(depression_severity_ini)) %>%
  # Calculate the proportions
  group_by(Cluster) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  # Create the plot
  ggplot(aes(x = Cluster, fill = depression_severity_ini)) +
  geom_bar(aes(y = ..prop.., group = depression_severity_ini), position = position_dodge()) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Percentage", title = "Percentage of Each Category of Depression Severity by Cluster") +
  theme_minimal()
p.severty.entero

#ggsave("./Result/severty.entero.pdf", p.severty.entero, height = 3.5, width = 10, dpi = 300)

############################################################
#Childhood Abuse
#####################################################
ACES <- read.csv("./Data/SchoolForTheWork-Aces_DATA_2023-10-24_1933 1.csv", header = T)
ACES

Sample_T2 <- left_join(Sample, ACES, by=c("record_id.x" = "record_id"))

p.ace.cluster <- Sample_T2 %>% dplyr::filter(Time.x == "T1") %>% ggplot(aes(x=Cluster, y=ace_total)) + geom_jitter() + geom_boxplot(alpha=0.8)
p.ace.cluster <- p.ace.cluster + facet_wrap(Time.x~depressed, ncol = 2)
p.ace.cluster

#ggsave("./ACE_Enterotype.pdf", p.ace.cluster, height = 4,width = 3, dpi = 300)

psy.ace <- dplyr::filter(Sample_T2, Time.x == "T1") %>% dplyr::select(bdi_total:psych_meds, ace_total,-total_score)
dim(psy.ace)

cor_matrix <- cor(psy.ace, method = "spearman", use = "pairwise.complete.obs")

cor_long <- melt(cor_matrix)

ggplot(cor_long, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() + # Creates a heatmap
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab", name="Spearman Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), axis.text.y = element_text(angle=0, vjust=0.5)) +
  labs(x = "", y = "")

p.ace <- pheatmap(cor_matrix, clustering_distance_rows = "correlation", clustering_method = "complete")
p.ace

ggsave(filename = "./Result/ACE_otherpsydata.pdf",p.ace, width = 10, height = 10, dpi = 300)

#################Call by different time point(Needs to change the time point for this)
pca.df.all <- filter(Sample, Time.x== "T5") %>% dplyr::select(all_of(phy.para), -total_score) 
pca.df.all <- pca.df.all[complete.cases(pca.df.all),]
zero_variance_cols <- apply(pca.df.all, 2, var) == 0

zero_variance_colnames <- names(pca.df.all)[zero_variance_cols]

pca.df.all <- pca.df.all[, !zero_variance_cols]
res.pca.all <- prcomp(pca.df.all, scale = TRUE)
fviz_eig(res.pca.all)
p4.p <- fviz_pca_ind(res.pca.all,
                     col.ind = Sample$depressed[match(rownames(pca.df.all), Sample$Row.names)],
                     geom = c("point"),
                     palette = c("#4DBBD5FF", "#DC0000FF"),
                     addEllipses = TRUE,
                     ellipse.type = "confidence",
                     legend.title = "Depression",
                     title ="48 Psychological Measurements, T4",
                     repel = TRUE)
p4.p

###########NMDS
pca.df <- Sample %>%  dplyr::select(all_of(phy.para), -total_score) 
pca.df <- pca.df[complete.cases(pca.df),]

pca.df_dist <- vegdist(pca.df, method = "bray")
pca.df_dist <- as.matrix(pca.df_dist, labels = T)
pca.df_dist_NMS <-metaMDS(pca.df_dist,
          distance = "bray",
          k = 4,
          maxit = 9999, 
          trymax = 500,
          wascores = TRUE)

data_scores <- as.data.frame(scores(pca.df_dist_NMS,display="site"))
data_scores$group <- Sample$Cluster[match(rownames(data_scores), Sample$Row.names)]
data_scores$depressed <- Sample$depressed[match(rownames(data_scores), Sample$Row.names)]
data_scores$time <- Sample$Time.x[match(rownames(data_scores), Sample$Row.names)]

p <- ggplot(data_scores, aes(x=NMDS1, y=NMDS2, color=as.factor(depressed))) + geom_point() + stat_ellipse(level = 0.50)
p <- p + facet_grid(. ~ time)
p

p <- ggplot(data_scores, aes(x=NMDS1, y=NMDS2, color=as.factor(depressed))) + geom_point() + stat_ellipse(level = 0.50)
p <- p + facet_grid(. ~ group)
p

table(Sample$Cluster, Sample$depressed)

###seperate cluster and run permanova
pca.df <- filter(Sample, Cluster == 3) %>% dplyr::select(all_of(phy.para), -total_score) 
pca.df <- pca.df[complete.cases(pca.df),]
dim(pca.df)
group.ent <- Sample$Cluster[match(rownames(pca.df), Sample$Row.names)]
group.dep <- Sample$depressed[match(rownames(pca.df), Sample$Row.names)]

pca.df_dist <- vegdist(pca.df, method = "bray")

adonis2(pca.df_dist ~ group.dep, method="bray",perm=9999)

#######################################################################################################
#Cytokine Analysis  
Cytokine.df <- merge(Cytokine, genus.meta, by.x=c("id", "Time"), by.y = c("id.x", "Time.x"))
cytokine.list <- colnames(Cytokine.df)[8:87]

Cytokine.df$Cluster <- as.character(Cytokine.df$Cluster)

# my_comparisons <- list( c("1", "2"), c("1", "3"), c("2", "3"))
# for (i in 1:80){
#   cytokine.id <- cytokine.list[i]
#   print(cytokine.id)
#   p.dp <- filter(Cytokine.df, depressed == "Depressed") %>% 
#     ggplot(aes_string(x= "Cluster", y=cytokine.id,group="Cluster")) + 
#     geom_jitter() + geom_boxplot(outlier.alpha = 0, alpha=0.5) + base_theme + xlab("Depressed") +
#     stat_compare_means(comparisons = my_comparisons,color="red")
#   
#   p.np <- filter(Cytokine.df, depressed != "Depressed") %>% 
#     ggplot(aes_string(x= "Cluster", y=cytokine.id,group="Cluster")) + 
#     geom_jitter() + geom_boxplot(outlier.alpha = 0, alpha=0.5) + base_theme + xlab("Not Depressed") +
#     stat_compare_means(comparisons = my_comparisons,color="red")
#   
#   p.para2 <- p.dp + p.np
#   p.para2 <- p.para2 & ggtitle(cytokine.id)
#   path2 <- paste0("./Result/cytokine.by.entero/", cytokine.id, ".pdf")
#   ggsave(filename = path2, p.para2, width = 5, height = 5, dpi = 300)
# }
#  
#  my_comparisons_time <- list(c("T1", "T3"), c("T1", "T4"), c("T1", "T5"), c("T3", "T4"), c("T3","T5"), c("T4","T5"))
# for (i in 1:80){
#   cytokine.id <- cytokine.list[i]
#   print(cytokine.id)
# 
#   p2.dp <- filter(Cytokine.df, depressed == "Depressed") %>%
#   ggplot(aes_string(x= "Time", y=cytokine.id,group="Time")) +
#   geom_jitter() + geom_boxplot(outlier.alpha = 0, alpha=0.5) + base_theme + xlab("Depressed") +
#   stat_compare_means(comparisons = my_comparisons_time,color="red")
# 
#   p2.np <- filter(Cytokine.df, depressed != "Depressed") %>%
#   ggplot(aes_string(x= "Time", y=cytokine.id,group="Time")) +
#   geom_jitter() + geom_boxplot(outlier.alpha = 0, alpha=0.5) + base_theme + xlab("Not Depressed") +
#   stat_compare_means(comparisons = my_comparisons_time,color="red")
# 
#   p2.para2 <- p2.dp + p2.np
#   p2.para2 <- p2.para2 & ggtitle(cytokine.id)
#   path3 <- paste0("./Result/cytokine.by.time/", cytokine.id, ".pdf")
#   ggsave(filename = path3, p2.para2, width = 5, height = 5, dpi = 300)
# }

#change cytokine to same grid
my_comparisons_time <- list(c("T1", "T3"), c("T1", "T4"), c("T1", "T5"), c("T3", "T4"), c("T3","T5"), c("T4","T5"))
for (i in 1:80){
  cytokine.id <- cytokine.list[i]
  print(cytokine.id)
  
  marker <- cytokine.id
  
  p.compare.grid <- Cytokine.df %>% mutate(depressed = ifelse(depressed == "Not Depressed", "Not_Depre", depressed)) %>% 
    ggplot(aes_string(x="Cluster", y=marker)) + geom_jitter(size=0.3) +
    geom_boxplot(aes(color = Cluster),outlier.alpha = 0, alpha=0.5) + xlab("Depressed") + 
    facet_grid(.~depressed) + scale_color_manual(values = entero_color_1) + xlab(NULL)+  
    stat_compare_means(label = "p.signif",comparisons = my_comparisons,color="red", label.y = c(1.1,1.5,1.3))+ coord_cartesian(ylim = c(-1, 2))+
    scale_x_discrete(labels=c("1" = "Ba", "2" = "Ru", "3" = "Pr")) + base_theme + theme(legend.position = "none")
  p.compare.grid
  path <- paste0("./Result/cytokine.by.entero/grid/", marker, "_grid.pdf")
  ggsave(filename = path, p.compare.grid, width = 2, height = 2.5, dpi = 300)
}


#############################################################
for (i in 1:80){
cytokine.id <- cytokine.list[i]
print(cytokine.id)

marker <- cytokine.id

p.compare.grid <- Cytokine.df %>% mutate(depressed = ifelse(depressed == "Not Depressed", "Not_Depre", depressed)) %>% 
  ggplot(aes_string(x="Cluster", y=marker)) + geom_jitter(size=0.3) +
  geom_boxplot(aes(color = Cluster),outlier.alpha = 0, alpha=0.5) + xlab("Depressed") + 
  facet_grid(.~depressed) + scale_color_manual(values = entero_color_1) + xlab(NULL)+
  scale_x_discrete(labels=c("1" = "Ba", "2" = "Ru", "3" = "Pr")) + base_theme + theme(legend.position = "none")
p.compare.grid
path <- paste0("./Result/cytokine.by.entero/Nosig/", marker, "_grid.pdf")
ggsave(filename = path, p.compare.grid, width = 2, height = 2.5, dpi = 300)
}


cytokine.id <- "ENA78.CXCL5"
print(cytokine.id)

marker <- cytokine.id

p.compare.grid <- Cytokine.df %>% mutate(depressed = ifelse(depressed == "Not Depressed", "Not_Depre", depressed)) %>% 
  ggplot(aes_string(x="Cluster", y=marker)) + geom_jitter(size=0.3) +
  geom_boxplot(aes(color = Cluster),outlier.alpha = 0, alpha=0.5) + xlab("Depressed") + 
  facet_grid(.~depressed) + scale_color_manual(values = entero_color_1) + xlab(NULL)+
  scale_x_discrete(labels=c("1" = "Ba", "2" = "Ru", "3" = "Pr")) + base_theme + theme(legend.position = "none") +
  coord_cartesian(ylim = c(-1, 2))
p.compare.grid
path <- paste0("./Result/cytokine.by.entero/Nosig/", marker, "_grid.pdf")
ggsave(filename = path, p.compare.grid, width = 2, height = 2.5, dpi = 300)


#by time
for (i in 1:80){
  cytokine.id <- cytokine.list[i]
  print(cytokine.id)

  marker <- cytokine.id

  p.compare.grid <- Cytokine.df %>%
    ggplot(aes_string(x="Time", y=marker)) + geom_jitter(size=0.3) +
    geom_boxplot(outlier.alpha = 0, alpha=0.5) + base_theme + xlab("Depressed") +
    stat_compare_means(comparisons = my_comparisons_time,color="red") + facet_grid(.~depressed)
  p.compare.grid
  path <- paste0("./Result/cytokine.by.time/", marker, "_grid.pdf")
  ggsave(filename = path, p.compare.grid, width = 3, height = 3, dpi = 300)
}


my_comparisons_T3 <- list(c("T1", "T3"), c("T3", "T4"), c("T3","T5"))
#plot single cytokine by time
p2.p <- filter(Cytokine.df) %>%
      ggplot(aes_string(x= "Time", y="TARC.CCL17",group="Time")) +
      geom_jitter() + geom_boxplot(outlier.alpha = 0, alpha=0.5) + base_theme + xlab("Not Depressed") 
p2.p <- p2.p+ facet_grid(.~depressed) + stat_compare_means(comparisons = my_comparisons_T3, color="red") #+ ylim(-1,10) 
p2.p
    

#############For figure 3, y axis needs to be customized by each cytokine
#CCL17
#T1-T2 or T2-T3, T2-T4 significantly different
list.3 <- c( "GROA",  "PDGFAA", "ENA78.CXCL5","SCD40L","MCP2.CCL8", "IL21", "IL7","IL5","I309.CCL1", "FAS")

#T1-T3 or T1-T4
list.4 <- c("MDC.CCL22","RANTES.CCL5", "MIP1A.CCL3","BCA1.CXCL13", "TRAIL.TNFSF10","CTACK.CCL27","PAI.1.SERPIN.E")

#Maybe
list.5 <- c("TARC.CCL17", "FLT3L","FAS","MIG.CXCL9", "RESISTIN", "IL33.NFHEV.MATURE.", "SCF", "I309.CCL1", "MCP4.CCL13","MCP2.CCL8","SDF1A.B.CXCL12")

#T3 differ
T3.Differ <- c("IL17A.CTLA8", "IFNG", "MDC.CCL22","FLT3L")
T4.Differ <- c("TRAIL.TNFSF10", "MDC.CCL22","FAS")

cytokine.list

cytokine_plot <- "FAS"

p.cytokine <- filter(Cytokine.df) %>%
  ggplot(aes_string(x= "Time", y=cytokine_plot, group="depressed")) +
  geom_jitter(aes(color=depressed), alpha = 0.2) +
  stat_summary(fun.data = mean_se, geom="errorbar", width= 0.2, aes(color=depressed)) +
  stat_summary(fun = mean, geom="point", shape=18, size=4, aes(color=depressed)) +
  stat_summary(fun = mean, geom="line", aes(color=depressed), linewidth = 0.8) + scale_color_manual(values = depress_color)+
  base_theme + theme(legend.position = "none") + xlab(NULL) +
  coord_cartesian(ylim = c(-1, 2)) + scale_x_discrete(labels=c("T1", "T2", "T3", "T4"))

p.cytokine

#Two Way Anova
fit <- aov(FAS  ~ depressed * Time, data = Cytokine.df)
summary(fit)

#path <- paste0("./Figures/Figure 3/Cytokine_By_Time/Big_text",cytokine_plot, ".pdf" )
#ggsave(filename = path, p.cytokine, width = 4, height = 2.5, dpi = 300)


################Perm Cytokine
cyto_PERMANOVA <- dplyr::select(Cytokine.df, EGF:SVCAM1)
df = dplyr::select(Cytokine.df,id.y, Time, Cluster, depressed)
d = vegdist(cyto_PERMANOVA,"euclidean")

physeq.ST.adonis = adonis2(d ~ Time + Cluster + depressed, by = "terms", df, permutations= 9999)
physeq.ST.adonis

physeq.ST.adonis = adonis2(d ~ Cluster + Time + depressed, by = "terms", df, permutations= 9999)
physeq.ST.adonis

physeq.ST.adonis = adonis2(d ~ depressed + Cluster + Time , by = "terms", df, permutations= 9999)
physeq.ST.adonis

physeq.ST.adonis = adonis2(d ~ Time + Cluster + depressed, by = "margin", df, permutations= 9999)
physeq.ST.adonis

physeq.ST.adonis = adonis2(d ~ Cluster + Time + depressed, by = "margin", df, permutations= 9999)
physeq.ST.adonis

physeq.ST.adonis = adonis2(d ~ depressed + Cluster + Time , by = "margin", df, permutations= 9999)
physeq.ST.adonis


################Heatmap
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(dendsort)

ND.cytokine <- Cytokine.df %>% filter(depressed != "Depressed") %>% dplyr::select(Cluster, EGF:SVCAM1)
DP.cytokine <- Cytokine.df %>% filter(depressed == "Depressed") %>% dplyr::select(Cluster, EGF:SVCAM1)

ND.cytokine.heat <- aggregate(ND.cytokine, by=list(ND.cytokine$Cluster),FUN=mean) %>% dplyr::select(-Cluster)
ND.cytokine.heat$Depressed <- "Not Depressed"
DP.cytokine.heat <- aggregate(DP.cytokine, by=list(DP.cytokine$Cluster),FUN=mean) %>% dplyr::select(-Cluster)
DP.cytokine.heat$Depressed <- "Depressed"
# write.csv(file = "./Result/Heatmap/ND.cytokine.heat.csv", ND.cytokine.heat)
# write.csv(file = "./Result/Heatmap/DP.cytokine.heat.csv", DP.cytokine.heat)

combined.cytokine <- rbind(ND.cytokine.heat,DP.cytokine.heat)
dim(combined.cytokine)
combined.cytokine <- t(dplyr::select(combined.cytokine, -Group.1, -Depressed))
colnames(combined.cytokine) <- c("ND.B", "ND.R", "ND.P", "D.B", "D.R", "D.P")

combined.cytokine.heat <- as_tibble(combined.cytokine) %>% dplyr::select(ND.B, ND.R, ND.P,D.P, D.B, D.R)
rownames(combined.cytokine.heat) <- rownames(combined.cytokine)
#set color
col_fun = colorRamp2(c(-0.5, 1), c("white", "dark red"))

circos.par(gap.after =  30)
circos.heatmap(combined.cytokine.heat, 
               col = col_fun, 
               rownames.side = "outside",
               dend.side = "inside",
               )

circos.clear()

#Hclust without heatmap

# Calculate the distance matrix
cor_matrix <- cor(combined.cytokine)
dist_matrix <- as.dist(1 - cor_matrix)
hc <- hclust(dist_matrix, method = "average")
plot(hc, labels = colnames(combined.cytokine), main = "Hierarchical Clustering of combined.cytokine", xlab = "Cytokines", ylab = "Distance (1 - Pearson Correlation)")

Hclust.cytokine <- Cytokine.df %>% dplyr::select(EGF:SVCAM1) %>% t() 

##############K Means Cluster by sample############
#https://www.statology.org/k-means-clustering-in-r/
Kmean.cytokine <- Cytokine.df %>% dplyr::select(EGF:SVCAM1)
rownames(Kmean.cytokine) <- Cytokine.df$sampleID
#no need to scale as the cytokine data is already scaled
#Kmean.cytokine_scaled <- scale(Kmean.cytokine) %>% data.frame()
Kmean.cytokine_scaled <- Kmean.cytokine 
class(Kmean.cytokine_scaled$EOTAXIN.CCL11)

fviz_nbclust(Kmean.cytokine_scaled, kmeans, method = "wss")

gap_stat <- clusGap(Kmean.cytokine_scaled,
                    FUN = kmeans,
                    nstart = 25,
                    K.max = 10,
                    B = 50)
fviz_gap_stat(gap_stat)

set.seed(17)
km.cytokine <- kmeans(Kmean.cytokine_scaled, 2, nstart = 25)
km.cytokine

fviz_cluster(km.cytokine, data = Kmean.cytokine_scaled)

Cytokine.df$cytokine.km <- km.cytokine$cluster

table(Cytokine.df$cytokine.km, Cytokine.df$Cluster)
table(Cytokine.df$cytokine.km, Cytokine.df$Time)
table(Cytokine.df$cytokine.km, Cytokine.df$depressed)

Cytokine.df %>% ggplot(aes(x=as.factor(cytokine.km), y=IL1RA)) + geom_jitter() + geom_boxplot()
Cytokine.df %>% ggplot(aes(x=as.factor(cytokine.km), y=IL1B)) + geom_jitter() + geom_boxplot()

Cytokine.df %>% ggplot(aes(x=as.factor(cytokine.km), y=GROA)) + geom_jitter() + geom_boxplot()

##############K Means Cluster by subject############
#Kmean.cytokine_bysub <- Cytokine.df %>% dplyr::select(EGF:SVCAM1) %>% 

Kmean.cytokine_bysub <- Cytokine.df %>% dplyr::group_by(id) %>% summarise(across(EGF:SVCAM1, ~ mean(.x, na.rm = TRUE))) %>% column_to_rownames(var = "id")
Kmean.cytokine_bysub %>% dim()
Kmean.cytokine_bysub_scaled <- scale(Kmean.cytokine_bysub) %>% data.frame()

fviz_nbclust(Kmean.cytokine_bysub_scaled, kmeans, method = "wss")
gap_stat <- clusGap(Kmean.cytokine_bysub_scaled,
                    FUN = kmeans,
                    nstart = 25,
                    K.max = 10,
                    B = 50)
fviz_gap_stat(gap_stat)

set.seed(17)
km.cytokine.bysub <- kmeans(Kmean.cytokine_bysub_scaled, 2, nstart = 25)
km.cytokine.bysub

fviz_cluster(km.cytokine.bysub, data = Kmean.cytokine_bysub_scaled)
Cytokine.df.km <- merge(Cytokine.df, data.frame(km.cytokine.bysub$cluster), by.x = "id", by.y="row.names")

table(Cytokine.df.km$km.cytokine.bysub.cluster, Cytokine.df.km$Cluster)
table(Cytokine.df.km$km.cytokine.bysub.cluster, Cytokine.df.km$Time)
table(Cytokine.df.km$km.cytokine.bysub.cluster, Cytokine.df.km$depressed)

##############
#Prepare BRMS data
# genus_count <- otu_table(physeqGenus) %>% data.frame()
# taxa <- tax_table(physeq) %>% data.frame()
# taxa
# colnames(genus_count) <- taxa$genus[match(colnames(genus_count), taxa$ASV)]
# sample_count <- sample_data(physeqGenus) %>% data.frame()
# 
# sample_count2 <- dplyr::select(sample_count, Row.names,record_id.x,Time.x,id.x,depressed)
# 
# Cytokine.table <- merge(sample_count2, Cytokine.df, by.x = c("id.x", "Time.x"), by.y = c("id", "Time"))
# Cytokine.table <- dplyr::select(Cytokine.table, -bdi_total.x, -safe.x, -Wells, -plate)
# 
# write.csv(file = "./Data/Genus.count.table.csv", sample_count)
# write.csv(file = "./Data/Cytokine.table.csv", Cytokine.table)

sample_count <- read.csv(file = "./Data/BRMSmodel/Genus.count.table.csv", header = T, row.names = 1)
Cytokine.table <- read.csv(file = "./Data/BRMSmodel/Cytokine.table.csv", header = T, row.names = 1)

rownames(sample_count)
rownames(Cytokine.table)

sample_count$Acetobacter
sample_count$Bacteroides

df.cytokine.microbiome <- merge(sample_count,Cytokine.table, by="row.names") %>% column_to_rownames("Row.names")

df.cytokine.microbiome
colnames(df.cytokine.microbiome)[65]

#extract list of bacteria and cytokine
bacteria.genus.list <- dplyr::select(df.cytokine.microbiome, Bacteroides:Acetobacter) %>% colnames()
length(bacteria.genus.list)
cytokines.list <- dplyr::select(df.cytokine.microbiome, CHEX4:SVCAM1) %>% colnames()
length(cytokines.list)

bacteria.genus.list_10 <- bacteria.genus.list[colSums(dplyr::select(df.cytokine.microbiome, Bacteroides:Acetobacter) != 0) > 12] 

bacteria.genus.list[colSums(dplyr::select(df.cytokine.microbiome, Bacteroides:Acetobacter) != 0) <= 12] %>% length()

bacteria <- bacteria.genus.list[3]
bacteria

cytokines <- cytokines.list[3]
cytokines
# 
# #start the loop
# Cytokine.pvalue <- data.frame()
# Depression.pvalue <- data.frame()
# Cytokine.depression.pvalue <- data.frame()
# 
# Cytokine.pvalue_F <- data.frame()
# Depression.pvalue_F <- data.frame()
# Cytokine.depression.pvalue_F <- data.frame()
# 
# for (j in 10:198){
#   print(j)
#   bacteria <- bacteria.genus.list_10[j]
#   print(bacteria)
# 
#   for (i in 1:77){
#     print(i)
#     cytokines <- cytokines.list[i]
#     f <- as.formula(paste(bacteria,"~", cytokines, "* depressed.x", sep = ""))
#     print(f)
# 
#     fm <- mixed_model(fixed = f , random = ~ 1 | id.x, data = df.cytokine.microbiome,
#                       family = zi.negative.binomial(), zi_fixed = ~ 1,iter_EM = 0, max_coef_value = 1000)
# 
#     x <- summary(fm)
# 
#     Cytokine.pvalue[i,1] <- cytokines
#     Cytokine.pvalue[i,2] <- x$coef_table[14]
#     Cytokine.pvalue[i,3] <- x$coef_table[2]
#     Cytokine.pvalue$bacteria <- bacteria
# 
#     Depression.pvalue[i,1] <- cytokines
#     Depression.pvalue[i,2] <- x$coef_table[15]
#     Depression.pvalue[i,3] <- x$coef_table[3]
#     Depression.pvalue$bacteria <- bacteria
# 
#     Cytokine.depression.pvalue[i,1] <- cytokines
#     Cytokine.depression.pvalue[i,2] <- x$coef_table[16]
#     Cytokine.depression.pvalue[i,3] <- x$coef_table[4]
#     Cytokine.depression.pvalue$bacteria <- bacteria
#     }
# 
#   Cytokine.pvalue_F <- rbind(Cytokine.pvalue_F,Cytokine.pvalue)
#   Depression.pvalue_F <- rbind(Depression.pvalue_F,Depression.pvalue)
#   Cytokine.depression.pvalue_F <- rbind(Cytokine.depression.pvalue_F,Cytokine.depression.pvalue)
#   }
# 
# colnames(Cytokine.pvalue) <- c("Cytokine", "P-Value", "Estimate", "Genus")
# colnames(Depression.pvalue) <- c("Cytokine", "P-Value", "Estimate", "Genus")
# colnames(Cytokine.depression.pvalue) <- c("Cytokine", "P-Value", "Estimate", "Genus")


####Check Single bacteria and Cytokine
bacteria <- "Prevotella"
cytokines <- "FLT3L"

f <- as.formula(paste(bacteria,"~", cytokines, "* depressed.x", sep = ""))

fm <- mixed_model(fixed = f , random = ~ 1 | id.x, data = df.cytokine.microbiome,
                  family = zi.negative.binomial(), zi_fixed = ~ 1,iter_EM = 0)

x <- summary(fm)
x

############loop cytokine through the data
cytokine.time.lmer <- data.frame()
y <- data.frame()

for (i in 1:77){
  print(i)
  cytokines <- cytokines.list[i]
  print(cytokines)
  f <- as.formula(paste(cytokines,"~ Time.x.x * depressed.x + (1 | id.x)", sep = ""))
  print(f)
  lmer.df <-lmerTest::lmer(f, data = df.cytokine.microbiome)
  x <- summary(lmer.df)
  y <- x$coefficients %>% data.frame()
  y$cytokine <- cytokines
  cytokine.time.lmer <- rbind(cytokine.time.lmer, y)
  } 

cytokine.time.lmer

cytokine.time.lmer$cater <- row.names(cytokine.time.lmer)
cytokine.time.lmer <- cytokine.time.lmer %>% filter(cater != "(Intercept)")
cytokine.time.lmer
cytokine.time.lmer$clean_rowname <- sub("(T[3-5])\\d+", "\\1", cytokine.time.lmer$cater)
cytokine.time.lmer$clean_rowname <- sub("(?<![T[3-5]])\\d+$", "", cytokine.time.lmer$clean_rowname, perl=TRUE)

#write.csv(file = "./Result/cytokine.by.time/table.csv",cytokine.time.lmer)

filter(cytokine.time.lmer, Pr...t..  < 0.05) %>% filter(clean_rowname == "depressed.xNot Depressed")

filter(cytokine.time.lmer, Pr...t..  < 0.05) %>% filter(clean_rowname == "Time.x.xT")



#https://cran.r-project.org/web/packages/sjPlot/vignettes/tab_mixed.html

cytokines <- "BCA1.CXCL13"
print(cytokines)
f <- as.formula(paste(cytokines,"~ Time.x.x * depressed.x + (1 | id.x)", sep = ""))
print(f)
lmer.df <-lmerTest::lmer(f, data = df.cytokine.microbiome)
tab_model(lmer.df)
plot_model(lmer.df, show.values = TRUE, value.offset = .3)
lmer.df

################Entero_Switch
En_Sw_df <- dplyr::select(Cytokine.df,id, Cluster, Time, depressed)

###T1 to T3
D_1_3 <- En_Sw_df %>% 
  filter(depressed == "Depressed") %>%
  filter(Time %in% c("T1", "T3"))

D_1_3 <- En_Sw_df %>% 
  filter(depressed == "Depressed") %>%
  filter(Time %in% c("T1", "T3")) %>% 
  filter(id %in% D_1_3$id[duplicated(D_1_3$id)])

D_1_3_wd <- dcast(D_1_3, id ~Time, value.var = "Cluster")
D_1_3_wd$ES_13 <- paste(D_1_3_wd$T1,D_1_3_wd$T3, sep = "_")
table(D_1_3_wd$ES_13)

ND_1_3 <- En_Sw_df %>% 
  filter(depressed != "Depressed") %>%
  filter(Time %in% c("T1", "T3"))

ND_1_3 <- En_Sw_df %>% 
  filter(depressed != "Depressed") %>%
  filter(Time %in% c("T1", "T3")) %>% 
  filter(id %in% ND_1_3$id[duplicated(ND_1_3$id)])

ND_1_3_wd <- dcast(ND_1_3, id ~Time, value.var = "Cluster")
ND_1_3_wd$ES_13 <- paste(ND_1_3_wd$T1,ND_1_3_wd$T3, sep = "_")
table(ND_1_3_wd$ES_13)

###T3 to T4
D_3_4 <- En_Sw_df %>% 
  filter(depressed == "Depressed") %>%
  filter(Time %in% c("T3", "T4"))

D_3_4 <- En_Sw_df %>% 
  filter(depressed == "Depressed") %>%
  filter(Time %in% c("T3", "T4")) %>% 
  filter(id %in% D_3_4$id[duplicated(D_3_4$id)])

D_3_4_wd <- dcast(D_3_4, id ~Time, value.var = "Cluster")
D_3_4_wd$ES_34 <- paste(D_3_4_wd$T3,D_3_4_wd$T4, sep = "_")
table(D_3_4_wd$ES_34)

ND_3_4 <- En_Sw_df %>% 
  filter(depressed != "Depressed") %>%
  filter(Time %in% c("T3", "T4"))

ND_3_4 <- En_Sw_df %>% 
  filter(depressed != "Depressed") %>%
  filter(Time %in% c("T3", "T4")) %>% 
  filter(id %in% ND_3_4$id[duplicated(ND_3_4$id)])

ND_3_4_wd <- dcast(ND_3_4, id ~Time, value.var = "Cluster")
ND_3_4_wd$ES_34 <- paste(ND_3_4_wd$T3,ND_3_4_wd$T4, sep = "_")
table(ND_3_4_wd$ES_34)

###T4 to T5
D_4_5 <- En_Sw_df %>% 
  filter(depressed == "Depressed") %>%
  filter(Time %in% c("T4", "T5"))

D_4_5 <- En_Sw_df %>% 
  filter(depressed == "Depressed") %>%
  filter(Time %in% c("T4", "T5")) %>% 
  filter(id %in% D_4_5$id[duplicated(D_4_5$id)])

D_4_5_wd <- dcast(D_4_5, id ~Time, value.var = "Cluster")
D_4_5_wd$ES_45 <- paste(D_4_5_wd$T4,D_4_5_wd$T5, sep = "_")
table(D_4_5_wd$ES_45)

ND_4_5 <- En_Sw_df %>% 
  filter(depressed != "Depressed") %>%
  filter(Time %in% c("T4", "T5"))

ND_4_5 <- En_Sw_df %>% 
  filter(depressed != "Depressed") %>%
  filter(Time %in% c("T4", "T5")) %>% 
  filter(id %in% ND_4_5$id[duplicated(ND_4_5$id)])

ND_4_5_wd <- dcast(ND_4_5, id ~Time, value.var = "Cluster")
ND_4_5_wd$ES_45 <- paste(ND_4_5_wd$T4,ND_4_5_wd$T5, sep = "_")
table(ND_4_5_wd$ES_45)

###summa
table(D_1_3_wd$ES_13)
table(D_3_4_wd$ES_34)
table(D_4_5_wd$ES_45)

sum(table(ND_1_3_wd$ES_13))
sum(table(ND_3_4_wd$ES_34))
sum(table(ND_4_5_wd$ES_45))

##############Run RPART
### Prepare the data
# rpart.data <- dplyr::select(Cytokine.df, Bacteroides:Acetobacter) 
# 
# rpart.data <- sapply(rpart.data, as.numeric, simplify=T)
# rpart.data2 <- rpart.data * 100000
# rpart.data2 <- round(rpart.data2, 0)
# rownames(rpart.data2) <- Cytokine.df$Sample
# 
# 
# data(saliva)
# data(throat)
# data(tonsils)
# 
# ### Create some covariates for our data set
# site <- c(rep("Saliva", nrow(saliva)), rep("Throat", nrow(throat)), 
#           rep("Tonsils", nrow(tonsils)))
# covars <- data.frame(Group=site)
# 
# ### Combine our data into a single object
# data <- rbind(saliva, throat, tonsils)
# 
# ### For a single rpart tree
# numCV <- 0
# numCon <- 0
# rpartRes <- DM.Rpart(data, covars, numCV=numCV, numCon=numCon)
# 
# rpartRes$fullTree$where

###anova
response_variables <- c("Observed", "Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")

data=filter(richness_meta, Cluster == "2")

for (var in response_variables) {
  formula <- as.formula(paste(var, "~ depressed * days"))
  result <- aov(formula, data=data)
  # Print the response variable name
  cat("\nResults for", var, ":\n")
  print(summary(result))
}



##########################
#read Dan's model

brms.df <- read.csv("D:/OneDrive - Stanford/Manuscript/mentalhealthmicrobiome/BBI R1/BRMS model/model-summaries_microbe-cluster-depression.csv", header = T)

brms.depressed.interaction <- filter(brms.df, term == "depressedDepressed:Time.x")

table(brms.depressed.interaction$cross0)

brms.depressed.interaction.signi <- brms.depressed.interaction %>% filter(cross0 == "FALSE") %>% filter(beatsNull == "TRUE")
brms.depressed.interaction.signi.background <- brms.depressed.interaction %>% filter(cytokine=="CHEX4")
hist(brms.depressed.interaction.signi.background$Estimate)

sig.estimates <- brms.depressed.interaction.signi.background$Estimate

background_estimates <- brms.depressed.interaction.signi.background$Estimate
target_estimates <- brms.depressed.interaction.signi$Estimate

background_estimates <- na.omit(background_estimates)
target_estimates <- na.omit(target_estimates)

alculate_p_values <- function(background, targets) {
  ecdf_bg <- ecdf(background)
  p_vals <- 2 * pmin((1 - ecdf_bg(targets)), ecdf_bg(targets))
  p_vals[p_vals > 1] <- 1
  return(p_vals)
}

# Use the function to calculate p-values
results_df <- data.frame(
  Estimate = target_estimates,
  P_Value = calculate_p_values(background_estimates, target_estimates)
)

head(results_df)
results_df$P_Value_FDR <- p.adjust(results_df$P_Value, method = "fdr")

brms.combined <- cbind(brms.depressed.interaction.signi,results_df)
hist(brms.combined$P_Value_FDR)

colnames(brms.combined)[4] <- "Ori.Esstimate"

brms.combined.nochex <- dplyr::filter(brms.combined,!cytokine %in% c("CHEX1", "CHEX2", "CHEX3", "CHEX4"))
dim(brms.combined.nochex)

table(brms.combined.nochex$cytokine)
brms.combined.nochex

brms.combined.nochex <- filter(brms.combined.nochex, P_Value_FDR < 0.05)

write.csv(file = "./Result/revision/interactionterm.csv", brms.combined.nochex)

filter(brms.combined.nochex, cytokine == "GROA")                 

# Load libraries
library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(scales)

# Calculate P_Value_FDR if not already present
brms.combined.nochex <- brms.combined.nochex %>%
  mutate(P_Value_FDR = p.adjust(2 * pmin(pnorm(Ori.Esstimate, lower.tail = FALSE), pnorm(Ori.Esstimate, lower.tail = TRUE)), method = "fdr"))

# Create a combined list of cytokines and microbiomes with their types
nodes_cytokines <- brms.combined.nochex %>%
  dplyr::select(name = cytokine) %>%
  mutate(type = "Cytokine")

nodes_microbiomes <- brms.combined.nochex %>%
  dplyr::select(name = microbe) %>%
  mutate(type = "Microbiome")

# Combine and calculate frequency
nodes <- bind_rows(nodes_cytokines, nodes_microbiomes) %>%
  group_by(name, type) %>%
  summarise(freq = n(), .groups = "drop")

# Prepare edges data frame
edges <- brms.combined.nochex %>%
  dplyr::select(cytokine, microbe, Ori.Esstimate, P_Value_FDR) %>%
  rename(from = cytokine, to = microbe)

# Create igraph object
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
# Assign node sizes based on frequency
V(g)$size <- V(g)$freq

# Assign node colors based on type
V(g)$type <- V(g)$type  # Already assigned in the nodes data frame

# Assign edge attributes
E(g)$P_Value_FDR <- E(g)$P_Value_FDR
E(g)$Ori_Estimate <- E(g)$Ori.Esstimate

# Normalize P_Value_FDR for color scaling
# Lower p-values (more significant) will be darker/redder
E(g)$color <- rescale(E(g)$P_Value_FDR, to = c(0, 1))

# Assign edge width based on absolute estimate
E(g)$width <- abs(E(g)$Ori.Esstimate) * 5  # Scaling factor for better visibility

# Define a color gradient function
edge_color_gradient <- colorRampPalette(c("red", "orange", "yellow", "green", "blue"))

# Assign colors based on P_Value_FDR
E(g)$color <- edge_color_gradient(100)[as.numeric(cut(E(g)$P_Value_FDR, breaks = 100))]

# Plot the network using ggraph with Kamada-Kawai layout
network_plot <- ggraph(g, layout = "kk") +  # Kamada-Kawai layout
  # Draw edges with color based on P_Value_FDR and width based on estimate
  geom_edge_link(aes(color = P_Value_FDR, width = Ori_Estimate), alpha = 0.8) +
  
  # Draw nodes with size based on frequency and color based on type
  geom_node_point(aes(size = freq, color = type)) +
  
  # Add node labels conditionally
  geom_node_text(
    aes(label = ifelse(type == "Cytokine", name, NA)),  # Label only cytokines
    repel = TRUE, 
    size = 4, 
    na.rm = TRUE
  ) +
  
  # Scale edge colors
  scale_edge_color_gradient(low = "grey", high = "blue", name = "P_Value_FDR") +
  
  # Scale edge widths
  scale_edge_width(range = c(0.2, 1), name = "Estimate") +
  
  # Scale node sizes
  scale_size_continuous(range = c(3, 10), name = "Frequency") +
  
  # Scale node colors manually
  scale_color_manual(
    values = c("Cytokine" = "orange", "Microbiome" = "skyblue"), 
    name = "Node Type"
  ) +
  
  # Theme adjustments
  theme_minimal() +
  
  # Titles and labels
  labs(
    title = "Cytokine-Microbiome Interaction Network",
    subtitle = "Node size: Frequency | Edge color: P_Value_FDR | Edge width: Estimate"
  ) +
  
  theme(legend.position = "right")

# Display the plot
# pdf(file = "Result/revision/network.interaction.pdf", width = 6, height = 6)
# print(network_plot)
# dev.off()
