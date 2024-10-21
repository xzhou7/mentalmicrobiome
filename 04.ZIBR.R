
#04.ZIBR
#Xin Zhou
#################ZIBR
library(ZIBR)
library(dplyr)
library(tidyr)

#mac
setwd("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXE12_Mental_Microbiome/")

#windows
setwd(dir = "C:/Users/zhoux/Box/XinZhouFiles/Projects/ZXE12_Mental_Microbiome/")

sample.df <- read.csv(file = "./Result/Result.table/Sample.table.csv", header = T, row.names = 1)
sample.df$Time <- 1
sample.df$Time[sample.df$Time.x == "T3"] <- 2
sample.df$Time[sample.df$Time.x == "T4"] <- 3
sample.df$Time[sample.df$Time.x == "T5"] <- 3

sample_freq <- table(sample.df$id.x) %>% data.frame()
sample.df.3 <- filter(sample.df,sample.df$id.x %in% sample_freq[sample_freq$Freq > 2,]$Var1 )
table(sample.df.3$depressed)

dplyr::select(sample.df.3, id.x, depressed) %>% unique()

sample.microbe <-dplyr::select(sample.df.3, Bacteroides:Acetobacter, Time, id.x,depressed)
ZIBR.df <- aggregate(.~Time + id.x + depressed, data=sample.microbe, FUN = mean)
table(ZIBR.df$id.x)

sample_freq2 <- table(ZIBR.df$id.x) %>% data.frame()
ZIBR.df <- filter(ZIBR.df,ZIBR.df$id.x %in% sample_freq2[sample_freq2$Freq > 2,]$Var1 )
table(ZIBR.df$id.x)
ZIBR.df$depressed_n <- 0
ZIBR.df$depressed_n[ZIBR.df$depressed=="Depressed"] <- 1
table(ZIBR.df$depressed_n)

dim(ZIBR.df)

trim.list <- data.frame(colSums(ZIBR.df[4:368]) == 0)
colnames(trim.list) <- "Col0"
trim.list <- filter(trim.list, Col0) #remove 65 taxa
ZIBR.df <- dplyr::select(ZIBR.df, -all_of(rownames(trim.list)))
dim(ZIBR.df)

zerofreq <- data.frame(colSums(ZIBR.df[4:303] == 0) < 4)
colnames(zerofreq) <- "Col0"
zerofreq <- filter(zerofreq, Col0)
ZIBR.df <- dplyr::select(ZIBR.df, -all_of(rownames(zerofreq)))
dim(ZIBR.df)

# zerofreq2 <- data.frame(colSums(ZIBR.df[4:272] != 0) < 3)
# colnames(zerofreq2) <- "Col0"
# zerofreq2 <- filter(zerofreq2, Col0)
# ZIBR.df <- dplyr::select(ZIBR.df, -all_of(rownames(zerofreq2)))
# dim(ZIBR.df)

#loop starts here
ZIBR_result <- data.frame()
for (i in 4:289){
  print(i)
  taxa <- colnames(ZIBR.df)[i]
  print(taxa)
  zibr.fit <- zibr(logistic.cov = data.frame(Treatment=ZIBR.df$depressed_n), 
                   beta.cov = data.frame(Treatment=ZIBR.df$depressed_n), 
                   Y = ZIBR.df[,i], subject.ind = ZIBR.df$id.x,
                   time.ind = ZIBR.df$Time)
  
  ZIBR_result[i,1] <- taxa
  #zero part estimation
  ZIBR_result[i,2] <- zibr.fit$logistic.est.table[2]
  ZIBR_result[i,3] <- zibr.fit$logistic.est.table[4]
  #none-zero part
  ZIBR_result[i,4] <- zibr.fit$beta.est.table[2]
  ZIBR_result[i,5] <- zibr.fit$beta.est.table[4]
  ZIBR_result[i,6] <- zibr.fit$joint.p}

colnames(ZIBR_result) <- c("Genus", "Logistic_Beta", "Logistic_p", "Beta_Beta", "Beta_p", "Joint_p")

ZIBR_result$p.adj <- p.adjust(ZIBR_result$Joint_p, method = "BH", n=length(ZIBR_result$Joint_p))

ZIBR_result %>% filter(p.adj < 0.15)

#write.csv(file = "./Result/Result.table/ZIBR.csv",ZIBR_result)

#plot single value (if necessary)
zibr.fit <- zibr(logistic.cov = data.frame(Treatment=ZIBR.df$depressed_n), 
                 beta.cov = data.frame(Treatment=ZIBR.df$depressed_n), 
                 Y = ZIBR.df$Leuconostoc, subject.ind = ZIBR.df$id.x,
                 time.ind = ZIBR.df$Time)

zibr.fit

#######################
#check single bacteria 
ZIBR.df %>% filter(depressed != "Depressed") %>% dplyr::select(Unclassified_Bacillales)

table(ZIBR.df$depressed)

# Function to compute non-zero mean
non_zero_mean <- function(x) {
  non_zero_values <- x[x != 0]
  return(mean(non_zero_values, na.rm = TRUE))
}

# Function to compute non-zero standard deviation
non_zero_sd <- function(x) {
  non_zero_values <- x[x != 0]
  return(sd(non_zero_values, na.rm = TRUE))
}

# Assuming the dataframe is named ZIBR.df
# Extract bacteria columns and calculate non-zero mean and SD
result <- ZIBR.df %>%
  dplyr::select(Prevotella:Acetobacter) %>%
  gather(key = "Bacteria", value = "Value") %>%
  group_by(Bacteria) %>%
  summarise(mean = non_zero_mean(Value),
            sd = non_zero_sd(Value))

result <- ZIBR.df %>%
  dplyr::select(depressed, Prevotella:Acetobacter) %>%  # Include the 'depressed' column
  gather(key = "Bacteria", value = "Value", -depressed) %>%  # Don't gather 'depressed'
  group_by(depressed, Bacteria) %>%  # Group by both 'depressed' and 'Bacteria'
  summarise(mean = non_zero_mean(Value),
            sd = non_zero_sd(Value))

de_result <- result %>% filter(depressed  == "Depressed") %>% ungroup() %>% dplyr::select(-depressed)
colnames(de_result) <- c("Genus",  "mean_dp", "sd_dp") 

nd_result <- result %>% filter(depressed  != "Depressed") %>% ungroup() %>% dplyr::select(-depressed)
colnames(nd_result) <- c("Genus",  "mean_nd", "sd_nd") 

ZIBR_result_meansd <- left_join(de_result,nd_result, by = c("Genus" = "Genus" )) %>% left_join(ZIBR_result,  by = c("Genus" = "Genus" ))
dim(ZIBR_result_meansd)

ZIBR_result_meansd %>% filter(p.adj < 0.15)

write.csv(file = "./Result/Result.table/ZIBR_with_mean.csv",ZIBR_result_meansd)

