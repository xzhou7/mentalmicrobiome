#Mediation Analysis 
#Author: Andre 
#Date Created: March 05, 2023
#Date Last Updated: August 1, 2023

library(mediation)
library(dplyr)
library(tidyverse)
library(compositions)

setwd(dir = "C:/Users/zhoux/Box/XinZhouFiles/Projects/ZXE12_Mental_Microbiome/")

#Load in data sets
cytokinesData <- read.csv("./Data/corrected_cytokines_data.csv")
psychometricData <- read.csv("./Data/Psychometric Data.csv")
genusData <- read.csv("./Data/Genus.table.csv", header = T, row.names = 1)
sampleData <- read.csv("./Data/SAMPLE_with_entero.csv", header = T, row.names = 1)

#Transforming genus data set flipping x and y axis for intersection
newGenusData <- genusData %>% t() %>% as.data.frame()

#Combining newGenusData and sampleData with matching "row names"
mergedGenusSample <- cbind(sampleData, newGenusData)

#merging cytokinesData and psychometricData by Time and id and getting rid of unnecessary columns
mergedCytokinesPsychometric <- inner_join(cytokinesData, psychometricData, by = c("Time", "id"))

#merging mergedGenusSample and mergedCytokinesPsychometric by Time and id getting rid of unnecessary columns
combinedData <- inner_join(mergedCytokinesPsychometric, mergedGenusSample, by = c("Time", "id"))

#create seperated data sets for X, M, Y
colnames(combinedData)
X <- combinedData %>% select(Bacteroides:Acetobacter)
M <- combinedData %>% select(EGF:SVCAM1)
Y <- combinedData %>% select(bdi_total.y:depressed)

#Applying CLR to Microbiome data
X_cls = X %>% 
  purrr::map(function(x) {
    x = compositions::clr(x) %>% 
      as.numeric()
    x
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

X_trim <- X_cls[, colSums(X>0) > 26]
dim(X_trim)

Y <- Y %>%
  select(-starts_with("item")) %>% select(-total_score, -psych_meds) %>% 
  mutate(depressed = ifelse(depressed == "Depressed", 1, 0))

# Testing with limited columns of X, M, and Y
# X <- X[, 1:5] #15
# X_cls <- X_cls[, 1:5]
# M <- M[, 1:2] #15
# Y <- Y[, 6:7] #5

X <- X_trim
# Create an empty data frame to store results
results <- data.frame()

# Loop over all combinations of X, M, and Y
for (x in colnames(X)) {
  for (m in colnames(M)) {
    for (y in colnames(Y)) {
      print(x)
      # Create a data frame for this combination of X, M, and Y
      data <- data.frame(X = X[[x]], M = M[[m]], Y = Y[[y]])
      
      # Exclude rows where current X,M,Y, have missing values
      data <- data[complete.cases(data$X, data$M, data$Y), ]
      
      if (nrow(data) > 0) {
      # Run the mediation analysis
      med.out <- mediation::mediate(model.m = glm(M ~ X, family = gaussian, data = data), 
                                    model.y = glm(Y ~ X + M,family = gaussian, data = data), 
                                    treat = "X", mediator = "M",
                                    boot = TRUE, sims = 100) #might want to change the number of simulations when using more powerful computer
      
      # Create a data frame for the results
      result_df <- data.frame(
        psychology = y,
        mediator = m,
        treat = x,
        ACME = med.out$d1,
        ACME_ci_lower = med.out$d1.ci[1],
        ACME_ci_upper = med.out$d1.ci[2],
        ACME_p = med.out$d1.p,
        ADE = med.out$z1,
        ADE_ci_lower = med.out$z1.ci[1],
        ADE_ci_upper = med.out$z1.ci[2],
        ADE_p = med.out$z1.p,
        total_effect = med.out$tau.coef,
        total_effect_ci_lower = med.out$tau.ci[1],
        total_effect_ci_upper = med.out$tau.ci[2],
        total_effect_p = med.out$tau.p,
        prop_mediate = med.out$n1,
        prop_mediate_ci_lower = med.out$n1.ci[1],
        prop_mediate_ci_upper = med.out$n1.ci[2],
        prop_mediate_p = med.out$n1.p
      )
      
      # Append the results to the larger data frame
      results <- rbind(results, result_df)
      }
    }
  }
}

#storing Before & After Log Transformation to compare
CLResults <- results
#Cleaning up row names for results tables
rownames(PreCLResults) <- NULL
rownames(CLResults) <- NULL

#storing results for presentation GROA & prevotella 
signifResultsCLR <- data.frame()

write.csv(results, file = "./mediation_results.csv", row.names = FALSE)
