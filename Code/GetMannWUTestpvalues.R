library(dplyr)
library(reshape)
library("pheatmap")
library(ggplot2)
library(reshape2)
#library(Hmisc)
library(stats)
library(RColorBrewer)


options(dplyr.summarise.inform = FALSE)
setwd("~/Dropbox/Team-Phylo/")

clusterDF <- read.csv(file = "Output/ClusterSizeDF_Nov2022.csv", header = TRUE) # read in cluster df 
DrugsToKeep <- c('Ciprofloxacin', "Gentamicin",  "Cefuroxime", "Ceftazidime", "AmoxiClav", "PipTaz")
clusterDF <- clusterDF[clusterDF$Drug %in% DrugsToKeep, ]
clusterDF$Drug<-factor(clusterDF$Drug, levels = c("PipTaz", "Cefuroxime", "Gentamicin" , "AmoxiClav", "Ceftazidime" , "Ciprofloxacin"))
levels(clusterDF$Drug)

#df = data.frame(row.names = levels(clusterDF$Drug))
#View(df)
set.seed(1)

pvaluelist<-c()
drug_1<-c()
drug_2<-c()

for (i in 1:6){
  drug1 = levels(clusterDF$Drug)[i]
  #  pvaluelist<-c()
  for (j in 1:6){
    drug2 = levels(clusterDF$Drug)[j]
    print(paste(i,j))
    if (i < j){
      clusterDF_Pair = clusterDF[clusterDF$Drug %in% c(drug1,drug2),]
      #print(paste(i,j))
      #print(wilcox.test(clustersize ~ Drug, data=clusterDF_Pair)$p.value)
      pvaluelist<-c(pvaluelist,wilcox.test(clustersize ~ Drug, data=clusterDF_Pair)$p.value)
    }
    if (i >=j) pvaluelist<-c(pvaluelist,NA)
    drug_1 <- c(drug_1, drug1) 
    drug_2 <- c(drug_2, drug2)
  }
}

df <- data.frame(drug_1, drug_2, pvaluelist)
names(df)[3] <- "p_values"

#check to see what data looks like
View(df)

df<-df[!is.na(df$p_values),]
df$drug_1<-factor(df$drug_1, levels = c("PipTaz", "Cefuroxime", "Gentamicin" , "AmoxiClav", "Ceftazidime","Ciprofloxacin" ))
df$drug_2<-factor(df$drug_2, levels = c("Ciprofloxacin", "Ceftazidime", "AmoxiClav", "Gentamicin", "Cefuroxime", "PipTaz"))

write.csv(x = df, file = "Output/MannWUtestpvalues.csv", row.names = FALSE)