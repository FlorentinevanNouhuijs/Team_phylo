library(dplyr)
#packages
install.packages("reshape")
library(reshape)
library(ggplot2)
library("pheatmap")


options(dplyr.summarise.inform = FALSE)

setwd("~/Dropbox/Team-Phylo/")
clusterDF <- read.csv(file = "Output/ClusterSizeDF_Nov2022.csv", header = TRUE) # read in cluster df 
DrugsToKeep <- c('Ciprofloxacin', "Gentamicin",  "Cefuroxime", "Ceftazidime", "AmoxiClav", "PipTaz")
clusterDF <- clusterDF[clusterDF$Drug %in% DrugsToKeep, ]
clusterDF$Drug<-factor(clusterDF$Drug, levels = c("PipTaz", "Cefuroxime", "Gentamicin" , "AmoxiClav", "Ceftazidime" , "Ciprofloxacin"))
levels(clusterDF$Drug)

df = data.frame(row.names = levels(clusterDF$Drug))

for (i in 1:6){
  drug1 = levels(clusterDF$Drug)[i]
  pvaluelist<-c()
  for (j in 1:6){
    drug2 = levels(clusterDF$Drug)[j]
    print(paste(i,j))
    if (i < j){
      clusterDF_Pair = clusterDF[clusterDF$Drug %in% c(drug1,drug2),]
      #print(paste(i,j))
      #print(wilcox.test(clustersize ~ Drug, data=clusterDF_Pair)$p.value)
      pvaluelist<-c(pvaluelist,wilcox.test(clustersize ~ Drug, data=clusterDF_Pair)$p.value)
    }
    if (i >=j) pvaluelist<-c(pvaluelist,1)
  }
df[i]<-pvaluelist
}
names(df)<-levels(clusterDF$Drug)

#check to see what data looks like
View(df)

pheatmap(df, color = c("red", "yellow", "blue", "black"), breaks = c(0,0.05, 0.1, 0.99,1),
         cutree_rows = 1, treeheight_col = 0, treeheight_row = 0, fontsize = 10)








