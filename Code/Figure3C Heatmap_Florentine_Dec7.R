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
    if (i >=j) pvaluelist<-c(pvaluelist,NA)
  }
  df[i]<-pvaluelist
}
names(df)<-levels(clusterDF$Drug)

#check to see what data looks like
View(df)

#df2 <- melt(df, id.vars = NULL)
#View(df2)

#delete empty row and empty column, Cipro and piptaz
df2 <- df[-1,]
drop <- c("Ciprofloxacin")
df2 = df2[,!(names(df2) %in% drop)]

#round digits down
df2<-round(df2,2)

#old version
#heatmaply(
#  df2,
#  scale_fill_gradient_fun = gradient_col, Colv=NA, Rowv=NA, cellnote=df2, breaks = myBreaks
#)

final_heatmap <- heatmaply(
  df2,
  colors = viridis(n = 256,  option = "magma", direction = -1), Colv=NA, Rowv=NA, cellnote = df2
)

#saving the image
png("Output/Figure_3c_Heatmap.png", width = 6, height = 4, units = "in", res = 200)
png("Output/x.png")
final_heatmap 
dev.off()


#plotly::orca(p = final_heatmap, #the graph to export
 #              file = "Output/graph 1.png") #the name and type of file (can be .png, .jpeg, etc.)
#save_image(p = final_heatmap, file = "Output/graph 1.png")


