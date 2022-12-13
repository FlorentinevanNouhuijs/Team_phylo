
#set working directory to get the clustersize df
#setwd("~/Documents/Faye_code")
setwd("~/Dropbox/Team-Phylo")

cluster_df <- read.csv("Output/ClusterSizeDF_Nov2022.csv")
#cluster_df <- read.csv("ClusterSizeDF_Nov2022.csv")
raw_data <- read.csv("Data/S1_Phylogroup.csv")
#raw_data <- read.csv("S1_Phylogroup.csv")


#loading packages
library(tidyverse)
library(ggpubr)
library(dbplyr)
library(rstatix)
library(phytools)
library(ape)
library(geiger)
library(readr)
library(dplyr)
library(ggplot2)

#viewing the data to select relevant columns
#View(raw_data)

#I need column 1, 8, 20-29, 32-42 to count the percentage that is unsusceptible
count_data <- raw_data[, c(1, 8, 20:29, 32:42)]
#View(count_data)
count_data$Phylogroup[count_data$Phylogroup=="F"]<-"D" ###MERGE D AND F

#creating avg dataframe 
avg_dataframe_drugs <- mutate(cluster_df)%>%
  group_by(Drug, phylogroup)%>%  #here I am grouping them by drug 
  summarise_at(vars(clustersize), list(avg_clade_size = mean))

#View(avg_dataframe_drugs) #check

OverviewData<-as.data.frame(avg_dataframe_drugs)
OverviewData$FractionRes = 0
for (Drug in unique(OverviewData$Drug)){
  for (phylogroup in unique(OverviewData$phylogroup)){
    SIR_List<-count_data[count_data$Phylogroup == phylogroup, Drug] 
    fractionRes = length(which(SIR_List == "R"))/length(SIR_List)
    #print(fractionRes)
    #print(length(SIR_List))
    OverviewData$FractionRes[OverviewData$phylogroup==phylogroup & OverviewData$Drug==Drug]<-fractionRes
  }
}

#reorder levels
OverviewData$Drug <- factor(OverviewData$Drug, levels = c("PipTaz", "Cefuroxime", "Gentamicin" , "AmoxiClav", "Ceftazidime" , "Ciprofloxacin"))

#trying this out with drug group 
plot <- ggplot(OverviewData,aes(x=avg_clade_size,y=FractionRes,col=Drug, shape = phylogroup))+
  labs(x = "Average cluster size", y = "Fraction resistant samples")+
  geom_point(size = 3.5, alpha = 0.8)+
  scale_x_continuous(trans='log2')+
  theme_bw()

png("Output/FractionRes_CladeSize_Fig4_nov2022.png", width = 6, height = 4, units = "in", res = 200)
plot 
dev.off()


#############################
####Linear model, now in results section. 
#############################

LM<-lm(FractionRes ~  avg_clade_size + Drug,  data = OverviewData) #ADD Drug class to the model 
summary(LM)

modelDP<-glm(clustersize ~  phylogroup + Drug,  data = cluster_df, family = "quasipoisson")
modelDtP<-glm(clustersize ~  phylogroup * Drug,  data = cluster_df, family = "quasipoisson")
modelP<-glm(clustersize ~  phylogroup,  data = cluster_df, family = "quasipoisson")
modelD<-glm(clustersize ~  Drug,  data = cluster_df, family = "quasipoisson")

summary(modelDP)
summary(modelD)
summary(modelP)

anova(modelDP, modelD, test = "F")
anova(modelDP, modelP, test = "F")
anova(modelD, modelP, test = "F")
anova(modelP, modelD, test = "F")
anova(modelDP, modelDtP, test = "F")

#############################
####Next, let's plot number of orgins / clusters vs average cluster size
###ACTUALLY, this doesn't make a lot of sense because, of course, the big subtrees are going to have more clusters (like B2)
#############################

if (FALSE){
OverviewData$NumClusters = 0
for (Drug in unique(OverviewData$Drug)){
  for (phylogroup in unique(OverviewData$phylogroup)){
    Cluster_List<-cluster_df$clustersize[cluster_df$Drug==Drug & cluster_df$phylogroup == phylogroup]
    OverviewData$NumClusters[OverviewData$phylogroup==phylogroup & OverviewData$Drug==Drug]<-length(Cluster_List)
  }
}

#trying this out with drug group 
plot <- ggplot(OverviewData,aes(x=NumClusters,y=FractionRes,col=Drug, shape = phylogroup))+
  labs(x = "Number of clusters", y = "Fraction resistant samples")+
  geom_point(size = 3.5, alpha = 0.7)+
  scale_x_continuous(trans='log2')+
  theme_bw()

png("Output/FractionRes_NumClusters_Fig4B_nov2022.png", width = 6, height = 4, units = "in", res = 200)
plot 
dev.off()
}
