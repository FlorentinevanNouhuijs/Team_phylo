#load packages 
library(dplyr)
library(reshape)
library("pheatmap")
library(ggplot2)
library(reshape2)
#library(Hmisc)
library(stats)
library(RColorBrewer)

setwd("~/Documents/GitHub/Team_phylo")

df <- read.csv(file = "Data/MannWUtestpvalues.csv")
df$drug_1<-factor(df$drug_1, levels = c("PipTaz", "Cefuroxime", "Gentamicin" , "AmoxiClav", "Ceftazidime","Ciprofloxacin" ))
df$drug_2<-factor(df$drug_2, levels = c("Ciprofloxacin", "Ceftazidime", "AmoxiClav", "Gentamicin", "Cefuroxime", "PipTaz"))

df$p_valuelabels <- round(df$p_values,3)
df$p_valuelabels[df$p_valuelabels<0.001] <- "<0.001"


# create a heatmap 
heatmap <- ggplot(df, aes(drug_1, drug_2)) +                           
  geom_tile(aes(fill = p_values)) +                                         
  geom_text(aes(label = p_valuelabels), color = "white") + # uncomment to have all actual p values shown
  scale_x_discrete(limits = c("PipTaz", "Cefuroxime", "Gentamicin", "AmoxiClav", "Ceftazidime")) +
  theme_bw() + # remove gray background
  labs(y = "Drugs", x = "Drugs") +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.border = element_blank()
  ) 
heatmap
#saving the image
png("Output/Figure_4c_Heatmap.png", width = 6, height = 5, units = "in", res = 200)
heatmap
dev.off()

