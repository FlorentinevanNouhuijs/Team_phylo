# CODE BELOW CREATES FIGURE 5 PHYLOCLUST HEATMAP
setwd("~/Documents/GitHub/Team_phylo")

# import necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

phyloclustDF <- read.csv(file = 'Output/phyloclust_pvalues_Dec2022.csv', header = TRUE) # read in phyloclust df

# for labeling pvalues on heatmap
phyloclustDF$pval_label <- round(phyloclustDF$pval ,digit=3) # Round off the pval column to 3 decimals
phyloclustDF$pval_label[phyloclustDF$pval == 0] <- "<10e-5" # for p val of 0 show it on the heatmap as "<10e-5"

phyloclustDF$pval_label[is.na(phyloclustDF$pval)] <- "NA" # for p val of 0 show it on the heatmap as "<10e-5"

phyloclustDF$stars[phyloclustDF$pval < 0.05] <- "*" # create star column where p values < 0.05 get a star



# create a heatmap 
heatmap <- ggplot(phyloclustDF, aes(Drug, Phylogroup)) +                           
  geom_tile(aes(fill = pval)) +                                         
  #geom_text(aes(label = stars), color = "white") + # uncomment to show just asterisks for p-values < 0.05
  geom_text(aes(label = pval_label), color = "white") + # uncomment to show actual p-values
  scale_x_discrete(limits = c("PipTaz", "Cefuroxime", "Gentamicin", "AmoxiClav", "Ceftazidime", "Ciprofloxacin")) +
  theme_bw() + # remove gray background
  theme( # remove grid
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) 
heatmap

# saving heatmap as pdf
#pdf(file = "Output/Fig5_phyloclust_heatmap_Dec2022.pdf")
png("Output/Fig5_phyloclust_heatmap_Dec2022.png", width = 6, height = 5, units = "in", res = 200)
heatmap
dev.off()
