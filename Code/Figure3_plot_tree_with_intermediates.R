#### Code that plots a tree for each Phylogroup 
#### Tree includes resistant isolates, intermediate isolates and susceptible isolates
#### This code uses the TreeFunctions.R code as a source

# remove all objects from Global environment
#### WARNING ONLY RUN THIS CODE WHEN YOU HAVE SAVED ALL PREVIOUS CODE
rm(list = ls())

# load packages
library(phytools)
library(ape)
library(ggplot2)
library(dplyr)
library("treeio")
library("ggtree")

#outgroup a: 'ERR435048'
#outgroup d: 'ERR434558'

#setwd to github working directory 
setwd("~/Documents/GitHub/Team_phylo")

# FOR CONTREE FILE UNCOMMENT LINES BELOW
#treefile = paste0("Data/cleaned_core_alignments_copy.contree") # read tree from contree file
CompleteTree <- ape::read.tree("Data/cleaned_core_alignments.contree") # load tree
DrugList = c("Cefuroxime", "Gentamicin", "AmoxiClav", "Ceftazidime", "Ciprofloxacin", "PipTaz")
PhyloGroupList = c("A", "B1", "B2", "D") ###NOV 2022 removed F. Will be merged with D automatically. 

source("Code/TreeFunctions.R") #using our script that has all the tree functions

### RUN ALL THE FUNCTIONS FOR ALL DRUGS AND ALL PHYLOGROUPS
if (FALSE){
  for (Phylogroup in PhyloGroupList){
    pdf(file = paste0("Output/Trees_PhyloGroup",Phylogroup , ".pdf"), width = 10, height = 15)
    par(mfrow=c(3,2))
    #par(mfrow=c(1,1))
    for (Drug in DrugList){
      ResData <- make.resdata(Phylogroup) # call function to create phylogroup A ResData
      #ResData <- add.outgroup(ResData,Phylogroup) # can it work without outgroup?
      myTree <- make.subtree.no_outgroup(ResData,Phylogroup, CompleteTree)
      y_list <- make.y(ResData, myTree, Drug)
      plot.tree(myTree, y_list, Phylogroup, Drug)
    }
    
    dev.off()
  }
}

### RUN ALL THE FUNCTIONS FOR ALL DRUGS AND ALL PHYLOGROUPS
### DIFFERENT TREE FUNCTION
#for (Phylogroup in PhyloGroupList){ # "D"
Phylogroup = "B2"; Drug = "Ciprofloxacin"
#for (Drug in DrugList){  ##c("Gentamicin", "Ciprofloxacin")){
    wh = 8; if (Phylogroup =="B2") wh = 10;
    cx = 0.2; if (Phylogroup =="B2") cx = 0.1;
    fs = 0.5; if (Phylogroup =="B2") fs = 0.2;
    png(file = paste0("Output/Trees/AllTrees_PhyloGroup",Phylogroup, Drug, ".png"), width = wh, height = wh, units = "in", res = 300)
    par(mfrow=c(1,1))
    ResData <- make.resdata(Phylogroup) # call function to create phylogroup A ResData
    #ResData <- add.outgroup(ResData,Phylogroup) # can it work without outgroup?
    myTree <- make.subtree.no_outgroup(ResData,Phylogroup, CompleteTree)
    y_list <- make.y(ResData, myTree, Drug)
    y_list<-as.factor(y_list)
    myTree<-as.phylo(myTree)
    myTree$tip.label<-paste0("    ", myTree$tip.label, "    ")
    cols<-c("white", "red", "pink")
    ## then plot it:
    plotTree(myTree,type = "fan", ftype="i",fsize=fs,color="darkgrey", offset=.5)
    myTree$tip.label<-substr(myTree$tip.label, 5,13)
    tiplabels(pie=to.matrix(y_list[myTree$tip.label],
                            levels(y_list)),piecol=cols,cex=cx)
    legend(x="topright",c("susceptible", "intermediate", "resistant"),pch=22,
           pt.bg=cols[c(1,3,2)],pt.cex=1.5,bty="n",cex=0.9)
    if (Phylogroup != "D") title(main=paste0("E. coli, Phylogroup ", Phylogroup, ", resistance to ", Drug),font.main=3,line=-2)
    if (Phylogroup == "D")title(main=paste0("E. coli, Phylogroup D and F, resistance to ", Drug),font.main=3,line=-2)
    dev.off()
  #}
#}

#### for loop that calculates how much data is missing for each drug and each phylogroup? 
    
DrugList = c("Cefuroxime", "Gentamicin", "AmoxiClav", "Ceftazidime", "Ciprofloxacin", "PipTaz")
PhyloGroupList = c("A", "B1", "B2", "D") 
    
# Outer loop for phylogroups
for (phylogroup in PhyloGroupList) {
  cat("Phylogroup:", phylogroup, "\n")
      
  # Inner loop for drugs
  for (drug in DrugList) {
    cat("\tDrug:", drug, "\n")
        
    # Calculate missing data and intermediates
    num_missing <- sum(is.na(ResData[ResData$Phylogroup == phylogroup, drug]) | ResData[ResData$Phylogroup == phylogroup, drug] == "")
    num_intermediate <- sum(ResData[ResData$Phylogroup == phylogroup, drug] == "I")
        
    # Print results
    cat("\tNumber of missing data:", num_missing, "\n")
    cat("\tNumber of intermediate:", num_intermediate, "\n")
      }
    }

### old code 
print (paste0(Phylogroup, " ", Drug))
print(paste0("number missing data: ", length(which(is.na(ResData[,Drug]) | ResData[,Drug] == ""))))
print(paste0("number intermediate: ",length(which(ResData[,Drug] == "I"))))



