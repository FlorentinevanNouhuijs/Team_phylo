### Code that runs phyloclust function and extracts p-values
# import necessary libraries
library(phytools)
library(ape)
library(RRphylo) # library that contains phyloclust function
library(tidyr)
library(dplyr)

setwd("~/Documents/GitHub/Team_phylo")

source("Scripts/Script_for_paper/TreeFunctions.R") # importing functions from a separate R script
treefile = paste0("Data/cleaned_core_alignments.contree") # read tree from the contree file
myTree <- ape::read.tree(treefile) # load tree

phylogroup_list = c("A","B1","B2","D")
drug_list = c("PipTaz", "Cefuroxime", "Gentamicin", "AmoxiClav", "Ceftazidime", "Ciprofloxacin")

# empty lists for df columns
phylogroup_col <- c() 
drug_col <- c()
pval_col <- c()

# for loop that will run the phyloclust function for each phylogroup + drug combination
for (phylogroup in phylogroup_list) {
  print(phylogroup)
  ResData <- make.resdata(phylogroup) # set up ResData df each time loop runs
  ResData <- add.outgroup(ResData, phylogroup) # add a row at the bottom of current phylogroup df from another phylogroup's row, we will use this as our outgroup
  subtree <- make.subtree(ResData, phylogroup, myTree) # set up phylogroup tree each time loop runs
  for (Drug in drug_list) {
    y <- make.y(ResData,subtree,Drug) # get y
    if (sum(y==2)>1) { # only runs the phyloclust if the y has more than 1 resistant phenotypes
      output <- phyloclust(subtree,y,focal=2,nsim=100000) # get output of phyloclust
      p_val <- output[[1]] # extract p-value from output
      print(paste0("Drug: ",Drug,", p-val: ",p_val))
      # saving values 
      phylogroup_col <- c(phylogroup_col, phylogroup) 
      drug_col <- c(drug_col, Drug) 
      pval_col <- c(pval_col, p_val)
    }
    else { # set p-value to NA if the drug/phylogroup has 1 or less resistant phenotypes (phyloclust will not run if this is the case)
      p_val <- "NA" 
      print(paste0("Drug: ",Drug,", p-val: ",p_val))
      # saving values 
      phylogroup_col <- c(phylogroup_col, phylogroup) 
      drug_col <- c(drug_col, Drug) 
      pval_col <- c(pval_col, p_val)}
  } 
}  

# Convert lists into dataframe columns
df <- data.frame(unlist(phylogroup_col), unlist(drug_col), unlist(pval_col))
# Names of columns of dataframe
names(df) <- c("Phylogroup", "Drug", "pval")
# saving df as a csv 
write.csv(x = df, file = "Output/phyloclust_pvalues_Dec2022.csv", row.names = FALSE)


plot.tree <- function(myTree,y_list, Phylogroup, Drug) { # function for plotting trees
  #if (sum(y==3)>=1){ # checks if there is an intermediate in the y variable, if TRUE then plot with pink
  plotTree(myTree,type="phylogram",fsize=0.35,ftype="i", branch.length = "none")# uncomment if you want fans style tree
  tips_palette <- c("white", "red", "pink") # palette for the tip labels
  mycols<-setNames(tips_palette[1:length(unique(y_list))],sort(unique(y_list))) # match palette colors to phenotypes
  #print(mycols)
  tiplabels(pie=to.matrix(y_list,sort(unique(y_list))),piecol=tips_palette,cex=0.25)
  #legend_palette <- c("white", "pink", "red") # palette for the legend column 
  #my_pheno <- c("S", "I", "R") # phenotypes for the legend column
  #names(legend_palette) <-my_pheno # assign palette colors to the respective phenotypes
  #legend_col <- legend_palette # set this to a new variable name
  title(main=paste0(Phylogroup," ", Drug),font.main=3,
        line=-2)
} 
### DEBUGGING (The code below looks at why phyloclust will not run if there is less than 2 resistant phenotypes)
# Example: Ceftazidime + B1 
ResData_B1 <- make.resdata("B1") 
ResData_B1 <- add.outgroup(ResData_B1, "B1") 
subtree_B1 <- make.subtree(ResData_B1, "B1", myTree) 
CTZ_y <- make.y(ResData_B1,subtree_B1,"Ceftazidime")
phyloclust(subtree_B1,CTZ_y,focal=2,nsim=100000) # Error here
plot.tree(subtree_B1, CTZ_y, "B1", "Ceftazidime") # check number of resistant phenotypes

# check where I get the error (phyloclust source code found here: https://rdrr.io/cran/RRphylo/src/R/phyloclust.R)
state <- treedataMatch(subtree_B1, CTZ_y)[[1]][,1]
2->st
cophenetic.phylo(subtree_B1)->cop
cop[which(state==st),which(state==st)]->subcop
mean(subcop[upper.tri(subcop)])->mds
dim(subcop)[1]->sl

r.mds<-array()

for(e in 1:100000){
  sample(subtree_B1$tip.label,sl)->test.tip
  cop[test.tip,test.tip]->r.cop
  mean(r.cop[upper.tri(r.cop)])->r.mds[e]
}
pval=length(which(r.mds<mds))/nsim

# check each line in the for loop starting with e = 1
e <- 1
sample(subtree_B1$tip.label,sl)->test.tip # ERROR COMES FROM THIS
cop[test.tip,test.tip]->r.cop
mean(r.cop[upper.tri(r.cop)])->r.mds[e]

# check variables
state
cop
subcop # 0
mds # NA value
sl # null value

# check why I dont get error for PipTaz
PTZ_y <- make.y(ResData_B1,subtree_B1,"PipTaz")
state <- treedataMatch(subtree_B1, PTZ_y)[[1]][,1]
2->st
cophenetic.phylo(subtree_B1)->cop
cop
cop[which(state==st),which(state==st)]->subcop
subcop # returns all of the resistant accessions and their cophenetic distances from each other
mean(subcop[upper.tri(subcop)])->mds
mds
dim(subcop)[1]->sl
sl # 3 resistant phenotypes

# Conclusion: we get an error if there aren't at least 2 resistant samples in a phylogroup + drug combo
# This is because phyloclust needs at least 2 samples to calculate their cophenetic distances



