#load libraries
library(adegenet); library(ape); 
library(cowplot); library(seqinr)
library(gridExtra); library(ggpubr); library(stringr); library("geiger")
library(phangorn); library("phytools"); library("lmtest"); library("treedater")
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(gridGraphics)
library(stringr)

#clear all data from environment in R 
rm(list = ls())

####load and subset Norway Data ###

setwd("~/Github/Team_Phylo")

#Read the meta data 
NorwayMetaData = read.csv(file = "Norway_data.csv", stringsAsFactors = T)

###View data to check 
View(NorwayMetaData)

# View all column names
colnames(NorwayMetaData) 

#check if there is phylogenetic data in the meta data: 
phylo_columns <- grepl("phylo", names(NorwayMetaData))
any_phylo_column <- any(phylo_columns)

if (any_phylo_column) {
  print("There is at least one column with 'phylo' in its name.")
} else {
  print("There are no columns with 'phylo' in their names.")
}
#[1] "There are no columns with 'phylo' in their names."

# check dimensions of data
dim(NorwayMetaData) #check both number. of rows and number of columns
#[1] 3254   63

# columns to keep: Drugs that are also in the Kallonnen data set and their accession numbers
# Drugs in the Kallonnen data set: PipTaz, Cefuroxime, Gentamicin, Ceftazidime, Cipro
# We only need the SIR columns: SIR: Susceptible, Intermediate, Resistance 
#NS: Non susceptible/Suscept
ColumnsToKeep_Norway <- c("Run.accession", 
                   "Ceftazidim_SIR", "Cefuroxim_SIR", "Ciprofloxacin_SIR",
                   "Gentamicin_SIR","Piperacillin_tazobactam_SIR")

#subset data 
NorwayMetaData_drugs_accession <- NorwayMetaData[, which(names(NorwayMetaData)%in%ColumnsToKeep_Norway)]
#Visual Check 
#View(NorwayMetaData_drugs_accession)

#for loops that prints the table for each of the drugs 
columns <- c("Ceftazidim_SIR", "Cefuroxim_SIR", "Ciprofloxacin_SIR", "Gentamicin_SIR", "Piperacillin_tazobactam_SIR")

# Loop through each column
for (col in columns) {
  cat(paste0(col, " Table:\n"))
  print(table(NorwayMetaData_drugs_accession[[col]]))
}


###### Now do the same for the kallonnen data ### 

#Read the meta data 
Kallonnen_MetaData = read.csv(file = "S1_Phylogroup.csv", stringsAsFactors = T)

###View data to check 
View(Kallonnen_MetaData)

# check dimensions of data
dim(Kallonnen_MetaData) #check both number. of rows and number of columns
#[1] 1509  407

# View all column names
colnames(Kallonnen_MetaData) 

#subset our drugs of interest in our dataset: PipTaz, Cefuroxime, Gentamicin, Ceftazidime, Cipro, Amoxiclav
ColumnsToKeep_Kallonnen <- c("ERR.accession", 
                          "AmoxiClav", "Ceftazidime", "Ciprofloxacin",
                          "Gentamicin","PipTaz", "Cefuroxime")
Kallonnen_drugs_accession <- Kallonnen_MetaData[, which(names(Kallonnen_MetaData)%in%ColumnsToKeep_Kallonnen)]
#Visual Check 
View(Kallonnen_drugs_accession)

#for loops that prints the table for each of the drugs 
columns_kallonnen <- c("AmoxiClav", "Ceftazidime", "Ciprofloxacin",
                          "Gentamicin","PipTaz", "Cefuroxime")

# Loop through each column
for (col in columns_kallonnen) {
  cat(paste0(col, " Table:\n"))
  print(table(Kallonnen_drugs_accession[[col]]))
}



