#### ALL TREE FUNCTIONS 
#### NEEDED FOR CountResDescendentsFromTips.R
#### NEEDED FOR plot_tree_with_intermediates_save_pdf_Nov_2022_PleuniVersion

#### MAKE RES DATA FUNCTION
make.resdata <- function(phylogroup) { # function to make a ResData df for any phylogroup
  Resfile = paste0("Data/E_coli_phylogroup_",phylogroup,".tsv") 
  ResData <- read.csv(Resfile, header = TRUE, sep="\t", stringsAsFactors = FALSE)
  if (phylogroup == "D"){
    ResfileF = paste0("Data/E_coli_phylogroup_","F",".tsv") 
    ResDataF <- read.csv(ResfileF, header = TRUE, sep="\t", stringsAsFactors = FALSE)
    ResData = rbind(ResData, ResDataF)
  }
  aacColumns = which(substr(names(ResData),1,3)=="aac")
  #ResData <- read.csv(Resfile, header = TRUE, sep="\t", stringsAsFactors = TRUE, na.strings=c("","NA"))[,c(1:42,aacColumns)]
  ResData <- ResData[,c(1:42,aacColumns)]
  names(ResData)[1]<-"Label" # rename first column to "Label"
  return(ResData)
}

#### ADD OUTGROUP FUNCTION ## PLEUNI: I AM NOT 100% SURE WHY WE NEED THIS?
### IT ADDS AN OUTGROUP FROM PHYLOGROUP D, EXPCEPT FOR PHYLOGROUP D IT ADDS AN OUTGROUP FROM A. 
add.outgroup <- function(ResData, phylogroup) {
  D_ResData <- make.resdata("D") # make df for phylogroup D to grab an outgroup row from
  A_ResData <- make.resdata("A") # make df for phylogroup A to grab an outgroup row from
  D_row <- D_ResData[55,] # take a row from phylogroup D df to make as the outgroup for the other phylogroups
  A_row <- A_ResData[90,] # take a row from phylogroup A df to make as the outgroup for phylogroup D
  if (phylogroup == "D"){ResData <- rbind(ResData, A_row)} # if phylogroup is D add outgroup from A df
  else {ResData <- rbind(ResData, D_row)} # else for other phylogroups use the outgroup from D df
  return(ResData)
}

#### FUNCTION TO SELECT RELEVANT SUBTREE AND MAKE SURE IT IS BIFURCATING
make.subtree <- function(ResData, phylogroup, CompleteTree) { # function that filters a subtree for any phylogroup
  err_tips <- ResData$Label #isolate label vector to get ERR accessions
  err_tips <- paste(err_tips) #convert elements in the vector to strings
  sub_tree <- keep.tip(CompleteTree, err_tips) #only keep the tips
  if (phylogroup == "D"){sub_tree <- root(sub_tree,  outgroup='ERR434558', resolve.root = TRUE)}
  else {sub_tree <- root(sub_tree,  outgroup='ERR435048', resolve.root = TRUE)}
  sub_tree <- multi2di(sub_tree, tol=1e-10) # needed to avoid Error in .check.tree(tree, ...) : 'tree must be bifurcating (no polytomies or unbranched nodes)'
  return(sub_tree)
}

#### FUNCTION TO SELECT RELEVANT SUBTREE AND MAKE SURE IT IS BIFURCATING
make.subtree.no_outgroup <- function(ResData, phylogroup, CompleteTree) { # function that filters a subtree for any phylogroup
  err_tips <- ResData$Label #isolate label vector to get ERR accessions
  err_tips <- paste(err_tips) #convert elements in the vector to strings
  sub_tree <- keep.tip(CompleteTree, err_tips) #only keep the tips
  #if (phylogroup == "D"){sub_tree <- root(sub_tree,  outgroup='ERR434558', resolve.root = TRUE)}
  #else {sub_tree <- root(sub_tree,  outgroup='ERR435048', resolve.root = TRUE)}
  sub_tree <- multi2di(sub_tree, tol=1e-10) # needed to avoid Error in .check.tree(tree, ...) : 'tree must be bifurcating (no polytomies or unbranched nodes)'
  return(sub_tree)
}

#### FUNCTION THAT CREATES NAMED LIST OF TIP LABELS AND PHENOTYPES (y)
make.y <- function(ResData, myTree, Drug) { # function that creates y variable
  ResistanceArray <-c()
  for (l in myTree$tip.label){
    #find right place in ResDataDetails
    r = which(ResData$Label == l)
    if (length(which(ResData$Label == l)) == 1){
      if (is.na(ResData[r, Drug])) ResistanceArray <- c(ResistanceArray, 1) ###Missing data coded as susceptible 
      if (!is.na(ResData[r, Drug])){
        if (ResData[r, Drug] == "R") ResistanceArray <- c(ResistanceArray, 2)
        if (ResData[r, Drug] == "S") ResistanceArray <- c(ResistanceArray, 1)
        if (ResData[r, Drug] == "") ResistanceArray <- c(ResistanceArray, 1) ##NOTE: for now missing data coded as susceptible. Should be removed. 
        if (ResData[r, Drug] == "I") ResistanceArray <- c(ResistanceArray, 3) ##NOTE: intermediate is coded as intermediate here!
      }
    }
  }
  y<-setNames(ResistanceArray,myTree$tip.label)
  return(y)
}

#### FUNCTION THAT CREATES THE TREE
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
} # TIP: Set prompt=TRUE above if you want to choose where the legend goes (need to click area on plot for it to appear)
