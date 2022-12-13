## read tree from string
library(ape)
text.string<-
  "(
  (
  ((S, S),S)
,(((R,R),(R,R))
,(S,S))
),S);"
vert.tree<-read.tree(text=text.string)
plot(vert.tree,no.margin=TRUE,edge.width=2)

library(ggplot2)
library(phytools)

png (filename = "Output/simpletree_cluster.png", width = 50, height = 50, units='mm', res = 300)
plotTree(vert.tree,type="phylogram",fsize=0.5,ftype="i") # uncomment if you want fans style tree
#plotTree(myTree,fsize=0.5,ftype="off")
tips_palette <- c("black", "red", "pink") # palette for the tip labels
mycols<-tips_palette[c(1,1,1,2,2,2,2,1,1,1)] # match palette colors to phenotypes
print(mycols)
tiplabels(pch = 16,col=mycols,cex=1)
dev.off()


text.string_2<-  "((((S, R),S),(((S,R),(S,S)),(R,S))),R);"
vert.tree_2<-read.tree(text=text.string_2)
png (filename = "Output/simpletree_singletons.png", width = 50, height = 50, units='mm', res = 300)
plotTree(vert.tree_2,type="phylogram",fsize=0.5,ftype="i") # uncomment if you want fans style tree
#plotTree(myTree,fsize=0.5,ftype="off")
tips_palette <- c("black", "red", "pink") # palette for the tip labels
mycols<-tips_palette[c(1,2,1,1,2,1,1,2,1,2)] # match palette colors to phenotypes
print(mycols)
tiplabels(pch = 16,col=mycols,cex=1)
dev.off()
