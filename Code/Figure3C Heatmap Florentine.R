library(dplyr)
library(reshape)
library("pheatmap")
library(ggplot2)
library(reshape2)
library(Hmisc)
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
View(df)

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

df2 <- melt(df, id.vars = NULL)
View(df2)

#delete empty row and empty column, Cipro and piptaz
df2 <- df[-1,]
drop <- c("Ciprofloxacin")
df2 = df2[,!(names(df2) %in% drop)]

#round digits down
df2<-round(df2,2)

#old version
heatmaply(
  df2,
  scale_fill_gradient_fun = gradient_col, Colv=NA, Rowv=NA, cellnote=df2, breaks = myBreaks
)

final_heatmap <- heatmaply(
  df2,
  colors = viridis(n = 256,  option = "magma", direction = -1), Colv=NA, Rowv=NA, cellnote = df2
)

#saving the image
png("Output/Figure_3c_Heatmap.png", width = 6, height = 4, units = "in", res = 200)
png("Output/x.png")
final_heatmap 
dev.off()

#changing the names of the columns


#create numerical values
cor(df2[sapply(df, is.numeric)])

melt(df2, na.rm = FALSE, value.name = "p-values")

colnames(df2) <- paste("Col", 1:5)
rownames(df2) <- paste("Row", 1:5)

df <- melt(df2)

x <- c("PipTaz", "Cefuroxime", "Gentamicin" , "AmoxiClav", "Ceftazidime")
y <- c("PipTaz", "Cefuroxime", "Gentamicin" , "AmoxiClav", "Ceftazidime" , "Ciprofloxacin")

ggplot(data1, aes(x = x,
                  y = y,
                  fill = value))+geom_tile()

df <- melt(df2)
colnames(df) <- c("x", "y", "value")

ggplot(df2, aes(x = x, y = y)) +
  geom_tile(color = "black") +
# geom_text(aes(label = value), color = "white", size = 4) +
  coord_fixed()

#create empty dataframe
new_df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("drug_1", "drug_2", "p_value"))

#adding rows

new_df[nrow(new_df) + 1,] <- c("Cefuroxime", "PipTaz", "0.1671587969")
new_df[nrow(new_df) + 1,] <- c("Gentamicin", "PipTaz", "0.0466210378")
new_df[nrow(new_df) + 1,] <- c("AmoxiClav", "PipTaz", "0.0043484649")
new_df[nrow(new_df) + 1,] <- c("Ceftazidime", "PipTaz", "0.0152232900")
new_df[nrow(new_df) + 1,] <- c("Ciprofloxacine", "PipTaz", "0.0002201143")
new_df[nrow(new_df) + 1,] <- c("AmoxiClav", "Cefuroxime", "0.052912085")
new_df[nrow(new_df) + 1,] <- c("Ceftazidime", "Cefuroxime", "0.178310067")
new_df[nrow(new_df) + 1,] <- c("Ciprofloxacin", "Cefuroxime", "0.003281052")
new_df[nrow(new_df) + 1,] <- c("Ciprofloxacin", "Gentamicin", "0.04765365")
new_df[nrow(new_df) + 1,] <- c("Ceftazidime", "Gentamicin", "0.56362487")
new_df[nrow(new_df) + 1,] <- c("AmoxiClav", "Gentamicin", "0.40979513")
new_df[nrow(new_df) + 1,] <- c("Ceftazidime","AmoxiClav", "1.0000000")
new_df[nrow(new_df) + 1,] <- c("Ciprofloxacin","AmoxiClav", "0.09910249")
new_df[nrow(new_df) + 1,] <- c("Ciprofloxacin", "Ceftazidime", "0.2649557")

#print
print(new_df)

#convert column to numerical 
new_df$p_value <- as.numeric(as.character(new_df$p_value))  # Convert one variable to numeric
new_df %>% mutate_if(is.numeric, ~round(., 2)) #round down to 2 digits

plot1<-ggplot(new_df, aes(drug_1, drug_2, fill= p_value)) + 
  geom_tile()
plot1


