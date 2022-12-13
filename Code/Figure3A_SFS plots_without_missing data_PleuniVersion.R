
setwd("~/Dropbox/Team-Phylo/")
#setwd("~/Documents/Faye_Code/")

clusterDF <- read.csv(file = "Output/ClusterSizeDF_Nov2022.csv", header = TRUE) # read in cluster df 
#View(clusterDF) # view df in new tab
clusterDF$clustersize[clusterDF$clustersize > 10 ] <- 10 # sets anything over 10 as just 10 (so that we can have them all grouped as 10+)
clusterDF$Drug<-factor(clusterDF$Drug, levels = c("PipTaz", "Cefuroxime", "Gentamicin" , "AmoxiClav", "Ceftazidime" , "Ciprofloxacin"))


# check unique drugs
unique(clusterDF$Drug)

# Function that creates a new column of the normalized cluster counts per drug
make.countDF <- function(df) { 
  df <- df %>% # overwrite df with the new one we make from code below
    group_by(Drug) %>% # get cluster size counts for each drug
    count(clustersize) # makes a new column called n (number of observations)
  #normalize y axis (n)
  df <- df %>% # overwrite df with the new one we make from code below
    group_by(Drug) %>% 
    mutate(normalized_counts = round(n / sum(n), 3)) # makes a new column by calculating normalized counts per drug 
  return(df)
}

#creating data frames for the drugs
df_cipro <- (filter(clusterDF, Drug =='Ciprofloxacin'))
df_genta <- (filter(clusterDF, Drug == "Gentamicin"))
df_cefuroxime <- (filter(clusterDF, Drug == "Cefuroxime"))
df_ceftazidime <- (filter(clusterDF, Drug == "Ceftazidime"))
df_amoxc <- (filter(clusterDF, Drug == "AmoxiClav"))
df_pip <- (filter(clusterDF, Drug == "PipTaz"))

#binding the dataframes together
df_total_data <- rbind(df_amoxc, df_pip, df_cipro, df_genta, df_cefuroxime, df_ceftazidime) 
df_total_data <- make.countDF(df_total_data) # call function to calculate normalized cluster counts 

View(df_total_data)

#reorder levels
levels(df_total_data$Drug) <- c("PipTaz", "Cefuroxime", "Gentamicin" , "AmoxiClav", "Ceftazidime" , "Ciprofloxacin")

# SFS PLOT FOR DRUGS
# Function that creates a plot for specified df
make.SFS <- function(df, title) { 
  p <- ggplot(data=df, aes(x = clustersize, y = normalized_counts, fill=Drug)) + # fill controls bar color by the drug groups
    geom_bar(stat="identity", color="black", position=position_dodge2(preserve = "single")) + # position dodge allows different color bars to be beside each other rather than stacked on top
    coord_cartesian(xlim = c(1, 10)) +
    scale_x_continuous(breaks=c(seq(1, 10, by=1)),labels=c("1","2","3","4","5","6","7","8","9","10+")) +
    #ggtitle(title) +
    labs(y = "Cluster size counts") 
    p + theme_bw()
    #theme(plot.title = element_text(hjust = 0.5))
}

png("Output/SFS_Fig3B_nov2022.png", width = 6, height = 4.5, units = "in", res = 200)
make.SFS(df_total_data, title = "Clustersize counts") 
dev.off()

# To do list:
#normalize y axis - DONE
#last column is 10 or more  - DONE
# try this with drugs - DONE

# Aug 23 To do list
# get total number of resistant strains
# for each drug calc percent of resistance
# how does cluster size correlate with total num of res strains
# do common strains have large cluster sizes?
# see if transmission or de novo evolution is higher for a drug
# create table with overview of dataset

df_totalNumClusters<-df_total_data %>% 
  group_by(Drug) %>% 
  summarise(sum = sum(n))

PlotNumClusters <-ggplot(data=df_totalNumClusters, aes(x = Drug, y = sum, fill=Drug)) + # fill controls bar color by the drug groups
  geom_bar(stat="identity", color="black", position=position_dodge()) + # position dodge allows different color bars to be beside each other rather than stacked on top
  labs(y = "Estimated number of clusters") +
  geom_text(aes(label = sum), vjust = 2) +
  theme_bw()+
  theme(legend.position = "none")
  
png("Output/Figure3APlotNumClusters.png", width = 6, height = 4, units = "in", res = 200)
PlotNumClusters
dev.off()
