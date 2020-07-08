setwd("~/Google_Drive/L1_Project/Analysis/Additional_analysis/L1_target_motif_analysis/Distance_L1_target/")

###All de novo L1s
#Read distance values (in bp) into R 
distance_denovoL1_motif_insertions <- read.table("distance_denovoAll_Insertions_targetmotifALL.bed", sep="\t", header=FALSE)
dim(distance_denovoL1_motif_insertions)
distance_insertions <- distance_denovoL1_motif_insertions[,8]
head(distance_insertions)
save(distance_insertions,file='motif_Distance_all_insertions.RData')
#Load the distance data
load("motif_Distance_all_insertions.RData")
#Histogram of coverage distribution 
#No transformation
#hist(distance_insertions, breaks = 10000, xlim = c(1,7000), ylim = c(1,600),main = "Distance distribution between all de novo L1s and \n L1 target motifs analyzed in the study",  xlab="Distance", ylab="Count")
#Plot frequency instead of count
pdf('Distance_target_all_denovoL1s.pdf',width=7,height=10)
hist(distance_insertions, breaks = 10000, xlim = c(1,7000), ylim = c(0,0.005), freq=FALSE,main = "Distance distribution between all de novo L1s and \n L1 target motifs analyzed in the study", xlab="Distance", ylab="Frequency")
abline(v=1000, col="blue")
dev.off()
#Calculate percentage of distance values under 1000  
length(distance_insertions)
distance_less_1k_insertions <- which(distance_insertions<=1000)
length(distance_less_1k_insertions)
length(distance_less_1k_insertions)/length(distance_insertions)


###Filtered de novo L1s
#Read distance values (in bp) into R 
distance_denovoL1_motif_insertions_filtered <- read.table("distance_denovoFiltered_Insertions_targetmotifALL.bed", sep="\t", header=FALSE)
dim(distance_denovoL1_motif_insertions_filtered)
distance_insertions_filtered <- distance_denovoL1_motif_insertions_filtered[,8]
head(distance_insertions_filtered)
save(distance_insertions_filtered,file='motif_Distance_all_insertions_filtered.RData')

#Load the distance data
load("motif_Distance_all_insertions_filtered.RData")
#Histogram of coverage distribution 
#No transformation
#hist(distance_insertions_filtered, breaks = 10000, xlim = c(1,7000), ylim = c(1,600),main = "Distance distribution between filtered de novo L1s and \n L1 target motifs analyzed in the study", xlab="Distance", ylab="Count")
#Plot frequency instead of count
pdf('Distance_target_filtered_denovoL1s.pdf',width=7,height=10)
hist(distance_insertions_filtered, breaks = 10000, xlim = c(1,7000),ylim = c(0,0.005),freq=FALSE,main = "Distance distribution between filtered de novo L1s and \n L1 target motifs analyzed in the study", xlab="Distance", ylab="Frequency")
abline(v=1000, col="blue")
dev.off()
#Calculate percentage of distance values under 1000  
length(distance_insertions_filtered)
distance_less_1k_insertions_filtered <- which(distance_insertions_filtered<=1000)
length(distance_less_1k_insertions_filtered)
length(distance_less_1k_insertions_filtered)/length(distance_insertions_filtered)


