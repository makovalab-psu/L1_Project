setwd("~/Google_Drive/L1_Project/Analysis/Additional_analysis/Enzyme_site_analysis/Enzyme_site_distribution/")

#Read distance values (in bp) into R 
#distance_MSPI_TaqI <- read.table("distance_MSPI_TaqI.bed", sep="\t", header=FALSE)
# distance_MSPI_TaqI <- read.table("distance_ALL.bed", sep="\t", header=FALSE)
# dim(distance_MSPI_TaqI)
# distance <- distance_MSPI_TaqI[,13]
# head(distance)
# 
# save(distance,file='RS_Distance.RData')

#Load the distance data
load("RS_Distance.RData")

#Plot the distribution
#Histogram of coverage distribution 
#No transformation
#hist(distance, breaks = 10000, xlim = c(1,10000), ylim = c(1,800000),main = "", xlab="Distance")
#abline(v=150, col="green")
#main = "Distance Distribution Between MSPI and TaqI Sites"

#Log transformation
pdf('Distance_distribution_MSPI_TaqI.pdf')
hist(log(distance+0.0001), main = "", xlab="Log Distance")
abline(v=log(150+0.0001), col="green")
dev.off()

#Calculate percentage of distance values under 150 bp  
length(distance)
distance_less_150 <- which(distance<=150)
length(distance_less_150)
length(distance_less_150)/length(distance)
     