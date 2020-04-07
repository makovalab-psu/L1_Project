setwd("~/Google_Drive/L1_Project/Analysis/Additional_analysis/Restriction_sites/")

#Read distance values (in bp) into R (not shown due to file size limit)
#Distances were calculated for any pair of adjacent sites (either MspI or TaqI)
# distance_MSPI_TaqI <- read.table("distance_ALL.bed", sep="\t", header=FALSE)
# dim(distance_MSPI_TaqI)
# distance <- distance_MSPI_TaqI[,13]
# head(distance)
# save(distance,file='RS_Distance.RData')

#Load distance data
load('RS_Distance.RData')

#Plot the distribution
#Histogram of coverage distribution 
#hist(distance, xlim=c(0,10000))
pdf("Distance_distribution_MSPI_TaqI.pdf")
hist(log(distance+0.0001), main = "Log Distance Distribution Between MspI and TaqI Sites", xlab="Log Distance")
#add cutoff of 150 bp (average read length)
abline(v=log(150+0.0001), col="green")
dev.off()

#Calculate percentage of distance values under 150 bp (average read length)
length(distance)
distance_less_150 <- which(distance<=150)
length(distance_less_150)
length(distance_less_150)/length(distance)
#=0.8746
     