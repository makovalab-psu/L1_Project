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
#Violin plot of coverage distribution 
#install.packages("vioplot")
library(vioplot)
vioplot(log(distance+0.0001), names=c("Log Distance Distribution Between MSPI and TaqI Sites"), ylab="Log Distance")
abline(h=log(150+0.0001), col="green")
#ylim = c(0,0.02))

#Histogram of coverage distribution 
#hist(distance, xlim=c(0,10000))
hist(log(distance+0.0001), main = "Log Distance Distribution Between MSPI and TaqI Sites", xlab="Log Distance")
#add cutoff of 150 bp (average read length)
abline(v=log(150+0.0001), col="green")


#Calculate percentage of distance values under 150 bp (average read length)
length(distance)
distance_less_150 <- which(distance<=150)
length(distance_less_150)
length(distance_less_150)/length(distance)
#=0.8746
     