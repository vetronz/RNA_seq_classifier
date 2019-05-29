getwd()


setwd("C:\\Users\\Myrsini\\Desktop\\TWG\\mega")


#To install core packages
source("http://bioconductor.org/biocLite.R")



#install lumi package
biocLite("lumi")

#load the lumi package
library("lumi")


#read in clinical data
targets<-read.delim(file="sample_details.txt", header=T)



head(targets)

tail(targets)

dim(targets)



category<-as.factor(targets$category)

#read in array data
filename<-"Mega Expt all arrays 21_03_14 n=656_Sample_Probe_Profile_GX format 02-03-15 No BG subtraction.txt"
object<-lumiR(filename, detectionTh = 0.01, 
              checkDupId = TRUE, QC = TRUE, 
              columnNameGrepPattern = list('AVG_SIGNAL', se.exprs='BEAD_STD', detection='Detection', beadNum='Avg_NBEADS'), 
              verbose = TRUE)


#get the control data
controld<-getControlData("Control probe profile.txt", type = 'LumiBatch')
plotControlData (getControlData("Control probe profile.txt")[,1:100])



#add the control data
object.c<-addControlData2lumi("Control probe profile.txt",object)

rm(object)

#structure of the object
str(object.c)



#CLASS

#qc plots
#subsetting for genes and samples
example.lumi<-object.c[1:1000,1:50]


plot(example.lumi, what='density',legend="FALSE")
## plot the MAplot M (log ratio) and A (mean average)

example.lumi<-object[1:1000,1:5]
plot(example.lumi, what='MAplot')

#the CDF plot in reverse direction can better show the different 
#at the high and middle expression range among different samples
plotCDF(example.lumi, reverse=TRUE)

## density plot of coefficient of varience
#The coefficient of variation (CV) is defined as the ratio of the standard deviation to the mean
plot(example.lumi, what='cv',legend="FALSE")


##Eucledian clustering 
plot(example.lumi, what='sampleRelation')

plot(example.lumi, what='sampleRelation', method='mds', color=c("01", "02", "01", "02"))

plot(example.lumi, what='boxplot')


#plot boxplot before normalisation
pdf(file="boxplot_before_normalisation_nobck.pdf", height=8, width=90)
boxplot(object.c,cex=0.7)
dev.off()

#background adjustment
pdf(file="boxplot_before_normalisation_bck_bgadjust.pdf", height=8, width=90)
boxplot(lumiB(object.c,"bgAdjust"),cex=0.7)
dev.off()

#try force positive

##NORMALISATION WITH METHOD
lumi.x.N<-lumiN(lumiQ(lumiT(lumiB(object.c,"bgAdjust"),
                            method='log')),method= "quantile") 
#quantile, 
#SSN (Simple Scaling Normalization), 
#RSN (Robust Spline Normalization), 
#loess normalization  
#Rank Invariant Normalization