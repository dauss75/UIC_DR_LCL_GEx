library(lumi); library(gtools); library(limma); library(bioDist); library(calibrate)

################################
##  directory setup
################################
HOME <- "/Users/sjung/Project/GlobusGenomics/UIC/310/"
InputDir <- paste(HOME,'InputData/', sep='')
InputFiles = list.files(path=InputDir, pattern='txt')
ControlDir <- paste(HOME,'ControlData/', sep='')
ControlFiles = list.files(path=ControlDir, pattern='txt')
PhenotypeDir = paste(HOME,'PhenotypeInfo/', sep='')
Output=paste(HOME,'output/',sep=''); dir.create(Output, showWarnings = FALSE)
Figs = paste(HOME,'Figs/', sep=''); dir.create(Figs, showWarnings = FALSE)

load(file = paste(Output,"sortedData.RData",sep=""))
geneList<-read.csv(file="DeltaAnalysisResults6.15.16.csv")
colnames(geneList)[4]<-"pval.equal.var"
colnames(geneList)[5]<-"pval.treat.eff"
geneList_p1<-geneList[geneList$pval.equal.var<0.01,]
data_p1<-sortedData[geneList_p1$geneID,]

library(corrplot)

cor_plot <- function(a,b, data) {
  M <- cor(data[,a:b])
  png(paste(Figs,a,"_",b,"_cor_p1.png", sep=""))
  corrplot.mixed(M, lower="number", upper="circle")
  dev.off()
}

cor_plot(1,5,sortedData)
cor_plot(11,13,sortedData);cor_plot(17,19,sortedData);cor_plot(23,25,sortedData)
cor_plot(29,31,sortedData);cor_plot(35,37,sortedData);cor_plot(41,43,sortedData); 
cor_plot(47,49,sortedData); cor_plot(53,57,sortedData); cor_plot(63,65,sortedData);
cor_plot(69,71,sortedData); cor_plot(75,77,sortedData); cor_plot(81,83,sortedData); 
cor_plot(87,89,sortedData); cor_plot(93,97,sortedData); cor_plot(103,105,sortedData);
cor_plot(109,111,sortedData); cor_plot(115,117,sortedData); cor_plot(121,123,sortedData); 
cor_plot(127,129,sortedData); cor_plot(133,135,sortedData); cor_plot(139,141,sortedData)

cor_plot(6,10,sortedData)
cor_plot(14,16,sortedData);cor_plot(20,22,sortedData);cor_plot(26,28,sortedData)
cor_plot(32,34,sortedData);cor_plot(38,40,sortedData);cor_plot(44,46,sortedData); 
cor_plot(50,52,sortedData); cor_plot(58,62,sortedData); cor_plot(66,68,sortedData);
cor_plot(72,74,sortedData); cor_plot(78,80,sortedData); cor_plot(84,86,sortedData); 
cor_plot(90,92,sortedData); cor_plot(98,102,sortedData); cor_plot(106,108,sortedData);
cor_plot(112,114,sortedData); cor_plot(118,120,sortedData); cor_plot(124,126,sortedData); 
cor_plot(130,132,sortedData); cor_plot(136,138,sortedData); cor_plot(142,144,sortedData)

sortedData.T<-data.frame(t(sortedData))
sortedData.T$Name<-substr(colnames(sortedData),1,nchar(colnames(sortedData))-2)
meansortedData.T<-aggregate(sortedData.T[,-ncol(sortedData.T)],list(sortedData.T$Name),mean)
rownames(meansortedData.T)<-meansortedData.T$Group.1
meansortedData.T$Group.1<-NULL
meansortedData<-data.frame(t(meansortedData.T))

cor_plot(1,44,meansortedData)
cor_plot(1,11,meansortedData)
cor_plot(12,22,meansortedData)
cor_plot(23,33,meansortedData)
cor_plot(34,44,meansortedData)

# library("PerformanceAnalytics")
# 
# cor_plot2 <- function(a,b, data) {
#   png(paste(Figs,a,"_",b,"_cor2.png", sep=""))
#   chart.Correlation(data, histogram=FALSE, pch="18")
#   dev.off()
# }
# 
# cor_plot2(1,5,sortedData)
# cor_plot2(11,13,sortedData);cor_plot2(17,19,sortedData);cor_plot2(23,25,sortedData)
# cor_plot2(29,31,sortedData);cor_plot2(35,37,sortedData);cor_plot2(41,43,sortedData); 
# cor_plot2(47,49,sortedData); cor_plot2(53,57,sortedData); cor_plot2(63,65,sortedData);
# cor_plot2(69,71,sortedData); cor_plot2(75,77,sortedData); cor_plot2(81,83,sortedData); 
# cor_plot2(87,89,sortedData); cor_plot2(93,97,sortedData); cor_plot2(103,105,sortedData);
# cor_plot2(109,111,sortedData); cor_plot2(115,117,sortedData); cor_plot2(121,123,sortedData); 
# cor_plot2(127,129,sortedData); cor_plot2(133,135,sortedData); cor_plot2(139,141,sortedData)
# 
# cor_plot2(6,10,sortedData)
# cor_plot2(14,16,sortedData);cor_plot2(20,22,sortedData);cor_plot2(26,28,sortedData)
# cor_plot2(32,34,sortedData);cor_plot2(38,40,sortedData);cor_plot2(44,46,sortedData); 
# cor_plot2(50,52,sortedData); cor_plot2(58,62,sortedData); cor_plot2(66,68,sortedData);
# cor_plot2(72,74,sortedData); cor_plot2(78,80,sortedData); cor_plot2(84,86,sortedData); 
# cor_plot2(90,92,sortedData); cor_plot2(98,102,sortedData); cor_plot2(106,108,sortedData);
# cor_plot2(112,114,sortedData); cor_plot2(118,120,sortedData); cor_plot2(124,126,sortedData); 
# cor_plot2(130,132,sortedData); cor_plot2(136,138,sortedData); cor_plot2(142,144,sortedData)

# 
# png(paste(Output,"1_5_cor.png",sep=""))
# chart.Correlation(sortedData[,1:5], histogram=FALSE, pch="18")
# dev.off()

HG <-c(seq(1,5),seq(11,13), seq(17,19), seq(23,25), seq(29,31), seq(35,37), seq(41,43), seq(47,49), 
       seq(53,57), seq(63,65), seq(69,71), seq(75,77), seq(81,83), seq(87,89), seq(93,97), seq(103,105),
       seq(109,111), seq(115,117), seq(121,123), seq(127,129), seq(133,135), seq(139,141))

LG <-c(seq(6,10),seq(14,16), seq(20,22), seq(26,28), seq(32,34), seq(38,40), seq(44,46), seq(50,52), 
       seq(58,62), seq(66,68), seq(72,74), seq(78,80), seq(84,86), seq(90,92), seq(98,102), seq(106,108),
       seq(112,114), seq(118,120), seq(124,126), seq(130,132), seq(136,138), seq(142,144))
HGData<-sortedData[,HG]
LGData<-sortedData[,LG]
deltaData<-sweep(HGData,1,LGData,"-") # delta x
deltaData.T<-data.frame(t(deltaData))
deltaData.T$Name<-substr(colnames(deltaData),1,nchar(colnames(deltaData))-2)

delta_cor_plot <- function(a,b, data) {
  M <- cor(data[,a:b])
  png(paste(Figs,a,"_",b,"_delta_cor_p1.png", sep=""))
  corrplot.mixed(M, lower="number", upper="circle")
  dev.off()
}

delta_cor_plot(1,5,deltaData)
delta_cor_plot(6,8,deltaData);delta_cor_plot(9,11,deltaData);
delta_cor_plot(12,14,deltaData); delta_cor_plot(15,17,deltaData);
delta_cor_plot(18,20,deltaData);delta_cor_plot(21,23,deltaData); 
delta_cor_plot(24,26,deltaData); delta_cor_plot(27,31,deltaData); 
delta_cor_plot(32,34,deltaData); delta_cor_plot(35,37,deltaData);
delta_cor_plot(38,40,deltaData);delta_cor_plot(41,43,deltaData); 
delta_cor_plot(44,46,deltaData); delta_cor_plot(47,51,deltaData); 
delta_cor_plot(52,54,deltaData); delta_cor_plot(55,57,deltaData);
delta_cor_plot(58,60,deltaData);delta_cor_plot(61,63,deltaData); 
delta_cor_plot(64,66,deltaData); delta_cor_plot(67,69,deltaData); 
delta_cor_plot(70,72,deltaData); 

delta_cor_plot(1,72,deltaData); 

# library("PerformanceAnalytics")
# 
# delta_cor_plot2 <- function(a,b, data) {
#   png(paste(Figs,a,"_",b,"_delta_cor2.png", sep=""))
#   chart.Correlation(data, histogram=FALSE, pch="18")
#   dev.off()
# }
# 
# delta_cor_plot2(1,5,deltaData)
# delta_cor_plot2(6,8,deltaData);delta_cor_plot2(9,11,deltaData);
# delta_cor_plot2(12,14,deltaData); delta_cor_plot2(15,17,deltaData);
# delta_cor_plot2(18,20,deltaData);delta_cor_plot2(21,23,deltaData); 
# delta_cor_plot2(24,26,deltaData); delta_cor_plot2(27,31,deltaData); 
# delta_cor_plot2(32,34,deltaData); delta_cor_plot2(35,37,deltaData);
# delta_cor_plot2(38,40,deltaData);delta_cor_plot2(41,43,deltaData); 
# delta_cor_plot2(44,46,deltaData); delta_cor_plot2(47,51,deltaData); 
# delta_cor_plot2(52,54,deltaData); delta_cor_plot2(55,57,deltaData);
# delta_cor_plot2(58,60,deltaData);delta_cor_plot2(61,63,deltaData); 
# delta_cor_plot2(64,66,deltaData); delta_cor_plot2(67,69,deltaData); 
# delta_cor_plot2(70,72,deltaData); 
# 
# delta_cor_plot2(1,72,deltaData); 


meanDeltaData.T<-aggregate(deltaData.T[,-ncol(deltaData.T)],list(deltaData.T$Name),mean)

# anova for inter-individual
group1=seq(1,nrow(meanDeltaData.T))
group2=c(rep(0,8),rep(1,7),rep(2,7))
groupName=substr(meanDeltaData.T$Group.1,1,nchar(meanDeltaData.T$Group)-5); 
rownames(meanDeltaData.T)<-groupName
meanDeltaData.T$Group.1<-NULL

meanDeltaData <- t(meanDeltaData.T)
delta_cor_plot(1,22,meanDeltaData)

rankByMAD <- function(x) {
  rowWise<-1; colWise<-2
  x1<-apply(x,colWise,as.numeric)
  rownames(x1)<-rownames(x)
  x1<-data.frame(x1)
  x1$mad<-apply(x1,rowWise,mad)
  madScrClmn<-ncol(x1)
  z <- x1[order(-x1[,madScrClmn]),] # sort by mad score
  return(z)
}

sorteddeltaDataMAD<-rankByMAD(deltaData)

sorteddeltaDataMAD1000<-sorteddeltaDataMAD[1:100, -ncol(sorteddeltaDataMAD)]

delta_cor_plot(1,5,sorteddeltaDataMAD1000)
delta_cor_plot(6,8,sorteddeltaDataMAD1000);delta_cor_plot(9,11,sorteddeltaDataMAD1000);
delta_cor_plot(12,14,sorteddeltaDataMAD1000); delta_cor_plot(15,17,sorteddeltaDataMAD1000);
delta_cor_plot(18,20,sorteddeltaDataMAD1000);delta_cor_plot(21,23,sorteddeltaDataMAD1000); 
delta_cor_plot(24,26,sorteddeltaDataMAD1000); delta_cor_plot(27,31,sorteddeltaDataMAD1000); 
delta_cor_plot(32,34,sorteddeltaDataMAD1000); delta_cor_plot(35,37,sorteddeltaDataMAD1000);
delta_cor_plot(38,40,sorteddeltaDataMAD1000);delta_cor_plot(41,43,sorteddeltaDataMAD1000); 
delta_cor_plot(44,46,sorteddeltaDataMAD1000); delta_cor_plot(47,51,sorteddeltaDataMAD1000); 
delta_cor_plot(52,54,sorteddeltaDataMAD1000); delta_cor_plot(55,57,sorteddeltaDataMAD1000);
delta_cor_plot(58,60,sorteddeltaDataMAD1000);delta_cor_plot(61,63,sorteddeltaDataMAD1000); 
delta_cor_plot(64,66,sorteddeltaDataMAD1000); delta_cor_plot(67,69,sorteddeltaDataMAD1000); 
delta_cor_plot(70,72,sorteddeltaDataMAD1000); 

delta_cor_plot(1,72,sorteddeltaDataMAD1000); 

meanDeltaDataMAD<-rankByMAD(meanDeltaData)
meanDeltaDataMAD1000<-meanDeltaDataMAD[1:50, -ncol(meanDeltaDataMAD)]
delta_cor_plot(1,22,meanDeltaDataMAD1000)