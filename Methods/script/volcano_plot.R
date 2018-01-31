library(lumi); library(gtools); library(limma); library(bioDist); library(calibrate)

##---------------
##   functions
##---------------
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

##--------------------
##  directory setup
##--------------------
HOME <- "/Users/sjung/Project/GlobusGenomics/UIC/310/"
Output=paste(HOME,'output/',sep=''); dir.create(Output, showWarnings = FALSE)
Figs = paste(HOME,'Figs/', sep=''); dir.create(Figs, showWarnings = FALSE)

load(file = paste(Output,"meanDeltaData.RData",sep=""))
load(file = paste(Output,"sortedDataGeneAnno.RData",sep=""))

meanDeltaDataMAD<-rankByMAD(meanDeltaData)

library(qvalue)
HOME<-"/Users/sjung/Project/GlobusGenomics/UIC/"
setwd(HOME)
data<-read.csv(file="DeltaAnalysisResults6.15.16.csv")
colnames(data)[4]<-"pval.equal.var"
colnames(data)[5]<-"pval.treat.eff"
p1<-data$pval.treat.eff
qval.treat.eff<-qvalue(p1)
data$qval.treat.eff<-qval.treat.eff$qvalues

rownames(data)<-data$geneID
dataSortedbyMAD<-data[rownames(meanDeltaDataMAD),]
dataSortedbyMAD$mad <- meanDeltaDataMAD$mad
rm(data)

write.csv(dataSortedbyMAD,"dataSortedbyMAD.csv",row.names = FALSE)

pdf(file = paste(Figs,"test1.pdf",sep=""),width=18,height=12,pointsize=20)
with(dataSortedbyMAD, plot(mad, -log10(qval.treat.eff), pch=20, col="gray", main="q-value vs variation ", xlab="mean absolute deviation", ylab="-log10(q-value)", xlim=c(0,1.0)))
with(subset(dataSortedbyMAD, qval.treat.eff<.1 & mad > .2), points(mad, -log10(qval.treat.eff), pch=20, col="red"), cex=1.1)
with(subset(dataSortedbyMAD, qval.treat.eff<.1 & mad > .2), textxy(mad, -log10(qval.treat.eff), labs=symbol, cex=.8))
dev.off()
# 
# pdf(file = paste(Figs,"DwoC.Insulin_DwoC.Con_Pval0.05.pdf",sep=""),width=10,height=10,pointsize=20)
# with(DwoC.Insulin_DwoC.Con, plot(logFC, -log10(P.Value), pch=20, col="gray", main="Diabets w/o complication (P.val <0.05)", xlab="log2 fold change", ylab="-log10(p-value)", xlim=c(-0.5,0.5)))
# with(subset(DwoC.Insulin_DwoC.Con, P.Value<.05), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
# with(subset(DwoC.Insulin_DwoC.Con, P.Value<.05), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
# dev.off()
# 
# pdf(file = paste(Figs,"Normal.Insulin_Normal.Con_Pval0.05.pdf",sep=""),width=18,height=12, pointsize=20)
# with(Normal.Insulin_Normal.Con, plot(logFC, -log10(P.Value), pch=20, col="gray", main="No disease (P.val <0.05)", xlab="log2 fold change", ylab="-log10(p-value)", xlim=c(-0.5,0.5)))
# with(subset(Normal.Insulin_Normal.Con, P.Value<.05), points(logFC, -log10(P.Value), pch=20, col="blue"), cex=1.1)
# with(subset(Normal.Insulin_Normal.Con, P.Value<.05), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
# dev.off()
