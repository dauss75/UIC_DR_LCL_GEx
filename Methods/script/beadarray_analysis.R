library(beadarray); library(gtools); library(limma)

## function
# MAD: mean absolute deviation
RANKbyMAD <- function(x) {
  x1<-apply(x,2,as.numeric)
  x1<-data.frame(x1)
  x1$mad<-apply(x1,1,mad)
  # x1$genes<-y[,1]
  rownames(x1)<-rownames(x)
  madScrClm<-ncol(x1)
  #   print(madScrClmn)
  z <- x1[order(-x1[,madScrClm]),] # sort by mad score
  return(z[,1:ncol(x1)-1])
}

## CV:coefficient of variation
CV <- function(x,na.rm=TRUE) {
  x1<-apply(x,2,as.numeric)
  x1_SD<-apply(x1,1,sd)
  x1_median<-apply(x1,1,median)  
  y1<-x1_SD/x1_median
  return(y1)
}

HOME <- "/Users/sjung/Project/GlobusGenomics/UIC-310/"
InputDir <- paste(HOME,'InputData/', sep='')
InputFiles = list.files(path=InputDir, pattern='txt')
ControlDir <- paste(HOME,'ControlData/', sep='')
ControlFiles = list.files(path=ControlDir, pattern='txt')
PhenotypeDir = paste(HOME,'PhenotypeInfo/', sep='')
Figs = paste(HOME,'Figs/', sep='')
dir.create(Figs, showWarnings = FALSE)

TEMP<-tempdir()
if (file.exists(TEMP)) {
  dir.create(file.path(TEMP, "/InputData"),showWarnings = FALSE)
  dir.create(file.path(TEMP, "/ControlData"), showWarnings = FALSE)
}

TempDir <- paste(TEMP,'/InputData/', sep='')
TempControlDir <- paste(TEMP,'/ControlData/', sep='')

for (i in 1:length(InputFiles)) {
  GSTextFile <- paste(InputDir,InputFiles[i], sep='')
  GSContProbeFile <- paste(ControlDir, ControlFiles[i], sep='')
  Text <- paste(TempDir,InputFiles[i], sep='')
  ContProbe <- paste(TempControlDir,ControlFiles[i], sep='')
  system(paste('sed /^\\s*$/d', GSTextFile, '>', Text, sep= ' '))
  system(paste('sed /^\\s*$/d', GSContProbeFile, '>', ContProbe, sep= ' '))
}

## import text file and control probe
InputFileList = list.files(path=TempDir, pattern='txt', full.names = TRUE)
ControlFileList = list.files(path=TempControlDir, pattern='txt', full.names = TRUE)

## read data
eset <- readBeadSummaryData(dataFile=InputFileList[1], qcFile=ControlFileList[1],
                            ProbeID="ProbeID", controlID="ProbeID",
                            skip=0, qc.skip=0,
                            annoCols=c("SYMBOL", "DEFINITION", "SYNONYMS", "CHROMOSOME", "ILMN_GENE", "SEARCH_KEY"))



lumi.Ds1 <- lumiR(InputFileList[1]); lumi.Ds2 <- lumiR(InputFileList[2]); 
lumi.Ds3 <- lumiR(InputFileList[3]); lumi.Ds4 <- lumiR(InputFileList[4]);

## add control data
lumi.Ds1 <- addControlData2lumi(ControlFileList[1],lumi.Ds1)
lumi.Ds2 <- addControlData2lumi(ControlFileList[2],lumi.Ds2)
lumi.Ds3 <- addControlData2lumi(ControlFileList[3],lumi.Ds3)
lumi.Ds4 <- addControlData2lumi(ControlFileList[4],lumi.Ds4)

## add NuID
if (require(lumiHumanIDMapping)){
  lumi1.nuId <- addNuID2lumi(lumi.Ds1, lib.mapping='lumiHumanIDMapping')
  lumi2.nuId <- addNuID2lumi(lumi.Ds2, lib.mapping='lumiHumanIDMapping')
  lumi3.nuId <- addNuID2lumi(lumi.Ds3, lib.mapping='lumiHumanIDMapping')
  lumi4.nuId <- addNuID2lumi(lumi.Ds4, lib.mapping='lumiHumanIDMapping')
}

rm(lumi.Ds1, lumi.Ds2, lumi.Ds3, lumi.Ds4)

## data processing (background correction, variance stabilization, quantile normalization)
lumi.N.Q1 <-  lumiExpresso(lumi1.nuId, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), variance.stabilize = TRUE, 
                           varianceStabilize.param = list(), normalize = TRUE, normalize.param = list(method='rsn'), QC.evaluation = TRUE, 
                           QC.param = list(), verbose = TRUE)
lumi.N.Q2 <-  lumiExpresso(lumi2.nuId, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), variance.stabilize = TRUE, 
                           varianceStabilize.param = list(), normalize = TRUE, normalize.param = list(method='rsn'), QC.evaluation = TRUE, 
                           QC.param = list(), verbose = TRUE)
lumi.N.Q3 <-  lumiExpresso(lumi3.nuId, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), variance.stabilize = TRUE, 
                           varianceStabilize.param = list(), normalize = TRUE, normalize.param = list(method='rsn'), QC.evaluation = TRUE, 
                           QC.param = list(), verbose = TRUE)
lumi.N.Q4 <-  lumiExpresso(lumi4.nuId, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), variance.stabilize = TRUE, 
                           varianceStabilize.param = list(), normalize = TRUE, normalize.param = list(method='rsn'), QC.evaluation = TRUE, 
                           QC.param = list(), verbose = TRUE)

## try individual
CorrectBg1<-lumiB(lumi1.nuId ) #background correction
CorrectBg2<-lumiB(lumi2.nuId ) #background correction
CorrectBg3<-lumiB(lumi3.nuId ) #background correction
CorrectBg4<-lumiB(lumi4.nuId ) #background correction

lumi.T1 <- lumiT(CorrectBg1, method='log2') # Transfer the Illumina data to stabilize the variance
lumi.T2 <- lumiT(CorrectBg2, method='log2') # Transfer the Illumina data to stabilize the variance
lumi.T3 <- lumiT(CorrectBg3, method='log2') # Transfer the Illumina data to stabilize the variance
lumi.T4 <- lumiT(CorrectBg4, method='log2') # Transfer the Illumina data to stabilize the variance

lumi.N1 <- lumiN(lumi.T1) # quantile normalization
lumi.N2 <- lumiN(lumi.T2) # quantile normalization
lumi.N3 <- lumiN(lumi.T3) # quantile normalization
lumi.N4 <- lumiN(lumi.T4) # quantile normalization

lumi.N.Q1 <- lumiQ(lumi.N1) 
lumi.N.Q2 <- lumiQ(lumi.N2) 
lumi.N.Q3 <- lumiQ(lumi.N3) 
lumi.N.Q4 <- lumiQ(lumi.N4) 

rm(lumi1.nuId, lumi2.nuId, lumi3.nuId, lumi4.nuId)

plot(lumi.N.Q1, what='density')
plot(lumi.N.Q1, what='boxplot')
plot(lumi.N.Q2, what='density')
plot(lumi.N.Q2, what='boxplot')
plot(lumi.N.Q3, what='density')
plot(lumi.N.Q3, what='boxplot')
plot(lumi.N.Q4, what='density')
plot(lumi.N.Q4, what='boxplot')
# find common probes to merge data
exprs1<-exprs(lumi.N.Q1); exprs2<-exprs(lumi.N.Q2); 
exprs3<-exprs(lumi.N.Q3); exprs4<-exprs(lumi.N.Q4)
dataList <- list(exprs1,exprs2,exprs3,exprs4)
common_probe_names = Reduce(intersect, lapply(dataList, row.names))
rm(dataList, lumi.N.Q1, lumi.N.Q2, lumi.N.Q3, lumi.N.Q4)

#combine data
dataMatrix<-combine(exprs1[common_probe_names,],exprs2[common_probe_names,],exprs3[common_probe_names,],exprs4[common_probe_names,])
rm(exprs1,exprs2,exprs3,exprs4)

colnames(dataMatrix)=gsub("SC-", "", colnames(dataMatrix))
## replace sample name
sampleInfo<-read.table(paste(PhenotypeDir,"20160203.MAG.layout for H12 chips_updated.txt",sep=''))
sampleInfo<-sampleInfo[mixedorder(sampleInfo$V1),]
colnames(dataMatrix)<-sampleInfo$V2

## remove probes that do not have annotated genes
if (require(lumiHumanAll.db) & require(annotate)) {
  dataMatrix <- dataMatrix[!is.na(getSYMBOL(rownames(dataMatrix), 'lumiHumanAll.db')),]
}

probeList <- rownames(dataMatrix)
if (require(lumiHumanAll.db) & require(annotate)){
  geneSymbol<- getSYMBOL(probeList, 'lumiHumanAll.db')
  geneName<- sapply(lookUp(probeList, 'lumiHumanAll.db', 'GENENAME'), function(x) x[1])
}
FLAG=FALSE
sortedData<-dataMatrix[,mixedorder(colnames(dataMatrix), decreasing=TRUE)]
if (FLAG==TRUE){
  write.table(colnames(sortedData), file=paste(PhenotypeDir,'target.txt',sep=''), row.names=FALSE, col.names=FALSE)
}
## ---- WILL DO LATER ----
## checking data quality  
# dataMatrix.T<-t(dataMatrix)
# x<-apply(dataMatrix.T,2,as.numeric)
# x<-data.frame(x)
# rownames(x)<-rownames(dataMatrix.T)
# dataMatrix.T<-cbind(dataMatrix.T,label=rownames(dataMatrix.T))
# sortedData<-dataMatrix[,mixedorder(colnames(dataMatrix), decreasing=TRUE)]
# temp<-RANKbyMAD(sortedData)
# for (i in seq(100,1000,100)){
#   fileName<-paste(Figs,"Gene_",i,".png",sep='')
#   png(file=fileName, width = 1024, height = 768)
#   d<-temp[1:i,]
#   # geneList<-paste("Gene_Top",i,".txt",sep="")
#   # write.table(rownames(d),file=geneList)
#   d2<-apply(d,2,as.numeric)
#   rownames(d2)<-rownames(d)
#   heatmap(d2)
#   dev.off()  
# }

# png(file=fileName, width = 1024, height = 768)
# z <- cor(temp)
# require(lattice)
# levelplot(z)
# dev.off()  

# pc <- prcomp(sortedData[,1:ncol(sortedData)-1], scale=TRUE)

# differential expression using limma
targets<-readTargets(paste(PhenotypeDir,"target.txt", sep=''))
ct<-factor(paste(targets$Condition, targets$Treatment,sep="."))
design <- model.matrix(~0+ct)
colnames(design) <- levels(ct)
corfit<-duplicateCorrelation(sortedData,design,block=targets$Subject)
fit <- lmFit(sortedData, design, block=targets$Subject, correlation=corfit$consensus)  
# fit <- lmFit(sortedData, design) 
contrast.matrix <- makeContrasts(
  DwC.Glucose_Con = disease_comp.Glucose-disease_comp.control,
  DwoC.Glucose_Con = disease.Glucose-disease.control,
  Normal.Glucose_Con = normal.Glucose-normal.control,
  levels=design)    
fit2<-contrasts.fit(fit, contrast.matrix)
fit2<-eBayes(fit2);  
fit2$genes <-data.frame(ID=probeList, geneSymbol=geneSymbol, geneName=geneName, stringsAsFactors=FALSE)

## A: DwC.Glucose vs DwC.Con
DwC.Glucose_DwC.Con <- topTable(fit2, coef="DwC.Glucose_Con", number=nrow(fit2$genes), sort.by="logFC", genelist=fit2$genes, adjust.method="BH")
DwoC.Glucose_DwoC.Con <- topTable(fit2, coef="DwoC.Glucose_Con", number=nrow(fit2$genes), sort.by="logFC", genelist=fit2$genes, adjust.method="BH")
Normal.Glucose_Normal.Con <- topTable(fit2, coef="Normal.Glucose_Con", number=nrow(fit2$genes), sort.by="logFC", genelist=fit2$genes, adjust.method="BH")
