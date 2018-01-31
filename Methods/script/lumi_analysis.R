library(lumi); library(gtools); library(limma); library(bioDist); library(calibrate)

################################
##  directory setup
################################
HOME <- "/Users/sjung/Project/GlobusGenomics/UIC-310/"
InputDir <- paste(HOME,'InputData/', sep='')
InputFiles = list.files(path=InputDir, pattern='txt')
ControlDir <- paste(HOME,'ControlData/', sep='')
ControlFiles = list.files(path=ControlDir, pattern='txt')
PhenotypeDir = paste(HOME,'PhenotypeInfo/', sep='')
Output=paste(HOME,'output/',sep=''); dir.create(Output, showWarnings = FALSE)
Figs = paste(HOME,'Figs/', sep=''); dir.create(Figs, showWarnings = FALSE)

# temp directory setup
#---------------------
TEMP<-tempdir()
if (file.exists(TEMP)) {
  dir.create(file.path(TEMP, "/InputData"),showWarnings = FALSE)
  dir.create(file.path(TEMP, "/ControlData"), showWarnings = FALSE)
}

TempDir <- paste(TEMP,'/InputData/', sep='')
TempControlDir <- paste(TEMP,'/ControlData/', sep='')

source(paste(HOME,"analysis/func.R",sep=""))

################################
##      parameter setup
################################

# figures
#--------
boxplot=FALSE; density=FALSE; cv=FALSE; sampleRelation=FALSE; pca=FALSE
NumOfSample <- c(1,46,47,92,93,144)

# filtering
#----------
filgering=TRUE
normalization.m = "quantile"
filter.Th = 0.01; filter.dp = 0

################################
##    fix data file format
################################
for (i in 1:length(InputFiles)) {
  GSTextFile <- paste(InputDir,InputFiles[i], sep='')
  GSContProbeFile <- paste(ControlDir, ControlFiles[i], sep='')
  Text <- paste(TempDir,InputFiles[i], sep='')
  ContProbe <- paste(TempControlDir,ControlFiles[i], sep='')
  system(paste('sed /^\\s*$/d', GSTextFile, '>', Text, sep= ' '))
  system(paste('sed /^\\s*$/d', GSContProbeFile, '>', ContProbe, sep= ' '))
}

# import text file and control probe
#-----------------------------------
InputFileList = list.files(path=TempDir, pattern='txt', full.names = TRUE)
ControlFileList = list.files(path=TempControlDir, pattern='txt', full.names = TRUE)

# read data
#----------
lumi.Ds1 <- lumiR(InputFileList[1]); sampleNames(lumi.Ds1)<-make.names(sampleNames(lumi.Ds1))
lumi.Ds2 <- lumiR(InputFileList[2]); sampleNames(lumi.Ds2)<-make.names(sampleNames(lumi.Ds2))
lumi.Ds3 <- lumiR(InputFileList[3]); sampleNames(lumi.Ds3)<-make.names(sampleNames(lumi.Ds3))
lumi.Ds4 <- lumiR(InputFileList[4]); sampleNames(lumi.Ds4)<-make.names(sampleNames(lumi.Ds4))

# add control data
#-----------------
rawData.ctrl1 <- addControlData2lumi(ControlFileList[1],lumi.Ds1)
rawData.ctrl2 <- addControlData2lumi(ControlFileList[2],lumi.Ds2)
rawData.ctrl3 <- addControlData2lumi(ControlFileList[3],lumi.Ds3)
rawData.ctrl4 <- addControlData2lumi(ControlFileList[4],lumi.Ds4)

# add NuID
#---------
if (require(lumiHumanIDMapping)){
  rawData.ctrl1.nuId <- addNuID2lumi(rawData.ctrl1, lib.mapping='lumiHumanIDMapping')
  rawData.ctrl2.nuId <- addNuID2lumi(rawData.ctrl2, lib.mapping='lumiHumanIDMapping')
  rawData.ctrl3.nuId <- addNuID2lumi(rawData.ctrl3, lib.mapping='lumiHumanIDMapping')
  rawData.ctrl4.nuId <- addNuID2lumi(rawData.ctrl4, lib.mapping='lumiHumanIDMapping')
}

rm(lumi.Ds1, lumi.Ds2, lumi.Ds3, lumi.Ds4)

## data processing (background correction, variance stabilization, quantile normalization)
normData1 <-  lumiExpresso(rawData.ctrl1.nuId, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), variance.stabilize = TRUE, 
                           varianceStabilize.param = list(method='log2'), normalize = TRUE, normalize.param = list(method='quantile'), QC.evaluation = TRUE 
)
normData2 <-  lumiExpresso(rawData.ctrl2.nuId, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), variance.stabilize = TRUE, 
                           varianceStabilize.param = list(method='log2'), normalize = TRUE, normalize.param = list(method='quantile'), QC.evaluation = TRUE
)
normData3 <-  lumiExpresso(rawData.ctrl3.nuId, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), variance.stabilize = TRUE, 
                           varianceStabilize.param = list(method='log2'), normalize = TRUE, normalize.param = list(method='quantile'), QC.evaluation = TRUE 
)
normData4 <-  lumiExpresso(rawData.ctrl4.nuId, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), variance.stabilize = TRUE, 
                           varianceStabilize.param = list(method='log2'), normalize = TRUE, normalize.param = list(method='quantile'), QC.evaluation = TRUE
)

rm(rawData.ctrl1, rawData.ctrl2, rawData.ctrl3, rawData.ctrl4)

# rawData expression data
#---------------------------------
eset.rawData1 <- exprs(rawData.ctrl1.nuId);
eset.rawData2 <- exprs(rawData.ctrl2.nuId);
eset.rawData3 <- exprs(rawData.ctrl3.nuId);
eset.rawData4 <- exprs(rawData.ctrl4.nuId);

# normData expression data
#---------------------------------
eset.normData1 <- exprs(normData1);
eset.normData2 <- exprs(normData2);
eset.normData3 <- exprs(normData3);
eset.normData4 <- exprs(normData4);

# find common probes to merge data
#---------------------------------
dataList <- list(eset.rawData1,eset.rawData2,eset.rawData3,eset.rawData4)
common_probe_names = Reduce(intersect, lapply(dataList, row.names))
rm(dataList)

#combine data
#------------
options(warn=-1) #supress warning for the sample name discrepency in control data

rawData<-combine(rawData.ctrl1.nuId[common_probe_names,],rawData.ctrl2.nuId[common_probe_names,],
                 rawData.ctrl3.nuId[common_probe_names,],rawData.ctrl4.nuId[common_probe_names,])

normData<-combine(normData1[common_probe_names,],normData2[common_probe_names,],
                  normData3[common_probe_names,],normData4[common_probe_names,])

rm(rawData.ctrl1.nuId, rawData.ctrl2.nuId, rawData.ctrl3.nuId, rawData.ctrl4.nuId,
   eset.rawData1,eset.rawData2,eset.rawData3,eset.rawData4,normData1, normData2, 
   normData3, normData4,eset.normData1,eset.normData2,eset.normData3,eset.normData4)

# replace sample name
#--------------------
sampleNames(rawData)<-gsub("SC.", "", sampleNames(rawData))
sampleNames(normData)<-gsub("SC.", "", sampleNames(normData))

sampleInfo<-read.table(paste(PhenotypeDir,"20160203.MAG.layout for H12 chips_updated.txt",sep=''))
sampleInfo<-sampleInfo[mixedorder(sampleInfo$V1),]
sampleNames(rawData)<-sampleInfo$V2; sampleNames(normData)<-sampleInfo$V2

# sort data by sample name
#-------------------------
rawData.sorted<-rawData[,mixedorder(sampleNames(rawData), decreasing=FALSE)]
normData.sorted<-normData[,mixedorder(sampleNames(normData), decreasing=FALSE)]
rm(eset.rawData,eset.normData,rawData,normData)

# Create array groups, array names and plot variables
#----------------------------------------------------
experimentFactor <- as.factor(substr(colnames(rawData.sorted), 1, nchar(colnames(rawData.sorted))-2))
colList          <- colorsByFactor(experimentFactor)
plotColors       <- colList$plotColors
legendColors     <- colList$legendColors
rm(colList)

##########################################
##          Generate QC plots
##########################################

# 1. Create symbol set for the array groups !! FIX HERE
#------------------------------------------
plotSymbols <- 18-as.numeric(experimentFactor)
legendSymbols <- sort(unique(plotSymbols), decreasing=TRUE)

# Define display parameters for the images  		                   
#-----------------------------------------
WIDTH <- 1000
HEIGHT <- 1414
POINTSIZE <- 24

if(!exists("MAXARRAY")) MAXARRAY <- 55

#----------------------------
# Data Quality Control graphs
# box plot, density plot, cv plot, clustering plot, pca plot
#-----------------------------------------------------------
if(boxplot) {
  print ("Plot boxplot for intensities") 
  for (i in 1:as.integer(length(NumOfSample)/2)){
    rawFigName=paste(Figs, "Raw_", NumOfSample[2*i-1], "_", NumOfSample[2*i-1], sep="")
    normFigName=paste(Figs, "Norm_", NumOfSample[2*i-1], "_", NumOfSample[2*i-1], sep="")
    box.plot(rawData.sorted[,NumOfSample[2*i-1]:NumOfSample[2*i]], rawFigName, col=plotColors, maxArray=MAXARRAY)
    box.plot(normData.sorted[,NumOfSample[2*i-1]:NumOfSample[2*i]], normFigName, col=plotColors, maxArray=MAXARRAY)
  }
  
  box.plot(rawData.sorted, paste(Figs, "Raw_1_144", sep=""), col=plotColors, maxArray=144)
  box.plot(normData.sorted, paste(Figs, "Norm_1_144", sep=""), col=plotColors, maxArray=144)
}

if(density) {
  print ("Plot density histogram for intensities")
  
  # raw and normalized data
  for (i in 1:as.integer(length(NumOfSample)/2)){
    rawFigName=paste(Figs, "Raw_", NumOfSample[2*i-1], "_", NumOfSample[2*i-1], sep="")
    normFigName=paste(Figs, "Norm_", NumOfSample[2*i-1], "_", NumOfSample[2*i-1], sep="")
    density.plot(rawData.sorted[,NumOfSample[2*i-1]:NumOfSample[2*i]], rawFigName, col=plotColors, maxArray=MAXARRAY)
    density.plot(normData.sorted[,NumOfSample[2*i-1]:NumOfSample[2*i]], normFigName, col=plotColors, maxArray=MAXARRAY)
  }
  # all the raw and normalized data
  density.plot(rawData.sorted, paste(Figs, "Raw_1_144", sep=""), col=plotColors, maxArray=144)
  density.plot(normData.sorted, paste(Figs, "Norm_1_144", sep=""), col=plotColors, maxArray=144)
}

if(cv) {
  print ("Plot density for coefficient of varience for intensities")
  
  # raw and normalized data
  for (i in 1:as.integer(length(NumOfSample)/2)){
    rawFigName=paste(Figs, "Raw_", NumOfSample[2*i-1], "_", NumOfSample[2*i-1], sep="")
    normFigName=paste(Figs, "Norm_", NumOfSample[2*i-1], "_", NumOfSample[2*i-1], sep="")
    cv.plot(rawData.sorted[,NumOfSample[2*i-1]:NumOfSample[2*i]], rawFigName, col=plotColors, maxArray=MAXARRAY)
    cv.plot(normData.sorted[,NumOfSample[2*i-1]:NumOfSample[2*i]], normFigName, col=plotColors, maxArray=MAXARRAY)
  }
  
  cv.plot(rawData.sorted, paste(Figs, "Raw_1_144", sep=""), col=plotColors, maxArray=144)
  cv.plot(normData.sorted, paste(Figs, "Norm_1_144", sep=""), col=plotColors, maxArray=144)
}

if(sampleRelation) {
  print ("Hierarchical clustering of data") 
  clusterOption1 = "spearman"; clusterOption2 = "complete"
  
  # raw data
  for (i in 1:as.integer(length(NumOfSample)/2)){
    clusterFun(rawData.sorted, range1=NumOfSample[2*i-1], range2=NumOfSample[2*i], normalized=FALSE, 
               experimentFactor=experimentFactor, clusterOption1=clusterOption1, clusterOption2=clusterOption2,
               plotColors=plotColors, legendColors=legendColors, plotSymbols=plotSymbols, legendSymbols=legendSymbols,
               WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=MAXARRAY) 
  }
  
  # raw data alltogether
  clusterFun(rawData.sorted, range1=1, range2=ncol(rawData.sorted), normalized=FALSE,  experimentFactor=experimentFactor,
             clusterOption1=clusterOption1, clusterOption2=clusterOption2,
             plotColors=plotColors, legendColors=legendColors,
             plotSymbols=plotSymbols, legendSymbols=legendSymbols,
             WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=12,MAXARRAY=144) 
  
  # normalized data
  for (i in 1:as.integer(length(NumOfSample)/2)){
    clusterFun(normData.sorted, range1=NumOfSample[2*i-1], range2=NumOfSample[2*i], normalized=TRUE,  
               experimentFactor=experimentFactor, clusterOption1=clusterOption1, clusterOption2=clusterOption2,
               plotColors=plotColors, legendColors=legendColors,plotSymbols=plotSymbols, legendSymbols=legendSymbols,
               WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=MAXARRAY) 
  }
  
  # normalized data alltogether
  clusterFun(normData.sorted, range1=1, range2=ncol(normData.sorted), normalized=TRUE,  experimentFactor=experimentFactor,
             clusterOption1=clusterOption1, clusterOption2=clusterOption2,
             plotColors=plotColors, legendColors=legendColors,
             plotSymbols=plotSymbols, legendSymbols=legendSymbols,
             WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=12,MAXARRAY=144) 
  
  # normalized data without replicates
  uniqueNormalData<-normData.sorted[,!duplicated(substr(colnames(normData.sorted),1,nchar(colnames(normData.sorted))-2))]
  clusterFun(uniqueNormalData, range1=1, range2=ncol(uniqueNormalData), normalized=TRUE,  experimentFactor=experimentFactor,
             clusterOption1=clusterOption1, clusterOption2=clusterOption2,
             plotColors=plotColors, legendColors=legendColors,
             plotSymbols=plotSymbols, legendSymbols=legendSymbols,
             WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=MAXARRAY) 
}

if(pca) {  
  print("PCA graph for data")
  groupsInLegend =  !( length(unique(levels(experimentFactor))) ) >=10
  
  # raw data
  for (i in 1:as.integer(length(NumOfSample)/2)){
    pcaFun(rawData.sorted, range1=NumOfSample[2*i-1], range2=NumOfSample[2*i], normalized=FALSE ,experimentFactor=experimentFactor, 
           plotColors=plotColors, legendColors=legendColors, plotSymbols=plotSymbols, legendSymbols=legendSymbols, 
           groupsInLegend=groupsInLegend ,namesInPlot=((max(nchar(sampleNames(rawData.sorted)))<=10) &&
                                                         (length(sampleNames(rawData.sorted))<=(MAXARRAY/2))),WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
    
  }
  # raw data alltogether 
  pcaFun(rawData.sorted, range1=1, range2=ncol(rawData.sorted), normalized=FALSE ,experimentFactor=experimentFactor, 
         plotColors=plotColors, legendColors=legendColors, plotSymbols=plotSymbols, legendSymbols=legendSymbols, 
         groupsInLegend=groupsInLegend, namesInPlot=((max(nchar(sampleNames(rawData.sorted)))<=10) &&
                                                       (length(sampleNames(rawData.sorted))<=(MAXARRAY/2))),WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=12)
  
  # normalized data
  for (i in 1:as.integer(length(NumOfSample)/2)){
    pcaFun(normData.sorted, range1=NumOfSample[2*i-1], range2=NumOfSample[2*i], normalized=TRUE, experimentFactor=experimentFactor, 
           plotColors=plotColors, legendColors=legendColors, plotSymbols=plotSymbols,legendSymbols=legendSymbols, groupsInLegend=groupsInLegend,
           namesInPlot=((max(nchar(sampleNames(rawData.sorted)))<=10) && (length(sampleNames(rawData.sorted))<=(MAXARRAY/2))),
           WIDTH=WIDTH,HEIGHT=HEIGHT, POINTSIZE=POINTSIZE)
  }
  # normalized data alltogether 
  pcaFun(normData.sorted, range1=1, range2=ncol(normData.sorted), normalized=TRUE, experimentFactor=experimentFactor, 
         plotColors=plotColors, legendColors=legendColors, plotSymbols=plotSymbols,legendSymbols=legendSymbols, groupsInLegend=groupsInLegend, 
         namesInPlot=((max(nchar(sampleNames(rawData.sorted)))<=10) && (length(sampleNames(rawData.sorted))<=(MAXARRAY/2))), 
         WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=12)
}

#################################
##        data analysis
#################################

# filtering
#----------------
filtering=TRUE
if(filtering) {
  # will save a filtered.normData file in the working dir.
  eset.filtered.normData = filterFun(rawData.sorted, normData.sorted, filter.Th, filter.dp)
}

# remove probes that do not have annotated genes
#------------------------------------------------
if (require(lumiHumanAll.db) & require(annotate)) {
  dataMatrix <- eset.filtered.normData[!is.na(getSYMBOL(rownames(eset.filtered.normData), 'lumiHumanAll.db')),]
}

probeList <- rownames(dataMatrix)
if (require(lumiHumanAll.db) & require(annotate)){
  geneSymbol<- getSYMBOL(probeList, 'lumiHumanAll.db')
  geneName<- sapply(lookUp(probeList, 'lumiHumanAll.db', 'GENENAME'), function(x) x[1])
}

sortedData<-dataMatrix[,mixedorder(colnames(dataMatrix), decreasing=FALSE)]
rm(dataMatrix, eset.filtered.normData)
# create a design matrix - manual modification is required
FLAG=FALSE
if (FLAG){
  write.table(colnames(sortedData), file=paste(PhenotypeDir,'target.txt',sep=''), row.names=FALSE, col.names=FALSE)
}

# differential expression analysis
#---------------------------------
targets<-readTargets(paste(PhenotypeDir,"target.txt", sep=''))
ct<-factor(paste(targets$Condition, targets$Treatment,sep="."))
design <- model.matrix(~0+ct)
colnames(design) <- levels(ct)
corfit<-duplicateCorrelation(sortedData,design,block=targets$Subject)
fit <- lmFit(sortedData, design, block=targets$Subject, correlation=corfit$consensus) 

contrast.matrix <- makeContrasts(
  DwC.Insulin_Con = disease_comp.insulin-disease_comp.control,
  DwoC.Insulin_Con = disease.insulin-disease.control,
  Normal.Insulin_Con = normal.insulin-normal.control,
  levels=design)    
fit2<-contrasts.fit(fit, contrast.matrix)
fit2<-eBayes(fit2);  
fit2$genes <-data.frame(ID=probeList, geneSymbol=geneSymbol, geneName=geneName, stringsAsFactors=FALSE)

# 1. Disease with Complication Insulin vs Control
# 2. Disease without Complication Insulin vs Control
# 3. Normal Insulin vs Control
#---------------------------------------------------
DwC.Insulin_DwC.Con <- topTable(fit2, coef="DwC.Insulin_Con", number=nrow(fit2$genes), sort.by="logFC", resort.by="p", genelist=fit2$genes, adjust.method="BH")
DwoC.Insulin_DwoC.Con <- topTable(fit2, coef="DwoC.Insulin_Con", number=nrow(fit2$genes), sort.by="logFC", resort.by="p", genelist=fit2$genes, adjust.method="BH")
Normal.Insulin_Normal.Con <- topTable(fit2, coef="Normal.Insulin_Con", number=nrow(fit2$genes), sort.by="logFC", resort.by="p", genelist=fit2$genes, adjust.method="BH")

DwC.Insulin_DwC.Con_Pval0.05 <- DwC.Insulin_DwC.Con[DwC.Insulin_DwC.Con$P.Value <= 0.05,]
write.csv(DwC.Insulin_DwC.Con_Pval0.05, file=paste(Output,"DwC.Insulin_DwC.Con_Pval0.05.csv",sep=""), row.names=FALSE)

DwoC.Insulin_DwoC.Con_Pval0.05 <- DwoC.Insulin_DwoC.Con[DwoC.Insulin_DwoC.Con$P.Value <= 0.05,]
write.csv(DwoC.Insulin_DwoC.Con_Pval0.05, file=paste(Output,"DwoC.Insulin_DwoC.Con_Pval0.05.csv",sep=""), row.names=FALSE)

Normal.Insulin_Normal.Con_Pval0.05 <- Normal.Insulin_Normal.Con[Normal.Insulin_Normal.Con$P.Value <= 0.05,]
write.csv(Normal.Insulin_Normal.Con_Pval0.05, file=paste(Output,"Normal.Insulin_Normal.Con_Pval0.05.csv",sep=""), row.names=FALSE)

pdf(file = paste(Figs,"DwC.Insulin_DwC.Con_Pval0.05.pdf",sep=""),width=18,height=12,pointsize=20)
with(DwC.Insulin_DwC.Con, plot(logFC, -log10(P.Value), pch=20, col="gray", main="Diabets w/ complication (P.val <0.05)", xlab="log2 fold change", ylab="-log10(p-value)", xlim=c(-0.5,0.5)))
with(subset(DwC.Insulin_DwC.Con, P.Value<.05), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(DwC.Insulin_DwC.Con, P.Value<.05), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()

pdf(file = paste(Figs,"DwoC.Insulin_DwoC.Con_Pval0.05.pdf",sep=""),width=10,height=10,pointsize=20)
with(DwoC.Insulin_DwoC.Con, plot(logFC, -log10(P.Value), pch=20, col="gray", main="Diabets w/o complication (P.val <0.05)", xlab="log2 fold change", ylab="-log10(p-value)", xlim=c(-0.5,0.5)))
with(subset(DwoC.Insulin_DwoC.Con, P.Value<.05), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(DwoC.Insulin_DwoC.Con, P.Value<.05), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()

pdf(file = paste(Figs,"Normal.Insulin_Normal.Con_Pval0.05.pdf",sep=""),width=18,height=12, pointsize=20)
with(Normal.Insulin_Normal.Con, plot(logFC, -log10(P.Value), pch=20, col="gray", main="No disease (P.val <0.05)", xlab="log2 fold change", ylab="-log10(p-value)", xlim=c(-0.5,0.5)))
with(subset(Normal.Insulin_Normal.Con, P.Value<.05), points(logFC, -log10(P.Value), pch=20, col="blue"), cex=1.1)
with(subset(Normal.Insulin_Normal.Con, P.Value<.05), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()
