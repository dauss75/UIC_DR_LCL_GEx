#------------------------------------------------------------
#  Library
#------------------------------------------------------------
library(lumi); library(gtools); library(limma); 
library(bioDist); library(calibrate)

#------------------------------------------------------------
#  Directory Setup
#------------------------------------------------------------
HOME <- "/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/"  # change if needed
InputDir <- paste(HOME,'data/input/', sep='')
InputFiles = list.files(path=InputDir, pattern='txt')
ControlDir <- paste(HOME,'data/control/', sep='')
ControlFiles = list.files(path=ControlDir, pattern='txt')
PhenotypeDir = paste(HOME,'data/phenotype/', sep='')
RData=paste(HOME,'rdata/',sep=''); dir.create(RData, showWarnings = FALSE)
Figs = paste(HOME,'figs/log2_rsn/', sep=''); dir.create(Figs, showWarnings = FALSE)
Output = paste(HOME,'output/', sep=''); dir.create(Figs, showWarnings = FALSE)

source(paste(HOME,'script/Function.R', sep=''))

TEMP<-tempdir()
if (file.exists(TEMP)) {
  dir.create(file.path(TEMP, "/InputData"),showWarnings = FALSE)
  dir.create(file.path(TEMP, "/ControlData"), showWarnings = FALSE)
}

TempDir <- paste(TEMP,'/InputData/', sep='')
TempControlDir <- paste(TEMP,'/ControlData/', sep='')

#------------------------------------------------------------
#  Fix Data File Format
#------------------------------------------------------------
for (i in 1:length(InputFiles)) {
  GSTextFile <- paste(InputDir,InputFiles[i], sep='')
  GSContProbeFile <- paste(ControlDir, ControlFiles[i], sep='')
  Text <- paste(TempDir,InputFiles[i], sep='')
  ContProbe <- paste(TempControlDir,ControlFiles[i], sep='')
  system(paste('sed /^\\s*$/d', GSTextFile, '>', Text, sep= ' '))
  system(paste('sed /^\\s*$/d', GSContProbeFile, '>', ContProbe, sep= ' '))
}

#------------------------------------------------------------
#  Import Treated and Baseline Data
#------------------------------------------------------------
InputFileList = list.files(path=TempDir, pattern='txt', full.names = TRUE)
ControlFileList = list.files(path=TempControlDir, pattern='txt', full.names = TRUE)

#------------------------------------------------------------
#  Lumi: 1. Read Expression Data 
#------------------------------------------------------------
lumi.Ds1 <- lumiR(InputFileList[1]); 
sampleNames(lumi.Ds1)<-make.names(sampleNames(lumi.Ds1))
lumi.Ds2 <- lumiR(InputFileList[2]); 
sampleNames(lumi.Ds2)<-make.names(sampleNames(lumi.Ds2))
lumi.Ds3 <- lumiR(InputFileList[3]); 
sampleNames(lumi.Ds3)<-make.names(sampleNames(lumi.Ds3))
lumi.Ds4 <- lumiR(InputFileList[4]); 
sampleNames(lumi.Ds4)<-make.names(sampleNames(lumi.Ds4))

#------------------------------------------------------------
#  Lumi: 2. Add Control Data 
#------------------------------------------------------------
rawData.ctrl1 <- addControlData2lumi(ControlFileList[1],lumi.Ds1)
rawData.ctrl2 <- addControlData2lumi(ControlFileList[2],lumi.Ds2)
rawData.ctrl3 <- addControlData2lumi(ControlFileList[3],lumi.Ds3)
rawData.ctrl4 <- addControlData2lumi(ControlFileList[4],lumi.Ds4)

if (require(lumiHumanIDMapping)){
  rawData.ctrl1.nuId <- addNuID2lumi(rawData.ctrl1, lib.mapping='lumiHumanIDMapping')
  rawData.ctrl2.nuId <- addNuID2lumi(rawData.ctrl2, lib.mapping='lumiHumanIDMapping')
  rawData.ctrl3.nuId <- addNuID2lumi(rawData.ctrl3, lib.mapping='lumiHumanIDMapping')
  rawData.ctrl4.nuId <- addNuID2lumi(rawData.ctrl4, lib.mapping='lumiHumanIDMapping')
}

rm(lumi.Ds1, lumi.Ds2, lumi.Ds3, lumi.Ds4)

#------------------------------------------------------------
#  Lumi: 3. Data Processing 
#------------------------------------------------------------

# b1. background corrected

Batch1 <-  lumiExpresso(rawData.ctrl1.nuId, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), 
                        variance.stabilize = FALSE, normalize = FALSE)
Batch2 <-  lumiExpresso(rawData.ctrl2.nuId, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), 
                        variance.stabilize = FALSE, normalize = FALSE)
Batch3 <-  lumiExpresso(rawData.ctrl3.nuId, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), 
                        variance.stabilize = FALSE, normalize = FALSE)
Batch4 <-  lumiExpresso(rawData.ctrl4.nuId, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), 
                        variance.stabilize = FALSE, normalize = FALSE)

#------------------------------------------------------------
# Merge batch data using common probes 
#------------------------------------------------------------
exprs.Batch1 <- exprs(Batch1);
exprs.Batch2 <- exprs(Batch2);
exprs.Batch3 <- exprs(Batch3);
exprs.Batch4 <- exprs(Batch4);

#------------------------------------------------------------
# find common probes and then combine data
#------------------------------------------------------------
dataList <- list(exprs.Batch1,exprs.Batch2,exprs.Batch3,exprs.Batch4)
common_probe_names = Reduce(intersect, lapply(dataList, row.names))
rm(dataList)

options(warn=-1) #supress warning for the sample name discrepency in control data
rawBatch<-combine(rawData.ctrl1.nuId[common_probe_names,],rawData.ctrl2.nuId[common_probe_names,],
                  rawData.ctrl3.nuId[common_probe_names,],rawData.ctrl4.nuId[common_probe_names,])

Batch <- combine(Batch1[common_probe_names,],Batch2[common_probe_names,],
                 Batch3[common_probe_names,],Batch4[common_probe_names,])

sampleNames(rawBatch)<-gsub("SC.", "", sampleNames(rawBatch))
sampleNames(Batch)<-gsub("SC.", "", sampleNames(Batch))

rm(Batch1,Batch2,Batch3,Batch4, common_probe_names, rawData.ctrl1.nuId, rawData.ctrl2.nuId, rawData.ctrl3.nuId, rawData.ctrl4.nuId)
rm(exprs.Batch1, exprs.Batch2, exprs.Batch3, exprs.Batch4, rawData.ctrl1, rawData.ctrl2, rawData.ctrl3, rawData.ctrl4)
#------------------------------------------------------------
# assign informative sample name and sort by sample name
#------------------------------------------------------------
sampleInfo<-read.table(paste(PhenotypeDir,"20160203.MAG.layout for H12 chips_updated.txt",sep=''))
sampleInfo<-sampleInfo[mixedorder(sampleInfo$V1),]
sampleNames(Batch)<-sampleInfo$V2
rawBatch.sorted<-Batch[,mixedorder(sampleNames(rawBatch), decreasing=FALSE)]
Batch.sorted<-Batch[,mixedorder(sampleNames(Batch), decreasing=FALSE)]
exprs.Batch.sorted <- exprs(Batch.sorted)
se.exprs.Batch.sorted <- se.exprs(Batch.sorted)
colnames(exprs.Batch.sorted)<-substr(colnames(Batch.sorted), 1, nchar(colnames(Batch.sorted))-2)
colnames(se.exprs.Batch.sorted)<-substr(colnames(Batch.sorted), 1, nchar(colnames(Batch.sorted))-2)
avg.exprs.Batch<-avearrays(exprs.Batch.sorted)
se.avg.exprs.Batch<-avearrays(se.exprs.Batch.sorted)
avg.Batch<-new('LumiBatch', exprs = avg.exprs.Batch, se.exprs = se.avg.exprs.Batch)


tt <-  lumiExpresso(avg.Batch, bg.correct = FALSE, variance.stabilize = TRUE, 
                           varianceStabilize.param = list(method='log2'), normalize = TRUE, )
normalizedBatch.sorted <- lumiN(tt, method='rsn')
# rawBatch.sorted<-rawBatch[,mixedorder(sampleNames(rawBatch), decreasing=FALSE)]
# save(rawBatch.sorted, file = paste(RData,"rawBatch.sorted.RData",sep=""))
# normalizedBatch.sorted<-normalizedBatch[,mixedorder(sampleNames(normalizedBatch), decreasing=FALSE)]
# save(normalizedBatch.sorted, file = paste(RData,"normalizedBatch.sorted.RData",sep=""))
# rm(rawBatch, normalizedBatch)

#------------------------------------------------------------
# QC plots
#------------------------------------------------------------
boxplot=TRUE; density=TRUE; sampleRelation=TRUE; cv=FALSE;
normalization.m = "quantile"
NumOfSample <- c(1,46,47,92,93,144)

experimentFactor <- as.factor(substr(colnames(normalizedBatch.sorted), 1, nchar(colnames(normalizedBatch.sorted))-2))
colList          <- colorsByFactor(experimentFactor)
plotColors       <- colList$plotColors
legendColors     <- colList$legendColors
rm(colList)

# 1. create symbol set for the array groups
plotSymbols <- 18-as.numeric(experimentFactor)
legendSymbols <- sort(unique(plotSymbols), decreasing=TRUE)

# 2. define display parameters for the images  		                   
WIDTH <- 1000; HEIGHT <- 1414; POINTSIZE <- 24
if(!exists("MAXARRAY")) MAXARRAY <- 55

# 3. boxplot
if(boxplot) {
  print ("Plot boxplot for intensities") 
  for (i in 1:as.integer(length(NumOfSample)/2)){
    rawFigName=paste(Figs, "boxplot_raw_", NumOfSample[2*i-1], "_", NumOfSample[2*i], sep="")
    normFigName=paste(Figs, "boxplot_normalized_", NumOfSample[2*i-1], "_", NumOfSample[2*i], sep="")
    box.plot(rawBatch.sorted[,NumOfSample[2*i-1]:NumOfSample[2*i]], rawFigName, col=plotColors, maxArray=MAXARRAY)
    box.plot(normalizedBatch.sorted[,NumOfSample[2*i-1]:NumOfSample[2*i]], normFigName, col=plotColors, maxArray=MAXARRAY)
  }
  
  box.plot(rawBatch.sorted, paste(Figs, "boxplot_raw_1-144", sep=""), col=plotColors, maxArray=144)
  box.plot(normalizedBatch.sorted, paste(Figs, "boxplot_normalized_1-144", sep=""), col=plotColors, maxArray=144)
}

# 4. density
if(density) {
  print ("Plot density histogram for intensities")
  
  # raw and normalized data
  for (i in 1:as.integer(length(NumOfSample)/2)){
    rawFigName=paste(Figs, "densityplot_raw_", NumOfSample[2*i-1], "_", NumOfSample[2*i], sep="")
    normFigName=paste(Figs, "densityplot_normalized_", NumOfSample[2*i-1], "_", NumOfSample[2*i], sep="")
    density.plot(rawBatch.sorted[,NumOfSample[2*i-1]:NumOfSample[2*i]], rawFigName, col=plotColors, maxArray=MAXARRAY)
    density.plot(normalizedBatch.sorted[,NumOfSample[2*i-1]:NumOfSample[2*i]], normFigName, col=plotColors, maxArray=MAXARRAY)
  }
  # all the raw and normalized data
  density.plot(rawBatch.sorted, paste(Figs, "densityplot_raw_1-144", sep=""), col=plotColors, maxArray=144)
  density.plot(normalizedBatch.sorted, paste(Figs, "densityplot_normalized_1-144", sep=""), col=plotColors, maxArray=144)
}

# 5. coefficient of variation
if(cv) {
  print ("Plot density for coefficient of varience for intensities")
  for (i in 1:as.integer(length(NumOfSample)/2)){
    rawFigName=paste(Figs, "cvplot_raw_", NumOfSample[2*i-1], "_", NumOfSample[2*i], sep="")
    normFigName=paste(Figs, "cvplot_normalized_", NumOfSample[2*i-1], "_", NumOfSample[2*i], sep="")
    cv.plot(rawBatch.sorted[,NumOfSample[2*i-1]:NumOfSample[2*i]], rawFigName, col=plotColors, maxArray=MAXARRAY)
    cv.plot(normalizedBatch.sorted[,NumOfSample[2*i-1]:NumOfSample[2*i]], normFigName, col=plotColors, maxArray=MAXARRAY)
  }
  cv.plot(rawBatch.sorted, paste(Figs, "cvplot_raw_1_144", sep=""), col=plotColors, maxArray=144)
  cv.plot(normalizedBatch.sorted, paste(Figs, "cvplot_normalized_1_144", sep=""), col=plotColors, maxArray=144)
}

# 6. sample relation
if(sampleRelation) {
  print ("Hierarchical clustering of data") 
  clusterOption1 = "spearman"; clusterOption2 = "complete"
  
  # raw data
  for (i in 1:as.integer(length(NumOfSample)/2)){
    clusterFun(rawBatch.sorted, range1=NumOfSample[2*i-1], range2=NumOfSample[2*i], normalized=FALSE, 
               experimentFactor=experimentFactor, clusterOption1=clusterOption1, clusterOption2=clusterOption2,
               plotColors=plotColors, legendColors=legendColors, plotSymbols=plotSymbols, legendSymbols=legendSymbols,
               WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=MAXARRAY) 
  }
  
  clusterFun(rawBatch.sorted, range1=1, range2=ncol(rawBatch.sorted), normalized=FALSE, experimentFactor=experimentFactor,
             clusterOption1=clusterOption1, clusterOption2=clusterOption2,
             plotColors=plotColors, legendColors=legendColors,
             plotSymbols=plotSymbols, legendSymbols=legendSymbols,
             WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=12,MAXARRAY=144) 
  
  # normalized data
  for (i in 1:as.integer(length(NumOfSample)/2)){
    clusterFun(normalizedBatch.sorted, range1=NumOfSample[2*i-1], range2=NumOfSample[2*i], normalized=TRUE,  
               experimentFactor=experimentFactor, clusterOption1=clusterOption1, clusterOption2=clusterOption2,
               plotColors=plotColors, legendColors=legendColors,plotSymbols=plotSymbols, legendSymbols=legendSymbols,
               WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=MAXARRAY) 
  }
  
  clusterFun(normalizedBatch.sorted, range1=1, range2=ncol(normalizedBatch.sorted), normalized=TRUE,  
             experimentFactor=experimentFactor,clusterOption1=clusterOption1, clusterOption2=clusterOption2,
             plotColors=plotColors, legendColors=legendColors, plotSymbols=plotSymbols, legendSymbols=legendSymbols,
             WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=12,MAXARRAY=144) 
  
  # replicates are averaged for normalized data
  uniqueNormalData<-normalizedBatch.sorted[,!duplicated(substr(colnames(normalizedBatch.sorted),1,nchar(colnames(normalizedBatch.sorted))-2))]
  clusterFun(uniqueNormalData, range1=1, range2=ncol(uniqueNormalData), normalized=TRUE,  experimentFactor=experimentFactor,
             clusterOption1=clusterOption1, clusterOption2=clusterOption2,plotColors=plotColors, legendColors=legendColors,
             plotSymbols=plotSymbols, legendSymbols=legendSymbols,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=MAXARRAY) 
}

rm(rawBatch, sampleInfo, uniqueNormalData)

#------------------------------------------------------------
# data filtering
#------------------------------------------------------------


# parallization
library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

filtering=TRUE

if(filtering) {
  filter.Th = 0.01; filter.dp = 0
  # will save a filtered.normData file in the working dir.
  eset.filtered.normalizedData = filterFunMatrix(rawBatch.sorted, normalizedBatch.sorted, filter.Th, filter.dp)
}
stopCluster(cl)
#if(filtering) {
#  filter.Th = 0.01; filter.dp = 0
  # will save a filtered.normData file in the working dir.
#  eset.filtered.normalizedData = filterFun(rawBatch.sorted, normalizedBatch.sorted, filter.Th, filter.dp)
#}


# remove probes that do not have annotated genes --------------------------
if (require(lumiHumanAll.db) & require(annotate)) {
  normalizedDataMatrix <- eset.filtered.normalizedData[!is.na(getSYMBOL(rownames(eset.filtered.normalizedData), 'lumiHumanAll.db')),]
}

probeList <- rownames(normalizedDataMatrix)
if (require(lumiHumanAll.db) & require(annotate)){
  geneSymbol<- getSYMBOL(probeList, 'lumiHumanAll.db')
  geneName<- sapply(lookUp(probeList, 'lumiHumanAll.db', 'GENENAME'), function(x) x[1])
}

normalizedDataMatrix.sorted<-normalizedDataMatrix[,mixedorder(colnames(normalizedDataMatrix), decreasing=FALSE)]
save(normalizedDataMatrix.sorted, file = paste(RData,"collapsed_normalizedDataMatrix.sorted.RData",sep=""))
rm(normalizedDataMatrix, eset.filtered.normalizedData)

#------------------------------------------------------------------
# remove group
#------------------------------------------------------------------


#------------------------------------------------------------------
# limma: differential expression analysis
#------------------------------------------------------------------

# create a design matrix
targets<-readTargets(paste(PhenotypeDir,"target.txt", sep=''))

ct<-factor(targets$Treatment) 
design <- model.matrix(~0+ct)
colnames(design) <- levels(ct)

fit<-lmFit(normalizedDataMatrix.sorted, design) # linear model fit for each gene
contrast.matrix <- makeContrasts(HG_SG = glucose-control, levels=design) # create a contrast matrix    
fit2<-contrasts.fit(fit, contrast.matrix) # given a linear model fit to the data, compute estimated coefficients and std errors
fit2<-eBayes(fit2); #  Given a linear model fit, compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression 
#  by empirical Bayes moderation
fit2$genes <-data.frame(ID=probeList, geneSymbol=geneSymbol, geneName=geneName, stringsAsFactors=FALSE)

#  Volcano plot for High Glucose vs. Standard Glucose

pval.cutoff<-0.01
# HG_SG <- topTable(fit2, coef="HG_SG", number=nrow(fit2$genes), sort.by="logFC", genelist=fit2$genes, adjust.method="BH")
HG_SG <- topTable(fit2, coef="HG_SG", number=nrow(fit2$genes), sort.by="logFC", resort.by="p", genelist=fit2$genes, adjust.method="BH")
pdf(file = paste(Figs,"HG_SG_Pval0.01.pdf",sep=""),width=18,height=12,pointsize=20)
with(HG_SG, plot(logFC, -log10(P.Value), pch=20, col="gray", main="HG vs SG (P.val <0.01)", xlab="log2 fold change", ylab="-log10(p-value)", xlim=c(-0.5,0.5)))
with(subset(HG_SG, P.Value<pval.cutoff), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(HG_SG, P.Value<pval.cutoff), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()

write.csv(HG_SG, paste(Output,"HG_SG.csv",sep=""))

#------------------------------------------------------------
#  calculate delta and mean delta
#------------------------------------------------------------

if (REUSE){
  load(file = paste(RData,"normalizedDataMatrix.sorted.RData",sep=""))
}

if (REUSE == FALSE){
  HighGlucoseSamples <-c(seq(1,5),seq(11,13), seq(17,19), seq(23,25), seq(29,31), seq(35,37), seq(41,43), seq(47,49), 
                         seq(53,57), seq(63,65), seq(69,71), seq(75,77), seq(81,83), seq(87,89), seq(93,97), seq(103,105),
                         seq(109,111), seq(115,117), seq(121,123), seq(127,129), seq(133,135), seq(139,141))
  
  NoGlucoseSamples <-c(seq(6,10),seq(14,16), seq(20,22), seq(26,28), seq(32,34), seq(38,40), seq(44,46), seq(50,52), 
                       seq(58,62), seq(66,68), seq(72,74), seq(78,80), seq(84,86), seq(90,92), seq(98,102), seq(106,108),
                       seq(112,114), seq(118,120), seq(124,126), seq(130,132), seq(136,138), seq(142,144))
  
  HighGlucoseData<-normalizedDataMatrix.sorted[,HighGlucoseSamples]
  NoGlucoseData<-normalizedDataMatrix.sorted[,NoGlucoseSamples]
  deltaData<-sweep(HighGlucoseData,1,NoGlucoseData,"-") # delta x
  save(deltaData, file = paste(RData,"deltaData.RData",sep=""))
  write.csv(deltaData, file= paste(Output,"deltaData.csv", sep=""))
  deltaData.T<-data.frame(t(deltaData))
  deltaData.T$Name<-substr(colnames(deltaData),1,nchar(colnames(deltaData))-2)
  meanDeltaData.T<-aggregate(deltaData.T[,-ncol(deltaData.T)],list(deltaData.T$Name), mean)
  
  group1=seq(1,nrow(meanDeltaData.T))
  group2=c(rep(0,8),rep(1,7),rep(2,7))
  groupName=substr(meanDeltaData.T$Group.1,1,nchar(meanDeltaData.T$Group)-5); 
  rownames(meanDeltaData.T)<-groupName
  meanDeltaData.T$Group.1<-NULL
  
  meanDeltaData <- t(meanDeltaData.T)
  rownames(meanDeltaData)<-rownames(deltaData)
  save(meanDeltaData, file = paste(RData,"meanDeltaData.RData",sep=""))
  write.csv(meanDeltaData, file= paste(Output,"meanDeltaData.csv", sep=""))
  rm(deltaData.T, meanDeltaData.T, normalizedDataMatrix, normalizedDataMatrix.sorted, 
     eset.filtered.normalizedData, HighGlucoseData, NoGlucoseData, Batch, geneName, geneSymbol,
     normalizedBatch.sorted, probeList, rawBatch.sorted)
}

if (REUSE){
  load(file = paste(RData,"deltaData.RData",sep=""))
  load(file = paste(RData,"meanDeltaData.RData",sep=""))
}

#--------------------------------------------------------------------------
#  Todo: export gene expression data for the genes in Mike's plos one paper
#  IL-1B (IL1B) fold change = 2.11, pvalue = 0.02 
#  PKCB (PRKCB) fold change = 2.30, p-value = 0.01
#  NFKB-p50 (NFKB1) fold change = 2.05, p-value = 0.01
#  NFKB-p65 (RELA) fold change = 2.82, p-value = 0.003
#  CD18 (ITGB2) fold change = 2.59, p-value = 0.02
#--------------------------------------------------------------------------

MG_geneList = c("IL1B", "PRKCB", "NFKB1", "RELA", "ITGB2")

# delta data --------------------------------------------------------------
deltaData.genesymbol<-deltaData[rownames(HG_SG),]
colnames(deltaData.genesymbol)<-gsub("_30mM", "", colnames(deltaData.genesymbol))
rownames(deltaData.genesymbol)<-HG_SG$geneSymbol
write.csv(deltaData.genesymbol, file= paste(Output,"deltaData.symbol.csv", sep=""))
deltaData.genesymbol<-2^deltaData.genesymbol
output1<-matrix(, ,ncol(deltaData.genesymbol))
for (i in 1:length(MG_geneList)) {
  tmp=deltaData.genesymbol[which(MG_geneList[i]==rownames(deltaData.genesymbol)), , drop=FALSE]
  print(nrow(tmp))
  output1<-rbind(output1,tmp)
}
output1 <- na.omit(output1)
write.csv(output1, file= paste(Output,"deltaData.5genes.csv", sep=""))

# meandelta data --------------------------------------------------------------
meanDeltaData.genesymbol<-meanDeltaData[rownames(HG_SG),]
colnames(meanDeltaData.genesymbol)<-gsub("_30mM", "", colnames(meanDeltaData.genesymbol))
rownames(meanDeltaData.genesymbol)<-HG_SG$geneSymbol
write.csv(meanDeltaData.genesymbol, file= paste(Output,"meanDeltaData.symbol.csv", sep=""))
meanDeltaData.genesymbol<-2^meanDeltaData.genesymbol
output2<-matrix(, ,ncol(meanDeltaData.genesymbol))
for (i in 1:length(MG_geneList)) {
  tmp=meanDeltaData.genesymbol[which(MG_geneList[i]==rownames(meanDeltaData.genesymbol)), , drop=FALSE]
  output2<-rbind(output2,tmp)
}
output2 <- na.omit(output2)
write.csv(output2, file= paste(Output,"meanDeltaData.5genes.csv", sep=""))

rm(deltaData, meanDeltaData)
#--------------------------------------------------------------------------
#  find genes with > 2 in any samples
#--------------------------------------------------------------------------
deltaData.genesymbol<-data.frame(deltaData.genesymbol)
meanDeltaData.genesymbol<-data.frame(meanDeltaData.genesymbol)

deltaData.genesymbol$max<-apply(deltaData.genesymbol,1,max)
deltaData.genesymbol$min<-apply(deltaData.genesymbol,1,min)

meanDeltaData.genesymbol$max<-apply(meanDeltaData.genesymbol,1,max)
meanDeltaData.genesymbol$min<-apply(meanDeltaData.genesymbol,1,min)

delta.<-deltaData.genesymbol[which(deltaData.genesymbol$max >= 2),]


#-------------------------------
#  STOP HERE
#-------------------------------




library(foreach)
library(doParallel)
library(ggplot2)
library(limma)
library(gtools); library(bioDist); library(calibrate)
library(plyr); library(reshape2); library(scales)
library(ggfortify)

HOME <- "/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/"  # change if needed
InputDir <- paste(HOME,'data/input/', sep='')
InputFiles = list.files(path=InputDir, pattern='txt')
ControlDir <- paste(HOME,'data/control/', sep='')
ControlFiles = list.files(path=ControlDir, pattern='txt')
PhenotypeDir = paste(HOME,'data/phenotype/', sep='')
RData=paste(HOME,'rdata/',sep=''); dir.create(RData, showWarnings = FALSE)
Figs = paste(HOME,'Figs/', sep=''); dir.create(Figs, showWarnings = FALSE)
Output = paste(HOME,'output/', sep=''); dir.create(Figs, showWarnings = FALSE)

source(paste(HOME,"pca/pca.R",sep=""))
source(paste(HOME,"script/Function.R",sep=""))

date=Sys.Date()

# load data
dat<-local(get(load(file=paste(RData,"normalizedDataMatrix_filtered.RData",sep=""))))
# gene_annotation<-read.csv(paste(PhenotypeDir,"raw_delta_gene_anno.csv",sep=""))
# rownames(gene_annotation)<-gene_annotation$GeneCode

probeList <- rownames(dat)  # probeID (nuID) gene annotation
if (require(lumiHumanAll.db) & require(annotate)){
  geneSymbol<- getSYMBOL(probeList, 'lumiHumanAll.db')
  geneName<- sapply(lookUp(probeList, 'lumiHumanAll.db', 'GENENAME'), function(x) x[1])
}
genes <-data.frame(ID=probeList, geneSymbol=geneSymbol, geneName=geneName, stringsAsFactors=FALSE)

##----------------------------------
## group comparison (each hg and sg)
##----------------------------------
#1a. No_DM vs PDR + No_PDR
# subject="No_DM_vs_Diabetes"
subject="all_replicate"
targets<-readTargets(paste(PhenotypeDir,"hg_sg_",subject,"_target.txt", sep=''))
group <- factor(paste(targets$Group,targets$Treatment,sep="."))
design <- model.matrix(~0+group)
colnames(design)<-levels(group)
rownames(design)<-targets$SampleName
corfit <- duplicateCorrelation(dat, block = targets$Subject)
fit <- lmFit(dat, design, block = targets$Subject, cor = corfit$consensus.correlation)

cm <- makeContrasts( hg_NoD_DwoC=NoD.T - DwoC.T, hg_NoD_DwC=NoD.T - DwC.T, 
                     hg_NoD_Diabetes=NoD.T - (DwC.T+DwoC.T)/2, hg_DwC_DwoC = DwC.T-DwoC.T,
                     sg_NoD_DwoC=NoD.C - DwoC.C, sg_NoD_DwC=NoD.C - DwC.C,
                     sg_NoD_Diabetes=NoD.C - (DwC.C+DwoC.C)/2, sg_DwC_DwoC = DwC.C-DwoC.C,
                     levels=design) # create a contrast matrix
fit2=contrasts.fit(fit,cm)
fit2=eBayes(fit2)

setEPS()
postscript(file = paste(Figs,"hg_sg_all_group_pvalues_",date,".eps",sep=""))

par(mfrow=c(2,4))
for (i in 1:ncol(fit2$p.value)) {
  hist(fit2$p.value[,i], main=colnames(fit2$p.value)[i])
}
dev.off()

#----------------
# hg NoD vs DwoC
#----------------

qval.cutoff=0.05; FC.cutoff=1 # FC=2

x1=topTable(fit2, coef="hg_NoD_DwoC", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="hg_NoD_DwoC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"hg_NoD_DwoC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
# pdf(file = paste(Figs,"hg_NoD_DwoC_q",qval.cutoff,"_","FC",FC.cutoff,"_",date,".pdf",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in high glucose (No_DM vs. No_PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"hg_NoD_DwoC_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"hg_NoD_DwoC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#----------------
# hg NoD vs DwC
#----------------
x1=topTable(fit2, coef="hg_NoD_DwC", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="hg_NoD_DwC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"hg_NoD_DwC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
# pdf(file = paste(Figs,"hg_NoD_DwC_q",qval.cutoff,"_","FC",FC.cutoff,"_",date,".pdf",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in high glucose (No_DM vs. PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"hg_NoD_DwC_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"hg_NoD_DwC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#---------------------
# hg NoD vs DwC+DwoC
#---------------------
x1=topTable(fit2, coef="hg_NoD_Diabetes", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="hg_NoD_Diabetes", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"hg_NoD_Diabetes_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
# pdf(file = paste(Figs,"hg_NoD_Diabetes_q",qval.cutoff,"_","FC",FC.cutoff,"_",date,".pdf",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in high glucose (No_DM vs. PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"hg_NoD_Diabetes_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"hg_NoD_Diabetes_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#---------------------
# hg DwC vs DwoC
#---------------------
x1=topTable(fit2, coef="hg_DwC_DwoC", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="hg_DwC_DwoC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"hg_DwC_DwoC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
# pdf(file = paste(Figs,"hg_DwC_DwoC_q",qval.cutoff,"_","FC",FC.cutoff,"_",date,".pdf",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in high glucose (PDR vs. No_PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"hg_DwC_DwoC_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"hg_DwC_DwoC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#----------------
# sg NoD vs DwoC
#----------------

qval.cutoff=0.05; FC.cutoff=1 # FC=1.12

x1=topTable(fit2, coef="sg_NoD_DwoC", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="sg_NoD_DwoC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"sg_NoD_DwoC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
# pdf(file = paste(Figs,"sg_NoD_DwoC_q",qval.cutoff,"_","FC",FC.cutoff,"_",date,".pdf",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in standard glucose (No_DM vs. No_PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"sg_NoD_DwoC_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"sg_NoD_DwoC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#----------------
# sg NoD vs DwC
#----------------
x1=topTable(fit2, coef="sg_NoD_DwC", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="sg_NoD_DwC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"sg_NoD_DwC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
# pdf(file = paste(Figs,"sg_NoD_DwC_q",qval.cutoff,"_","FC",FC.cutoff,"_",date,".pdf",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in standard glucose (No_DM vs. PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"sg_NoD_DwC_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"sg_NoD_DwC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#-------------------
# sg NoD vs DwC+DwoC
#-------------------
x1=topTable(fit2, coef="sg_NoD_Diabetes", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="sg_NoD_Diabetes", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"sg_NoD_Diabetes_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
# pdf(file = paste(Figs,"sg_NoD_Diabetes_q",qval.cutoff,"_","FC",FC.cutoff,"_",date,".pdf",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in standard glucose (No_DM vs. No_PDR + PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"sg_NoD_Diabetes_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"sg_NoD_Diabetes_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#---------------------
# sg DwC vs DwoC
#---------------------
x1=topTable(fit2, coef="sg_DwC_DwoC", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="sg_DwC_DwoC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"sg_DwC_DwoC_q",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
# pdf(file = paste(Figs,"sg_DwC_DwoC_q",qval.cutoff,"_","FC",FC.cutoff,"_",date,".pdf",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in standard glucose (PDR vs. No_PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"sg_DwC_DwoC_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"sg_DwC_DwoC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)


##--------------------------------------------------
## treatment response high vs standard for each group
##--------------------------------------------------
#-----------------
# 1. No_DM (NoD)
#-----------------
subject="No_DM_replicate"
NoD<-dat[,grep("NoD",colnames(dat))] # grep all normal samples
targets<-readTargets(paste(PhenotypeDir,"hg_sg_",subject,"_target.txt", sep=''))
Treat <- factor(targets$Treatment,levels=c("C","T"))
Replicates <- factor(targets$rep)
design <- model.matrix(~Replicates+Treat)

corfit <- duplicateCorrelation(NoD, block = targets$Subject)
fit <-lmFit(NoD,design,block=targets$Subject,correlation=corfit$consensus.correlation)
fit<-eBayes(fit)

qval.cutoff=0.1; FC.cutoff=0.17 # FC=1.12

x1=topTable(fit, coef="TreatT", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="TreatT", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"NoD_treatment_effect_qval_",qval.cutoff,"_", date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("No_DM treatment effect (high vs. standard glucose)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-.7,.7)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

cat("number of genes (q-val <0.05): ",sum(x1$adj.P.Val<qval.cutoff))
write.csv(x1_subset, paste(Output,"NoD_treatment_effect_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1, paste(Output,"NoD_treatment_effect_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(y1, paste(Output,"NoD_treatment_effect_",date,".csv",sep=""),row.names = F)

setEPS()
postscript(file = paste(Figs,"NoD_treatment_effect_pvalue_distribution_",date,".eps",sep=""),width=18,height=12,pointsize=20)
par(mfrow=c(1,2))
hist(y1$P.Value, main="",xlab="p-value")
hist(y1$adj.P.Val,main="", xlab="adjusted p-value")
mtext("No_DM treatment effect (high vs. standard glucose)", side=3, outer=TRUE, line=-3)
dev.off()

# 2. No_PDR (DwoC)
subject="No_PDR_replicate"
DwoC<-dat[,grep("DwoC",colnames(dat))] # grep all normal samples
targets<-readTargets(paste(PhenotypeDir,"hg_sg_",subject,"_target.txt", sep=''))
Treat <- factor(targets$Treatment,levels=c("C","T"))
Replicates <- factor(targets$rep)
design <- model.matrix(~Replicates+Treat)

corfit <- duplicateCorrelation(DwoC, block = targets$Subject)
fit <-lmFit(DwoC,design,block=targets$Subject,correlation=corfit$consensus.correlation)
fit<-eBayes(fit)

qval.cutoff=0.1; FC.cutoff=0.17 # FC=1.12

x1=topTable(fit, coef="TreatT", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="TreatT", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"DwoC_treatment_effect_qval_",qval.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("No_PDR treatment effect (high vs. standard glucose)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-.5,.5)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

cat("number of genes (q-val <0.05): ",sum(x1$adj.P.Val<qval.cutoff))
write.csv(x1_subset, paste(Output,"DwoC_treatment_effect_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1, paste(Output,"DwoC_treatment_effect_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(y1, paste(Output,"DwoC_treatment_effect_",date,".csv",sep=""),row.names = F)

setEPS()
postscript(file = paste(Figs,"DwoC_treatment_effect_pvalue_distribution_",date,".eps",sep=""),width=18,height=12,pointsize=20)
par(mfrow=c(1,2))
hist(y1$P.Value, main="",xlab="p-value")
hist(y1$adj.P.Val,main="", xlab="adjusted p-value")
mtext("No_PDR treatment effect (high vs. standard glucose)", side=3, outer=TRUE, line=-3)
dev.off()

# 3. PDR (DwC)
subject="PDR_replicate"
DwC<-dat[,grep("DwC",colnames(dat))] # grep all normal samples
targets<-readTargets(paste(PhenotypeDir,"hg_sg_",subject,"_target.txt", sep=''))
Treat <- factor(targets$Treatment,levels=c("C","T"))
Replicates <- factor(targets$rep)
design <- model.matrix(~Replicates+Treat)

corfit <- duplicateCorrelation(DwC, block = targets$Subject)
fit <-lmFit(DwC,design,block=targets$Subject,correlation=corfit$consensus.correlation)
fit<-eBayes(fit)

qval.cutoff=0.1; FC.cutoff=0.17 # FC=1.12

x1=topTable(fit, coef="TreatT", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="TreatT", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"DwC_treatment_effect_q",qval.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("PDR treatment effect (high vs. standard glucose)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-.5,.5)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

cat("number of genes (q-val <0.1): ",sum(x1$adj.P.Val<qval.cutoff))
write.csv(x1_subset, paste(Output,"DwC_treatment_effect_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1, paste(Output,"DwC_treatment_effect_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(y1, paste(Output,"DwC_treatment_effect_",date,".csv",sep=""),row.names = F)

setEPS()
postscript(file = paste(Figs,"DwC_treatment_effect_pvalue_distribution_",date,".eps",sep=""),width=18,height=12,pointsize=20)
par(mfrow=c(1,2))
hist(y1$P.Value, main="",xlab="p-value")
hist(y1$adj.P.Val,main="", xlab="adjusted p-value")
mtext("PDR treatment effect (high vs. standard glucose)", side=3, outer=TRUE, line=-3)
dev.off()

##---------------------
##    Diabetes
##---------------------
subject="Diabetes"
Diabetes=cbind(DwC,DwoC)
targets<-readTargets(paste(PhenotypeDir,"hg_sg_",subject,"_target.txt", sep=''))
Treat <- factor(targets$Treatment,levels=c("C","T"))
Replicates <- factor(targets$rep)
design <- model.matrix(~Replicates+Treat)

corfit <- duplicateCorrelation(Diabetes, block = targets$Subject)
fit <-lmFit(Diabetes,design,block=targets$Subject,correlation=corfit$consensus.correlation)
fit<-eBayes(fit)

qval.cutoff=0.1; FC.cutoff=0.17 # FC=1.12

x1=topTable(fit, coef="TreatT", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="TreatT", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"Diabetes_treatment_effect_qval_",qval.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("No_PDR and PDR treatment effect (high vs. standard glucose)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-0.5,0.5)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

cat("number of genes (q-val <0.05): ",sum(x1$adj.P.Val<qval.cutoff))
write.csv(x1_subset, paste(Output,"Diabetes_treatment_effect_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1, paste(Output,"Diabetes_treatment_effect_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(y1, paste(Output,"Diabetes_treatment_effect_",date,".csv",sep=""),row.names = F)

setEPS()
postscript(file = paste(Figs,"Diabetes_treatment_effect_pvalue_distribution_",date,".eps",sep=""),width=18,height=12,pointsize=20)
par(mfrow=c(1,2))
hist(y1$P.Value, main="",xlab="p-value")
hist(y1$adj.P.Val,main="", xlab="adjusted p-value")
mtext("No_PDR and PDR treatment effect (high vs. standard glucose)", side=3, outer=TRUE, line=-3)
dev.off()

##---------------------------------------
##    all group
##---------------------------------------
subject="all_replicate"
targets<-readTargets(paste(PhenotypeDir,"hg_sg_",subject,"_target.txt", sep=''))
Treat <- factor(targets$Treatment,levels=c("C","T"))
Replicates <- factor(targets$rep)
design <- model.matrix(~Replicates+Treat)

corfit <- duplicateCorrelation(dat, block = targets$Subject)
fit <-lmFit(dat,design,block=targets$Subject,correlation=corfit$consensus.correlation)
fit<-eBayes(fit)

qval.cutoff=0.05; FC.cutoff=0.17 # FC=1.12

# x1=topTable(fit, coef="TreatT", n=nrow(genes), p.value=qval.cutoff,genelist=genes)
tmp=topTable(fit, coef="TreatT", n=nrow(genes),adjust.method="BH",genelist=genes)
x1=tmp[tmp$P.Value<qval.cutoff,]
y1=topTable(fit, coef="TreatT", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$P.Value<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"All_treatment_effect_pval_",qval.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(P.Value), pch=20, col="gray", main=paste("No_DM, No_PDR and PDR treatment effect (high vs. standard glucose)",sep=""), xlab="log2 fold change", ylab="-log10(P.Value)", xlim=c(-0.5,0.5)))
with(subset(y1, P.Value<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(y1, P.Value<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()

cat("number of genes (p-val <0.05): ",sum(x1$P.Value<qval.cutoff))
write.csv(x1_subset, paste(Output,"All_treatment_effect_pval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1, paste(Output,"All_treatment_effect_pval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(y1, paste(Output,"All_treatment_effect_",date,".csv",sep=""),row.names = F)

setEPS()
postscript(file = paste(Figs,"All_treatment_effect_pvalue_distribution_",date,".eps",sep=""),width=18,height=12,pointsize=20)
# par(mfrow=c(1,2))
hist(y1$P.Value, main="",xlab="p-value")
mtext("No_DM, No_PDR and PDR treatment effect (high vs. standard glucose)", side=3, outer=TRUE, line=-3)
dev.off()


# 
# 
# #---------------------------------------------------------------------------
# # PCA by group (HG, SG, delta)
# # 1. No_DM
# subject<-"No_DM"
# NoD_sg<-NoD[,grep("norm",colnames(NoD))]; NoD_hg<-NoD[,grep("30mM",colnames(NoD))];
# NoD_delta <- NoD_hg-NoD_sg
# pca_NoD_sg <- prcomp(NoD_sg); pca_NoD_hg <- prcomp(NoD_hg); pca_NoD_delta <- prcomp(NoD_delta)
# 
# # 1a. variance explained
# pdf(file = paste(Figs,subject,"_sg_pca_variance_explained_",date,".pdf",sep=""))
# plot_variance_explained(pca_NoD_sg) + theme_bw(base_size = 18)
# dev.off()
# 
# pdf(file = paste(Figs,subject,"_hg_pca_variance_explained_",date,".pdf",sep=""))
# plot_variance_explained(pca_NoD_hg) + theme_bw(base_size = 18)
# dev.off()
# 
# pdf(file = paste(Figs,subject,"_delta_pca_variance_explained_",date,".pdf",sep=""))
# plot_variance_explained(pca_NoD_delta) + theme_bw(base_size = 18)
# dev.off()
# 
# # 2. No_PDR
# subject<-"No_PDR"
# DwoC_sg<-DwoC[,grep("norm",colnames(DwoC))]; DwoC_hg<-DwoC[,grep("30mM",colnames(DwoC))];
# DwoC_delta <- DwoC_hg-DwoC_sg
# pca_DwoC_sg <- prcomp(DwoC_sg); pca_DwoC_hg <- prcomp(DwoC_hg); pca_DwoC_delta <- prcomp(DwoC_delta)
# 
# # 2a. variance explained
# pdf(file = paste(Figs,subject,"_sg_pca_variance_explained_",date,".pdf",sep=""))
# plot_variance_explained(pca_DwoC_sg) + theme_bw(base_size = 18)
# dev.off()
# 
# pdf(file = paste(Figs,subject,"_hg_pca_variance_explained_",date,".pdf",sep=""))
# plot_variance_explained(pca_DwoC_hg) + theme_bw(base_size = 18)
# dev.off()
# 
# pdf(file = paste(Figs,subject,"_delta_pca_variance_explained_",date,".pdf",sep=""))
# plot_variance_explained(pca_DwoC_delta) + theme_bw(base_size = 18)
# dev.off()
# 
# 
# # 3. PDR
# subject<-"PDR"
# DwC_sg<-DwC[,grep("norm",colnames(DwC))]; DwC_hg<-DwC[,grep("30mM",colnames(DwC))];
# DwC_delta <- DwC_hg-DwC_sg
# pca_DwC_sg <- prcomp(DwC_sg); pca_DwC_hg <- prcomp(DwC_hg); pca_DwC_delta <- prcomp(DwC_delta)
# 
# # 3a. variance explained
# pdf(file = paste(Figs,subject,"_sg_pca_variance_explained_",date,".pdf",sep=""))
# plot_variance_explained(pca_DwC_sg) + theme_bw(base_size = 18)
# dev.off()
# 
# pdf(file = paste(Figs,subject,"_hg_pca_variance_explained_",date,".pdf",sep=""))
# plot_variance_explained(pca_DwC_hg) + theme_bw(base_size = 18)
# dev.off()
# 
# pdf(file = paste(Figs,subject,"_delta_pca_variance_explained_",date,".pdf",sep=""))
# plot_variance_explained(pca_DwC_delta) + theme_bw(base_size = 18)
# dev.off()
# 
# # 4. data exploration
# data_explore_by_pca1(NoD_hg)
# data_explore_by_pca1(NoD_sg)
# data_explore_by_pca1(NoD_delta)
# data_explore_by_pca1(NoD)
# 
# data_explore_by_pca1(DwC_hg)
# data_explore_by_pca1(DwC_sg)
# data_explore_by_pca1(DwC_delta)
# data_explore_by_pca1(DwC)
# 
# data_explore_by_pca1(DwoC_hg)
# data_explore_by_pca1(DwoC_sg)
# data_explore_by_pca1(DwoC_delta)
# data_explore_by_pca1(DwoC)
# 
# #---------------------------------------------------------------------------
# # PCA and covariates (HG, SG, delta)
# # 1. find confounding factor using collasped data
# meta<-read.table(paste(PhenotypeDir, "meanDeltaSampleInfo3.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
# 
# # 1a. HG - collapse data
# Diabete_HG<-cbind(DwC_hg, DwoC_hg)
# Diabete_HG.T<-data.frame(t(Diabete_HG))
# Diabete_HG.T$Name<-substr(rownames(Diabete_HG.T),1,nchar(rownames(Diabete_HG.T))-7)
# Diabete_HG.T.agg<-aggregate(Diabete_HG.T[,-ncol(Diabete_HG.T)],list(Diabete_HG.T$Name), mean)
# rownames(Diabete_HG.T.agg)<-Diabete_HG.T.agg$Group.1; Diabete_HG.T.agg$Group.1<-NULL
# tmp<-apply(Diabete_HG.T.agg,2,as.numeric); rownames(tmp)<-rownames(Diabete_HG.T.agg)
# Diabete_HG_mean<-data.frame(t(tmp))
# Diabete_HG_mean_input<-Diabete_HG_mean[,rownames(meta)]
# rm(Diabete_HG_mean, Diabete_HG.T, Diabete_HG.T.agg)
# 
# # 1b. HG - PCA variance explained
# pca_Diabete_HG_mean_input <- prcomp(Diabete_HG_mean_input, scale = FALSE, center = TRUE)
# pdf(file = paste(Figs,"Diabete_HG_mean_pca_variance_explained_",date,".pdf",sep=""),width=18)
# plot_variance_explained(pca_Diabete_HG_mean_input, 10) + theme_bw(base_size = 18)
# dev.off()
# 
# # 1c. HG - correlation between PCA and covariates
# pca_Diabete_HG_mean_input_cor<-correlate_pcs(pca_Diabete_HG_mean_input, meta, npcs = 10, min.cor = 0)
# pca_Diabete_HG_mean_input_cor<-as.data.frame(pca_Diabete_HG_mean_input_cor);
# pca_Diabete_HG_mean_input_cor<-abs(pca_Diabete_HG_mean_input_cor)
# 
# tt<-cbind(rownames(pca_Diabete_HG_mean_input_cor),pca_Diabete_HG_mean_input_cor);
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "Correlation")
# pdf(file = paste(Figs,"Diabete_HG_mean_cor_pca_covariate_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = Correlation),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(Correlation, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))
# dev.off()
# rm(pca_Diabete_HG_mean_input_cor, tt,tt.m)
# 
# # 1d. HG - proportion of variance for PCA and covariate
# pca_Diabete_HG_mean_input_r2<-proportional_pcs(pca_Diabete_HG_mean_input, meta, npcs = 10, min.cor = 0)
# pca_Diabete_HG_mean_input_r2<-as.data.frame(pca_Diabete_HG_mean_input_r2); pca_Diabete_HG_mean_input_r2<-abs(pca_Diabete_HG_mean_input_r2)
# tt<-cbind(rownames(pca_Diabete_HG_mean_input_r2),pca_Diabete_HG_mean_input_r2);
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "adjusted.R2")
# pdf(file = paste(Figs,"Diabete_HG_mean_r2_pca_covariate_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = adjusted.R2),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(adjusted.R2, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))
# dev.off()
# rm(pca_Diabete_HG_mean_input_r2, tt,tt.m)
# 
# # 2a. SG - collapse data
# Diabete_SG<-cbind(DwC_sg, DwoC_sg)
# Diabete_SG.T<-data.frame(t(Diabete_SG))
# Diabete_SG.T$Name<-substr(rownames(Diabete_SG.T),1,nchar(rownames(Diabete_SG.T))-7)
# Diabete_SG.T.agg<-aggregate(Diabete_SG.T[,-ncol(Diabete_SG.T)],list(Diabete_SG.T$Name), mean)
# rownames(Diabete_SG.T.agg)<-Diabete_SG.T.agg$Group.1; Diabete_SG.T.agg$Group.1<-NULL
# tmp<-apply(Diabete_SG.T.agg,2,as.numeric); rownames(tmp)<-rownames(Diabete_SG.T.agg)
# Diabete_SG_mean<-data.frame(t(tmp))
# Diabete_SG_mean_input<-Diabete_SG_mean[,rownames(meta)]
# rm(Diabete_SG_mean, Diabete_SG.T, Diabete_SG.T.agg)
# 
# # 2b. SG - PCA variance explained
# pca_Diabete_SG_mean_input <- prcomp(Diabete_SG_mean_input, scale = FALSE, center = TRUE)
# pdf(file = paste(Figs,"Diabete_SG_mean_pca_variance_explained_",date,".pdf",sep=""),width=18)
# plot_variance_explained(pca_Diabete_SG_mean_input, 10) + theme_bw(base_size = 18)
# dev.off()
# 
# # 2c. SG - correlation between PCA and covariates
# pca_Diabete_SG_mean_input_cor<-correlate_pcs(pca_Diabete_SG_mean_input, meta, npcs = 10, min.cor = 0)
# pca_Diabete_SG_mean_input_cor<-as.data.frame(pca_Diabete_SG_mean_input_cor);
# pca_Diabete_SG_mean_input_cor<-abs(pca_Diabete_SG_mean_input_cor)
# 
# tt<-cbind(rownames(pca_Diabete_SG_mean_input_cor),pca_Diabete_SG_mean_input_cor);
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "Correlation")
# pdf(file = paste(Figs,"Diabete_SG_mean_cor_pca_covariate_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = Correlation),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(Correlation, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))
# dev.off()
# rm(pca_Diabete_SG_mean_input_cor, tt,tt.m)
# 
# # 2d. SG - proportion of variance for PCA and covariate
# pca_Diabete_SG_mean_input_r2<-proportional_pcs(pca_Diabete_SG_mean_input, meta, npcs = 10, min.cor = 0)
# pca_Diabete_SG_mean_input_r2<-as.data.frame(pca_Diabete_SG_mean_input_r2); pca_Diabete_SG_mean_input_r2<-abs(pca_Diabete_SG_mean_input_r2)
# tt<-cbind(rownames(pca_Diabete_SG_mean_input_r2),pca_Diabete_SG_mean_input_r2);
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "adjusted.R2")
# pdf(file = paste(Figs,"Diabete_SG_mean_r2_pca_covariate_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = adjusted.R2),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(adjusted.R2, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))
# dev.off()
# rm(pca_Diabete_SG_mean_input_r2, tt,tt.m)
# 
# # 3a. delta - collapse data
# Diabete_delta<-cbind(DwC_delta, DwoC_delta)
# Diabete_delta.T<-data.frame(t(Diabete_delta))
# Diabete_delta.T$Name<-substr(rownames(Diabete_delta.T),1,nchar(rownames(Diabete_delta.T))-7)
# Diabete_delta.T.agg<-aggregate(Diabete_delta.T[,-ncol(Diabete_delta.T)],list(Diabete_delta.T$Name), mean)
# rownames(Diabete_delta.T.agg)<-Diabete_delta.T.agg$Group.1; Diabete_delta.T.agg$Group.1<-NULL
# tmp<-apply(Diabete_delta.T.agg,2,as.numeric); rownames(tmp)<-rownames(Diabete_delta.T.agg)
# Diabete_delta_mean<-data.frame(t(tmp))
# Diabete_delta_mean_input<-Diabete_delta_mean[,rownames(meta)]
# rm(Diabete_delta_mean, Diabete_delta.T, Diabete_delta.T.agg)
# 
# # 3b. delta - PCA variance explained
# pca_Diabete_delta_mean_input <- prcomp(Diabete_delta_mean_input, scale = FALSE, center = TRUE)
# pdf(file = paste(Figs,"Diabete_delta_mean_pca_variance_explained_",date,".pdf",sep=""),width=18)
# plot_variance_explained(pca_Diabete_delta_mean_input, 15) + theme_bw(base_size = 18)
# dev.off()
# 
# # 3c. delta - correlation between PCA and covariates
# pca_Diabete_delta_mean_input_cor<-correlate_pcs(pca_Diabete_delta_mean_input, meta, npcs = 15, min.cor = 0)
# pca_Diabete_delta_mean_input_cor<-as.data.frame(pca_Diabete_delta_mean_input_cor);
# pca_Diabete_delta_mean_input_cor<-abs(pca_Diabete_delta_mean_input_cor)
# 
# tt<-cbind(rownames(pca_Diabete_delta_mean_input_cor),pca_Diabete_delta_mean_input_cor);
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "Correlation")
# pdf(file = paste(Figs,"Diabete_delta_mean_cor_pca_covariate_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = Correlation),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(Correlation, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))
# dev.off()
# rm(pca_Diabete_delta_mean_input_cor, tt,tt.m)
# 
# # 3d. delta - proportion of variance for PCA and covariate
# pca_Diabete_delta_mean_input_r2<-proportional_pcs(pca_Diabete_delta_mean_input, meta, npcs = 15, min.cor = 0)
# pca_Diabete_delta_mean_input_r2<-as.data.frame(pca_Diabete_delta_mean_input_r2); pca_Diabete_delta_mean_input_r2<-abs(pca_Diabete_delta_mean_input_r2)
# tt<-cbind(rownames(pca_Diabete_delta_mean_input_r2),pca_Diabete_delta_mean_input_r2);
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "adjusted.R2")
# pdf(file = paste(Figs,"Diabete_delta_mean_r2_pca_covariate_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = adjusted.R2),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(adjusted.R2, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))
# dev.off()
# rm(pca_Diabete_delta_mean_input_r2, tt,tt.m)
# 
# rm(pca_Diabete_delta_mean_input, pca_Diabete_HG_mean_input, pca_Diabete_SG_mean_input, pca_DwC_delta, pca_DwC_hg,
#    pca_DwC_sg, pca_DwoC_delta, pca_DwoC_hg, pca_DwoC_sg, pca_tmp)
# 
# #--------------------------------------------
# # covariate correlation
# meta.cor<-correlate_object(meta); meta.cor<-abs(meta.cor)
# tt.m<-melt(meta.cor); colnames(tt.m)<-c("Covariate1", "Covariate2", "Correlation")
# 
# pdf(file = paste(Figs,"covariate_correlation_",date,".pdf",sep=""),width=14,height=7)
# ggplot(tt.m, aes(Covariate1, Covariate2)) + geom_tile(aes(fill = Correlation),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(Correlation, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                            axis.text.x = element_text(colour="grey20",size=10,face="bold",angle = 90, hjust = 1),
#                                                                                                                                                                                                            axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))
# dev.off()
# 
# rm(meta.cor, tt.m)
# 
# #-----------------------------------------------
# ## tretment effect: DE analysis on Diabetes (DwC+DwoC)
# #----------------------------------------------
# meta<-read.table(paste(PhenotypeDir, "meanDeltaSampleInfo3.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
# Diabetes<-dat[,-grep("NoD",colnames(dat))] # grep all normal samples
# subject='Diabetes'
# targets<-readTargets(paste(PhenotypeDir,"hg_sg_",subject,"_target.txt", sep=''))
# blocks<-factor(targets$Subject)
# Treat <- factor(targets$Treatment,levels=c("C","T"))
# Replicates <- factor(targets$rep)
# design <- model.matrix(~Replicates+Treat)
# 
# corfit <- duplicateCorrelation(Diabetes, block = blocks)
# fit <- lmFit(Diabetes, design,block = blocks, cor = corfit$consensus.correlation)
# fit2=eBayes(fit)
# 
# qval.cutoff=0.2; FC.cutoff=0.17 # FC=1.2
# x5=topTable(fit2,coef="TreatT",n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
# y5=topTable(fit2, coef="TreatT",n=nrow(genes),adjust.method="BH",genelist=genes)
# x5_subset<-x5[x5$adj.P.Val<qval.cutoff & abs(x5$logFC) > FC.cutoff,]
# 
# pdf(file = paste(Figs,subject,"_q",qval.cutoff,"_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# with(y5, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Diabetes with and without complication subjects: high glucose vs. standard glucose (q.val <",qval.cutoff,"& FC >",FC.cutoff,sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-0.5,0.5)))
# with(subset(y5, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
# with(subset(y5, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
# dev.off()
# 
# cat("number of genes (q-val <0.2): ",sum(x5$adj.P.Val<qval.cutoff))
# write.csv(x5_subset, paste(Output,subject,"_",qval.cutoff,"_FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
# write.csv(x5, paste(Output,subject,"_",qval.cutoff,"_",date,".csv",sep=""),row.names=F)
# write.csv(y5, paste(Output,subject,"_",date,".csv",sep=""),row.names=F)
# 
# pdf(file = paste(Figs,subject,"_pvalue_distribution_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# par(mfrow=c(1,2))
# hist(y5$P.Value, main="",xlab="p-value")
# hist(y5$adj.P.Val,main="", xlab="adjusted p-value")
# mtext(paste(subject ," subject (treatment effect)",sep=""), side=3, outer=TRUE, line=-3)
# dev.off()
# 
# # #--------------------------------------------------------------------------
# # ## DE analysis on Diabetes (DwC+DwoC) after adjusting for growth rate
# growth_rate<-meta$GROWTH_RATE; gr <- factor(rep(growth_rate,each=6))
# design2 <- model.matrix(~Treat+gr)
# 
# corfit <- duplicateCorrelation(Diabetes, block = blocks)
# fit3 <- lmFit(Diabetes, design2,block = blocks, cor = corfit$consensus.correlation)
# fit4=eBayes(fit3)
# 
# qval.cutoff=0.2; FC.cutoff=0.17 # FC=1.2
# x5=topTable(fit4,coef="TreatT",n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
# y5=topTable(fit4, coef="TreatT",n=nrow(genes),adjust.method="BH",genelist=genes)
# x5_subset<-x5[x5$adj.P.Val<qval.cutoff & abs(x5$logFC) > FC.cutoff,]
# 
# pdf(file = paste(Figs,subject,"no_growth_rate_q",qval.cutoff,"_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# with(y5, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Diabetes with and without complication subjects: high glucose vs. standard glucose (q.val <",qval.cutoff,"& FC >",FC.cutoff,sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-0.5,0.5)))
# with(subset(y5, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
# with(subset(y5, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
# dev.off()
# 
# cat("number of genes (q-val <0.2): ",sum(x5$adj.P.Val<qval.cutoff))
# write.csv(x5_subset, paste(Output,subject,"_no_growth_rate_",qval.cutoff,"_FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
# write.csv(x5, paste(Output,subject,"_no_growth_rate_",qval.cutoff,"_",date,".csv",sep=""),row.names=F)
# write.csv(y5, paste(Output,subject,"_no_growth_rate_",date,".csv",sep=""),row.names=F)
# 
# pdf(file = paste(Figs,subject,"_no_growth_rate_pvalue_distribution_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# par(mfrow=c(1,2))
# hist(y5$P.Value, main="",xlab="p-value")
# hist(y5$adj.P.Val,main="", xlab="adjusted p-value")
# mtext(paste(subject ," subject (treatment effect) without growth rate",sep=""), side=3, outer=TRUE, line=-3)
# dev.off()
# #
# 
# # ## DE analysis on Diabetes (DwC+DwoC) after adjusting for DCCT_SBP_MEAN
# rm(fit3,fit4)
# dcct_sbp_mean<-meta$DCCT_SBP_MEAN; dsm <- factor(rep(dcct_sbp_mean,each=6))
# design2 <- model.matrix(~Treat+dsm)
# 
# corfit <- duplicateCorrelation(Diabetes, block = blocks)
# fit3 <- lmFit(Diabetes, design2,block = blocks, cor = corfit$consensus.correlation)
# fit4=eBayes(fit3)
# 
# qval.cutoff=0.2; FC.cutoff=0.17 # FC=1.2
# x5=topTable(fit4,coef="TreatT",n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
# y5=topTable(fit4, coef="TreatT",n=nrow(genes),adjust.method="BH",genelist=genes)
# x5_subset<-x5[x5$adj.P.Val<qval.cutoff & abs(x5$logFC) > FC.cutoff,]
# 
# pdf(file = paste(Figs,subject,"no_dcct_sbp_mean_q",qval.cutoff,"_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# with(y5, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Diabetes with and without complication subjects: high glucose vs. standard glucose (q.val <",qval.cutoff,"& FC >",FC.cutoff,sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-0.5,0.5)))
# with(subset(y5, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
# with(subset(y5, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
# dev.off()
# 
# cat("number of genes (q-val <0.2): ",sum(x5$adj.P.Val<qval.cutoff))
# write.csv(x5_subset, paste(Output,subject,"_no_dcct_sbp_mean_",qval.cutoff,"_FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
# write.csv(x5, paste(Output,subject,"_no_dcct_sbp_mean_",qval.cutoff,"_",date,".csv",sep=""),row.names=F)
# write.csv(y5, paste(Output,subject,"_no_dcct_sbp_mean_",date,".csv",sep=""),row.names=F)
# 
# pdf(file = paste(Figs,subject,"_no_dcct_sbp_mean_pvalue_distribution_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# par(mfrow=c(1,2))
# hist(y5$P.Value, main="",xlab="p-value")
# hist(y5$adj.P.Val,main="", xlab="adjusted p-value")
# mtext(paste(subject ," subject (treatment effect) without dcct_sbp_mean",sep=""), side=3, outer=TRUE, line=-3)
# dev.off()


# # All_DE_input2<-matrix(1,nrow(All_DE2),ncol(All_DE2))
# # growth_rate2<-rep(growth_rate,each=3)
# # griwth_rate2<-rep(meta$DCCT_SBP_MEAN,each=3)
#
# # for (i in 1:nrow(All_DE_input2)){
# # gene<-as.numeric(All_DE2[i,])
# # All_DE_input2[i,]<-summary(lm(gene~growth_rate2))$residual
# # }
# # rownames(All_DE_input2)<-rownames(All_DE2)
# # colnames(All_DE_input2)<-colnames(All_DE2)
#
# # design=model.matrix(~growth_rate2)
#
# corfit <- duplicateCorrelation(All_DE2, block = blocks)
# fitTrtMean <- lmFit(All_DE2, model.matrix(~as.factor(growth_rate2)), block = blocks, cor = corfit$consensus.correlation)
# efit.contrast=eBayes(fitTrtMean)
#
# qval.cutoff=0.2; FC.cutoff=0.17 # FC=1.2
# x5=topTable(efit.contrast,n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
# y5=topTable(efit.contrast, n=nrow(genes),adjust.method="BH",genelist=genes)
# x5_subset<-x5[x5$adj.P.Val<qval.cutoff & abs(x5$logFC) > FC.cutoff,]
#
# pdf(file = paste(Figs,"no_growth_rate_sg_",subject,"_q",qval.cutoff,"_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# with(y5, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Diabetes with and without complication subjects: high glucose vs. standard glucose (q.val <",qval.cutoff,"& FC >",FC.cutoff,sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-0.5,0.5)))
# with(subset(y5, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
# with(subset(y5, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
# dev.off()
#
# cat("number of genes (q-val <0.2): ",sum(x5$adj.P.Val<qval.cutoff))
# write.csv(x5_subset, paste(Output,"no_growth_rate_sg_",subject,"_",qval.cutoff,"_FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
# write.csv(x5, paste(Output,"no_growth_rate_sg_",subject,"_",qval.cutoff,"_",date,".csv",sep=""),row.names=F)
# write.csv(y5, paste(Output,"no_growth_rate_sg_",subject,"_",date,".csv",sep=""),row.names=F)
#
# pdf(file = paste(Figs,"no_growth_rate_sg_",subject,"_pvalue_distribution_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# par(mfrow=c(1,2))
# hist(y5$P.Value, main="",xlab="p-value")
# hist(y5$adj.P.Val,main="", xlab="adjusted p-value")
# mtext(paste("no growth rate ", subject ," subject (standard glucose)",sep=""), side=3, outer=TRUE, line=-3)
# dev.off()
#----------------------------------------------
#--------------------------------------------------------------------------
## DE analysis on Diabetes (DwC+DwoC) after adjusting for growth rate
# meta<-read.table(paste(PhenotypeDir, "meanDeltaSampleInfo3.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
# Diabetes<-dat[,-grep("NoD",colnames(dat))] # grep all normal samples
# subject="Diabetes"
# targets<-readTargets(paste(PhenotypeDir,"hg_sg_",subject,"_target.txt", sep=''))
# Treat <- factor(targets$Treatment,levels=c("C","T"))
# Replicates <- factor(targets$rep)
# design <- model.matrix(~Replicates+Treat)
#
# corfit <- duplicateCorrelation(Diabetes, block = targets$Subject)
# fit <-lmFit(Diabetes,design,block=targets$Subject,correlation=corfit$consensus.correlation)
# fit2<-eBayes(fit)
# tmp2<-topTable(fit2,  coef="TreatT", n=nrow(genes), adjust.method="BH",genelist=genes)
# png("test2.png")
# hist(tmp2$P.Value)
# dev.off()


# adjust covariate - growth rate

# growth_rate<-meta$GROWTH_RATE; gr <- factor(rep(growth_rate,each=6))
# design3 <- model.matrix(~Treat+gr)
# corfit3 <- duplicateCorrelation(Diabetes, block = targets$Subject)
# fit4 <-lmFit(Diabetes,design3,block=targets$Subject,correlation=corfit3$consensus.correlation)
# fit4<-eBayes(fit4)
# tmp4<-topTable(fit4,coef="TreatT", n=nrow(genes), adjust.method="BH",genelist=genes)
# png("test4.png")
# hist(tmp4$P.Value)
# dev.off()
#
#
#
#
# fit2<-contrasts.fit(fit,cm)
# fit2 <- eBayes(fit2)
#
# topTable(fit2, coef="DiseasedvsNormalForTissueA")
# fitTrtMean <- lmFit(All_DE2, model.matrix(~as.factor(growth_rate2)), block = blocks, cor = corfit$consensus.correlation)
# efit.contrast=eBayes(fitTrtMean)
#
# qval.cutoff=0.2; FC.cutoff=0.17 # FC=1.2
# x5=topTable(efit.contrast,n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
# y5=topTable(efit.contrast, n=nrow(genes),adjust.method="BH",genelist=genes)
# x5_subset<-x5[x5$adj.P.Val<qval.cutoff & abs(x5$logFC) > FC.cutoff,]
#
# pdf(file = paste(Figs,"no_growth_rate_sg_",subject,"_q",qval.cutoff,"_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# with(y5, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Diabetes with and without complication subjects: high glucose vs. standard glucose (q.val <",qval.cutoff,"& FC >",FC.cutoff,sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-0.5,0.5)))
# with(subset(y5, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
# with(subset(y5, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
# dev.off()
#
# cat("number of genes (q-val <0.2): ",sum(x5$adj.P.Val<qval.cutoff))
# write.csv(x5_subset, paste(Output,"no_growth_rate_sg_",subject,"_",qval.cutoff,"_FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
# write.csv(x5, paste(Output,"no_growth_rate_sg_",subject,"_",qval.cutoff,"_",date,".csv",sep=""),row.names=F)
# write.csv(y5, paste(Output,"no_growth_rate_sg_",subject,"_",date,".csv",sep=""),row.names=F)
#
# pdf(file = paste(Figs,"no_growth_rate_sg_",subject,"_pvalue_distribution_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# par(mfrow=c(1,2))
# hist(y5$P.Value, main="",xlab="p-value")
# hist(y5$adj.P.Val,main="", xlab="adjusted p-value")
# mtext(paste("no growth rate ", subject ," subject (standard glucose)",sep=""), side=3, outer=TRUE, line=-3)
# dev.off()
# #----------------------------------------------
