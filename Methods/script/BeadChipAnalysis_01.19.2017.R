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
Figs = paste(HOME,'Figs/', sep=''); dir.create(Figs, showWarnings = FALSE)
Output = paste(HOME,'output/', sep=''); dir.create(Figs, showWarnings = FALSE)

source(paste(HOME,'script/Function.R', sep=''))

REUSE = FALSE

if (REUSE == TRUE) {
  load(file = paste(RData,"rawBatch.sorted.RData",sep=""))
  load(file = paste(RData,"normalizedBatch.sorted.RData",sep=""))
} else {
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
  #  Lumi: 1. Read Expression Data from GenomeStudio
  #------------------------------------------------------------
  lumi.Ds1 <- lumiR(InputFileList[1]); 
  lumi.Ds2 <- lumiR(InputFileList[2]); 
  lumi.Ds3 <- lumiR(InputFileList[3]); 
  lumi.Ds4 <- lumiR(InputFileList[4]); 
  
  #------------------------------------------------------------
  #  Lumi: 2. Add Control Data 
  #------------------------------------------------------------
  rawData.ctrl1 <- addControlData2lumi(ControlFileList[1], lumi.Ds1)
  rawData.ctrl2 <- addControlData2lumi(ControlFileList[2], lumi.Ds2)
  rawData.ctrl3 <- addControlData2lumi(ControlFileList[3], lumi.Ds3)
  rawData.ctrl4 <- addControlData2lumi(ControlFileList[4], lumi.Ds4)
  
  if (require(lumiHumanIDMapping)){
    rawData.ctrl1.nuId <- addNuID2lumi(rawData.ctrl1, lib.mapping='lumiHumanIDMapping')
    rawData.ctrl2.nuId <- addNuID2lumi(rawData.ctrl2, lib.mapping='lumiHumanIDMapping')
    rawData.ctrl3.nuId <- addNuID2lumi(rawData.ctrl3, lib.mapping='lumiHumanIDMapping')
    rawData.ctrl4.nuId <- addNuID2lumi(rawData.ctrl4, lib.mapping='lumiHumanIDMapping')
  }
  
  rm(lumi.Ds1, lumi.Ds2, lumi.Ds3, lumi.Ds4, rawData.ctrl1, rawData.ctrl2, rawData.ctrl3, rawData.ctrl4)
  
  #------------------------------------------------------------
  #  Lumi: 3. Combine Data and Perform preprocessing 
  #  a. background correction; b. variance stabilization; 
  #  c. normalization
  #------------------------------------------------------------
  exprs1 <- exprs(rawData.ctrl1.nuId); exprs2 <- exprs(rawData.ctrl2.nuId);
  exprs3 <- exprs(rawData.ctrl3.nuId); exprs4 <- exprs(rawData.ctrl4.nuId);
  
  dataList <- list(exprs1, exprs2, exprs3, exprs4)
  common_probe_names = Reduce(intersect, lapply(dataList, row.names))
  
  rawBatch<-combine(rawData.ctrl1.nuId[common_probe_names,],rawData.ctrl2.nuId[common_probe_names,],
                    rawData.ctrl3.nuId[common_probe_names,],rawData.ctrl4.nuId[common_probe_names,])
  
  normalizedBatch <-  lumiExpresso(rawBatch, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), variance.stabilize = TRUE, list(method='log2'), 
                                   normalize = TRUE, list(method='rsn'))
  
  rm(exprs1, exprs2, exprs3, exprs4,rawData.ctrl1.nuId, rawData.ctrl2.nuId, rawData.ctrl3.nuId, rawData.ctrl4.nuId)
  rm(common_probe_names, dataList)
  
  #------------------------------------------------------------
  # assign informative sample name and sort by sample name
  #------------------------------------------------------------
  rawBatch.sorted<-rawBatch[,mixedorder(sampleNames(rawBatch), decreasing=TRUE)]
  # save(rawBatch.sorted, file = paste(RData,"rawBatch.sorted.RData",sep=""))
  normalizedBatch.sorted<-normalizedBatch[,mixedorder(sampleNames(normalizedBatch), decreasing=TRUE)]
  # save(normalizedBatch.sorted, file = paste(RData,"normalizedBatch.sorted.RData",sep=""))
  rm(rawBatch, normalizedBatch)
}

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
# rm(colList)

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

#------------------------------------------------------------
# data filtering
#------------------------------------------------------------
filtering=TRUE

if(filtering) {
  filter.Th = 0.01; filter.dp = round(ncol(normalizedBatch.sorted)*.1) # detection 20%
  # will save a filtered.normData file in the working dir.
  eset.filtered.normalizedData = filterFun(rawBatch.sorted, normalizedBatch.sorted, filter.Th, filter.dp)
}

# remove probes that do not have annotated genes --------------------------
if (require(lumiHumanAll.db) & require(annotate)) {
  normalizedDataMatrix <- eset.filtered.normalizedData[!is.na(getSYMBOL(rownames(eset.filtered.normalizedData), 'lumiHumanAll.db')),]
}

probeList <- rownames(normalizedDataMatrix)
if (require(lumiHumanAll.db) & require(annotate)){
  geneSymbol<- getSYMBOL(probeList, 'lumiHumanAll.db')
  geneName<- sapply(lookUp(probeList, 'lumiHumanAll.db', 'GENENAME'), function(x) x[1])
}

sampleInfo<-read.table(paste(PhenotypeDir,"20160203.MAG.layout for H12 chips_updated.txt",sep=''))
sampleInfo<-sampleInfo[mixedorder(sampleInfo$V1),]
colnames(normalizedDataMatrix)<-sampleInfo$V2 

normalizedDataMatrix.sorted<-normalizedDataMatrix[,mixedorder(colnames(normalizedDataMatrix), decreasing=FALSE)]
rm(normalizedDataMatrix, eset.filtered.normalizedData)

## Average technical replicates
S1<-c("DwC_1026_30mM.2", "DwC_1026_30mM.3","DwC_1026_30mM.4")
SN1<-c("DwC_1026_norm.2", "DwC_1026_norm.3","DwC_1026_norm.4")
S2<-c("DwoC_2318_30mM.3", "DwoC_2318_30mM.4","DwoC_2318_30mM.5")
SN2<-c("DwoC_2318_norm.3", "DwoC_2318_norm.4","DwoC_2318_norm.5")
S3<-c("DwoC_25224_30mM.1", "DwoC_25224_30mM.2","DwoC_25224_30mM.3")
SN3<-c("DwoC_25224_norm.1", "DwoC_25224_norm.2","DwoC_25224_norm.3")

S1.expr<-data.frame(apply(normalizedDataMatrix.sorted[,S1],1,mean));  colnames(S1.expr)<-S1[1]
SN1.expr<-data.frame(apply(normalizedDataMatrix.sorted[,SN1],1,mean));  colnames(SN1.expr)<-SN1[1]

S2.expr<-data.frame(apply(normalizedDataMatrix.sorted[,S2],1,mean));  colnames(S2.expr)<-S2[1]
SN2.expr<-data.frame(apply(normalizedDataMatrix.sorted[,SN2],1,mean));  colnames(SN2.expr)<-SN2[1]

S3.expr<-data.frame(apply(normalizedDataMatrix.sorted[,S3],1,mean));  colnames(S3.expr)<-S3[1]
SN3.expr<-data.frame(apply(normalizedDataMatrix.sorted[,SN3],1,mean));  colnames(SN3.expr)<-SN3[1]

drops=c(S1, S2, S3, SN1, SN2, SN3)
normalizedDataMatrix.sorted<-normalizedDataMatrix.sorted[ , !(colnames(normalizedDataMatrix.sorted) %in% drops)]

normalizedDataMatrix.sorted2<-cbind(normalizedDataMatrix.sorted, S1.expr, SN1.expr, S2.expr, SN2.expr, S3.expr, SN3.expr)
normalizedDataMatrix.sorted2<-normalizedDataMatrix.sorted2[,mixedorder(colnames(normalizedDataMatrix.sorted2), decreasing=FALSE)]
targets2<-readTargets(paste(PhenotypeDir,"target2.txt", sep=''))
colnames(normalizedDataMatrix.sorted2) <- targets2$Case_Number

save(normalizedDataMatrix.sorted2, file = paste(RData,"normalizedDataMatrix.sorted.RData",sep=""))
#--------------------------------------------------------------------------------------
# limma: differential expression analysis -  High Glucose (HG) vs Standard Glucose (SG)
#--------------------------------------------------------------------------------------

# create a design matrix
targets<-readTargets(paste(PhenotypeDir,"target2.txt", sep=''))

#--- test 1 (treatment)
ct<-factor(targets$Treatment) 
ct<-factor(paste(targets$Group, targets$Treatment,sep=".")) #testing including group
design <- model.matrix(~0+ct)
colnames(design) <- levels(ct)

fit<-lmFit(normalizedDataMatrix.sorted2, design) # linear model fit for each gene
contrast.matrix <- makeContrasts(HG_SG = glucose-control, levels=design) # create a contrast matrix    
fit2<-contrasts.fit(fit, contrast.matrix) # given a linear model fit to the data, compute estimated coefficients and std errors
fit2<-eBayes(fit2); #  Given a linear model fit, compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression 
#  by empirical Bayes moderation
fit2$genes <-data.frame(ID=probeList, geneSymbol=geneSymbol, geneName=geneName, stringsAsFactors=FALSE)

#--- test 2 (treatment + group)  MODIFY THIS PART

ct3<-factor(paste(targets$Group, targets$Treatment,sep=".")) #testing including group
design3 <- model.matrix(~0+ct3)
colnames(design3) <- levels(ct3)

fit3<-lmFit(normalizedDataMatrix.sorted2, design3) # linear model fit for each gene
contrast.matrix3 <- makeContrasts(HG_SG = glucose-control, levels=design) # create a contrast matrix    
fit4<-contrasts.fit(fit3, contrast.matrix3) # given a linear model fit to the data, compute estimated coefficients and std errors
fit4<-eBayes(fit4); #  Given a linear model fit, compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression 
#  by empirical Bayes moderation
fit4$genes <-data.frame(ID=probeList, geneSymbol=geneSymbol, geneName=geneName, stringsAsFactors=FALSE)


#  Volcano plot for High Glucose vs. Standard Glucose

pval.cutoff<-0.01
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
  HighGlucoseSamples <-c(seq(1,3),seq(7,9), seq(13,15), seq(19,21), seq(25,27), seq(31,33), seq(37,39), seq(43,45), 
                         seq(49,51), seq(55,57), seq(61,63), seq(67,69), seq(73,75), seq(79,81), seq(85,87), seq(91,93),
                         seq(97,99), seq(103,105), seq(109,111), seq(115,117), seq(121,123), seq(127,129))
  
  NoGlucoseSamples <-c(seq(4,6),seq(10,12), seq(16,18), seq(22,24), seq(28,30), seq(34,36), seq(40,42), seq(46,48), 
                       seq(52,54), seq(58,60), seq(64,66), seq(70,72), seq(76,78), seq(82,84), seq(88,90), seq(94,96),
                       seq(100,102), seq(106,108), seq(112,114), seq(118,120), seq(124,126), seq(130,132))
  
  HighGlucoseData<-normalizedDataMatrix.sorted2[,HighGlucoseSamples]  # data in log2 scale
  NoGlucoseData<-normalizedDataMatrix.sorted2[,NoGlucoseSamples]      # data in log2 scale
  deltaData <- HighGlucoseData-NoGlucoseData
  colnames(deltaData)<-gsub("_30mM", "", colnames(deltaData))
  save(deltaData, file = paste(RData,"deltaData_nuID.RData",sep=""))
  write.csv(deltaData, file= paste(Output,"deltaData_nuID.csv", sep=""))
  deltaData.T<-data.frame(t(deltaData))
  deltaData.T$Name<-substr(colnames(deltaData),1,nchar(colnames(deltaData))-2)
  meanDeltaData.T<-aggregate(deltaData.T[,-ncol(deltaData.T)],list(deltaData.T$Name), mean)
  
  group1=seq(1,nrow(meanDeltaData.T))
  group2=c(rep(0,8),rep(1,7),rep(2,7))
  rownames(meanDeltaData.T)<-meanDeltaData.T$Group.1
  meanDeltaData.T$Group.1<-NULL
  
  meanDeltaData <- t(meanDeltaData.T)
  rownames(meanDeltaData) <- rownames(deltaData)
  save(meanDeltaData, file = paste(RData,"meanDeltaData_nuID.RData",sep=""))
  write.csv(meanDeltaData, file= paste(Output,"meanDeltaData_nuID.csv", sep=""))
  rm(deltaData.T, meanDeltaData.T, normalizedDataMatrix.sorted, 
     HighGlucoseData, NoGlucoseData, geneName, geneSymbol,
     normalizedBatch.sorted, probeList, rawBatch.sorted)
}

if (REUSE){
  load(file = paste(RData,"deltaData_nuID.RData",sep=""))
  load(file = paste(RData,"meanDeltaData_nuID.RData",sep=""))
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
deltaData.genesymbol<-data.frame(deltaData[rownames(HG_SG),])
deltaData.genesymbol$geneSymbol<-HG_SG$geneSymbol
write.csv(deltaData.genesymbol, file= paste(Output,"deltaData.symbol.csv", sep=""))
deltaData.genesymbol[,-ncol(deltaData.genesymbol)] <- 2^deltaData.genesymbol[,-ncol(deltaData.genesymbol)]
output1<-matrix(,, ncol = ncol(deltaData.genesymbol))
colnames(output1)<-colnames(deltaData.genesymbol)
for (i in 1:length(MG_geneList)) {
  tmp=deltaData.genesymbol[which(MG_geneList[i]==deltaData.genesymbol$geneSymbol), , drop=FALSE]
  print(nrow(tmp))
  output1<-rbind(output1,tmp)
}
output1 <- na.omit(output1)
write.csv(output1, file= paste(Output,"deltaData.5genes.csv", sep=""))

# meandelta data --------------------------------------------------------------
meanDeltaData.genesymbol<-data.frame(meanDeltaData[rownames(HG_SG),])
meanDeltaData.genesymbol$geneSymbol<-HG_SG$geneSymbol
write.csv(meanDeltaData.genesymbol, file= paste(Output,"meanDeltaData.symbol.csv", sep=""))
meanDeltaData.genesymbol[,-ncol(meanDeltaData.genesymbol)] <- 2^meanDeltaData.genesymbol[,-ncol(meanDeltaData.genesymbol)]
output2<-matrix(, ,ncol(meanDeltaData.genesymbol))
colnames(output2)<-colnames(meanDeltaData.genesymbol)
for (i in 1:length(MG_geneList)) {
  tmp=meanDeltaData.genesymbol[which(MG_geneList[i]==meanDeltaData.genesymbol$geneSymbol), , drop=FALSE]
  output2<-rbind(output2,tmp)
}
output2 <- na.omit(output2)
write.csv(output2, file= paste(Output,"meanDeltaData.5genes.csv", sep=""))

rm(deltaData, meanDeltaData)

#--------------------------------------------------------------------------
#  find genes with > 2 in any samples
#--------------------------------------------------------------------------
deltaData.genesymbol$max<-apply(deltaData.genesymbol[,-ncol(deltaData.genesymbol)],1,max)
deltaData.genesymbol$min<-apply(deltaData.genesymbol[,-ncol(deltaData.genesymbol)],1,min)
deltaData.genesymbol$min<-as.numeric(deltaData.genesymbol$min)


meanDeltaData.genesymbol$max<-apply(meanDeltaData.genesymbol[,-ncol(meanDeltaData.genesymbol)],1,max)
meanDeltaData.genesymbol$min<-apply(meanDeltaData.genesymbol[,-ncol(meanDeltaData.genesymbol)],1,min)
meanDeltaData.genesymbol$min<-as.numeric(meanDeltaData.genesymbol$min)

## 2-fold up or down in at least 1 gene across samples 
delta.twofold.up<-deltaData.genesymbol[which(deltaData.genesymbol$max >= 2),]
delta.twofold.down<-deltaData.genesymbol[which(deltaData.genesymbol$min <= 0.5),]

meanDeltaData.twofold.up<-meanDeltaData.genesymbol[which(meanDeltaData.genesymbol$max >= 2),]
meanDeltaData.twofold.down<-meanDeltaData.genesymbol[which(meanDeltaData.genesymbol$min <= 0.5),]

write.csv(delta.twofold.up, file= paste(Output,"delta.twofold.up.csv", sep=""))
write.csv(delta.twofold.down, file= paste(Output,"delta.twofold.down.csv", sep=""))

write.csv(meanDeltaData.twofold.up, file= paste(Output,"meanDeltaData.twofold.up.csv", sep=""))
write.csv(meanDeltaData.twofold.down, file= paste(Output,"meanDeltaData.twofold.down.csv", sep=""))

