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
  
  rm(lumi.Ds1, lumi.Ds2, lumi.Ds3, lumi.Ds4, rawData.ctrl1, rawData.ctrl2, rawData.ctrl3, rawData.ctrl4)
  
  #------------------------------------------------------------
  #  Lumi: 3. Data Processing 
  #  a. background correction
  #  b. variance stabilization (1.log2, 2.vst)
  #  c. quantile normalization
  #------------------------------------------------------------
  
  LOG1=TRUE
  
  # b1. log2 & quantile
  if(LOG1){
    Batch1 <-  lumiExpresso(rawData.ctrl1.nuId, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), variance.stabilize = TRUE, 
                            varianceStabilize.param = list(method='log2'), normalize = FALSE)
    Batch2 <-  lumiExpresso(rawData.ctrl2.nuId, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), variance.stabilize = TRUE, 
                            varianceStabilize.param = list(method='log2'), normalize = FALSE)
    Batch3 <-  lumiExpresso(rawData.ctrl3.nuId, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), variance.stabilize = TRUE, 
                            varianceStabilize.param = list(method='log2'), normalize = FALSE)
    Batch4 <-  lumiExpresso(rawData.ctrl4.nuId, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), variance.stabilize = TRUE, 
                            varianceStabilize.param = list(method='log2'), normalize = FALSE)
  }
  
  
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
  
  normalizedBatch <- lumiN(Batch, method='rsn')
  # normalizedBatch <- lumiN(Batch, method='quantile')
  
  rm(rawData.ctrl1.nuId, rawData.ctrl2.nuId, rawData.ctrl3.nuId, rawData.ctrl4.nuId,Batch1, 
     Batch2, Batch3, Batch4, exprs.Batch1,exprs.Batch2, 
     exprs.Batch3, exprs.Batch4)
  
  sampleNames(rawBatch)<-gsub("SC.", "", sampleNames(rawBatch)) # replace sample name
  sampleNames(normalizedBatch)<-gsub("SC.", "", sampleNames(normalizedBatch)) # replace sample name
  rm(common_probe_names)
  #------------------------------------------------------------
  # assign informative sample name and sort by sample name
  #------------------------------------------------------------
  sampleInfo<-read.table(paste(PhenotypeDir,"20160203.MAG.layout for H12 chips_updated.txt",sep=''))
  sampleInfo<-sampleInfo[mixedorder(sampleInfo$V1),]
  sampleNames(rawBatch)<-sampleInfo$V2
  sampleNames(normalizedBatch)<-sampleInfo$V2
  
  rawBatch.sorted<-rawBatch[,mixedorder(sampleNames(rawBatch), decreasing=FALSE)]
  save(rawBatch.sorted, file = paste(RData,"rawBatch.sorted.RData",sep=""))
  normalizedBatch.sorted<-normalizedBatch[,mixedorder(sampleNames(normalizedBatch), decreasing=FALSE)]
  save(normalizedBatch.sorted, file = paste(RData,"normalizedBatch.sorted.RData",sep=""))
  rm(rawBatch, normalizedBatch)
  
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
} 
#------------------------------------------------------------
# data filtering
#------------------------------------------------------------


filtering=TRUE

if(filtering) {
  filter.Th = 0.01; filter.dp = 0
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

normalizedDataMatrix.sorted<-normalizedDataMatrix[,mixedorder(colnames(normalizedDataMatrix), decreasing=FALSE)]
save(normalizedDataMatrix.sorted, file = paste(RData,"normalizedDataMatrix.sorted.RData",sep=""))
rm(normalizedDataMatrix, eset.filtered.normalizedData)


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
