
#------------------------------------------------------------
#  Directory Setup
#------------------------------------------------------------
HOME <- "/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/"  # change if needed
RData=paste(HOME,'rdata/',sep=''); dir.create(RData, showWarnings = FALSE)
Figs = paste(HOME,'Figs/', sep=''); dir.create(Figs, showWarnings = FALSE)
Output = paste(HOME,'output/', sep=''); dir.create(Figs, showWarnings = FALSE)
PhenotypeDir = paste(HOME,'data/phenotype/', sep='')

# source(paste(HOME,'script/Function.R', sep=''))

#------------------------------------------------------------
# batch effect check
#------------------------------------------------------------
library(ExpressionNormalizationWorkflow);

#--- covarites for batch, patient group
# expression matrix: 21964 X 144
exprs<-read.table(paste(Output, "normalizedDataMatrix.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
covrts<-read.table(paste(PhenotypeDir, "sampleInfo.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
inpData<-expSetobj(exprs,covrts)
pct_thrsh <- 0.75
cvrts_eff_var<-c("batch","patient_group","treatment")
png(paste(Figs,"batch_effect.jpg",sep=""),width = 1024, height = 768)
pvcAnaly(inpData, pct_thrsh, cvrts_eff_var)
dev.off()

#--- covarites from GWU; this does not include normal samples
# expression matrix: 21964 X 102 (excluding normal samples)
covrts2<-read.csv(paste(PhenotypeDir, "sampleInfo2.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
exprs2<-exprs[,rownames(covrts2)]
inpData2<-expSetobj(exprs2,covrts2)
cvrts_eff_var2<-c("batch","group","treatment", "HAB1C","sex", "HDL3", "LDL3", "SBP3","DBP3","PULSE3", "BMI3", "AGE3")
png(paste(Figs,"PCA_covarites_all_samples.jpg",sep=""),width = 1024, height = 768)
pvcAnaly(inpData2, pct_thrsh, cvrts_eff_var2)
dev.off()

#--- with avergaed expression data
# expression matrix: 21964 X 30 (excluding normal samples)
cvrts3<-read.table(paste(PhenotypeDir, "meanSampleInfo.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
meanExprs<-read.table(paste(Output, "meanExprs.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
inpData3<-expSetobj(meanExprs,cvrts3)
cvrts_eff_var3<-c("group","treatment", "sex", "HDL3", "LDL3", "SBP3","DBP3","PULSE3", "BMI3", "AGE3")
png(paste(Figs,"PCA_covarites_all_averaged_samples.jpg",sep=""),width = 1024, height = 768,  pointsize = 14)
pvcAnaly(inpData3, pct_thrsh, cvrts_eff_var3)
dev.off()

#--- with delta data
cvrts4<-read.table(paste(PhenotypeDir, "meanDeltaSampleInfo.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
meanDeltaExprs<-read.table(paste(Output, "meanDeltaExprs.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
inpData4<-expSetobj(meanDeltaExprs,cvrts4)
cvrts_eff_var4<-c("group", "sex", "HDL3", "LDL3", "SBP3","DBP3","PULSE3", "BMI3", "AGE3")
png(paste(Figs,"PCA_covarites_averaged_delta_samples.jpg",sep=""),width = 1024, height = 768,  pointsize = 14)
pvcAnaly(inpData4, pct_thrsh, cvrts_eff_var4)
dev.off()
