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

probeList <- rownames(dat)  # probeID (nuID) gene annotation
if (require(lumiHumanAll.db) & require(annotate)){
  geneSymbol<- getSYMBOL(probeList, 'lumiHumanAll.db')
  geneName<- sapply(lookUp(probeList, 'lumiHumanAll.db', 'GENENAME'), function(x) x[1])
}
genes <-data.frame(ID=probeList, geneSymbol=geneSymbol, geneName=geneName, stringsAsFactors=FALSE)

# PCA by group (HG, SG, delta)
# 1. No_DM
subject<-"No_DM"
NoD<-dat[,grep("NoD",colnames(dat))] # grep all normal samples
NoD_sg<-NoD[,grep("norm",colnames(NoD))]; NoD_hg<-NoD[,grep("30mM",colnames(NoD))];
NoD_delta <- NoD_hg-NoD_sg

# 2. No_PDR
subject<-"No_PDR"
DwoC<-dat[,grep("DwoC",colnames(dat))] # grep all normal samples
DwoC_sg<-DwoC[,grep("norm",colnames(DwoC))]; DwoC_hg<-DwoC[,grep("30mM",colnames(DwoC))];
DwoC_delta <- DwoC_hg-DwoC_sg

# 3. PDR
DwC<-dat[,grep("DwC",colnames(dat))] # grep all normal samples
DwC_sg<-DwC[,grep("norm",colnames(DwC))]; DwC_hg<-DwC[,grep("30mM",colnames(DwC))];
DwC_delta <- DwC_hg-DwC_sg


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
corfit$consensus.correlation
fit <-lmFit(NoD,design,block=targets$Subject,correlation=corfit$consensus.correlation)
fit<-eBayes(fit)

qval.cutoff=0.1; FC.cutoff=0.17 # FC=1.12

x1=topTable(fit, coef="TreatT", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="TreatT", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]
NoD_pval<-y1[y1$P.Value<0.05,]
NoD_qval_10<-x1
NoD_qval_20<-topTable(fit, coef="TreatT", n=nrow(genes), p.value=0.2,adjust.method="BH",genelist=genes)
NoD_qval_50<-topTable(fit, coef="TreatT", n=nrow(genes), p.value=0.5,adjust.method="BH",genelist=genes)

setEPS()
postscript(file = paste(Figs,"NoD_treatment_effect_qval_",qval.cutoff,"_", date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("No_DM treatment effect (high vs. standard glucose)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-.7,.7)))
with(subset(y1, adj.P.Val<qval.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

cat("number of genes (q-val <0.1): ",sum(x1$adj.P.Val<qval.cutoff))
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
DwoC<-dat[,grep("DwoC",colnames(dat))] 
targets<-readTargets(paste(PhenotypeDir,"hg_sg_",subject,"_target.txt", sep=''))
Treat <- factor(targets$Treatment,levels=c("C","T"))
Replicates <- factor(targets$rep)
design <- model.matrix(~Replicates+Treat)

corfit <- duplicateCorrelation(DwoC, block = targets$Subject)
corfit$consensus.correlation
fit <-lmFit(DwoC,design,block=targets$Subject,correlation=corfit$consensus.correlation)
fit<-eBayes(fit)

qval.cutoff=0.1; FC.cutoff=0.17 # FC=1.12

x1=topTable(fit, coef="TreatT", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="TreatT", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]
DwoC_pval<-y1[y1$P.Value<0.05,]
DwoC_qval<-x1

DwoC_qval_10<-x1
DwoC_qval_20<-topTable(fit, coef="TreatT", n=nrow(genes), p.value=0.2,adjust.method="BH",genelist=genes)
DwoC_qval_50<-topTable(fit, coef="TreatT", n=nrow(genes), p.value=0.5,adjust.method="BH",genelist=genes)

setEPS()
postscript(file = paste(Figs,"DwoC_treatment_effect_qval_",qval.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("No_PDR treatment effect (high vs. standard glucose)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-.5,.5)))
with(subset(y1, adj.P.Val<qval.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()
# with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
# with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))

cat("number of genes (q-val <0.1): ",sum(x1$adj.P.Val<qval.cutoff))
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
DwC<-dat[,grep("DwC",colnames(dat))]
targets<-readTargets(paste(PhenotypeDir,"hg_sg_",subject,"_target.txt", sep=''))
Treat <- factor(targets$Treatment,levels=c("C","T"))
Replicates <- factor(targets$rep)
design <- model.matrix(~Replicates+Treat)

corfit <- duplicateCorrelation(DwC, block = targets$Subject)
corfit$consensus.correlation
fit <-lmFit(DwC,design,block=targets$Subject,correlation=corfit$consensus.correlation)
fit<-eBayes(fit)

qval.cutoff=0.1; FC.cutoff=0.17 # FC=1.12

x1=topTable(fit, coef="TreatT", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="TreatT", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]
DwC_pval<-y1[y1$P.Value<0.05,]
DwC_qval_10<-x1
DwC_qval_20<-topTable(fit, coef="TreatT", n=nrow(genes), p.value=0.2,adjust.method="BH",genelist=genes)
DwC_qval_50<-topTable(fit, coef="TreatT", n=nrow(genes), p.value=0.5,adjust.method="BH",genelist=genes)

setEPS()
postscript(file = paste(Figs,"DwC_treatment_effect_q",qval.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("PDR treatment effect (high vs. standard glucose)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-.5,.5)))
with(subset(y1, adj.P.Val<qval.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
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
corfit$consensus.correlation
fit <-lmFit(Diabetes,design,block=targets$Subject,correlation=corfit$consensus.correlation)
fit<-eBayes(fit)

qval.cutoff=0.1; FC.cutoff=0.17 # FC=1.12

x1=topTable(fit, coef="TreatT", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="TreatT", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]
Diabetes_pval<-y1[y1$P.Value<0.05,]
Diabetes_qval_10<-x1
Diabetes_qval_20<-topTable(fit, coef="TreatT", n=nrow(genes), p.value=0.2,adjust.method="BH",genelist=genes)
Diabetes_qval_50<-topTable(fit, coef="TreatT", n=nrow(genes), p.value=0.5,adjust.method="BH",genelist=genes)

setEPS()
postscript(file = paste(Figs,"Diabetes_treatment_effect_qval_",qval.cutoff,"_logFC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("No_PDR and PDR treatment effect (high vs. standard glucose)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-0.5,0.5)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

cat("number of genes (q-val <0.1): ",sum(x1$adj.P.Val<qval.cutoff))
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
corfit$consensus.correlation
fit <-lmFit(dat,design,block=targets$Subject,correlation=corfit$consensus.correlation)
fit<-eBayes(fit)

qval.cutoff=0.1; FC.cutoff=0.17 # FC=1.12

x1=topTable(fit, coef="TreatT", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="TreatT", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]
All_pval<-y1[y1$P.Value<0.05,]
All_qval_10<-x1
All_qval_20<-topTable(fit, coef="TreatT", n=nrow(genes), p.value=0.2,adjust.method="BH",genelist=genes)
All_qval_50<-topTable(fit, coef="TreatT", n=nrow(genes), p.value=0.5,adjust.method="BH",genelist=genes)

setEPS()
postscript(file = paste(Figs,"All_treatment_effect_qval_",qval.cutoff,"_logFC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("No_DM, No_PDR and PDR treatment effect (high vs. standard glucose)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-0.5,0.5)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

cat("number of genes (q-val <0.1): ",sum(x1$adj.P.Val<qval.cutoff))
write.csv(x1_subset, paste(Output,"Diabetes_treatment_effect_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1, paste(Output,"All_treatment_effect_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(y1, paste(Output,"All_treatment_effect_",date,".csv",sep=""),row.names = F)

setEPS()
postscript(file = paste(Figs,"All_treatment_effect_pvalue_distribution_",date,".eps",sep=""),width=18,height=12,pointsize=20)
par(mfrow=c(1,2))
hist(y1$P.Value, main="",xlab="p-value")
hist(y1$adj.P.Val,main="", xlab="adjusted p-value")
mtext("No_DM, No_PDR and PDR treatment effect (high vs. standard glucose)", side=3, outer=TRUE, line=-3)
dev.off()

##----------------------------------
## group comparison (each hg, sg, and delta)
##----------------------------------

#1a. delta
delta<-cbind(DwC_delta,DwoC_delta,NoD_delta)
colnames(delta)<-substr(colnames(delta), 1, nchar(colnames(delta))-7)
subject="all_replicate"
targets<-readTargets(paste(PhenotypeDir,"delta_",subject,"_target.txt", sep=''))
group <- factor(targets$Group)
design <- model.matrix(~0+group)
colnames(design)<-levels(group)
rownames(design)<-targets$SampleName
corfit <- duplicateCorrelation(delta, block = targets$Subject)
fit <- lmFit(delta, design, block = targets$Subject, cor = corfit$consensus.correlation)

cm <- makeContrasts( delta_NoD_DwoC=NoD - DwoC, delta_NoD_DwC=NoD - DwC, 
                     delta_NoD_Diabetes=NoD - (DwC+DwoC)/2, delta_DwC_DwoC = DwC-DwoC,
                     delta_DwoC_DwC = DwoC-DwC,
                     levels=design) # create a contrast matrix
fit2=contrasts.fit(fit,cm)
fit2=eBayes(fit2)

setEPS()
postscript(file = paste(Figs,"delta_all_group_pvalues_",date,".eps",sep=""))

par(mfrow=c(2,3))
for (i in 1:ncol(fit2$p.value)) {
  hist(fit2$p.value[,i], main=colnames(fit2$p.value)[i])
}
dev.off()

#----------------
# delta NoD vs DwoC
#----------------

pval.cutoff=0.05; FC.cutoff=0.263 # FC=1.2

x1=topTable(fit2, coef="delta_NoD_DwoC", n=nrow(genes), adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="delta_NoD_DwoC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$P.Value<pval.cutoff & abs(x1$logFC) > FC.cutoff,]
delta_NoD_DwoC_pval_FC<-x1_subset

setEPS()
postscript(file = paste(Figs,"delta_NoD_DwoC_pvalue_",pval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(P.Value), pch=20, col="gray", main=paste("Group difference in delta (No_DM vs. No_PDR)",sep=""), xlab="log2 fold change", ylab="-log10(P.Value)", xlim=c(-1.0,1.0)))
with(subset(y1, P.Value < pval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(y1, P.Value < pval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"delta_NoD_DwoC_pvalue_all_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"delta_NoD_DwoC_pvalue_",pval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#----------------
# delta NoD vs DwC
#----------------
x1=topTable(fit2, coef="delta_NoD_DwC", n=nrow(genes), adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="delta_NoD_DwC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$P.Value<pval.cutoff & abs(x1$logFC) > FC.cutoff,]
delta_NoD_DwC_pval_FC<-x1_subset

setEPS()
postscript(file = paste(Figs,"delta_NoD_DwC_pvalue_",pval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(P.Value), pch=20, col="gray", main=paste("Group difference in delta (No_DM vs. PDR)",sep=""), xlab="log2 fold change", ylab="-log10(P.Value)", xlim=c(-1.0,1.0)))
with(subset(y1, P.Value<pval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(y1, P.Value<pval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"delta_NoD_DwC_pvalue_all_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"delta_NoD_DwC_pvalue_",pval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#---------------------
# delta NoD vs DwC+DwoC
#---------------------
x1=topTable(fit2, coef="delta_NoD_Diabetes", n=nrow(genes), adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="delta_NoD_Diabetes", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$P.Value<pval.cutoff & abs(x1$logFC) > FC.cutoff,]
delta_NoD_Diabetes_pval_FC<-x1_subset

setEPS()
postscript(file = paste(Figs,"delta_NoD_Diabetes_pvalue_",pval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(P.Value), pch=20, col="gray", main=paste("Group difference in delta (No_DM vs. PDR + No_PDR)",sep=""), xlab="log2 fold change", ylab="-log10(P.Value)", xlim=c(-1.0,1.0)))
with(subset(y1, P.Value<pval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(y1, P.Value<pval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"delta_NoD_Diabetes_pvalue_all_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"delta_NoD_Diabetes_pvalue_",pval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#---------------------
# delta DwC vs DwoC
#---------------------
x1=topTable(fit2, coef="delta_DwC_DwoC", n=nrow(genes), adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="delta_DwC_DwoC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$P.Value<pval.cutoff & abs(x1$logFC) > FC.cutoff,]
delta_DwC_DwoC_pval_FC<-x1_subset

setEPS()
postscript(file = paste(Figs,"delta_DwC_DwoC_pvalue_",pval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(P.Value), pch=20, col="gray", main=paste("Group difference in delta (PDR vs. No_PDR)",sep=""), xlab="log2 fold change", ylab="-log10(P.Value)", xlim=c(-1,1)))
with(subset(y1, P.Value<pval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(y1, P.Value<pval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"delta_DwC_DwoC_pvalue_all_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"delta_DwC_DwoC_pvalue_",pval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#---------------------
# delta DwoC vs DwC
#---------------------
x1=topTable(fit2, coef="delta_DwoC_DwC", n=nrow(genes), adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="delta_DwoC_DwC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$P.Value<pval.cutoff & abs(x1$logFC) > FC.cutoff,]
delta_DwoC_DwC_pval_FC<-x1_subset

setEPS()
postscript(file = paste(Figs,"delta_DwoC_DwC_pvalue_",pval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(P.Value), pch=20, col="gray", main=paste("Group difference in delta (No_PDR vs. PDR)",sep=""), xlab="log2 fold change", ylab="-log10(P.Value)", xlim=c(-1,1)))
with(subset(y1, P.Value<pval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(y1, P.Value<pval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"delta_DwoC_DwC_pvalue_all_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"delta_DwoC_DwC_pvalue_",pval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)



library(VennDiagram)

# p-val <0.05
a1<-length(NoD_pval$ID)
a2<-length(DwoC_pval$ID)
a3<-length(DwC_pval$ID)
a4<-length(Diabetes_pval$ID)
a5<-length(All_pval$ID)

n12=length(intersect(NoD_pval$ID,DwoC_pval$ID))
n13=length(intersect(NoD_pval$ID,DwC_pval$ID))
n14=length(intersect(NoD_pval$ID,Diabetes_pval$ID))
n15=length(intersect(NoD_pval$ID,All_pval$ID))
n23=length(intersect(DwoC_pval$ID,DwC_pval$ID))
n24=length(intersect(DwoC_pval$ID,Diabetes_pval$ID))
n25=length(intersect(DwoC_pval$ID,All_pval$ID))
n34=length(intersect(DwC_pval$ID,Diabetes_pval$ID))
n35=length(intersect(DwC_pval$ID,All_pval$ID))
n45=length(intersect(Diabetes_pval$ID,All_pval$ID))

n123=length(intersect(intersect(NoD_pval$ID,DwoC_pval$ID),DwC_pval$ID))
n124=length(intersect(intersect(NoD_pval$ID,DwoC_pval$ID),Diabetes_pval$ID))
n125=length(intersect(intersect(NoD_pval$ID,DwoC_pval$ID),All_pval$ID))
n134=length(intersect(intersect(NoD_pval$ID,DwC_pval$ID),Diabetes_pval$ID))
n135=length(intersect(intersect(NoD_pval$ID,DwC_pval$ID),All_pval$ID))
n145=length(intersect(intersect(NoD_pval$ID,Diabetes_pval$ID),All_pval$ID))
n234=length(intersect(intersect(DwoC_pval$ID,DwC_pval$ID),Diabetes_pval$ID))
n235=length(intersect(intersect(DwoC_pval$ID,DwC_pval$ID),All_pval$ID))
n245=length(intersect(intersect(DwoC_pval$ID,Diabetes_pval$ID),All_pval$ID))
n345=length(intersect(intersect(DwC_pval$ID,Diabetes_pval$ID),All_pval$ID))
n1234=length(intersect(intersect(intersect(NoD_pval$ID,DwoC_pval$ID),DwC_pval$ID),Diabetes_pval$ID))
n1235=length(intersect(intersect(intersect(NoD_pval$ID,DwoC_pval$ID),DwC_pval$ID),All_pval$ID))
n1245=length(intersect(intersect(intersect(NoD_pval$ID,DwoC_pval$ID),Diabetes_pval$ID),All_pval$ID))
n1345=length(intersect(intersect(intersect(NoD_pval$ID,DwC_pval$ID),Diabetes_pval$ID),All_pval$ID))
n2345=length(intersect(intersect(intersect(DwoC_pval$ID,DwC_pval$ID),Diabetes_pval$ID),All_pval$ID))
n12345=n1234=length(intersect(intersect(intersect(intersect(NoD_pval$ID,DwoC_pval$ID),DwC_pval$ID),Diabetes_pval$ID),All_pval$ID))

plot.new()
pdf(paste(Figs,"NoD_DwC_DwoC_pval0.05.pdf",sep=""))
# grid.newpage()
draw.triple.venn(area1=a1, area2=a2, area3=a3,n12=n12, n13=n13, n23=n23,n123=n123, category = c("NoD", "DwoC", "DwC"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"NoD_Diabetes_pval0.05.pdf",sep=""))
draw.pairwise.venn(area1=a1, area2=a4, cross.area=n14, category = c("NoD", "Diabetes"),lty = "solid", fill = c("skyblue", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"NoD_DwC_DwoC_Diabetes_pval0.05.pdf",sep=""))
draw.quad.venn(area1=a1, area2=a2, area3=a3, area4=a4, n12, n13, n14, n23, n24, n34, n123, n124, n134,
                    n234, n1234, category = c("NoD", "DwoC", "DwC", "Diabetes"),
                    lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()

plot.new()
pdf(paste(Figs,"NoD_DwC_DwoC_Diabetes_All_pval0.05.pdf",sep=""),width=12)
draw.quintuple.venn(area1=a1, area2=a2, area3=a3, area4=a4, area5=a5, n12, n13, n14, n15, n23, n24, n25, n34, n35, n45, n123, n124, n125, n134,
                    n135, n145, n234, n235, n245, n345, n1234, n1235,n1245, n1345, n2345, n12345, category = c("NoD", "DwoC", "DwC", "Diabetes", "All"),
                    lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow","orange"))
dev.off()

# q-val <0.1
a1<-length(NoD_qval_10$ID)
a2<-length(DwoC_qval_10$ID)
a3<-length(DwC_qval_10$ID)
a4<-length(Diabetes_qval_10$ID)
a5<-length(All_qval_10$ID)
n12=length(intersect(NoD_qval_10$ID,DwoC_qval_10$ID))
n13=length(intersect(NoD_qval_10$ID,DwC_qval_10$ID))
n14=length(intersect(NoD_qval_10$ID,Diabetes_qval_10$ID))
n15=length(intersect(NoD_qval_10$ID,All_qval_10$ID))
n23=length(intersect(DwoC_qval_10$ID,DwC_qval_10$ID))
n24=length(intersect(DwoC_qval_10$ID,Diabetes_qval_10$ID))
n25=length(intersect(DwoC_qval_10$ID,All_qval_10$ID))
n34=length(intersect(DwC_qval_10$ID,Diabetes_qval_10$ID))
n35=length(intersect(DwC_qval_10$ID,All_qval_10$ID))
n45=length(intersect(Diabetes_qval_10$ID,All_qval_10$ID))

n123=length(intersect(intersect(NoD_qval_10$ID,DwoC_qval_10$ID),DwC_qval_10$ID))
n124=length(intersect(intersect(NoD_qval_10$ID,DwoC_qval_10$ID),Diabetes_qval_10$ID))
n125=length(intersect(intersect(NoD_qval_10$ID,DwoC_qval_10$ID),All_qval_10$ID))
n134=length(intersect(intersect(NoD_qval_10$ID,DwC_qval_10$ID),Diabetes_qval_10$ID))
n135=length(intersect(intersect(NoD_qval_10$ID,DwC_qval_10$ID),All_qval_10$ID))
n145=length(intersect(intersect(NoD_qval_10$ID,Diabetes_qval_10$ID),All_qval_10$ID))
n234=length(intersect(intersect(DwoC_qval_10$ID,DwC_qval_10$ID),Diabetes_qval_10$ID))
n235=length(intersect(intersect(DwoC_qval_10$ID,DwC_qval_10$ID),All_qval_10$ID))
n245=length(intersect(intersect(DwoC_qval_10$ID,Diabetes_qval_10$ID),All_qval_10$ID))
n345=length(intersect(intersect(DwC_qval_10$ID,Diabetes_qval_10$ID),All_qval_10$ID))
n1234=length(intersect(intersect(intersect(NoD_qval_10$ID,DwoC_qval_10$ID),DwC_qval_10$ID),Diabetes_qval_10$ID))
n1235=length(intersect(intersect(intersect(NoD_qval_10$ID,DwoC_qval_10$ID),DwC_qval_10$ID),All_qval_10$ID))
n1245=length(intersect(intersect(intersect(NoD_qval_10$ID,DwoC_qval_10$ID),Diabetes_qval_10$ID),All_qval_10$ID))
n1345=length(intersect(intersect(intersect(NoD_qval_10$ID,DwC_qval_10$ID),Diabetes_qval_10$ID),All_qval_10$ID))
n2345=length(intersect(intersect(intersect(DwoC_qval_10$ID,DwC_qval_10$ID),Diabetes_qval_10$ID),All_qval_10$ID))
n12345=n1234=length(intersect(intersect(intersect(intersect(NoD_qval_10$ID,DwoC_qval_10$ID),DwC_qval_10$ID),Diabetes_qval_10$ID),All_qval_10$ID))

# grid.newpage()
plot.new()
pdf(paste(Figs,"NoD_DwC_DwoC_qval0.1.pdf",sep=""))
draw.triple.venn(area1=a1, area2=a2, area3=a3,n12=n12, n13=n13, n23=n23,n123=n123, category = c("NoD", "DwoC", "DwC"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"NoD_Diabetes_qval0.1.pdf",sep=""))
draw.pairwise.venn(area1=a1, area2=a4, cross.area=n14, category = c("NoD", "Diabetes"),lty = "solid", fill = c("skyblue", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"NoD_DwC_DwoC_Diabetes_qval0.1.pdf",sep=""))
draw.quad.venn(area1=a1, area2=a2, area3=a3, area4=a4, n12, n13, n14, n23, n24, n34, n123, n124, n134,
                    n234, n1234, category = c("NoD", "DwoC", "DwC", "Diabetes"),
                    lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()

# grid.newpage()
plot.new()
pdf(paste(Figs,"NoD_DwC_DwoC_Diabetes_All_qval0.1.pdf",sep=""),width=12)
draw.quintuple.venn(area1=a1, area2=a2, area3=a3, area4=a4, area5=a5, n12, n13, n14, n15, n23, n24, n25, n34, n35, n45, n123, n124, n125, n134,
                    n135, n145, n234, n235, n245, n345, n1234, n1235,n1245, n1345, n2345, n12345, category = c("NoD", "DwoC", "DwC", "Diabetes", "All"),
                    lty = "blank", fill = c("skyblue", "pink1", "mediumorchid", "yellow","orange"))
dev.off()

# q-val <0.2
a1<-length(NoD_qval_20$ID)
a2<-length(DwoC_qval_20$ID)
a3<-length(DwC_qval_20$ID)
a4<-length(Diabetes_qval_20$ID)
a5<-length(All_qval_20$ID)
n12=length(intersect(NoD_qval_20$ID,DwoC_qval_20$ID))
n13=length(intersect(NoD_qval_20$ID,DwC_qval_20$ID))
n14=length(intersect(NoD_qval_20$ID,Diabetes_qval_20$ID))
n15=length(intersect(NoD_qval_20$ID,All_qval_20$ID))
n23=length(intersect(DwoC_qval_20$ID,DwC_qval_20$ID))
n24=length(intersect(DwoC_qval_20$ID,Diabetes_qval_20$ID))
n25=length(intersect(DwoC_qval_20$ID,All_qval_20$ID))
n34=length(intersect(DwC_qval_20$ID,Diabetes_qval_20$ID))
n35=length(intersect(DwC_qval_20$ID,All_qval_20$ID))
n45=length(intersect(Diabetes_qval_20$ID,All_qval_20$ID))

n123=length(intersect(intersect(NoD_qval_20$ID,DwoC_qval_20$ID),DwC_qval_20$ID))
n124=length(intersect(intersect(NoD_qval_20$ID,DwoC_qval_20$ID),Diabetes_qval_20$ID))
n125=length(intersect(intersect(NoD_qval_20$ID,DwoC_qval_20$ID),All_qval_20$ID))
n134=length(intersect(intersect(NoD_qval_20$ID,DwC_qval_20$ID),Diabetes_qval_20$ID))
n135=length(intersect(intersect(NoD_qval_20$ID,DwC_qval_20$ID),All_qval_20$ID))
n145=length(intersect(intersect(NoD_qval_20$ID,Diabetes_qval_20$ID),All_qval_20$ID))
n234=length(intersect(intersect(DwoC_qval_20$ID,DwC_qval_20$ID),Diabetes_qval_20$ID))
n235=length(intersect(intersect(DwoC_qval_20$ID,DwC_qval_20$ID),All_qval_20$ID))
n245=length(intersect(intersect(DwoC_qval_20$ID,Diabetes_qval_20$ID),All_qval_20$ID))
n345=length(intersect(intersect(DwC_qval_20$ID,Diabetes_qval_20$ID),All_qval_20$ID))
n1234=length(intersect(intersect(intersect(NoD_qval_20$ID,DwoC_qval_20$ID),DwC_qval_20$ID),Diabetes_qval_20$ID))
n1235=length(intersect(intersect(intersect(NoD_qval_20$ID,DwoC_qval_20$ID),DwC_qval_20$ID),All_qval_20$ID))
n1245=length(intersect(intersect(intersect(NoD_qval_20$ID,DwoC_qval_20$ID),Diabetes_qval_20$ID),All_qval_20$ID))
n1345=length(intersect(intersect(intersect(NoD_qval_20$ID,DwC_qval_20$ID),Diabetes_qval_20$ID),All_qval_20$ID))
n2345=length(intersect(intersect(intersect(DwoC_qval_20$ID,DwC_qval_20$ID),Diabetes_qval_20$ID),All_qval_20$ID))
n12345=n1234=length(intersect(intersect(intersect(intersect(NoD_qval_20$ID,DwoC_qval_20$ID),DwC_qval_20$ID),Diabetes_qval_20$ID),All_qval_20$ID))

# grid.newpage()
plot.new()
pdf(paste(Figs,"NoD_DwC_DwoC_qval0.2.pdf",sep=""))
draw.triple.venn(area1=a1, area2=a2, area3=a3,n12=n12, n13=n13, n23=n23,n123=n123, category = c("NoD", "DwoC", "DwC"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"NoD_Diabetes_qval0.2.pdf",sep=""))
draw.pairwise.venn(area1=a1, area2=a4, cross.area=n14, category = c("NoD", "Diabetes"),lty = "solid", fill = c("skyblue", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"NoD_DwC_DwoC_Diabetes_qval0.2.pdf",sep=""))
draw.quad.venn(area1=a1, area2=a2, area3=a3, area4=a4, n12, n13, n14, n23, n24, n34, n123, n124, n134,
                    n234, n1234, category = c("NoD", "DwoC", "DwC", "Diabetes"),
                    lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()

# grid.newpage()
plot.new()
pdf(paste(Figs,"NoD_DwC_DwoC_Diabetes_All_qval0.2.pdf",sep=""),width=12)
draw.quintuple.venn(area1=a1, area2=a2, area3=a3, area4=a4, area5=a5, n12, n13, n14, n15, n23, n24, n25, n34, n35, n45, n123, n124, n125, n134,
                    n135, n145, n234, n235, n245, n345, n1234, n1235,n1245, n1345, n2345, n12345, category = c("NoD", "DwoC", "DwC", "Diabetes", "All"),
                    lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow","orange"))
dev.off()

# q-val <0.5
a1<-length(NoD_qval_50$ID)
a2<-length(DwoC_qval_50$ID)
a3<-length(DwC_qval_50$ID)
a4<-length(Diabetes_qval_50$ID)
a5<-length(All_qval_50$ID)
n12=length(intersect(NoD_qval_50$ID,DwoC_qval_50$ID))
n13=length(intersect(NoD_qval_50$ID,DwC_qval_50$ID))
n14=length(intersect(NoD_qval_50$ID,Diabetes_qval_50$ID))
n15=length(intersect(NoD_qval_50$ID,All_qval_50$ID))
n23=length(intersect(DwoC_qval_50$ID,DwC_qval_50$ID))
n24=length(intersect(DwoC_qval_50$ID,Diabetes_qval_50$ID))
n25=length(intersect(DwoC_qval_50$ID,All_qval_50$ID))
n34=length(intersect(DwC_qval_50$ID,Diabetes_qval_50$ID))
n35=length(intersect(DwC_qval_50$ID,All_qval_50$ID))
n45=length(intersect(Diabetes_qval_50$ID,All_qval_50$ID))

n123=length(intersect(intersect(NoD_qval_50$ID,DwoC_qval_50$ID),DwC_qval_50$ID))
n124=length(intersect(intersect(NoD_qval_50$ID,DwoC_qval_50$ID),Diabetes_qval_50$ID))
n125=length(intersect(intersect(NoD_qval_50$ID,DwoC_qval_50$ID),All_qval_50$ID))
n134=length(intersect(intersect(NoD_qval_50$ID,DwC_qval_50$ID),Diabetes_qval_50$ID))
n135=length(intersect(intersect(NoD_qval_50$ID,DwC_qval_50$ID),All_qval_50$ID))
n145=length(intersect(intersect(NoD_qval_50$ID,Diabetes_qval_50$ID),All_qval_50$ID))
n234=length(intersect(intersect(DwoC_qval_50$ID,DwC_qval_50$ID),Diabetes_qval_50$ID))
n235=length(intersect(intersect(DwoC_qval_50$ID,DwC_qval_50$ID),All_qval_50$ID))
n245=length(intersect(intersect(DwoC_qval_50$ID,Diabetes_qval_50$ID),All_qval_50$ID))
n345=length(intersect(intersect(DwC_qval_50$ID,Diabetes_qval_50$ID),All_qval_50$ID))
n1234=length(intersect(intersect(intersect(NoD_qval_50$ID,DwoC_qval_50$ID),DwC_qval_50$ID),Diabetes_qval_50$ID))
n1235=length(intersect(intersect(intersect(NoD_qval_50$ID,DwoC_qval_50$ID),DwC_qval_50$ID),All_qval_50$ID))
n1245=length(intersect(intersect(intersect(NoD_qval_50$ID,DwoC_qval_50$ID),Diabetes_qval_50$ID),All_qval_50$ID))
n1345=length(intersect(intersect(intersect(NoD_qval_50$ID,DwC_qval_50$ID),Diabetes_qval_50$ID),All_qval_50$ID))
n2345=length(intersect(intersect(intersect(DwoC_qval_50$ID,DwC_qval_50$ID),Diabetes_qval_50$ID),All_qval_50$ID))
n12345=n1234=length(intersect(intersect(intersect(intersect(NoD_qval_50$ID,DwoC_qval_50$ID),DwC_qval_50$ID),Diabetes_qval_50$ID),All_qval_50$ID))

plot.new()
pdf(paste(Figs,"NoD_DwC_DwoC_qval0.5.pdf",sep=""))
draw.triple.venn(area1=a1, area2=a2, area3=a3,n12=n12, n13=n13, n23=n23,n123=n123, category = c("NoD", "DwoC", "DwC"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"NoD_Diabetes_qval0.5.pdf",sep=""))
draw.pairwise.venn(area1=a1, area2=a4, cross.area=n14, category = c("NoD", "Diabetes"),lty = "solid", fill = c("skyblue", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"NoD_DwC_DwoC_Diabetes_qval0.5.pdf",sep=""))
draw.quad.venn(area1=a1, area2=a2, area3=a3, area4=a4, n12, n13, n14, n23, n24, n34, n123, n124, n134,
                    n234, n1234, category = c("NoD", "DwoC", "DwC", "Diabetes"),
                    lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()

plot.new()
pdf(paste(Figs,"NoD_DwC_DwoC_Diabetes_All_qval0.5.pdf",sep=""),width=12)
draw.quintuple.venn(area1=a1, area2=a2, area3=a3, area4=a4, area5=a5, n12, n13, n14, n15, n23, n24, n25, n34, n35, n45, n123, n124, n125, n134,
                    n135, n145, n234, n235, n245, n345, n1234, n1235,n1245, n1345, n2345, n12345, category = c("NoD", "DwoC", "DwC", "Diabetes", "All"),
                    lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow","orange"))
dev.off()

# delta pval=0.5 & FC.cutoff=1.2
a1<-length(delta_NoD_DwoC_pval_FC$ID)
a2<-length(delta_NoD_DwC_pval_FC$ID)
a3<-length(delta_NoD_Diabetes_pval_FC$ID)
a4<-length(delta_DwC_DwoC_pval_FC$ID)
n12=length(intersect(delta_NoD_DwoC_pval_FC$ID,delta_NoD_DwC_pval_FC$ID))
n13=length(intersect(delta_NoD_DwoC_pval_FC$ID,delta_NoD_Diabetes_pval_FC$ID))
n14=length(intersect(delta_NoD_DwoC_pval_FC$ID,delta_DwC_DwoC_pval_FC$ID))
n23=length(intersect(delta_NoD_DwC_pval_FC$ID,delta_NoD_Diabetes_pval_FC$ID))
n24=length(intersect(delta_NoD_DwC_pval_FC$ID,delta_DwC_DwoC_pval_FC$ID))
n34=length(intersect(delta_NoD_Diabetes_pval_FC$ID,delta_DwC_DwoC_pval_FC$ID))

n123=length(intersect(intersect(delta_NoD_DwoC_pval_FC$ID,delta_NoD_DwC_pval_FC$ID),delta_NoD_Diabetes_pval_FC$ID))
n124=length(intersect(intersect(delta_NoD_DwoC_pval_FC$ID,delta_NoD_DwC_pval_FC$ID),delta_DwC_DwoC_pval_FC$ID))
n125=length(intersect(intersect(delta_NoD_DwoC_pval_FC$ID,delta_NoD_DwC_pval_FC$ID),All_pval$ID))
n134=length(intersect(intersect(delta_NoD_DwoC_pval_FC$ID,delta_NoD_Diabetes_pval_FC$ID),delta_DwC_DwoC_pval_FC$ID))
n234=length(intersect(intersect(delta_NoD_DwC_pval_FC$ID,delta_NoD_Diabetes_pval_FC$ID),delta_DwC_DwoC_pval_FC$ID))
n1234=length(intersect(intersect(intersect(delta_NoD_DwoC_pval_FC$ID,delta_NoD_DwC_pval_FC$ID),delta_NoD_Diabetes_pval_FC$ID),delta_DwC_DwoC_pval_FC$ID))


plot.new()
pdf(paste(Figs,"delta_NoD_DwoC_and_delta_NoD_DwC_pval0.05_FC1.2.pdf",sep=""))
draw.pairwise.venn(area1=a1, area2=a2, cross.area=n12, category = c("NoD_vs_DwoC", "NoD_vs_DwC"),lty = "solid", fill = c("skyblue", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"delta_NoD_DwoC_and_delta_NoD_DwC_and_NoD_Diabetes_pval0.05_FC1.2.pdf",sep=""))
draw.triple.venn(area1=a1, area2=a2, area3=a3,n12=n12, n13=n13, n23=n23,n123=n123, category = c("NoD_vs_DwoC", "NoD_vs_DwC", "NoD_vs_Diabetes"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()


plot.new()
pdf(paste(Figs,"delta_all_group_comparison_pval0.05_FC1.2.pdf",sep=""),width=12)

draw.quad.venn(area1=a1, area2=a2, area3=a3, area4=a4,  n12=n12, n13=n13, n14=n14, n23=n23, n24=n24,
               n34=n34, n123=n123, n124=n124, n134=n134, n234=n234, n1234=n1234, category = c("NoD_vs_DwoC", "NoD_vs_DwC", "NoD_vs_Diabetes", "DwC_vs_DwoC"),
                 lty = "blank", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()
