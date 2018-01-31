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
source(paste(HOME,"script/andrew_Func.R",sep=""))

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

if (require(lumiHumanIDMapping)) {
  probeID<-nuID2probeID(rownames(dat), lib.mapping = "lumiHumanIDMapping")
}


genes <-data.frame(ID=probeList, probeID=probeID,geneSymbol=geneSymbol, geneName=geneName, stringsAsFactors=FALSE)

# PCA by group (HG, SG, delta)
# 1. No_DM
subject<-"No_DM"
NoD<-dat[,grep("NoD",colnames(dat))] # grep all normal samples
NoD_sg<-NoD[,grep("norm",colnames(NoD))]; 
NoD_hg<-NoD[,grep("30mM",colnames(NoD))];
NoD_delta <- NoD_hg-NoD_sg

NoD_sg_avg<-prep_data(NoD_sg)
NoD_hg_avg<-prep_data(NoD_hg)
NoD_delta_avg<-prep_data(NoD_delta)
rm(NoD_hg, NoD_sg, NoD_delta)

# 2. No_PDR
subject<-"No_PDR"
DwoC<-dat[,grep("DwoC",colnames(dat))] # grep all normal samples
DwoC_sg<-DwoC[,grep("norm",colnames(DwoC))]; DwoC_hg<-DwoC[,grep("30mM",colnames(DwoC))];
DwoC_delta <- DwoC_hg-DwoC_sg

DwoC_sg_avg<-prep_data(DwoC_sg)
DwoC_hg_avg<-prep_data(DwoC_hg)
DwoC_delta_avg<-prep_data(DwoC_delta)
rm(DwoC_hg, DwoC_sg, DwoC_delta)

# 3. PDR
DwC<-dat[,grep("DwC",colnames(dat))] # grep all normal samples
DwC_sg<-DwC[,grep("norm",colnames(DwC))]; DwC_hg<-DwC[,grep("30mM",colnames(DwC))];
DwC_delta <- DwC_hg-DwC_sg

DwC_sg_avg<-prep_data(DwC_sg)
DwC_hg_avg<-prep_data(DwC_hg)
DwC_delta_avg<-prep_data(DwC_delta)
rm(DwC_hg, DwC_sg, DwC_delta)


##----------------------------------
## group comparison (each hg, sg, and delta)
##----------------------------------

#1a. delta
delta<-cbind(DwC_delta_avg,DwoC_delta_avg,NoD_delta_avg)
subject="all_replicate"
targets<-readTargets(paste(PhenotypeDir,"avg_delta_",subject,"_target.txt", sep=''))
group <- factor(targets$Group)
design <- model.matrix(~0+group)
colnames(design)<-levels(group)
rownames(design)<-targets$SampleName
fit <- lmFit(delta, design)

cm <- makeContrasts( delta_NoD_DwoC=NoD - DwoC, delta_NoD_DwC=NoD - DwC, 
                     delta_NoD_Diabetes=NoD - (DwC+DwoC)/2, delta_DwoC_DwC = DwoC-DwC,
                     levels=design) # create a contrast matrix
fit2=contrasts.fit(fit,cm)
fit2=eBayes(fit2)

setEPS()
postscript(file = paste(Figs,"avg_delta_all_group_pvalues_",date,".eps",sep=""))

par(mfrow=c(2,2))
for (i in 1:ncol(fit2$p.value)) {
  hist(fit2$p.value[,i], main=colnames(fit2$p.value)[i])
}
dev.off()

#----------------
# delta NoD vs DwoC
#----------------

qval.cutoff=0.05; FC.cutoff=0.263 # FC=1.2

x1=topTable(fit2, coef="delta_NoD_DwoC", n=nrow(genes), adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="delta_NoD_DwoC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$P.Value<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"avg_delta_NoD_DwoC_pvalue_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(P.Value), pch=20, col="gray", main=paste("Group difference in delta (No_DM vs. No_PDR)",sep=""), xlab="log2 fold change", ylab="-log10(P.Value)", xlim=c(-1.0,1.0)))
with(subset(y1, P.Value < qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(y1, P.Value < qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()
write.csv(y1, paste(Output,"avg_delta_NoD_DwoC_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"avg_delta_NoD_DwoC_pvalue_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#----------------
# delta NoD vs DwC
#----------------
x1=topTable(fit2, coef="delta_NoD_DwC", n=nrow(genes), adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="delta_NoD_DwC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$P.Value<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"avg_delta_NoD_DwC_pvalue_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(P.Value), pch=20, col="gray", main=paste("Group difference in delta (No_DM vs. PDR)",sep=""), xlab="log2 fold change", ylab="-log10(P.Value)", xlim=c(-1.0,1.0)))
with(subset(y1, P.Value<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(y1, P.Value<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()

write.csv(y1, paste(Output,"avg_delta_NoD_DwC_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"avg_delta_NoD_DwC_pvalue_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#---------------------
# delta NoD vs DwC+DwoC
#---------------------
x1=topTable(fit2, coef="delta_NoD_Diabetes", n=nrow(genes), adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="delta_NoD_Diabetes", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$P.Value<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"avg_delta_NoD_Diabetes_pvalue_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(P.Value), pch=20, col="gray", main=paste("Group difference in delta (No_DM vs. PDR + No_PDR)",sep=""), xlab="log2 fold change", ylab="-log10(P.Value)", xlim=c(-1.0,1.0)))
with(subset(y1, P.Value<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(y1, P.Value<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()

write.csv(y1, paste(Output,"avg_delta_NoD_Diabetes_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"avg_delta_NoD_Diabetes_pvalue_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#---------------------
# delta DwC vs DwoC
#---------------------
x1=topTable(fit2, coef="delta_DwoC_DwC", n=nrow(genes), adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="delta_DwoC_DwC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$P.Value<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"avg_delta_DwC_DwoC_pvalue_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(P.Value), pch=20, col="gray", main=paste("Group difference in delta (PDR vs. No_PDR)",sep=""), xlab="log2 fold change", ylab="-log10(P.Value)", xlim=c(-1,1)))
with(subset(y1, P.Value<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(y1, P.Value<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()

write.csv(y1, paste(Output,"avg_delta_DwC_DwoC_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"avg_delta_DwC_DwoC_pvalue_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)


##-----------------------------------------
#  high glucose
# ##-----------------------------------------
#1a. No_DM vs PDR + No_PDR
subject="No_DM_vs_Diabetes"
subject="all_replicate"
targets<-readTargets(paste(PhenotypeDir,"hg_sg_",subject,"_target.txt", sep=''))
group <- factor(paste(targets$Group,targets$Treatment,sep="."))
design <- model.matrix(~0+group)
colnames(design)<-levels(group)
rownames(design)<-targets$SampleName

hg<-cbind(DwC_hg_avg,DwoC_hg_avg,NoD_hg_avg)
fit <- lmFit(hg, design)

hg_cm <- makeContrasts( hg_NoD_DwoC=NoD - DwoC, hg_NoD_DwC=NoD - DwC,
                     hg_NoD_Diabetes=NoD - (DwC+DwoC)/2, hg_DwoC_DwC = DwoC-DwC,
                     levels=design) # create a contrast matrix

cm <- makeContrasts( hg_NoD_DwoC=NoD.T - DwoC.T, hg_NoD_DwC=NoD.T - DwC.T,
                     hg_NoD_Diabetes=NoD.T - (DwC.T+DwoC.T)/2, hg_DwC_DwoC = DwC.T-DwoC.T,
                     sg_NoD_DwoC=NoD.C - DwoC.C, sg_NoD_DwC=NoD.C - DwC.C,
                     sg_NoD_Diabetes=NoD.C - (DwC.C+DwoC.C)/2, sg_DwC_DwoC = DwC.C-DwoC.C,
                     levels=design) # create a contrast matrix
fit2=contrasts.fit(fit,hg_cm)
fit2=eBayes(fit2)

setEPS()
postscript(file = paste(Figs,"avg_hg_group_pvalues_",date,".eps",sep=""))

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
postscript(file = paste(Figs,"avg_hg_NoD_DwoC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
# pdf(file = paste(Figs,"hg_NoD_DwoC_q",qval.cutoff,"_","FC",FC.cutoff,"_",date,".pdf",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in high glucose (No_DM vs. No_PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"avg_hg_NoD_DwoC_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"avg_hg_NoD_DwoC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#----------------
# hg NoD vs DwC
#----------------
x1=topTable(fit2, coef="hg_NoD_DwC", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="hg_NoD_DwC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"avg_hg_NoD_DwC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
# pdf(file = paste(Figs,"hg_NoD_DwC_q",qval.cutoff,"_","FC",FC.cutoff,"_",date,".pdf",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in high glucose (No_DM vs. PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"avg_hg_NoD_DwC_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"avg_hg_NoD_DwC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#---------------------
# hg NoD vs DwC+DwoC
#---------------------
x1=topTable(fit2, coef="hg_NoD_Diabetes", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="hg_NoD_Diabetes", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"avg_hg_NoD_Diabetes_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
# pdf(file = paste(Figs,"hg_NoD_Diabetes_q",qval.cutoff,"_","FC",FC.cutoff,"_",date,".pdf",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in high glucose (No_DM vs. PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"avg_hg_NoD_Diabetes_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"avg_hg_NoD_Diabetes_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#---------------------
# hg DwC vs DwoC
#---------------------
x1=topTable(fit2, coef="hg_DwC_DwoC", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="hg_DwC_DwoC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"avg_hg_DwC_DwoC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
# pdf(file = paste(Figs,"hg_DwC_DwoC_q",qval.cutoff,"_","FC",FC.cutoff,"_",date,".pdf",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in high glucose (PDR vs. No_PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"avg_hg_DwC_DwoC_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"avg_hg_DwC_DwoC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

##-----------------------------------------
#  standard glucose
##-----------------------------------------
sg<-cbind(DwC_sg_avg,DwoC_sg_avg,NoD_sg_avg)

fit <- lmFit(sg, design)

sg_cm <- makeContrasts( sg_NoD_DwoC=NoD - DwoC, sg_NoD_DwC=NoD - DwC,
                        sg_NoD_Diabetes=NoD - (DwC+DwoC)/2, sg_DwoC_DwC = DwoC-DwC,
                        levels=design) # create a contrast matrix

fit2=contrasts.fit(fit,sg_cm)
fit2=eBayes(fit2)

setEPS()
postscript(file = paste(Figs,"avg_sg_group_pvalues_",date,".eps",sep=""))

par(mfrow=c(2,4))
for (i in 1:ncol(fit2$p.value)) {
  hist(fit2$p.value[,i], main=colnames(fit2$p.value)[i])
}
dev.off()
#----------------
# sg NoD vs DwoC
#----------------

qval.cutoff=0.05; FC.cutoff=1 # FC=1.12

x1=topTable(fit2, coef="sg_NoD_DwoC", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="sg_NoD_DwoC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"avg_sg_NoD_DwoC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in standard glucose (No_DM vs. No_PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"avg_sg_NoD_DwoC_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"avg_sg_NoD_DwoC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#----------------
# sg NoD vs DwC
#----------------
x1=topTable(fit2, coef="sg_NoD_DwC", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="sg_NoD_DwC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"avg_sg_NoD_DwC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in standard glucose (No_DM vs. PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"avg_sg_NoD_DwC_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"avg_sg_NoD_DwC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#-------------------
# sg NoD vs DwC+DwoC
#-------------------
x1=topTable(fit2, coef="sg_NoD_Diabetes", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="sg_NoD_Diabetes", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"avg_sg_NoD_Diabetes_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in standard glucose (No_DM vs. No_PDR + PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"avg_sg_NoD_Diabetes_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"avg_sg_NoD_Diabetes_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#---------------------
# sg DwC vs DwoC
#---------------------
x1=topTable(fit2, coef="sg_DwC_DwoC", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit2, coef="sg_DwC_DwoC", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"avg_sg_DwC_DwoC_q",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in standard glucose (PDR vs. No_PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(x1, paste(Output,"avg_sg_DwC_DwoC_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"avg_sg_DwC_DwoC_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)


##--------------------------------------------------
## treatment response high vs standard for each group
##--------------------------------------------------
#-----------------
# 1. No_DM (NoD)
#-----------------
avg_NoD<-cbind(NoD_hg_avg,NoD_sg_avg)
subject="No_DM_replicate"
targets<-readTargets(paste(PhenotypeDir,"avg_hg_sg_",subject,"_target.txt", sep=''))
Paired <- factor(targets$paired)
Treat <- factor(targets$Treatment)
design <- model.matrix(~Paired+Treat)
fit <- lmFit(avg_NoD, design)
fit<- eBayes(fit)

qval.cutoff=0.1; FC.cutoff=0.17 # FC=1.12

x1=topTable(fit, coef="TreatT", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="TreatT", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"avg_NoD_treatment_effect_pval_",qval.cutoff,"_logFC_",FC.cutoff ,date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("No_DM treatment effect (high vs. standard glucose)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-.7,.7)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

cat("number of genes (q-val <0.1): ",sum(x1$adj.P.Val<qval.cutoff))
write.csv(x1_subset, paste(Output,"avg_NoD_treatment_effect_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1, paste(Output,"avg_NoD_treatment_effect_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(y1, paste(Output,"avg_NoD_treatment_effect_",date,".csv",sep=""),row.names = F)

setEPS()
postscript(file = paste(Figs,"avg_NoD_treatment_effect_pvalue_distribution_",date,".eps",sep=""),width=18,height=12,pointsize=20)
par(mfrow=c(1,2))
hist(y1$P.Value, main="",xlab="p-value")
hist(y1$adj.P.Val,main="", xlab="adjusted p-value")
mtext("avg_No_DM treatment effect (high vs. standard glucose)", side=3, outer=TRUE, line=-3)
dev.off()

# 2. No_PDR (DwoC)
subject="No_PDR_replicate"
avg_DwoC<-cbind(DwoC_hg_avg,DwoC_sg_avg) 
targets<-readTargets(paste(PhenotypeDir,"avg_hg_sg_",subject,"_target.txt", sep=''))
Paired <- factor(targets$paired)
Treat <- factor(targets$Treatment)
design <- model.matrix(~Paired+Treat)
fit <- lmFit(avg_DwoC, design)
fit<- eBayes(fit)

qval.cutoff=0.4; FC.cutoff=0.17 # FC=1.12

x1=topTable(fit, coef="TreatT", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="TreatT", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

setEPS()
postscript(file = paste(Figs,"avg_DwoC_treatment_effect_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("No_PDR treatment effect (high vs. standard glucose)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-.5,.5)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()
# with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
# with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))

cat("number of genes (q-val <0.1): ",sum(x1$adj.P.Val<qval.cutoff))
write.csv(x1_subset, paste(Output,"avg_DwoC_treatment_effect_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1, paste(Output,"avg_DwoC_treatment_effect_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(y1, paste(Output,"avg_DwoC_treatment_effect_",date,".csv",sep=""),row.names = F)

setEPS()
postscript(file = paste(Figs,"avg_DwoC_treatment_effect_pvalue_distribution_",date,".eps",sep=""),width=18,height=12,pointsize=20)
par(mfrow=c(1,2))
hist(y1$P.Value, main="",xlab="p-value")
hist(y1$adj.P.Val,main="", xlab="adjusted p-value")
mtext("avg_No_PDR treatment effect (high vs. standard glucose)", side=3, outer=TRUE, line=-3)
dev.off()


#-----------------
#  STOP HERE #
#-----------------




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

setEPS()
postscript(file = paste(Figs,"Diabetes_treatment_effect_qval_",qval.cutoff,"_logFC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("No_PDR and PDR treatment effect (high vs. standard glucose)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-0.5,0.5)))
# with(subset(y1, adj.P.Val<qval.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
# with(subset(y1, adj.P.Val<qval.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))

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

setEPS()
postscript(file = paste(Figs,"All_treatment_effect_qval_",qval.cutoff,"_logFC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("No_DM, No_PDR and PDR treatment effect (high vs. standard glucose)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-0.5,0.5)))
# with(subset(y1, adj.P.Val<qval.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
# with(subset(y1, adj.P.Val<qval.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))

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
