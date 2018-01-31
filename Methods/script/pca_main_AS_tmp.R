library(ggplot2)
library(reshape2);

HOME <- "/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/"
PhenotypeDir = paste(HOME,'data/phenotype/', sep='')
RData=paste(HOME,'rdata/',sep=''); dir.create(RData, showWarnings = FALSE)
Figs = paste(HOME,'Figs/', sep=''); dir.create(Figs, showWarnings = FALSE)

setwd(HOME)
source(paste(HOME,"script/pca.R",sep=""))
source(paste(HOME,"script/Function.R",sep=""))
source(paste(HOME,"script/andrew_Func.R",sep=""))

library("ggplot2")
library("gridExtra")
library("cowplot")
library("pscl")

date=Sys.Date()

# load data
dat<-local(get(load(file=paste(RData,"collapsed_normalizedDataMatrix.sorted.RData",sep=""))))
#---------------------------------------------------------------------------
# PCA by group (HG, SG, delta)
# 1. No_DM (NoD)
subject<-"No_DM"
NoD<-dat[,grep("NoD",colnames(dat))] # grep all normal samples
NoD_sg<-NoD[,grep("norm",colnames(NoD))]; NoD_sg<-apply(NoD_sg,2,as.numeric);rownames(NoD_sg)<-rownames(dat);
NoD_hg<-NoD[,grep("30mM",colnames(NoD))]; NoD_hg<-apply(NoD_hg,2,as.numeric);rownames(NoD_hg)<-rownames(dat);
NoD_delta <- NoD_hg-NoD_sg

# 2. No_PDR (DwoC)
subject<-"No_PDR"
DwoC<-dat[,grep("DwoC",colnames(dat))] # grep all normal samples
DwoC_sg<-DwoC[,grep("norm",colnames(DwoC))]; DwoC_sg<-apply(DwoC_sg,2,as.numeric);rownames(DwoC_sg)<-rownames(dat);
DwoC_hg<-DwoC[,grep("30mM",colnames(DwoC))]; DwoC_hg<-apply(DwoC_hg,2,as.numeric);rownames(DwoC_hg)<-rownames(dat);
DwoC_delta <- DwoC_hg-DwoC_sg

# 3. PDR (DwC)
DwC<-dat[,grep("DwC",colnames(dat))] # grep all normal samples
DwC_sg<-DwC[,grep("norm",colnames(DwC))]; DwC_sg<-apply(DwC_sg,2,as.numeric);rownames(DwC_sg)<-rownames(dat);
DwC_hg<-DwC[,grep("30mM",colnames(DwC))]; DwC_hg<-apply(DwC_hg,2,as.numeric);rownames(DwC_hg)<-rownames(dat);
DwC_delta <- DwC_hg-DwC_sg

######
meta<-read.table(paste(PhenotypeDir, "meanDeltaSampleInfo3.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
group<-meta$GROUP

meta1<-read.table(paste(PhenotypeDir, "sample_group.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
group_all<-meta1$GROUP
## ------------------------------------------------

# 1a. HG - group removed
All_HG <- cbind(DwC_hg, DwoC_hg, NoD_hg)

data_input = All_HG

data_input2<-matrix(1,nrow(data_input),ncol(data_input))

for (i in 1:nrow(data_input)){
  gene<-data_input[i,]
  data_input2[i,]<-summary(lm(gene~group_all))$residual
}

rownames(data_input2)<-rownames(data_input)
colnames(data_input2)<-colnames(data_input)

data_input<-data_input2
colnames(data_input)<-substr(colnames(data_input),1,nchar(colnames(data_input))-5)
rm(data_input2)

# 1b. HG - PCA variance explained
pca_hg <- prcomp(data_input, scale = FALSE, center = TRUE)
pca_hg = pca_hg$rotation
pca_hg = as.data.frame(pca_hg)

plot_pca(pca_hg, "PCA_HG")

#--------------------------------------
## DE analysis of high glucose
meta1<-read.table(paste(PhenotypeDir, "sample_group_NoD_vs_Diabetes.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
group <-as.factor(meta1$GROUP)
pc<-pca_hg$PC1
fit <- eBayes(lmFit(data_input, model.matrix(~group+pc)))
y <- topTable(fit, coef=2, number=Inf)


# subject='delta_PDR_vs_No_PDR'
# targets<-readTargets(paste(PhenotypeDir,subject,"_target.txt", sep=''))
# blocks<-factor(targets$Subject)

corfit <- duplicateCorrelation(All_DE2, block = blocks)
fit <- lmFit(All_DE2, block = blocks, cor = corfit$consensus.correlation)
fit2=eBayes(fit)

qval.cutoff=0.1; FC.cutoff=0.17 # FC=1.2
x5=topTable(fit2,n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y5=topTable(fit2, n=nrow(genes),adjust.method="BH",genelist=genes)
x5_subset<-x5[x5$adj.P.Val<qval.cutoff & abs(x5$logFC) > FC.cutoff,]

pdf(file = paste(Figs,"no_growth_rate_",subject,"_q",qval.cutoff,"_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
with(y5, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Diabetes with and without complication subjects: high glucose vs. standard glucose (q.val <",qval.cutoff,"& FC >",FC.cutoff,sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-0.5,0.5)))
with(subset(y5, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y5, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

cat("number of genes (q-val <0.2): ",sum(x5$adj.P.Val<qval.cutoff))
write.csv(x5_subset, paste(Output,"no_growth_rate_",subject,"_",qval.cutoff,"_FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x5, paste(Output,"no_growth_rate_",subject,"_",qval.cutoff,"_",date,".csv",sep=""),row.names=F)
write.csv(y5, paste(Output,"no_growth_rate_",subject,"_",date,".csv",sep=""),row.names=F)

pdf(file = paste(Figs,"no_growth_rate_",subject,"_pvalue_distribution_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
par(mfrow=c(1,2))
hist(y5$P.Value, main="",xlab="p-value")
hist(y5$adj.P.Val,main="", xlab="adjusted p-value")
mtext(paste("no growth rate ", subject ," subject (standard glucose)",sep=""), side=3, outer=TRUE, line=-3)
dev.off()

#Normal vs DwoC
#Normal vs DwC
#Normal vs DwC and DwoC
#--------------------------------------

# 1c. HG - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_hg)
mtch = match(meta.names, pca.names)
pca_hg = pca_hg[mtch,]
pca_hg = as.data.frame(pca_hg)

pca_hg_cor<-correlate_pcs(pca_hg, meta, npcs = 10, min.cor = 0)
pca_hg_cor<-as.data.frame(pca_hg_cor);

plot_cor_pca_cov(pca_hg_cor, "pca_hg_cov_corr_all_samples")

pca_hg_p = pca.meta.regress(pca_hg, meta, no.pcs = 6)
plot_P_pca_cov(pca_hg_p$Ps, "pca_hg_p_all_samples")

write.table(file = "Predict_cov_hg.txt", pca_hg_p$Pred, quote=F, col.names = TRUE, row.names=F)

######### REPEAT ABOVE BUT EXCLUDING NO_HG 

Diabete_HG <- cbind(DwC_hg, DwoC_hg)
data_input=Diabete_HG
colnames(data_input)<-substr(colnames(data_input),1,nchar(colnames(data_input))-5)

data_input2<-matrix(1,nrow(data_input),ncol(data_input))

# group removed
for (i in 1:nrow(data_input)){
  gene<-data_input[i,]
  data_input2[i,]<-summary(lm(gene~group))$residual
}

rownames(data_input2)<-rownames(data_input)
colnames(data_input2)<-colnames(data_input)

data_input<-data_input2
# 1b. HG - PCA variance explained
pca_hg_no_norm <- prcomp(data_input, scale = FALSE, center = TRUE)

# 1c. HG - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_hg_no_norm$rotation)
mtch = match(meta.names, pca.names)
pca_hg_no_norm = pca_hg_no_norm$rotation[mtch,]
pca_hg_no_norm = as.data.frame(pca_hg_no_norm)

pca_hg_no_norm_cor<-correlate_pcs(pca_hg_no_norm, meta, npcs = 10, min.cor = 0)
pca_hg_no_norm_cor<-as.data.frame(pca_hg_no_norm_cor);

plot_cor_pca_cov(pca_hg_no_norm_cor, "PCA_meta_cor_NO_NoD")

pca_hg_p_no_norm = pca.meta.regress(pca_hg_no_norm, meta, no.pcs = 6)
plot_P_pca_cov(pca_hg_p_no_norm$Ps, "pca_hg_p_NO_NoD")

write.table(file = "Predict_cov_hg_No_NoD.txt", pca_hg_p_no_norm$Pred, quote=F, col.names = TRUE, row.names=F)

### COMPARE PCAs OF DIABETES WHEN NORMAL IS INCLUDED AND EXCLUDED.
PC_cor = calc_pc_cors(pca_hg, pca_hg_no_norm)
write.table(file = "PC_cor_hg_allSamps_vs_noNoD.txt", PC_cor, quote=F, col.names=TRUE, row.names=T)

var.select = pca.meta.regress.2var(pca_hg_no_norm, meta)
write.table(file = "Variable_select_hg_no_NoD.txt", var.select, quote=FALSE, col.names=TRUE, row.names=FALSE)


#####---------------########
# 2a. SG - collapse data
Diabete_SG<-cbind(DwC_sg, DwoC_sg, NoD_sg)
data_input = Diabete_SG
colnames(data_input)<-substr(colnames(data_input),1,nchar(colnames(data_input))-5)

data_input2<-matrix(1,nrow(data_input),ncol(data_input))
for (i in 1:nrow(data_input)){
  gene<-data_input[i,]
  data_input2[i,]<-summary(lm(gene~group_all))$residual
}

rownames(data_input2)<-rownames(data_input)
colnames(data_input2)<-colnames(data_input)

data_input<-data_input2

# 2b. SG - PCA variance explained
pca_sg <- prcomp(data_input, scale = FALSE, center = TRUE)
pca_sg = pca_sg$rotation
pca_sg = as.data.frame(pca_sg)

plot_pca(pca_sg, "PCA_SG")

# 2c. SG - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_sg)
mtch = match(meta.names, pca.names)
pca_sg = pca_sg[mtch,]
pca_sg = as.data.frame(pca_sg)

pca_sg_cor<-correlate_pcs(pca_sg, meta, npcs = 10, min.cor = 0)
pca_sg_cor<-as.data.frame(pca_sg_cor);

plot_cor_pca_cov(pca_sg_cor, "pca_sg_cov_corr_all_samples")

pca_sg_p = pca.meta.regress(pca_sg, meta)
plot_P_pca_cov(pca_sg_p$Ps, "pca_sg_p_all_samples")

write.table(file = "Predict_cov_sg.txt", pca_sg_p$Pred, quote=F, col.names = TRUE, row.names=F)


######### REPEAT ABOVE BUT EXCLUDING NO_NoD
Diabete_SG <- cbind(DwC_sg, DwoC_sg)
data_input<-Diabete_SG
colnames(data_input)<-substr(colnames(data_input),1,nchar(colnames(data_input))-5)

data_input2<-matrix(1,nrow(data_input),ncol(data_input))
for (i in 1:nrow(data_input)){
  gene<-data_input[i,]
  data_input2[i,]<-summary(lm(gene~group))$residual
}

rownames(data_input2)<-rownames(data_input)
colnames(data_input2)<-colnames(data_input)

data_input<-data_input2

# 2b-2. HG - PCA variance explained
pca_sg_no_norm <- prcomp(data_input, scale = FALSE, center = TRUE)

# 2c. HG - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_sg_no_norm$rotation)
mtch = match(meta.names, pca.names)
pca_sg_no_norm = pca_sg_no_norm$rotation[mtch,]
pca_sg_no_norm = as.data.frame(pca_sg_no_norm)

pca_sg_no_norm_cor<-correlate_pcs(pca_sg_no_norm, meta, npcs = 10, min.cor = 0)
pca_sg_no_norm_cor<-as.data.frame(pca_sg_no_norm_cor);

plot_cor_pca_cov(pca_sg_no_norm_cor, "PCA_SG_meta_cor_NO_NoD")

pca_sg_p_no_norm = pca.meta.regress(pca_sg_no_norm, meta, no.pcs = 6)
plot_P_pca_cov(pca_sg_p_no_norm$Ps, "pca_sg_p_No_NoDm")

write.table(file = "Predict_cov_sg_NO_NoD.txt", pca_sg_p_no_norm$Pred, quote=F, col.names = TRUE, row.names=F)

### COMPARE PCAS OF DIABETES WHEN NORMAL IS INCLUDED AND EXCLUDED.
PC_cor = calc_pc_cors(pca_sg, pca_sg_no_norm)
write.table(file = "PC_cor_sg_allSamps_vs_noNorm.txt", PC_cor, quote=F, col.names=TRUE, row.names=T)

var.select = pca.meta.regress.2var(pca_sg_no_norm, meta)
write.table(file = "Variable_select_sg_no_NoD.txt", var.select, quote=F, row.names=F, col.names=T)

###------------------------------------------------###
# 3a. Delta
Diabete_delta<-cbind(DwC_delta, DwoC_delta, NoD_delta)
data_input = Diabete_delta
colnames(data_input)<-substr(colnames(data_input),1,nchar(colnames(data_input))-5)

data_input2<-matrix(1,nrow(data_input),ncol(data_input))
for (i in 1:nrow(data_input)){
  gene<-data_input[i,]
  data_input2[i,]<-summary(lm(gene~group_all))$residual
}

rownames(data_input2)<-rownames(data_input)
colnames(data_input2)<-colnames(data_input)

data_input<-data_input2

# PCA variance explained
pca_delta <- prcomp(data_input, scale = FALSE, center = TRUE)
pca_delta = pca_delta$rotation
pca_delta = as.data.frame(pca_delta)

plot_pca(pca_delta, "PCA_delta")

# 2c. SG - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_delta)
mtch = match(meta.names, pca.names)
pca_delta = pca_delta[mtch,]
pca_delta = as.data.frame(pca_delta)

pca_delta_cor<-correlate_pcs(pca_delta, meta, npcs = 10, min.cor = 0)
pca_delta_cor<-as.data.frame(pca_delta_cor);

plot_cor_pca_cov(pca_delta_cor, "pca_delta_cov_corr_all_samples")

pca_delta_p = pca.meta.regress(pca_delta, meta)
plot_P_pca_cov(pca_delta_p$Ps, "pca_delta_p_all_samples")

write.table(file = "Predict_cov_delta.txt", pca_delta_p$Pred, quote=F, col.names = TRUE, row.names=F)

######### REPEAT ABOVE BUT EXCLUDING NO_NoD
Diabete_delta <- cbind(DwC_delta, DwoC_delta)
data_input = Diabete_delta
colnames(data_input)<-substr(colnames(data_input),1,nchar(colnames(data_input))-5)

data_input2<-matrix(1,nrow(data_input),ncol(data_input))
for (i in 1:nrow(data_input)){
  gene<-data_input[i,]
  data_input2[i,]<-summary(lm(gene~group))$residual
}

rownames(data_input2)<-rownames(data_input)
colnames(data_input2)<-colnames(data_input)

data_input<-data_input2

# 2b-2. HG - PCA variance explained
pca_delta_no_norm <- prcomp(data_input, scale = FALSE, center = TRUE)

# 1c. HG - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_delta_no_norm$rotation)
mtch = match(meta.names, pca.names)
pca_delta_no_norm = pca_delta_no_norm$rotation[mtch,]
pca_delta_no_norm = as.data.frame(pca_delta_no_norm)


pca_delta_no_norm_cor<-correlate_pcs(pca_delta_no_norm, meta, npcs = 10, min.cor = 0)
pca_delta_no_norm_cor<-as.data.frame(pca_delta_no_norm_cor);

plot_cor_pca_cov(pca_delta_no_norm_cor, "PCA_delta_meta_cor_NO_NoD")

pca_delta_p_noNorm = pca.meta.regress(pca_delta_no_norm, meta, no.pcs = 6)
plot_P_pca_cov(pca_delta_p_noNorm$Ps, "pca_delta_p_NO_NoD")

write.table(file = "Predict_cov_delta_NO_NoD.txt", pca_delta_p_noNorm$Pred, quote=F, col.names = TRUE, row.names=F)

### COMPARE PCAS OF DIABETES WHEN NORMAL IS INCLUDED AND EXCLUDED.
PC_cor = calc_pc_cors(pca_delta, pca_delta_no_norm)
write.table(file = "PC_cor_delta_allSamps_NO_NoD.txt", PC_cor, quote=F, col.names=TRUE, row.names=T)

var.select = pca.meta.regress.2var(pca_delta_no_norm, meta)
write.table(file = "Variable_select_delta.txt", var.select, quote=F, row.names=F, col.names=T)

#Test for differential expression (DE) between
#Normal vs DwoC
#Normal vs DwC
#Normal vs DwC and DwoC
#DwoC vs DwC

#These analyses will be performed under both 1) Standard glucose conditions, and 2) High glucose conditions.
#Covariates: All analyses will include PC1 (wtih group removed) as a covariate, and the last analysis will also include growth rate.