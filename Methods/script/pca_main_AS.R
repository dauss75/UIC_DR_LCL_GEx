# PCA and Covariates
# Between Group Analysis

library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)
library(pscl)
library(limma)
library(calibrate)

HOME <- "/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/"
PhenotypeDir = paste(HOME,'data/phenotype/', sep='')
RData=paste(HOME,'rdata/',sep=''); dir.create(RData, showWarnings = FALSE)
Figs = paste(HOME,'Figs/', sep=''); dir.create(Figs, showWarnings = FALSE)
Output = paste(HOME,'output/', sep=''); dir.create(Figs, showWarnings = FALSE)

setwd(HOME)
source(paste(HOME,"script/pca.R",sep=""))
source(paste(HOME,"script/Function.R",sep=""))
source(paste(HOME,"script/andrew_Func.R",sep=""))

date=Sys.Date()

# load data
dat<-local(get(load(file=paste(RData,"normalizedDataMatrix_filtered.RData",sep=""))))

probeList <- rownames(dat)  # probeID (nuID) gene annotation
if (require(lumiHumanAll.db) & require(annotate)){
  geneSymbol<- getSYMBOL(probeList, 'lumiHumanAll.db')
  geneName<- sapply(lookUp(probeList, 'lumiHumanAll.db', 'GENENAME'), function(x) x[1])
}
genes <-data.frame(ID=probeList, geneSymbol=geneSymbol, geneName=geneName, stringsAsFactors=FALSE)
#---------------------------------------------------------------------------
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

######

meta<-read.table(paste(PhenotypeDir, "meanDeltaSampleInfo3.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
group<-meta$GROUP


meta1<-read.table(paste(PhenotypeDir, "sample_group.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
group_all<-meta1$GROUP
## ------------------------------------------------

# 1a. HG - collapse data
All_HG <- cbind(DwC_hg, DwoC_hg, NoD_hg)

data_input = prep_data(All_HG)

# 1b. HG - PCA variance explained
pca_hg <- prcomp(data_input, scale = FALSE, center = TRUE)

setEPS()
postscript(file = paste(Figs,"pca_hg_variance_explained_",date,".eps",sep=""),width=11,height=7)
plot_variance_explained(pca_hg, 10) + theme_bw(base_size = 18)
dev.off()

pca_hg = pca_hg$rotation
pca_hg = as.data.frame(pca_hg)

file_name=paste(Figs,"PCA_HG",sep="")
plot_pca(pca_hg, file_name)

# 1c. HG - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_hg)
mtch = match(meta.names, pca.names)
pca_hg = pca_hg[mtch,]
pca_hg = as.data.frame(pca_hg)

pca_hg_cor<-correlate_pcs(pca_hg, meta, npcs = 10, min.cor = 0)
pca_hg_cor<-as.data.frame(pca_hg_cor);

file_name=paste(Figs,"pca_hg_cov_corr_allSamples",sep="")
plot_cor_pca_cov(pca_hg_cor, file_name)

pca_hg_p = pca.meta.regress(pca_hg, meta, no.pcs = 6)
file_name=paste(Figs,"pca_hg_pval_allSamples",sep="")
plot_P_pca_cov(pca_hg_p$Ps, file_name)

write.table(file = paste(Output,"Predict_cov_hg.txt",sep=""), pca_hg_p$Pred, quote=F, col.names = TRUE, row.names=F)

# remove group and then repeat the analysis
data_input2<-matrix(1,nrow(data_input),ncol(data_input))
for (i in 1:nrow(data_input)){
  gene<-as.numeric(data_input[i,])
  data_input2[i,]<-summary(lm(gene~group_all))$residual
}

rownames(data_input2)<-rownames(data_input)
colnames(data_input2)<-colnames(data_input)

# 1b. HG - PCA variance explained
pca_hg2 <- prcomp(data_input2, scale = FALSE, center = TRUE)

setEPS()
postscript(file = paste(Figs,"pca_hg_variance_explained_group_removed_",date,".eps",sep=""),width=11,height=7)
plot_variance_explained(pca_hg2, 10) + theme_bw(base_size = 18)
dev.off()

pca_hg2 = pca_hg2$rotation
pca_hg2 = as.data.frame(pca_hg2)

file_name=paste(Figs,"PCA_HG_group_removed",sep="")
plot_pca(pca_hg2, file_name)

# 1c. HG - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_hg2)
mtch = match(meta.names, pca.names)
pca_hg2 = pca_hg2[mtch,]
pca_hg2 = as.data.frame(pca_hg2)

pca_hg_cor2<-correlate_pcs(pca_hg2, meta, npcs = 10, min.cor = 0)
pca_hg_cor2<-as.data.frame(pca_hg_cor2);

file_name=paste(Figs,"pca_hg_cov_corr_allSamples_group_removed",sep="")
plot_cor_pca_cov(pca_hg_cor2, file_name)

pca_hg_p2 = pca.meta.regress(pca_hg2, meta, no.pcs = 6)
file_name=paste(Figs,"pca_hg_pval_allSamples_group_removed",sep="")
plot_P_pca_cov(pca_hg_p2$Ps, file_name)

write.table(file = paste(Output,"Predict_cov_hg_group_removed.txt",sep=""), pca_hg_p2$Pred, quote=F, col.names = TRUE, row.names=F)


######### REPEAT ABOVE BUT EXCLUDING NO_HG

Diabete_HG <- cbind(DwC_hg, DwoC_hg)

data_input = prep_data(Diabete_HG)

# 1b. HG - PCA variance explained
pca_hg_no_norm <- prcomp(data_input, scale = FALSE, center = TRUE)

setEPS()
postscript(file = paste(Figs,"pca_hg_variance_explained_NoD_excluded_",date,".eps",sep=""),width=11,height=7)
plot_variance_explained(pca_hg_no_norm, 10) + theme_bw(base_size = 18)
dev.off()

# 1c. HG - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_hg_no_norm$rotation)
mtch = match(meta.names, pca.names)
pca_hg_no_norm = pca_hg_no_norm$rotation[mtch,]
pca_hg_no_norm = as.data.frame(pca_hg_no_norm)

pca_hg_no_norm_cor<-correlate_pcs(pca_hg_no_norm, meta, npcs = 10, min.cor = 0)
pca_hg_no_norm_cor<-as.data.frame(pca_hg_no_norm_cor);

file_name=paste(Figs,"PCA_hg_meta_cor_no_norm",sep="")
plot_cor_pca_cov(pca_hg_no_norm_cor, file_name)


pca_hg_p_no_norm = pca.meta.regress(pca_hg_no_norm, meta, no.pcs = 6)
file_name=paste(Figs,"pca_hg_p_no_norm",sep="")
plot_P_pca_cov(pca_hg_p_no_norm$Ps, file_name)

write.table(file = paste(Output,"Predict_cov_hg_nonorm.txt",sep=""), pca_hg_p_no_norm$Pred, quote=F, col.names = TRUE, row.names=F)


data_input2<-matrix(1,nrow(data_input),ncol(data_input))
for (i in 1:nrow(data_input)){
  gene<-as.numeric(data_input[i,])
  data_input2[i,]<-summary(lm(gene~group))$residual
}

rownames(data_input2)<-rownames(data_input)
colnames(data_input2)<-colnames(data_input)

# 1b. HG - PCA variance explained
pca_hg_no_norm2 <- prcomp(data_input2, scale = FALSE, center = TRUE)

setEPS()
postscript(file = paste(Figs,"pca_pca_hg_variance_explained_NoD_excluded_group_removed_",date,".eps",sep=""),width=11,height=7)
plot_variance_explained(pca_hg_no_norm2, 10) + theme_bw(base_size = 18)
dev.off()

# 1c. HG - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_hg_no_norm2$rotation)
mtch = match(meta.names, pca.names)
pca_hg_no_norm2 = pca_hg_no_norm2$rotation[mtch,]
pca_hg_no_norm2 = as.data.frame(pca_hg_no_norm2)

pca_hg_no_norm_cor2<-correlate_pcs(pca_hg_no_norm2, meta, npcs = 10, min.cor = 0)
pca_hg_no_norm_cor2<-as.data.frame(pca_hg_no_norm_cor2);

file_name=paste(Figs,"PCA_hg_meta_cor_no_norm_group_removed",sep="")
plot_cor_pca_cov(pca_hg_no_norm_cor2, file_name)

pca_hg_p_no_norm2 = pca.meta.regress(pca_hg_no_norm2, meta, no.pcs = 6)
file_name=paste(Figs,"pca_hg_p_no_norm_group_removed",sep="")
plot_P_pca_cov(pca_hg_p_no_norm2$Ps, file_name)
# plot_P_pca_cov(pca_hg_p_no_norm$Ps, "pca_hg_p_no_norm")

write.table(file = paste(Output,"Predict_cov_hg_nonorm_group_removed.txt",sep=""), pca_hg_p_no_norm2$Pred, quote=F, col.names = TRUE, row.names=F)

### COMPARE PCAS OF DIABETES WHEN NORMAL IS INCLUDED AND EXCLUDED.
PC_cor = calc_pc_cors(pca_hg, pca_hg_no_norm)
write.table(file = paste(Output,"PC_cor_hg_allSamps_noNorm.txt",sep=""), PC_cor, quote=F, col.names=TRUE, row.names=T)

# var.select = pca.meta.regress.2var(pca_hg_no_norm, meta)
# write.table(file = paste(Output,"Variable_select_hg.txt",sep=""), var.select, quote=FALSE, col.names=TRUE, row.names=FALSE)

# group removed
PC_cor2 = calc_pc_cors(pca_hg2, pca_hg_no_norm2)
write.table(file = paste(Output,"PC_cor_hg_allSamps_noNorm_group_removed.txt",sep=""), PC_cor2, quote=F, col.names=TRUE, row.names=T)

# var.select2 = pca.meta.regress.2var(pca_hg_no_norm2, meta)
# write.table(file = paste(Output,"Variable_select_hg_group_removed.txt",sep=""), var.select2, quote=FALSE, col.names=TRUE, row.names=FALSE)

#####---------------########
# 2a. SG - collapse data
All_SG<-cbind(DwC_sg, DwoC_sg, NoD_sg)
data_input = prep_data(All_SG)

# 2b. SG - PCA variance explained
pca_sg <- prcomp(data_input, scale = FALSE, center = TRUE)

setEPS()
postscript(file = paste(Figs,"pca_sg_variance_explained_",date,".eps",sep=""),width=11,height=7)
plot_variance_explained(pca_sg, 10) + theme_bw(base_size = 18)
dev.off()

pca_sg = pca_sg$rotation
pca_sg = as.data.frame(pca_sg)

file_name=paste(Figs,"PCA_SG",sep="")
plot_pca(pca_sg, file_name)

# 2c. SG - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_sg)
mtch = match(meta.names, pca.names)
pca_sg = pca_sg[mtch,]
pca_sg = as.data.frame(pca_sg)

pca_sg_cor<-correlate_pcs(pca_sg, meta, npcs = 10, min.cor = 0)
pca_sg_cor<-as.data.frame(pca_sg_cor);

file_name=paste(Figs,"pca_sg_cov_corr_allSamples",sep="")
plot_cor_pca_cov(pca_sg_cor, file_name)

pca_sg_p = pca.meta.regress(pca_sg, meta)
file_name=paste(Figs,"pca_sg_p_allSamples",sep="")
plot_P_pca_cov(pca_sg_p$Ps, file_name)

write.table(file = paste(Output,"Predict_cov_sg.txt",sep=""), pca_sg_p$Pred, quote=F, col.names = TRUE, row.names=F)

# repeat after group removed

data_input2<-matrix(1,nrow(data_input),ncol(data_input))
for (i in 1:nrow(data_input)){
  gene<-as.numeric(data_input[i,])
  data_input2[i,]<-summary(lm(gene~group_all))$residual
}

rownames(data_input2)<-rownames(data_input)
colnames(data_input2)<-colnames(data_input)

# 2b. SG - PCA variance explained
pca_sg2 <- prcomp(data_input2, scale = FALSE, center = TRUE)

setEPS()
postscript(file = paste(Figs,"pca_sg_variance_explained_group_removed_",date,".eps",sep=""),width=11,height=7)
plot_variance_explained(pca_sg2, 10) + theme_bw(base_size = 18)
dev.off()

pca_sg2 = pca_sg2$rotation
pca_sg2 = as.data.frame(pca_sg2)

file_name=paste(Figs,"PCA_SG_group_removed",sep="")
plot_pca(pca_sg2, file_name)

# 2c. SG - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_sg2)
mtch = match(meta.names, pca.names)
pca_sg2 = pca_sg2[mtch,]
pca_sg2 = as.data.frame(pca_sg2)

pca_sg_cor2<-correlate_pcs(pca_sg2, meta, npcs = 10, min.cor = 0)
pca_sg_cor2<-as.data.frame(pca_sg_cor2);

file_name=paste(Figs,"pca_sg_cov_corr_allSamples_group_removed",sep="")
plot_cor_pca_cov(pca_sg_cor2, file_name)

pca_sg_p2 = pca.meta.regress(pca_sg2, meta)
file_name=paste(Figs,"pca_sg_p_allSamples_group_removed",sep="")
plot_P_pca_cov(pca_sg_p2$Ps, file_name)

write.table(file = paste(Output,"Predict_cov_sg_group_removed.txt",sep=""), pca_sg_p2$Pred, quote=F, col.names = TRUE, row.names=F)

######### REPEAT ABOVE BUT EXCLUDING NO_HG
Diabete_SG <- cbind(DwC_sg, DwoC_sg)
data_input = prep_data(Diabete_SG)

# 2b-2. SG - PCA variance explained
pca_sg_no_norm <- prcomp(data_input, scale = FALSE, center = TRUE)

setEPS()
postscript(file = paste(Figs,"pca_sg_variance_explained_NoD_excluded_",date,".eps",sep=""),width=11,height=7)
plot_variance_explained(pca_sg_no_norm, 10) + theme_bw(base_size = 18)
dev.off()

# 2c. HG - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_sg_no_norm$rotation)
mtch = match(meta.names, pca.names)
pca_sg_no_norm = pca_sg_no_norm$rotation[mtch,]
pca_sg_no_norm = as.data.frame(pca_sg_no_norm)

pca_sg_no_norm_cor<-correlate_pcs(pca_sg_no_norm, meta, npcs = 10, min.cor = 0)
pca_sg_no_norm_cor<-as.data.frame(pca_sg_no_norm_cor);

file_name=paste(Figs,"PCA_SG_meta_cor_no_norm",sep="")
plot_cor_pca_cov(pca_sg_no_norm_cor, file_name)

pca_sg_p_no_norm = pca.meta.regress(pca_sg_no_norm, meta, no.pcs = 6)
file_name=paste(Figs,"pca_sg_p_noNorm",sep="")
plot_P_pca_cov(pca_sg_p_no_norm$Ps, file_name)

write.table(file = paste(Output,"Predict_cov_sg_noNorm.txt",sep=""), pca_sg_p_no_norm$Pred, quote=F, col.names = TRUE, row.names=F)

data_input2<-matrix(1,nrow(data_input),ncol(data_input))
for (i in 1:nrow(data_input)){
  gene<-as.numeric(data_input[i,])
  data_input2[i,]<-summary(lm(gene~group))$residual
}

rownames(data_input2)<-rownames(data_input)
colnames(data_input2)<-colnames(data_input)


# 2b-2. HG - PCA variance explained
pca_sg_no_norm2 <- prcomp(data_input2, scale = FALSE, center = TRUE)

setEPS()
postscript(file = paste(Figs,"pca_sg_variance_explained_NoD_excluded_group_removed_",date,".eps",sep=""),width=11,height=7)
plot_variance_explained(pca_sg_no_norm2, 10) + theme_bw(base_size = 18)
dev.off()

# 2c. HG - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_sg_no_norm2$rotation)
mtch = match(meta.names, pca.names)
pca_sg_no_norm2 = pca_sg_no_norm2$rotation[mtch,]
pca_sg_no_norm2 = as.data.frame(pca_sg_no_norm2)

pca_sg_no_norm_cor2<-correlate_pcs(pca_sg_no_norm2, meta, npcs = 10, min.cor = 0)
pca_sg_no_norm_cor2<-as.data.frame(pca_sg_no_norm_cor2);

file_name=paste(Figs,"PCA_SG_meta_cor_no_norm_group_removed",sep="")
plot_cor_pca_cov(pca_sg_no_norm_cor2, file_name)

pca_sg_p_no_norm2 = pca.meta.regress(pca_sg_no_norm2, meta, no.pcs = 6)
file_name=paste(Figs,"pca_sg_p_noNorm_group_removed",sep="")
plot_P_pca_cov(pca_sg_p_no_norm2$Ps, file_name)

write.table(file = paste(Output,"Predict_cov_sg_noNorm_group_removed.txt",sep=""), pca_sg_p_no_norm2$Pred, quote=F, col.names = TRUE, row.names=F)


### COMPARE PCAS OF DIABETES WHEN NORMAL IS INCLUDED AND EXCLUDED.
PC_cor = calc_pc_cors(pca_sg, pca_sg_no_norm)
write.table(file = paste(Output,"PC_cor_sg_allSamps_noNorm.txt",sep=""), PC_cor, quote=F, col.names=TRUE, row.names=T)

var.select = pca.meta.regress.2var(pca_sg_no_norm, meta)
write.table(file = paste(Output,"Variable_select_sg.txt",sep=""), var.select, quote=F, row.names=F, col.names=T)

PC_cor2 = calc_pc_cors(pca_sg2, pca_sg_no_norm2)
write.table(file = paste(Output,"PC_cor_sg_allSamps_noNorm_group_removed.txt",sep=""), PC_cor, quote=F, col.names=TRUE, row.names=T)

var.select2 = pca.meta.regress.2var(pca_sg_no_norm2, meta)
write.table(file = paste(Output,"Variable_select_sg_group_removed.txt",sep=""), var.select2, quote=F, row.names=F, col.names=T)



###------------------------------------------------###
# 3a. Delta - collapse data
All_delta<-cbind(DwC_delta, DwoC_delta, NoD_delta)
data_input = prep_data(All_delta)

# 3b. Delta - PCA variance explained
pca_delta <- prcomp(data_input, scale = FALSE, center = TRUE)

setEPS()
postscript(file = paste(Figs,"pca_delta_variance_explained_",date,".eps",sep=""),width=11,height=7)
plot_variance_explained(pca_delta, 10) + theme_bw(base_size = 18)
dev.off()

pca_delta = pca_delta$rotation
pca_delta = as.data.frame(pca_delta)

file_name=paste(Figs,"PCA_delta",sep="")
plot_pca(pca_delta, file_name)

# 3c. Delta - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_delta)
mtch = match(meta.names, pca.names)
pca_delta = pca_delta[mtch,]
pca_delta = as.data.frame(pca_delta)

pca_delta_cor<-correlate_pcs(pca_delta, meta, npcs = 10, min.cor = 0)
pca_delta_cor<-as.data.frame(pca_delta_cor);

file_name=paste(Figs,"pca_delta_cov_corr_allSamples",sep="")
plot_cor_pca_cov(pca_delta_cor, file_name)

pca_delta_p = pca.meta.regress(pca_delta, meta)
file_name=paste(Figs,"pca_delta_p_allSamples",sep="")
plot_P_pca_cov(pca_delta_p$Ps, file_name)

write.table(file = paste(Output,"Predict_cov_delta.txt",sep=""), pca_delta_p$Pred, quote=F, col.names = TRUE, row.names=F)

data_input2<-matrix(1,nrow(data_input),ncol(data_input))
for (i in 1:nrow(data_input)){
  gene<-as.numeric(data_input[i,])
  data_input2[i,]<-summary(lm(gene~group_all))$residual
}

rownames(data_input2)<-rownames(data_input)
colnames(data_input2)<-colnames(data_input)

# data_input<-data_input2

# 3b. Delta - PCA variance explained
pca_delta2 <- prcomp(data_input2, scale = FALSE, center = TRUE)

setEPS()
postscript(file = paste(Figs,"pca_delta_variance_explained_group_removed_",date,".eps",sep=""),width=11,height=7)
plot_variance_explained(pca_delta2, 10) + theme_bw(base_size = 18)
dev.off()

pca_delta2 = pca_delta2$rotation
pca_delta2 = as.data.frame(pca_delta2)

file_name=paste(Figs,"PCA_delta_group_removed",sep="")
plot_pca(pca_delta2, file_name)

# 3c. Delta - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_delta2)
mtch = match(meta.names, pca.names)
pca_delta2 = pca_delta2[mtch,]
pca_delta2 = as.data.frame(pca_delta2)

pca_delta_cor2<-correlate_pcs(pca_delta2, meta, npcs = 10, min.cor = 0)
pca_delta_cor2<-as.data.frame(pca_delta_cor2);

file_name=paste(Figs,"pca_delta_cov_corr_allSamples_group_removed",sep="")
plot_cor_pca_cov(pca_delta_cor2, file_name)

pca_delta_p2 = pca.meta.regress(pca_delta2, meta)
file_name=paste(Figs,"pca_delta_p_allSamples_group_removed",sep="")
plot_P_pca_cov(pca_delta_p2$Ps, file_name)

write.table(file = paste(Output,"Predict_cov_delta_group_removed.txt",sep=""), pca_delta_p2$Pred, quote=F, col.names = TRUE, row.names=F)


######### REPEAT ABOVE BUT EXCLUDING NO_HG
Diabete_delta <- cbind(DwC_delta, DwoC_delta)
data_input = prep_data(Diabete_delta)

pca_delta_no_norm <- prcomp(data_input, scale = FALSE, center = TRUE)

setEPS()
postscript(file = paste(Figs,"pca_delta_variance_explained_NoD_excluded_",date,".eps",sep=""),width=11,height=7)
plot_variance_explained(pca_delta_no_norm, 10) + theme_bw(base_size = 18)
dev.off()

# 3c. Delta - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_delta_no_norm$rotation)
mtch = match(meta.names, pca.names)
pca_delta_no_norm = pca_delta_no_norm$rotation[mtch,]
pca_delta_no_norm = as.data.frame(pca_delta_no_norm)

pca_delta_no_norm_cor<-correlate_pcs(pca_delta_no_norm, meta, npcs = 10, min.cor = 0)
pca_delta_no_norm_cor<-as.data.frame(pca_delta_no_norm_cor);

file_name=paste(Figs,"PCA_delta_meta_cor_no_norm",sep="")
plot_cor_pca_cov(pca_delta_no_norm_cor, file_name)

pca_delta_p_noNorm = pca.meta.regress(pca_delta_no_norm, meta, no.pcs = 6)
file_name=paste(Figs,"pca_delta_p_noNorm",sep="")
plot_P_pca_cov(pca_delta_p_noNorm$Ps, file_name)

write.table(file = paste(Output,"Predict_cov_delta_noNorm.txt",sep=""), pca_delta_p_noNorm$Pred, quote=F, col.names = TRUE, row.names=F)

### COMPARE PCAS OF DIABETES WHEN NORMAL IS INCLUDED AND EXCLUDED.
PC_cor = calc_pc_cors(pca_delta, pca_delta_no_norm)
write.table(file = paste(Output,"PC_cor_delta_allSamps_noNorm.txt",sep=""), PC_cor, quote=F, col.names=TRUE, row.names=T)

var.select = pca.meta.regress.2var(pca_delta_no_norm, meta)
write.table(file = paste(Output,"Variable_select_delta.txt",sep=""), var.select, quote=F, row.names=F, col.names=T)

data_input2<-matrix(1,nrow(data_input),ncol(data_input))
for (i in 1:nrow(data_input)){
  gene<-as.numeric(data_input[i,])
  data_input2[i,]<-summary(lm(gene~group))$residual
}

rownames(data_input2)<-rownames(data_input)
colnames(data_input2)<-colnames(data_input)

# 3b-2. Delta - PCA variance explained
pca_delta_no_norm2 <- prcomp(data_input2, scale = FALSE, center = TRUE)

setEPS()
postscript(file = paste(Figs,"pca_delta_variance_explained_NoD_excluded_group_removed_",date,".eps",sep=""),width=11,height=7)
plot_variance_explained(pca_delta_no_norm2, 10) + theme_bw(base_size = 18)
dev.off()

# 3c. Delta - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_delta_no_norm2$rotation)
mtch = match(meta.names, pca.names)
pca_delta_no_norm2 = pca_delta_no_norm2$rotation[mtch,]
pca_delta_no_norm2 = as.data.frame(pca_delta_no_norm2)

pca_delta_no_norm_cor2<-correlate_pcs(pca_delta_no_norm2, meta, npcs = 10, min.cor = 0)
pca_delta_no_norm_cor2<-as.data.frame(pca_delta_no_norm_cor2);

file_name=paste(Figs,"PCA_delta_meta_cor_no_norm_group_removed",sep="")
plot_cor_pca_cov(pca_delta_no_norm_cor2, file_name)

pca_delta_p_noNorm2 = pca.meta.regress(pca_delta_no_norm2, meta, no.pcs = 6)
file_name=paste(Figs,"pca_delta_p_noNorm_group_removed",sep="")
plot_P_pca_cov(pca_delta_p_noNorm2$Ps, file_name)

write.table(file = paste(Output,"Predict_cov_delta_noNorm_group_removed.txt",sep=""), pca_delta_p_noNorm2$Pred, quote=F, col.names = TRUE, row.names=F)

### COMPARE PCAS OF DIABETES WHEN NORMAL IS INCLUDED AND EXCLUDED.
PC_cor2 = calc_pc_cors(pca_delta2, pca_delta_no_norm2)
write.table(file = paste(Output,"PC_cor_delta_allSamps_noNorm_group_removed.txt",sep=""), PC_cor2, quote=F, col.names=TRUE, row.names=T)

var.select = pca.meta.regress.2var(pca_delta_no_norm2, meta)
write.table(file = paste(Output,"Variable_select_delta_group_removed.txt",sep=""), var.select, quote=F, row.names=F, col.names=T)










#--------------------------------------
## DE analysis for high glucose
#--------------------------------------

All_HG <- cbind(DwC_hg, DwoC_hg, NoD_hg)

#collapse data
data_input = prep_data(All_HG)

data_input2<-matrix(1,nrow(data_input),ncol(data_input))
#regress group out
for (i in 1:nrow(data_input)){
  gene<-as.numeric(data_input[i,])
  data_input2[i,]<-summary(lm(gene~group_all))$residual
}

rownames(data_input2)<-rownames(data_input)
colnames(data_input2)<-colnames(data_input)

#--------------------------------------
# HG NoD vs DwoC
#--------------------------------------

# meta_sorted<-na.omit(meta_sorted)
meta_tmp<-read.table(paste(PhenotypeDir, "sample_group_NoD_vs_DwoC.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
group <-as.factor(meta_tmp$GROUP)
growth_rate<-meta[rownames(meta_tmp),"GROWTH_RATE"]
data_input_tmp<-data_input[,rownames(meta_tmp)]; data_input_tmp<-apply(data_input_tmp,2,as.numeric); rownames(data_input_tmp)<-rownames(data_input);

# PC1 from data after group regressed out
data_input2_tmp<-data_input2[,rownames(meta_tmp)]; data_input2_tmp<-apply(data_input2_tmp,2,as.numeric); rownames(data_input2_tmp)<-rownames(data_input2);

pca_NoD_DwoC <- prcomp(data_input2_tmp, scale = FALSE, center = TRUE)
pca_NoD_DwoC = pca_NoD_DwoC$rotation
pca_NoD_DwoC = as.data.frame(pca_NoD_DwoC)

pc<-pca_NoD_DwoC$PC1
design<-model.matrix(~pc+growth_rate+group)
fit <- eBayes(lmFit(data_input_tmp, design))

qval.cutoff=0.05; FC.cutoff=1 # FC=2
x1=topTable(fit, coef="group2", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="group2", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

FC.cutoff1.5=0.5849625
hg_NoD_DwoC_qval5_FC2<-x1_subset
hg_NoD_DwoC_qval10_FC2<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff,]
hg_NoD_DwoC_qval20_FC2<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff,]
hg_NoD_DwoC_qval5_FC1.5<-x1[x1$adj.P.Val<0.05 & abs(x1$logFC) > FC.cutoff1.5,]
hg_NoD_DwoC_qval10_FC1.5<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff1.5,]
hg_NoD_DwoC_qval20_FC1.5<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff1.5,]

setEPS()
postscript(file = paste(Figs,"hg_NoD_DwoC_PC1_GrowthRate_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in high glucose (No_DM vs. No_PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(y1, paste(Output,"hg_NoD_DwoC_PC1_GrowthRate_all_",date,".csv",sep=""),row.names = F)
write.csv(x1, paste(Output,"hg_NoD_DwoC_PC1_GrowthRate_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"hg_NoD_DwoC_PC1_GrowthRate_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#--------------------------------------
# HG NoD vs DwC
#--------------------------------------
meta_tmp<-read.table(paste(PhenotypeDir, "sample_group_NoD_vs_DwC.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
group <-as.factor(meta_tmp$GROUP)
growth_rate<-meta[rownames(meta_tmp),"GROWTH_RATE"]
data_input_tmp<-data_input[,rownames(meta_tmp)]; data_input_tmp<-apply(data_input_tmp,2,as.numeric); rownames(data_input_tmp)<-rownames(data_input);

# PC1 from data after group regressed out
data_input2_tmp<-data_input2[,rownames(meta_tmp)]; data_input2_tmp<-apply(data_input2_tmp,2,as.numeric); rownames(data_input2_tmp)<-rownames(data_input2);

pca_NoD_DwC <- prcomp(data_input2_tmp, scale = FALSE, center = TRUE)
pca_NoD_DwC = pca_NoD_DwC$rotation
pca_NoD_DwC = as.data.frame(pca_NoD_DwC)

pc<-pca_NoD_DwC$PC1
design<-model.matrix(~pc+growth_rate+group)
fit <- eBayes(lmFit(data_input_tmp, design))

qval.cutoff=0.05; FC.cutoff=1 # FC=2
x1=topTable(fit, coef="group2", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="group2", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

hg_NoD_DwC_qval5_FC2<-x1_subset
hg_NoD_DwC_qval10_FC2<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff,]
hg_NoD_DwC_qval20_FC2<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff,]
hg_NoD_DwC_qval5_FC1.5<-x1[x1$adj.P.Val<0.05 & abs(x1$logFC) > FC.cutoff1.5,]
hg_NoD_DwC_qval10_FC1.5<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff1.5,]
hg_NoD_DwC_qval20_FC1.5<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff1.5,]

setEPS()
postscript(file = paste(Figs,"hg_NoD_DwC_PC1_GrowthRate_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in high glucose (No_DM vs. PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(y1, paste(Output,"hg_NoD_DwC_PC1_GrowthRate_all_",date,".csv",sep=""),row.names = F)
write.csv(x1, paste(Output,"hg_NoD_DwC_PC1_GrowthRate_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"hg_NoD_DwC_PC1_GrowthRate_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)


#--------------------------------------
# HG NoD vs Diabetes
#--------------------------------------
meta_tmp<-read.table(paste(PhenotypeDir, "sample_group_NoD_vs_Diabetes.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
group <-as.factor(meta_tmp$GROUP)
growth_rate<-meta[rownames(meta_tmp),"GROWTH_RATE"]
data_input_tmp<-data_input[,rownames(meta_tmp)]; data_input_tmp<-apply(data_input_tmp,2,as.numeric); rownames(data_input_tmp)<-rownames(data_input);

# PC1 from data after group regressed out
data_input2_tmp<-data_input2[,rownames(meta_tmp)]; data_input2_tmp<-apply(data_input2_tmp,2,as.numeric); rownames(data_input2_tmp)<-rownames(data_input2);

pca_NoD_Diabetes <- prcomp(data_input2_tmp, scale = FALSE, center = TRUE)
pca_NoD_Diabetes = pca_NoD_Diabetes$rotation
pca_NoD_Diabetes = as.data.frame(pca_NoD_Diabetes)

pc<-pca_NoD_Diabetes$PC1
design<-model.matrix(~pc+growth_rate+group)
fit <- eBayes(lmFit(data_input_tmp, design))

qval.cutoff=0.05; FC.cutoff=1 # FC=2
x1=topTable(fit, coef="group2", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="group2", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

hg_NoD_Diabetes_qval5_FC2<-x1_subset
hg_NoD_Diabetes_qval10_FC2<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff,]
hg_NoD_Diabetes_qval20_FC2<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff,]
hg_NoD_Diabetes_qval5_FC1.5<-x1[x1$adj.P.Val<0.05 & abs(x1$logFC) > FC.cutoff1.5,]
hg_NoD_Diabetes_qval10_FC1.5<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff1.5,]
hg_NoD_Diabetes_qval20_FC1.5<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff1.5,]

setEPS()
postscript(file = paste(Figs,"hg_NoD_Diabetes_PC1_GrowthRate_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in high glucose (No_DM vs. No_PDR + PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(y1, paste(Output,"hg_NoD_Diabetes_PC1_GrowthRate_all_",date,".csv",sep=""),row.names = F)
write.csv(x1, paste(Output,"hg_NoD_Diabetes_PC1_GrowthRate_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"hg_NoD_Diabetes_PC1_GrowthRate_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#--------------------------------------
# HG DwC vs DwoC
#--------------------------------------
meta_tmp<-read.table(paste(PhenotypeDir, "sample_group_DwC_vs_DwoC.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
group <-as.factor(meta_tmp$GROUP)
growth_rate<-meta[rownames(meta_tmp),"GROWTH_RATE"]
data_input_tmp<-data_input[,rownames(meta_tmp)]; data_input_tmp<-apply(data_input_tmp,2,as.numeric); rownames(data_input_tmp)<-rownames(data_input);

# PC1 from data after group regressed out
data_input2_tmp<-data_input2[,rownames(meta_tmp)]; data_input2_tmp<-apply(data_input2_tmp,2,as.numeric); rownames(data_input2_tmp)<-rownames(data_input2);

pca_DwC_DwoC <- prcomp(data_input2_tmp, scale = FALSE, center = TRUE)
pca_DwC_DwoC = pca_DwC_DwoC$rotation
pca_DwC_DwoC = as.data.frame(pca_DwC_DwoC)

pc<-pca_DwC_DwoC$PC1
design<-model.matrix(~pc+growth_rate+group)
fit <- eBayes(lmFit(data_input_tmp, design))

pval.cutoff=0.05; FC.cutoff=1 # FC=2
x1=topTable(fit, coef="group1", n=nrow(genes), adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="group1", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$P.Value<pval.cutoff & abs(x1$logFC) > FC.cutoff,]

hg_DwC_DwoC_qval5_FC2<-x1_subset
hg_DwC_DwoC_qval10_FC2<-x1[x1$P.Value<0.1 & abs(x1$logFC) > FC.cutoff,]
hg_DwC_DwoC_qval20_FC2<-x1[x1$P.Value<0.2 & abs(x1$logFC) > FC.cutoff,]
hg_DwC_DwoC_qval5_FC1.5<-x1[x1$P.Value<0.05 & abs(x1$logFC) > FC.cutoff1.5,]
hg_DwC_DwoC_qval10_FC1.5<-x1[x1$P.Value<0.1 & abs(x1$logFC) > FC.cutoff1.5,]
hg_DwC_DwoC_qval20_FC1.5<-x1[x1$P.Value<0.2 & abs(x1$logFC) > FC.cutoff1.5,]

setEPS()
postscript(file = paste(Figs,"hg_DwC_DwoC_PC1_GrowthRate_pvalue_",pval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(P.Value), pch=20, col="gray", main=paste("Group difference in high glucose (PDR vs. nPDR)",sep=""), xlab="log2 fold change", ylab="-log10(P.Value)", xlim=c(-3,3)))
with(subset(y1, P.Value<pval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(y1, P.Value<pval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()

write.csv(y1, paste(Output,"hg_DwC_DwoC_PC1_GrowthRate_all_",date,".csv",sep=""),row.names = F)
write.csv(x1, paste(Output,"hg_DwC_DwoC_PC1_GrowthRate_pvalue_",pval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"hg_DwC_DwoC_PC1_GrowthRate_pvalue_",pval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

## find common genes
library(VennDiagram)

# q-val <0.05 & FC > 2
hg_a1<-length(hg_NoD_DwoC_qval5_FC2$ID)
hg_a2<-length(hg_NoD_DwC_qval5_FC2$ID)
hg_a3<-length(hg_NoD_Diabetes_qval5_FC2$ID)
hg_a4<-length(hg_DwC_DwoC_qval5_FC2$ID)

hg_n12=length(intersect(hg_NoD_DwoC_qval5_FC2$ID, hg_NoD_DwC_qval5_FC2$ID))
hg_n13=length(intersect(hg_NoD_DwoC_qval5_FC2$ID, hg_NoD_Diabetes_qval5_FC2$ID))
hg_n14=length(intersect(hg_NoD_DwoC_qval5_FC2$ID, hg_DwC_DwoC_qval5_FC2$ID))
hg_n23=length(intersect(hg_NoD_DwC_qval5_FC2$ID, hg_NoD_Diabetes_qval5_FC2$ID))
hg_n24=length(intersect(hg_NoD_DwC_qval5_FC2$ID, hg_DwC_DwoC_qval5_FC2$ID))
hg_n34=length(intersect(hg_NoD_Diabetes_qval5_FC2$ID, hg_DwC_DwoC_qval5_FC2$ID))
hg_n123=length(intersect(intersect(hg_NoD_DwoC_qval5_FC2$ID, hg_NoD_DwC_qval5_FC2$ID), hg_NoD_Diabetes_qval5_FC2$ID))
hg_n124=length(intersect(intersect(hg_NoD_DwoC_qval5_FC2$ID, hg_NoD_DwC_qval5_FC2$ID), hg_DwC_DwoC_qval5_FC2$ID))
hg_n134=length(intersect(intersect(hg_NoD_DwoC_qval5_FC2$ID, hg_NoD_Diabetes_qval5_FC2$ID), hg_DwC_DwoC_qval5_FC2$ID))
hg_n234=length(intersect(intersect(hg_NoD_DwC_qval5_FC2$ID, hg_NoD_Diabetes_qval5_FC2$ID), hg_DwC_DwoC_qval5_FC2$ID))
hg_n1234=length(intersect(intersect(intersect(hg_NoD_DwoC_qval5_FC2$ID, hg_NoD_DwC_qval5_FC2$ID), hg_NoD_Diabetes_qval5_FC2$ID), hg_DwC_DwoC_qval5_FC2$ID))

plot.new()
pdf(paste(Figs,"HG_three_Group_Comparison_qval5_FC2.pdf",sep=""),width=12)
draw.triple.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, hg_n12, hg_n13, hg_n23, hg_n123,
                 category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"HG_four_Group_Comparison_qval5_FC2.pdf",sep=""),width=12)
draw.quad.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, area4=hg_a4, hg_n12, hg_n13, hg_n14, hg_n23, 
               hg_n24, hg_n34, hg_n123, hg_n124, hg_n134,
               hg_n234, hg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwC_DwoC"),
               lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()


# q-val <0.1 & FC > 2
hg_a1<-length(hg_NoD_DwoC_qval10_FC2$ID)
hg_a2<-length(hg_NoD_DwC_qval10_FC2$ID)
hg_a3<-length(hg_NoD_Diabetes_qval10_FC2$ID)
hg_a4<-length(hg_DwC_DwoC_qval10_FC2$ID)

hg_n12=length(intersect(hg_NoD_DwoC_qval10_FC2$ID, hg_NoD_DwC_qval10_FC2$ID))
hg_n13=length(intersect(hg_NoD_DwoC_qval10_FC2$ID, hg_NoD_Diabetes_qval10_FC2$ID))
hg_n14=length(intersect(hg_NoD_DwoC_qval10_FC2$ID, hg_DwC_DwoC_qval10_FC2$ID))
hg_n23=length(intersect(hg_NoD_DwC_qval10_FC2$ID, hg_NoD_Diabetes_qval10_FC2$ID))
hg_n24=length(intersect(hg_NoD_DwC_qval10_FC2$ID, hg_DwC_DwoC_qval10_FC2$ID))
hg_n34=length(intersect(hg_NoD_Diabetes_qval10_FC2$ID, hg_DwC_DwoC_qval10_FC2$ID))
hg_n123=length(intersect(intersect(hg_NoD_DwoC_qval10_FC2$ID, hg_NoD_DwC_qval10_FC2$ID), hg_NoD_Diabetes_qval10_FC2$ID))
hg_n124=length(intersect(intersect(hg_NoD_DwoC_qval10_FC2$ID, hg_NoD_DwC_qval10_FC2$ID), hg_DwC_DwoC_qval10_FC2$ID))
hg_n134=length(intersect(intersect(hg_NoD_DwoC_qval10_FC2$ID, hg_NoD_Diabetes_qval10_FC2$ID),hg_DwC_DwoC_qval10_FC2$ID))
hg_n234=length(intersect(intersect(hg_NoD_DwC_qval10_FC2$ID, hg_NoD_Diabetes_qval10_FC2$ID), hg_DwC_DwoC_qval10_FC2$ID))
hg_n1234=length(intersect(intersect(intersect(hg_NoD_DwoC_qval10_FC2$ID, hg_NoD_DwC_qval10_FC2$ID), hg_NoD_Diabetes_qval10_FC2$ID), hg_DwC_DwoC_qval10_FC2$ID))

plot.new()
pdf(paste(Figs,"HG_three_Group_Comparison_qval10_FC2.pdf",sep=""),width=12)
draw.triple.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, hg_n12, hg_n13, hg_n23, hg_n123,
                 category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"HG_four_Group_Comparison_qval10_FC2.pdf",sep=""),width=12)
draw.quad.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, area4=hg_a4, hg_n12, hg_n13, hg_n14, hg_n23, hg_n24, hg_n34, hg_n123, hg_n124, hg_n134,
               hg_n234, hg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwC_DwoC"),
               lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()

# q-val <0.2 & FC > 2
hg_a1<-length(hg_NoD_DwoC_qval20_FC2$ID)
hg_a2<-length(hg_NoD_DwC_qval20_FC2$ID)
hg_a3<-length(hg_NoD_Diabetes_qval20_FC2$ID)
hg_a4<-length(hg_DwC_DwoC_qval20_FC2$ID)

hg_n12=length(intersect(hg_NoD_DwoC_qval20_FC2$ID, hg_NoD_DwC_qval20_FC2$ID))
hg_n13=length(intersect(hg_NoD_DwoC_qval20_FC2$ID, hg_NoD_Diabetes_qval20_FC2$ID))
hg_n14=length(intersect(hg_NoD_DwoC_qval20_FC2$ID, hg_DwC_DwoC_qval20_FC2$ID))
hg_n23=length(intersect(hg_NoD_DwC_qval20_FC2$ID, hg_NoD_Diabetes_qval20_FC2$ID))
hg_n24=length(intersect(hg_NoD_DwC_qval20_FC2$ID, hg_DwC_DwoC_qval20_FC2$ID))
hg_n34=length(intersect(hg_NoD_Diabetes_qval20_FC2$ID, hg_DwC_DwoC_qval20_FC2$ID))
hg_n123=length(intersect(intersect(hg_NoD_DwoC_qval20_FC2$ID, hg_NoD_DwC_qval20_FC2$ID), hg_NoD_Diabetes_qval20_FC2$ID))
hg_n124=length(intersect(intersect(hg_NoD_DwoC_qval20_FC2$ID, hg_NoD_DwC_qval20_FC2$ID), hg_DwC_DwoC_qval20_FC2$ID))
hg_n134=length(intersect(intersect(hg_NoD_DwoC_qval20_FC2$ID, hg_NoD_Diabetes_qval20_FC2$ID), hg_DwC_DwoC_qval20_FC2$ID))
hg_n234=length(intersect(intersect(hg_NoD_DwC_qval20_FC2$ID, hg_NoD_Diabetes_qval20_FC2$ID), hg_DwC_DwoC_qval20_FC2$ID))
hg_n1234=length(intersect(intersect(intersect(hg_NoD_DwoC_qval20_FC2$ID, hg_NoD_DwC_qval20_FC2$ID), hg_NoD_Diabetes_qval20_FC2$ID), hg_DwC_DwoC_qval20_FC2$ID))

plot.new()
pdf(paste(Figs,"HG_three_Group_Comparison_qval20_FC2.pdf",sep=""),width=12)
draw.triple.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, hg_n12, hg_n13, hg_n23, hg_n123,
                 category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"HG_four_Group_Comparison_qval20_FC2.pdf",sep=""),width=12)
draw.quad.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, area4=hg_a4, hg_n12, hg_n13, hg_n14, hg_n23, hg_n24, hg_n34,hg_n123, hg_n124, hg_n134,
               hg_n234, hg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwC_DwoC"),
               lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()

# q-val <0.05 & FC > 1.5
hg_a1<-length(hg_NoD_DwoC_qval5_FC1.5$ID)
hg_a2<-length(hg_NoD_DwC_qval5_FC1.5$ID)
hg_a3<-length(hg_NoD_Diabetes_qval5_FC1.5$ID)
hg_a4<-length(hg_DwC_DwoC_qval5_FC1.5$ID)

hg_n12=length(intersect(hg_NoD_DwoC_qval5_FC1.5$ID, hg_NoD_DwC_qval5_FC1.5$ID))
hg_n13=length(intersect(hg_NoD_DwoC_qval5_FC1.5$ID, hg_NoD_Diabetes_qval5_FC1.5$ID))
hg_n14=length(intersect(hg_NoD_DwoC_qval5_FC1.5$ID, hg_DwC_DwoC_qval5_FC1.5$ID))
hg_n23=length(intersect(hg_NoD_DwC_qval5_FC1.5$ID, hg_NoD_Diabetes_qval5_FC1.5$ID))
hg_n24=length(intersect(hg_NoD_DwC_qval5_FC1.5$ID, hg_DwC_DwoC_qval5_FC1.5$ID))
hg_n34=length(intersect(hg_NoD_Diabetes_qval5_FC1.5$ID, hg_DwC_DwoC_qval5_FC1.5$ID))
hg_n123=length(intersect(intersect(hg_NoD_DwoC_qval5_FC1.5$ID, hg_NoD_DwC_qval5_FC1.5$ID), hg_NoD_Diabetes_qval5_FC1.5$ID))
hg_n124=length(intersect(intersect(hg_NoD_DwoC_qval5_FC1.5$ID, hg_NoD_DwC_qval5_FC1.5$ID), hg_DwC_DwoC_qval5_FC1.5$ID))
hg_n134=length(intersect(intersect(hg_NoD_DwoC_qval5_FC1.5$ID, hg_NoD_Diabetes_qval5_FC1.5$ID), hg_DwC_DwoC_qval5_FC1.5$ID))
hg_n234=length(intersect(intersect(hg_NoD_DwC_qval5_FC1.5$ID, hg_NoD_Diabetes_qval5_FC1.5$ID), hg_DwC_DwoC_qval5_FC1.5$ID))
hg_n1234=length(intersect(intersect(intersect(hg_NoD_DwoC_qval5_FC1.5$ID, hg_NoD_DwC_qval5_FC1.5$ID), hg_NoD_Diabetes_qval5_FC1.5$ID), hg_DwC_DwoC_qval5_FC1.5$ID))

plot.new()
pdf(paste(Figs,"HG_three_Group_Comparison_qval5_FC1.5.pdf",sep=""),width=12)
draw.triple.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, hg_n12, hg_n13, hg_n23, hg_n123,
                 category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"HG_four_Group_Comparison_qval5_FC1.5.pdf",sep=""),width=12)
draw.quad.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, area4=hg_a4, hg_n12, hg_n13, hg_n14, hg_n23, hg_n24, hg_n34, hg_n123, hg_n124, hg_n134,
               hg_n234, hg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwC_DwoC"),
               lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()


# q-val <0.1 & FC > 1.5
hg_a1<-length(hg_NoD_DwoC_qval10_FC1.5$ID)
hg_a2<-length(hg_NoD_DwC_qval10_FC1.5$ID)
hg_a3<-length(hg_NoD_Diabetes_qval10_FC1.5$ID)
hg_a4<-length(hg_DwC_DwoC_qval10_FC1.5$ID)

hg_n12=length(intersect(hg_NoD_DwoC_qval10_FC1.5$ID, hg_NoD_DwC_qval10_FC1.5$ID))
hg_n13=length(intersect(hg_NoD_DwoC_qval10_FC1.5$ID, hg_NoD_Diabetes_qval10_FC1.5$ID))
hg_n14=length(intersect(hg_NoD_DwoC_qval10_FC1.5$ID, hg_DwC_DwoC_qval10_FC1.5$ID))
hg_n23=length(intersect(hg_NoD_DwC_qval10_FC1.5$ID, hg_NoD_Diabetes_qval10_FC1.5$ID))
hg_n24=length(intersect(hg_NoD_DwC_qval10_FC1.5$ID, hg_DwC_DwoC_qval10_FC1.5$ID))
hg_n34=length(intersect(hg_NoD_Diabetes_qval10_FC1.5$ID, hg_DwC_DwoC_qval10_FC1.5$ID))
hg_n123=length(intersect(intersect(hg_NoD_DwoC_qval10_FC1.5$ID, hg_NoD_DwC_qval10_FC1.5$ID), hg_NoD_Diabetes_qval10_FC1.5$ID))
hg_n124=length(intersect(intersect(hg_NoD_DwoC_qval10_FC1.5$ID, hg_NoD_DwC_qval10_FC1.5$ID), hg_DwC_DwoC_qval10_FC1.5$ID))
hg_n134=length(intersect(intersect(hg_NoD_DwoC_qval10_FC1.5$ID, hg_NoD_Diabetes_qval10_FC1.5$ID), hg_DwC_DwoC_qval10_FC1.5$ID))
hg_n234=length(intersect(intersect(hg_NoD_DwC_qval10_FC1.5$ID, hg_NoD_Diabetes_qval10_FC1.5$ID), hg_DwC_DwoC_qval10_FC1.5$ID))
hg_n1234=length(intersect(intersect(intersect(hg_NoD_DwoC_qval10_FC1.5$ID, hg_NoD_DwC_qval10_FC1.5$ID), hg_NoD_Diabetes_qval10_FC1.5$ID), hg_DwC_DwoC_qval10_FC1.5$ID))

plot.new()
pdf(paste(Figs,"HG_three_Group_Comparison_qval10_FC1.5.pdf",sep=""),width=12)
draw.triple.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, hg_n12, hg_n13, hg_n23, hg_n123,
                 category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"HG_four_Group_Comparison_qval10_FC1.5.pdf",sep=""),width=12)
draw.quad.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, area4=hg_a4, hg_n12, hg_n13, hg_n14, hg_n23, hg_n24, hg_n34, hg_n123, hg_n124, hg_n134,
               hg_n234, hg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwC_DwoC"),
               lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()

# q-val <0.2 & FC > 1.5
hg_a1<-length(hg_NoD_DwoC_qval20_FC1.5$ID)
hg_a2<-length(hg_NoD_DwC_qval20_FC1.5$ID)
hg_a3<-length(hg_NoD_Diabetes_qval20_FC1.5$ID)
hg_a4<-length(hg_DwC_DwoC_qval20_FC1.5$ID)

hg_n12=length(intersect(hg_NoD_DwoC_qval20_FC1.5$ID, hg_NoD_DwC_qval20_FC1.5$ID))
hg_n13=length(intersect(hg_NoD_DwoC_qval20_FC1.5$ID, hg_NoD_Diabetes_qval20_FC1.5$ID))
hg_n14=length(intersect(hg_NoD_DwoC_qval20_FC1.5$ID, hg_DwC_DwoC_qval20_FC1.5$ID))
hg_n23=length(intersect(hg_NoD_DwC_qval20_FC1.5$ID, hg_NoD_Diabetes_qval20_FC1.5$ID))
hg_n24=length(intersect(hg_NoD_DwC_qval20_FC1.5$ID, hg_DwC_DwoC_qval20_FC1.5$ID))
hg_n34=length(intersect(hg_NoD_Diabetes_qval20_FC1.5$ID, hg_DwC_DwoC_qval20_FC1.5$ID))
hg_n123=length(intersect(intersect(hg_NoD_DwoC_qval20_FC1.5$ID, hg_NoD_DwC_qval20_FC1.5$ID), hg_NoD_Diabetes_qval20_FC1.5$ID))
hg_n124=length(intersect(intersect(hg_NoD_DwoC_qval20_FC1.5$ID, hg_NoD_DwC_qval20_FC1.5$ID), hg_DwC_DwoC_qval20_FC1.5$ID))
hg_n134=length(intersect(intersect(hg_NoD_DwoC_qval20_FC1.5$ID, hg_NoD_Diabetes_qval20_FC1.5$ID), hg_DwC_DwoC_qval20_FC1.5$ID))
hg_n234=length(intersect(intersect(hg_NoD_DwC_qval20_FC1.5$ID, hg_NoD_Diabetes_qval20_FC1.5$ID), hg_DwC_DwoC_qval20_FC1.5$ID))
hg_n1234=length(intersect(intersect(intersect(hg_NoD_DwoC_qval20_FC1.5$ID, hg_NoD_DwC_qval20_FC1.5$ID), hg_NoD_Diabetes_qval20_FC1.5$ID), hg_DwC_DwoC_qval20_FC1.5$ID))

plot.new()
pdf(paste(Figs,"HG_three_Group_Comparison_qval20_FC1.5.pdf",sep=""),width=12)
draw.triple.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, hg_n12, hg_n13, hg_n23, hg_n123,
                 category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"HG_four_Group_Comparison_qval20_FC1.5.pdf",sep=""),width=12)
draw.quad.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, area4=hg_a4, hg_n12, hg_n13, hg_n14, hg_n23, hg_n24, hg_n34, hg_n123, hg_n124, hg_n134,
               hg_n234, hg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwC_DwoC"),
               lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()

#--------------------------------------
## DE analysis of standard glucose
#--------------------------------------
Diabete_SG<-cbind(DwC_sg, DwoC_sg, NoD_sg)
data_input = prep_data(Diabete_SG)

data_input2<-matrix(1,nrow(data_input),ncol(data_input))
for (i in 1:nrow(data_input)){
  gene<-as.numeric(data_input[i,])
  data_input2[i,]<-summary(lm(gene~group_all))$residual
}

rownames(data_input2)<-rownames(data_input)
colnames(data_input2)<-colnames(data_input)

#--------------------------------------
# SG NoD vs DwoC
#--------------------------------------
meta_tmp<-read.table(paste(PhenotypeDir, "sample_group_NoD_vs_DwoC.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
group <-as.factor(meta_tmp$GROUP)
growth_rate<-meta[rownames(meta_tmp),"GROWTH_RATE"]
data_input_tmp<-data_input[,rownames(meta_tmp)]; data_input_tmp<-apply(data_input_tmp,2,as.numeric); rownames(data_input_tmp)<-rownames(data_input);

# PC1 from data after group regressed out
data_input2_tmp<-data_input2[,rownames(meta_tmp)]; data_input2_tmp<-apply(data_input2_tmp,2,as.numeric); rownames(data_input2_tmp)<-rownames(data_input2);

pca_NoD_DwoC <- prcomp(data_input2_tmp, scale = FALSE, center = TRUE)
pca_NoD_DwoC = pca_NoD_DwoC$rotation
pca_NoD_DwoC = as.data.frame(pca_NoD_DwoC)

pc<-pca_NoD_DwoC$PC1
design<-model.matrix(~pc+growth_rate+group)
fit <- eBayes(lmFit(data_input_tmp, design))

qval.cutoff=0.05; FC.cutoff=1 # FC=2
x1=topTable(fit, coef="group2", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="group2", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

FC.cutoff1.5=0.5849625
sg_NoD_DwoC_qval5_FC2<-x1_subset
sg_NoD_DwoC_qval10_FC2<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff,]
sg_NoD_DwoC_qval20_FC2<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff,]
sg_NoD_DwoC_qval5_FC1.5<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff1.5,]
sg_NoD_DwoC_qval10_FC1.5<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff1.5,]
sg_NoD_DwoC_qval20_FC1.5<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff1.5,]

setEPS()
postscript(file = paste(Figs,"sg_NoD_DwoC_PC1_GrowthRate_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in standard glucose (No_DM vs. No_PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(y1, paste(Output,"sg_NoD_DwoC_PC1_GrowthRate_all_",date,".csv",sep=""),row.names = F)
write.csv(x1, paste(Output,"sg_NoD_DwoC_PC1_GrowthRate_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"sg_NoD_DwoC_PC1_GrowthRate_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#--------------------------------------
# SG NoD vs DwC
#--------------------------------------
meta_tmp<-read.table(paste(PhenotypeDir, "sample_group_NoD_vs_DwC.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
group <-as.factor(meta_tmp$GROUP)
growth_rate<-meta[rownames(meta_tmp),"GROWTH_RATE"]
data_input_tmp<-data_input[,rownames(meta_tmp)]; data_input_tmp<-apply(data_input_tmp,2,as.numeric); rownames(data_input_tmp)<-rownames(data_input);

# PC1 from data after group regressed out
data_input2_tmp<-data_input2[,rownames(meta_tmp)]; data_input2_tmp<-apply(data_input2_tmp,2,as.numeric); rownames(data_input2_tmp)<-rownames(data_input2);

pca_NoD_DwC <- prcomp(data_input2_tmp, scale = FALSE, center = TRUE)
pca_NoD_DwC = pca_NoD_DwC$rotation
pca_NoD_DwC = as.data.frame(pca_NoD_DwC)

pc<-pca_NoD_DwC$PC1
design<-model.matrix(~pc+growth_rate+group)
fit <- eBayes(lmFit(data_input_tmp, design))

qval.cutoff=0.05; FC.cutoff=1 # FC=2
x1=topTable(fit, coef="group2", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="group2", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

sg_NoD_DwC_qval5_FC2<-x1_subset
sg_NoD_DwC_qval10_FC2<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff,]
sg_NoD_DwC_qval20_FC2<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff,]
sg_NoD_DwC_qval5_FC1.5<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff1.5,]
sg_NoD_DwC_qval10_FC1.5<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff1.5,]
sg_NoD_DwC_qval20_FC1.5<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff1.5,]

setEPS()
postscript(file = paste(Figs,"sg_NoD_DwC_PC1_GrowthRate_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in standard glucose (No_DM vs. PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(y1, paste(Output,"sg_NoD_DwC_PC1_GrowthRate_all_",date,".csv",sep=""),row.names = F)
write.csv(x1, paste(Output,"sg_NoD_DwC_PC1_GrowthRate_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"sg_NoD_DwC_PC1_GrowthRate_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)


#--------------------------------------
# SG NoD vs Diabetes
#--------------------------------------
meta_tmp<-read.table(paste(PhenotypeDir, "sample_group_NoD_vs_Diabetes.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
group <-as.factor(meta_tmp$GROUP)
growth_rate<-meta[rownames(meta_tmp),"GROWTH_RATE"]
data_input_tmp<-data_input[,rownames(meta_tmp)]; data_input_tmp<-apply(data_input_tmp,2,as.numeric); rownames(data_input_tmp)<-rownames(data_input);

# PC1 from data after group regressed out
data_input2_tmp<-data_input2[,rownames(meta_tmp)]; data_input2_tmp<-apply(data_input2_tmp,2,as.numeric); rownames(data_input2_tmp)<-rownames(data_input2);

pca_NoD_Diabetes <- prcomp(data_input2_tmp, scale = FALSE, center = TRUE)
pca_NoD_Diabetes = pca_NoD_Diabetes$rotation
pca_NoD_Diabetes = as.data.frame(pca_NoD_Diabetes)

pc<-pca_NoD_Diabetes$PC1
design<-model.matrix(~pc+growth_rate+group)
fit <- eBayes(lmFit(data_input_tmp, design))

qval.cutoff=0.05; FC.cutoff=1 # FC=2
x1=topTable(fit, coef="group2", n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="group2", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]

sg_NoD_Diabetes_qval5_FC2<-x1_subset
sg_NoD_Diabetes_qval10_FC2<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff,]
sg_NoD_Diabetes_qval20_FC2<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff,]
sg_NoD_Diabetes_qval5_FC1.5<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff1.5,]
sg_NoD_Diabetes_qval10_FC1.5<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff1.5,]
sg_NoD_Diabetes_qval20_FC1.5<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff1.5,]

setEPS()
postscript(file = paste(Figs,"sg_NoD_Diabetes_PC1_GrowthRate_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Group difference in standard glucose (No_DM vs. No_PDR + PDR)",sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-3,3)))
with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()

write.csv(y1, paste(Output,"sg_NoD_Diabetes_PC1_GrowthRate_all_",date,".csv",sep=""),row.names = F)
write.csv(x1, paste(Output,"sg_NoD_Diabetes_PC1_GrowthRate_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"sg_NoD_Diabetes_PC1_GrowthRate_qval_",qval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#--------------------------------------
# SG DwC vs DwoC
#--------------------------------------
meta_tmp<-read.table(paste(PhenotypeDir, "sample_group_DwC_vs_DwoC.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
group <-as.factor(meta_tmp$GROUP)
growth_rate<-meta[rownames(meta_tmp),"GROWTH_RATE"]
data_input_tmp<-data_input[,rownames(meta_tmp)]; data_input_tmp<-apply(data_input_tmp,2,as.numeric); rownames(data_input_tmp)<-rownames(data_input);

# PC1 from data after group regressed out
data_input2_tmp<-data_input2[,rownames(meta_tmp)]; data_input2_tmp<-apply(data_input2_tmp,2,as.numeric); rownames(data_input2_tmp)<-rownames(data_input2);

pca_DwC_DwoC <- prcomp(data_input2_tmp, scale = FALSE, center = TRUE)
pca_DwC_DwoC = pca_DwC_DwoC$rotation
pca_DwC_DwoC = as.data.frame(pca_DwC_DwoC)

pc<-pca_DwC_DwoC$PC1
design<-model.matrix(~pc+growth_rate+group)
fit <- eBayes(lmFit(data_input_tmp, design))

pval.cutoff=0.05; FC.cutoff=1 # FC=2
x1=topTable(fit, coef="group1", n=nrow(genes), adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="group1", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$P.Value<pval.cutoff & abs(x1$logFC) > FC.cutoff,]

sg_DwC_DwoC_qval5_FC2<-x1_subset
sg_DwC_DwoC_qval10_FC2<-x1[x1$P.Value<0.1 & abs(x1$logFC) > FC.cutoff,]
sg_DwC_DwoC_qval20_FC2<-x1[x1$P.Value<0.2 & abs(x1$logFC) > FC.cutoff,]
sg_DwC_DwoC_qval5_FC1.5<-x1[x1$P.Value<qval.cutoff & abs(x1$logFC) > FC.cutoff1.5,]
sg_DwC_DwoC_qval10_FC1.5<-x1[x1$P.Value<0.1 & abs(x1$logFC) > FC.cutoff1.5,]
sg_DwC_DwoC_qval20_FC1.5<-x1[x1$P.Value<0.2 & abs(x1$logFC) > FC.cutoff1.5,]

# q-val <0.05 & FC > 2
sg_a1<-length(NoD_DwoC_qval5_FC2$ID)
sg_a2<-length(sg_NoD_DwC_qval5_FC2$ID)
sg_a3<-length(sg_NoD_Diabetes_qval5_FC2$ID)
sg_a4<-length(sg_DwC_DwoC_qval5_FC2$ID)

sg_n12=length(intersect(sg_NoD_DwoC_qval5_FC2$ID, sg_NoD_DwC_qval5_FC2$ID))
sg_n13=length(intersect(sg_NoD_DwoC_qval5_FC2$ID, sg_NoD_Diabetes_qval5_FC2$ID))
sg_n14=length(intersect(sg_NoD_DwoC_qval5_FC2$ID, sg_DwC_DwoC_qval5_FC2$ID))
sg_n23=length(intersect(sg_NoD_DwC_qval5_FC2$ID, sg_NoD_Diabetes_qval5_FC2$ID))
sg_n24=length(intersect(sg_NoD_DwC_qval5_FC2$ID, sg_DwC_DwoC_qval5_FC2$ID))
sg_n34=length(intersect(sg_NoD_Diabetes_qval5_FC2$ID, sg_DwC_DwoC_qval5_FC2$ID))
sg_n123=length(intersect(intersect(sg_NoD_DwoC_qval5_FC2$ID, sg_NoD_DwC_qval5_FC2$ID), sg_NoD_Diabetes_qval5_FC2$ID))
sg_n124=length(intersect(intersect(sg_NoD_DwoC_qval5_FC2$ID, sg_NoD_DwC_qval5_FC2$ID), sg_DwC_DwoC_qval5_FC2$ID))
sg_n134=length(intersect(intersect(sg_NoD_DwoC_qval5_FC2$ID, sg_NoD_Diabetes_qval5_FC2$ID), sg_DwC_DwoC_qval5_FC2$ID))
sg_n234=length(intersect(intersect(sg_NoD_DwC_qval5_FC2$ID, sg_NoD_Diabetes_qval5_FC2$ID), sg_DwC_DwoC_qval5_FC2$ID))
sg_n1234=length(intersect(intersect(intersect(sg_NoD_DwoC_qval5_FC2$ID, sg_NoD_DwC_qval5_FC2$ID), sg_NoD_Diabetes_qval5_FC2$ID), sg_DwC_DwoC_qval5_FC2$ID))

plot.new()
pdf(paste(Figs,"SG_three_Group_Comparison_qval5_FC2.pdf",sep=""),width=12)
draw.triple.venn(area1=a1, area2=a2, area3=a3, n12, n13, n23, n123,
                 category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"SG_four_Group_Comparison_qval5_FC2.pdf",sep=""),width=12)
draw.quad.venn(area1=sg_a1, area2=sg_a2, area3=sg_a3, area4=sg_a4, sg_n12, sg_n13, sg_n14, sg_n23, sg_n24, sg_n34, sg_n123, sg_n124, sg_n134,
               sg_n234, sg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwC_DwoC"),
               lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()


# q-val <0.1 & FC > 2
a1<-length(sg_NoD_DwoC_qval10_FC2$ID)
a2<-length(sg_NoD_DwC_qval10_FC2$ID)
a3<-length(sg_NoD_Diabetes_qval10_FC2$ID)
a4<-length(sg_DwC_DwoC_qval10_FC2$ID)

n12=length(intersect(sg_NoD_DwoC_qval10_FC2$ID, sg_NoD_DwC_qval10_FC2$ID))
n13=length(intersect(sg_NoD_DwoC_qval10_FC2$ID, sg_NoD_Diabetes_qval10_FC2$ID))
n14=length(intersect(sg_NoD_DwoC_qval10_FC2$ID, sg_DwC_DwoC_qval10_FC2$ID))
n23=length(intersect(NoD_DwC_qval10_FC2$ID, sg_NoD_Diabetes_qval10_FC2$ID))
n24=length(intersect(NoD_DwC_qval10_FC2$ID, sg_DwC_DwoC_qval10_FC2$ID))
n34=length(intersect(NoD_Diabetes_qval10_FC2$ID, sg_DwC_DwoC_qval10_FC2$ID))
n123=length(intersect(intersect(sg_NoD_DwoC_qval10_FC2$ID, sg_NoD_DwC_qval10_FC2$ID), sg_NoD_Diabetes_qval10_FC2$ID))
n124=length(intersect(intersect(sg_NoD_DwoC_qval10_FC2$ID, sg_NoD_DwC_qval10_FC2$ID), sg_DwC_DwoC_qval10_FC2$ID))
n134=length(intersect(intersect(sg_NoD_DwoC_qval10_FC2$ID, sg_NoD_Diabetes_qval10_FC2$ID), sg_DwC_DwoC_qval10_FC2$ID))
n234=length(intersect(intersect(NoD_DwC_qval10_FC2$ID, sg_NoD_Diabetes_qval10_FC2$ID), sg_DwC_DwoC_qval10_FC2$ID))
n1234=length(intersect(intersect(intersect(sg_NoD_DwoC_qval10_FC2$ID, sg_NoD_DwC_qval10_FC2$ID), sg_NoD_Diabetes_qval10_FC2$ID), sg_DwC_DwoC_qval10_FC2$ID))

plot.new()
pdf(paste(Figs,"SG_three_Group_Comparison_qval10_FC2.pdf",sep=""),width=12)
draw.triple.venn(area1=a1, area2=a2, area3=a3, n12, n13, n23, n123,
                 category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"SG_four_Group_Comparison_qval10_FC2.pdf",sep=""),width=12)
draw.quad.venn(area1=a1, area2=a2, area3=a3, area4=a4, n12, n13, n14, n23, n24, n34, n123, n124, n134,
               n234, n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwC_DwoC"),
               lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()

# q-val <0.2 & FC > 2
a1<-length(sg_NoD_DwoC_qval20_FC2$ID)
a2<-length(sg_NoD_DwC_qval20_FC2$ID)
a3<-length(sg_NoD_Diabetes_qval20_FC2$ID)
a4<-length(sg_DwC_DwoC_qval20_FC2$ID)

n12=length(intersect(sg_NoD_DwoC_qval20_FC2$ID, sg_NoD_DwC_qval20_FC2$ID))
n13=length(intersect(sg_NoD_DwoC_qval20_FC2$ID, sg_NoD_Diabetes_qval20_FC2$ID))
n14=length(intersect(sg_NoD_DwoC_qval20_FC2$ID, sg_DwC_DwoC_qval20_FC2$ID))
n23=length(intersect(NoD_DwC_qval20_FC2$ID, sg_NoD_Diabetes_qval20_FC2$ID))
n24=length(intersect(NoD_DwC_qval20_FC2$ID, sg_DwC_DwoC_qval20_FC2$ID))
n34=length(intersect(NoD_Diabetes_qval20_FC2$ID, sg_DwC_DwoC_qval20_FC2$ID))
n123=length(intersect(intersect(sg_NoD_DwoC_qval20_FC2$ID, sg_NoD_DwC_qval20_FC2$ID), sg_NoD_Diabetes_qval20_FC2$ID))
n124=length(intersect(intersect(sg_NoD_DwoC_qval20_FC2$ID, sg_NoD_DwC_qval20_FC2$ID), sg_DwC_DwoC_qval20_FC2$ID))
n134=length(intersect(intersect(sg_NoD_DwoC_qval20_FC2$ID, sg_NoD_Diabetes_qval20_FC2$ID), sg_DwC_DwoC_qval20_FC2$ID))
n234=length(intersect(intersect(NoD_DwC_qval20_FC2$ID, sg_NoD_Diabetes_qval20_FC2$ID), sg_DwC_DwoC_qval20_FC2$ID))
n1234=length(intersect(intersect(intersect(sg_NoD_DwoC_qval20_FC2$ID, sg_NoD_DwC_qval20_FC2$ID), sg_NoD_Diabetes_qval20_FC2$ID), sg_DwC_DwoC_qval20_FC2$ID))

plot.new()
pdf(paste(Figs,"SG_three_Group_Comparison_qval20_FC2.pdf",sep=""),width=12)
draw.triple.venn(area1=sg_a1, area2=sg_a2, area3=sg_a3, sg_n12, sg_n13, sg_n23, sg_n123,
                 category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"SG_four_Group_Comparison_qval20_FC2.pdf",sep=""),width=12)
draw.quad.venn(area1=sg_a1, area2=sg_a2, area3=sg_a3, area4=sg_a4, sg_n12, 
               sg_n13, sg_n14, sg_n23, sg_n24, sg_n34, sg_n123, sg_n124, sg_n134,
               sg_n234, sg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwC_DwoC"),
               lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()

# q-val <0.05 & FC > 1.5
sg_a1<-length(sg_NoD_DwoC_qval5_FC1.5$ID)
sg_a2<-length(sg_NoD_DwC_qval5_FC1.5$ID)
sg_a3<-length(sg_NoD_Diabetes_qval5_FC1.5$ID)
sg_a4<-length(sg_DwC_DwoC_qval5_FC1.5$ID)

sg_n12=length(intersect(sg_NoD_DwoC_qval5_FC1.5$ID, sg_NoD_DwC_qval5_FC1.5$ID))
sg_n13=length(intersect(sg_NoD_DwoC_qval5_FC1.5$ID, sg_NoD_Diabetes_qval5_FC1.5$ID))
sg_n14=length(intersect(sg_NoD_DwoC_qval5_FC1.5$ID, sg_DwC_DwoC_qval5_FC1.5$ID))
sg_n23=length(intersect(sg_NoD_DwC_qval5_FC1.5$ID, sg_NoD_Diabetes_qval5_FC1.5$ID))
sg_n24=length(intersect(sg_NoD_DwC_qval5_FC1.5$ID, sg_DwC_DwoC_qval5_FC1.5$ID))
sg_n34=length(intersect(sg_NoD_Diabetes_qval5_FC1.5$ID, sg_DwC_DwoC_qval5_FC1.5$ID))
sg_n123=length(intersect(intersect(sg_NoD_DwoC_qval5_FC1.5$ID, sg_NoD_DwC_qval5_FC1.5$ID), sg_NoD_Diabetes_qval5_FC1.5$ID))
sg_n124=length(intersect(intersect(sg_NoD_DwoC_qval5_FC1.5$ID, sg_NoD_DwC_qval5_FC1.5$ID), sg_DwC_DwoC_qval5_FC1.5$ID))
sg_n134=length(intersect(intersect(sg_NoD_DwoC_qval5_FC1.5$ID, sg_NoD_Diabetes_qval5_FC1.5$ID), sg_DwC_DwoC_qval5_FC1.5$ID))
sg_n234=length(intersect(intersect(sg_NoD_DwC_qval5_FC1.5$ID, sg_NoD_Diabetes_qval5_FC1.5$ID), sg_DwC_DwoC_qval5_FC1.5$ID))
sg_n1234=length(intersect(intersect(intersect(sg_NoD_DwoC_qval5_FC1.5$ID, sg_NoD_DwC_qval5_FC1.5$ID), sg_NoD_Diabetes_qval5_FC1.5$ID), sg_DwC_DwoC_qval5_FC1.5$ID))

plot.new()
pdf(paste(Figs,"SG_three_Group_Comparison_qval5_FC1.5.pdf",sep=""),width=12)
draw.triple.venn(area1=sg_a1, area2=sg_a2, area3=sg_a3, sg_n12, sg_n13, sg_n23, sg_n123,
                 category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"SG_four_Group_Comparison_qval5_FC1.5.pdf",sep=""),width=12)
draw.quad.venn(area1=sg_a1, area2=sg_a2, area3=sg_a3, area4=sg_a4, sg_n12, sg_n13, sg_n14, sg_n23, sg_n24, sg_n34, sg_n123, sg_n124, sg_n134,
               sg_n234, sg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwC_DwoC"),
               lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()


# q-val <0.1 & FC > 1.5
sg_a1<-length(sg_NoD_DwoC_qval10_FC1.5$ID)
sg_a2<-length(sg_NoD_DwC_qval10_FC1.5$ID)
sg_a3<-length(sg_NoD_Diabetes_qval10_FC1.5$ID)
sg_a4<-length(sg_DwC_DwoC_qval10_FC1.5$ID)

sg_n12=length(intersect(sg_NoD_DwoC_qval10_FC1.5$ID, sg_NoD_DwC_qval10_FC1.5$ID))
sg_n13=length(intersect(sg_NoD_DwoC_qval10_FC1.5$ID, sg_NoD_Diabetes_qval10_FC1.5$ID))
sg_n14=length(intersect(sg_NoD_DwoC_qval10_FC1.5$ID, sg_DwC_DwoC_qval10_FC1.5$ID))
sg_n23=length(intersect(NoD_DwC_qval10_FC1.5$ID, sg_NoD_Diabetes_qval10_FC1.5$ID))
sg_n24=length(intersect(NoD_DwC_qval10_FC1.5$ID, sg_DwC_DwoC_qval10_FC1.5$ID))
sg_n34=length(intersect(NoD_Diabetes_qval10_FC1.5$ID, sg_DwC_DwoC_qval10_FC1.5$ID))
sg_n123=length(intersect(intersect(sg_NoD_DwoC_qval10_FC1.5$ID, sg_NoD_DwC_qval10_FC1.5$ID), sg_NoD_Diabetes_qval10_FC1.5$ID))
sg_n124=length(intersect(intersect(sg_NoD_DwoC_qval10_FC1.5$ID, sg_NoD_DwC_qval10_FC1.5$ID), sg_DwC_DwoC_qval10_FC1.5$ID))
sg_n134=length(intersect(intersect(sg_NoD_DwoC_qval10_FC1.5$ID, sg_NoD_Diabetes_qval10_FC1.5$ID), sg_DwC_DwoC_qval10_FC1.5$ID))
sg_n234=length(intersect(intersect(NoD_DwC_qval10_FC1.5$ID, sg_NoD_Diabetes_qval10_FC1.5$ID), sg_DwC_DwoC_qval10_FC1.5$ID))
sg_n1234=length(intersect(intersect(intersect(sg_NoD_DwoC_qval10_FC1.5$ID, sg_NoD_DwC_qval10_FC1.5$ID), sg_NoD_Diabetes_qval10_FC1.5$ID), sg_DwC_DwoC_qval10_FC1.5$ID))

plot.new()
pdf(paste(Figs,"SG_three_Group_Comparison_qval10_FC1.5.pdf",sep=""),width=12)
draw.triple.venn(area1=a1, area2=a2, area3=a3, n12, n13, n23, n123,
                 category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"SG_four_Group_Comparison_qval10_FC1.5.pdf",sep=""),width=12)
draw.quad.venn(area1=a1, area2=a2, area3=a3, area4=a4, n12, n13, n14, n23, n24, n34, n123, n124, n134,
               n234, n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwC_DwoC"),
               lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()

# q-val <0.2 & FC > 1.5
sg_a1<-length(sg_NoD_DwoC_qval20_FC1.5$ID)
sg_a2<-length(sg_NoD_DwC_qval20_FC1.5$ID)
sg_a3<-length(sg_NoD_Diabetes_qval20_FC1.5$ID)
sg_a4<-length(sg_DwC_DwoC_qval20_FC1.5$ID)

sg_n12=length(intersect(sg_NoD_DwoC_qval20_FC1.5$ID, sg_NoD_DwC_qval20_FC1.5$ID))
sg_n13=length(intersect(sg_NoD_DwoC_qval20_FC1.5$ID, sg_NoD_Diabetes_qval20_FC1.5$ID))
sg_n14=length(intersect(sg_NoD_DwoC_qval20_FC1.5$ID, sg_DwC_DwoC_qval20_FC1.5$ID))
sg_n23=length(intersect(sg_NoD_DwC_qval20_FC1.5$ID, sg_NoD_Diabetes_qval20_FC1.5$ID))
sg_n24=length(intersect(sg_NoD_DwC_qval20_FC1.5$ID, sg_DwC_DwoC_qval20_FC1.5$ID))
sg_n34=length(intersect(sg_NoD_Diabetes_qval20_FC1.5$ID, sg_DwC_DwoC_qval20_FC1.5$ID))
sg_n123=length(intersect(intersect(sg_NoD_DwoC_qval20_FC1.5$ID, sg_NoD_DwC_qval20_FC1.5$ID), sg_NoD_Diabetes_qval20_FC1.5$ID))
sg_n124=length(intersect(intersect(sg_NoD_DwoC_qval20_FC1.5$ID, sg_NoD_DwC_qval20_FC1.5$ID), sg_DwC_DwoC_qval20_FC1.5$ID))
sg_n134=length(intersect(intersect(sg_NoD_DwoC_qval20_FC1.5$ID, sg_NoD_Diabetes_qval20_FC1.5$ID), sg_DwC_DwoC_qval20_FC1.5$ID))
sg_n234=length(intersect(intersect(sg_NoD_DwC_qval20_FC1.5$ID, sg_NoD_Diabetes_qval20_FC1.5$ID), sg_DwC_DwoC_qval20_FC1.5$ID))
sg_n1234=length(intersect(intersect(intersect(sg_NoD_DwoC_qval20_FC1.5$ID, sg_NoD_DwC_qval20_FC1.5$ID), sg_NoD_Diabetes_qval20_FC1.5$ID), sg_DwC_DwoC_qval20_FC1.5$ID))

plot.new()
pdf(paste(Figs,"SG_three_Group_Comparison_qval20_FC1.5.pdf",sep=""),width=12)
draw.triple.venn(area1=a1, area2=a2, area3=a3, n12, n13, n23, n123,
                 category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"SG_four_Group_Comparison_qval20_FC1.5.pdf",sep=""),width=12)
draw.quad.venn(area1=a1, area2=a2, area3=a3, area4=a4, n12, n13, n14, n23, n24, n34, n123, n124, n134,
               n234, n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwC_DwoC"),
               lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()


setEPS()
postscript(file = paste(Figs,"sg_DwC_DwoC_PC1_GrowthRate_pvalue_",pval.cutoff,"_log2FC_",FC.cutoff,"_",date,".eps",sep=""))
with(y1, plot(logFC, -log10(P.Value), pch=20, col="gray", main=paste("Group difference in standard glucose (No_PDR vs. PDR)",sep=""), xlab="log2 fold change", ylab="-log10(P.Value)", xlim=c(-3,3)))
with(subset(y1, P.Value<pval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(y1, P.Value<pval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()

write.csv(y1, paste(Output,"sg_DwC_DwoC_PC1_GrowthRate_all_",date,".csv",sep=""),row.names = F)
write.csv(x1, paste(Output,"sg_DwC_DwoC_PC1_GrowthRate_pvalue_",pval.cutoff,"_",date,".csv",sep=""),row.names = F)
write.csv(x1_subset, paste(Output,"sg_DwC_DwoC_PC1_GrowthRate_pvalue_",pval.cutoff,"_log2FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)

#--------------------------------------
