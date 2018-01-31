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
dat<-local(get(load(file=paste(RData,"normalizedDataMatrix_filtered.RData",sep=""))))

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
## ------------------------------------------------

# 1a. HG - collapse data
Diabete_HG <- cbind(DwC_hg, DwoC_hg, NoD_hg)

data_input = prep_data(Diabete_HG)

# 1b. HG - PCA variance explained
pca_hg <- prcomp(data_input, scale = FALSE, center = TRUE)
pca_hg = pca_hg$rotation
pca_hg = as.data.frame(pca_hg)

plot_pca(pca_hg, "PCA_HG")

# 1c. HG - correlation between PCA and covariates
## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
meta.names = rownames(meta)
pca.names = rownames(pca_hg)
mtch = match(meta.names, pca.names)
pca_hg = pca_hg[mtch,]
pca_hg = as.data.frame(pca_hg)

pca_hg_cor<-correlate_pcs(pca_hg, meta, npcs = 10, min.cor = 0)
pca_hg_cor<-as.data.frame(pca_hg_cor);

plot_cor_pca_cov(pca_hg_cor, "pca_hg_cov_corr_allSamples")

pca_hg_p = pca.meta.regress(pca_hg, meta, no.pcs = 6)
plot_P_pca_cov(pca_hg_p$Ps, "pca_gh_p_allSamples")

write.table(file = "Predict_cov_hg.txt", pca_hg_p$Pred, quote=F, col.names = TRUE, row.names=F)
##
######### REPEAT ABOVE BUT EXCLUDING NO_HG
Diabete_HG <- cbind(DwC_hg, DwoC_hg)

data_input = prep_data(Diabete_HG)

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

plot_cor_pca_cov(pca_hg_no_norm_cor, "PCA_meta_cor_no_norm")

pca_hg_p_no_norm = pca.meta.regress(pca_hg_no_norm, meta, no.pcs = 6)
plot_P_pca_cov(pca_hg_p_no_norm$Ps, "pca_hg_p_no_norm")

write.table(file = "Predict_cov_hg_nonorm.txt", pca_hg_p_no_norm$Pred, quote=F, col.names = TRUE, row.names=F)

### COMPARE PCAS OF DIABETES WHEN NORMAL IS INCLUDED AND EXCLUDED.
PC_cor = calc_pc_cors(pca_hg, pca_hg_no_norm)
write.table(file = "PC_cor_hg_allSamps_noNorm.txt", PC_cor, quote=F, col.names=TRUE, row.names=T)

var.select = pca.meta.regress.2var(pca_hg_no_norm, meta)
write.table(file = "Variable_select_hg.txt", var.select, quote=FALSE, col.names=TRUE, row.names=FALSE)

#####---------------########
# 2a. SG - collapse data
Diabete_SG<-cbind(DwC_sg, DwoC_sg, NoD_sg)
data_input = prep_data(Diabete_SG)

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

plot_cor_pca_cov(pca_sg_cor, "pca_sg_cov_corr_allSamples")

pca_sg_p = pca.meta.regress(pca_sg, meta)
plot_P_pca_cov(pca_sg_p$Ps, "pca_sg_p_allSamples")

write.table(file = "Predict_cov_sg.txt", pca_sg_p$Pred, quote=F, col.names = TRUE, row.names=F)


###
######### REPEAT ABOVE BUT EXCLUDING NO_HG
Diabete_SG <- cbind(DwC_sg, DwoC_sg)
data_input = prep_data(Diabete_SG)

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

plot_cor_pca_cov(pca_sg_no_norm_cor, "PCA_SG_meta_cor_no_norm")

pca_sg_p_no_norm = pca.meta.regress(pca_sg_no_norm, meta, no.pcs = 6)
plot_P_pca_cov(pca_sg_p_no_norm$Ps, "pca_sg_p_noNorm")

write.table(file = "Predict_cov_sg_noNorm.txt", pca_sg_p_no_norm$Pred, quote=F, col.names = TRUE, row.names=F)

### COMPARE PCAS OF DIABETES WHEN NORMAL IS INCLUDED AND EXCLUDED.
PC_cor = calc_pc_cors(pca_sg, pca_sg_no_norm)
write.table(file = "PC_cor_sg_allSamps_noNorm.txt", PC_cor, quote=F, col.names=TRUE, row.names=T)

var.select = pca.meta.regress.2var(pca_sg_no_norm, meta)
write.table(file = "Variable_select_sg.txt", var.select, quote=F, row.names=F, col.names=T)

###------------------------------------------------###
# 3a. SG - collapse data
Diabete_delta<-cbind(DwC_delta, DwoC_delta, NoD_delta)
data_input = prep_data(Diabete_delta)

# 2b. SG - PCA variance explained
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

plot_cor_pca_cov(pca_delta_cor, "pca_delta_cov_corr_allSamples")

pca_delta_p = pca.meta.regress(pca_delta, meta)
plot_P_pca_cov(pca_delta_p$Ps, "pca_delta_p_allSamples")

write.table(file = "Predict_cov_delta.txt", pca_delta_p$Pred, quote=F, col.names = TRUE, row.names=F)

###
######### REPEAT ABOVE BUT EXCLUDING NO_HG
Diabete_delta <- cbind(DwC_delta, DwoC_delta)
data_input = prep_data(Diabete_delta)

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

plot_cor_pca_cov(pca_delta_no_norm_cor, "PCA_delta_meta_cor_no_norm")

pca_delta_p_noNorm = pca.meta.regress(pca_delta_no_norm, meta, no.pcs = 6)
plot_P_pca_cov(pca_delta_p_noNorm$Ps, "pca_delta_p_noNorm")

write.table(file = "Predict_cov_delta_noNorm.txt", pca_delta_p_noNorm$Pred, quote=F, col.names = TRUE, row.names=F)

### COMPARE PCAS OF DIABETES WHEN NORMAL IS INCLUDED AND EXCLUDED.
PC_cor = calc_pc_cors(pca_delta, pca_delta_no_norm)
write.table(file = "PC_cor_delta_allSamps_noNorm.txt", PC_cor, quote=F, col.names=TRUE, row.names=T)

var.select = pca.meta.regress.2var(pca_delta_no_norm, meta)
write.table(file = "Variable_select_delta.txt", var.select, quote=F, row.names=F, col.names=T)