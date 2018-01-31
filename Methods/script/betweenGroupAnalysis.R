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
# hg_NoD_DwoC_qval5_FC2<-x1_subset
# hg_NoD_DwoC_qval10_FC2<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff,]
# hg_NoD_DwoC_qval20_FC2<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff,]
hg_NoD_DwoC_qval5_FC1.5<-x1[x1$adj.P.Val<0.05 & abs(x1$logFC) > FC.cutoff1.5,]
write.csv(hg_NoD_DwoC_qval5_FC1.5, paste(Output,"hg_NoD_DwoC_qval5_FC1.5_",date,".csv",sep=""),row.names = F)
# hg_NoD_DwoC_qval10_FC1.5<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff1.5,]
# hg_NoD_DwoC_qval20_FC1.5<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff1.5,]

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

# hg_NoD_DwC_qval5_FC2<-x1_subset
# hg_NoD_DwC_qval10_FC2<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff,]
# hg_NoD_DwC_qval20_FC2<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff,]
hg_NoD_DwC_qval5_FC1.5<-x1[x1$adj.P.Val<0.05 & abs(x1$logFC) > FC.cutoff1.5,]
write.csv(hg_NoD_DwC_qval5_FC1.5, paste(Output,"hg_NoD_DwC_qval5_FC1.5_",date,".csv",sep=""),row.names = F)
# hg_NoD_DwC_qval10_FC1.5<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff1.5,]
# hg_NoD_DwC_qval20_FC1.5<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff1.5,]

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

# hg_NoD_Diabetes_qval5_FC2<-x1_subset
# hg_NoD_Diabetes_qval10_FC2<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff,]
# hg_NoD_Diabetes_qval20_FC2<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff,]
hg_NoD_Diabetes_qval5_FC1.5<-x1[x1$adj.P.Val<0.05 & abs(x1$logFC) > FC.cutoff1.5,]
write.csv(hg_NoD_Diabetes_qval5_FC1.5, paste(Output,"hg_NoD_Diabetes_qval5_FC1.5_",date,".csv",sep=""),row.names = F)
# hg_NoD_Diabetes_qval10_FC1.5<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff1.5,]
# hg_NoD_Diabetes_qval20_FC1.5<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff1.5,]

#--------------------------------------
# HG DwoC vs DwC
#--------------------------------------
meta_tmp<-read.table(paste(PhenotypeDir, "sample_group_DwoC_vs_DwC.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
group <-as.factor(meta_tmp$GROUP)
growth_rate<-meta[rownames(meta_tmp),"GROWTH_RATE"]
data_input_tmp<-data_input[,rownames(meta_tmp)]; data_input_tmp<-apply(data_input_tmp,2,as.numeric); rownames(data_input_tmp)<-rownames(data_input);

# PC1 from data after group regressed out
data_input2_tmp<-data_input2[,rownames(meta_tmp)]; data_input2_tmp<-apply(data_input2_tmp,2,as.numeric); rownames(data_input2_tmp)<-rownames(data_input2);

pca_DwoC_DwC <- prcomp(data_input2_tmp, scale = FALSE, center = TRUE)
pca_DwoC_DwC = pca_DwoC_DwC$rotation
pca_DwoC_DwC = as.data.frame(pca_DwoC_DwC)

pc<-pca_DwoC_DwC$PC1
design<-model.matrix(~pc+growth_rate+group)
fit <- eBayes(lmFit(data_input_tmp, design))

pval.cutoff=0.05; FC.cutoff=1 # FC=2
x1=topTable(fit, coef="group1", n=nrow(genes), adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="group1", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$P.Value<pval.cutoff & abs(x1$logFC) > FC.cutoff,]

# hg_DwoC_DwC_qval5_FC2<-x1_subset
# hg_DwoC_DwC_qval10_FC2<-x1[x1$P.Value<0.1 & abs(x1$logFC) > FC.cutoff,]
# hg_DwoC_DwC_qval20_FC2<-x1[x1$P.Value<0.2 & abs(x1$logFC) > FC.cutoff,]
hg_DwoC_DwC_qval5_FC1.5<-x1[x1$P.Value<0.05 & abs(x1$logFC) > FC.cutoff1.5,]
write.csv(hg_DwoC_DwC_qval5_FC1.5, paste(Output,"hg_DwoC_DwC_qval5_FC1.5_",date,".csv",sep=""),row.names = F)
# hg_DwoC_DwC_qval10_FC1.5<-x1[x1$P.Value<0.1 & abs(x1$logFC) > FC.cutoff1.5,]
# hg_DwoC_DwC_qval20_FC1.5<-x1[x1$P.Value<0.2 & abs(x1$logFC) > FC.cutoff1.5,]

## find common genes
library(VennDiagram)

# q-val <0.05 & FC > 2
# hg_a1<-length(hg_NoD_DwoC_qval5_FC2$ID)
# hg_a2<-length(hg_NoD_DwC_qval5_FC2$ID)
# hg_a3<-length(hg_NoD_Diabetes_qval5_FC2$ID)
# hg_a4<-length(hg_DwoC_DwC_qval5_FC2$ID)
# 
# hg_n12=length(intersect(hg_NoD_DwoC_qval5_FC2$ID, hg_NoD_DwC_qval5_FC2$ID))
# hg_n13=length(intersect(hg_NoD_DwoC_qval5_FC2$ID, hg_NoD_Diabetes_qval5_FC2$ID))
# hg_n14=length(intersect(hg_NoD_DwoC_qval5_FC2$ID, hg_DwoC_DwC_qval5_FC2$ID))
# hg_n23=length(intersect(hg_NoD_DwC_qval5_FC2$ID, hg_NoD_Diabetes_qval5_FC2$ID))
# hg_n24=length(intersect(hg_NoD_DwC_qval5_FC2$ID, hg_DwoC_DwC_qval5_FC2$ID))
# hg_n34=length(intersect(hg_NoD_Diabetes_qval5_FC2$ID, hg_DwoC_DwC_qval5_FC2$ID))
# hg_n123=length(intersect(intersect(hg_NoD_DwoC_qval5_FC2$ID, hg_NoD_DwC_qval5_FC2$ID), hg_NoD_Diabetes_qval5_FC2$ID))
# hg_n124=length(intersect(intersect(hg_NoD_DwoC_qval5_FC2$ID, hg_NoD_DwC_qval5_FC2$ID), hg_DwoC_DwC_qval5_FC2$ID))
# hg_n134=length(intersect(intersect(hg_NoD_DwoC_qval5_FC2$ID, hg_NoD_Diabetes_qval5_FC2$ID), hg_DwoC_DwC_qval5_FC2$ID))
# hg_n234=length(intersect(intersect(hg_NoD_DwC_qval5_FC2$ID, hg_NoD_Diabetes_qval5_FC2$ID), hg_DwoC_DwC_qval5_FC2$ID))
# hg_n1234=length(intersect(intersect(intersect(hg_NoD_DwoC_qval5_FC2$ID, hg_NoD_DwC_qval5_FC2$ID), hg_NoD_Diabetes_qval5_FC2$ID), hg_DwoC_DwC_qval5_FC2$ID))
# 
# plot.new()
# pdf(paste(Figs,"HG_three_Group_Comparison_qval5_FC2.pdf",sep=""),width=12)
# draw.triple.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, hg_n12, hg_n13, hg_n23, hg_n123,
#                  category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
#                  lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
# dev.off()
# 
# plot.new()
# pdf(paste(Figs,"HG_four_Group_Comparison_qval5_FC2.pdf",sep=""),width=12)
# draw.quad.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, area4=hg_a4, hg_n12, hg_n13, hg_n14, hg_n23, 
#                hg_n24, hg_n34, hg_n123, hg_n124, hg_n134,
#                hg_n234, hg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwoC_DwC"),
#                lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
# dev.off()
# 
# 
# # q-val <0.1 & FC > 2
# hg_a1<-length(hg_NoD_DwoC_qval10_FC2$ID)
# hg_a2<-length(hg_NoD_DwC_qval10_FC2$ID)
# hg_a3<-length(hg_NoD_Diabetes_qval10_FC2$ID)
# hg_a4<-length(hg_DwoC_DwC_qval10_FC2$ID)
# 
# hg_n12=length(intersect(hg_NoD_DwoC_qval10_FC2$ID, hg_NoD_DwC_qval10_FC2$ID))
# hg_n13=length(intersect(hg_NoD_DwoC_qval10_FC2$ID, hg_NoD_Diabetes_qval10_FC2$ID))
# hg_n14=length(intersect(hg_NoD_DwoC_qval10_FC2$ID, hg_DwoC_DwC_qval10_FC2$ID))
# hg_n23=length(intersect(hg_NoD_DwC_qval10_FC2$ID, hg_NoD_Diabetes_qval10_FC2$ID))
# hg_n24=length(intersect(hg_NoD_DwC_qval10_FC2$ID, hg_DwoC_DwC_qval10_FC2$ID))
# hg_n34=length(intersect(hg_NoD_Diabetes_qval10_FC2$ID, hg_DwoC_DwC_qval10_FC2$ID))
# hg_n123=length(intersect(intersect(hg_NoD_DwoC_qval10_FC2$ID, hg_NoD_DwC_qval10_FC2$ID), hg_NoD_Diabetes_qval10_FC2$ID))
# hg_n124=length(intersect(intersect(hg_NoD_DwoC_qval10_FC2$ID, hg_NoD_DwC_qval10_FC2$ID), hg_DwoC_DwC_qval10_FC2$ID))
# hg_n134=length(intersect(intersect(hg_NoD_DwoC_qval10_FC2$ID, hg_NoD_Diabetes_qval10_FC2$ID),hg_DwoC_DwC_qval10_FC2$ID))
# hg_n234=length(intersect(intersect(hg_NoD_DwC_qval10_FC2$ID, hg_NoD_Diabetes_qval10_FC2$ID), hg_DwoC_DwC_qval10_FC2$ID))
# hg_n1234=length(intersect(intersect(intersect(hg_NoD_DwoC_qval10_FC2$ID, hg_NoD_DwC_qval10_FC2$ID), hg_NoD_Diabetes_qval10_FC2$ID), hg_DwoC_DwC_qval10_FC2$ID))
# 
# plot.new()
# pdf(paste(Figs,"HG_three_Group_Comparison_qval10_FC2.pdf",sep=""),width=12)
# draw.triple.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, hg_n12, hg_n13, hg_n23, hg_n123,
#                  category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
#                  lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
# dev.off()
# 
# plot.new()
# pdf(paste(Figs,"HG_four_Group_Comparison_qval10_FC2.pdf",sep=""),width=12)
# draw.quad.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, area4=hg_a4, hg_n12, hg_n13, hg_n14, hg_n23, hg_n24, hg_n34, hg_n123, hg_n124, hg_n134,
#                hg_n234, hg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwoC_DwC"),
#                lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
# dev.off()
# 
# # q-val <0.2 & FC > 2
# hg_a1<-length(hg_NoD_DwoC_qval20_FC2$ID)
# hg_a2<-length(hg_NoD_DwC_qval20_FC2$ID)
# hg_a3<-length(hg_NoD_Diabetes_qval20_FC2$ID)
# hg_a4<-length(hg_DwoC_DwC_qval20_FC2$ID)
# 
# hg_n12=length(intersect(hg_NoD_DwoC_qval20_FC2$ID, hg_NoD_DwC_qval20_FC2$ID))
# hg_n13=length(intersect(hg_NoD_DwoC_qval20_FC2$ID, hg_NoD_Diabetes_qval20_FC2$ID))
# hg_n14=length(intersect(hg_NoD_DwoC_qval20_FC2$ID, hg_DwoC_DwC_qval20_FC2$ID))
# hg_n23=length(intersect(hg_NoD_DwC_qval20_FC2$ID, hg_NoD_Diabetes_qval20_FC2$ID))
# hg_n24=length(intersect(hg_NoD_DwC_qval20_FC2$ID, hg_DwoC_DwC_qval20_FC2$ID))
# hg_n34=length(intersect(hg_NoD_Diabetes_qval20_FC2$ID, hg_DwoC_DwC_qval20_FC2$ID))
# hg_n123=length(intersect(intersect(hg_NoD_DwoC_qval20_FC2$ID, hg_NoD_DwC_qval20_FC2$ID), hg_NoD_Diabetes_qval20_FC2$ID))
# hg_n124=length(intersect(intersect(hg_NoD_DwoC_qval20_FC2$ID, hg_NoD_DwC_qval20_FC2$ID), hg_DwoC_DwC_qval20_FC2$ID))
# hg_n134=length(intersect(intersect(hg_NoD_DwoC_qval20_FC2$ID, hg_NoD_Diabetes_qval20_FC2$ID), hg_DwoC_DwC_qval20_FC2$ID))
# hg_n234=length(intersect(intersect(hg_NoD_DwC_qval20_FC2$ID, hg_NoD_Diabetes_qval20_FC2$ID), hg_DwoC_DwC_qval20_FC2$ID))
# hg_n1234=length(intersect(intersect(intersect(hg_NoD_DwoC_qval20_FC2$ID, hg_NoD_DwC_qval20_FC2$ID), hg_NoD_Diabetes_qval20_FC2$ID), hg_DwoC_DwC_qval20_FC2$ID))
# 
# plot.new()
# pdf(paste(Figs,"HG_three_Group_Comparison_qval20_FC2.pdf",sep=""),width=12)
# draw.triple.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, hg_n12, hg_n13, hg_n23, hg_n123,
#                  category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
#                  lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
# dev.off()
# 
# plot.new()
# pdf(paste(Figs,"HG_four_Group_Comparison_qval20_FC2.pdf",sep=""),width=12)
# draw.quad.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, area4=hg_a4, hg_n12, hg_n13, hg_n14, hg_n23, hg_n24, hg_n34,hg_n123, hg_n124, hg_n134,
#                hg_n234, hg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwoC_DwC"),
#                lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
# dev.off()

# q-val <0.05 & FC > 1.5
hg_a1<-length(hg_NoD_DwoC_qval5_FC1.5$ID)
hg_a2<-length(hg_NoD_DwC_qval5_FC1.5$ID)
hg_a3<-length(hg_NoD_Diabetes_qval5_FC1.5$ID)
hg_a4<-length(hg_DwoC_DwC_qval5_FC1.5$ID)

hg_n12=length(intersect(hg_NoD_DwoC_qval5_FC1.5$ID, hg_NoD_DwC_qval5_FC1.5$ID))
hg_n13=length(intersect(hg_NoD_DwoC_qval5_FC1.5$ID, hg_NoD_Diabetes_qval5_FC1.5$ID))
hg_n14=length(intersect(hg_NoD_DwoC_qval5_FC1.5$ID, hg_DwoC_DwC_qval5_FC1.5$ID))
hg_n23=length(intersect(hg_NoD_DwC_qval5_FC1.5$ID, hg_NoD_Diabetes_qval5_FC1.5$ID))
hg_n24=length(intersect(hg_NoD_DwC_qval5_FC1.5$ID, hg_DwoC_DwC_qval5_FC1.5$ID))
hg_n34=length(intersect(hg_NoD_Diabetes_qval5_FC1.5$ID, hg_DwoC_DwC_qval5_FC1.5$ID))
hg_n123=length(intersect(intersect(hg_NoD_DwoC_qval5_FC1.5$ID, hg_NoD_DwC_qval5_FC1.5$ID), hg_NoD_Diabetes_qval5_FC1.5$ID))
hg_n124=length(intersect(intersect(hg_NoD_DwoC_qval5_FC1.5$ID, hg_NoD_DwC_qval5_FC1.5$ID), hg_DwoC_DwC_qval5_FC1.5$ID))
hg_n134=length(intersect(intersect(hg_NoD_DwoC_qval5_FC1.5$ID, hg_NoD_Diabetes_qval5_FC1.5$ID), hg_DwoC_DwC_qval5_FC1.5$ID))
hg_n234=length(intersect(intersect(hg_NoD_DwC_qval5_FC1.5$ID, hg_NoD_Diabetes_qval5_FC1.5$ID), hg_DwoC_DwC_qval5_FC1.5$ID))
hg_n1234=length(intersect(intersect(intersect(hg_NoD_DwoC_qval5_FC1.5$ID, hg_NoD_DwC_qval5_FC1.5$ID), hg_NoD_Diabetes_qval5_FC1.5$ID), hg_DwoC_DwC_qval5_FC1.5$ID))

plot.new()
pdf(paste(Figs,"HG_three_Group_Comparison_qval5_FC1.5.pdf",sep=""),width=12)
draw.triple.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, hg_n12, hg_n13, hg_n23, hg_n123,
                 category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"HG_four_Group_Comparison_qval5_FC1.5.pdf",sep=""),width=12)
draw.quad.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, area4=hg_a4, hg_n12, hg_n13, hg_n14, hg_n23, hg_n24, hg_n34, hg_n123, hg_n124, hg_n134,
               hg_n234, hg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwoC_DwC"),
               lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()


# # q-val <0.1 & FC > 1.5
# hg_a1<-length(hg_NoD_DwoC_qval10_FC1.5$ID)
# hg_a2<-length(hg_NoD_DwC_qval10_FC1.5$ID)
# hg_a3<-length(hg_NoD_Diabetes_qval10_FC1.5$ID)
# hg_a4<-length(hg_DwoC_DwC_qval10_FC1.5$ID)
# 
# hg_n12=length(intersect(hg_NoD_DwoC_qval10_FC1.5$ID, hg_NoD_DwC_qval10_FC1.5$ID))
# hg_n13=length(intersect(hg_NoD_DwoC_qval10_FC1.5$ID, hg_NoD_Diabetes_qval10_FC1.5$ID))
# hg_n14=length(intersect(hg_NoD_DwoC_qval10_FC1.5$ID, hg_DwoC_DwC_qval10_FC1.5$ID))
# hg_n23=length(intersect(hg_NoD_DwC_qval10_FC1.5$ID, hg_NoD_Diabetes_qval10_FC1.5$ID))
# hg_n24=length(intersect(hg_NoD_DwC_qval10_FC1.5$ID, hg_DwoC_DwC_qval10_FC1.5$ID))
# hg_n34=length(intersect(hg_NoD_Diabetes_qval10_FC1.5$ID, hg_DwoC_DwC_qval10_FC1.5$ID))
# hg_n123=length(intersect(intersect(hg_NoD_DwoC_qval10_FC1.5$ID, hg_NoD_DwC_qval10_FC1.5$ID), hg_NoD_Diabetes_qval10_FC1.5$ID))
# hg_n124=length(intersect(intersect(hg_NoD_DwoC_qval10_FC1.5$ID, hg_NoD_DwC_qval10_FC1.5$ID), hg_DwoC_DwC_qval10_FC1.5$ID))
# hg_n134=length(intersect(intersect(hg_NoD_DwoC_qval10_FC1.5$ID, hg_NoD_Diabetes_qval10_FC1.5$ID), hg_DwoC_DwC_qval10_FC1.5$ID))
# hg_n234=length(intersect(intersect(hg_NoD_DwC_qval10_FC1.5$ID, hg_NoD_Diabetes_qval10_FC1.5$ID), hg_DwoC_DwC_qval10_FC1.5$ID))
# hg_n1234=length(intersect(intersect(intersect(hg_NoD_DwoC_qval10_FC1.5$ID, hg_NoD_DwC_qval10_FC1.5$ID), hg_NoD_Diabetes_qval10_FC1.5$ID), hg_DwoC_DwC_qval10_FC1.5$ID))
# 
# plot.new()
# pdf(paste(Figs,"HG_three_Group_Comparison_qval10_FC1.5.pdf",sep=""),width=12)
# draw.triple.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, hg_n12, hg_n13, hg_n23, hg_n123,
#                  category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
#                  lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
# dev.off()
# 
# plot.new()
# pdf(paste(Figs,"HG_four_Group_Comparison_qval10_FC1.5.pdf",sep=""),width=12)
# draw.quad.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, area4=hg_a4, hg_n12, hg_n13, hg_n14, hg_n23, hg_n24, hg_n34, hg_n123, hg_n124, hg_n134,
#                hg_n234, hg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwoC_DwC"),
#                lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
# dev.off()
# 
# # q-val <0.2 & FC > 1.5
# hg_a1<-length(hg_NoD_DwoC_qval20_FC1.5$ID)
# hg_a2<-length(hg_NoD_DwC_qval20_FC1.5$ID)
# hg_a3<-length(hg_NoD_Diabetes_qval20_FC1.5$ID)
# hg_a4<-length(hg_DwoC_DwC_qval20_FC1.5$ID)
# 
# hg_n12=length(intersect(hg_NoD_DwoC_qval20_FC1.5$ID, hg_NoD_DwC_qval20_FC1.5$ID))
# hg_n13=length(intersect(hg_NoD_DwoC_qval20_FC1.5$ID, hg_NoD_Diabetes_qval20_FC1.5$ID))
# hg_n14=length(intersect(hg_NoD_DwoC_qval20_FC1.5$ID, hg_DwoC_DwC_qval20_FC1.5$ID))
# hg_n23=length(intersect(hg_NoD_DwC_qval20_FC1.5$ID, hg_NoD_Diabetes_qval20_FC1.5$ID))
# hg_n24=length(intersect(hg_NoD_DwC_qval20_FC1.5$ID, hg_DwoC_DwC_qval20_FC1.5$ID))
# hg_n34=length(intersect(hg_NoD_Diabetes_qval20_FC1.5$ID, hg_DwoC_DwC_qval20_FC1.5$ID))
# hg_n123=length(intersect(intersect(hg_NoD_DwoC_qval20_FC1.5$ID, hg_NoD_DwC_qval20_FC1.5$ID), hg_NoD_Diabetes_qval20_FC1.5$ID))
# hg_n124=length(intersect(intersect(hg_NoD_DwoC_qval20_FC1.5$ID, hg_NoD_DwC_qval20_FC1.5$ID), hg_DwoC_DwC_qval20_FC1.5$ID))
# hg_n134=length(intersect(intersect(hg_NoD_DwoC_qval20_FC1.5$ID, hg_NoD_Diabetes_qval20_FC1.5$ID), hg_DwoC_DwC_qval20_FC1.5$ID))
# hg_n234=length(intersect(intersect(hg_NoD_DwC_qval20_FC1.5$ID, hg_NoD_Diabetes_qval20_FC1.5$ID), hg_DwoC_DwC_qval20_FC1.5$ID))
# hg_n1234=length(intersect(intersect(intersect(hg_NoD_DwoC_qval20_FC1.5$ID, hg_NoD_DwC_qval20_FC1.5$ID), hg_NoD_Diabetes_qval20_FC1.5$ID), hg_DwoC_DwC_qval20_FC1.5$ID))
# 
# plot.new()
# pdf(paste(Figs,"HG_three_Group_Comparison_qval20_FC1.5.pdf",sep=""),width=12)
# draw.triple.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, hg_n12, hg_n13, hg_n23, hg_n123,
#                  category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
#                  lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
# dev.off()
# 
# plot.new()
# pdf(paste(Figs,"HG_four_Group_Comparison_qval20_FC1.5.pdf",sep=""),width=12)
# draw.quad.venn(area1=hg_a1, area2=hg_a2, area3=hg_a3, area4=hg_a4, hg_n12, hg_n13, hg_n14, hg_n23, hg_n24, hg_n34, hg_n123, hg_n124, hg_n134,
#                hg_n234, hg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwoC_DwC"),
#                lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
# dev.off()

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
# sg_NoD_DwoC_qval5_FC2<-x1_subset
# sg_NoD_DwoC_qval10_FC2<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff,]
# sg_NoD_DwoC_qval20_FC2<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff,]
sg_NoD_DwoC_qval5_FC1.5<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff1.5,]
write.csv(sg_NoD_DwoC_qval5_FC1.5, paste(Output,"sg_NoD_DwoC_qval5_FC1.5_",date,".csv",sep=""),row.names = F)
# sg_NoD_DwoC_qval10_FC1.5<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff1.5,]
# sg_NoD_DwoC_qval20_FC1.5<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff1.5,]

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

# sg_NoD_DwC_qval5_FC2<-x1_subset
# sg_NoD_DwC_qval10_FC2<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff,]
# sg_NoD_DwC_qval20_FC2<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff,]
sg_NoD_DwC_qval5_FC1.5<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff1.5,]
write.csv(sg_NoD_DwC_qval5_FC1.5, paste(Output,"sg_NoD_DwC_qval5_FC1.5_",date,".csv",sep=""),row.names = F)
# sg_NoD_DwC_qval10_FC1.5<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff1.5,]
# sg_NoD_DwC_qval20_FC1.5<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff1.5,]

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

# sg_NoD_Diabetes_qval5_FC2<-x1_subset
# sg_NoD_Diabetes_qval10_FC2<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff,]
# sg_NoD_Diabetes_qval20_FC2<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff,]
sg_NoD_Diabetes_qval5_FC1.5<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff1.5,]
write.csv(sg_NoD_Diabetes_qval5_FC1.5, paste(Output,"sg_NoD_Diabetes_qval5_FC1.5_",date,".csv",sep=""),row.names = F)
# sg_NoD_Diabetes_qval10_FC1.5<-x1[x1$adj.P.Val<0.1 & abs(x1$logFC) > FC.cutoff1.5,]
# sg_NoD_Diabetes_qval20_FC1.5<-x1[x1$adj.P.Val<0.2 & abs(x1$logFC) > FC.cutoff1.5,]

#--------------------------------------
# SG DwoC vs DwC
#--------------------------------------
meta_tmp<-read.table(paste(PhenotypeDir, "sample_group_DwoC_vs_DwC.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
group <-as.factor(meta_tmp$GROUP)
growth_rate<-meta[rownames(meta_tmp),"GROWTH_RATE"]
data_input_tmp<-data_input[,rownames(meta_tmp)]; data_input_tmp<-apply(data_input_tmp,2,as.numeric); rownames(data_input_tmp)<-rownames(data_input);

# PC1 from data after group regressed out
data_input2_tmp<-data_input2[,rownames(meta_tmp)]; data_input2_tmp<-apply(data_input2_tmp,2,as.numeric); rownames(data_input2_tmp)<-rownames(data_input2);

pca_DwoC_DwC <- prcomp(data_input2_tmp, scale = FALSE, center = TRUE)
pca_DwoC_DwC = pca_DwoC_DwC$rotation
pca_DwoC_DwC = as.data.frame(pca_DwoC_DwC)

pc<-pca_DwoC_DwC$PC1
design<-model.matrix(~pc+growth_rate+group)
fit <- eBayes(lmFit(data_input_tmp, design))

pval.cutoff=0.05; FC.cutoff=1 # FC=2
x1=topTable(fit, coef="group1", n=nrow(genes), adjust.method="BH",genelist=genes)
y1=topTable(fit, coef="group1", n=nrow(genes),adjust.method="BH",genelist=genes)
x1_subset<-x1[x1$P.Value<pval.cutoff & abs(x1$logFC) > FC.cutoff,]

# sg_DwoC_DwC_qval5_FC2<-x1_subset
# sg_DwoC_DwC_qval10_FC2<-x1[x1$P.Value<0.1 & abs(x1$logFC) > FC.cutoff,]
# sg_DwoC_DwC_qval20_FC2<-x1[x1$P.Value<0.2 & abs(x1$logFC) > FC.cutoff,]
sg_DwoC_DwC_qval5_FC1.5<-x1[x1$P.Value<qval.cutoff & abs(x1$logFC) > FC.cutoff1.5,]
write.csv(sg_DwoC_DwC_qval5_FC1.5, paste(Output,"sg_DwoC_DwC_qval5_FC1.5_",date,".csv",sep=""),row.names = F)
# sg_DwoC_DwC_qval10_FC1.5<-x1[x1$P.Value<0.1 & abs(x1$logFC) > FC.cutoff1.5,]
# sg_DwoC_DwC_qval20_FC1.5<-x1[x1$P.Value<0.2 & abs(x1$logFC) > FC.cutoff1.5,]

# # q-val <0.05 & FC > 2
# sg_a1<-length(sg_NoD_DwoC_qval5_FC2$ID)
# sg_a2<-length(sg_NoD_DwC_qval5_FC2$ID)
# sg_a3<-length(sg_NoD_Diabetes_qval5_FC2$ID)
# sg_a4<-length(sg_DwoC_DwC_qval5_FC2$ID)
# 
# sg_n12=length(intersect(sg_NoD_DwoC_qval5_FC2$ID, sg_NoD_DwC_qval5_FC2$ID))
# sg_n13=length(intersect(sg_NoD_DwoC_qval5_FC2$ID, sg_NoD_Diabetes_qval5_FC2$ID))
# sg_n14=length(intersect(sg_NoD_DwoC_qval5_FC2$ID, sg_DwoC_DwC_qval5_FC2$ID))
# sg_n23=length(intersect(sg_NoD_DwC_qval5_FC2$ID, sg_NoD_Diabetes_qval5_FC2$ID))
# sg_n24=length(intersect(sg_NoD_DwC_qval5_FC2$ID, sg_DwoC_DwC_qval5_FC2$ID))
# sg_n34=length(intersect(sg_NoD_Diabetes_qval5_FC2$ID, sg_DwoC_DwC_qval5_FC2$ID))
# sg_n123=length(intersect(intersect(sg_NoD_DwoC_qval5_FC2$ID, sg_NoD_DwC_qval5_FC2$ID), sg_NoD_Diabetes_qval5_FC2$ID))
# sg_n124=length(intersect(intersect(sg_NoD_DwoC_qval5_FC2$ID, sg_NoD_DwC_qval5_FC2$ID), sg_DwoC_DwC_qval5_FC2$ID))
# sg_n134=length(intersect(intersect(sg_NoD_DwoC_qval5_FC2$ID, sg_NoD_Diabetes_qval5_FC2$ID), sg_DwoC_DwC_qval5_FC2$ID))
# sg_n234=length(intersect(intersect(sg_NoD_DwC_qval5_FC2$ID, sg_NoD_Diabetes_qval5_FC2$ID), sg_DwoC_DwC_qval5_FC2$ID))
# sg_n1234=length(intersect(intersect(intersect(sg_NoD_DwoC_qval5_FC2$ID, sg_NoD_DwC_qval5_FC2$ID), sg_NoD_Diabetes_qval5_FC2$ID), sg_DwoC_DwC_qval5_FC2$ID))
# 
# plot.new()
# pdf(paste(Figs,"SG_three_Group_Comparison_qval5_FC2.pdf",sep=""),width=12)
# draw.triple.venn(area1=a1, area2=a2, area3=a3, n12, n13, n23, n123,
#                  category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
#                  lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
# dev.off()
# 
# plot.new()
# pdf(paste(Figs,"SG_four_Group_Comparison_qval5_FC2.pdf",sep=""),width=12)
# draw.quad.venn(area1=sg_a1, area2=sg_a2, area3=sg_a3, area4=sg_a4, sg_n12, sg_n13, sg_n14, sg_n23, sg_n24, sg_n34, sg_n123, sg_n124, sg_n134,
#                sg_n234, sg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwoC_DwC"),
#                lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
# dev.off()
# 
# 
# # q-val <0.1 & FC > 2
# a1<-length(sg_NoD_DwoC_qval10_FC2$ID)
# a2<-length(sg_NoD_DwC_qval10_FC2$ID)
# a3<-length(sg_NoD_Diabetes_qval10_FC2$ID)
# a4<-length(sg_DwoC_DwC_qval10_FC2$ID)
# 
# n12=length(intersect(sg_NoD_DwoC_qval10_FC2$ID, sg_NoD_DwC_qval10_FC2$ID))
# n13=length(intersect(sg_NoD_DwoC_qval10_FC2$ID, sg_NoD_Diabetes_qval10_FC2$ID))
# n14=length(intersect(sg_NoD_DwoC_qval10_FC2$ID, sg_DwoC_DwC_qval10_FC2$ID))
# n23=length(intersect(NoD_DwC_qval10_FC2$ID, sg_NoD_Diabetes_qval10_FC2$ID))
# n24=length(intersect(NoD_DwC_qval10_FC2$ID, sg_DwoC_DwC_qval10_FC2$ID))
# n34=length(intersect(NoD_Diabetes_qval10_FC2$ID, sg_DwoC_DwC_qval10_FC2$ID))
# n123=length(intersect(intersect(sg_NoD_DwoC_qval10_FC2$ID, sg_NoD_DwC_qval10_FC2$ID), sg_NoD_Diabetes_qval10_FC2$ID))
# n124=length(intersect(intersect(sg_NoD_DwoC_qval10_FC2$ID, sg_NoD_DwC_qval10_FC2$ID), sg_DwoC_DwC_qval10_FC2$ID))
# n134=length(intersect(intersect(sg_NoD_DwoC_qval10_FC2$ID, sg_NoD_Diabetes_qval10_FC2$ID), sg_DwoC_DwC_qval10_FC2$ID))
# n234=length(intersect(intersect(NoD_DwC_qval10_FC2$ID, sg_NoD_Diabetes_qval10_FC2$ID), sg_DwoC_DwC_qval10_FC2$ID))
# n1234=length(intersect(intersect(intersect(sg_NoD_DwoC_qval10_FC2$ID, sg_NoD_DwC_qval10_FC2$ID), sg_NoD_Diabetes_qval10_FC2$ID), sg_DwoC_DwC_qval10_FC2$ID))
# 
# plot.new()
# pdf(paste(Figs,"SG_three_Group_Comparison_qval10_FC2.pdf",sep=""),width=12)
# draw.triple.venn(area1=a1, area2=a2, area3=a3, n12, n13, n23, n123,
#                  category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
#                  lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
# dev.off()
# 
# plot.new()
# pdf(paste(Figs,"SG_four_Group_Comparison_qval10_FC2.pdf",sep=""),width=12)
# draw.quad.venn(area1=a1, area2=a2, area3=a3, area4=a4, n12, n13, n14, n23, n24, n34, n123, n124, n134,
#                n234, n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwoC_DwC"),
#                lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
# dev.off()
# 
# # q-val <0.2 & FC > 2
# a1<-length(sg_NoD_DwoC_qval20_FC2$ID)
# a2<-length(sg_NoD_DwC_qval20_FC2$ID)
# a3<-length(sg_NoD_Diabetes_qval20_FC2$ID)
# a4<-length(sg_DwoC_DwC_qval20_FC2$ID)
# 
# n12=length(intersect(sg_NoD_DwoC_qval20_FC2$ID, sg_NoD_DwC_qval20_FC2$ID))
# n13=length(intersect(sg_NoD_DwoC_qval20_FC2$ID, sg_NoD_Diabetes_qval20_FC2$ID))
# n14=length(intersect(sg_NoD_DwoC_qval20_FC2$ID, sg_DwoC_DwC_qval20_FC2$ID))
# n23=length(intersect(NoD_DwC_qval20_FC2$ID, sg_NoD_Diabetes_qval20_FC2$ID))
# n24=length(intersect(NoD_DwC_qval20_FC2$ID, sg_DwoC_DwC_qval20_FC2$ID))
# n34=length(intersect(NoD_Diabetes_qval20_FC2$ID, sg_DwoC_DwC_qval20_FC2$ID))
# n123=length(intersect(intersect(sg_NoD_DwoC_qval20_FC2$ID, sg_NoD_DwC_qval20_FC2$ID), sg_NoD_Diabetes_qval20_FC2$ID))
# n124=length(intersect(intersect(sg_NoD_DwoC_qval20_FC2$ID, sg_NoD_DwC_qval20_FC2$ID), sg_DwoC_DwC_qval20_FC2$ID))
# n134=length(intersect(intersect(sg_NoD_DwoC_qval20_FC2$ID, sg_NoD_Diabetes_qval20_FC2$ID), sg_DwoC_DwC_qval20_FC2$ID))
# n234=length(intersect(intersect(NoD_DwC_qval20_FC2$ID, sg_NoD_Diabetes_qval20_FC2$ID), sg_DwoC_DwC_qval20_FC2$ID))
# n1234=length(intersect(intersect(intersect(sg_NoD_DwoC_qval20_FC2$ID, sg_NoD_DwC_qval20_FC2$ID), sg_NoD_Diabetes_qval20_FC2$ID), sg_DwoC_DwC_qval20_FC2$ID))
# 
# plot.new()
# pdf(paste(Figs,"SG_three_Group_Comparison_qval20_FC2.pdf",sep=""),width=12)
# draw.triple.venn(area1=sg_a1, area2=sg_a2, area3=sg_a3, sg_n12, sg_n13, sg_n23, sg_n123,
#                  category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
#                  lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
# dev.off()
# 
# plot.new()
# pdf(paste(Figs,"SG_four_Group_Comparison_qval20_FC2.pdf",sep=""),width=12)
# draw.quad.venn(area1=sg_a1, area2=sg_a2, area3=sg_a3, area4=sg_a4, sg_n12, 
#                sg_n13, sg_n14, sg_n23, sg_n24, sg_n34, sg_n123, sg_n124, sg_n134,
#                sg_n234, sg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwoC_DwC"),
#                lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
# dev.off()

# q-val <0.05 & FC > 1.5
sg_a1<-length(sg_NoD_DwoC_qval5_FC1.5$ID)
sg_a2<-length(sg_NoD_DwC_qval5_FC1.5$ID)
sg_a3<-length(sg_NoD_Diabetes_qval5_FC1.5$ID)
sg_a4<-length(sg_DwoC_DwC_qval5_FC1.5$ID)

sg_n12=length(intersect(sg_NoD_DwoC_qval5_FC1.5$ID, sg_NoD_DwC_qval5_FC1.5$ID))
sg_n13=length(intersect(sg_NoD_DwoC_qval5_FC1.5$ID, sg_NoD_Diabetes_qval5_FC1.5$ID))
sg_n14=length(intersect(sg_NoD_DwoC_qval5_FC1.5$ID, sg_DwoC_DwC_qval5_FC1.5$ID))
sg_n23=length(intersect(sg_NoD_DwC_qval5_FC1.5$ID, sg_NoD_Diabetes_qval5_FC1.5$ID))
sg_n24=length(intersect(sg_NoD_DwC_qval5_FC1.5$ID, sg_DwoC_DwC_qval5_FC1.5$ID))
sg_n34=length(intersect(sg_NoD_Diabetes_qval5_FC1.5$ID, sg_DwoC_DwC_qval5_FC1.5$ID))
sg_n123=length(intersect(intersect(sg_NoD_DwoC_qval5_FC1.5$ID, sg_NoD_DwC_qval5_FC1.5$ID), sg_NoD_Diabetes_qval5_FC1.5$ID))
sg_n124=length(intersect(intersect(sg_NoD_DwoC_qval5_FC1.5$ID, sg_NoD_DwC_qval5_FC1.5$ID), sg_DwoC_DwC_qval5_FC1.5$ID))
sg_n134=length(intersect(intersect(sg_NoD_DwoC_qval5_FC1.5$ID, sg_NoD_Diabetes_qval5_FC1.5$ID), sg_DwoC_DwC_qval5_FC1.5$ID))
sg_n234=length(intersect(intersect(sg_NoD_DwC_qval5_FC1.5$ID, sg_NoD_Diabetes_qval5_FC1.5$ID), sg_DwoC_DwC_qval5_FC1.5$ID))
sg_n1234=length(intersect(intersect(intersect(sg_NoD_DwoC_qval5_FC1.5$ID, sg_NoD_DwC_qval5_FC1.5$ID), sg_NoD_Diabetes_qval5_FC1.5$ID), sg_DwoC_DwC_qval5_FC1.5$ID))

plot.new()
pdf(paste(Figs,"SG_three_Group_Comparison_qval5_FC1.5.pdf",sep=""),width=12)
draw.triple.venn(area1=sg_a1, area2=sg_a2, area3=sg_a3, sg_n12, sg_n13, sg_n23, sg_n123,
                 category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
                 lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

plot.new()
pdf(paste(Figs,"SG_four_Group_Comparison_qval5_FC1.5.pdf",sep=""),width=12)
draw.quad.venn(area1=sg_a1, area2=sg_a2, area3=sg_a3, area4=sg_a4, sg_n12, sg_n13, sg_n14, sg_n23, sg_n24, sg_n34, sg_n123, sg_n124, sg_n134,
               sg_n234, sg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwoC_DwC"),
               lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
dev.off()


# # q-val <0.1 & FC > 1.5
# sg_a1<-length(sg_NoD_DwoC_qval10_FC1.5$ID)
# sg_a2<-length(sg_NoD_DwC_qval10_FC1.5$ID)
# sg_a3<-length(sg_NoD_Diabetes_qval10_FC1.5$ID)
# sg_a4<-length(sg_DwoC_DwC_qval10_FC1.5$ID)
# 
# sg_n12=length(intersect(sg_NoD_DwoC_qval10_FC1.5$ID, sg_NoD_DwC_qval10_FC1.5$ID))
# sg_n13=length(intersect(sg_NoD_DwoC_qval10_FC1.5$ID, sg_NoD_Diabetes_qval10_FC1.5$ID))
# sg_n14=length(intersect(sg_NoD_DwoC_qval10_FC1.5$ID, sg_DwoC_DwC_qval10_FC1.5$ID))
# sg_n23=length(intersect(NoD_DwC_qval10_FC1.5$ID, sg_NoD_Diabetes_qval10_FC1.5$ID))
# sg_n24=length(intersect(NoD_DwC_qval10_FC1.5$ID, sg_DwoC_DwC_qval10_FC1.5$ID))
# sg_n34=length(intersect(NoD_Diabetes_qval10_FC1.5$ID, sg_DwoC_DwC_qval10_FC1.5$ID))
# sg_n123=length(intersect(intersect(sg_NoD_DwoC_qval10_FC1.5$ID, sg_NoD_DwC_qval10_FC1.5$ID), sg_NoD_Diabetes_qval10_FC1.5$ID))
# sg_n124=length(intersect(intersect(sg_NoD_DwoC_qval10_FC1.5$ID, sg_NoD_DwC_qval10_FC1.5$ID), sg_DwoC_DwC_qval10_FC1.5$ID))
# sg_n134=length(intersect(intersect(sg_NoD_DwoC_qval10_FC1.5$ID, sg_NoD_Diabetes_qval10_FC1.5$ID), sg_DwoC_DwC_qval10_FC1.5$ID))
# sg_n234=length(intersect(intersect(NoD_DwC_qval10_FC1.5$ID, sg_NoD_Diabetes_qval10_FC1.5$ID), sg_DwoC_DwC_qval10_FC1.5$ID))
# sg_n1234=length(intersect(intersect(intersect(sg_NoD_DwoC_qval10_FC1.5$ID, sg_NoD_DwC_qval10_FC1.5$ID), sg_NoD_Diabetes_qval10_FC1.5$ID), sg_DwoC_DwC_qval10_FC1.5$ID))
# 
# plot.new()
# pdf(paste(Figs,"SG_three_Group_Comparison_qval10_FC1.5.pdf",sep=""),width=12)
# draw.triple.venn(area1=a1, area2=a2, area3=a3, n12, n13, n23, n123,
#                  category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
#                  lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
# dev.off()
# 
# plot.new()
# pdf(paste(Figs,"SG_four_Group_Comparison_qval10_FC1.5.pdf",sep=""),width=12)
# draw.quad.venn(area1=a1, area2=a2, area3=a3, area4=a4, n12, n13, n14, n23, n24, n34, n123, n124, n134,
#                n234, n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwoC_DwC"),
#                lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
# dev.off()
# 
# # q-val <0.2 & FC > 1.5
# sg_a1<-length(sg_NoD_DwoC_qval20_FC1.5$ID)
# sg_a2<-length(sg_NoD_DwC_qval20_FC1.5$ID)
# sg_a3<-length(sg_NoD_Diabetes_qval20_FC1.5$ID)
# sg_a4<-length(sg_DwoC_DwC_qval20_FC1.5$ID)
# 
# sg_n12=length(intersect(sg_NoD_DwoC_qval20_FC1.5$ID, sg_NoD_DwC_qval20_FC1.5$ID))
# sg_n13=length(intersect(sg_NoD_DwoC_qval20_FC1.5$ID, sg_NoD_Diabetes_qval20_FC1.5$ID))
# sg_n14=length(intersect(sg_NoD_DwoC_qval20_FC1.5$ID, sg_DwoC_DwC_qval20_FC1.5$ID))
# sg_n23=length(intersect(sg_NoD_DwC_qval20_FC1.5$ID, sg_NoD_Diabetes_qval20_FC1.5$ID))
# sg_n24=length(intersect(sg_NoD_DwC_qval20_FC1.5$ID, sg_DwoC_DwC_qval20_FC1.5$ID))
# sg_n34=length(intersect(sg_NoD_Diabetes_qval20_FC1.5$ID, sg_DwoC_DwC_qval20_FC1.5$ID))
# sg_n123=length(intersect(intersect(sg_NoD_DwoC_qval20_FC1.5$ID, sg_NoD_DwC_qval20_FC1.5$ID), sg_NoD_Diabetes_qval20_FC1.5$ID))
# sg_n124=length(intersect(intersect(sg_NoD_DwoC_qval20_FC1.5$ID, sg_NoD_DwC_qval20_FC1.5$ID), sg_DwoC_DwC_qval20_FC1.5$ID))
# sg_n134=length(intersect(intersect(sg_NoD_DwoC_qval20_FC1.5$ID, sg_NoD_Diabetes_qval20_FC1.5$ID), sg_DwoC_DwC_qval20_FC1.5$ID))
# sg_n234=length(intersect(intersect(sg_NoD_DwC_qval20_FC1.5$ID, sg_NoD_Diabetes_qval20_FC1.5$ID), sg_DwoC_DwC_qval20_FC1.5$ID))
# sg_n1234=length(intersect(intersect(intersect(sg_NoD_DwoC_qval20_FC1.5$ID, sg_NoD_DwC_qval20_FC1.5$ID), sg_NoD_Diabetes_qval20_FC1.5$ID), sg_DwoC_DwC_qval20_FC1.5$ID))
# 
# plot.new()
# pdf(paste(Figs,"SG_three_Group_Comparison_qval20_FC1.5.pdf",sep=""),width=12)
# draw.triple.venn(area1=sg_a1, area2=sg_a2, area3=sg_a3, sg_n12, sg_n13, sg_n23, sg_n123,
#                  category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes"),
#                  lty = "solid", fill = c("skyblue", "pink1", "mediumorchid"))
# dev.off()
# 
# plot.new()
# pdf(paste(Figs,"SG_four_Group_Comparison_qval20_FC1.5.pdf",sep=""),width=12)
# draw.quad.venn(area1=sg_a1, area2=sg_a2, area3=sg_a3, area4=sg_a4, sg_n12, sg_n13, sg_n14, sg_n23, sg_n24, 
#                sg_n34, sg_n123, sg_n124, sg_n134, sg_n234, sg_n1234, category = c("NoD_DwoC", "NoD_DwC", "NoD_Diabetes", "DwoC_DwC"),
#                lty = "solid", fill = c("skyblue", "pink1", "mediumorchid", "yellow"))
# dev.off()

rownames(hg_NoD_DwC_qval5_FC1.5)<-hg_NoD_DwC_qval5_FC1.5$ID
rownames(sg_NoD_DwC_qval5_FC1.5)<-sg_NoD_DwC_qval5_FC1.5$ID
a1<-length(hg_NoD_DwC_qval5_FC1.5$ID)
a2<-length(sg_NoD_DwC_qval5_FC1.5$ID)

n12=length(intersect(hg_NoD_DwC_qval5_FC1.5$ID, sg_NoD_DwC_qval5_FC1.5$ID))

hg_sg_NoD_DwC_qval5_FC1.5<-hg_NoD_DwC_qval5_FC1.5[intersect(hg_NoD_DwC_qval5_FC1.5$ID, sg_NoD_DwC_qval5_FC1.5$ID),]
write.csv(hg_sg_NoD_DwC_qval5_FC1.5, paste(Output,"hg_sg_NoD_DwC_qval5_FC1.5_",date,".csv",sep=""),row.names = F)


plot.new()
pdf(paste(Figs,"HG_SG_common_gene_NoD_DwC_qval5_FC1.5_qval5_FC1.5.pdf",sep=""),width=12)
draw.pairwise.venn(area1=a1, area2=a2, n12,
                 category = c("hg_NoD_DwC", "sg_NoD_DwC"),
                 lty = "solid", fill = c("skyblue", "pink"))
dev.off()

rownames(hg_NoD_DwoC_qval5_FC1.5)<-hg_NoD_DwoC_qval5_FC1.5$ID
rownames(sg_NoD_DwoC_qval5_FC1.5)<-sg_NoD_DwoC_qval5_FC1.5$ID
a1<-length(hg_NoD_DwoC_qval5_FC1.5$ID)
a2<-length(sg_NoD_DwoC_qval5_FC1.5$ID)

n12=length(intersect(hg_NoD_DwoC_qval5_FC1.5$ID, sg_NoD_DwoC_qval5_FC1.5$ID))

hg_sg_NoD_DwoC_qval5_FC1.5<-hg_NoD_DwoC_qval5_FC1.5[intersect(hg_NoD_DwoC_qval5_FC1.5$ID, sg_NoD_DwoC_qval5_FC1.5$ID),]
write.csv(hg_sg_NoD_DwoC_qval5_FC1.5, paste(Output,"hg_sg_NoD_DwoC_qval5_FC1.5_",date,".csv",sep=""),row.names = F)



plot.new()
pdf(paste(Figs,"HG_SG_common_gene_NoD_DwoC_qval5_FC1.5_qval5_FC1.5.pdf",sep=""),width=12)
draw.pairwise.venn(area1=a1, area2=a2, n12,
                   category = c("hg_NoD_DwoC", "sg_NoD_DwoC"),
                   lty = "solid", fill = c("skyblue", "pink"))
dev.off()

rownames(hg_NoD_Diabetes_qval5_FC1.5)<-hg_NoD_Diabetes_qval5_FC1.5$ID
rownames(sg_NoD_Diabetes_qval5_FC1.5)<-sg_NoD_Diabetes_qval5_FC1.5$ID
a1<-length(hg_NoD_Diabetes_qval5_FC1.5$ID)
a2<-length(sg_NoD_Diabetes_qval5_FC1.5$ID)

n12=length(intersect(hg_NoD_Diabetes_qval5_FC1.5$ID, sg_NoD_Diabetes_qval5_FC1.5$ID))

hg_sg_NoD_Diabetes_qval5_FC1.5<-hg_NoD_Diabetes_qval5_FC1.5[intersect(hg_NoD_Diabetes_qval5_FC1.5$ID, sg_NoD_Diabetes_qval5_FC1.5$ID),]
write.csv(hg_sg_NoD_Diabetes_qval5_FC1.5, paste(Output,"hg_sg_NoD_Diabetes_qval5_FC1.5_",date,".csv",sep=""),row.names = F)

plot.new()
pdf(paste(Figs,"HG_SG_common_gene_NoD_Diabetes_qval5_FC1.5_qval5_FC1.5.pdf",sep=""),width=12)
draw.pairwise.venn(area1=a1, area2=a2, n12,
                   category = c("hg_NoD_Diabetes", "sg_NoD_Diabetes"),
                   lty = "solid", fill = c("skyblue", "pink"))
dev.off()

rownames(hg_DwoC_DwC_qval5_FC1.5)<-hg_DwoC_DwC_qval5_FC1.5$ID
rownames(sg_DwoC_DwC_qval5_FC1.5)<-sg_DwoC_DwC_qval5_FC1.5$ID
a1<-length(hg_DwoC_DwC_qval5_FC1.5$ID)
a2<-length(sg_DwoC_DwC_qval5_FC1.5$ID)

n12=length(intersect(hg_DwoC_DwC_qval5_FC1.5$ID, sg_DwoC_DwC_qval5_FC1.5$ID))

hg_sg_DwoC_DwC_qval5_FC1.5<-hg_DwoC_DwC_qval5_FC1.5[intersect(hg_DwoC_DwC_qval5_FC1.5$ID, sg_DwoC_DwC_qval5_FC1.5$ID),]
write.csv(hg_sg_DwoC_DwC_qval5_FC1.5, paste(Output,"hg_sg_DwoC_DwC_qval5_FC1.5_",date,".csv",sep=""),row.names = F)


plot.new()
pdf(paste(Figs,"HG_SG_common_gene_DwoC_DwC_qval5_FC1.5_qval5_FC1.5.pdf",sep=""),width=12)
draw.pairwise.venn(area1=a1, area2=a2, n12,
                   category = c("hg_DwoC_DwC", "sg_DwoC_DwC"),
                   lty = "solid", fill = c("skyblue", "pink"))
dev.off()
