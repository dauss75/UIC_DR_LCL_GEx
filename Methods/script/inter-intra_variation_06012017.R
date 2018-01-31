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

# parallization
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

dat<-local(get(load(file=paste(RData,"normalizedDataMatrix.RData",sep=""))))
# dat<-local(get(load(file=paste(RData,"normalizedDataMatrix_filtered.RData",sep=""))))
gene_annotation<-read.csv(paste(PhenotypeDir,"raw_delta_gene_anno.csv",sep=""))
rownames(gene_annotation)<-gene_annotation$GeneCode

# for normal samples
NoD<-dat[,grep("NoD",colnames(dat))] # grep all normal samples
DwC<-dat[,grep("DwC",colnames(dat))] # grep all DwC samples
DwoC<-dat[,grep("DwoC",colnames(dat))] # grep all DwoC samples

# separate high glucose and standard glucose
NoD_norm<-NoD[,grep("norm",colnames(NoD))]; NoD_30mM<-NoD[,grep("30mM",colnames(NoD))]
NoD_DE<-NoD_30mM-NoD_norm; colnames(NoD_DE)<-gsub("_30mM","",colnames(NoD_DE))

DwC_norm<-DwC[,grep("norm",colnames(DwC))]; DwC_30mM<-DwC[,grep("30mM",colnames(DwC))]
DwC_DE<-DwC_30mM-DwC_norm; colnames(DwC_DE)<-gsub("_30mM","",colnames(DwC_DE))

DwoC_norm<-DwoC[,grep("norm",colnames(DwoC))]; DwoC_30mM<-DwoC[,grep("30mM",colnames(DwoC))]
DwoC_DE<-DwoC_30mM-DwoC_norm; colnames(DwoC_DE)<-gsub("_30mM","",colnames(DwoC_DE))

sg<-dat[,grep("norm",colnames(dat))]; hg<-dat[,grep("30mM",colnames(dat))]
All_DE<-hg-sg; colnames(All_DE)<-gsub("_30mM","",colnames(All_DE))
Diabetes<-cbind(DwC_DE,DwoC_DE)

#------normal samples with standard glucose (sg)
NoD_sg_output<-var_proc(NoD_norm,gene_annotation)
write.csv(NoD_sg_output, paste(Output,"sg_NoD_variation_",date,".csv",sep=""))
NoD_sg_output[,2:ncol(NoD_sg_output)]<-data.frame(apply(NoD_sg_output[,2:ncol(NoD_sg_output)],2,as.numeric)); 

inter<-log2(NoD_sg_output$`inter-individual variance`)
intra<-log2(NoD_sg_output$`intra-individual variance`)
data<-as.data.frame(cbind(inter,intra))

tmp<-cbind(log2(NoD_sg_output$`inter-individual variance`),log2(NoD_sg_output$`intra-individual variance`))
rownames(tmp)<-rownames(NoD_sg_output); colnames(tmp)<-c('inter-individual variance', 'intra-individual variance')

setEPS()
postscript(file = paste(Figs,"sg_NoD_inter_vs_intra-individual_",date,".eps",sep=""))

ggplot(tmp, aes(y = `inter-individual variance`, x = `intra-individual variance`)) + 
  geom_point(shape=1,size=0.1) +theme(text = element_text(size=20)) + theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-10, 7.5) + ylim(-10, 7.5)+geom_abline(intercept=0, slope=1)+xlab("log2(intra-individual variance)") + 
  ylab("log2(inter-individual variance)") + ggtitle("No diabetes (standard glucose)")
dev.off()

setEPS()
postscript(file = paste(Figs,"sg_NoD_inter_vs_intra-individual_density_",date,".eps",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="No diabetes (standard glucose)")
dev.off()

setEPS()
postscript(file = paste(Figs,"sg_NoD_inter_vs_intra-individual_QQplot_",date,".eps",sep=""))
ggd.qqplot(NoD_sg_output$`p-value`,"No diabetes (standard glucose)")
dev.off()

NoD_sg_output_pval<-NoD_sg_output$`p-value`
#------normal samples with high glucose (hg)
NoD_hg_output<-var_proc(NoD_30mM,gene_annotation)

write.csv(NoD_hg_output, paste(Output,"hg_NoD_glucose_variation_",date,".csv",sep=""))
NoD_hg_output[,2:ncol(NoD_hg_output)]<-data.frame(apply(NoD_hg_output[,2:ncol(NoD_hg_output)],2,as.numeric)); 

inter<-log2(NoD_hg_output$`inter-individual variance`)
intra<-log2(NoD_hg_output$`intra-individual variance`)
data<-as.data.frame(cbind(inter,intra))

tmp<-cbind(log2(NoD_hg_output$`inter-individual variance`),log2(NoD_hg_output$`intra-individual variance`))
rownames(tmp)<-rownames(NoD_hg_output); colnames(tmp)<-c('inter-individual variance', 'intra-individual variance')

setEPS()
postscript(file = paste(Figs,"hg_NoD_inter_vs_intra-individual_",date,".eps",sep=""))
ggplot(tmp, aes(y = `inter-individual variance`, x = `intra-individual variance`)) + 
  geom_point(shape=1,size=0.1) +theme(text = element_text(size=20)) + theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-10, 7.5) + ylim(-10, 7.5)+geom_abline(intercept=0, slope=1)+xlab("log2(intra-individual variance)") + 
  ylab("log2(inter-individual variance)") + ggtitle("No diabetes (high glucose)")
dev.off()

setEPS()
postscript(file = paste(Figs,"hg_NoD_inter_vs_intra-individual_density_",date,".eps",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="No diabetes (high glucose)")
dev.off()

setEPS()
postscript(file = paste(Figs,"hg_NoD_inter_vs_intra-individual_QQplot_",date,".eps",sep=""))
ggd.qqplot(NoD_hg_output$`p-value`,"No diabetes (high glucose)")
dev.off()

NoD_hg_output_pval<-NoD_hg_output$`p-value`

# find outliers
NoD_sg_outlier<-NoD_sg_output[which(NoD_sg_output$`intra-individual variance` >NoD_sg_output$`inter-individual variance`),]
NoD_hg_outlier<-NoD_hg_output[which(NoD_hg_output$`intra-individual variance` >NoD_hg_output$`inter-individual variance`),]

NoD_common_outlier<-intersect(rownames(NoD_sg_outlier),rownames(NoD_hg_outlier))
nrow(NoD_sg_outlier);nrow(NoD_hg_outlier);length(NoD_common_outlier)

#------normal samples delta
NoD_delta<-NoD_30mM-NoD_norm
NoD_delta_output<-var_proc(NoD_delta,gene_annotation)

write.csv(NoD_delta_output, paste(Output,"delta_NoD_variation_",date,".csv",sep=""))
NoD_delta_output[,2:ncol(NoD_delta_output)]<-data.frame(apply(NoD_delta_output[,2:ncol(NoD_delta_output)],2,as.numeric)); 

inter<-log2(NoD_delta_output$`inter-individual variance`)
intra<-log2(NoD_delta_output$`intra-individual variance`)
data<-as.data.frame(cbind(inter,intra))

tmp<-cbind(log2(NoD_delta_output$`inter-individual variance`),log2(NoD_delta_output$`intra-individual variance`))
rownames(tmp)<-rownames(NoD_delta_output); colnames(tmp)<-c('inter-individual variance', 'intra-individual variance')

setEPS()
postscript(file = paste(Figs,"delta_NoD_inter_vs_intra-individual_",date,".eps",sep=""))
ggplot(tmp, aes(y = `inter-individual variance`, x = `intra-individual variance`)) + 
  geom_point(shape=1,size=0.1) +theme(text = element_text(size=20)) + theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-10, 7.5) + ylim(-10, 7.5)+geom_abline(intercept=0, slope=1)+xlab("log2(intra-individual variance)") + 
  ylab("log2(inter-individual variance)") + ggtitle("No diabetes (high vs. standard glucose)")
dev.off()

setEPS()
postscript(file = paste(Figs,"delta_NoD_inter_vs_intra-individual_density_",date,".eps",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="No diabetes (high vs. standard glucose)")
dev.off()

setEPS()
postscript(file = paste(Figs,"delta_NoD_inter_vs_intra-individual_QQplot_",date,".eps",sep=""))
ggd.qqplot(NoD_delta_output$`p-value`,"No diabetes (high vs. standard glucose)")
dev.off()

NoD_delta_output_pval<-NoD_delta_output$`p-value`

# find outliers
NoD_delta_outlier<-NoD_delta_output[which(NoD_delta_output$`intra-individual variance` >NoD_delta_output$`inter-individual variance`),]
dim(NoD_delta_outlier)

#------diabets with complication samples with standard glucose (sg)
DwC_sg_output<-var_proc(DwC_norm,gene_annotation)
write.csv(DwC_sg_output, paste(Output,"sg_DwC_variation_",date,".csv",sep=""))
DwC_sg_output[,2:ncol(DwC_sg_output)]<-data.frame(apply(DwC_sg_output[,2:ncol(DwC_sg_output)],2,as.numeric)); 

inter<-log2(DwC_sg_output$`inter-individual variance`)
intra<-log2(DwC_sg_output$`intra-individual variance`)
data<-as.data.frame(cbind(inter,intra))

setEPS()
postscript(file = paste(Figs,"sg_DwC_inter_vs_intra-individual_density_",date,".eps",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="Proliferative diabetic retinopathy (standard glucose)")
dev.off()

tmp<-cbind(log2(DwC_sg_output$`inter-individual variance`),log2(DwC_sg_output$`intra-individual variance`))
rownames(tmp)<-rownames(DwC_sg_output); colnames(tmp)<-c('inter-individual variance', 'intra-individual variance')

setEPS()
postscript(file = paste(Figs,"sg_DwC_inter_vs_intra-individual_",date,".eps",sep=""))
ggplot(tmp, aes(y = `inter-individual variance`, x = `intra-individual variance`)) + 
  geom_point(shape=1,size=0.1) +theme(text = element_text(size=20)) + theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-10, 7.5) + ylim(-10, 7.5)+geom_abline(intercept=0, slope=1)+xlab("log2(intra-individual variance)") + 
  ylab("log2(inter-individual variance)") + ggtitle("PDR (standard glucose)")
dev.off()

setEPS()
postscript(file = paste(Figs,"sg_DwC_inter_vs_intra-individual_QQplot_",date,".eps",sep=""))
ggd.qqplot(DwC_sg_output$`p-value`,"PDR (standard glucose)")
dev.off()

DwC_sg_output_pval<-DwC_sg_output$`p-value`

#------diabets with complication samples with high glucose (hg)
DwC_hg_output<-var_proc(DwC_30mM,gene_annotation)
write.csv(DwC_hg_output, paste(Output,"hg_DwC_variation_",date,".csv",sep=""))
DwC_hg_output[,2:ncol(DwC_hg_output)]<-data.frame(apply(DwC_hg_output[,2:ncol(DwC_hg_output)],2,as.numeric)); 

inter<-log2(DwC_hg_output$`inter-individual variance`)
intra<-log2(DwC_hg_output$`intra-individual variance`)
data<-as.data.frame(cbind(inter,intra))


setEPS()
postscript(file = paste(Figs,"hg_DwC_inter_vs_intra-individual_density_",date,".eps",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="Proliferative diabetic retinopathy  (high glucose)")
dev.off()

tmp<-cbind(log2(DwC_hg_output$`inter-individual variance`),log2(DwC_hg_output$`intra-individual variance`))
rownames(tmp)<-rownames(DwC_hg_output); colnames(tmp)<-c('inter-individual variance', 'intra-individual variance')

setEPS()
postscript(file = paste(Figs,"hg_DwC_inter_vs_intra-individual_",date,".eps",sep=""))
ggplot(tmp, aes(y = `inter-individual variance`, x = `intra-individual variance`)) + 
  geom_point(shape=1,size=0.1) +theme(text = element_text(size=20)) + theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-10, 7.5) + ylim(-10, 7.5)+geom_abline(intercept=0, slope=1)+xlab("log2(intra-individual variance)") + 
  ylab("log2(inter-individual variance)") + ggtitle("PDR (high glucose)")
dev.off()

setEPS()
postscript(file = paste(Figs,"hg_DwC_inter_vs_intra-individual_QQplot_",date,".eps",sep=""))
ggd.qqplot(DwC_hg_output$`p-value`,"PDR (high glucose)")
dev.off()

DwC_hg_output_pval<-DwC_hg_output$`p-value`

# find outliers
DwC_sg_outlier<-DwC_sg_output[which(DwC_sg_output$`intra-individual variance` >DwC_sg_output$`inter-individual variance`),]
DwC_hg_outlier<-DwC_hg_output[which(DwC_hg_output$`intra-individual variance` >DwC_hg_output$`inter-individual variance`),]

DwC_common_outlier<-intersect(rownames(DwC_sg_outlier),rownames(DwC_hg_outlier))
nrow(DwC_sg_outlier);nrow(DwC_hg_outlier);length(DwC_common_outlier)

#------diabets with complication samples delta
DwC_delta<-DwC_30mM-DwC_norm
DwC_delta_output<-var_proc(DwC_delta,gene_annotation)

write.csv(DwC_delta_output, paste(Output,"delta_DwC_variation_",date,".csv",sep=""))
DwC_delta_output[,2:ncol(DwC_delta_output)]<-data.frame(apply(DwC_delta_output[,2:ncol(DwC_delta_output)],2,as.numeric)); 

inter<-log2(DwC_delta_output$`inter-individual variance`)
intra<-log2(DwC_delta_output$`intra-individual variance`)
data<-as.data.frame(cbind(inter,intra))


setEPS()
postscript(file = paste(Figs,"delta_DwC_inter_vs_intra-individual_density_",date,".eps",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="Proliferative diabetic retinopathy (high vs. standard glucose)")
dev.off()

tmp<-cbind(log2(DwC_delta_output$`inter-individual variance`),log2(DwC_delta_output$`intra-individual variance`))
rownames(tmp)<-rownames(DwC_delta_output); colnames(tmp)<-c('inter-individual variance', 'intra-individual variance')
setEPS()
postscript(file = paste(Figs,"delta_DwC_inter_vs_intra-individual_",date,".eps",sep=""))
ggplot(tmp, aes(y = `inter-individual variance`, x = `intra-individual variance`)) + 
  geom_point(shape=1,size=0.1) +theme(text = element_text(size=20)) + theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-10, 7.5) + ylim(-10, 7.5)+geom_abline(intercept=0, slope=1)+xlab("log2(intra-individual variance)") + 
  ylab("log2(inter-individual variance)") + ggtitle("PDR (high vs. standard glucose)")
dev.off()

setEPS()
postscript(file = paste(Figs,"delta_DwC_inter_vs_intra-individual_QQplot_",date,".eps",sep=""))
ggd.qqplot(DwC_delta_output$`p-value`,"PDR (high vs. standard glucose)")
dev.off()

DwC_delta_output_pval<-DwC_delta_output$`p-value`

# find outliers
DwC_delta_outlier<-DwC_delta_output[which(DwC_delta_output$`intra-individual variance` >DwC_delta_output$`inter-individual variance`),]
dim(DwC_delta_outlier)

#------diabets without complication samples with standard glucose (sg)
DwoC_sg_output<-var_proc(DwoC_norm,gene_annotation)
write.csv(DwoC_sg_output, paste(Output,"sg_DwoC_variation_",date,".csv",sep=""))
DwoC_sg_output[,2:ncol(DwoC_sg_output)]<-data.frame(apply(DwoC_sg_output[,2:ncol(DwoC_sg_output)],2,as.numeric)); 

inter<-log2(DwoC_sg_output$`inter-individual variance`)
intra<-log2(DwoC_sg_output$`intra-individual variance`)
data<-as.data.frame(cbind(inter,intra))

setEPS()
postscript(file = paste(Figs,"sg_DwoC_inter_vs_intra-individual_density_",date,".eps",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="Diabetes without retinopathy (standard glucose)")
dev.off()

tmp<-cbind(log2(DwoC_sg_output$`inter-individual variance`),log2(DwoC_sg_output$`intra-individual variance`))
rownames(tmp)<-rownames(DwoC_sg_output); colnames(tmp)<-c('inter-individual variance', 'intra-individual variance')

setEPS()
postscript(file = paste(Figs,"sg_DwoC_inter_vs_intra-individual_",date,".eps",sep=""))
ggplot(tmp, aes(y = `inter-individual variance`, x = `intra-individual variance`)) + 
  geom_point(shape=1,size=0.1) +theme(text = element_text(size=20)) + theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-10, 7.5) + ylim(-10, 7.5)+geom_abline(intercept=0, slope=1)+xlab("log2(intra-individual variance)") + 
  ylab("log2(inter-individual variance)") + ggtitle("No_PDR (standard glucose)")
dev.off()

setEPS()
postscript(file = paste(Figs,"sg_DwoC_inter_vs_intra-individual_QQplot_",date,".eps",sep=""))
ggd.qqplot(DwoC_sg_output$`p-value`,"No_PDR (standard glucose)")
dev.off()

DwoC_sg_output_pval<-DwoC_sg_output$`p-value`

#------diabets without complication samples with high glucose (hg)
DwoC_hg_output<-var_proc(DwoC_30mM,gene_annotation)
write.csv(DwoC_hg_output, paste(Output,"hg_DwoC_variation_",date,".csv",sep=""))
DwoC_hg_output[,2:ncol(DwoC_hg_output)]<-data.frame(apply(DwoC_hg_output[,2:ncol(DwoC_hg_output)],2,as.numeric)); 

inter<-log2(DwoC_hg_output$`inter-individual variance`)
intra<-log2(DwoC_hg_output$`intra-individual variance`)
data<-as.data.frame(cbind(inter,intra))

setEPS()
postscript(file = paste(Figs,"hg_DwoC_inter_vs_intra-individual_density_",date,".eps",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="Diabetes without retinopathy (high glucose)")
dev.off()

tmp<-cbind(log2(DwoC_hg_output$`inter-individual variance`),log2(DwoC_hg_output$`intra-individual variance`))
rownames(tmp)<-rownames(DwoC_sg_output); colnames(tmp)<-c('inter-individual variance', 'intra-individual variance')

setEPS()
postscript(file = paste(Figs,"hg_DwoC_inter_vs_intra-individual_",date,".eps",sep=""))
ggplot(tmp, aes(y = `inter-individual variance`, x = `intra-individual variance`)) + 
  geom_point(shape=1,size=0.1) +theme(text = element_text(size=20)) + theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-10, 7.5) + ylim(-10, 7.5)+geom_abline(intercept=0, slope=1)+xlab("log2(intra-individual variance)") + 
  ylab("log2(inter-individual variance)") + ggtitle("No_PDR (high glucose)")
dev.off()

setEPS()
postscript(file = paste(Figs,"hg_DwoC_inter_vs_intra-individual_QQplot_",date,".eps",sep=""))
ggd.qqplot(DwoC_hg_output$`p-value`,"No_PDR (high glucose)")
dev.off()

DwoC_hg_output_pval<-DwoC_hg_output$`p-value`

# find outliers
DwoC_sg_outlier<-DwoC_sg_output[which(DwoC_sg_output$`intra-individual variance` >DwoC_sg_output$`inter-individual variance`),]
DwoC_hg_outlier<-DwoC_hg_output[which(DwoC_hg_output$`intra-individual variance` >DwoC_hg_output$`inter-individual variance`),]

DwoC_common_outlier<-intersect(rownames(DwoC_sg_outlier),rownames(DwoC_hg_outlier))

nrow(DwoC_sg_outlier);nrow(DwoC_hg_outlier);length(DwoC_common_outlier)

#------diabets without complication samples delta
DwoC_delta<-DwoC_30mM-DwoC_norm
DwoC_delta_output<-var_proc(DwoC_delta,gene_annotation)

write.csv(DwoC_delta_output, paste(Output,"delta_DwoC_variation_",date,".csv",sep=""))
DwoC_delta_output[,2:ncol(DwoC_delta_output)]<-data.frame(apply(DwoC_delta_output[,2:ncol(DwoC_delta_output)],2,as.numeric)); 

inter<-log2(DwoC_delta_output$`inter-individual variance`)
intra<-log2(DwoC_delta_output$`intra-individual variance`)
data<-as.data.frame(cbind(inter,intra))

setEPS()
postscript(file = paste(Figs,"delta_DwoC_inter_vs_intra-individual_density_",date,".eps",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="Diabetes without retinopathy  (high vs. standard glucose)")
dev.off()

tmp<-cbind(log2(DwoC_delta_output$`inter-individual variance`),log2(DwoC_delta_output$`intra-individual variance`))
rownames(tmp)<-rownames(DwoC_sg_output); colnames(tmp)<-c('inter-individual variance', 'intra-individual variance')

setEPS()
postscript(file = paste(Figs,"delta_DwoC_inter_vs_intra-individual_",date,".eps",sep=""))
ggplot(tmp, aes(y = `inter-individual variance`, x = `intra-individual variance`)) + 
  geom_point(shape=1,size=0.1) +theme(text = element_text(size=20)) + theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-10, 7.5) + ylim(-10, 7.5)+geom_abline(intercept=0, slope=1)+xlab("log2(intra-individual variance)") + 
  ylab("log2(inter-individual variance)") + ggtitle("No_PDR (high vs. standard glucose)")
dev.off()

setEPS()
postscript(file = paste(Figs,"delta_DwoC_inter_vs_intra-individual_QQplot_",date,".eps",sep=""))
ggd.qqplot(DwoC_delta_output$`p-value`,"No_PDR (high vs. standard glucose)")
dev.off()

DwoC_delta_output_pval<-DwoC_delta_output$`p-value`

# find outliers
DwoC_delta_outlier<-DwoC_delta_output[which(DwoC_delta_output$`intra-individual variance` >DwoC_delta_output$`inter-individual variance`),]
dim(DwoC_delta_outlier)

tmp_outlier<-intersect(NoD_common_outlier,DwC_common_outlier); 
common_outlier<-intersect(tmp_outlier,DwoC_common_outlier)  # 23 outliers
save(common_outlier, file = paste(RData,"common_outlier.RData",sep=""))

dat2<-dat[-which(rownames(dat) %in% common_outlier),]
save(dat2,file = paste(RData,"normalizedDataMatrix_filtered.RData",sep=""))
#---------------------------------------------------------------------------------
##------inter and intra-individual variability for standard glucose, high glucose, all-together
#------all standard glucose

sg_output<-var_proc(sg,gene_annotation)
sg_output<-sg_output[order(sg_output$`p-value`,decreasing = TRUE),]
write.csv(sg_output, paste(Output,"sg_all_inter_intra_individual_variation_",date,".csv",sep=""))
sg_output[,2:ncol(sg_output)]<-data.frame(apply(sg_output[,2:ncol(sg_output)],2,as.numeric));

inter<-log2(sg_output$`inter-individual variance`)
intra<-log2(sg_output$`intra-individual variance`)
data<-as.data.frame(cbind(inter,intra))

setEPS()
postscript(file = paste(Figs,"sg_all_inter_vs_intra-individual_density_",date,".eps",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="standard glucose samples inter and intra-individual variability")
dev.off()

tmp<-cbind(log2(sg_output$`inter-individual variance`),log2(sg_output$`intra-individual variance`))
rownames(tmp)<-rownames(sg_output); colnames(tmp)<-c('inter-individual variance', 'intra-individual variance')


setEPS()
postscript(file = paste(Figs,"sg_all_inter_vs_intra-individual_",date,".eps",sep=""))
ggplot(tmp, aes(y = `inter-individual variance`, x = `intra-individual variance`)) + 
  geom_point(shape=1,size=0.1) +theme(text = element_text(size=20)) + theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-10, 7.5) + ylim(-10, 7.5)+geom_abline(intercept=0, slope=1)+xlab("log2(intra-individual variance)") + 
  ylab("log2(inter-individual variance)") + ggtitle("standard glucose")
dev.off()

setEPS()
postscript(file = paste(Figs,"sg_all_inter_vs_intra-individual_QQplot_",date,".eps",sep=""))
ggd.qqplot(sg_output$`p-value`,"standard glucose")
dev.off()

sg_output_pval<-sg_output$`p-value`

#------all high glucose
hg_output<-var_proc(hg,gene_annotation)
hg_output<-hg_output[order(hg_output$`p-value`,decreasing = TRUE),]
write.csv(hg_output, paste(Output,"hg_all_inter_intra_individual_variation_",date,".csv",sep=""))

hg_output[,2:ncol(hg_output)]<-data.frame(apply(hg_output[,2:ncol(hg_output)],2,as.numeric));

inter<-log2(hg_output$`inter-individual variance`)
intra<-log2(hg_output$`intra-individual variance`)
data<-as.data.frame(cbind(inter,intra))

setEPS()
postscript(file = paste(Figs,"hg_all_inter_vs_intra-individual_density_",date,".eps",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="high glucose samples inter and intra-individual variability")
dev.off()

tmp<-cbind(log2(hg_output$`inter-individual variance`),log2(hg_output$`intra-individual variance`))
rownames(tmp)<-rownames(hg_output); colnames(tmp)<-c('inter-individual variance', 'intra-individual variance')


setEPS()
postscript(file = paste(Figs,"hg_all_inter_vs_intra-individual_",date,".eps",sep=""))
ggplot(tmp, aes(y = `inter-individual variance`, x = `intra-individual variance`)) + 
  geom_point(shape=1,size=0.1) +theme(text = element_text(size=20)) + theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-10, 7.5) + ylim(-10, 7.5)+geom_abline(intercept=0, slope=1)+xlab("log2(intra-individual variance)") + 
  ylab("log2(inter-individual variance)") + ggtitle("high glucose")
dev.off()

setEPS()
postscript(file = paste(Figs,"hg_all_inter_vs_intra-individual_QQplot_",date,".eps",sep=""))
ggd.qqplot(hg_output$`p-value`,"high glucose")
dev.off()

hg_output_pval<-hg_output$`p-value`

#------standard + high glucose
all_output<-var_proc(dat,gene_annotation)
all_output<-all_output[order(all_output$`p-value`,decreasing = TRUE),]
write.csv(all_output, paste(Output,"all(hg+sg)_samples_inter_intra_individual_variation_",date,".csv",sep=""))

all_output[,2:ncol(all_output)]<-data.frame(apply(all_output[,2:ncol(all_output)],2,as.numeric));

inter<-log2(all_output$`inter-individual variance`)
intra<-log2(all_output$`intra-individual variance`)
data<-as.data.frame(cbind(inter,intra))

setEPS()
postscript(file = paste(Figs,"all(hg+sg)_inter_vs_intra-individual_density_",date,".eps",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="all samples inter and intra-individual variability")
dev.off()

tmp<-cbind(log2(all_output$`inter-individual variance`),log2(all_output$`intra-individual variance`))
rownames(tmp)<-rownames(all_output); colnames(tmp)<-c('inter-individual variance', 'intra-individual variance')


setEPS()
postscript(file = paste(Figs,"all(hg+sg)_inter_vs_intra-individual_",date,".eps",sep=""))
ggplot(tmp, aes(y = `inter-individual variance`, x = `intra-individual variance`)) + 
  geom_point(shape=1,size=0.1) +theme(text = element_text(size=20)) + theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-10, 7.5) + ylim(-10, 7.5)+geom_abline(intercept=0, slope=1)+xlab("log2(intra-individual variance)") + 
  ylab("log2(inter-individual variance)") + ggtitle("all (high + standard) samples")
dev.off()

setEPS()
postscript(file = paste(Figs,"all(hg+sg)_inter_vs_intra-individual_QQplot_",date,".eps",sep=""))
ggd.qqplot(all_output$`p-value`,"all (high and standard glucose) samples")
dev.off()

all_output_pval<-all_output$`p-value`


# delta
de<-hg-sg
de_output<-var_proc(de,gene_annotation)
de_output<-de_output[order(de_output$`p-value`,decreasing = TRUE),]
write.csv(de_output, paste(Output,"delta_all_samples_inter_intra_individual_variation_",date,".csv",sep=""))

de_output[,2:ncol(de_output)]<-data.frame(apply(de_output[,2:ncol(de_output)],2,as.numeric));

inter<-log2(de_output$`inter-individual variance`)
intra<-log2(de_output$`intra-individual variance`)
data<-as.data.frame(cbind(inter,intra))

setEPS()
postscript(file = paste(Figs,"delta_all_inter_vs_intra-individual_density_",date,".eps",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="delta all samples inter and intra-individual variability")
dev.off()

tmp<-cbind(log2(de_output$`inter-individual variance`),log2(de_output$`intra-individual variance`))
rownames(tmp)<-rownames(de_output); colnames(tmp)<-c('inter-individual variance', 'intra-individual variance')


setEPS()
postscript(file = paste(Figs,"delta_all_inter_vs_intra-individual_",date,".eps",sep=""))
ggplot(tmp, aes(y = `inter-individual variance`, x = `intra-individual variance`)) + 
  geom_point(shape=1,size=0.1) +theme(text = element_text(size=20)) + theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-10, 7.5) + ylim(-10, 7.5)+geom_abline(intercept=0, slope=1)+xlab("log2(intra-individual variance)") + 
  ylab("log2(inter-individual variance)") + ggtitle("delta all samples")
dev.off()

setEPS()
postscript(file = paste(Figs,"delta_all_inter_vs_intra-individual_QQplot_",date,".eps",sep=""))
ggd.qqplot(de_output$`p-value`,"delta (high vs standard glucose)")
dev.off()

de_output_pval<-de_output$`p-value`

stopCluster(cl)

rm(NoD,DwC,DwoC,NoD_norm,NoD_30mM,DwC_norm,DwC_30mM,DwoC_norm,DwoC_30mM,sg,hg,NoD_sg_output,
   inter,intra,data,NoD_hg_output,DwC_hg_output,DwoC_hg_output,DwoC_sg_output,hg_output,sg_output,DwC_sg_output)


#---------------------
# p-val distribution
pval<-cbind(NoD_sg_output_pval,NoD_hg_output_pval,DwC_sg_output_pval,DwC_hg_output_pval,
            DwoC_sg_output_pval,DwoC_hg_output_pval,sg_output_pval,hg_output_pval,all_output_pval)
colnames(pval)<-c("NoD_sg","NoD_hg","DwC_sg","DwC_hg","DwoC_sg","DwoC_hg","sg_all","hg_all","all(hg+sg)")

setEPS()
postscript(paste(Figs,"inter_intra_variation_pvalue_0.01_distribution.eps",sep=""))
par(mfrow=c(3,3))
for (i in 1:ncol(pval)) {
  hist(pval[,i], main=paste(colnames(pval)[i],": P>0.01=",sum(pval[,i]>0.01), "(", signif(sum(pval[,i]>0.01)/length(pval[,i])*100,digits=3),"%)",sep=""))
}
dev.off()

setEPS()
postscript(file = paste(Figs,"inter_intra_variation_pvalue_0.05_distribution.eps",sep=""))
par(mfrow=c(3,3))
for (i in 1:ncol(pval)) {
  hist(pval[,i], main=paste(colnames(pval)[i],": P>0.05=",sum(pval[,i]>0.05), "(", signif(sum(pval[,i]>0.05)/length(pval[,i])*100,digits=3),"%)",sep=""))
}
dev.off()

setEPS()
postscript(file = paste(Figs,"inter_intra_variation_pvalue_0.1_distribution.eps",sep=""))
par(mfrow=c(3,3))
for (i in 1:ncol(pval)) {
  hist(pval[,i], main=paste(colnames(pval)[i],": P>0.1=",sum(pval[,i]>0.1), " (", signif(sum(pval[,i]>0.1)/length(pval[,i])*100,digits=3),"%)",sep=""))
}
dev.off()

# #--------------------------------------------------------------------------------------
# # exploring PCA analysis on all samples (132) which include sg_normal, sg_DwC, sg_DwoC
# # and hg_normal, hg_DwC, hg_DwoC
# #--------------------------------------------------------------------------------------
# 
# #----> standard glucose
# pca_tmp <- prcomp(sg)
# title<-deparse(substitute(sg))
# pdf(file = paste(Figs,title,"_pca_variance_explained_",date,".pdf",sep=""))
# plot_variance_explained(pca_tmp) + theme_bw(base_size = 18)
# dev.off()
# 
# data_explore_by_pca(sg)
# 
# #----> high glucose
# pca_tmp <- prcomp(hg)
# title<-deparse(substitute(hg))
# pdf(file = paste(Figs,title,"_pca_variance_explained_",date,".pdf",sep=""))
# plot_variance_explained(pca_tmp) + theme_bw(base_size = 18)
# dev.off()
# 
# data_explore_by_pca(hg)
# 
# #----> standard + high glucose
# pca_tmp <- prcomp(dat)
# title<-deparse(substitute(dat))
# pdf(file = paste(Figs,title,"_pca_variance_explained_",date,".pdf",sep=""))
# plot_variance_explained(pca_tmp) + theme_bw(base_size = 18)
# dev.off()
# 
# data_explore_by_pca(dat)
# 
# #------ PCA and Covariates 1: sex including normal
# # meta<-read.table(paste(PhenotypeDir, "meanDeltaSampleInfo-age-sex.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
# # genes<-5000
# # 
# # hg.T<-data.frame(t(hg)); 
# # hg.T$Name<-substr(rownames(hg.T),1,nchar(rownames(hg.T))-2)
# # hgagg<-aggregate(hg.T[,-ncol(hg.T)],list(hg.T$Name), mean)
# # rownames(hgagg)<-hgagg$Group.1; hgagg$Group.1<-NULL
# # hgmean<-apply(hgagg,2,as.numeric); rownames(hgmean)<-rownames(hgagg)
# # hg_mean<-data.frame(t(hgmean)); colnames(hg_mean)<-substr(colnames(hg_mean),1,nchar(colnames(hg_mean))-5)
# # hg_mean$CV<-CV(hg_mean); hg_mean$mean<-rowMeans(hg_mean[,-ncol(hg_mean)])
# # hg_mean_sorted<-hg_mean[order(hg_mean$mean,decreasing = TRUE),]
# # hg_input<-hg_mean_sorted[,rownames(meta)]
# # hg_input_5000<-hg_mean_sorted[1:genes,rownames(meta)]
# # 
# # sg.T<-data.frame(t(sg)); 
# # sg.T$Name<-substr(rownames(sg.T),1,nchar(rownames(sg.T))-2)
# # sgagg<-aggregate(sg.T[,-ncol(sg.T)],list(sg.T$Name), mean)
# # rownames(sgagg)<-sgagg$Group.1; sgagg$Group.1<-NULL
# # sgmean<-apply(sgagg,2,as.numeric); rownames(sgmean)<-rownames(sgagg)
# # sg_mean<-data.frame(t(sgmean)); colnames(sg_mean)<-substr(colnames(sg_mean),1,nchar(colnames(sg_mean))-5)
# # sg_mean$CV<-CV(sg_mean); sg_mean$mean<-rowMeans(sg_mean[,-ncol(sg_mean)])
# # sg_mean_sorted<-sg_mean[order(sg_mean$mean,decreasing = TRUE),]
# # sg_input<-sg_mean_sorted[,rownames(meta)]
# # sg_input_5000<-sg_mean_sorted[1:genes,rownames(meta)]
# # 
# # All_DE.T<-data.frame(t(All_DE))
# # All_DE.T$Name<-substr(rownames(All_DE.T),1,nchar(rownames(All_DE.T))-2)
# # All_DEagg<-aggregate(All_DE.T[,-ncol(All_DE.T)],list(All_DE.T$Name), mean)
# # rownames(All_DEagg)<-All_DEagg$Group.1; All_DEagg$Group.1<-NULL
# # meandelta<-apply(All_DEagg,2,as.numeric); rownames(meandelta)<-rownames(All_DEagg)
# # delta_mean<-data.frame(t(meandelta)); delta_mean$CV<-CV(delta_mean); delta_mean$mean<-rowMeans(delta_mean[,-ncol(delta_mean)]) 
# # delta_mean_sorted<-delta_mean[order(delta_mean$mean,decreasing = TRUE),]
# # delta_input<-delta_mean_sorted[,rownames(meta)]
# # delta_input_5000<-delta_mean_sorted[1:genes,rownames(meta)]
# # 
# # rm(hg_mean, sg_mean, delta_mean)
# # 
# # subject1="sg_mean";subject2="hg_mean";subject3="delta_mean"
# # 
# # 
# # #-----------------------------------------------
# # ## standard glucose
# # #-----------------------------------------------
# # pca_sg_mean <- prcomp(sg_input, scale = FALSE, center = TRUE)
# # pdf(file = paste(Figs,subject1,"_pca_variance_explained_",date,".pdf",sep=""))
# # plot_variance_explained(pca_sg_mean) + theme_bw(base_size = 18)
# # dev.off()
# # 
# # pca_sg_mean_5000 <- prcomp(sg_input_5000, scale = FALSE, center = TRUE)
# # pdf(file = paste(Figs,subject1,"_pca_variance_explained_",genes,"genes_",date,".pdf",sep=""))
# # plot_variance_explained(pca_sg_mean_5000) + theme_bw(base_size = 18)
# # dev.off()
# # 
# # # data_explore_by_pca2(sg_input)
# # # data_explore_by_pca2(sg_input_5000)
# # 
# # # correlate pcs and proportional pcs (edit the function in pca.R)
# # sg_cor_pcs1<-correlate_pcs(pca_sg_mean, meta$SEX, npcs = 10, min.cor = 0)
# # sg_cor_pcs1<-as.data.frame(sg_cor_pcs1); sg_cor_pcs1<-abs(sg_cor_pcs1)
# # 
# # tt<-cbind(rownames(sg_cor_pcs1),sg_cor_pcs1); 
# # tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "Correlation")
# # pdf(file = paste(Figs,subject1,"_pca_covariate_",date,".pdf",sep=""),width=11,height=7)
# # ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = Correlation),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(Correlation, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# # dev.off()
# # 
# # sg_cor_pcs2<-correlate_pcs(pca_sg_mean_5000, meta$SEX, npcs = 10, min.cor = 0)
# # sg_cor_pcs2<-as.data.frame(sg_cor_pcs2); sg_cor_pcs2<-abs(sg_cor_pcs2)
# # tt<-cbind(rownames(sg_cor_pcs2),sg_cor_pcs2); 
# # tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "Correlation")
# # pdf(file = paste(Figs,subject1,"_pca_covariate_",genes,"genes_",date,".pdf",sep=""),width=11,height=7)
# # ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = Correlation),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(Correlation, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# # dev.off()
# # 
# # sg_r2_pcs1<-proportional_pcs(pca_sg_mean, meta$SEX, npcs = 10, min.cor = 0)
# # sg_r2_pcs1<-as.data.frame(sg_r2_pcs1); 
# # tt<-cbind(rownames(sg_r2_pcs1),sg_r2_pcs1); 
# # tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "adjusted.R2")
# # pdf(file = paste(Figs,subject1,"_pca_covariate_r2_",date,".pdf",sep=""),width=11,height=7)
# # ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = adjusted.R2),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(adjusted.R2, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# # dev.off()
# # 
# # sg_r2_pcs2<-proportional_pcs(pca_sg_mean_5000, meta$SEX, npcs = 10, min.cor = 0)
# # sg_r2_pcs2<-as.data.frame(sg_r2_pcs2); sg_r2_pcs2<-abs(sg_r2_pcs2)
# # tt<-cbind(rownames(sg_r2_pcs2),sg_r2_pcs2); 
# # tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "adjusted.R2")
# # pdf(file = paste(Figs,subject1,"_pca_covariate_r2_",genes,"genes_",date,".pdf",sep=""),width=11,height=7)
# # ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = adjusted.R2),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(adjusted.R2, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# # dev.off()
# # 
# # #-----------------------------------------------
# # ## high glucose
# # #-----------------------------------------------
# # pca_hg_mean <- prcomp(hg_input)
# # pdf(file = paste(Figs,subject2,"_pca_variance_explained_",date,".pdf",sep=""))
# # plot_variance_explained(pca_hg_mean) + theme_bw(base_size = 18)
# # dev.off()
# # 
# # pca_hg_mean_5000 <- prcomp(hg_input_5000, scale = FALSE, center = TRUE)
# # pdf(file = paste(Figs,subject2,"_DwC_and_DwoC_pca_variance_explained_",genes,"genes_",date,".pdf",sep=""))
# # plot_variance_explained(pca_hg_mean_5000) + theme_bw(base_size = 18)
# # dev.off()
# # 
# # data_explore_by_pca2(hg_input)
# # data_explore_by_pca2(hg_input_5000)
# # 
# # hg_cor_pcs1<-correlate_pcs(pca_hg_mean, meta, npcs = 10, min.cor = 0)
# # hg_cor_pcs1<-as.data.frame(hg_cor_pcs1); hg_cor_pcs1<-abs(hg_cor_pcs1)
# # tt<-cbind(rownames(hg_cor_pcs1),hg_cor_pcs1); 
# # tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "Correlation")
# # pdf(file = paste(Figs,subject2,"_pca_covariate_",date,".pdf",sep=""),width=11,height=7)
# # ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = Correlation),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(Correlation, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# # dev.off()
# # 
# # hg_cor_pcs2<-correlate_pcs(pca_hg_mean_5000, meta, npcs = 10, min.cor = 0)
# # hg_cor_pcs2<-as.data.frame(hg_cor_pcs2); hg_cor_pcs2<-abs(hg_cor_pcs2)
# # tt<-cbind(rownames(hg_cor_pcs2),hg_cor_pcs2); 
# # tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "Correlation")
# # pdf(file = paste(Figs,subject2,"_pca_covariate_",genes,"genes_",date,".pdf",sep=""),width=11,height=7)
# # ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = Correlation),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(Correlation, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# # dev.off()
# # 
# # 
# # hg_r2_pcs1<-proportional_pcs(pca_hg_mean, meta, npcs = 10, min.cor = 0)
# # hg_r2_pcs1<-as.data.frame(hg_r2_pcs1); hg_r2_pcs1<-abs(hg_r2_pcs1)
# # tt<-cbind(rownames(hg_r2_pcs1),hg_r2_pcs1); 
# # tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "adjusted.R2")
# # pdf(file = paste(Figs,subject2,"_pca_covariate_r2_",date,".pdf",sep=""),width=11,height=7)
# # ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = adjusted.R2),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(adjusted.R2, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# # dev.off()
# # 
# # hg_r2_pcs2<-proportional_pcs(pca_hg_mean_5000, meta, npcs = 10, min.cor = 0)
# # hg_r2_pcs2<-as.data.frame(hg_r2_pcs2); hg_r2_pcs2<-abs(hg_r2_pcs2)
# # tt<-cbind(rownames(hg_r2_pcs2),hg_r2_pcs2); 
# # tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "adjusted.R2")
# # pdf(file = paste(Figs,subject2,"_pca_covariate_r2_",genes,"genes_",date,".pdf",sep=""),width=11,height=7)
# # ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = adjusted.R2),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(adjusted.R2, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# # dev.off()
# # 
# # #-----------------------------------------------
# # ## delta
# # #-----------------------------------------------
# # pca_delta_mean <- prcomp(delta_input, scale = FALSE, center = TRUE)
# # pdf(file = paste(Figs,subject3,"_DwC_and_DwoC_pca_variance_explained_",date,".pdf",sep=""))
# # plot_variance_explained(pca_delta_mean, 15) + theme_bw(base_size = 18)
# # dev.off()
# # 
# # data_explore_by_pca2(delta_input)
# # 
# # delta_cor_pcs1<-correlate_pcs(pca_delta_mean, meta, npcs = 15, min.cor = 0)
# # delta_cor_pcs1<-as.data.frame(delta_cor_pcs1); delta_cor_pcs1<-abs(delta_cor_pcs1)
# # 
# # tt<-cbind(rownames(delta_cor_pcs1),delta_cor_pcs1); 
# # tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "Correlation")
# # pdf(file = paste(Figs,subject3,"_pca_covariate_",date,".pdf",sep=""),width=11,height=7)
# # ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = Correlation),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(Correlation, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# # dev.off()
# # 
# # delta_r2_pcs1<-proportional_pcs(pca_delta_mean, meta, npcs = 15, min.cor = 0)
# # delta_r2_pcs1<-as.data.frame(delta_r2_pcs1); delta_r2_pcs1<-abs(delta_r2_pcs1)
# # tt<-cbind(rownames(delta_r2_pcs1),delta_r2_pcs1); 
# # tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "adjusted.R2")
# # pdf(file = paste(Figs,subject3,"_pca_covariate_r2_",date,".pdf",sep=""),width=11,height=7)
# # ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = adjusted.R2),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(adjusted.R2, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
# #                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# # dev.off()
# 
# 
# 
# #------ PCA and Covariates
# meta<-read.table(paste(PhenotypeDir, "meanDeltaSampleInfo3.csv",sep=""), header=TRUE, sep=",", row.names=1, as.is=TRUE)
# genes<-5000
# 
# hg.T<-data.frame(t(hg)); 
# hg.T$Name<-substr(rownames(hg.T),1,nchar(rownames(hg.T))-2)
# hgagg<-aggregate(hg.T[,-ncol(hg.T)],list(hg.T$Name), mean)
# rownames(hgagg)<-hgagg$Group.1; hgagg$Group.1<-NULL
# hgmean<-apply(hgagg,2,as.numeric); rownames(hgmean)<-rownames(hgagg)
# hg_mean<-data.frame(t(hgmean)); colnames(hg_mean)<-substr(colnames(hg_mean),1,nchar(colnames(hg_mean))-5)
# hg_mean$CV<-CV(hg_mean); hg_mean$mean<-rowMeans(hg_mean[,-ncol(hg_mean)])
# hg_mean_sorted<-hg_mean[order(hg_mean$mean,decreasing = TRUE),]
# hg_input<-hg_mean_sorted[,rownames(meta)]
# hg_input_5000<-hg_mean_sorted[1:genes,rownames(meta)]
# 
# sg.T<-data.frame(t(sg)); 
# sg.T$Name<-substr(rownames(sg.T),1,nchar(rownames(sg.T))-2)
# sgagg<-aggregate(sg.T[,-ncol(sg.T)],list(sg.T$Name), mean)
# rownames(sgagg)<-sgagg$Group.1; sgagg$Group.1<-NULL
# sgmean<-apply(sgagg,2,as.numeric); rownames(sgmean)<-rownames(sgagg)
# sg_mean<-data.frame(t(sgmean)); colnames(sg_mean)<-substr(colnames(sg_mean),1,nchar(colnames(sg_mean))-5)
# sg_mean$CV<-CV(sg_mean); sg_mean$mean<-rowMeans(sg_mean[,-ncol(sg_mean)])
# sg_mean_sorted<-sg_mean[order(sg_mean$mean,decreasing = TRUE),]
# sg_input<-sg_mean_sorted[,rownames(meta)]
# sg_input_5000<-sg_mean_sorted[1:genes,rownames(meta)]
# 
# All_DE.T<-data.frame(t(All_DE))
# All_DE.T$Name<-substr(rownames(All_DE.T),1,nchar(rownames(All_DE.T))-2)
# All_DEagg<-aggregate(All_DE.T[,-ncol(All_DE.T)],list(All_DE.T$Name), mean)
# rownames(All_DEagg)<-All_DEagg$Group.1; All_DEagg$Group.1<-NULL
# meandelta<-apply(All_DEagg,2,as.numeric); rownames(meandelta)<-rownames(All_DEagg)
# delta_mean<-data.frame(t(meandelta)); delta_mean$CV<-CV(delta_mean); delta_mean$mean<-rowMeans(delta_mean[,-ncol(delta_mean)]) 
# delta_mean_sorted<-delta_mean[order(delta_mean$mean,decreasing = TRUE),]
# delta_input<-delta_mean_sorted[,rownames(meta)]
# delta_input_5000<-delta_mean_sorted[1:genes,rownames(meta)]
# 
# rm(hg_mean, sg_mean, delta_mean)
# 
# subject1="sg_mean";subject2="hg_mean";subject3="delta_mean"
# 
# #-----------------------------------------------
# ## standard glucose
# #-----------------------------------------------
# pca_sg_mean <- prcomp(sg_input, scale = FALSE, center = TRUE)
# pdf(file = paste(Figs,subject1,"_DwC_and_DwoC_pca_variance_explained_",date,".pdf",sep=""))
# plot_variance_explained(pca_sg_mean) + theme_bw(base_size = 18)
# dev.off()
# 
# pca_sg_mean_5000 <- prcomp(sg_input_5000, scale = FALSE, center = TRUE)
# pdf(file = paste(Figs,subject1,"_DwC_and_DwoC_pca_variance_explained_",genes,"genes_",date,".pdf",sep=""))
# plot_variance_explained(pca_sg_mean_5000) + theme_bw(base_size = 18)
# dev.off()
# 
# data_explore_by_pca2(sg_input)
# data_explore_by_pca2(sg_input_5000)
# 
# # correlate pcs and proportional pcs (edit the function in pca.R)
# sg_cor_pcs1<-correlate_pcs(pca_sg_mean, meta, npcs = 10, min.cor = 0)
# sg_cor_pcs1<-as.data.frame(sg_cor_pcs1); sg_cor_pcs1<-abs(sg_cor_pcs1)
# 
# tt<-cbind(rownames(sg_cor_pcs1),sg_cor_pcs1); 
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "Correlation")
# pdf(file = paste(Figs,subject1,"_pca_covariate_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = Correlation),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(Correlation, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# dev.off()
# 
# sg_cor_pcs2<-correlate_pcs(pca_sg_mean_5000, meta, npcs = 10, min.cor = 0)
# sg_cor_pcs2<-as.data.frame(sg_cor_pcs2); sg_cor_pcs2<-abs(sg_cor_pcs2)
# tt<-cbind(rownames(sg_cor_pcs2),sg_cor_pcs2); 
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "Correlation")
# pdf(file = paste(Figs,subject1,"_pca_covariate_",genes,"genes_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = Correlation),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(Correlation, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# dev.off()
# 
# sg_r2_pcs1<-proportional_pcs(pca_sg_mean, meta, npcs = 10, min.cor = 0)
# sg_r2_pcs1<-as.data.frame(sg_r2_pcs1); 
# tt<-cbind(rownames(sg_r2_pcs1),sg_r2_pcs1); 
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "adjusted.R2")
# pdf(file = paste(Figs,subject1,"_pca_covariate_r2_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = adjusted.R2),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(adjusted.R2, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# dev.off()
# 
# sg_r2_pcs2<-proportional_pcs(pca_sg_mean_5000, meta, npcs = 10, min.cor = 0)
# sg_r2_pcs2<-as.data.frame(sg_r2_pcs2); sg_r2_pcs2<-abs(sg_r2_pcs2)
# tt<-cbind(rownames(sg_r2_pcs2),sg_r2_pcs2); 
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "adjusted.R2")
# pdf(file = paste(Figs,subject1,"_pca_covariate_r2_",genes,"genes_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = adjusted.R2),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(adjusted.R2, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# dev.off()
# 
# #-----------------------------------------------
# ## high glucose
# #-----------------------------------------------
# pca_hg_mean <- prcomp(hg_input)
# pdf(file = paste(Figs,subject2,"_DwC_and_DwoC_pca_variance_explained_",date,".pdf",sep=""))
# plot_variance_explained(pca_hg_mean) + theme_bw(base_size = 18)
# dev.off()
# 
# pca_hg_mean_5000 <- prcomp(hg_input_5000, scale = FALSE, center = TRUE)
# pdf(file = paste(Figs,subject2,"_DwC_and_DwoC_pca_variance_explained_",genes,"genes_",date,".pdf",sep=""))
# plot_variance_explained(pca_hg_mean_5000) + theme_bw(base_size = 18)
# dev.off()
# 
# data_explore_by_pca2(hg_input)
# data_explore_by_pca2(hg_input_5000)
# 
# hg_cor_pcs1<-correlate_pcs(pca_hg_mean, meta, npcs = 10, min.cor = 0)
# hg_cor_pcs1<-as.data.frame(hg_cor_pcs1); hg_cor_pcs1<-abs(hg_cor_pcs1)
# tt<-cbind(rownames(hg_cor_pcs1),hg_cor_pcs1); 
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "Correlation")
# pdf(file = paste(Figs,subject2,"_pca_covariate_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = Correlation),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(Correlation, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# dev.off()
# 
# hg_cor_pcs2<-correlate_pcs(pca_hg_mean_5000, meta, npcs = 10, min.cor = 0)
# hg_cor_pcs2<-as.data.frame(hg_cor_pcs2); hg_cor_pcs2<-abs(hg_cor_pcs2)
# tt<-cbind(rownames(hg_cor_pcs2),hg_cor_pcs2); 
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "Correlation")
# pdf(file = paste(Figs,subject2,"_pca_covariate_",genes,"genes_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = Correlation),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(Correlation, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# dev.off()
# 
# 
# hg_r2_pcs1<-proportional_pcs(pca_hg_mean, meta, npcs = 10, min.cor = 0)
# hg_r2_pcs1<-as.data.frame(hg_r2_pcs1); hg_r2_pcs1<-abs(hg_r2_pcs1)
# tt<-cbind(rownames(hg_r2_pcs1),hg_r2_pcs1); 
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "adjusted.R2")
# pdf(file = paste(Figs,subject2,"_pca_covariate_r2_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = adjusted.R2),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(adjusted.R2, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# dev.off()
# 
# hg_r2_pcs2<-proportional_pcs(pca_hg_mean_5000, meta, npcs = 10, min.cor = 0)
# hg_r2_pcs2<-as.data.frame(hg_r2_pcs2); hg_r2_pcs2<-abs(hg_r2_pcs2)
# tt<-cbind(rownames(hg_r2_pcs2),hg_r2_pcs2); 
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "adjusted.R2")
# pdf(file = paste(Figs,subject2,"_pca_covariate_r2_",genes,"genes_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = adjusted.R2),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(adjusted.R2, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# dev.off()
# 
# #-----------------------------------------------
# ## delta
# #-----------------------------------------------
# pca_delta_mean <- prcomp(delta_input, scale = FALSE, center = TRUE)
# pdf(file = paste(Figs,subject3,"_DwC_and_DwoC_pca_variance_explained_",date,".pdf",sep=""))
# plot_variance_explained(pca_delta_mean, 15) + theme_bw(base_size = 18)
# dev.off()
# 
# data_explore_by_pca2(delta_input)
# 
# delta_cor_pcs1<-correlate_pcs(pca_delta_mean, meta, npcs = 15, min.cor = 0)
# delta_cor_pcs1<-as.data.frame(delta_cor_pcs1); delta_cor_pcs1<-abs(delta_cor_pcs1)
# 
# tt<-cbind(rownames(delta_cor_pcs1),delta_cor_pcs1); 
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "Correlation")
# pdf(file = paste(Figs,subject3,"_pca_covariate_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = Correlation),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(Correlation, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# dev.off()
# 
# delta_r2_pcs1<-proportional_pcs(pca_delta_mean, meta, npcs = 15, min.cor = 0)
# delta_r2_pcs1<-as.data.frame(delta_r2_pcs1); delta_r2_pcs1<-abs(delta_r2_pcs1)
# tt<-cbind(rownames(delta_r2_pcs1),delta_r2_pcs1); 
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "adjusted.R2")
# pdf(file = paste(Figs,subject3,"_pca_covariate_r2_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = adjusted.R2),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(adjusted.R2, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# dev.off()
# 
# # ----------------------------------------------------------------
# ##  remove confounder from meanDelta(growth rate, DCCT_SBP_MEAN, DCCT_DBP_MEAN)
# # ----------------------------------------------------------------
# 
# gr<-meta$GROWTH_RATE; dcct_sbp<-meta$DCCT_SBP_MEAN; dcct_dbp<-meta$DCCT_DBP_MEAN
# delta_input2<-matrix(1,nrow(delta_input),ncol(delta_input))
# 
# for (i in 1:nrow(delta_input)){
#   gene<-as.numeric(delta_input[i,])
#   delta_input2[i,]<-summary(lm(gene~gr+dcct_sbp+dcct_dbp))$residual
# }
# rownames(delta_input2)<-rownames(delta_input)
# colnames(delta_input2)<-colnames(delta_input)
# 
# pca_delta_mean2 <- prcomp(delta_input2, scale = FALSE, center = TRUE)
# pdf(file = paste(Figs,subject3,"_DwC_and_DwoC_no_confounder_pca_variance_explained_",date,".pdf",sep=""))
# plot_variance_explained(pca_delta_mean2, 15) + theme_bw(base_size = 15)
# dev.off()
# 
# data_explore_by_pca2(delta_input2)
# 
# delta_cor_pcs2<-correlate_pcs(pca_delta_mean2, meta, npcs = 11, min.cor = 0)
# delta_cor_pcs2<-as.data.frame(delta_cor_pcs2); delta_cor_pcs2<-abs(delta_cor_pcs2)
# 
# tt<-cbind(rownames(delta_cor_pcs2),delta_cor_pcs2); 
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "Correlation")
# pdf(file = paste(Figs,subject3,"_pca_no_confounder_covariate_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = Correlation),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(Correlation, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# dev.off()
# 
# delta_r2_pcs2<-proportional_pcs(pca_delta_mean2, meta, npcs = 11, min.cor = 0)
# delta_r2_pcs2<-as.data.frame(delta_r2_pcs2); 
# tt<-cbind(rownames(delta_r2_pcs2),delta_r2_pcs2); 
# tt.m<-melt(tt); colnames(tt.m)<-c("Covariate", "PC", "adjusted.R2")
# pdf(file = paste(Figs,subject3,"_pca_no_confounder_covariate_r2_",date,".pdf",sep=""),width=11,height=7)
# ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = adjusted.R2),colour = "white") + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(adjusted.R2, 2))) + theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.x = element_text(colour="grey20",size=10,face="bold"),
#                                                                                                                                                                                                   axis.text.y = element_text(colour="grey20",size=10,face="bold"),                                                                                                                                                                                                   axis.title.x = element_text(colour="grey20",size=10,face="bold"))  
# dev.off()
# 
# # ----------------------------------------------------------------
# ##  remove confounder from Delta (growth rate, DCCT_SBP_MEAN, DCCT_DBP_MEAN)
# # ----------------------------------------------------------------
# 
# gr<-rep(meta$GROWTH_RATE, each=3); dcct_sbp<-rep(meta$DCCT_SBP_MEAN, each=3); dcct_dbp<-rep(meta$DCCT_DBP_MEAN, each=3)
# adj_Diabetes<-matrix(1,nrow(Diabetes),ncol(Diabetes))
# 
# for (i in 1:nrow(Diabetes)){
#   gene<-as.numeric(Diabetes[i,])
#   adj_Diabetes[i,]<-summary(lm(gene~gr+dcct_sbp+dcct_dbp))$residual
# }
# rownames(adj_Diabetes)<-rownames(Diabetes)
# colnames(adj_Diabetes)<-colnames(Diabetes)
# 
# #--------------------------------------------------------------------------------------
# # limma: differential expression analysis for treatment effect
# # based on All_DE t-test
# #--------------------------------------------------------------------------------------
# 
# probeList <- rownames(All_DE)  # probeID (nuID) gene annotation
# if (require(lumiHumanAll.db) & require(annotate)){
#   geneSymbol<- getSYMBOL(probeList, 'lumiHumanAll.db')
#   geneName<- sapply(lookUp(probeList, 'lumiHumanAll.db', 'GENENAME'), function(x) x[1])
# }
# genes <-data.frame(ID=probeList, geneSymbol=geneSymbol, geneName=geneName, stringsAsFactors=FALSE)
# 
# 
# #----------------------------------
# # DE after confounding factors removed
# #----------------------------------
# 
# # subject='DwC_DwoC'
# # targets<-readTargets(paste(PhenotypeDir,"hg_sg_",subject,"_target.txt", sep=''))
# # blocks<-factor(targets$Subject)
# # colnames(adj_Diabetes)<-substr(colnames(adj_Diabetes),1,nchar(colnames(adj_Diabetes))-2)
# # 
# # corfit <- duplicateCorrelation(adj_Diabetes, block = blocks)
# # fitTrtMean <- lmFit(adj_Diabetes, block = blocks, cor = corfit$consensus.correlation)
# # efit.contrast=eBayes(fitTrtMean)
# # 
# # # qval.cutoff=0.1; FC.cutoff=0.17 # FC=1.2
# #  qval.cutoff=0; FC.cutoff=0 # FC=1.2
# # x5=topTable(efit.contrast,n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
# # y5=topTable(efit.contrast, n=nrow(genes),adjust.method="BH",genelist=genes)
# # x5_subset<-x5[x5$adj.P.Val<qval.cutoff & abs(x5$logFC) > FC.cutoff,]
# # 
# # pdf(file = paste(Figs,"treatment_effect_",subject,"_q",qval.cutoff,"_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# # with(y5, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Diabetes with and without complication subjects: high glucose vs. standard glucose (q.val <",qval.cutoff,"& FC >",FC.cutoff,sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-0.5,0.5)))
# # with(subset(y5, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
# # with(subset(y5, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
# # dev.off()
# # 
# # cat("number of genes (q-val <0.05): ",sum(x2$adj.P.Val<qval.cutoff))
# # write.csv(x5_subset, paste(Output,"treatment_effect_",subject,"_",qval.cutoff,"_FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
# # write.csv(x5, paste(Output,"treatment_effect_",subject,"_",qval.cutoff,"_",date,".csv",sep=""),row.names=F)
# # write.csv(y5, paste(Output,"treatment_effect_",subject,"_",date,".csv",sep=""),row.names=F)
# # 
# # pdf(file = paste(Figs,"treatment_effect_",subject,"_pvalue_distribution",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# # par(mfrow=c(1,2))
# # hist(y5$P.Value, main="",xlab="p-value")
# # hist(y5$adj.P.Val,main="", xlab="adjusted p-value")
# # mtext(paste("treatment effect ", subject ," subject (high glucose vs. standard glucose)",sep=""), side=3, outer=TRUE, line=-3)
# # dev.off()
# 
# #--------------------------------------------------------------------------------------
# # High Glucose (HG) vs Standard Glucose (SG)
# # based on All_DE t-test
# #--------------------------------------------------------------------------------------
# 
# # create a design matrix
# all_targets<-readTargets(paste(PhenotypeDir,"hg_sg_all_target.txt", sep=''))
# 
# colnames(All_DE)<-substr(colnames(All_DE),1,nchar(colnames(All_DE))-2)
# blocks<-factor(all_targets$Subject)
# 
# corfit <- duplicateCorrelation(All_DE, block = blocks)
# fitTrtMean <- lmFit(All_DE, block = blocks, cor = corfit$consensus.correlation)
# efit.contrast=eBayes(fitTrtMean)
# 
# qval.cutoff=0.1; FC.cutoff=0.17 # FC=1.2
# 
# x1=topTable(efit.contrast,n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
# y1=topTable(efit.contrast, n=nrow(genes),adjust.method="BH",genelist=genes)
# x1_subset<-x1[x1$adj.P.Val<qval.cutoff & abs(x1$logFC) > FC.cutoff,]
# 
# pdf(file = paste(Figs,"treatment_effect_q",qval.cutoff,"_",date,".pdf",sep=""))
# with(y1, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("high glucose vs. standard glucose (q.val <",qval.cutoff," & FC >)",FC.cutoff,sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-0.5,0.5)))
# with(subset(y1, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
# with(subset(y1, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
# dev.off()
# 
# 
# cat("number of genes (q-val <0.05): ",sum(x2$adj.P.Val<qval.cutoff))
# write.csv(x1_subset, paste(Output,"treatment_effect_qval_",qval.cutoff,"_FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
# write.csv(x1, paste(Output,"treatment_effect_qval_",qval.cutoff,"_",date,".csv",sep=""),row.names = F)
# write.csv(y1, paste(Output,"treatment_effect_",date,".csv",sep=""),row.names = F)
# 
# pdf(file = paste(Figs,"treatment_effect_pvalue_distribution_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# par(mfrow=c(1,2))
# hist(y1$P.Value, main="",xlab="p-value")
# hist(y1$adj.P.Val,main="", xlab="adjusted p-value")
# mtext("treatment effect (high glucose vs. standard glucose)", side=3, outer=TRUE, line=-3)
# dev.off()
# 
# #------normal subject only
# norm_targets<-readTargets(paste(PhenotypeDir,"hg_sg_normal_target.txt", sep=''))
# blocks<-factor(norm_targets$Subject)
# 
# # colnames(NoD)<-substr(colnames(NoD),1,nchar(colnames(NoD))-3)
# corfit <- duplicateCorrelation(NoD_DE, block = blocks)
# fitTrtMean <- lmFit(NoD_DE, block = blocks, cor = corfit$consensus.correlation)
# efit.contrast=eBayes(fitTrtMean)
# 
# qval.cutoff=0.2; FC.cutoff=0.17 # FC=1.2
# 
# # standard glucose Normal vs Diabetes without Complication
# x2=topTable(efit.contrast,n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
# y2=topTable(efit.contrast, n=nrow(genes),adjust.method="BH",genelist=genes)
# x2_subset<-x2[x2$adj.P.Val<qval.cutoff & abs(x2$logFC) > FC.cutoff,]
# 
# pdf(file = paste(Figs,"treatment_effect_normal_qval_",qval_cutoff,"_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# with(y2, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("normal subject high glucose vs. standard glucose (q.val <",qval.cutoff,"& FC >",FC.cutoff,sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-0.5,0.5)))
# with(subset(y2, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
# with(subset(y2, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
# dev.off()
# 
# cat("number of genes (q-val <0.1): ",sum(x2$adj.P.Val<qval.cutoff))
# write.csv(x2_subset, paste(Output,"treatment_effect_normal_q",qval.cutoff,"_FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
# write.csv(x2, paste(Output,"treatment_effect_normal_q",qval.cutoff,"_",date,".csv",sep=""),row.names=F)
# write.csv(y2, paste(Output,"treatment_effect_normal_",date,".csv",sep=""),row.names=F)
# 
# pdf(file = paste(Figs,"treatment_effect_normal_pvalue_distribution_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# par(mfrow=c(1,2))
# hist(y2$P.Value, main="",xlab="p-value")
# hist(y2$adj.P.Val,main="", xlab="adjusted p-value")
# mtext("treatment effect normal subject (high glucose vs. standard glucose)", side=3, outer=TRUE, line=-3)
# dev.off()
# 
# #------DwC subject only
# subject='DwC'
# dwc_targets<-readTargets(paste(PhenotypeDir,"hg_sg_",subject,"_target.txt", sep=''))
# blocks<-factor(dwc_targets$Subject)
# # colnames(DwC_DE)<-substr(colnames(DwC_DE),1,nchar(colnames(DwC_DE))-2)
# 
# corfit <- duplicateCorrelation(DwC_DE, block = blocks)
# fitTrtMean <- lmFit(DwC_DE, block = blocks, cor = corfit$consensus.correlation)
# efit.contrast=eBayes(fitTrtMean)
# 
# qval.cutoff=0.1; FC.cutoff=0.17 # FC=1.2
# x3=topTable(efit.contrast,n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
# y3=topTable(efit.contrast, n=nrow(genes),adjust.method="BH",genelist=genes)
# x3_subset<-x3[x3$adj.P.Val<qval.cutoff & abs(x3$logFC) > FC.cutoff,]
# 
# pdf(file = paste(Figs,"treatment_effect_",subject,"_q",qval.cutoff,"_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# with(y3, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Diabetes with complication subjects: high glucose vs. standard glucose (q.val <",qval.cutoff,"& FC >",FC.cutoff,sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-0.5,0.5)))
# with(subset(y3, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
# with(subset(y3, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
# dev.off()
# 
# cat("number of genes (q-val <0.05): ",sum(x2$adj.P.Val<qval.cutoff))
# 
# write.csv(x3_subset, paste(Output,"treatment_effect_DwC_q",qval.cutoff,"_FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
# write.csv(x3, paste(Output,"treatment_effect_",subject,"_q",qval.cutoff,"_",date,".csv",sep=""),row.names=F)
# write.csv(y3, paste(Output,"treatment_effect_",subject,"_",date,".csv",sep=""),row.names=F)
# 
# pdf(file = paste(Figs,"treatment_effect_",subject,"_pvalue_distribution",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# par(mfrow=c(1,2))
# hist(y3$P.Value, main="",xlab="p-value")
# hist(y3$adj.P.Val,main="", xlab="adjusted p-value")
# mtext(paste("treatment effect ", subject ," subject (high glucose vs. standard glucose)",sep=""), side=3, outer=TRUE, line=-3)
# dev.off()
# 
# #------DwoC subject only
# subject='DwoC'
# targets<-readTargets(paste(PhenotypeDir,"hg_sg_",subject,"_target.txt", sep=''))
# blocks<-factor(targets$Subject)
# 
# corfit <- duplicateCorrelation(DwoC_DE, block = blocks)
# fitTrtMean <- lmFit(DwoC_DE, block = blocks, cor = corfit$consensus.correlation)
# efit.contrast=eBayes(fitTrtMean)
# 
# qval.cutoff=0.1; FC.cutoff=0.17 # FC=1.2
# x4=topTable(efit.contrast,n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
# y4=topTable(efit.contrast, n=nrow(genes),adjust.method="BH",genelist=genes)
# x4_subset<-x4[x4$adj.P.Val<qval.cutoff & abs(x4$logFC) > FC.cutoff,]
# 
# pdf(file = paste(Figs,"treatment_effect_",subject,"_q",qval.cutoff,"_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# with(y4, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Diabetes with complication subjects: high glucose vs. standard glucose (q.val <",qval.cutoff,"& FC >",FC.cutoff,sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-0.5,0.5)))
# with(subset(y4, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
# with(subset(y4, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
# dev.off()
# 
# cat("number of genes (q-val <0.05): ",sum(x2$adj.P.Val<qval.cutoff))
# write.csv(x4_subset, paste(Output,"treatment_effect_DwoC_q",qval.cutoff,"_FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
# write.csv(x4, paste(Output,"treatment_effect_",subject,"_",qval.cutoff,"_",date,".csv",sep=""),row.names=F)
# write.csv(y4, paste(Output,"treatment_effect_",subject,"_",date,".csv",sep=""),row.names=F)
# 
# pdf(file = paste(Figs,"treatment_effect_",subject,"_pvalue_distribution",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# par(mfrow=c(1,2))
# hist(y4$P.Value, main="",xlab="p-value")
# hist(y4$adj.P.Val,main="", xlab="adjusted p-value")
# mtext(paste("treatment effect ", subject ," subject (high glucose vs. standard glucose)",sep=""), side=3, outer=TRUE, line=-3)
# dev.off()
# 
# #------Dwc + DwoC subjects
# subject='DwC_DwoC'
# targets<-readTargets(paste(PhenotypeDir,"hg_sg_",subject,"_target.txt", sep=''))
# blocks<-factor(targets$Subject)
# 
# corfit <- duplicateCorrelation(Diabetes, block = blocks)
# fitTrtMean <- lmFit(Diabetes, block = blocks, cor = corfit$consensus.correlation)
# efit.contrast=eBayes(fitTrtMean)
# 
# qval.cutoff=0.1; FC.cutoff=0.17 # FC=1.2
# x5=topTable(efit.contrast,n=nrow(genes), p.value=qval.cutoff,adjust.method="BH",genelist=genes)
# y5=topTable(efit.contrast, n=nrow(genes),adjust.method="BH",genelist=genes)
# x5_subset<-x5[x5$adj.P.Val<qval.cutoff & abs(x5$logFC) > FC.cutoff,]
# 
# pdf(file = paste(Figs,"treatment_effect_",subject,"_q",qval.cutoff,"_",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# with(y5, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main=paste("Diabetes with and without complication subjects: high glucose vs. standard glucose (q.val <",qval.cutoff,"& FC >",FC.cutoff,sep=""), xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-0.5,0.5)))
# with(subset(y5, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
# with(subset(y5, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
# dev.off()
# 
# cat("number of genes (q-val <0.05): ",sum(x2$adj.P.Val<qval.cutoff))
# write.csv(x5_subset, paste(Output,"treatment_effect_",subject,"_",qval.cutoff,"_FC_",FC.cutoff,"_",date,".csv",sep=""),row.names = F)
# write.csv(x5, paste(Output,"treatment_effect_",subject,"_",qval.cutoff,"_",date,".csv",sep=""),row.names=F)
# write.csv(y5, paste(Output,"treatment_effect_",subject,"_",date,".csv",sep=""),row.names=F)
# 
# pdf(file = paste(Figs,"treatment_effect_",subject,"_pvalue_distribution",date,".pdf",sep=""),width=18,height=12,pointsize=20)
# par(mfrow=c(1,2))
# hist(y5$P.Value, main="",xlab="p-value")
# hist(y5$adj.P.Val,main="", xlab="adjusted p-value")
# mtext(paste("treatment effect ", subject ," subject (high glucose vs. standard glucose)",sep=""), side=3, outer=TRUE, line=-3)
# dev.off()
