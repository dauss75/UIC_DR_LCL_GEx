# Inter_Intra_variance2 <- function(x){
#   output<-matrix(, nrow = nrow(x), ncol = 5)
#   colnames(output)<-c("geneID","inter-individual variance", "intra-individual variance", "F-value", "p-value")
#   for (i in 1:nrow(x)) {
#     data = data.frame(gene = t(x[i,]), group = factor(substr(colnames(x),1,nchar(colnames(x))-2))) 
#     colnames(data)[1]<-"gene"
#     out = anova(lm(gene~group,data))
#     inter_val = out["group","Mean Sq"]
#     intra_val = out["Residuals","Mean Sq"]
#     Fval = out["group","F value"]
#     Pval = out["group","Pr(>F)"]
#     output[i,]<-c(row.names(dat)[i],inter_val, intra_val, Fval, Pval)
#   }
#   return(output)
# }

Inter_Intra_variance <- function(x){
  finalMatrix <- foreach(i=1:nrow(x), .combine=rbind) %dopar% {
    data = data.frame(gene = t(x[i,]), group = factor(substr(colnames(x),1,nchar(colnames(x))-2))) 
    colnames(data)[1]<-"gene"
    out = anova(lm(gene~group,data))
    inter_val = out["group","Mean Sq"]
    intra_val = out["Residuals","Mean Sq"]
    Fval = out["group","F value"]
    Pval = out["group","Pr(>F)"]
    
    c(row.names(x)[i],inter_val, intra_val, Fval, Pval)
  }
  colnames(finalMatrix)<-c("nuID","inter-individual variance", "intra-individual variance", "F-value", "p-value")
  return(finalMatrix)
}

CV <- function(x, na.rm=TRUE) {
  #row: genes, column: samples
  x1<-apply(x,2,as.numeric)
  x1_SD<-apply(x1,1,sd)
  x1_median<-apply(x1,1,median)  
  y1<-x1_SD/x1_median
  return(y1)
}

plot.multi.dens <- function(s)
{
  junk.x = NULL
  junk.y = NULL
  for(i in 1:length(s)) {
    junk.x = c(junk.x, density(s[[i]])$x)
    junk.y = c(junk.y, density(s[[i]])$y)
  }
  xr <- range(junk.x)
  yr <- range(junk.y)
  plot(density(s[[1]]), xlim = xr, ylim = yr, main = "", xlab="log2(variance)")
  for(i in 1:length(s)) {
    lines(density(s[[i]]), xlim = xr, ylim = yr, col = i)
  }
}

var_proc <-function(x,y){
  #x: data, y: gene annotation
  var<-Inter_Intra_variance(x)
  cv <-CV(x)
  variation<-as.data.frame(cbind(var,cv))
  colnames(variation)[6]<-"CV"
  rownames(variation)<-variation$nuID; variation$nuID<-NULL
  final_output<-cbind(y[rownames(variation),2],variation)
  colnames(final_output)[1]<-"geneID"
  final_output_sorted<-final_output[order(final_output$`p-value`),]
  return(final_output_sorted)
}

library(foreach)
library(doParallel)
library(ggplot2)
library(limma)
library(gtools); library(bioDist); library(calibrate)

HOME <- "/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/"  # change if needed
InputDir <- paste(HOME,'data/input/', sep='')
InputFiles = list.files(path=InputDir, pattern='txt')
ControlDir <- paste(HOME,'data/control/', sep='')
ControlFiles = list.files(path=ControlDir, pattern='txt')
PhenotypeDir = paste(HOME,'data/phenotype/', sep='')
RData=paste(HOME,'rdata/',sep=''); dir.create(RData, showWarnings = FALSE)
Figs = paste(HOME,'Figs/', sep=''); dir.create(Figs, showWarnings = FALSE)
Output = paste(HOME,'output/', sep=''); dir.create(Figs, showWarnings = FALSE)

dat<-local(get(load(file=paste(RData,"normalizedDataMatrix.sorted.RData",sep=""))))
gene_annotation<-read.csv(paste(PhenotypeDir,"raw_delta_gene_anno.csv",sep=""))
rownames(gene_annotation)<-gene_annotation$GeneCode

# parallization
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

# for normal samples
NoD<-dat[,grep("NoD",colnames(dat))] # grep all normal samples
DwC<-dat[,grep("DwC",colnames(dat))] # grep all DwC samples
DwoC<-dat[,grep("DwoC",colnames(dat))] # grep all DwoC samples

# separate high glucose and standard glucose
NoD_norm<-NoD[,grep("norm",colnames(NoD))]
NoD_30mM<-NoD[,grep("30mM",colnames(NoD))]

DwC_norm<-DwC[,grep("norm",colnames(DwC))]
DwC_30mM<-DwC[,grep("30mM",colnames(DwC))]

DwoC_norm<-DwoC[,grep("norm",colnames(DwoC))]
DwoC_30mM<-DwoC[,grep("30mM",colnames(DwoC))]

sg<-dat[,grep("norm",colnames(dat))]
hg<-dat[,grep("30mM",colnames(dat))]

# normal samples with standard glucose (sg)
NoD_sg_output<-var_proc(NoD_norm,gene_annotation)
write.csv(NoD_sg_output, paste(Output,"Normal_standard_glucose_variation_01.19.2017.csv",sep=""))
subset
tmp<-NoD_sg_output[,c(2,3)]; colnames(tmp)<-c("inter","intra")
ggplot(tmp, aes(y = inter, x = intra)) + geom_point()
# NoD_sg_output_top5k<-NoD_sg_output[order(NoD_sg_output$CV,decreasing = TRUE)[1:5000],]

#---density plot
inter<-log2(as.numeric(as.character(NoD_sg_output_top5k$`inter-individual variance`)))
intra<-log2(as.numeric(as.character(NoD_sg_output_top5k$`intra-individual variance`)))
data<-as.data.frame(cbind(inter,intra))
png(paste(Figs,"normal_standard_glucose_01.13.2017.png",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="Normal samples with standard glucose")
dev.off()
#---

# normal samples with high glucose (hg)
NoD_hg_output<-var_proc(NoD_30mM,gene_annotation)
write.csv(NoD_hg_output, paste(Output,"Normal_high_glucose_variation_01.13.2017.csv",sep=""))
NoD_hg_output_top5k<-NoD_hg_output[order(NoD_hg_output$CV,decreasing = TRUE)[1:5000],]

#---density plot
inter<-log2(as.numeric(as.character(NoD_hg_output_top5k$`inter-individual variance`)))
intra<-log2(as.numeric(as.character(NoD_hg_output_top5k$`intra-individual variance`)))
data<-as.data.frame(cbind(inter,intra))

png(paste(Figs,"normal_high_glucose_01.13.2017.png",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="Normal samples with high glucose")
dev.off()
#---

# diabets with complication samples with standard glucose (sg)
DwC_sg_output<-var_proc(DwC_norm,gene_annotation)
write.csv(DwC_sg_output, paste(Output,"Diabetes_with_complication_standard_glucose_variation_01.13.2017.csv",sep=""))
DwC_sg_output_top5k<-DwC_sg_output[order(DwC_sg_output$CV,decreasing = TRUE)[1:5000],]

#---density plot
inter<-log2(as.numeric(as.character(DwC_sg_output_top5k$`inter-individual variance`)))
intra<-log2(as.numeric(as.character(DwC_sg_output_top5k$`intra-individual variance`)))
data<-as.data.frame(cbind(inter,intra))
png(paste(Figs,"DwC_standard_glucose_01.13.2017.png",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="diabets with complication samples with standard glucose")
dev.off()
#---

# diabets with complication samples with high glucose (hg)
DwC_hg_output<-var_proc(DwC_30mM,gene_annotation)
write.csv(DwC_hg_output, paste(Output,"DwC_high_glucose_variation_01.13.2017.csv",sep=""))
DwC_hg_output_top5k<-DwC_hg_output[order(DwC_hg_output$CV,decreasing = TRUE)[1:5000],]

#---density plot
inter<-log2(as.numeric(as.character(DwC_hg_output_top5k$`inter-individual variance`)))
intra<-log2(as.numeric(as.character(DwC_hg_output_top5k$`intra-individual variance`)))
data<-as.data.frame(cbind(inter,intra))

png(paste(Figs,"DwC_high_glucose_01.13.2017.png",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="diabets with complication samples with high glucose")
dev.off()
#---

# diabets without complication samples with standard glucose (sg)
DwoC_sg_output<-var_proc(DwoC_norm,gene_annotation)
write.csv(DwoC_sg_output, paste(Output,"Diabetes_without_complication_standard_glucose_variation_01.13.2017.csv",sep=""))
DwoC_sg_output_top5k<-DwoC_sg_output[order(DwoC_sg_output$CV,decreasing = TRUE)[1:5000],]

#---density plot
inter<-log2(as.numeric(as.character(DwoC_sg_output_top5k$`inter-individual variance`)))
intra<-log2(as.numeric(as.character(DwoC_sg_output_top5k$`intra-individual variance`)))
data<-as.data.frame(cbind(inter,intra))
png(paste(Figs,"DwoC_standard_glucose_01.13.2017.png",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="diabets without complication samples with standard glucose")
dev.off()
#---

# diabets without complication samples with high glucose (hg)
DwoC_hg_output<-var_proc(DwoC_30mM,gene_annotation)
write.csv(DwoC_hg_output, paste(Output,"Diabetes_without_complication_high_glucose_variation_01.13.2017.csv",sep=""))
DwoC_hg_output_top5k<-DwoC_hg_output[order(DwoC_hg_output$CV,decreasing = TRUE)[1:5000],]

#---density plot
inter<-log2(as.numeric(as.character(DwoC_hg_output_top5k$`inter-individual variance`)))
intra<-log2(as.numeric(as.character(DwoC_hg_output_top5k$`intra-individual variance`)))
data<-as.data.frame(cbind(inter,intra))

png(paste(Figs,"DwoC_high_glucose_01.13.2017.png",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="diabets without complication samples with high glucose")
dev.off()
#---

# inter and intra-individual variability for standard glucose, high glucose, all-together
#all standard glucose
sg_output<-var_proc(sg,gene_annotation)
sg_output<-sg_output[order(sg_output$`p-value`,decreasing = TRUE),]
write.csv(sg_output, paste(Output,"standard_glucose_samples_inter_intra_individual_variation_01.13.2017.csv",sep=""))
sg_output_top5k<-sg_output[order(sg_output$CV,decreasing = TRUE)[1:5000],]

#---density plot
inter<-log2(as.numeric(as.character(sg_output_top5k$`inter-individual variance`)))
intra<-log2(as.numeric(as.character(sg_output_top5k$`intra-individual variance`)))
data<-as.data.frame(cbind(inter,intra))

png(paste(Figs,"sg_output_01.13.2017.png",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="standard glucose samples inter and intra-individual variability")
dev.off()

#all hg
hg_output<-var_proc(hg,gene_annotation)
hg_output<-hg_output[order(hg_output$`p-value`,decreasing = TRUE),]
write.csv(hg_output, paste(Output,"high_glucose_samples_inter_intra_individual_variation_01.13.2017.csv",sep=""))
hg_output_top5k<-hg_output[order(hg_output$CV,decreasing = TRUE)[1:5000],]

#---density plot
inter<-log2(as.numeric(as.character(hg_output_top5k$`inter-individual variance`)))
intra<-log2(as.numeric(as.character(hg_output_top5k$`intra-individual variance`)))
data<-as.data.frame(cbind(inter,intra))

png(paste(Figs,"hg_output_01.13.2017.png",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="high glucose samples inter and intra-individual variability")
dev.off()
#all
all_output<-var_proc(dat,gene_annotation)
all_output<-all_output[order(all_output$`p-value`,decreasing = TRUE),]
write.csv(all_output, paste(Output,"all_samples_inter_intra_individual_variation_01.13.2017.csv",sep=""))
all_output_top5k<-all_output[order(all_output$CV,decreasing = TRUE)[1:5000],]

#---density plot
inter<-log2(as.numeric(as.character(all_output_top5k$`inter-individual variance`)))
intra<-log2(as.numeric(as.character(all_output_top5k$`intra-individual variance`)))
data<-as.data.frame(cbind(inter,intra))

png(paste(Figs,"all_output_01.13.2017.png",sep=""))
plot.multi.dens( list(inter, intra))
legend("topright", c("inter-individual","intra-individual"),col = c("black", "red"), lty = c(1, 1))
title(main="all samples inter and intra-individual variability")
dev.off()
#---
stopCluster(cl)

rm(NoD,DwC,DwoC,NoD_norm,NoD_30mM,DwC_norm,DwC_30mM,DwoC_norm,DwoC_30mM,sg,hg,NoD_sg_output,NoD_sg_output_top5k,
   inter,intra,data,NoD_hg_output,NoD_hg_output_top5k,DwC_hg_output_top5k,DwC_hg_output,DwoC_sg_output_top5k,
   DwoC_hg_output_top5k,DwoC_hg_output,DwoC_sg_output,hg_output,hg_output_top5k,sg_output,sg_output_top5k,
   DwC_sg_output,DwC_sg_output_top5k)

# rows: genes, columns: samples
#--------------------------------------------------------------------------------------
# limma: differential expression analysis -  High Glucose (HG) vs Standard Glucose (SG)
# based on ANOVA
#--------------------------------------------------------------------------------------

# create a design matrix
targets<-readTargets(paste(PhenotypeDir,"target4.txt", sep=''))

#--- test 1 (treatment)
trts<-factor(targets$Group) 
blocks<-factor(targets$Subject)

design.trt=model.matrix(~0+trts)
corfit <- duplicateCorrelation(dat, design.trt, block = blocks)
fitTrtMean <- lmFit(dat, design.trt, block = blocks, cor = corfit$consensus.correlation)

contrast.matrix <- makeContrasts( sg_NoDVsDwoC=trtsNoD_sg - trtsDwoC_sg, sg_DwCVsDwoC=trtsDwC_sg-trtsDwoC_sg, 
                                  hg_NoDVsDwoC=trtsNoD_hg - trtsDwoC_hg, hg_DwCVsDwoC=trtsDwC_hg-trtsDwoC_hg,
                                  NoD=trtsNoD_hg - trtsNoD_sg, DwC=trtsDwC_hg - trtsDwC_sg,
                                  DwoC=trtsDwoC_hg - trtsDwoC_sg,
                                  hg_sg=(trtsNoD_hg+trtsDwC_hg+trtsDwoC_hg)/3 - (trtsNoD_sg+trtsDwC_sg+trtsDwoC_sg)/3,
                                  levels=design.trt) # create a contrast matrix    
fit.contrast=contrasts.fit(fitTrtMean,contrast.matrix)
efit.contrast=eBayes(fit.contrast)
par(mfrow=c(2,4))
for (i in 1:ncol(efit.contrast$p.value)) {
  hist(efit.contrast$p.value[,i], main=colnames(efit.contrast$p.value)[i])
}


probeList <- rownames(dat)  # probeID (nuID) gene annotation
if (require(lumiHumanAll.db) & require(annotate)){
  geneSymbol<- getSYMBOL(probeList, 'lumiHumanAll.db')
  geneName<- sapply(lookUp(probeList, 'lumiHumanAll.db', 'GENENAME'), function(x) x[1])
}

genes <-data.frame(ID=probeList, geneSymbol=geneSymbol, geneName=geneName, stringsAsFactors=FALSE)
qval.cutoff=0.01; FC.cutoff=1
# standard glucose Normal vs Diabetes without Complication
sg_NoDVsDwoC_qval1=topTable(efit.contrast,coef="sg_NoDVsDwoC", n=5000, p.value=qval.cutoff,adjust.method="BH",genelist=genes)
sg_NoDVsDwoC=topTable(efit.contrast,coef="sg_NoDVsDwoC", n=nrow(genes),adjust.method="BH",genelist=genes)
pdf(file = paste(Figs,"Normal Vs Diabetes without Complication (standard glucose) qval0.01.pdf",sep=""),width=18,height=12,pointsize=20)
with(sg_NoDVsDwoC, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main="Normal Vs Diabets without Complication for standard glucose (q.val <0.01 & FC >2)", xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-2.5,2.5)))
with(subset(sg_NoDVsDwoC, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(sg_NoDVsDwoC, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()
cat("number of genes (q-val <0.01 & FC >2): ",sum(sg_NoDVsDwoC$adj.P.Val<qval.cutoff & abs(sg_NoDVsDwoC$logFC)>FC.cutoff))
write.csv(sg_NoDVsDwoC_qval1, paste(Output,"Normal Vs Diabetes without Complication (standard glucose) qval0.01.csv",sep=""))
rm(sg_NoDVsDwoC,sg_NoDVsDwoC_qval1)

# high glucose Normal vs Diabetes without Complication
hg_NoDVsDwoC_qval1=topTable(efit.contrast,coef="hg_NoDVsDwoC", n=5000, p.value=qval.cutoff,adjust.method="BH",genelist=genes)
hg_NoDVsDwoC=topTable(efit.contrast,coef="hg_NoDVsDwoC", n=nrow(genes),adjust.method="BH",genelist=genes)
pdf(file = paste(Figs,"Normal Vs Diabetes without Complication (high glucose) qval0.01.pdf",sep=""),width=18,height=12,pointsize=20)
with(hg_NoDVsDwoC, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main="Normal Vs Diabets without Complication for high glucose (q.val <0.01 & FC >2)", xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-2.5,2.5)))
with(subset(hg_NoDVsDwoC, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(hg_NoDVsDwoC, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()
cat("number of genes (q-val <0.01 & FC >2): ",sum(hg_NoDVsDwoC$adj.P.Val<qval.cutoff & abs(hg_NoDVsDwoC$logFC)>FC.cutoff))
write.csv(hg_NoDVsDwoC_qval1, paste(Output,"Normal Vs Diabetes without Complication (high glucose) qval0.01.csv",sep=""))
rm(hg_NoDVsDwoC,hg_NoDVsDwoC_qval1)

qval.cutoff=0.05; FC.cutoff=0.5849
# standard glucose Diabetes with Complication vs Diabetes without Complication
sg_DwCVsDwoC_qval5=topTable(efit.contrast,coef="sg_DwCVsDwoC", n=5000, p.value=qval.cutoff,adjust.method="BH",genelist=genes)
sg_DwCVsDwoC=topTable(efit.contrast,coef="sg_DwCVsDwoC", n=nrow(genes),adjust.method="BH",genelist=genes)
pdf(file = paste(Figs,"Diabetes with Complication Vs Diabetes without Complication (standard glucose) qval0.05.pdf",sep=""),width=18,height=12,pointsize=20)
with(sg_DwCVsDwoC, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main="Diabetes with Complication Vs Diabets without Complication for standard glucose (q.val <0.05 & FC >1.5)", xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-2.5,2.5)))
with(subset(sg_DwCVsDwoC, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(sg_DwCVsDwoC, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()
cat("number of genes (q-val <0.05 & FC >1.5): ",sum(sg_DwCVsDwoC$adj.P.Val<qval.cutoff & abs(sg_DwCVsDwoC$logFC)>FC.cutoff))
write.csv(sg_DwCVsDwoC_qval5, paste(Output,"Diabetes with Complication Vs Diabetes without Complication (standard glucose) qval0.05.csv",sep=""))
rm(sg_DwCVsDwoC_qval5, sg_DwCVsDwoC)

# high glucose Diabetes with Complication vs Diabetes without Complication
hg_DwCVsDwoC_qval5=topTable(efit.contrast,coef="hg_DwCVsDwoC", n=5000, p.value=qval.cutoff,adjust.method="BH",genelist=genes)
hg_DwCVsDwoC=topTable(efit.contrast,coef="hg_DwCVsDwoC", n=nrow(genes),adjust.method="BH",genelist=genes)
pdf(file = paste(Figs,"Diabetes with Complication Vs Diabetes without Complication (high glucose) qval0.05.pdf",sep=""),width=18,height=12,pointsize=20)
with(hg_DwCVsDwoC, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main="Diabetes with Complication Vs Diabets without Complication for high glucose (q.val <0.05 & FC >1.5)", xlab="log2 fold change", ylab="-log10(adj.P.Val)", xlim=c(-2.5,2.5)))
with(subset(hg_DwCVsDwoC, adj.P.Val<qval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(subset(hg_DwCVsDwoC, adj.P.Val<qval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
dev.off()
cat("number of genes (q-val <0.05 & FC >1.5): ",sum(hg_DwCVsDwoC$adj.P.Val<qval.cutoff & abs(hg_DwCVsDwoC$logFC)>FC.cutoff))

write.csv(hg_DwCVsDwoC_qval5, paste(Output,"Diabetes with Complication Vs Diabetes without Complication (high glucose) qval0.05.csv",sep=""))
rm(hg_DwCVsDwoC_qval5, hg_DwCVsDwoC)

# extra: treatment effect
pval.cutoff=0.05; FC.cutoff=0.5849

# normal subject (high vs standard glucose)
NoD=topTable(efit.contrast,coef="NoD", n=nrow(genes),adjust.method="BH",genelist=genes)
NoD_pval5=NoD[NoD$P.Value<pval.cutoff,]
pdf(file = paste(Figs,"Normal subject (high vs standard glucose).pdf",sep=""),width=18,height=12,pointsize=20)
with(NoD, plot(logFC, -log10(P.Value), pch=20, col="gray", main="Normal Subject (high vs. standard glucose) (p.value <0.05)", xlab="log2 fold change", ylab="-log10(P.val)", xlim=c(-1,1)))
with(subset(NoD, P.Value<pval.cutoff ), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(NoD, P.Value<pval.cutoff) , textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()
cat("number of genes (p-val <0.05): ",sum(NoD$P.Value<pval.cutoff))
write.csv(NoD_pval5, paste(Output,"Normal Subject (high vs standard glucose) pval0.05.csv",sep=""))
rm(NoD, NoD_pval5)

# diabetes with complication subject (high vs standard glucose)
DwC=topTable(efit.contrast,coef="DwC", n=nrow(genes),adjust.method="BH",genelist=genes)
DwC_pval5=DwC[DwC$P.Value<pval.cutoff,]
pdf(file = paste(Figs,"Diabetes with complication subject (high vs standard glucose).pdf",sep=""),width=18,height=12,pointsize=20)
with(DwC, plot(logFC, -log10(P.Value), pch=20, col="gray", main="Diabetes with Complication (high vs. standard glucose) (p.value <0.05)", xlab="log2 fold change", ylab="-log10(P.val)", xlim=c(-1,1)))
with(subset(DwC, P.Value<pval.cutoff ), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(DwC, P.Value<pval.cutoff) , textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()
cat("number of genes (p-val <0.05): ",sum(DwC$P.Value<pval.cutoff))
write.csv(DwC_pval5, paste(Output,"Diabetes with complication Subject (high vs standard glucose) pval0.05.csv",sep=""))
rm(DwC, DwC_pval5)

# diabetes without complication subject (high vs standard glucose)
DwoC=topTable(efit.contrast,coef="DwoC", n=nrow(genes),adjust.method="BH",genelist=genes)
DwoC_pval5=DwoC[DwoC$P.Value<pval.cutoff,]
pdf(file = paste(Figs,"Diabetes without complication subject (high vs standard glucose).pdf",sep=""),width=18,height=12,pointsize=20)
with(DwoC, plot(logFC, -log10(P.Value), pch=20, col="gray", main="Diabetes without Complication (high vs. standard glucose) (p.value <0.05)", xlab="log2 fold change", ylab="-log10(P.val)", xlim=c(-1,1)))
with(subset(DwoC, P.Value<pval.cutoff ), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(DwoC, P.Value<pval.cutoff) , textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()
cat("number of genes (p-val <0.05): ",sum(DwoC$P.Value<pval.cutoff))
write.csv(DwoC_pval5, paste(Output,"Diabetes without complication Subject (high vs standard glucose) pval0.05.csv",sep=""))
rm(DwoC, DwoC_pval5)

# all subject (high vs standard glucose)
hg_sg=topTable(efit.contrast,coef="hg_sg", n=nrow(genes),adjust.method="BH",genelist=genes)
hg_sg_pval5=hg_sg[hg_sg$P.Value<pval.cutoff,]
pdf(file = paste(Figs,"All subject (high vs standard glucose).pdf",sep=""),width=18,height=12,pointsize=20)
with(hg_sg, plot(logFC, -log10(P.Value), pch=20, col="gray", main="All subject (high vs. standard glucose) (p.value <0.05)", xlab="log2 fold change", ylab="-log10(P.val)", xlim=c(-1,1)))
with(subset(hg_sg, P.Value<pval.cutoff ), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(hg_sg, P.Value<pval.cutoff) , textxy(logFC, -log10(P.Value), labs=geneSymbol, cex=.8))
dev.off()
cat("number of genes (p-val <0.05): ",sum(hg_sg$P.Value<pval.cutoff))
write.csv(hg_sg_pval5, paste(Output,"All Subject (high vs standard glucose) pval0.05.csv",sep=""))
rm(hg_sg, hg_sg_pval5)

# treatment effect based on t-test
hg.T<-data.frame(t(hg)); sg.T<-data.frame(t(sg))
hg.T$Name<-substr(rownames(hg.T),1,nchar(rownames(hg.T))-2); sg.T$Name<-substr(rownames(sg.T),1,nchar(rownames(sg.T))-2)
hgagg<-aggregate(hg.T[,-ncol(hg.T)],list(hg.T$Name), mean); sgagg<-aggregate(sg.T[,-ncol(sg.T)],list(sg.T$Name), mean)
rownames(hgagg)<-hgagg$Group.1; hgagg$Group.1<-NULL
hgmean<-apply(hgagg,2,as.numeric); rownames(hgmean)<-rownames(hgagg)
rownames(sgagg)<-sgagg$Group.1; sgagg$Group.1<-NULL
sgmean<-apply(sgagg,2,as.numeric); rownames(sgmean)<-rownames(sgagg)

pairedttest_pval <- foreach(i=1:ncol(hgmean), .combine=rbind) %dopar% {
  pvalue<-t.test(hgmean[,i],sgmean[,i], paired=TRUE)$p.value
}




# PCA and covariates

