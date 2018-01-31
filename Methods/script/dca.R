library(DGCA); library(limma)
data(darmanis); data(design_mat)

HOME <- "/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/"  # change if needed
PhenotypeDir = paste(HOME,'data/phenotype/', sep='')
RData=paste(HOME,'rdata/',sep=''); 

source(paste(HOME,"pca/pca.R",sep=""))
source(paste(HOME,"script/Function.R",sep=""))

dat<-local(get(load(file=paste(RData,"normalizedDataMatrix_filtered.RData",sep=""))))
probeList <- rownames(dat)  # probeID (nuID) gene annotation
if (require(lumiHumanAll.db) & require(annotate)){
  geneSymbol<- getSYMBOL(probeList, 'lumiHumanAll.db')
  geneName<- sapply(lookUp(probeList, 'lumiHumanAll.db', 'GENENAME'), function(x) x[1])
}
genes <-data.frame(ID=probeList, geneSymbol=geneSymbol, geneName=geneName, stringsAsFactors=FALSE)

subject="delta_PDR_vs_No_PDR"
targets<-readTargets(paste(PhenotypeDir,subject,"_target.txt", sep=''))
group<-as.factor(targets$Group)
design_mat <- model.matrix(~0+group)
colnames(design_mat)<-levels(group)


# 2. No_PDR
subject<-"No_PDR"
DwoC<-dat[,grep("DwoC",colnames(dat))] # grep all normal samples
DwoC_sg<-DwoC[,grep("norm",colnames(DwoC))]; DwoC_hg<-DwoC[,grep("30mM",colnames(DwoC))];
DwoC_delta <- DwoC_hg-DwoC_sg

# 3. PDR
DwC<-dat[,grep("DwC",colnames(dat))] # grep all normal samples
DwC_sg<-DwC[,grep("norm",colnames(DwC))]; DwC_hg<-DwC[,grep("30mM",colnames(DwC))];
DwC_delta <- DwC_hg-DwC_sg

exprs<-cbind(DwC_delta,DwoC_delta)
rownames(exprs)<-make.names(genes$geneSymbol, unique=TRUE)

#filter the genes
library(matrixStats, quietly = TRUE)
exprs_mean_filtered = filterGenes(exprs, filterTypes = "central", filterCentralType = "median", filterCentralPercentile = 0.3)
exprs_cv_filtered = filterGenes(exprs_mean_filtered, filterTypes = "dispersion", filterDispersionType = "cv", filterDispersionPercentile = 0.3)

cor_res = getCors(inputMat = exprs_cv_filtered, design = design_mat) 
dcPairs_res = pairwiseDCor(cor_res, compare = c("DwC", "DwoC"))
str(dcPairs_res)
dd_pairs = dcTopPairs(dcPairs_res, nPairs = 1000, classify = TRUE, adjust = "BH")
write.csv(dd_pairs,"/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/DCA/DCA_top1000pairs.csv",row.names = F)


# ddcor_res1 = ddcorAll(inputMat = exprs, design = design_mat, compare = c("DwC", "DwoC"), adjust = "perm", nPerm = 10, splitSet = "IFI27L1")
# head(ddcor_res1)
# write.csv(ddcor_res1,"/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/DCA/DCA_IFI27L1.csv",row.names = F)

library("Cairo")
library("ggplot2")
# setEPS()
# postscript(file = "/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/DCA/NUP98_IFI27L1_Correlation.eps")
ggsave("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/DCA/top1_NUP98_IFI27L1.eps", device=cairo_ps)
plotCors(inputMat = exprs, design = design_mat, compare = c("DwC", "DwoC"), geneA = "NUP98", geneB = "IFI27L1")
dev.off()

ggsave("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/DCA/top2_PALLD_IFI27L1.eps", device=cairo_ps)
plotCors(inputMat = exprs, design = design_mat, compare = c("DwC", "DwoC"), geneA = "PALLD", geneB = "IFI27L1")
dev.off()

ggsave("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/DCA/top3_ANAPC11_IFI27L1.eps", device=cairo_ps)
plotCors(inputMat = exprs, design = design_mat, compare = c("DwC", "DwoC"), geneA = "ANAPC11", geneB = "IFI27L1")
dev.off()


plotVals(inputMat = exprs, design = design_mat, compare = c("DwC", "DwoC"), gene = "IFI27L1")

ddcor_res2 = ddcorAll(inputMat = exprs, design = design_mat, compare = c("DwC", "DwoC"), adjust = "perm", nPerm = 10, splitSet = "IFI27L1.1")
head(ddcor_res2)
write.csv(ddcor_res2,"/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/DCA/DCA_IFI27L1.1.csv",row.names = F)

ggsave("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/DCA/top1_DAGLB_IFI27L1.1.eps", device=cairo_ps)
plotCors(inputMat = exprs, design = design_mat, compare = c("DwC", "DwoC"), geneA = "DAGLB", geneB = "IFI27L1.1")
dev.off()

ggsave("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/DCA/top2_SLC27A3_IFI27L1.1.eps", device=cairo_ps)
plotCors(inputMat = exprs, design = design_mat, compare = c("DwC", "DwoC"), geneA = "SLC27A3", geneB = "IFI27L1.1")
dev.off()

ggsave("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/DCA/top3_AGBL5_IFI27L1.1.eps", device=cairo_ps)
plotCors(inputMat = exprs, design = design_mat, compare = c("DwC", "DwoC"), geneA = "AGBL5", geneB = "IFI27L1.1")
dev.off()

## RG_all
rownames(dat)<-make.names(genes$geneSymbol, unique=TRUE)
subject="all_replicate"
targets<-readTargets(paste(PhenotypeDir,"hg_sg_",subject,"_target.txt", sep=''))
Treat <- factor(targets$Treatment,levels=c("C","T"))
# Replicates <- factor(targets$rep)
design <- model.matrix(~0+Treat)

#filter the genes
library(matrixStats, quietly = TRUE)
exprs_mean_filtered = filterGenes(dat, filterTypes = "central", filterCentralType = "median", filterCentralPercentile = 0.3)
exprs_cv_filtered = filterGenes(exprs_mean_filtered, filterTypes = "dispersion", filterDispersionType = "cv", filterDispersionPercentile = 0.3)

cor_res = getCors(inputMat = exprs_cv_filtered, design = design) 
dcPairs_res = pairwiseDCor(cor_res, compare = c("TreatT", "TreatC"))
str(dcPairs_res)
dd_pairs = dcTopPairs(dcPairs_res, nPairs=10000, classify = TRUE, adjust = "BH")
write.csv(dd_pairs,"/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/DCA/DCA_RGAll_top10000pairs.csv",row.names = F)


