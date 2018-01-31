library(data.table)

sampleid<-fread("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/Match.csv")
rownames(sampleid)<-sampleid$GWAS_ID
# prepare expression data
exprs<-fread("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/edic_predixcan_lymphocyte_predicted_expression.txt", header=T)
exprs<-as.data.frame(exprs)
rownames(exprs)<-exprs$FID; exprs$FID<-NULL; exprs$IID<-NULL
exprs.t<-data.frame(t(exprs))
for (i in 1:length(rownames(exprs.t))){
  rownames(exprs.t)[i]<-unlist(strsplit(as.character(rownames(exprs.t)[i]),".",fixed=T))[1]
}

IFI27L1.lymphocyte.exprs<-data.frame(exprs.t[match("ENSG00000165948",rownames(exprs.t)),rownames(sampleid)])

library(splitstackshape) ## added
# sampleid<-expandRows(sampleid, count = 3, count.is.col = FALSE)  ## added

# genotypes
chr14<-fread("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/EDIC/Predix_chr14.dose", header=F)
chr14<-data.frame(chr14)
rownames(chr14)<-chr14[,2]
chr14[,1:6]<-NULL

header<-fread("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/EDIC/LASER_nohead.ped", header=F)
colnames(chr14)<-header$V2
chr14.variant<-chr14[,rownames(sampleid)]
# rm(chr14, chr14.variant)

library(ggplot2)
chr14.variant.rs1122361<-data.frame(t(chr14.variant["rs1122361",]))
chr14.variant.rs1122361$ID<-row.names(chr14.variant.rs1122361)
chr14.variant.rs1122361<-data.frame(expandRows(chr14.variant.rs1122361, count = 3, count.is.col = FALSE))  ## added
 
# qplot(t(chr14.variant.rs1122361),t(IFI27L1.lymphocyte.exprs),col=sampleid$type)

chr14.variant.rs12587898<-data.frame(t(chr14.variant["rs12587898",]))
chr14.variant.rs12587898$ID<-row.names(chr14.variant.rs12587898)
chr14.variant.rs12587898<-data.frame(expandRows(chr14.variant.rs12587898, count = 3, count.is.col = FALSE))  ## added

# qplot(t(chr14.variant.rs12587898),t(IFI27L1.lymphocyte.exprs),col=sampleid$type)

chr14.variant.rs748762<-data.frame(t(chr14.variant["rs748762",]))
chr14.variant.rs748762$ID<-row.names(chr14.variant.rs748762)
chr14.variant.rs748762<-data.frame(expandRows(chr14.variant.rs748762, count = 3, count.is.col = FALSE))  ## added

# qplot(t(chr14.variant.rs748762),t(IFI27L1.lymphocyte.exprs),col=sampleid$type)

chr14.variant.rs3814821<-data.frame(t(chr14.variant["rs3814821",]))
chr14.variant.rs3814821$ID<-row.names(chr14.variant.rs3814821)
chr14.variant.rs3814821<-data.frame(expandRows(chr14.variant.rs3814821, count = 3, count.is.col = FALSE))  ## added

# qplot(t(chr14.variant.rs3814821),t(IFI27L1.lymphocyte.exprs),col=sampleid$type)

## read RG_PDR_nPDR
HOME <- "/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/"  # change if needed
PhenotypeDir = paste(HOME,'data/phenotype/', sep='')
RData=paste(HOME,'rdata/',sep=''); 

source(paste(HOME,"script/Function.R",sep=""))
source(paste(HOME,"script/andrew_Func.R",sep=""))

dat<-local(get(load(file=paste(RData,"normalizedDataMatrix_filtered.RData",sep=""))))
probeList <- rownames(dat)  # probeID (nuID) gene annotation
if (require(lumiHumanAll.db) & require(annotate)){
  geneSymbol<- getSYMBOL(probeList, 'lumiHumanAll.db')
  geneName<- sapply(lookUp(probeList, 'lumiHumanAll.db', 'GENENAME'), function(x) x[1])
}
genes <-data.frame(ID=probeList, geneSymbol=geneSymbol, geneName=geneName, stringsAsFactors=FALSE)

# 2. No_PDR
subject<-"No_PDR"
DwoC<-dat[,grep("DwoC",colnames(dat))] # grep all normal samples
DwoC_sg<-DwoC[,grep("norm",colnames(DwoC))]; DwoC_hg<-DwoC[,grep("30mM",colnames(DwoC))];
DwoC_delta <- DwoC_hg-DwoC_sg
# DwoC_delta_avg<-prep_data(DwoC_delta)
# rm(DwoC, DwoC_sg, DwoC_hg, DwoC_delta)

# 3. PDR
DwC<-dat[,grep("DwC",colnames(dat))] # grep all normal samples
DwC_sg<-DwC[,grep("norm",colnames(DwC))]; DwC_hg<-DwC[,grep("30mM",colnames(DwC))];
DwC_delta <- DwC_hg-DwC_sg
# DwC_delta_avg<-prep_data(DwC_delta)
# rm(DwC, DwC_sg, DwC_hg, DwC_delta)
RG_PDR_nPDR<-cbind(DwC_delta,DwoC_delta)  ## added
# RG_PDR_nPDR<-cbind(DwC_delta_avg,DwoC_delta_avg)
IFI27L1.PDR_nPDR<-genes[match("IFI27L1",genes$geneSymbol),]
IFI27L1.PDR_nPDR.exprs<-RG_PDR_nPDR[IFI27L1.PDR_nPDR$ID,]
for (i in 1:length(colnames(IFI27L1.PDR_nPDR.exprs))){
  colnames(IFI27L1.PDR_nPDR.exprs)[i]<-unlist(strsplit(as.character(colnames(IFI27L1.PDR_nPDR.exprs)[i]),"_",fixed=T))[2]
}
IFI27L1.PDR_nPDR.exprs<-IFI27L1.PDR_nPDR.exprs[,order(as.double(colnames(IFI27L1.PDR_nPDR.exprs)))]
rownames(IFI27L1.PDR_nPDR.exprs)<-"IFI27L1"
# rownames(sampleid)<-sampleid$PATIENT

# library(splitstackshape) ## added
sampleid<-expandRows(sampleid, count = 3, count.is.col = FALSE)  ## added


# IFI27L1.PDR_nPDR.exprs<-IFI27L1.PDR_nPDR.exprs[,sampleid$PATIENT]
colnames(IFI27L1.PDR_nPDR.exprs)<-sampleid$GWAS_ID
# cor(t(chr14.variant.rs1122361),t(IFI27L1.PDR_nPDR.exprs))
library(ggplot2)
setEPS()
postscript(file = "/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/rs1122361_all.eps")
qplot(chr14.variant.rs1122361[,1],t(IFI27L1.PDR_nPDR.exprs),col=sampleid$type, main="IFI27L1 vs rs11222361", ylab="Gene Expression (RG_PDR_nPDR)", xlab="Imputed genotype (rs1122361)") +labs(col='Group') +theme(plot.title = element_text(hjust = 0.5))
dev.off()

setEPS()
postscript(file = "/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/rs12587898_all.eps")
qplot(chr14.variant.rs12587898[,1],t(IFI27L1.PDR_nPDR.exprs),col=sampleid$type, main="IFI27L1 vs rs12587898", ylab="Gene Expression (RG_PDR_nPDR)", xlab="Imputed genotype (rs12587898)") +labs(col='Group') +theme(plot.title = element_text(hjust = 0.5))
dev.off()

setEPS()
postscript(file = "/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/rs748762_all.eps")
qplot(chr14.variant.rs748762[,1],t(IFI27L1.PDR_nPDR.exprs),col=sampleid$type, main="IFI27L1 vs rs748762", ylab="Gene Expression (RG_PDR_nPDR)", xlab="Imputed genotype (rs748762)") +labs(col='Group') +theme(plot.title = element_text(hjust = 0.5))
dev.off()

setEPS()
postscript(file = "/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/rs3814821_all.eps")
qplot(chr14.variant.rs3814821[,1],t(IFI27L1.PDR_nPDR.exprs),col=sampleid$type, main="IFI27L1 vs rs3814821", ylab="Gene Expression (RG_PDR_nPDR)", xlab="Imputed genotype (rs3814821)") +labs(col='Group') +theme(plot.title = element_text(hjust = 0.5))
dev.off()
