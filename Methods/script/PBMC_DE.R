library(oligo); library(limma)
library(annotate)

data="/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/data/02-22-2017/PBMC_Microarray/20100127.PBMC.csv"
pbmc<-read.csv(data)
pbmc$X<-NULL
genesymbol<-data.frame(pbmc$GENE_SYMBOL); pbmc$GENE_SYMBOL<-NULL;

pheno<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/data/02-22-2017/pbmc.pheno.csv")
treatments<-factor(pheno$Group)

design<-model.matrix(~0+treatments)
colnames(design)<-levels(treatments)
cont.matrix<-makeContrasts(levels=colnames(design), DE=case-control)
fit<-lmFit(pbmc,design)
fit2<-contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2)

# annodb <- "hugene10sttranscriptcluster.db"
# ID<-row.names(eb.ls)

# genes<-getSYMBOL(ID,"hugene10sttranscriptcluster.db")


x1=topTable(fit2, coef="DE", genelist=genesymbol,n=nrow(genesymbol))
x1<-x1[complete.cases(x1), ]
write.csv(x1,"/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/5.indenepdnet_data/PBMC_DE.csv", row.names=F)
