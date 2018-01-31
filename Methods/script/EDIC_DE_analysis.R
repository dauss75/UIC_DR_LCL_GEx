library(oligo); library(limma)
library(annotate)

data="/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/data/02-22-2017/EDIC_lymphocyte_affy_data"
celFiles <- list.celfiles(data, full.names=TRUE)
affy.data <- read.celfiles(celFiles)
eset <- rma(affy.data)
# exprSet<-exprs(eset.rma)

# library(hgu95av2.db)
# library(hugene10sttranscriptcluster.db)

ID <- featureNames(eset)
# Symbol <- getSYMBOL(ID,"hgu95av2.db")
Symbol <- getSYMBOL(ID,"hugene10sttranscriptcluster.db")
fData(eset) <- data.frame(Symbol=Symbol)

pheno<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/data/02-22-2017/lymphocyte.pheno.data.corrected.csv")
treatments<-factor(pheno$Group)

design<-model.matrix(~0+treatments)
colnames(design)<-levels(treatments)
cont.matrix<-makeContrasts(levels=colnames(design), DE=case-control)
fit<-lmFit(eset,design)
fit2<-contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2)

# annodb <- "hugene10sttranscriptcluster.db"
# ID<-row.names(eb.ls)

# genes<-getSYMBOL(ID,"hugene10sttranscriptcluster.db")


x1=topTable(fit2, coef="DE", n=length(ID))
x1<-x1[complete.cases(x1), ]
write.csv(x1,"/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/5.indenepdnet_data/EDIC_lymphocyte_affy_DE.csv", row.names=F)
