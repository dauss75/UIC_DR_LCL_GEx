library(limma)
set.seed(1234)
HOME<-"/Users/sja517/Documents/results/IFN/MountSinaicolaboration2/"
OUTPUT_T <- "/Users/sja517/Documents/results/IFN/MountSinaicolaboration2/output_T/"
source(paste(HOME,"R-script/idat2lumibatch.R",sep=""))
RawData<-paste(HOME,"raw/",sep="")
setwd(RawData)

# read bead intensity values
idatFiles = dir(pattern="idat")
lumibatchOut<-idat2lumibatch(idatFiles)

if (require(lumiHumanIDMapping)){
    lumi.nuId <- addNuID2lumi(lumibatchOut, lib.mapping='lumiHumanIDMapping')
}

# dataMatrix1 <- exprs(lumibatchOut)
# plot(lumibatchOut, what='density')
# plot(lumibatchOut, what='boxplot', title=F)

idatBgFixed<-lumiB(lumi.nuId ) #background correction

lumi.T <- lumiT(idatBgFixed, method='log2') # Transfer the Illumina data to stabilize the variance
lumi.N <- lumiN(lumi.T) # quantile normalization
lumi.N.Q <- lumiQ(lumi.N) #quality control estimation after normalization

plot(lumi.N.Q, what='density')
plot(lumi.N.Q, what='boxplot')

# filtering probes
dataMatrix <- exprs(lumi.N.Q)

if (require(lumiHumanAll.db) & require(annotate)) {
#     dataMatrix <- dataMatrix[!is.na(getSYMBOL(rownames(dataMatrix), 'f')),]
    dataMatrix <- dataMatrix[!is.na(getSYMBOL(rownames(dataMatrix), 'lumiHumanAll.db')),]
}

probeList <- rownames(dataMatrix)
if (require(lumiHumanAll.db) & require(annotate)){
    geneSymbol<- getSYMBOL(probeList, 'lumiHumanAll.db')
    geneName<- sapply(lookUp(probeList, 'lumiHumanAll.db', 'GENENAME'), function(x) x[1])
}

## limma data analysis
## sort data by target table
targetID<-readTargets("Targets.txt")
targetID$ID<-paste(targetID$Sentrix_ID,"_",targetID$Sentrix_Position,sep="")
rownames(targetID)<-targetID$ID
dataMatrix_sorted<-dataMatrix[,targetID$ID]

## specify sample type by modified target table
targets<-readTargets("Targets2.txt")
ct<-factor(paste(targets$Condition, targets$Treatment,sep="."))
design <- model.matrix(~0+ct)
colnames(design) <- levels(ct)
corfit<-duplicateCorrelation(dataMatrix_sorted,design,block=targets$Subject)
fit <- lmFit(dataMatrix_sorted, design, block=targets$Subject, correlation=corfit$consensus)  
contrast.matrix <- makeContrasts(
                                 PV.IFN_Con = Diseased.IFN-Diseased.Con,
                                 PV.IFN_Normal.IFN   = Diseased.IFN-Normal.IFN,
                                 PV.Con_Normal.Con   = Diseased.Con-Normal.Con,
                                 Normal.IFN_Normal.Con     = Normal.IFN-Normal.Con,                                 
                                 levels=design)    
fit2<-contrasts.fit(fit, contrast.matrix)
fit2<-eBayes(fit2);  
fit2$genes <-data.frame(ID=probeList, geneSymbol=geneSymbol, geneName=geneName, stringsAsFactors=FALSE)

FDR<-8; 
pco5 <- 0.05 ## adj p-val cutoff
pco10 <- 0.1 ## adj p-val cutoff
lfc <- 1 ## absolute logFC cutoff


# PN_IFN <- PN_IFN[PN_IFN$adj.P.Val <= pco & abs(PN_IFN$logFC) >= lfc ,]
# PN_Con <- PN_Con[PN_Con$adj.P.Val <= pco & abs(PN_Con$logFC) >= lfc ,]
# up<-PN_IFN[order(-PN_IFN$logFC),]
# down<-PN_IFN[order(PN_IFN$logFC),]

## A: PV.IFN vs Normal.IFN
# PV.IFN_Normal.IFN <- topTable(fit2, coef="PV.IFN_Normal.IFN", number=nrow(fit2$genes), sort.by="logFC", genelist=fit2$genes, adjust.method="BH")
# PV.IFN_Normal.IFN_pval5 <- PV.IFN_Normal.IFN[PV.IFN_Normal.IFN$P.Value <= 0.05 & abs(PV.IFN_Normal.IFN$logFC) >= lfc ,]
# # write.csv(PV.IFN_Normal.IFN_pval5, file=paste(OUTPUT_T,"PV.IFN_Normal.IFN_pval5.csv",sep=""))
# PV.IFN_Normal.IFN_qval10 <- PV.IFN_Normal.IFN[PV.IFN_Normal.IFN$adj.P.Val <= pco10  & abs(PV.IFN_Normal.IFN$logFC) >= lfc ,]
# print(length(unique(PV.IFN_Normal.IFN_qval10$geneSymbol)))
# write.csv(PV.IFN_Normal.IFN_qval10, file=paste(OUTPUT_T,"PV.IFN_Normal.IFN_qval10.csv",sep=""))
# rm(PV.IFN_Normal.IFN, PV.IFN_Normal.IFN_pval5)

## B: PV.Con vs Normal.Con
PV.Con_Normal.Con <- topTable(fit2, coef="PV.Con_Normal.Con", number=nrow(fit2$genes), sort.by="logFC", genelist=fit2$genes, adjust.method="BH")
PV.Con_Normal.Con_pval5 <- PV.Con_Normal.Con[PV.Con_Normal.Con$P.Value <= 0.05 & abs(PV.Con_Normal.Con$logFC) >= lfc ,]
write.csv(PV.Con_Normal.Con_pval5, file=paste(OUTPUT_T,"PV.Con_Normal.Con_pval5.csv",sep=""))
PV.Con_Normal.Con_qval10 <- PV.Con_Normal.Con[PV.Con_Normal.Con$adj.P.Val <= pco10 & abs(PV.Con_Normal.Con$logFC) >= lfc ,]
write.csv(PV.Con_Normal.Con_qval10, file=paste(OUTPUT_T,"PV.Con_Normal.Con_qval10.csv",sep=""))
rm(PV.Con_Normal.Con, PV.Con_Normal.Con_pval5)

tmp3<-PV.Con_Normal.Con[!duplicated(PV.Con_Normal.Con$geneSymbol),]
with(tmp3, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main="Volcano plot", xlab="log2 fold change", ylab="-log10(adjusted p-value)", xlim=c(-4,4.5)))
with(subset(tmp3, adj.P.Val<.1 & abs(logFC)>1), points(logFC, -log10(adj.P.Val), pch=20, col="blue"), cex=1.1)

# Label points with the textxy function from the calibrate plot
library(calibrate)

# with(subset(tmp3, adj.P.Val<.2 & logFC > 1), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
# with(subset(tmp3, adj.P.Val<.2 & logFC < -1), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
with(subset(tmp3, adj.P.Val<.1 & logFC > 1), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
with(subset(tmp3, adj.P.Val<.1 & logFC < -1), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))

## C: PV.IFN vs PV.Con
PV.IFN_Con <- topTable(fit2, coef="PV.IFN_Con", number=nrow(fit2$genes), sort.by="logFC", genelist=fit2$genes, adjust.method="BH")
# PV.IFN_Con_pval5 <- PV.IFN_Con[PV.IFN_Con$P.Value <= 0.05 & abs(PV.IFN_Con$logFC) >= lfc ,]
# write.csv(PV.IFN_Con_pval5, file=paste(OUTPUT_T,"PV.IFN_Con_pval5.csv",sep=""))
PV.IFN_Con_qval10 <- PV.IFN_Con[PV.IFN_Con$adj.P.Val <= pco10 & abs(PV.IFN_Con$logFC) >= lfc ,]
write.csv(PV.IFN_Con_qval10, file=paste(OUTPUT_T,"PV.IFN_Con_qval10.csv",sep=""))
rm(PV.IFN_Con, PV.IFN_Con_pval5)

# tmp3<-PV.IFN_Con[!duplicated(PV.IFN_Con$geneSymbol),]
# with(tmp3, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main="Volcano plot", xlab="log2 fold change", ylab="-log10(adjusted p-value)", xlim=c(-4,4.5)))
# with(subset(tmp3, adj.P.Val<.1 & abs(logFC)>1), points(logFC, -log10(adj.P.Val), pch=20, col="blue"), cex=1.1)

# Label points with the textxy function from the calibrate plot
library(calibrate)

# with(subset(tmp3, adj.P.Val<.2 & logFC > 1), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
# with(subset(tmp3, adj.P.Val<.2 & logFC < -1), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))
with(subset(tmp3, adj.P.Val<.1 & logFC > 1), textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))


## D: Normal.IFN vs Normal.Con
Normal.IFN_Con <- topTable(fit2, coef="Normal.IFN_Normal.Con", number=nrow(fit2$genes), sort.by="logFC", genelist=fit2$genes, adjust.method="BH")
Normal.IFN_Con_pval5 <- Normal.IFN_Con[Normal.IFN_Con$P.Value <= 0.05 & abs(Normal.IFN_Con$logFC) >= lfc ,]
# write.csv(Normal.IFN_Con_pval5, file=paste(OUTPUT_T,"Normal.IFN_Con_pval5.csv",sep=""))
Normal.IFN_Con_qval10 <- Normal.IFN_Con[Normal.IFN_Con$adj.P.Val <= pco10 & abs(Normal.IFN_Con$logFC) >= lfc ,]
write.csv(Normal.IFN_Con_pval5, file=paste(OUTPUT_T,"Normal.IFN_Con_qval0.csv",sep=""))
rm(Normal.IFN_Con, Normal.IFN_Con_pval5)

## A - B (PV.IFN vs Normal.IFN) - (PV.Con vs Normal.Con)
# A_B<-intersect(PV.IFN_Normal.IFN_pval5$ID, PV.Con_Normal.Con_pval5$ID)
# A_B_result<-PV.IFN_Normal.IFN_pval5[-which(PV.IFN_Normal.IFN_pval5$ID %in% A_B),]

## A - C (PV.IFN vs Normal.IFN) - (PV.IFN vs PV.Con)
# A_C<-intersect(PV.IFN_Normal.IFN_pval5$ID, PV.IFN_Con_pval5$ID)
# A_C_result<-PV.IFN_Normal.IFN_pval5[-which(PV.IFN_Normal.IFN_pval5$ID %in% A_C),]

## A- D (PV.IFN vs Normal.IFN) - (Normal.IFN vs Normal.Con)
# A_D<-intersect(PV.IFN_Normal.IFN_pval5$ID, Normal.IFN_Con_pval5$ID)
# A_D_result<-PV.IFN_Normal.IFN_pval5[-which(PV.IFN_Normal.IFN_pval5$ID %in% A_D),]

## C - D (PV.IFN vs PV.Con) - (Normal.IFN vs Normal.Con)
tmp3<-PV.IFN_Con
with(tmp3, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main="Volcano plot", xlab="log2 fold change", ylab="-log10(adjusted p-value)", xlim=c(-4,4.5)))
with(subset(tmp3, adj.P.Val<.1 & abs(logFC)>1), points(logFC, -log10(adj.P.Val), pch=20, col="blue"), cex=1.1)

C_D<-intersect(PV.IFN_Con_qval10$geneSymbol, Normal.IFN_Con_qval10$geneSymbol)
C_D_result<-PV.IFN_Con_qval10[-which(PV.IFN_Con_qval10$geneSymbol %in% C_D),]
with(C_D_result, points(logFC, -log10(adj.P.Val), pch=20, col="red"), cex=1.1)
with(C_D_result, textxy(logFC, -log10(adj.P.Val), labs=geneSymbol, cex=.8))

write.csv(C_D_result, file=paste(OUTPUT_T,"PV.IFN_Con_Normal.IFN_Con_qval10.csv",sep=""))
## C - D - B (PV.IFN vs PV.Con) - (Normal.IFN vs Normal.Con) - (PV.Con vs Normal.Con)


rownames(PV.Con_Normal.Con_pval5)<-PV.Con_Normal.Con_pval5[,2]

temp<-read.csv("/Users/sja517/Documents/results/IFN/MountSinaicolaboration2/output_T/up-regulated_gene_list.csv")

## volcano plot
# tmp<-patientvsnormal_IFN[!(patientvsnormal_IFN$ID %in% sigNormal_IFN_Con$ID),]
# tmp2<-tmp[!(tmp$ID %in% sigPatient_IFN_Con$ID),]
# tmp3<-tmp2[!(tmp2$ID %in% sigPatient_vs_Normal_Con$ID),]
# up_tmp3<- tmp3[tmp3$logFC >0,]

with(tmp3, plot(logFC, -log10(adj.P.Val), pch=20, col="gray", main="Volcano plot", xlab="log2 fold change", ylab="-log10(adjusted p-value)", xlim=c(-4,4.5)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
# with(subset(tmp3, adj.P.Val < .05 ), points(logFC, pch=20, col="red"), cex=1.1)
# with(subset(tmp3, abs(logFC)>1), points(logFC, -log10(adj.P.Val), pch=20, col="orange"), cex=1.1)
with(subset(tmp3, adj.P.Val<.05 & abs(logFC)>1), points(logFC, -log10(adj.P.Val), pch=20, col="blue"), cex=1.1)


# Label points with the textxy function from the calibrate plot
library(calibrate)

with(subset(tmp3, adj.P.Val<.05 & logFC > 2.5), textxy(logFC, -log10(adj.P.Val), labs=Gene, cex=.8))
with(subset(tmp3, adj.P.Val<.05 & logFC < -1.5), textxy(logFC, -log10(adj.P.Val), labs=Gene, cex=.8))
with(subset(tmp3, adj.P.Val<.01 & logFC > 1), textxy(logFC, -log10(adj.P.Val), labs=Gene, cex=.8))