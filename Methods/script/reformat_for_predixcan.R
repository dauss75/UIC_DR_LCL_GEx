library(data.table)

HOME="/Users/sjung/Project/GlobusGenomics/UIC/predixcan/EDIC/"
dose_files = list.files(path=HOME, pattern='dose', full.names=TRUE)
for (i in 1:length(dose_files)){
  tt<-fread(dose_files[i], header=F)
  # tt<-na.omit(tt)
  # tt<-as.data.frame(append(tt, list( MAF = 0.05), after = 5))
  fwrite(tt,paste(HOME,"chr",i,".new.txt",sep=""),sep=" ", row.names=F, col.names = F)
} 

sample<-fread("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/genotype/samples.txt")
header<-fread("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/samples/samples.txt")
phenotype<-na.omit(sample)
colnames(phenotype)<-colnames(header)
fwrite(phenotype,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/phenotype/phenotype.txt", sep="\t")

# sample<-fread("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/GOKIND/LASER2.ped")
sample<-fread("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/EDIC/LASER.ped")

## post-processing
edic_expr<-fread("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/edic_predixcan_lymphocytes_predicted_expression.txt")
edic_assoc<-fread("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/edic_predixcan_lymphocytes_association.txt")
edic_assoc<-na.omit(edic_assoc)

library(qvalue)
qobj <- qvalue(p = edic_assoc$p)
summary(qobj)

#           <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1   <1
#p-value        1      3    27     82   182  353 3420
#q-value        0      0     1      1     1    1 3420
#local FDR      0      0     1      1     1    1    2

library("biomaRt")

mart <- useDataset(dataset = "hsapiens_gene_ensembl",
                   mart = useMart("ensembl",
                                  host    = "www.ensembl.org"))

attribs = c("ensembl_gene_id", "entrezgene", "hgnc_symbol",
            "chromosome_name", "start_position","end_position")


gene_set<-NULL
for (i in 1:length(edic_assoc$gene)){
  gene_set[i]<-unlist(strsplit(as.character(edic_assoc$gene[i]),".",fixed=T))[1]
}
edic_assoc$gene<-gene_set

resultTable <- getBM(attributes = attribs, filters = "ensembl_gene_id",
                     values = edic_assoc$gene, mart = mart)

resultTable2<-resultTable[!resultTable$hgnc_symbol=="",]
tt<-merge(edic_assoc,resultTable2,by.x="gene",by.y="ensembl_gene_id")
tt2<-tt[order(tt$p),]
tt2$entrezgene<-NULL
edic_pred_assoc_result<-tt2[!duplicated(tt2),]

write.csv(edic_pred_assoc_result,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/edic_predixcan_lymphocyte_assoc_result.csv",row.names=F)

## find overlap with GTeX
# GTex_Home<-"/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/results/"
# 
# GTex_HG_PDR_nPDR<-read.csv(paste(GTex_Home,"486snp_gene_in_HG_PDR_nPDR.csv",sep=""))
# GTex_SG_PDR_nPDR<-read.csv(paste(GTex_Home,"486snp_gene_in_SG_PDR_nPDR.csv",sep=""))
# GTex_RG_PDR_nPDR<-read.csv(paste(GTex_Home,"486snp_gene_in_RG_PDR_nPDR.csv",sep=""))
# 
# GTex_HG_PDR_nPDR_predixcan<-merge(GTex_HG_PDR_nPDR,edic_pred_assoc_result, by.x="geneSymbol", by.y="hgnc_symbol")

# intersect - edic vs (hg, sg, rg) pdr_npdr
edic_assoc_result<-fread("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/edic_predixcan_lymphocyte_assoc_result.csv") #3338
edic_assoc_result_p5<-edic_assoc_result[edic_assoc_result$p<=0.05,]; #178
edic_assoc_result_p5$beta<-NULL;edic_assoc_result_p5$t<-NULL; edic_assoc_result_p5$`se(beta)`<-NULL

RG_PDR_nPDR<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/4.Glucose_Response_analyses/difference_between_group/delta_DwC_DwoC_pvalue_all_2017-09-13.csv")
RG_PDR_nPDR$ID<-NULL; RG_PDR_nPDR$t<-NULL; RG_PDR_nPDR$adj.P.Val<-NULL; RG_PDR_nPDR$B<-NULL
RG_PDR_nPDR_p5<-RG_PDR_nPDR[RG_PDR_nPDR$P.Value<=0.05,] #654
predixcan_p5_RG_PDR_nPDR_p5<-merge(RG_PDR_nPDR_p5, edic_assoc_result_p5, by.x="geneSymbol", by.y="hgnc_symbol")
predixcan_p5_RG_PDR_nPDR_p5<-predixcan_p5_RG_PDR_nPDR_p5[order(predixcan_p5_RG_PDR_nPDR_p5$P.Value),]
predixcan_p5_RG_PDR_nPDR_p5<-predixcan_p5_RG_PDR_nPDR_p5[!duplicated(predixcan_p5_RG_PDR_nPDR_p5$geneSymbol),]


HG_PDR_nPDR<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/3.Between_Group_Analysis/HG/DwC_vs_DwoC/hg_DwC_DwoC_PC1_GrowthRate_all_2017-09-13.csv")
HG_PDR_nPDR$ID<-NULL; HG_PDR_nPDR$t<-NULL; HG_PDR_nPDR$adj.P.Val<-NULL; HG_PDR_nPDR$B<-NULL
HG_PDR_nPDR_p5<-HG_PDR_nPDR[HG_PDR_nPDR$P.Value<=0.05,] #677
predixcan_p5_HG_PDR_nPDR_p5<-merge(HG_PDR_nPDR_p5, edic_assoc_result_p5, by.x="geneSymbol", by.y="hgnc_symbol") #8
predixcan_p5_HG_PDR_nPDR_p5<-predixcan_p5_HG_PDR_nPDR_p5[order(predixcan_p5_HG_PDR_nPDR_p5$P.Value),]
predixcan_p5_HG_PDR_nPDR_p5<-predixcan_p5_HG_PDR_nPDR_p5[!duplicated(predixcan_p5_HG_PDR_nPDR_p5$geneSymbol),]

SG_PDR_nPDR<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/3.Between_Group_Analysis/SG/DwC_vs_DwoC/sg_DwC_DwoC_PC1_GrowthRate_all_2017-09-13.csv")
SG_PDR_nPDR$ID<-NULL; SG_PDR_nPDR$t<-NULL; SG_PDR_nPDR$adj.P.Val<-NULL; SG_PDR_nPDR$B<-NULL
SG_PDR_nPDR_p5<-SG_PDR_nPDR[HG_PDR_nPDR$P.Value<=0.05,] #677
predixcan_p5_SG_PDR_nPDR_p5<-merge(SG_PDR_nPDR_p5, edic_assoc_result_p5, by.x="geneSymbol", by.y="hgnc_symbol") #6
predixcan_p5_SG_PDR_nPDR_p5<-predixcan_p5_SG_PDR_nPDR_p5[order(predixcan_p5_SG_PDR_nPDR_p5$P.Value),]
predixcan_p5_SG_PDR_nPDR_p5<-predixcan_p5_SG_PDR_nPDR_p5[!duplicated(predixcan_p5_SG_PDR_nPDR_p5$geneSymbol),]


write.csv(predixcan_p5_RG_PDR_nPDR_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/predixcan_lymphocyte_p5_RG_PDR_nPDR_p5.csv",row.names = FALSE)
write.csv(predixcan_p5_HG_PDR_nPDR_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/predixcan_lymphocyte_p5_HG_PDR_nPDR_p5.csv",row.names = FALSE)
write.csv(predixcan_p5_SG_PDR_nPDR_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/predixcan_lymphocyte_p5_SG_PDR_nPDR_p5.csv",row.names = FALSE)


# intersect - 487 gtex, edic, (hg, sg, rg) pdr_npdr
predixcan_p5_RG_PDR_nPDR_p5<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/predixcan_lymphocyte_p5_RG_PDR_nPDR_p5.csv")
predixcan_p5_HG_PDR_nPDR_p5<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/predixcan_lymphocyte_p5_HG_PDR_nPDR_p5.csv")
predixcan_p5_SG_PDR_nPDR_p5<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/predixcan_lymphocyte_p5_SG_PDR_nPDR_p5.csv")
gtex_487<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/486snp_vs_GTex_nonredundant.csv")

gtex_predixcan_p5_RG_PDR_nPDR_p5<-merge(predixcan_p5_RG_PDR_nPDR_p5, gtex_487, by.x="geneSymbol", by.y="hgnc_symbol") 
gtex_predixcan_p5_HG_PDR_nPDR_p5<-merge(predixcan_p5_HG_PDR_nPDR_p5, gtex_487, by.x="geneSymbol", by.y="hgnc_symbol") 
gtex_predixcan_p5_SG_PDR_nPDR_p5<-merge(predixcan_p5_SG_PDR_nPDR_p5, gtex_487, by.x="geneSymbol", by.y="hgnc_symbol") 

write.csv(gtex_predixcan_p5_RG_PDR_nPDR_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/gtex487snps_edic_predixcan_lymphocyte_p5_RG_PDR_nPDR_p5.csv",row.names = FALSE)
write.csv(gtex_predixcan_p5_HG_PDR_nPDR_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/gtex487snps_edic_predixcan_lymphocyte_p5_HG_PDR_nPDR_p5.csv",row.names = FALSE)
write.csv(gtex_predixcan_p5_SG_PDR_nPDR_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/gtex487snps_edic_predixcan_lymphocyte_p5_SG_PDR_nPDR_p5.csv",row.names = FALSE)

# intersect - gsea edic vs (hg, sg, rg) pdr_npdr, IPA (done)

# DCA

# check with Andrew 


