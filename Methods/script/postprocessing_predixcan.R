library(data.table)

tissue="lymphocyte"  #cerebellum, heart_atrial_appendage, heart_left_ventricle, liver, lymphocyte, whole_blood

edic_expr<-fread(paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/edic_predixcan_",tissue,"_predicted_expression.txt",sep=""))
edic_assoc<-fread(paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/edic_predixcan_",tissue,"_association.txt",sep=""))
edic_assoc<-na.omit(edic_assoc)

library(qvalue)
qobj <- qvalue(p = edic_assoc$p)
sink(file=paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/edic_predixcan_",tissue,"_association_pval_summary.txt",sep=""))
summary(qobj)
sink()

qobj_data<-data.frame(gene=edic_assoc$gene,pval=qobj$pvalues, FDR=qobj$qvalues)

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
fdr_summary<-merge(qobj_data,resultTable2,by.x="gene",by.y="ensembl_gene_id")
fdr_summary2<-data.frame(gene=fdr_summary$hgnc_symbol, pval=fdr_summary$pval,qval=fdr_summary$FDR)
tt<-merge(edic_assoc,resultTable2,by.x="gene",by.y="ensembl_gene_id")
tt2<-tt[order(tt$p),]
tt2$entrezgene<-NULL
edic_pred_assoc_result<-tt2[!duplicated(tt2),]

write.csv(edic_pred_assoc_result,paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/edic_predixcan_",tissue,"_assoc_result.csv",sep=""),row.names=F)
write.csv(fdr_summary2, paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/edic_predixcan_",tissue,"_assoc_fdr_result.csv",sep=""),row.names=F)



edic_assoc_result<-fread(paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/edic_predixcan_",tissue,"_assoc_result.csv",sep="")) #3338
edic_assoc_result_p5<-edic_assoc_result[edic_assoc_result$p<=0.05,]; #130
edic_assoc_result_p5$beta<-NULL;edic_assoc_result_p5$t<-NULL; edic_assoc_result_p5$`se(beta)`<-NULL

#------------------------------------------------------------#
#      intersect - edic vs PBMC, lympocyte, monocyte         #
#------------------------------------------------------------#
EDIC_lymphocyte<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/5.indenepdnet_data/EDIC_lymphocyte_affy_DE.csv")
EDIC_lymphocyte_p5<-EDIC_lymphocyte[EDIC_lymphocyte$P.Value<=0.05,];
EDIC_monocyte<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/5.indenepdnet_data/EDIC_monocyte_affy_DE.csv")
EDIC_monocyte_p5<-EDIC_monocyte[EDIC_monocyte$P.Value<=0.05,];
PBMC<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/5.indenepdnet_data/PBMC_DE.csv")
PBMC_p5<-PBMC[PBMC$P.Value<=0.05,];

predixcan_p5_EDIC_lymphocyte_p5<-merge(EDIC_lymphocyte_p5, edic_assoc_result_p5, by.x="Symbol", by.y="hgnc_symbol")
predixcan_p5_EDIC_monocyte_p5<-merge(EDIC_monocyte_p5, edic_assoc_result_p5, by.x="Symbol", by.y="hgnc_symbol")
predixcan_p5_PBMC_p5<-merge(PBMC_p5, edic_assoc_result_p5, by.x="pbmc.GENE_SYMBOL", by.y="hgnc_symbol")

write.csv(predixcan_p5_EDIC_lymphocyte_p5,paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/intersect_edic_predixcan_",tissue,"_p5_EDIC_lymphocyte_p5.csv",sep=""),row.names = FALSE)
write.csv(predixcan_p5_EDIC_monocyte_p5,paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/intersect_edic_predixcan_",tissue,"_p5_EDIC_monocyte_p5.csv",sep=""),row.names = FALSE)
write.csv(predixcan_p5_PBMC_p5,paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/intersect_edic_predixcan_",tissue,"_p5_PBMC_p5.csv",sep=""),row.names = FALSE)

#------------------------------------------------------------#
#           intersect - edic vs (hg, sg, rg) pdr_npdr        #
#------------------------------------------------------------#
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


write.csv(predixcan_p5_RG_PDR_nPDR_p5,paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/intersect_edic_predixcan_",tissue,"_p5_RG_PDR_nPDR_p5.csv",sep=""),row.names = FALSE)
write.csv(predixcan_p5_HG_PDR_nPDR_p5,paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/intersect_edic_predixcan_",tissue,"_p5_HG_PDR_nPDR_p5.csv",sep=""),row.names = FALSE)
write.csv(predixcan_p5_SG_PDR_nPDR_p5,paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/intersect_edic_predixcan_",tissue,"_p5_SG_PDR_nPDR_p5.csv",sep=""),row.names = FALSE)


# intersect - 487 gtex, edic, (hg, sg, rg) pdr_npdr
# predixcan_p5_RG_PDR_nPDR_p5<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/predixcan_lymphocyte_p5_RG_PDR_nPDR_p5.csv")
# predixcan_p5_HG_PDR_nPDR_p5<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/predixcan_lymphocyte_p5_HG_PDR_nPDR_p5.csv")
# predixcan_p5_SG_PDR_nPDR_p5<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/predixcan_lymphocyte_p5_SG_PDR_nPDR_p5.csv")
#------------------------------
# gtex v6
#------------------------------
gtex_487<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/486snp_vs_GTex_nonredundant.csv")

gtex_predixcan_p5_RG_PDR_nPDR_p5<-merge(predixcan_p5_RG_PDR_nPDR_p5, gtex_487, by.x="geneSymbol", by.y="hgnc_symbol") 
gtex_predixcan_p5_HG_PDR_nPDR_p5<-merge(predixcan_p5_HG_PDR_nPDR_p5, gtex_487, by.x="geneSymbol", by.y="hgnc_symbol") 
gtex_predixcan_p5_SG_PDR_nPDR_p5<-merge(predixcan_p5_SG_PDR_nPDR_p5, gtex_487, by.x="geneSymbol", by.y="hgnc_symbol") 

write.csv(gtex_predixcan_p5_RG_PDR_nPDR_p5,paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/intersect_gtex487snps_edic_predixcan_",tissue,"_p5_RG_PDR_nPDR_p5.csv",sep=""),row.names = FALSE)
write.csv(gtex_predixcan_p5_HG_PDR_nPDR_p5,paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/intersect_gtex487snps_edic_predixcan_",tissue,"_p5_HG_PDR_nPDR_p5.csv",sep=""),row.names = FALSE)
write.csv(gtex_predixcan_p5_SG_PDR_nPDR_p5,paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/intersect_gtex487snps_edic_predixcan_",tissue,"_p5_SG_PDR_nPDR_p5.csv",sep=""),row.names = FALSE)

predixcan_cerebellum<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/cerebellum/edic_predixcan_cerebellum_assoc_result.csv")
predixcan_cerebellum_p5<-predixcan_cerebellum[predixcan_cerebellum$p<=0.05,]; rm(predixcan_cerebellum)

predixcan_heart_atrial_appendage<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/heart_atrial_appendage/edic_predixcan_heart_atrial_appendage_assoc_result.csv")
predixcan_heart_atrial_appendage_p5<-predixcan_heart_atrial_appendage[predixcan_heart_atrial_appendage$p<=0.05,]; rm(predixcan_heart_atrial_appendage)

predixcan_heart_left_ventricle<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/heart_left_ventricle/edic_predixcan_heart_left_ventricle_assoc_result.csv")
predixcan_heart_left_ventricle_p5<-predixcan_heart_left_ventricle[predixcan_heart_left_ventricle$p<=0.05,]; rm(predixcan_heart_left_ventricle)

predixcan_liver<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/liver/edic_predixcan_liver_assoc_result.csv")
predixcan_liver_p5<-predixcan_liver[predixcan_liver$p<=0.05,]; rm(predixcan_liver)

predixcan_lymphocyte<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/edic_predixcan_lymphocyte_assoc_result.csv")
predixcan_lymphocyte_p5<-predixcan_lymphocyte[predixcan_lymphocyte$p<=0.05,]; rm(predixcan_lymphocyte)

predixcan_whole_blood<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/whole_blood/edic_predixcan_whole_blood_assoc_result.csv")
predixcan_whole_blood_p5<-predixcan_whole_blood[predixcan_whole_blood$p<=0.05,]; rm(predixcan_whole_blood)

gtex487_predixcan_cerebellum_p5 <-merge(predixcan_cerebellum_p5, gtex_487, by="hgnc_symbol"); rm(predixcan_cerebellum_p5) 
gtex487_predixcan_heart_atrial_appendage_p5 <-merge(predixcan_heart_atrial_appendage_p5, gtex_487, by="hgnc_symbol"); rm(predixcan_heart_atrial_appendage_p5) 
gtex487_predixcan_heart_left_ventricle_p5 <-merge(predixcan_heart_left_ventricle_p5, gtex_487, by="hgnc_symbol"); rm(predixcan_heart_left_ventricle_p5) 
gtex487_predixcan_liver_p5 <-merge(predixcan_liver_p5, gtex_487, by="hgnc_symbol"); rm(predixcan_liver_p5) 
gtex487_predixcan_lymphocyte_p5 <-merge(predixcan_lymphocyte_p5, gtex_487, by="hgnc_symbol"); rm(predixcan_lymphocyte_p5) 
gtex487_predixcan_whole_blood_p5 <-merge(predixcan_whole_blood_p5, gtex_487, by="hgnc_symbol"); rm(predixcan_whole_blood_p5) 

write.csv(gtex487_predixcan_cerebellum_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/cerebellum/intersect_gtex487_predixcan_cerebellum_p5.csv",row.names = FALSE)
write.csv(gtex487_predixcan_heart_atrial_appendage_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/heart_atrial_appendage/intersect_gtex487_predixcan_heart_atrial_appendage_p5.csv",row.names = FALSE)
write.csv(gtex487_predixcan_heart_left_ventricle_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/heart_left_ventricle/intersect_gtex487_predixcan_heart_left_ventricle_p5.csv",row.names = FALSE)
write.csv(gtex487_predixcan_liver_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/liver/intersect_gtex487_predixcan_liver_p5.csv",row.names = FALSE)
write.csv(gtex487_predixcan_lymphocyte_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/intersect_gtex487_predixcan_lymphocyte_p5.csv",row.names = FALSE)
write.csv(gtex487_predixcan_whole_blood_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/whole_blood/intersect_gtex487_predixcan_whole_blood_p5.csv",row.names = FALSE)

#------------------------------
## gtex v7
#------------------------------

gtex_487<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/486snp_vs_GTex_v7_nonredundant.csv")

gtex_predixcan_p5_RG_PDR_nPDR_p5<-merge(predixcan_p5_RG_PDR_nPDR_p5, gtex_487, by.x="geneSymbol", by.y="hgnc_symbol") 
gtex_predixcan_p5_HG_PDR_nPDR_p5<-merge(predixcan_p5_HG_PDR_nPDR_p5, gtex_487, by.x="geneSymbol", by.y="hgnc_symbol") 
gtex_predixcan_p5_SG_PDR_nPDR_p5<-merge(predixcan_p5_SG_PDR_nPDR_p5, gtex_487, by.x="geneSymbol", by.y="hgnc_symbol") 

write.csv(gtex_predixcan_p5_RG_PDR_nPDR_p5,paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/intersect_gtex_v7_487snps_edic_predixcan_",tissue,"_p5_RG_PDR_nPDR_p5.csv",sep=""),row.names = FALSE)
write.csv(gtex_predixcan_p5_HG_PDR_nPDR_p5,paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/intersect_gtex_v7_487snps_edic_predixcan_",tissue,"_p5_HG_PDR_nPDR_p5.csv",sep=""),row.names = FALSE)
write.csv(gtex_predixcan_p5_SG_PDR_nPDR_p5,paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/intersect_gtex_v7_487snps_edic_predixcan_",tissue,"_p5_SG_PDR_nPDR_p5.csv",sep=""),row.names = FALSE)

predixcan_cerebellum<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/cerebellum/edic_predixcan_cerebellum_assoc_result.csv")
predixcan_cerebellum_p5<-predixcan_cerebellum[predixcan_cerebellum$p<=0.05,]; rm(predixcan_cerebellum)

predixcan_heart_atrial_appendage<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/heart_atrial_appendage/edic_predixcan_heart_atrial_appendage_assoc_result.csv")
predixcan_heart_atrial_appendage_p5<-predixcan_heart_atrial_appendage[predixcan_heart_atrial_appendage$p<=0.05,]; rm(predixcan_heart_atrial_appendage)

predixcan_heart_left_ventricle<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/heart_left_ventricle/edic_predixcan_heart_left_ventricle_assoc_result.csv")
predixcan_heart_left_ventricle_p5<-predixcan_heart_left_ventricle[predixcan_heart_left_ventricle$p<=0.05,]; rm(predixcan_heart_left_ventricle)

predixcan_liver<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/liver/edic_predixcan_liver_assoc_result.csv")
predixcan_liver_p5<-predixcan_liver[predixcan_liver$p<=0.05,]; rm(predixcan_liver)

predixcan_lymphocyte<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/edic_predixcan_lymphocyte_assoc_result.csv")
predixcan_lymphocyte_p5<-predixcan_lymphocyte[predixcan_lymphocyte$p<=0.05,]; rm(predixcan_lymphocyte)

predixcan_whole_blood<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/whole_blood/edic_predixcan_whole_blood_assoc_result.csv")
predixcan_whole_blood_p5<-predixcan_whole_blood[predixcan_whole_blood$p<=0.05,]; rm(predixcan_whole_blood)

gtex487_predixcan_cerebellum_p5 <-merge(predixcan_cerebellum_p5, gtex_487, by="hgnc_symbol"); rm(predixcan_cerebellum_p5) 
gtex487_predixcan_heart_atrial_appendage_p5 <-merge(predixcan_heart_atrial_appendage_p5, gtex_487, by="hgnc_symbol"); rm(predixcan_heart_atrial_appendage_p5) 
gtex487_predixcan_heart_left_ventricle_p5 <-merge(predixcan_heart_left_ventricle_p5, gtex_487, by="hgnc_symbol"); rm(predixcan_heart_left_ventricle_p5) 
gtex487_predixcan_liver_p5 <-merge(predixcan_liver_p5, gtex_487, by="hgnc_symbol"); rm(predixcan_liver_p5) 
gtex487_predixcan_lymphocyte_p5 <-merge(predixcan_lymphocyte_p5, gtex_487, by="hgnc_symbol"); rm(predixcan_lymphocyte_p5) 
gtex487_predixcan_whole_blood_p5 <-merge(predixcan_whole_blood_p5, gtex_487, by="hgnc_symbol"); rm(predixcan_whole_blood_p5) 

write.csv(gtex487_predixcan_cerebellum_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/cerebellum/intersect_gtex_v7_487_predixcan_cerebellum_p5.csv",row.names = FALSE)
write.csv(gtex487_predixcan_heart_atrial_appendage_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/heart_atrial_appendage/intersect_gtex_v7_487_predixcan_heart_atrial_appendage_p5.csv",row.names = FALSE)
write.csv(gtex487_predixcan_heart_left_ventricle_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/heart_left_ventricle/intersect_gtex_v7_487_predixcan_heart_left_ventricle_p5.csv",row.names = FALSE)
write.csv(gtex487_predixcan_liver_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/liver/intersect_gtex_v7_487_predixcan_liver_p5.csv",row.names = FALSE)
write.csv(gtex487_predixcan_lymphocyte_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/intersect_gtex_v7_487_predixcan_lymphocyte_p5.csv",row.names = FALSE)
write.csv(gtex487_predixcan_whole_blood_p5,"/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/whole_blood/intersect_gtex_v7_487_predixcan_whole_blood_p5.csv",row.names = FALSE)

#------------------------------
# predixcan fold change
#------------------------------
rm(list=ls())
library(data.table)
tissue="lymphocyte"  #cerebellum, heart_atrial_appendage, heart_left_ventricle, liver, lymphocyte, whole_blood
library(limma)

predixcan<-fread(paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/edic_predixcan_",tissue,"_predicted_expression.txt",sep=""))
predixcan<-data.frame(predixcan)
rownames(predixcan)<-predixcan[,1]
predixcan[,1:2]<-NULL
predixcan.t<-data.frame(t(predixcan))


# find gene symbol
library("biomaRt")

mart <- useDataset(dataset = "hsapiens_gene_ensembl",
                   mart = useMart("ensembl",
                                  host    = "www.ensembl.org"))

# attribs = c("ensembl_gene_id", "entrezgene", "hgnc_symbol",
            # "chromosome_name", "start_position","end_position")
attribs = c("ensembl_gene_id", "hgnc_symbol")


gene_set<-NULL
for (i in 1:length(rownames(predixcan.t))){
  ensl<-unlist(strsplit(as.character(rownames(predixcan.t)[i]),".",fixed=T))[1]
  gene_set[i]<-ensl
  rownames(predixcan.t)[i]<-ensl
}

genesymbol <- getBM(attributes = attribs, filters = "ensembl_gene_id",
                     values = gene_set, mart = mart)

genesymbol<-genesymbol[!genesymbol$hgnc_symbol=="",]
rownames(genesymbol)<-genesymbol$ensembl_gene_id; genesymbol$ensembl_gene_id<-NULL

pheno<-fread("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/phenotype/igrowth.txt")
pheno$FID<-NULL; pheno$pid<-NULL; pheno$mid<-NULL; pheno$sex<-NULL
group<-factor(pheno$igrowth)
extrs<-predixcan.t[,pheno$IID]


design<-model.matrix(~0+group)
colnames(design)<-c("case","control")
cont.matrix<-makeContrasts(levels=colnames(design), DE=case-control)
fit<-lmFit(extrs,design)
fit$gene<-genesymbol
fit2<-contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2)

x1=topTable(fit2, coef="DE", n=nrow(genesymbol ))

# x1=topTable(fit2, coef="DE",genelist=genesymbol, n=nrow(genesymbol))
x1<-na.omit(x1)
x2<-x1[rownames(genesymbol),]
x2<-cbind(genesymbol$hgnc_symbol,x2)
x2<-na.omit(x2); x2$t<-NULL;x2$B<-NULL
write.csv(x2,paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/",tissue,"/predixcan_",tissue,"_FC.csv",sep=""), row.names=F)

x3<-x2[order(x2$P.Value),]

library(ggplot2)
library(limma)
library(gtools); library(bioDist); library(calibrate)
library(plyr); library(reshape2); library(scales)
library(ggfortify)

pval.cutoff<-0.01
# FC.cutoff <-0.17
geneSymbol<-as.character(x3$`genesymbol$hgnc_symbol`)
setEPS()
postscript(file = paste("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/pval.cutoff_", pval.cutoff, ".eps",sep=""))
with(x3, plot(logFC, -log10(P.Value), pch=20, col="gray", main=paste("Group difference (PDR vs. nPDR) with P-value <0.01",sep=""), xlab="log2 fold change", ylab="-log10(P.Value)", xlim=c(-0.2,0.2)))
with(subset(x3, P.Value < pval.cutoff), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
with(subset(x3, P.Value < pval.cutoff), textxy(logFC, -log10(P.Value), labs=geneSymbol[1:10], cex=.8))

# with(subset(x3, P.Value < pval.cutoff & abs(logFC) > FC.cutoff), points(logFC, -log10(P.Value), pch=20, col="red"), cex=1.1)
# with(subset(x3, P.Value < pval.cutoff  & abs(logFC) > FC.cutoff), textxy(logFC, -log10(P.Value), labs=genesymbol, cex=.8))
dev.off()
# DCA

# check with Andrew 