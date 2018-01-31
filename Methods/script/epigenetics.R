final_result <- function(file1, file2){
  m1<-merge(file1,file2,by.x="geneSymbol",by.y="Symbol")
  write.csv(m1,paste("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/Epigenome/common_genes_between_",deparse(substitute(file1)),"_",deparse(substitute(file2)),".csv",sep=""), row.names=F)
}

Hyper_MLs_Monos<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/epigenome/Hyper-MLs_Monos.csv",header = T)
Hypo_MLs_Monos<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/epigenome/Hypo-MLs_Monos.csv",header = T)
Hyper_MLs_WB_DNA<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/epigenome/Hyper-MLs_WB_DNA.csv",header = T)
Hypo_MLs_WB_DNAs<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/epigenome/Hypo-MLs_WB_DNA.csv",header = T)

RG_PDR_nPDR<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/4.Glucose_Response_analyses/difference_between_group/delta_DwC_DwoC_pvalue_all_2017-09-13.csv")
RG_PDR_nPDR$ID<-NULL; RG_PDR_nPDR$t<-NULL; RG_PDR_nPDR$adj.P.Val<-NULL; RG_PDR_nPDR$B<-NULL
HG_PDR_nPDR<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/3.Between_Group_Analysis/HG/DwC_vs_DwoC/hg_DwC_DwoC_PC1_GrowthRate_all_2017-09-13.csv")
HG_PDR_nPDR$ID<-NULL; HG_PDR_nPDR$t<-NULL; HG_PDR_nPDR$adj.P.Val<-NULL; HG_PDR_nPDR$B<-NULL
SG_PDR_nPDR<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/3.Between_Group_Analysis/SG/DwC_vs_DwoC/sg_DwC_DwoC_PC1_GrowthRate_all_2017-09-13.csv")
SG_PDR_nPDR$ID<-NULL; SG_PDR_nPDR$t<-NULL; SG_PDR_nPDR$adj.P.Val<-NULL; SG_PDR_nPDR$B<-NULL

final_result(RG_PDR_nPDR,Hyper_MLs_Monos)
final_result(HG_PDR_nPDR,Hyper_MLs_Monos)
final_result(SG_PDR_nPDR,Hyper_MLs_Monos)

final_result(RG_PDR_nPDR,Hypo_MLs_Monos)
final_result(HG_PDR_nPDR,Hypo_MLs_Monos)
final_result(SG_PDR_nPDR,Hypo_MLs_Monos)

final_result(RG_PDR_nPDR,Hyper_MLs_WB_DNA)
final_result(HG_PDR_nPDR,Hyper_MLs_WB_DNA)
final_result(SG_PDR_nPDR,Hyper_MLs_WB_DNA)

final_result(RG_PDR_nPDR,Hypo_MLs_WB_DNAs)
final_result(HG_PDR_nPDR,Hypo_MLs_WB_DNAs)
final_result(SG_PDR_nPDR,Hypo_MLs_WB_DNAs)


