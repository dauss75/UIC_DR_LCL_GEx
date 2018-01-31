snp<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/gwas/MG_gwas_data.csv")
# gtex_files = list.files(path="/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/GTEx_Analysis_v6p_eQTL", pattern='txt', full.names=TRUE)

snp<-data.frame(snp[!duplicated(snp),])
snp_gene<-data.frame(matrix(data=NA, nrow=1000, ncol = 2))
x<-read.table("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/GTEx_Analysis_v6p_eQTL/Cells_EBV-transformed_lymphocytes_Analysis.v6p.egenes.txt",header=T,row.names=1)


j=1
for (i in 1:nrow(snp)){
  print(i)
  # for (j in 1:length(gtex_files)) {
      # x <-read.table(gtex_files[j],header=T,row.names=1)

      if (nrow(x[grep(snp[i,1], x$rs_id_dbSNP142_GRCh37p13, value = FALSE),])!=0){
        print(j)
        print(x[grep(snp[i,1], x$rs_id_dbSNP142_GRCh37p13, value = FALSE),c(1,16)])
        snp_gene[j,]<-unlist(x[grep(snp[i,1], x$rs_id_dbSNP142_GRCh37p13, value = FALSE),c(1,16)])
        j=j+1
      }
  # }
}

snp_gene<-na.omit(snp_gene)
write.csv(snp_gene,"/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/Cells_EBV-transformed_lymphocytes_snp_gene.csv", row.names=F)
      
snp_gene<-data.frame(matrix(data=NA, nrow=1000, ncol = 2))
y<-read.table("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/Brain_Cerebellum_Analysis.v6p.egenes.txt",header=T,row.names=1)

j=1
for (i in 1:nrow(snp)){
  print(i)
  # for (j in 1:length(gtex_files)) {
  # x <-read.table(gtex_files[j],header=T,row.names=1)
  
  if (nrow(y[grep(snp[i,1], y$rs_id_dbSNP142_GRCh37p13, value = FALSE),])!=0){
    print(j)
    print(y[grep(snp[i,1], y$rs_id_dbSNP142_GRCh37p13, value = FALSE),c(1,16)])
    snp_gene[j,]<-unlist(y[grep(snp[i,1], y$rs_id_dbSNP142_GRCh37p13, value = FALSE),c(1,16)])
    j=j+1
  }
  # }
}

snp_gene<-na.omit(snp_gene)
write.csv(snp_gene,"/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/Brain_Cerebellum_Analysis_snp_gene.csv", row.names=F)

rg_nod_vs_d<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/4.Glucose_Response_analyses/difference_between_group/delta_NoD_Diabetes_pvalue_all_2017-06-13.csv")
hg_nod_vs_d<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/3.Between_Group_Analysis/HG/NoD_vs_Diabetes/hg_NoD_Diabetes_PC1_GrowthRate_all_2017-06-12.csv")
sg_nod_vs_d<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/3.Between_Group_Analysis/SG/NoD_vs_Diabetes/sg_NoD_Diabetes_PC1_GrowthRate_all_2017-06-12.csv")


