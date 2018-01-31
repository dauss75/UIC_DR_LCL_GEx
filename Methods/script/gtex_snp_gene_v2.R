library("biomaRt")

mart <- useDataset(dataset = "hsapiens_gene_ensembl",
                   mart = useMart("ensembl",
                                  host    = "www.ensembl.org"))

attribs = c("ensembl_gene_id", "entrezgene", "hgnc_symbol",
            "chromosome_name", "start_position","end_position")

# final_result <- function(no_snp, expr_data){
#   gtex<-read.csv(paste("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/",no_snp,"snp_vs_GTex_nonredundant.csv",sep=""))
#   gtex<-gtex[!duplicated(gtex$gene),]
#   tmp<-expr_data[which(expr_data$geneSymbol %in% gtex$gene),]
#   tmp<-tmp[!duplicated(tmp$geneSymbol),]
#   write.csv(tmp,paste("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/",no_snp,"snp_gene_in_",deparse(substitute(expr_data)),".csv",sep=""), row.names=F)
# }

final_result <- function(no_snp, expr_data){
  # gtex<-read.csv(paste("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/",no_snp,"snp_vs_GTex_nonredundant.csv",sep=""))
  gtex<-read.csv(paste("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/",no_snp,"snp_vs_GTex_v7_nonredundant.csv",sep=""))
  # gtex<-gtex[!duplicated(gtex$gene),]
  m1<-merge(expr_data,gtex,by.x="geneSymbol",by.y="hgnc_symbol")
  m1<-m1[order(m1$P.Value),]
  m1<-m1[!duplicated(m1$geneSymbol),]
  write.csv(m1,paste("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/",no_snp,"snp_gene_in_",deparse(substitute(expr_data)),"gtex_v7.csv",sep=""), row.names=F)
}

convert_snpid<- function(ref, snp){
  m1<-merge(ref,snp,by.x="rs_id_dbSNP147_GRCh37p13",by.y="dbsnp")
  return(m1)
}

i=2  ## 1 for Mike data and 2 for GWAS data

library(data.table)
ref_snp<-fread("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/snp_lookup_table.csv", header=T)
ref_snp<-data.frame(ref_snp)



if (i==1){
  ## Mike data
  snp<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/gwas/MG_gwas_data.csv", header=F)
  colnames(snp)<-"dbsnp"
  snp<-snp[!duplicated(snp),,drop=FALSE]
  new_snpid<-convert_snpid(ref_snp,snp)
}
if (i==2){
  ## GWAS data
  snp<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/gwas/gwas_t1d_snp.csv")
  colnames(snp)<-"dbsnp"
  snp<-snp[!duplicated(snp),,drop=FALSE]
  new_snpid<-convert_snpid(ref_snp,snp)
}

# gtex_files = list.files(path="/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/GTEx_Analysis_v6p_eQTL", pattern='txt', full.names=TRUE)
gtex_files = list.files(path="/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/GTEx_Analysis_v7_eQTL", pattern='txt', full.names=TRUE)

snp<-new_snpid


library(foreach);library(doParallel)
no_cores <- detectCores() - 1
cl<-makeCluster(no_cores)
registerDoParallel(cl)
all<-NULL; ls<-NULL
for (j in 1:length(gtex_files)){
  x <-fread(gtex_files[j],header=T)
  ls<-foreach(i=1:nrow(snp), .combine=rbind) %dopar% {
    
    if (nrow(x[grep(snp[i,2], x$variant_id, value = FALSE),])!=0){
      tt<-x[grep(snp[i,2], x$variant_id, value = FALSE),]
      tissue<-unlist(strsplit(basename(gtex_files[j]),".",fixed=TRUE))[1]
      cbind(as.character(tissue),as.character(snp[i,1]),as.character(tt$variant_id),as.character(tt$gene_id), as.character(tt$pval_nominal))
      # cbind(as.character(tissue),as.character(tt[1]),as.character(tt[2]),as.character(tt[3]))
    }
  }
  all<-rbind(all,ls)
  print(gtex_files[j])
  print(nrow(ls))
}

stopCluster(cl)

no_snp<-nrow(snp)

colnames(all)<-c("tissue","dbsnp","variant_id","gene", "pval_nominal")
write.csv(all,paste("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/",no_snp,"snp_vs_GTex_v7_Full.csv",sep=""), row.names=F)

nonredundant<-data.frame(all[!duplicated(all[,c(2,4)]),])
write.csv(nonredundant,paste("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/",no_snp,"snp_vs_GTex_v7_nonredundant.csv",sep=""), row.names=F)

gene_set<-NULL
for (i in 1:length(nonredundant$gene)){
  gene_set[i]<-unlist(strsplit(as.character(nonredundant$gene[i]),".",fixed=T))[1]
}
nonredundant$gene<-gene_set
resultTable <- getBM(attributes = attribs, filters = "ensembl_gene_id",
                     values = gene_set, mart = mart)
resultTable2<-resultTable[!resultTable$hgnc_symbol=="",]
tt<-merge(nonredundant,resultTable2,by.x="gene",by.y="ensembl_gene_id")
tt<-tt[!duplicated(tt[,c(1,2,3)]),]


write.csv(tt,paste("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/",no_snp,"snp_vs_GTex_v7_nonredundant.csv",sep=""), row.names=F)
# convert dbsnp to GTex variant id
# snp_lookup_table<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/snp_lookup_table.csv")


## intersect with our expression data
RG_NoD_vs_D<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/4.Glucose_Response_analyses/difference_between_group/delta_NoD_Diabetes_pvalue_all_2017-06-13.csv")
RG_NoD_vs_D$ID<-NULL; RG_NoD_vs_D$t<-NULL; RG_NoD_vs_D$B<-NULL
HG_NoD_vs_D<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/3.Between_Group_Analysis/HG/NoD_vs_Diabetes/hg_NoD_Diabetes_PC1_GrowthRate_all_2017-06-12.csv")
HG_NoD_vs_D$ID<-NULL; HG_NoD_vs_D$t<-NULL; HG_NoD_vs_D$B<-NULL
SG_NoD_vs_D<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/3.Between_Group_Analysis/SG/NoD_vs_Diabetes/sg_NoD_Diabetes_PC1_GrowthRate_all_2017-06-12.csv")
SG_NoD_vs_D$ID<-NULL; SG_NoD_vs_D$t<-NULL; SG_NoD_vs_D$B<-NULL

# no_snp=78

final_result(no_snp,RG_NoD_vs_D)
final_result(no_snp,HG_NoD_vs_D)
final_result(no_snp,SG_NoD_vs_D)



RG_PDR_nPDR<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/4.Glucose_Response_analyses/difference_between_group/delta_DwC_DwoC_pvalue_all_2017-09-13.csv")
RG_PDR_nPDR$ID<-NULL; RG_PDR_nPDR$t<-NULL; RG_PDR_nPDR$adj.P.Val<-NULL; RG_PDR_nPDR$B<-NULL
HG_PDR_nPDR<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/3.Between_Group_Analysis/HG/DwC_vs_DwoC/hg_DwC_DwoC_PC1_GrowthRate_all_2017-09-13.csv")
HG_PDR_nPDR$ID<-NULL; HG_PDR_nPDR$t<-NULL; HG_PDR_nPDR$adj.P.Val<-NULL; HG_PDR_nPDR$B<-NULL
SG_PDR_nPDR<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/3.Between_Group_Analysis/SG/DwC_vs_DwoC/sg_DwC_DwoC_PC1_GrowthRate_all_2017-09-13.csv")
SG_PDR_nPDR$ID<-NULL; SG_PDR_nPDR$t<-NULL; SG_PDR_nPDR$adj.P.Val<-NULL; SG_PDR_nPDR$B<-NULL

# no_snp<-487
final_result(no_snp,RG_PDR_nPDR)
final_result(no_snp,HG_PDR_nPDR)
final_result(no_snp,SG_PDR_nPDR)
