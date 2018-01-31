snp_gtex<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/data/edic_pval/486snp_vs_GTex_nonredundant.csv")
edic<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/data/edic_pval/edic_pval.csv")

for (i in 1:nrow(snp_gtex)) {
  print(i)
  print(match(snp_gtex$dbsnp[i],edic$dbsnp))
  snp_gtex$EDIC.p.val[i]<-edic[match(snp_gtex$dbsnp[i],edic$dbsnp),2]
}

write.csv(snp_gtex,"/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/data/edic_pval/486snp_vs_GTex_nonredundant2.csv",row.names = FALSE)
