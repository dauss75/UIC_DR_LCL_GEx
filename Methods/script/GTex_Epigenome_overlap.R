GTex_Home<-"/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/GTex/results/"
Epi_Home<-"/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/epigenome/results/"

GTex_HG_PDR_nPDR<-read.csv(paste(GTex_Home,"486snp_gene_in_HG_PDR_nPDR.csv",sep=""))
GTex_SG_PDR_nPDR<-read.csv(paste(GTex_Home,"486snp_gene_in_SG_PDR_nPDR.csv",sep=""))
GTex_RG_PDR_nPDR<-read.csv(paste(GTex_Home,"486snp_gene_in_RG_PDR_nPDR.csv",sep=""))

Epi_HG_PDR_nPDR_Hyper_Monos<-read.csv(paste(Epi_Home,"common_genes_between_HG_PDR_nPDR_Hyper_MLs_Monos.csv",sep=""))
Epi_HG_PDR_nPDR_Hyper_WB<-read.csv(paste(Epi_Home,"common_genes_between_HG_PDR_nPDR_Hyper_MLs_WB_DNA.csv",sep=""))
Epi_HG_PDR_nPDR_Hypo_Monos<-read.csv(paste(Epi_Home,"common_genes_between_HG_PDR_nPDR_Hypo_MLs_Monos.csv",sep=""))
Epi_HG_PDR_nPDR_Hypo_WB<-read.csv(paste(Epi_Home,"common_genes_between_HG_PDR_nPDR_Hypo_MLs_WB_DNA.csv",sep=""))

Epi_SG_PDR_nPDR_Hyper_Monos<-read.csv(paste(Epi_Home,"common_genes_between_SG_PDR_nPDR_Hyper_MLs_Monos.csv",sep=""))
Epi_SG_PDR_nPDR_Hyper_WB<-read.csv(paste(Epi_Home,"common_genes_between_SG_PDR_nPDR_Hyper_MLs_WB_DNA.csv",sep=""))
Epi_SG_PDR_nPDR_Hypo_Monos<-read.csv(paste(Epi_Home,"common_genes_between_SG_PDR_nPDR_Hypo_MLs_Monos.csv",sep=""))
Epi_SG_PDR_nPDR_Hypo_WB<-read.csv(paste(Epi_Home,"common_genes_between_SG_PDR_nPDR_Hypo_MLs_WB_DNA.csv",sep=""))

Epi_RG_PDR_nPDR_Hyper_Monos<-read.csv(paste(Epi_Home,"common_genes_between_RG_PDR_nPDR_Hyper_MLs_Monos.csv",sep=""))
Epi_RG_PDR_nPDR_Hyper_WB<-read.csv(paste(Epi_Home,"common_genes_between_RG_PDR_nPDR_Hyper_MLs_WB_DNA.csv",sep=""))
Epi_RG_PDR_nPDR_Hypo_Monos<-read.csv(paste(Epi_Home,"common_genes_between_RG_PDR_nPDR_Hypo_MLs_Monos.csv",sep=""))
Epi_RG_PDR_nPDR_Hypo_WB<-read.csv(paste(Epi_Home,"common_genes_between_RG_PDR_nPDR_Hypo_MLs_WB_DNA.csv",sep=""))

GTex_Epi_HG_PDR_nPDR_Hyper_Monos<-merge(GTex_HG_PDR_nPDR,Epi_HG_PDR_nPDR_Hyper_Monos, by.x="geneSymbol", by.y="geneSymbol")
GTex_Epi_HG_PDR_nPDR_Hypo_Monos<-merge(GTex_HG_PDR_nPDR,Epi_HG_PDR_nPDR_Hypo_Monos, by.x="geneSymbol", by.y="geneSymbol")
GTex_Epi_HG_PDR_nPDR_Hyper_WB<-merge(GTex_HG_PDR_nPDR,Epi_HG_PDR_nPDR_Hyper_WB, by.x="geneSymbol", by.y="geneSymbol")
GTex_Epi_HG_PDR_nPDR_Hypo_WB<-merge(GTex_HG_PDR_nPDR,Epi_HG_PDR_nPDR_Hypo_WB, by.x="geneSymbol", by.y="geneSymbol")

GTex_Epi_SG_PDR_nPDR_Hyper_Monos<-merge(GTex_SG_PDR_nPDR,Epi_SG_PDR_nPDR_Hyper_Monos, by.x="geneSymbol", by.y="geneSymbol")
GTex_Epi_SG_PDR_nPDR_Hypo_Monos<-merge(GTex_SG_PDR_nPDR,Epi_SG_PDR_nPDR_Hypo_Monos, by.x="geneSymbol", by.y="geneSymbol")
GTex_Epi_SG_PDR_nPDR_Hyper_WB<-merge(GTex_SG_PDR_nPDR,Epi_SG_PDR_nPDR_Hyper_WB, by.x="geneSymbol", by.y="geneSymbol")
GTex_Epi_SG_PDR_nPDR_Hypo_WB<-merge(GTex_SG_PDR_nPDR,Epi_SG_PDR_nPDR_Hypo_WB, by.x="geneSymbol", by.y="geneSymbol")

## HDAC4 histone deacetylase 4 0.20210299 4.618714 0.04190599 ENSG00000068024 Esophagus_Gastroesophageal_Junction_Analysis rs3828222

GTex_Epi_RG_PDR_nPDR_Hyper_Monos<-merge(GTex_RG_PDR_nPDR,Epi_RG_PDR_nPDR_Hyper_Monos, by.x="geneSymbol", by.y="geneSymbol")
GTex_Epi_RG_PDR_nPDR_Hypo_Monos<-merge(GTex_RG_PDR_nPDR,Epi_RG_PDR_nPDR_Hypo_Monos, by.x="geneSymbol", by.y="geneSymbol")
GTex_Epi_RG_PDR_nPDR_Hyper_WB<-merge(GTex_RG_PDR_nPDR,Epi_RG_PDR_nPDR_Hyper_WB, by.x="geneSymbol", by.y="geneSymbol")
GTex_Epi_RG_PDR_nPDR_Hypo_WB<-merge(GTex_RG_PDR_nPDR,Epi_RG_PDR_nPDR_Hypo_WB, by.x="geneSymbol", by.y="geneSymbol")

