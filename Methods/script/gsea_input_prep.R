#------------------------------------------------------------------------
# 1. RG_All
# setwd("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/4.Glucose_Response_analyses/RG")
RGAll<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/4.Glucose_Response_analyses/RG/All_treatment_effect_2017-06-13.csv")
# prank<-sign(t1$logFC)*(1/t1$P.Value)
prank<-sign(RGAll$logFC)*(-log10(RGAll$P.Value))
geneSymbol<-as.character(RGAll$geneSymbol)
RGAll2<-cbind(geneSymbol,prank)
RGAll3<-RGAll2[!duplicated(RGAll2[,1]),]
RGAll4<-RGAll3[order(as.numeric(RGAll3[,2]),decreasing = TRUE),]
write.table(RGAll4,file="/Users/sjung/Desktop/test/RGAll.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 1a. RG_All_FDR
# setwd("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/4.Glucose_Response_analyses/RG")
RGAll<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/4.Glucose_Response_analyses/RG/All_treatment_effect_2017-06-13.csv")
# prank<-sign(t1$logFC)*(1/t1$P.Value)
prank<-sign(RGAll$logFC)*(-log10(RGAll$adj.P.Val))
geneSymbol<-as.character(RGAll$geneSymbol)
RGAll2<-cbind(geneSymbol,prank)
RGAll3<-RGAll2[!duplicated(RGAll2[,1]),]
RGAll4<-RGAll3[order(as.numeric(RGAll3[,2]),decreasing = TRUE),]
write.table(RGAll4,file="/Users/sjung/Desktop/test/RGAll_FDR.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)


# 1. RG_nDM
# setwd("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/4.Glucose_Response_analyses/RG")
RGnDM<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/4.Glucose_Response_analyses/RG/NoD_treatment_effect_2017-06-13.csv")
# prank<-sign(t1$logFC)*(1/t1$P.Value)
prank<-sign(RGnDM$logFC)*(-log10(RGnDM$P.Value))
geneSymbol<-as.character(RGnDM$geneSymbol)
RGnDM2<-cbind(geneSymbol,prank)
RGnDM3<-RGnDM2[!duplicated(RGnDM2[,1]),]
RGnDM4<-RGnDM3[order(as.numeric(RGnDM3[,2]),decreasing = TRUE),]
write.table(RGnDM4,file="/Users/sjung/Desktop/test/RG_nDM.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 2. RG_DM
RGPDR<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/4.Glucose_Response_analyses/RG/DwC_treatment_effect_2017-06-13.csv")
RGnPDR<-read.csv("DwoC_treatment_effect_2017-06-13.csv")
RGDM<-rbind(RGPDR, RGnPDR)
prank<-sign(RGDM$logFC)*(-log10(RGDM$P.Value))
geneSymbol<-as.character(RGDM$geneSymbol)
RGDM2<-cbind(geneSymbol,prank)
RGDM3<-RGDM2[!duplicated(RGDM2[,1]),]
RGDM4<-RGDM3[order(as.numeric(RGDM3[,2]),decreasing = TRUE),]
write.table(RGDM4,file="/Users/sjung/Desktop/test/RG_DM.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 3. RG_PDR
RGPDR<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/4.Glucose_Response_analyses/RG/DwC_treatment_effect_2017-06-13.csv")
prank<-sign(RGPDR$logFC)*(-log10(RGPDR$P.Value))
geneSymbol<-as.character(RGPDR$geneSymbol)
RGPDR2<-cbind(geneSymbol,prank)
RGPDR3<-RGPDR2[!duplicated(RGPDR2[,1]),]
RGPDR4<-RGPDR3[order(as.numeric(RGPDR3[,2]),decreasing = TRUE),]
write.table(RGPDR4,file="/Users/sjung/Desktop/test/RG_PDR.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 3. RG_nPDR
RGnPDR<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/4.Glucose_Response_analyses/RG/DwoC_treatment_effect_2017-06-13.csv")
prank<-sign(RGnPDR$logFC)*(-log10(RGnPDR$P.Value))
geneSymbol<-as.character(RGnPDR$geneSymbol)
RGnPDR2<-cbind(geneSymbol,prank)
RGnPDR3<-RGnPDR2[!duplicated(RGnPDR2[,1]),]
RGnPDR4<-RGnPDR3[order(as.numeric(RGnPDR3[,2]),decreasing = TRUE),]
write.table(RGnPDR4,file="/Users/sjung/Desktop/test/RG_nPDR.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 4. RG_DM_nDM
RG_DM_nDM<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/4.Glucose_Response_analyses/difference_between_group/delta_NoD_Diabetes_pvalue_all_2017-06-13.csv")
prank<-sign(RG_DM_nDM$logFC)*(-log10(RG_DM_nDM$P.Value))
geneSymbol<-as.character(RG_DM_nDM$geneSymbol)
RG_DM_nDM2<-cbind(geneSymbol,prank)
RG_DM_nDM3<-RG_DM_nDM2[!duplicated(RG_DM_nDM2[,1]),]
RG_DM_nDM4<-RG_DM_nDM3[order(as.numeric(RG_DM_nDM3[,2]),decreasing = TRUE),]
write.table(RG_DM_nDM4,file="/Users/sjung/Desktop/test/RG_DM_nDM.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 5. RG_PDR_nPDR
RG_PDR_nPDR<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/4.Glucose_Response_analyses/difference_between_group/delta_DwC_DwoC_pvalue_all_2017-09-13.csv")
prank<-sign(RG_PDR_nPDR$logFC)*(-log10(RG_PDR_nPDR$P.Value))
geneSymbol<-as.character(RG_PDR_nPDR$geneSymbol)
RG_PDR_nPDR2<-cbind(geneSymbol,prank)
RG_PDR_nPDR3<-RG_PDR_nPDR2[!duplicated(RG_PDR_nPDR2[,1]),]
RG_PDR_nPDR4<-RG_PDR_nPDR3[order(as.numeric(RG_PDR_nPDR3[,2]),decreasing = TRUE),]
write.table(RG_PDR_nPDR4,file="/Users/sjung/Desktop/test/RG_PDR_nPDR.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 6. BG_SG_nDM_DM
BG_SG_nDM_DM<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/3.Between_Group_Analysis/SG/NoD_vs_Diabetes/sg_NoD_Diabetes_PC1_GrowthRate_all_2017-06-12.csv")
prank<-sign(BG_SG_nDM_DM$logFC)*(-log10(BG_SG_nDM_DM$P.Value))
geneSymbol<-as.character(BG_SG_nDM_DM$geneSymbol)
BG_SG_nDM_DM2<-cbind(geneSymbol,prank)
BG_SG_nDM_DM3<-BG_SG_nDM_DM2[!duplicated(BG_SG_nDM_DM2[,1]),]
BG_SG_nDM_DM4<-BG_SG_nDM_DM3[order(as.numeric(BG_SG_nDM_DM3[,2]),decreasing = TRUE),]
write.table(BG_SG_nDM_DM4,file="/Users/sjung/Desktop/test/BG_SG_nDM_DM.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 7. BG_SG_PDR_nPDR
BG_SG_PDR_nPDR<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/3.Between_Group_Analysis/SG/DwC_vs_DwoC/sg_DwoC_DwC_PC1_GrowthRate_all_2017-06-12.csv")
prank<-sign(BG_SG_PDR_nPDR$logFC)*(-log10(BG_SG_PDR_nPDR$P.Value))
geneSymbol<-as.character(BG_SG_PDR_nPDR$geneSymbol)
BG_SG_PDR_nPDR2<-cbind(geneSymbol,prank)
BG_SG_PDR_nPDR3<-BG_SG_PDR_nPDR2[!duplicated(BG_SG_PDR_nPDR2[,1]),]
BG_SG_PDR_nPDR4<-BG_SG_PDR_nPDR3[order(as.numeric(BG_SG_PDR_nPDR3[,2]),decreasing = TRUE),]
write.table(BG_SG_PDR_nPDR4,file="/Users/sjung/Desktop/test/BG_SG_PDR_nPDR.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 8. BG_HG_nDM_DM
BG_HG_nDM_DM<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/3.Between_Group_Analysis/HG/NoD_vs_Diabetes/hg_NoD_Diabetes_PC1_GrowthRate_all_2017-06-12.csv")
prank<-sign(BG_HG_nDM_DM$logFC)*(-log10(BG_HG_nDM_DM$P.Value))
geneSymbol<-as.character(BG_HG_nDM_DM$geneSymbol)
BG_HG_nDM_DM2<-cbind(geneSymbol,prank)
BG_HG_nDM_DM3<-BG_HG_nDM_DM2[!duplicated(BG_HG_nDM_DM2[,1]),]
BG_HG_nDM_DM4<-BG_HG_nDM_DM3[order(as.numeric(BG_HG_nDM_DM3[,2]),decreasing = TRUE),]
write.table(BG_HG_nDM_DM4,file="/Users/sjung/Desktop/test/BG_HG_nDM_DM.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 9. BG_HG_PDR_nPDR
BG_HG_PDR_nPDR<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/3.Between_Group_Analysis/HG/DwC_vs_DwoC/hg_DwoC_DwC_PC1_GrowthRate_all_2017-06-12.csv")
prank<-sign(BG_HG_PDR_nPDR$logFC)*(-log10(BG_HG_PDR_nPDR$P.Value))
geneSymbol<-as.character(BG_HG_PDR_nPDR$geneSymbol)
BG_HG_PDR_nPDR2<-cbind(geneSymbol,prank)
BG_HG_PDR_nPDR3<-BG_HG_PDR_nPDR2[!duplicated(BG_HG_PDR_nPDR2[,1]),]
BG_HG_PDR_nPDR4<-BG_HG_PDR_nPDR3[order(as.numeric(BG_HG_PDR_nPDR3[,2]),decreasing = TRUE),]
write.table(BG_HG_PDR_nPDR4,file="/Users/sjung/Desktop/test/BG_HG_PDR_nPDR.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 10. EDIC_lymphocyte
EDIC_lymphocyte<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/5.indenepdnet_data/EDIC_lymphocyte_affy_DE.csv")
prank<-sign(EDIC_lymphocyte$logFC)*(-log10(EDIC_lymphocyte$P.Value))
geneSymbol<-as.character(EDIC_lymphocyte$Symbol)
EDIC_lymphocyte2<-cbind(geneSymbol,prank)
EDIC_lymphocyte3<-EDIC_lymphocyte2[!duplicated(EDIC_lymphocyte2[,1]),]
EDIC_lymphocyte4<-EDIC_lymphocyte3[order(as.numeric(EDIC_lymphocyte3[,2]),decreasing = TRUE),]
write.table(EDIC_lymphocyte4,file="/Users/sjung/Desktop/test/EDIC_lymphocyte.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 11. PBMC
pbmc<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/5.indenepdnet_data/EDIC_lymphocyte_affy_DE.csv")
prank<-sign(pbmc$logFC)*(-log10(pbmc$P.Value))
geneSymbol<-as.character(pbmc$Symbol)
pbmc2<-cbind(geneSymbol,prank)
pbmc3<-pbmc2[!duplicated(pbmc2[,1]),]
pbmc4<-pbmc3[order(as.numeric(pbmc3[,2]),decreasing = TRUE),]
write.table(pbmc4,file="/Users/sjung/Desktop/test/pbmc.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 10. EDIC_monocyte
EDIC_monocyte<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/MG_BS/output/5.indenepdnet_data/EDIC_monocyte_affy_DE.csv")
prank<-sign(EDIC_monocyte$logFC)*(-log10(EDIC_monocyte$P.Value))
geneSymbol<-as.character(EDIC_monocyte$Symbol)
EDIC_monocyte2<-cbind(geneSymbol,prank)
EDIC_monocyte3<-EDIC_monocyte2[!duplicated(EDIC_monocyte2[,1]),]
EDIC_monocyte4<-EDIC_monocyte3[order(as.numeric(EDIC_monocyte3[,2]),decreasing = TRUE),]
write.table(EDIC_monocyte4,file="/Users/sjung/Desktop/test/EDIC_monocyte.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 11. Predixcan_lymphoblast
Predixcan_lymphoblast<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/lymphocyte/edic_predixcan_lymphocyte_assoc_result.csv")
prank <- -log10(Predixcan_lymphoblast$p)
geneSymbol<-as.character(Predixcan_lymphoblast$hgnc_symbol)
Predixcan_lymphoblast2<-cbind(geneSymbol,prank)
Predixcan_lymphoblast3<-Predixcan_lymphoblast2[!duplicated(Predixcan_lymphoblast2[,1]),]
Predixcan_lymphoblast4<-Predixcan_lymphoblast3[order(as.numeric(Predixcan_lymphoblast3[,2]),decreasing = TRUE),]
write.table(Predixcan_lymphoblast4,file="/Users/sjung/Desktop/test/Predixcan_lymphoblast.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 12. Predixcan_cerebellum
Predixcan_cerebellum<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/cerebellum/edic_predixcan_cerebellum_assoc_result.csv")
prank <- -log10(Predixcan_cerebellum$p)
geneSymbol<-as.character(Predixcan_cerebellum$hgnc_symbol)
Predixcan_cerebellum2<-cbind(geneSymbol,prank)
Predixcan_cerebellum3<-Predixcan_cerebellum2[!duplicated(Predixcan_cerebellum2[,1]),]
Predixcan_cerebellum4<-Predixcan_cerebellum3[order(as.numeric(Predixcan_cerebellum3[,2]),decreasing = TRUE),]
write.table(Predixcan_cerebellum4,file="/Users/sjung/Desktop/test/Predixcan_cerebellum.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 13. Predixcan_heart_left_ventricle
Predixcan_heart_left_ventricle<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/heart_left_ventricle/edic_predixcan_heart_left_ventricle_assoc_result.csv")
prank <- -log10(Predixcan_heart_left_ventricle$p)
geneSymbol<-as.character(Predixcan_heart_left_ventricle$hgnc_symbol)
Predixcan_heart_left_ventricle2<-cbind(geneSymbol,prank)
Predixcan_heart_left_ventricle3<-Predixcan_heart_left_ventricle2[!duplicated(Predixcan_heart_left_ventricle2[,1]),]
Predixcan_heart_left_ventricle4<-Predixcan_heart_left_ventricle3[order(as.numeric(Predixcan_heart_left_ventricle3[,2]),decreasing = TRUE),]
write.table(Predixcan_heart_left_ventricle4,file="/Users/sjung/Desktop/test/Predixcan_heart_left_ventricle.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 14. Predixcan_heart_atrial_appendage
Predixcan_heart_atrial_appendage<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/heart_atrial_appendage/edic_predixcan_heart_atrial_appendage_assoc_result.csv")
prank <- -log10(Predixcan_heart_atrial_appendage$p)
geneSymbol<-as.character(Predixcan_heart_atrial_appendage$hgnc_symbol)
Predixcan_heart_atrial_appendage2<-cbind(geneSymbol,prank)
Predixcan_heart_atrial_appendage3<-Predixcan_heart_atrial_appendage2[!duplicated(Predixcan_heart_atrial_appendage2[,1]),]
Predixcan_heart_atrial_appendage4<-Predixcan_heart_atrial_appendage3[order(as.numeric(Predixcan_heart_atrial_appendage3[,2]),decreasing = TRUE),]
write.table(Predixcan_heart_atrial_appendage4,file="/Users/sjung/Desktop/test/Predixcan_heart_atrial_appendage.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 15. Predixcan_whole_blood
Predixcan_whole_blood<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/whole_blood/edic_predixcan_whole_blood_assoc_result.csv")
prank <- -log10(Predixcan_whole_blood$p)
geneSymbol<-as.character(Predixcan_whole_blood$hgnc_symbol)
Predixcan_whole_blood2<-cbind(geneSymbol,prank)
Predixcan_whole_blood3<-Predixcan_whole_blood2[!duplicated(Predixcan_whole_blood2[,1]),]
Predixcan_whole_blood4<-Predixcan_whole_blood3[order(as.numeric(Predixcan_whole_blood3[,2]),decreasing = TRUE),]
write.table(Predixcan_whole_blood4,file="/Users/sjung/Desktop/test/Predixcan_whole_blood.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)

# 16. Predixcan_liver
Predixcan_liver<-read.csv("/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/results/liver/edic_predixcan_liver_assoc_result.csv")
prank <- -log10(Predixcan_liver$p)
geneSymbol<-as.character(Predixcan_liver$hgnc_symbol)
Predixcan_liver2<-cbind(geneSymbol,prank)
Predixcan_liver3<-Predixcan_liver2[!duplicated(Predixcan_liver2[,1]),]
Predixcan_liver4<-Predixcan_liver3[order(as.numeric(Predixcan_liver3[,2]),decreasing = TRUE),]
write.table(Predixcan_liver4,file="/Users/sjung/Desktop/test/Predixcan_liver.rnk", sep="\t", col.names = FALSE, row.names = FALSE,  quote = FALSE)
