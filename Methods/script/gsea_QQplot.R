## gsea up and down-regulated pathway QQplot

library(ggplot2)
library(gtools); library(bioDist); library(calibrate)
library(plyr); library(reshape2); library(scales)
library(ggfortify)
library(qqman)

up<-read.csv("/Users/sjung/gsea_home/output/c5/RG_All.GseaPreranked.1500567875906/gsea_report_for_na_pos_1500567875906.csv")
down<-read.csv("/Users/sjung/gsea_home/output/c5/RG_All.GseaPreranked.1500567875906/gsea_report_for_na_neg_1500567875906.csv")

fdr<-0.01 # 1%
up_fdr1<-up[up$FDR.q.val<=fdr,]; 
up_fdr1$group<-2
down_fdr1<-down[down$FDR.q.val<=fdr,]; down_fdr1$group<-4
fdr1<-rbind(up_fdr1,down_fdr1)
fdr1_sorted<-fdr1[order(fdr1$FDR.q.val),]
fdr1_sorted$FDR.q.val2<-fdr1_sorted$FDR.q.val+1e-10
rm(down,fdr1)

setEPS()
postscript(file = "/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/qq_plot.eps")

# qq(fdr1_sorted$FDR.q.val2, main = "RG_All", xlim = c(0, 3), ylim = c(0, 4), pch = 17, col=fdr1_sorted$group, cex = 1, las = 1)
qq(fdr1_sorted$FDR.q.val2, main = "RG_All", pch = 17, col=fdr1_sorted$group, cex = 1, las = 1)

dev.off()

fdr2<-rbind(up,down)
fdr2_sorted<-fdr2[order(fdr2$FDR.q.val),]
fdr2_sorted$FDR.q.val2<-fdr2_sorted$FDR.q.val+1e-10
rm(down,fdr1)

setEPS()
postscript(file = "/Users/sjung/Project/GlobusGenomics/UIC/predixcan/edic_run/qq_plot2.eps")

# qq(fdr1_sorted$FDR.q.val2, main = "RG_All", xlim = c(0, 3), ylim = c(0, 4), pch = 17, col=fdr1_sorted$group, cex = 1, las = 1)
qq(fdr2_sorted$FDR.q.val2, main = "RG_All", pch = 17, col=fdr1_sorted$group, cex = 1, las = 1)

dev.off()
