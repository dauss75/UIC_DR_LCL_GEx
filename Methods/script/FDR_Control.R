library(qvalue)
HOME<-"/Users/sjung/Project/GlobusGenomics/UIC/"
setwd(HOME)


# pre-processing data for equal variance and treat effect
#--------------------------------------------------------
data<-read.csv(file="DeltaAnalysisResults10.5.16.csv")
colnames(data)[which(colnames(data)=='pvalue..test.equal.variance.')]<-"pval.equal.var"
colnames(data)[which(colnames(data)=='p_t..test.treatment.effect.')]<-"pval.treat.eff"

p1<-data$pval.equal.var
p2<-data$pval.treat.eff

qval.equal.var<-qvalue(p1)
qval.treat.eff<-qvalue(p2)

data$qval.equal.var<-qval.equal.var$qvalues
data$qval.treat.eff<-qval.treat.eff$qvalues

data.equal.var.p0.01<-data[which(data$pval.equal.var<0.01),]
data.treat.eff.q0.05<-data[which(data$qval.treat.eff<0.05),]

# plot for equal variance and treat effect
#--------------------------------------------------------
library(ggplot2)
# 1. pval.equal.var distribution
ggplot(data,aes(x=pval.equal.var))+ geom_histogram(aes(y=..density..), binwidth=.025, colour="black", fill="white")+     geom_density(alpha=.2, color="blue", linetype="dashed", size=0.5)+geom_vline(aes(xintercept=0.01), color="red", linetype="dashed", size=1) +theme_bw()
ggplot(data,aes(x=qval.equal.var))+ geom_histogram(aes(y=..density..), binwidth=.025, colour="black", fill="white")+     geom_density(alpha=.2, color="blue", linetype="dashed", size=0.5) +theme_bw()


# 2. pval.treat.eff distribution
ggplot(data,aes(x=pval.treat.eff))+ geom_histogram(aes(y=..density..), binwidth=.025, colour="black", fill="white")+     geom_density(alpha=.2, color="blue", linetype="dashed", size=0.5)+theme_bw()
ggplot(data,aes(x=qval.treat.eff))+ geom_histogram(aes(y=..density..), binwidth=.025, colour="black", fill="white")+     geom_density(alpha=.2, color="blue", linetype="dashed", size=0.5)+geom_vline(aes(xintercept=0.1), color="red", linetype="dashed", size=1) +theme_bw()


geneExp<-read.csv(file=paste(HOME,"310/output/mean_delta.csv",sep=""))


# data$qval.treat.eff<-qval.treat.eff$qvalues
data.q0.05<-data[which(data$qval.treat.eff<0.05),]
write.csv(data.q0.05, paste(HOME,"310/output/q5%.csv",sep=""))



data.p0.01q0.05<-intersect(data.q0.1$symbol,data.p0.01$symbol)

rownames(geneExp)<-geneExp$X
rownames(data.q0.05)<-data.q0.05$geneID
tt<-merge(data.q0.05,geneExp,by="row.names", all=T)

