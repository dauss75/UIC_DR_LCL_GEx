HOME<-"/Users/sjung/Project/GlobusGenomics/UIC/"
setwd(HOME)
delta<-read.csv(file=paste(HOME,"310/output/raw_delta.csv", sep=""))
rownames(delta)<-delta$X
delta$X<-NULL

delta.T<-data.frame(t(delta))
delta.T$Name<-substr(rownames(delta.T),1,nchar(rownames(delta.T))-2)
delta.T$Name2<-substr(rownames(delta.T),1,4)
delta.pca<-prcomp(delta.T[,c(-ncol(delta.T)+1,-ncol(delta.T))], scale=TRUE, center=TRUE)
plot(delta.pca, type="l", main="delta PCA")

# plots for PCs
library(ggbiplot)
ggbiplot(delta.pca, choices = 1:2 , var.axes =FALSE, scale=TRUE, var.scale=1, circle=TRUE)
ggbiplot(delta.pca, choices = 2:3 , var.axes =FALSE, scale=TRUE, var.scale=1, circle=TRUE)
ggbiplot(delta.pca, choices = 3:4 , var.axes =FALSE, scale=TRUE, var.scale=1, circle=TRUE)
ggbiplot(delta.pca, choices = 4:5 , var.axes =FALSE, scale=TRUE, var.scale=1, circle=TRUE)
ggbiplot(delta.pca, choices = 5:6 , var.axes =FALSE, scale=TRUE, var.scale=1, circle=TRUE)



library(ggfortify)
autoplot(delta.pca, data=delta.T, colour="Name", label = FALSE, label.size = 4)
autoplot(delta.pca, data=delta.T, colour="Name", label = TRUE, label.size = 4)
autoplot(delta.pca, data=delta.T, colour="Name2", label = FALSE, label.size = 4)
autoplot(delta.pca, data=delta.T, colour="Name2", label = TRUE, label.size = 4)

autoplot(delta.pca, data=delta.T, frame.colour="Name", frame=TRUE)
autoplot(delta.pca, data=delta.T, frame.colour="Name", frame=TRUE)
autoplot(delta.pca, data=delta.T, frame.colour="Name2", frame=TRUE)




