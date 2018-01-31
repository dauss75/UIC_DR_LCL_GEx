#' Get the variance explained by principal components
#' A list with class "prcomp" containing: sdev, rotation, x.
variance_explained <- function(plist) {
  rotation <- as.data.frame(plist$rotation)
  variance <- plist$sdev ^ 2
  data.frame(
    Component = factor(1:ncol(rotation), levels = 1:ncol(rotation)),
    Variance = variance / sum(variance),
    CumulativeVariance = ( cumsum(variance) / sum(variance) )
  )
}

#' Plot the variance explained by the first few principal components.
plot_variance_explained <- function(plist, n = 10, cumulative = FALSE) {
  # Get the variance explained as a dataframe.
  options(scipen = 999)
  dat <- head(variance_explained(plist), n)
  dat$Variance<-as.numeric(substr(dat$Variance,1,5))

  if (cumulative) {
    # Cumulative line plot.
    ggplot(data = dat) +
      geom_line(alpha = 0.8,
                aes(x = Component, y = CumulativeVariance, group = 1)) +
      geom_point(size = 3, aes(x = Component, y = CumulativeVariance)) +
      xlab("PC") +
      ylab("Cumulative Fraction of Variance Explained") +
      geom_abline(intercept = 0, slope = 100 * 1 / ncol(rotation),
                  alpha = 0.25) +
      scale_y_continuous(limits = c(0, max(dat$CumulativeVariance)))
  } else {
    ggplot(data = dat, aes(x = Component, y = Variance, fill=Component)) +
      geom_bar(position = 'dodge', stat="identity") +
      ylab("Fraction of Variance Explained") +
      xlab("Principal Component") +
      geom_text(aes(label=Variance), vjust=-0.2)
  }
}

#' Correlate principal components with factors in another data.frame.
correlate_object <- function(df) {
  options(scipen = 999)
  # pca.r = as.data.frame(pca$rotation)[ , 1:npcs]
  df = as.data.frame(df)
  # Make all of the columns numeric, so we can run cor().
  for (col in colnames(df)) {
    if (class(df[ , col]) != "numeric" && class(df[ , col]) != "integer") {
      df[ , col] = factorToNumeric(df[ , col])
    }
  }
  cor(df)
}

#' Correlate principal components with factors in another data.frame.
correlate_pcs <- function(pca, df, npcs = 5, min.cor = 0.5) {
  options(scipen = 999)
  pca.r = as.data.frame(pca)[ , 1:npcs]
  df = as.data.frame(df)
  # Make all of the columns numeric, so we can run cor().
  for (col in colnames(df)) {
    if (class(df[ , col]) != "numeric" && class(df[ , col]) != "integer") {
      df[ , col] = factorToNumeric(df[ , col])
    }
  }

  pca.cor = cor(cbind(pca.r, df))
  result = list()
  # create a list for each PC and covariates
  for (col in colnames(pca.r)) {
    res=pca.cor[col,-grep("PC",names(pca.cor[col,]))]
    result[[col]] = res
  }
  result
}

#' Correlate principal components with factors in another data.frame.
proportional_pcs <- function(pca, df, npcs = 15, min.cor = 0) {
  options(scipen = 999)
  pca.r = as.data.frame(pca$rotation)[ , 1:npcs]
  df = as.data.frame(df)
  # Make all of the columns numeric, so we can run cor().
  for (col in colnames(df)) {
    if (class(df[ , col]) != "numeric" && class(df[ , col]) != "integer") {
      df[ , col] = factorToNumeric(df[ , col])
    }
  }
  ###instead of getting correlation, do the linear regression and get r.sequred value which which
  ###is proportional

  R<-matrix(0,ncol(df),ncol(pca.r))
  for (i in 1:ncol(df)){
     for (j in 1:ncol(pca.r)){
       R[i,j] <- max(0, summary(lm(df[,i]~pca.r[,j]))$adj.r.squared)
     }
  }

  rownames(R)<-colnames(df);colnames(R)<-colnames(pca.r)

  R
}

#' List the genes with the greatest absolute values of loading factors.
loading_values = function(pca, n_pcs = 2, n_items = 30) {
  pc_scores = data.frame(pca$x)
  result = list()
  for (i in 1:n_pcs) {
    code = paste0("PC", i)
    sx = pc_scores[order(pc_scores[, code]), ]
    result[[code]] = c(
      head(sx[, code, drop = TRUE], n_items / 2),
      tail(sx[, code, drop = TRUE], n_items / 2)
    )
  }
  result
}

#' Return a numeric vector with the levels of a factor.
#'
#' The purpose of this function is to convert a non-numeric factor to a numeric
#' factor. This is useful when you want to compute correlation with non-numeric
#' vectors.
factorToNumeric = function(xs) {
  uniqs = unique(xs)
  values = 1:length(uniqs)
  as.numeric(values[match(xs, uniqs)])
}

data_explore_by_pca <-function(x){
  options(scipen = 999)
  group<-as.factor(gsub("_","",substr(colnames(x),1,nchar(colnames(x))-12)))

  pca_tmp <- prcomp(x)
  title<-deparse(substitute(x))

  prinComp<-as.data.frame(pca_tmp$rotation); prinComp$group<-group

  p1<-ggplot(prinComp, aes(x = PC1, y = PC2, colour = group, shape=group)) + geom_point(size = 3, show.legend = FALSE) + ggtitle("PC1 vs PC2")+theme(text = element_text(size=20))
  p2<-ggplot(prinComp, aes(x = PC2, y = PC3, colour = group, shape=group)) + geom_point(size = 3, show.legend = FALSE) + ggtitle("PC2 vs PC3")+theme(text = element_text(size=20))
  p3<-ggplot(prinComp, aes(x = PC1, y = PC3, colour = group, shape=group)) + geom_point(size = 3) + ggtitle("PC1 vs PC3")+theme(text = element_text(size=20))
  p4<-ggplot(prinComp, aes(x = PC3, y = PC4, colour = group, shape=group)) + geom_point(size = 3) + ggtitle("PC3 vs PC4")+theme(text = element_text(size=20))
  pdf(file = paste(Figs,title,"_pca_explore_by_group_",date,".pdf",sep=""),width=18,height=12)
  multiplot(p1, p2, p3, p4, cols=2)
  dev.off()
}

data_explore_by_pca1 <-function(x){
  options(scipen = 999)
  group<-as.factor(gsub("_",".",substr(colnames(x),5,nchar(colnames(x))-7)))

  pca_tmp <- prcomp(x)
  title<-deparse(substitute(x))

  prinComp<-as.data.frame(pca_tmp$rotation); prinComp$group<-group

  p1<-ggplot(prinComp, aes(x = PC1, y = PC2, colour = group, shape=group)) + scale_shape_manual(values=1:nlevels(prinComp$group)) + geom_point(size = 3, show.legend = FALSE, stroke=2) + ggtitle("PC1 vs PC2")+theme(text = element_text(size=20))
  p2<-ggplot(prinComp, aes(x = PC2, y = PC3, colour = group, shape=group)) + scale_shape_manual(values=1:nlevels(prinComp$group)) + geom_point(size = 3, show.legend = FALSE, stroke=2) + ggtitle("PC2 vs PC3")+theme(text = element_text(size=20))
  p3<-ggplot(prinComp, aes(x = PC1, y = PC3, colour = group, shape=group)) + scale_shape_manual(values=1:nlevels(prinComp$group)) + geom_point(size = 3, stroke=2) + ggtitle("PC1 vs PC3")+theme(text = element_text(size=20))
  p4<-ggplot(prinComp, aes(x = PC3, y = PC4, colour = group, shape=group)) + scale_shape_manual(values=1:nlevels(prinComp$group)) + geom_point(size = 3, stroke=2) + ggtitle("PC3 vs PC4")+theme(text = element_text(size=20))
  pdf(file = paste(Figs,title,"_pca_explore_by_group_",date,".pdf",sep=""),width=18,height=12)
  multiplot(p1, p2, p3, p4, cols=2)
  dev.off()
}

data_explore_by_pca2 <-function(x){
  options(scipen = 999)
  group<-as.factor(gsub("_","",substr(colnames(x),1,nchar(colnames(x))-5)))

  pca_tmp <- prcomp(x)
  title<-deparse(substitute(x))

  prinComp<-as.data.frame(pca_tmp$rotation); prinComp$group<-group

  p1<-ggplot(prinComp, aes(x = PC1, y = PC2, colour = group, shape=group)) + geom_point(size = 3, show.legend = FALSE) + ggtitle("PC1 vs PC2")+theme(text = element_text(size=20))
  p2<-ggplot(prinComp, aes(x = PC2, y = PC3, colour = group, shape=group)) + geom_point(size = 3, show.legend = FALSE) + ggtitle("PC2 vs PC3")+theme(text = element_text(size=20))
  p3<-ggplot(prinComp, aes(x = PC1, y = PC3, colour = group, shape=group)) + geom_point(size = 3) + ggtitle("PC1 vs PC3")+theme(text = element_text(size=20))
  p4<-ggplot(prinComp, aes(x = PC3, y = PC4, colour = group, shape=group)) + geom_point(size = 3) + ggtitle("PC3 vs PC4")+theme(text = element_text(size=20))
  pdf(file = paste(Figs,title,"_pca_Dwc_DwoC_explore_by_group_",date,".pdf",sep=""),width=18,height=12)
  multiplot(p1, p2, p3, p4, cols=2)
  dev.off()
}
