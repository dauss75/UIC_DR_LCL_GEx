
prep_data <- function(data){
  
  data <- data.frame(t(data))
  data$Name <- substr(rownames(data), 1, nchar(rownames(data))-7)
  data.agg <- aggregate(data[, -ncol(data)], list(data$Name), mean)
  rownames(data.agg) <- data.agg$Group.1;
  data.agg$Group.1 <- NULL
  tmp <- apply(data.agg, 2, as.numeric);
  rownames(tmp) <- rownames(data.agg)
  data_mean <- data.frame(t(tmp))
  #Diabete_HG_mean_input<-Diabete_HG_mean[,rownames(meta)]
  data_mean_input <- data_mean
  
  return(data_mean_input)
}



process_corr <- function(data, meta){
  
  # 1a. HG - collapse data
  #data = cbind(DwC_hg, DwoC_hg, NoD_hg)
  data <- data.frame(t(data))
  data$Name<-substr(rownames(data), 1, nchar(rownames(data))-7)
  data.agg<-aggregate(data[,-ncol(data)], list(data$Name), mean)
  rownames(data.agg) <- data.agg$Group.1;
  data.agg$Group.1 <- NULL
  tmp <- apply(data.agg, 2, as.numeric);
  rownames(tmp) <- rownames(data.agg)
  data_mean <- data.frame(t(tmp))
  #Diabete_HG_mean_input<-Diabete_HG_mean[,rownames(meta)]
  data_mean_input <- data_mean
  rm(data_mean, data.agg)
  
  # 1b. HG - PCA variance explained
  pca <- prcomp(data_mean_input, scale = FALSE, center = TRUE)
  
  
  # 1c. HG - correlation between PCA and covariates
  ## MAKE SURE META AND PC INFO HAVE THE SAME INDIVIDIAULS ##
  meta.names = rownames(meta)
  pca.names = rownames(pca_Diabete_HG_all_mean_input$rotation)
  mtch = match(meta.names, pca.names)
  pca$rotation = pca$rotation[mtch,]
  
  pca_cor<-correlate_pcs(pca, meta, npcs = 10, min.cor = 0)
  pca_cor<-as.data.frame(pca_cor);
  
  return(pca_cor)
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

calc_pc_cors <- function(pc1, pc2){
  
  rnames2 = rownames(pc2)
  rnames1 = rownames(pc1)
  
  mtch = match(rnames2, rnames1)
  
  pc1 = pc1[mtch,]
  
  pc_cor = matrix(0, 10,10)
  for (i in 1:10){
    for (j in 1:10){
      pc_cor[i,j] = cor(pc1[,i], pc2[,j])
    }
  }
  rownames(pc_cor) = colnames(pc_cor) = paste("PC",1:10, sep=""                                          )
  return(pc_cor)
}


plot_cor_pca_cov <- function(pca_cor, figName){
  date=Sys.Date()
  tt<-cbind(rownames(pca_cor), pca_cor);
  tt.m<-melt(tt);
  colnames(tt.m)<-c("Covariate", "PC", "Correlation")
  tt.m$Corr_signed = tt.m$Correlation
  tt.m$Correlation = abs(tt.m$Correlation)
  
  setEPS()
  postscript(file = paste(figName,"_",date,".eps",sep=""),width=11,height=7)
  
  # pdf(file = paste(figName,".pdf",sep=""),width=11,height=7)
  plt = ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = Correlation),colour = "white") +
    scale_fill_gradient(low = "white", high = "red") +
    geom_text(aes(label = round(Corr_signed, 2))) +
    theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
          axis.text.x = element_text(colour="grey20",size=10,face="bold"),
          axis.text.y = element_text(colour="grey20",size=10,face="bold"),
          axis.title.x = element_text(colour="grey20",size=10,face="bold"))
  print(plt)
  dev.off()
  print(paste("Wrote figure", figName))
  
}

plot_P_pca_cov <- function(pca_p, figName){
  
  tt<-cbind(pca_p);
  tt.m<-melt(tt);
  colnames(tt.m)<-c("Covariate", "PC", "P")
  tt.m$logP = -log10(tt.m$P)
  
  setEPS()
  postscript(file = paste(figName,"_",date,".eps",sep=""),width=11,height=7)
  
  # pdf(file = paste(figName,".pdf",sep=""),width=11,height=7)
  plt = ggplot(tt.m, aes(PC, Covariate)) + geom_tile(aes(fill = logP),colour = "white") +
    scale_fill_gradient(low = "white", high = "red") +
    geom_text(aes(label = round(logP, 1))) +
    theme(axis.title.y = element_text(colour="grey20",size=10,face="bold"),
          axis.text.x = element_text(colour="grey20",size=10,face="bold"),
          axis.text.y = element_text(colour="grey20",size=10,face="bold"),
          axis.title.x = element_text(colour="grey20",size=10,face="bold"))
  print(plt)
  dev.off()
  print(paste("Wrote figure", figName))
  
}


plot_pca <- function(pca, figName){
  
  date=Sys.Date()
  group = rep(0,nrow(pca))
  group[grep("DwC", rownames(pca))] = "DwC"
  group[grep("DwoC", rownames(pca))] = "DwoC"
  group[grep("NoD", rownames(pca))] = "NoD"
  
  plots = list()
  cnt = 1
  for (i in 1:4){
    for( j in 1:4){
      plots[[cnt]] = ggplot(pca, aes_string(x=paste0("PC",i), y=paste0("PC",j))) + geom_point(aes(colour = group, size = 1.5)) + ggtitle(paste0("PC",i," vs PC",j))
      cnt = cnt+1
    }
  }
  
  png(file = paste0(figName,"_",date, ".png"), width = 18, height = 18, units="in", res=200)
  print(plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
                  plots[[5]], plots[[6]], plots[[7]], plots[[8]],
                  plots[[9]], plots[[10]], plots[[11]], plots[[12]],
                  plots[[13]], plots[[14]], plots[[15]], plots[[16]], ncol = 4, nrow = 4))
  dev.off()
  
  print(paste("Wrote figure", figName))
}

pca.meta.regress <- function(pcs, meta, no.pcs = 6){
  ## CALCULATE THE FOR THE PAIRWISE REGRESSIONS
  ## CALCULATE HOW WELL THE PCS CAN PREDICT THE CLINICAL INFO
  mtch = match(rownames(meta), rownames(pcs))
  pcs = pcs[mtch,]
  Ps = c()
  Pred = c()
  
  for (i in 1:ncol(meta)){
    
    # print(i)
    ## DETERMINE IF THE VARIABLE (I) IS A FACTOR AND CHANGE IT TO A BINARY VARIABLE
    var = meta[, i]
    
    is.binary = FALSE
    if (is.character(var)){ var = as.numeric(as.factor(var))-1; is.binary=TRUE; print(is.binary)}
    
    data.full = as.data.frame(cbind(var, pcs))
    
    ## CALCULATE THE P-VALUE OF THE ASSOCATION BETWEEN CURRENT VARIABLES (I)
    ## AND THE FIRST 10 PCS
    Ps.tmp = c()
    
    for (j in 1:no.pcs){
      # print(paste0("J = ",j))
      if(is.binary){
        tmp.glm = glm(data.full$var ~ data.full[,j+1], family=binomial(link = "logit"))
        # summary(tmp.glm)
        tmp = anova(tmp.glm, test = "Chisq")
        Ps.tmp = rbind(Ps.tmp, tmp$Pr[2])
      }else{
        tmp = summary(lm(data.full$var ~ data.full[,j+1]))
        Ps.tmp = rbind(Ps.tmp, tmp$coefficients[2,4])
      }
    }
    
    Ps = rbind(Ps, c(names(meta)[i], Ps.tmp))
    
    ## ADJUST DATA TO ONLY HAVE THE SIGNIFICANT PCS (P<.05) ##
    ind = which(Ps.tmp <= 0.05)
    if (length(ind)==0){
      Pred = rbind(Pred, c(names(meta)[i], NA))
      next
    }
    data = cbind(var, pcs[,ind])
    data = as.data.frame(data)
    data$var = as.numeric(var)
    if (is.binary){
      
      tmp = glm(var ~ ., data=as.data.frame(data), family=binomial(link = "logit"))
      r2 = pR2(tmp)[5]    # access the model fit
      Pred = rbind(Pred, c(names(meta)[i], r2))
    }else{
      tmp = summary(lm(var ~ ., data=data))
      Pred = rbind(Pred, c(names(meta)[i], tmp$adj.r.squared))
    }
  }
  print("Pred")
  Pred = as.data.frame(Pred)
  names(Pred) = c("factor", "R.adj")
  Ps = as.data.frame(Ps)
  names(Ps) = c("Factor", names(pcs)[1:no.pcs])
  print("apply")
  Ps[,-1] = apply(Ps[,-1], 2, as.numeric)
  print("About to return")
  return(list(Ps = Ps, Pred = Pred))
}


pca.meta.regress.2var <- function(pcs, meta){
  
  ### FOR EACH PC, DETERMINE WHICH VARIABLES ARE SIGNIFICANT AT A
  ### GIVEN LEVEL, AND THEN ALL PAIR-WISE CONDITIONAL REGRESSIONS
  ### AND ALSO REPORT THE CORRELATION BETWEEN THE VARIABLES
  
  mtch = match(rownames(meta), rownames(pcs))
  pcs = pcs[mtch,]
  Ps = c()
  
  covs = meta
  for (i in 1:ncol(covs)){
    if (is.character(covs[,i])){
      covs[,i] = as.numeric(as.factor(covs[,i]))-1
    }
  }
  for (i in 1:ncol(pcs)){
    pc = pcs[,i]
    marg.sign = c()
    
    for (j in 1:ncol(covs)){
      data = as.data.frame(cbind(pc, covs[,j]))
      p = summary(lm(pc~., data=data))
      if (p$coefficient[2,4] <= .05){
        marg.sign = c(marg.sign, j)
        Ps = rbind(Ps, c(paste0("PC",i), names(meta)[j], names(meta)[j], p$coefficient[2,4], p$r.squared))
      }
    }
    
    for (j in marg.sign){
      for (k in marg.sign){
        
        if (j == k){ next }
        
        data = as.data.frame(cbind(pc, covs[,c(j,k)]))
        if (abs(cor(data[,2:3])[1,2])==1){
          Ps = rbind(Ps, c(paste0("PC",i), names(meta)[k], names(meta)[j], NA ,NA))
          next
        }
        
        p = summary(lm(pc~., data=data))
        Ps = rbind(Ps, c(paste0("PC",i), names(meta)[k], names(meta)[j], p$coefficients[3,4] ,p$r.squared))
      }
    }
  }
  
  Ps = as.data.frame(Ps)
  names(Ps) = c("PC", "Tested.Cov", "Conditioned.Cov", "P", "R2")
  
  return(Ps)
}


cor_meta <- function(meta){
  
  ## convert factors or characters to 0/1
  for (i in 1:ncol(meta)){
    
    if (is.character(meta[,i])){
      meta[,i] = as.numeric(as.factor(var))-1
    }
    
    cor.meta = cor(meta)
    
  }
}