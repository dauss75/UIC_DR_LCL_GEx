## CV:coefficient of variation
CV <- function(x,na.rm=TRUE) {

  x1<-apply(x,2,as.numeric)
  x1_SD<-apply(x1,1,sd)
  x1_median<-apply(x1,1,median)  
  y1<-x1_SD/x1_median
  return(y1)
}

# 1. load data


# 2. run CV to pick top 5000 genes



# 3. run PCA with covariates



# 4. plot