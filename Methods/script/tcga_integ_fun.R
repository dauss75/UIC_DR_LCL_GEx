#::: rank by MAD or CV :::
rankByMAD <- function(x) {
    ## x is a matrix and the last column is gene list 
    ## that will be used later for removing correlated transcript
    data<-x[,-ncol(x)]
    genes<-x[,ncol(x)]
    rowWise<-1; colWise<-2
    x1<-apply(data,colWise,as.numeric)
    rownames(x1)<-rownames(x)
    x1<-data.frame(x1)
    x1$mad<-apply(x1,rowWise,mad)
    x1$genes<-genes
    #     rownames(x1)<-rownames(x)
    madScrClmn<-ncol(x1)-1
    z <- x1[order(-x1[,madScrClmn]),] # sort by mad score
    # z$mad<-NULL
    return(z)
}

rankByMAD_prad <- function(x) {
    ## x is a matrix where row is for genes and column is for samples
    rowWise<-1; colWise<-2
    x1<-apply(x,colWise,as.numeric)
    rownames(x1)<-rownames(x)
    x1<-data.frame(x1)
    x1$mad<-apply(x1,rowWise,mad)
#     x1$genes<-genes
    #     rownames(x1)<-rownames(x)
#     madScrClmn<-ncol(x1)
    z <- x1[order(-x1[,ncol(x1)]),] # sort by mad score
    return(z)
}

# ::: remove isoforms :::
rmIsoform <-function(x){
    ## remove isoforms that are correlated
    ## the input matrix x contains genes and mad scores 
    ## sort first by gene, remove correlated isoforms, 
    gene.sorted<-x[with(x, order(genes)), ] # sort by genes
    genLen = nrow(gene.sorted) ## number of genes
    mad_gene_clm<-2
    smplsz = ncol(gene.sorted)-mad_gene_clm   ## number of samples
    gene.sorted$keep<-rep(1,genLen)
    
    indxMAD<-smplsz + 1 ## the index of mad column
    indxKp<-smplsz + 3  ## the index of keep column
    gnName<-gene.sorted$genes[1]
    
    # gnName<-mad.ranked.Log10Mtrx$genes[1]
    LeadGnIndx<-1
    StrtGnIndex<-1
    
    while (LeadGnIndx<=genLen) {
        if (gene.sorted$genes[LeadGnIndx]==gnName) {
            LeadGnIndx=LeadGnIndx+1
        }
        else {
            if (StrtGnIndex!=LeadGnIndx) {
                crrMtrx<-cor(t(gene.sorted[StrtGnIndex:(LeadGnIndx-1),1:smplsz])) 
                for ( k in 1:(nrow(crrMtrx)-1) ) {
                    crrIsofrms<-colnames(crrMtrx)[crrMtrx[k,k:nrow(crrMtrx)] > 0.8]
                    if (length(crrIsofrms)>1) {
                        crrIsofrms<-crrIsofrms[-which.max(gene.sorted[crrIsofrms,indxMAD])]
                        gene.sorted[crrIsofrms,indxKp]=-1
                    }
                }
            }
            gnName<-gene.sorted$genes[LeadGnIndx]
            StrtGnIndex<-LeadGnIndx
            LeadGnIndx<-LeadGnIndx+1    
        }  
    }
    sum(gene.sorted$keep==1) 
    
    fltrd.data <-gene.sorted[gene.sorted$keep==1,]
    madScrClmn<-ncol(fltrd.data)-mad_gene_clm
    fltrd.data <- fltrd.data[order(-fltrd.data[,madScrClmn]),] # sort by mad score
    fltrd.data<-fltrd.data[,1:smplsz]
    return(fltrd.data)
}

extract_survival_info <-function(x){
    # x is the TCGA clinical_patient_cancer.txt file
    # this step will remove samples that do not meet criteria such as days_to_deat
    cln_list<-c("bcr_patient_barcode","vital_status", "days_to_last_followup","days_to_death", "pathologic_stage")
    colnames(x)<-x[1,]; x<-x[-1:-2,cln_list]
    status_1<-c("[Not Applicable]","[Not Available]")
    y<-x[!((x$days_to_last_followup %in% status_1) & (x$days_to_death %in% status_1)),] # remove when both columns match to status_1
    y$days_to_death[y$days_to_death %in% status_1] <-y$days_to_last_followup[y$days_to_death %in% status_1]         # replace days_to_death by days_to_last_followup
    y$vital_status[y$vital_status == "Alive"]<-FALSE
    y$vital_status[y$vital_status == "Dead"]<-TRUE
    y<-y[y$days_to_death > 0,] 
    rownames(y)<-y$bcr_patient_barcode
    y$days_to_death <- as.numeric(y$days_to_death); y$vital_status <- as.logical(y$vital_status)
    return(y)
}

extract_prad_survival_info <-function(x){
    library(Hmisc)
    # x is the TCGA clinical_patient_cancer.txt file
    # this step will remove samples that do not meet criteria such as days_to_deat
    cln_list<-c("bcr_patient_barcode","vital_status", "days_to_last_followup","days_to_death", "primary_pattern", "secondary_pattern", "tertiary_pattern", "gleason_score", "psa_value", "clinical_M", "clinical_N", "clinical_T", "pathologic_M", "pathologic_N", "pathologic_T")
    colnames(x)<-x[1,]; 
    x<-x[-1:-2,cln_list]
    status_1<-c("[Not Applicable]","[Not Available]", "[Unknown]")
    y<-x[!((x$days_to_last_followup %in% status_1) & (x$days_to_death %in% status_1)),] # remove when both columns match to status_1
    y$days_to_death[y$days_to_death %in% status_1] <-y$days_to_last_followup[y$days_to_death %in% status_1]         # replace days_to_death by days_to_last_followup
    y$vital_status[y$vital_status == "Alive"]<-FALSE
    y$vital_status[y$vital_status == "Dead"]<-TRUE
    y<-y[y$days_to_death > 0,] 
    rownames(y)<-y$bcr_patient_barcode
    y$days_to_death <- as.numeric(y$days_to_death); y$vital_status <- as.logical(y$vital_status)

    for (i in 1:nrow(y)){
        if ((y$pathologic_M[i] %in% status_1) & (y$clinical_M[i] %nin% status_1)) {
            y$pathologic_M[i] <-y$clinical_M[i]
        }
        if ((y$pathologic_N[i] %in% status_1) & (y$clinical_N[i] %nin% status_1)) {
            y$pathologic_N[i] <-y$clinical_N[i]
        }
        if ((y$pathologic_T[i] %in% status_1) & (y$clinical_T[i] %nin% status_1)) {
            y$pathologic_T[i] <-y$clinical_T[i]
        }
    }
    y$clinical_M<-NULL; y$clinical_N<-NULL; y$clinical_T<-NULL
    return(y)
}

extract_prad_cqcf_survival_info <-function(x){
    # x is the TCGA clinical_patient_cancer.txt file
    # this step will remove samples that do not meet criteria such as days_to_deat
    cln_list<-c("bcr_patient_barcode","gleason_score_primary", "gleason_score_secondary","gleason_score_combined", "race")
    colnames(x)<-x[1,]; 
    x<-x[-1:-2,cln_list]
    #     status_1<-c("[Not Applicable]","[Not Available]")
    #     y<-x[!((x$days_to_last_followup %in% status_1) & (x$days_to_death %in% status_1)),] # remove when both columns match to status_1
    #     y$days_to_death[y$days_to_death %in% status_1] <-y$days_to_last_followup[y$days_to_death %in% status_1]         # replace days_to_death by days_to_last_followup
    #     y$vital_status[y$vital_status == "Alive"]<-FALSE
    #     y$vital_status[y$vital_status == "Dead"]<-TRUE
    #     y<-y[y$days_to_death > 0,] 
    rownames(x)<-x$bcr_patient_barcode
    #     y$days_to_death <- as.numeric(y$days_to_death); y$vital_status <- as.logical(y$vital_status)
    return(x)
}

coxph_p_value <-function(x, y){
    # last two columns of the matrix x are vital_status and days_to_death and y is a cutoff value
    p_values<-NULL
    p_value_clm<-3
    print(paste("process", ncol(x)-2, "number of lines", sep=" "))
    for (i in 1:(ncol(x)-2)){   
        form <- as.formula(paste("Surv(days_to_death,vital_status) ~ ", colnames(x)[i], sep = ""))
        res<-summary(coxph(form, data=x))
        p_values<-c(p_values,res$logtest[p_value_clm]) # test, df, pvalue
        print(paste("line", i, sep= " "))
    }
    print("coxph is done")
    names(p_values)<-colnames(x)[1:(ncol(x)-2)] # -2 is exclude vital_status and days_to_death
    p_sorted<-sort(p_values,decreasing=FALSE)
    write.table(p_sorted,file="LUAD_p_values_sorted_.csv",sep="\t")
    p_sorted_with_cutoff<-p_sorted[p_sorted < y]
    return(p_sorted_with_cutoff)
    #     z<-t(as.matrix(x[,names(p_sorted_with_cutoff)]))
    #     return(z)
}