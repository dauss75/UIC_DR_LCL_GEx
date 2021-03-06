# Methods

## Table of Contents

- [PCA and Covariates](#pca-and-covariates)
- [DE models with and without collapsing data](#de-models-with-and-without-collapsing-data)
- [DE models with clinical covariates controlled](#de-models-with-clinical-covariates-controlled)
- [gene set enrichment analysis (GSEA)](#gene-set-enrichment-analysis-gsea)
- [Predixcan](#predixcan)

## PCA and Covariates

To identify important covariates to condition on the differential expression (DE) analyses, we performed principal components analysis (PCA) and association testing between covariates and PCs .
We describe here how the effect of group (i.e. PDR, nPDR) was removed prior to performing the PCA.

   1. run PCA using the PCA function [prcomp](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/prcomp.html) in R with metadata by: (1) using all the gene expression data and (2) after the patients group removed such that we can identify which covariate is the confounding factor.

      `pca<-prcomp(input_data, scale = FALSE, center = TRUE)`

   2. as a sanity check, explore the variance explained for each PC to look at the distribution.

   ![Screenshot](figure/figure1.png "an example of variance explained for the high glucose data")

   3. compute Pearson's correlation coefficient (r) between eigenvector and metadata. The p-value of 0.05 (1.3 in log scale) is considered as significant.

      `pca_meta_cor<-cor(pca$rotation, metadata)`

    ![Screenshot](figure/figure2.png)

   4. Now, remove the group by retrieving the residual data from applying the linear regression model that yields the difference between the gene expression data of the dependent variable group (g) and the fitted values (g').

      `input_data_group_removed<-summary(lm(input_data~group))$residual`

   5. Then, repeat the step in ii and iii.
   ![Screenshot](figure/figure3.png)
   ![Screenshot](figure/figure4.png)

## DE models with and without collapsing data

We evaluated the gene expression data in two different ways for DE analysis using (1) collapsing into a mean and (2) all replicates that fitted into a mixed model taking into account the correlation between repeated measures. We show an example with PDR that compares the gene expression between high glucose and stand glucose.

   1. with a mean data
      - build a disign matrix using the [model.matrix](https://www.rdocumentation.org/packages/stats/versions/3.4.3/topics/model.matrix) function for paired treatment.
         ```
         subject="avg_PDR"
         avg_PDR<-cbind(avg_PDR_hg,avg_PDR_sg)
         targets<-readTargets(paste("avg_hg_sg_",subject,"_target.txt", sep=''))
         Paired <- factor(targets$paired)
         Treat <- factor(targets$Treatment)
         design <- model.matrix(~Paired+Treat)
         fit <- lmFit(avg_PDR, design)
         fit<- eBayes(fit)
         ```
         ![Screenshot](figure/Figure9.png "design matrix")

         ![Screenshot](figure/figure10.png "overall p-value distribution")

         ![Screenshot](figure/figure11.png "volcano plot")

   2. with replicates
      - build a design matrix using the [model.matrix](https://www.rdocumentation.org/packages/stats/versions/3.4.3/topics/model.matrix) function with both replicates (purple box) and treatment (red box).
      Note that "C" and "T" in treatment stand for high glucose and standard glucose, respectively.

         ```
         subject="PDR_replicate”
         PDR<-input_data[,grep("PDR",colnames(input_data))]
         targets<-readTargets(paste("hg_sg_",subject,"_target.txt", sep=''))
         Treat <- factor(targets$Treatment,levels=c("C","T"))
         Replicates <- factor(targets$rep)
         design <- model.matrix(~Replicates+Treat)
         ```
        ![Screenshot](figure/figure6.png "design matrix")

        Then, use the [duplicateCorrelation](http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/dupcor.html) function to estimate the correlation between technical replicates using a mixed linear model that returns a consensus correlation, a robust average of the individual correlation. We use the the value to [lmFit](http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/lmFit.html).

         ```
         corfit <- duplicateCorrelation(PDR, block = targets$Subject)
         fit <-lmFit(PDR, design, block=targets$Subject, correlation=corfit$consensus.correlation)
         fit<-eBayes(fit)
         ```

         ![Screenshot](figure/figure7.png "overall p- and q-value distribution")


         ![Screenshot](figure/figure8.png "Volcano plot for DE genes")

- __conclusion:  we find that the gene rank from both analysis is about the same except the statistical power is much stronger with the data using all replicates.__


## DE models with clinical covariates controlled

Once the covariates (PC1 and growth rate) are identified, we eliminate them as a necessary predictor in DE models.
   - get PC1 after group regressed out as well as growth rate from metadata.

      ```
      pca_tmp<-pca$rotation
      pc1<-pca_tmp$PC1
      growth_rate<-meta[,"GROWTH_RATE"]
      ```

   - Using limma, the Bioconductor R pcakage:
      - build a design matrix.

         `design<-model.matrix(~pc1+growth_rate+group)`

      - Given a linear model fit, compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical Bayes moderation of the standard errors towards a common value.


         `fit <- eBayes(lmFit(input_data, design))`

      - Then, apply FDR < 5% and log2(fold change) > 1 (or fold change > 2).

         ```
         qval.cutoff=0.05;
         FC.cutoff=2
         x1=topTable(fit, coef="group2", n=nrow(genes), p.value=qval.cutoff, adjust.method="BH", genelist=genes)
         ```

   - This gives 25 genes for difference between nDM vs. nPDR in high glucose.

   ![Screenshot](figure/figure5.png)

## gene set enrichment analysis (GSEA)

This GSEA analysis was performed by using preranked gene list in the [javaGSEA Desktop Application](https://github.com/GSEA-MSigDB/gsea-desktop) tool, available at http://software.broadinstitute.org/gsea/downloads.jsp, along with the c2 (curated gene sets) and c5 (gene ontology (GO) gene sets) gene set collections from the [Molecular Signature Database](http://software.broadinstitute.org/gsea/msigdb/index.jsp). Specific steps are as follows:

- Input: ranked genes by this formula: sign(Fold Change) x - log10(p-value).  Full gene list is provided and duplicated genes are removed resulting 11579 genes
- Parameter setting
   - permutations are done by gene set
   - The runs are based on weighted (ranked).
   - Gene sets database
      - c2.all.v6.0, curated from various sources such as online pathway databases, the biomedical literature, and knowledge of domain experts
      - c5.all.v6.0,  based on GO terms in the collection belong to one of three GO ontologies: molecular function (MF), cellular component (CC) or biological process (BP)
   - Min and max gene set size: 15 and 500

   ![Screenshot](figure/figure12.png)

-------------
## Predixcan

- We download the predict DB data ([v7](GTEx-V7_HapMap-2017-11-29.tar.gz)) from [here](http://predictdb.hakyimlab.org/).
- Then, we run predixcan for prediction and association using the following command:

```
./PrediXcan.py --predict --assoc --linear \
               --weights weights/file_name.db \
               --dosages genotype \
               --samples samples.txt \
               --pheno phenotype/igrowth.txt \
               --output_prefix results/file_name
```

Here is a real example we used:

   '''
   ./PrediXcan.py --predict --assoc --linear --weights weights_v7/gtex_v7_Adipose_Subcutaneous_imputed_europeans_tw_0.5_signif.db --dosages genotype --samples sample.txt  --pheno phenotype/igrowth.txt --output_prefix results_v7/gtex_v7_Adipose_Subcutaneous
   '''
