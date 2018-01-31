# Methods

- principal components analysis (PCA) and association testing between covariates
and PCs for identifying important covariates to condition on for differential expression (DE) analyses.
Describe how the effect of group was removed prior to performing PCA.

   - we run PCA using the R PCA function [prcomp](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/prcomp.html) with metadata by: (1) using all the data and (2) after the group removed such that we can identify which covariate is the confounding factor.

      `prcomp(data_input, scale = FALSE, center = TRUE)`

         - 1. as a sanity check, explore the variance explained for each PC to look at the distribution.
         - 2. compute Pearson's correlation coefficient (r) between eigenvector and metadata. The p-value of 0.05 is considered as significant.
            `cor(pca$rotation, metadata)`
         - 3. Now, remove the group by retrieving the residual data from applying the linear regression model that yields the difference between the gene expression data of the dependent variable group (g) and the fitted values (g')
            `summary(lm(data_input~group))$residual`
         - 4. Then, repeat the step 2





- Describe assessing correlation between covariates as well as performing conditional regression
to determine if any clinical covariates can be eliminated as a necessary predictor in DE models.
Covariates: All analyses will include PC1 (with group removed) as a covariate and growth rate.

- Describe the model and test that will be used to test for DE genes. Describe the threshold that
will be used determine significance (FDR levels). Address how the repeated measures are being
treated for analysis. Collapsing into a mean vs. fitted into a mixed model that takes into account
the correlation between repeated measures.

For genes found to be interesting based on the significance criteria above, we will investigate
if any of the association might be explained by a covariate that was not included in the model.
These include: A1C, Sex, Age, Duration, EBV copy number, BMI, LDL, HDL, Pulse, SBP, DBP, Triglycerides

Describe how results will be summarized. E.g. volcano plots, etc.

Summarize significant genes found in each group and how to determine which are shared and unique.
a.    Consider more liberal thresholds for allowing genes to be DE in different groups.
 i.         Require same direction?
ii.         Bonferroni threshold.
iii.         Calculate FDR in significant genes.

Describe gene set enrichment analysis (GSEA)
