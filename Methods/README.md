# Methods

- Describe principal components analysis and association testing between covariates
and PCs for the purpose of identifying important covariates to condition on for DE analyses.
Describe how the effect of group was removed prior to performing PCA.





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
