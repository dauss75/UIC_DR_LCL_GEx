cd "/Users/dcao/Box Sync/Current Consulting Projects/Mike Grassi/Genomics"

***********************************
* FDR 
***********************************

use treat_test, clear
sort geneID
save treat_test, replace

use anova_result3, clear
sort geneID
merge geneID using treat_test

tab _merge

drop _merge

replace df_t = 21

foreach i in bonferroni sidak holm holland liu1 liu2 hochberg simes yekutieli krieger {
	multproc, pvalue(p_t) method(`i')
}

foreach i in bonferroni sidak holm holland liu1 liu2 hochberg rom simes yekutieli krieger {
	multproc, pvalue(p1) method(`i')
}
foreach i in bonferroni sidak holm holland liu1 liu2 hochberg rom simes yekutieli krieger {
	multproc, pvalue(p2) method(`i')
}
foreach i in bonferroni sidak holm holland liu1 liu2 hochberg simes yekutieli krieger {
	multproc, pvalue(pvalue) method(`i')
}

gen mse1 = ss1/df1
gen mse2 = ss2/df2

list geneID symbol mse1 mse2 pvalue p1 p2 if pvalue < 0.05
