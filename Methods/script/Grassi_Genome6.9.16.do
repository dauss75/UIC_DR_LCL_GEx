cd "/Users/dcao/Box Sync/Current Consulting Projects/Mike Grassi/Genomics"
set more off
insheet using "raw_delta_gene_anno.csv", clear
drop if _n == 1
ren v1 genecode
ren v2 symbol
ren v3 name

sort genecode
save GeneCodeSymbol, replace

insheet using "raw_delta.csv", clear
sort genecode
merge genecode using GeneCodeSymbol
tab _merge
drop _merge

encode genecode, gen(geneID)

foreach i in 1026  {
	foreach j in 1 2 3 4 5{
		ren dwc_`i'_30mm`j' S`i'D`j'
	}
}

foreach i in 6009 11167 13206 17075 21127 26008 27141 {
	foreach j in 1 2 3 {
		ren dwc_`i'_30mm`j' S`i'D`j'
	}
}

foreach i in 2318  25224{
	foreach j in 1 2 3 4 5{
		ren dwoc_`i'_30mm`j' S`i'D`j'
	}
}

foreach i in 3395 6154 16362 18296 21183 {
	foreach j in 1 2 3 {
		ren dwoc_`i'_30mm`j' S`i'D`j'
	}
}

foreach i in 7012 7344 11985 14381 14520 14569 14581 {
	foreach j in 1 2 3 {
		ren nod_`i'_30mm`j' S`i'D`j'
	}
}

reshape long S1026D S6009D S11167D S13206D S17075D S21127D S26008D S27141D S2318D S25224D S3395D S6154D S16362D S18296D S21183D S7012D S7344D S11985D S14381D S14520D S14569D S14581D, i(geneID) j(rep)

foreach i in 1026 6009 11167 13206 17075 21127 26008 27141 2318  25224 3395 6154 16362 18296 21183 7012 7344 11985 14381 14520 14569 14581 {
	ren S`i'D Delta`i'
}

reshape long Delta, i(geneID rep) j(subID)

save Deltalong, replace

use Deltalong, clear

preserve

	bysort geneID subID: egen mDelta = mean(Delta)
	bysort geneID subID: keep if _n == 1
	keep geneID subID mDelta
	save mDelta, replace
	gen df_t = .
	gen p_t = .
	gen mu_t = .
	gen se_t = .
	gen t = .
	forvalues i = 1/21975 {
		dis "ttest gene `i'"
		quietly ttest mDelta == 0 if geneID == `i'
		quietly replace mu_t = `r(mu_1)' if geneID == `i'
		quietly replace se_t = `r(se)' if geneID == `i'
		quietly replace t = `r(t)' if geneID == `i'
		quietly replace df_t = `r(df_t)' if geneID == `i'
		quietly replace p_t = `r(p)' if geneID == `i'
	}
	bysort geneID: keep if _n == 1
	keep geneID mu_t se_t t df_t p_t
	order geneID mu_t se_t t df_t p_t
	save treat_test, replace
restore

gen ss1 = .
gen ss2 = .
gen df1 = .
gen df2 = .

gen F1 = .
gen F2 = .

gen rss = .
gen df_r = .
drop if rep > 3
forvalues i = 1/21975 {
	dis "anova gene `i'"
	quietly anova Delta subID rep if geneID == `i', repeated(rep)
	quietly replace ss1 = `e(ss_1)' if geneID == `i'
	quietly replace ss2 = `e(ss_2)' if geneID == `i'
	quietly replace rss = `e(rss)' if geneID == `i'
	quietly replace df1 = `e(df_1)' if geneID == `i'
	quietly replace df2 = `e(df_2)' if geneID == `i'
	quietly replace df_r = `e(df_r)' if geneID == `i'

	quietly replace F1 = `e(F_1)' if geneID == `i'
	quietly replace F2 = `e(F_2)' if geneID == `i'

}
bysort geneID: keep if _n == 1

gen p1 = 1-F(df1,df_r, F1)
gen p2 = 1-F(df2,df_r, F2)
gen fvalue = (ss1/df1)/(ss2/df2)
gen pvalue = (1-F(df1, df2, fvalue))

keep  geneID genecode symbol name ss1 ss2 rss df1 df2 df_r  F1 F2 fvalue p1 p2 pvalue 
order geneID genecode symbol name ss1 ss2 rss df1 df2 df_r  F1 p1 F2 p2 fvalue pvalue 
save anova_result3, replace


