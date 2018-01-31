gr<-read.csv(paste(PhenotypeDir,"cell_count.csv", sep=''))
rownames(gr)<-gr$t;gr$t<-NULL
gr_log<-log(gr)
gr_log.t<-as.data.frame(t(gr_log))
gr_log.t$day<-c(0,3,5,7)

varlist<-names(gr_log.t)[1:ncol(gr_log.t)-1]
models <- lapply(varlist, function(x) {
  lm(substitute(i ~ day, list(i = as.name(x))), data = gr_log.t)
})

gr_log$rate<-NA
gr_log$r2<-NA
gr_log$pval<-NA
for (i in 1:length(models)){
  gr_log$rate[i]<-models[[i]]$coefficients[2]
  gr_log$r2[i]<-summary(models[[i]])[8]
  gr_log$pval[i]<-summary(models[[i]])[[4]][2,4]
}

my.df <- data.frame(lapply(gr_log, as.numeric), stringsAsFactors=FALSE)
rownames(my.df)<-rownames(gr_log)
write.csv(my.df, paste(PhenotypeDir,"growth_rate.csv",sep=""))

png(paste(Figs,"dummy1.png",sep=""),res=100)
