options(stringsAsFactors = FALSE)
source('code/r.functions/paper.figures.F.R')
source('code/r.functions/load.all.data.F.R')
library(SAJR)

anns = readRDS('Rdata/anns.Rdata')
all.anns = readRDS('Rdata/all.anns.Rdata')

for(s in names(anns)){
	print(s)
	t = readRDS(paste0('Rdata/',s,'.as.u.all.Rdata'))
	m = cbind(t$seg,cnt=apply(t$i+t$e,1,sum))
	rm(t)
	gc()
	m = split(m,paste(m$chr_id,m$start,m$stop))
	print(table(sapply(m,nrow)))
	m = m[sapply(m,nrow)==2]
	dupl = do.call(rbind,lapply(m,function(x){x=x[order(x$cnt),];data.frame(dupl.rate=x$cnt[1]/x$cnt[2],sid=rownames(x)[1])}))
	all.anns[[s]]$antisense.dupl.rate = NA
	all.anns[[s]][dupl$sid,'antisense.dupl.rate'] = dupl$dupl.rate
	anns[[s]] = all.anns[[s]][rownames(anns[[s]]),]
}

# saveRDS(anns,'Rdata/anns.Rdata')
# saveRDS(all.anns,'Rdata/all.anns.Rdata')

