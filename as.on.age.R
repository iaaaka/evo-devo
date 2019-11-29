setwd('~/skoltech/projects/evo.devo/')
options(stringsAsFactors = FALSE)
source('code/r.functions/load.all.data.F.R')
source('code/r.functions/as.on.age.F.R')
library(SAJR)

species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
anns = readRDS('Rdata/anns.Rdata') 
psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
border.stages = readRDS('Rdata/border.stages.Rdata')

#complexity on age


col = unique(meta[,c('tissue','col')])
col = setNames(col$col,substr(col$tissue,1,1))
col = c(col,const='black',bh='orange','2ts'='#444444','3ts'='#777777','4ts'='#AAAAAA','5ts'='#DDDDDD')


pdf('figures/AS.ad.complexity.on.age.pdf',w=12,h=7)
par(mfrow=c(2,4),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(4,2,1.5,0),oma=c(0,0,1.5,1))
plotSpeciesADAltProp('mouse',main='All cassette',col=col,psi.thr = 0.1)
plotSpeciesADAltProp('mouse',TRUE,main='Cassettes with dPSI > 0.5',col=col,psi.thr=0.1)
mtext('mouse',3,outer = T)
plotSpeciesADAltProp('human',main='All cassette',col=col)
plotSpeciesADAltProp('human',TRUE,main='Cassettes with dPSI > 0.5',col=col)
mtext('human',3,outer = T)
dev.off()


getAlternativityOnAge = function(psi,m,tissues=unique(m$tissue),rm.na = c('no','obs','seg')[1],uniq=FALSE,ret.segments=FALSE){
	m = m[rownames(m) %in% colnames(psi) & m$tissue %in% tissues,]
	psi = psi[,rownames(m)]
	p = ifelse(is.na(psi),'-',ifelse(psi<0.1,0,ifelse(psi>0.9,1,'a')))
	s = unique(m[colnames(psi),c('days','stage')])
	s = s[order(s$days),'stage']
	sapply(s,function(x){
		p = p[,m$stage==x]
		if(rm.na == 'seg')
			p = p[apply(p=='-',1,sum)==0,]
		if(uniq)
			t = apply(p,1,function(y)paste(sort(unique(y)),collapse=''))
		else
			t = apply(p,1,paste,collapse='')
		if(rm.na == 'obs')
			t = gsub('-','',t,fixed=TRUE)
		if(ret.segments)
			t
		else
			table(t)
		})
}

m = getAgeASchanges(psi.tsm,meta.tsm,0.5,border.stages,'mouse')
mm = apply(m,2,function(x){x %in% c('u','d') & anns$mouse$sites=='ad'})

z = getAlternativityOnAge(psi.tsm$mouse[anns$mouse$sites=='ad' & anns$mouse$cod !='n'& apply(mm,1,sum)>0,],meta.tsm,c('brain','heart','liver','ovary','testis'),'seg',uniq=T)
barplot(sapply(z,'[',c('00000','aaaaa','11111')))
t = sapply(z,sum) - sapply(lapply(z,'[',c('0','a','1')),sum)
barplot(rbind(t,sapply(z,'[',c('0','a','1'))))


barplot(sapply(z,'[',c('a0000','0a000','00a00','000a0','0000a')),col=params$tissue.col[c('brain','heart','liver','ovary','testis')])



z = getAlternativityOnAge(psi.tsm$mouse[anns$mouse$sites=='ad' & anns$mouse$cod !='n' & apply(mm,1,sum)>0,],meta.tsm,c('brain','heart','liver','ovary','testis'),'no',uniq=F,ret.segments = T)
sort(table(z[z[,1]=='11111',14]))
sort(table(z[z[,1]=='00000',14]))
sort(table(z[z[,1]=='aaaaa',14]))
