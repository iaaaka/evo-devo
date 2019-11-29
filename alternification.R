setwd('~/skoltech/projects/evo.devo/')
options(stringsAsFactors = FALSE)
source('code/r.functions/load.all.data.F.R')
library(SAJR)

species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
#orth.seg.ad.all = readRDS('Rdata/orth.seg.ad.all.Rdata')
orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
my.ge.cod = readRDS('Rdata/my.ge.cod.Rdata')
#all.anns = readRDS('Rdata/all.anns.Rdata')
#my.ge.cod = readRDS('Rdata/my.ge.cod.Rdata')
orth.seg.ad.all.tsm = readRDS('Rdata/orth.seg.ad.all.tsm.Rdata')

params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)


# check on  PSI thr ####
psi = list()
psi$mean = sapply(orth.seg.ad.all.tsm,function(x)apply(x,1,mean,na.rm=T))
psi$min = sapply(orth.seg.ad.all.tsm,function(x)apply(x,1,min,na.rm=T))
psi$min4 = sapply(orth.seg.ad.all.tsm,function(x)apply(x,1,function(z)sort(z)[4]))
psi$ann = sapply(orth.seg.ad.all,function(x)x$seg$type)

boxplot(psi$min4[,1] ~psi$ann[,1])
table(psi$min4[,1]<0.9 & psi$min4[,1]>0.1,psi$ann[,1])



d = psi$min4

min4.spsp = sapply(c(1,0.98,0.95,0.9,0.85,0.8),function(t){
	f = apply(is.na(d) | is.infinite(d),1,sum)==0
	apply(d[f,]<t,1,function(s)paste(species$short[s],collapse = ''))
	})
u = unique(as.vector(min4.spsp))
t = apply(min4.spsp,2,function(c)table(factor(c,levels = u)))
t3 = t[nchar(rownames(t))==3,]
t3[order(t3[,1],decreasing = T),]

plot(t3['mrb',]/apply(t3[rownames(t3) != 'mrb',],2,mean))

plot(table(nchar(min4.spsp[,6]))[-1])

t2 = t[nchar(rownames(t))==2,]
t2[order(t2[,1],decreasing = T),]

f = rownames(min4.spsp)[min4.spsp[,6]=='hqmrboc']
table(f,apply(psi$ann=='ALT',1,sum))

p = do.call(cbind,lapply(orth.seg.ad.all.tsm, function(x)x[f,]))
c = cor(p,u='p')
mds = cmdscale(1-c,k=2)
m = meta.tsm[rownames(mds),]
plot(mds,col=m$col)
plot(mds,col=params$species.col[m$species])

f =apply(psi$ann=='ALT',1,sum)==7
table(f)
p1 = do.call(cbind,lapply(orth.seg.ad.all.tsm, function(x)x[f,]))
c1 = cor(p1,u='p')
mds1 = cmdscale(1-c1,k=2)
m = meta.tsm[rownames(mds),]
plot(mds1,col=m$col)
plot(mds1,col=params$species.col[m$species])

# look on 6th exon of FAS ####
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')
exon.birth.one = readRDS('Rdata/exon.birth.one.Rdata')
seg2ens = readRDS('Rdata/seg2ens.Rdata')
psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
anns = readRDS('Rdata/anns.Rdata')
h = orth.seg.ad$human$seg
h = anns$human
h[h$chr_id=='10' & h$start<=90770572 & h$stop >= 90770510, ]
plotTissueAgeProile(psi.tsm$human['hum.8192.s10',],meta.tsm)
seg2ens$human['hum.8192.s10'] # correct, ENSG00000026103 is FAS
# there are no ortholog in opossum according to ensembl
'hum.8192.s10' %in% do.call(rbind,exon.birth.one)$seg_id #so, this exon is not newborn
hqmrb = read.table('processed/orth.segs/only.ad/hqmrb.0.0.orth.ad.segs')
sids = unlist(hqmrb[hqmrb$V2=='hum.8192.s10',2:6]) # rabit have bad overlap (compare to hqmrb.0.6.orth.ad.segs)
per.tissue.age.qv$human[sids[1],]
per.tissue.age.qv$macaque[sids[2],]
r = list()
for(i in 1:5){
	print(i)
	r[[i]] = readRDS(paste0('Rdata/',rownames(species)[i],'.as.u.all.Rdata'))[sids[i],]
	gc()
	#
}
#saveRDS(r,'Rdata/hum.8192.s10.FAS.orth.ir.Rdata')
r = readRDS('Rdata/hum.8192.s10.FAS.orth.ir.Rdata')

pdf('figures/FAS.6th.exon.hum.8192.s10.pdf',w=10,h=6)
par(mfrow=c(2,3),tck=-0.02,mgp=c(1.3,0.4,0),mar=c(5,3,1.5,0),oma=c(0,0,2,1))
for(i in 1:5)
	plotTissueAgeProile(r[[i]]$ir[1,],meta,ylim=c(0,1),main=rownames(species)[i],age.axis = 'rank',ylab='PSI')
dev.off()

# old ####
orth.seg.ad.all.tsm = lapply(orth.seg.ad.all,function(x){x = x[,colnames(x$ir) %in% rownames(meta)];m=meta[colnames(x$ir),];calcMeanCols(x$ir,paste(m$species,m$tissue,m$stage))})
#saveRDS(orth.seg.ad.all.tsm,'Rdata/orth.seg.ad.all.tsm.Rdata')
orth.seg.ad.all.tsm = readRDS('Rdata/orth.seg.ad.all.tsm.Rdata')
alt.thrs = c(0.9,0.9)
min.psi = sapply(orth.seg.ad.all.tsm,function(x)apply(x,1,min,na.rm=T))
max.na.fraq = sapply(orth.seg.ad.all.tsm,function(x)apply(is.na(x),1,mean))
hist(apply(max.na.fraq,1,max))
f = apply(max.na.fraq,1,max)<0.5
table(f)
hist(min.psi[f,1],200)

alt.sp = apply(min.psi,1,function(x){
	r = '-'
	if(sum(is.na(x)) == 0 && sum(x>alt.thrs[1] & x<alt.thrs[2])==0){
		r=paste(species$short[x<alt.thrs[1]],collapse='')
	}
	r
	})
table(alt.sp!='-',f)
sort(table(alt.sp[f]))

orth.seg.ad.all.ids = sapply(orth.seg.ad.all,function(x)rownames(x$seg))
orth.seg.ad.all.ids = cbind(as.data.frame(orth.seg.ad.all.ids),species=alt.sp)


pdf('figures/alternification.exon.PSI.on.evol.ages[0.8,0.9].pdf',w=12,h=5)
par(mfcol=c(2,5),tck=-0.01,mgp=c(1.1,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,2,1))
m=plotExpAndPsiForNewExons('mouse',c('m','mr','mrb','hqmrb','hqmrbo'),orth.seg.ad.all,orth.seg.ad.all.ids[f,],my.ge.cod,center = T,scale = T,max.na.prop = 1,ylab.fun.name = 'z-score')
mtext('mouse',3,outer = TRUE)
m=plotExpAndPsiForNewExons('human',c('h','hq','hqmrb','hqmrbo'),orth.seg.ad.all,orth.seg.ad.all.ids[f,],my.ge.cod,center = T,scale = T,max.na.prop = 1,ylab.fun.name = 'z-score')
mtext('human',3,outer = TRUE)

par(mfcol=c(2,5),tck=-0.01,mgp=c(1.1,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,2,1))
m=plotExpAndPsiForNewExons('mouse',c('m','mr','mrb','hqmrb','hqmrbo'),orth.seg.ad.all,orth.seg.ad.all.ids[f,],my.ge.cod,center = F,scale = F,max.na.prop = 1,ylab.fun.name = '')
mtext('mouse',3,outer = TRUE)
m=plotExpAndPsiForNewExons('human',c('h','hq','hqmrb','hqmrbo'),orth.seg.ad.all,orth.seg.ad.all.ids[f,],my.ge.cod,center = F,scale = F,max.na.prop = 1,ylab.fun.name = '')
mtext('human',3,outer = TRUE)
dev.off()
