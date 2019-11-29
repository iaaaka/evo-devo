options(stringsAsFactors = FALSE)
source('code/r.functions/paper.figures.F.R')
# source('code/r.functions/load.all.data.F.R')
# source('code/r.functions/paper.figures.F.R')
# source('code/r.functions/ad.on.ge.F.R')
library(SAJR)
# library(ape)
# library(GO.db)
# library(GenomicRanges)

species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
orth.seg.ad.tsm = readRDS('Rdata/orth.seg.ad.tsm.Rdata')
params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)

##
orth.per.tissue.age.qv = readRDS('Rdata/orth.per.tissue.age.qv.Rdata')
t = do.call(cbind,orth.per.tissue.age.qv)
orth.dev = apply(t<0.05,1,sum,na.rm=T)>0

ancient = apply(sapply(orth.seg.ad,function(x)x$seg$type)=='ALT',1,sum)==7
all.psi = do.call(cbind,lapply(orth.seg.ad,'[[','ir'))
all.psi = all.psi[,rownames(meta)]
#orth.seg.cor = list()

plot(sapply(split(orth.dev,apply(sapply(orth.seg.ad,function(x)x$seg$type)=='ALT',1,sum)),mean))
#aLl/aNcient; Pearson/Spearman; All/Dev,Not-age
# orth.seg.cor$lma = as.matrix(dist(t(all.psi),m='manh'))
# orth.seg.cor$npa = cor(all.psi[ancient,],use = 'pair',method = 'p')
# orth.seg.cor$nsa = cor(all.psi[ancient,],use = 'pair',method = 's')
# orth.seg.cor$lpa = cor(all.psi,          use = 'pair',method = 'p')
# orth.seg.cor$lsa = cor(all.psi,          use = 'pair',method = 's')
# 
# orth.seg.cor$npd = cor(all.psi[ancient & orth.dev,],use = 'pair',method = 'p')
# orth.seg.cor$nsd = cor(all.psi[ancient & orth.dev,],use = 'pair',method = 's')
# orth.seg.cor$lpd = cor(all.psi[orth.dev,          ],use = 'pair',method = 'p')
# orth.seg.cor$lsd = cor(all.psi[orth.dev,          ],use = 'pair',method = 's')
# 
# orth.seg.cor$npn = cor(all.psi[ancient & !orth.dev,],use = 'pair',method = 'p')
# orth.seg.cor$nsn = cor(all.psi[ancient & !orth.dev,],use = 'pair',method = 's')
# orth.seg.cor$lpn = cor(all.psi[!orth.dev,          ],use = 'pair',method = 'p')
# orth.seg.cor$lsn = cor(all.psi[!orth.dev,          ],use = 'pair',method = 's')
#saveRDS(orth.seg.cor,'Rdata/orth.seg.cor.Rdata')

orth.seg.cor = readRDS('Rdata/orth.seg.cor.Rdata')

ts = unique(meta$tissue)
tissues = c(list(ts,ts[3:7],ts[4:7],ts[4:6]),ts)

caclMDSFromCor = function(c,m,tissues,k=2){
	m = m[colnames(c),]
	r = lapply(tissues,function(t){cmdscale(1-c[m$tissue %in% t,m$tissue %in% t],k=k)})
	names(r) = sapply(tissues,function(x){paste(substr(x,1,1),collapse='')})
	r
}

orth.seg.mds = lapply(orth.seg.cor[-9:-10],function(x)caclMDSFromCor(x,meta,tissues,k=2))

page.names=c('ancient, Pearson (1590)','ancient Spearman (1590)','all, Pearson (46210)','all, Spearman (46210)',
						 'ancient, devAS, Pearson (1517)' ,'ancient, devAS, Spearman (1517)' ,'all, devAS, Pearson (5953)' ,'all, devAS, Spearman (5953)',
						 'ancient, !devAS, Pearson (73)','ancient, !devAS, Spearman (73)','all, !devAS, Pearson (40257)','all, !devAS, Spearman (40257)')[-9:-10]

pdf('figures/MDSs/orth.mds.pdf',w=36,h=6)
par(mfcol=c(2,12),tck=-0.01,mgp=c(2,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,1,1))
for(i in 1:length(orth.seg.mds)){
	for(j in 1:length(orth.seg.mds[[i]])){
		x = orth.seg.mds[[i]][[j]]
		m = meta[rownames(x),]
		cs = m$cex
		cs = (cs-min(cs))/(max(cs)-min(cs))*2+0.1
		plot(x,pch=19,cex=cs,col=m$col,xlab='Dimension 1',ylab='Dimension 2',main=paste0('tissues: ',names(orth.seg.mds[[i]])[j]))
		plot(x,pch=19,cex=m$cex,col=params$species.col[m$species],xlab='Dimension 1',ylab='Dimension 2',main=paste0('tissues: ',names(orth.seg.mds[[i]])[j]))
	}
	mtext(page.names[i],3,outer = T)
	plot.new()
	legend('topleft',fill=params$tissue.col,legend=names(params$tissue.col))
	plot.new()
	legend('topleft',fill=params$species.col,legend=names(params$species.col))
}
dev.off()

# 4D MDSs
lpa4=cmdscale(1-orth.seg.cor$lpa,k=4)
npa4=cmdscale(1-orth.seg.cor$npa,k=4)

pdf('figures/MDSs/orth.mds.4D.pdf',w=9,h=6)
par(mfrow=c(2,3),tck=-0.01,mgp=c(2,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,1,1))
m = meta[rownames(lpa4),]
cs = m$cex
cs = (cs-min(cs))/(max(cs)-min(cs))*2+0.1
plot(lpa4[,1:2],pch=19,cex=cs,col=m$col,xlab='Dimension 1',ylab='Dimension 2',main='All Pearson (46210)')
plot(lpa4[,3:4],pch=19,cex=cs,col=m$col,xlab='Dimension 3',ylab='Dimension 3',main='All Pearson (46210)')
plot.new()
legend('topleft',fill=params$tissue.col,legend=names(params$tissue.col))
plot(lpa4[,1:2],pch=19,cex=cs,col=params$species.col[m$species],xlab='Dimension 1',ylab='Dimension 2',main='All Pearson (46210)')
plot(lpa4[,3:4],pch=19,cex=cs,col=params$species.col[m$species],xlab='Dimension 3',ylab='Dimension 3',main='All Pearson (46210)')
plot.new()
legend('topleft',fill=params$species.col,legend=names(params$species.col))


m = meta[rownames(npa4),]
cs = m$cex
cs = (cs-min(cs))/(max(cs)-min(cs))*2+0.1
plot(npa4[,1:2],pch=19,cex=cs,col=m$col,xlab='Dimension 1',ylab='Dimension 2',main='ancient Pearson (1590)')
plot(npa4[,3:4],pch=19,cex=cs,col=m$col,xlab='Dimension 3',ylab='Dimension 3',main='ancient Pearson (1590)')
plot.new()
legend('topleft',fill=params$tissue.col,legend=names(params$tissue.col))
plot(npa4[,1:2],pch=19,cex=cs,col=params$species.col[m$species],xlab='Dimension 1',ylab='Dimension 2',main='ancient Pearson (1590)')
plot(npa4[,3:4],pch=19,cex=cs,col=params$species.col[m$species],xlab='Dimension 3',ylab='Dimension 3',main='ancient Pearson (1590)')
plot.new()
legend('topleft',fill=params$species.col,legend=names(params$species.col))
dev.off()

## MDSs for devAS per tissue
#tissue; ancient/all; pearson/sperman
orth.seg.cor.tis = list()
orth.seg.cor.cnt = c()
for(t in unique(meta$tissue)){
	print(t)
	sgn=apply(sapply(orth.per.tissue.age.qv,function(x)x[,t])<0.05,1,sum,na.rm=T)>0
	psil = all.psi[sgn,meta$tissue==t]
	psin = all.psi[sgn & ancient,meta$tissue==t]
	orth.seg.cor.cnt[paste0(substr(t,1,1),'l')] = sum(sgn)
	orth.seg.cor.cnt[paste0(substr(t,1,1),'n')] = sum(sgn & ancient)
	orth.seg.cor.tis[[paste0(substr(t,1,1),'ls')]] = cor(psil,use='pair',method = 'spearman')
	orth.seg.cor.tis[[paste0(substr(t,1,1),'lp')]] = cor(psil,use='pair',method = 'pearson')
	
	orth.seg.cor.tis[[paste0(substr(t,1,1),'ns')]] = cor(psin,use='pair',method = 'spearman')
	orth.seg.cor.tis[[paste0(substr(t,1,1),'np')]] = cor(psin,use='pair',method = 'pearson')
}

pdf('figures/MDSs/per.tissue.devAS.MDS.pdf',w=12,h=24)
par(mfrow=c(8,4),tck=-0.01,mgp=c(2,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,1,1))
ns = c(np='ancient, Pearson',ns='ancient, Spearman',lp='all, Pearson',ls='all, Spearman')
for(t in unique(meta$tissue)){
	for(i in 1:length(ns)){
		x = cmdscale(1-orth.seg.cor.tis[[paste0(substr(t,1,1),names(ns)[i])]])
		c = orth.seg.cor.cnt[paste0(substr(t,1,1),substr(names(ns)[i],1,1))]
		cs=meta[rownames(x),'cex']
		cs = (cs-min(cs))/(max(cs)-min(cs))*2+0.1
		plot(x,cex=cs,pch=19,col=params$species.col[meta[rownames(x),'species']],xlab='Dimension 1',ylab='Dimension 2',main=paste0(t,', ',ns[i],' (',c,')'))
	}
}
plot.new()
legend('topleft',fill=params$species.col,legend=names(params$species.col))
dev.off()

# conservation for devAS and rest
age.al.i = read.table('input/margarida.age.alignment.final.2017.05.17.improved.tab',sep='\t',header = T)[,c(rownames(species),'to.remove')]




pdf('figures/devAS/species.cor.dev-vs-not.dev.pdf',w=6,h=12)
par(mfrow=c(4,2),tck=-0.01,mgp=c(2,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,1,1))

plotSpCor2Cor(orth.seg.cor$lpd,orth.seg.cor$lpn,age.al.i[,c('mouse','rat')],meta,main='Mouse-rat, Pearson',xlab='devAS (5953)',ylab='not devAS (40257)')
plotSpCor2Cor(orth.seg.cor$lsd,orth.seg.cor$lsn,age.al.i[,c('mouse','rat')],meta,main='Mouse-rat, Spearman',xlab='devAS (5953)',ylab='not devAS (40257)')

plotSpCor2Cor(orth.seg.cor$lpd,orth.seg.cor$lpn,age.al.i[,c('mouse','human')],meta,main='Mouse-human, Pearson',xlab='devAS (5953)',ylab='not devAS (40257)')
plotSpCor2Cor(orth.seg.cor$lsd,orth.seg.cor$lsn,age.al.i[,c('mouse','human')],meta,main='Mouse-human, Spearman',xlab='devAS (5953)',ylab='not devAS (40257)')

plotSpCor2Cor(orth.seg.cor$npd,orth.seg.cor$nsd,age.al.i[,c('mouse','rat')],meta,main='Mouse-rat',xlab='ancient devAS (1517), Pearson',ylab='ancient devAS (1517), Spearman')
plotSpCor2Cor(orth.seg.cor$lpd,orth.seg.cor$lsd,age.al.i[,c('mouse','rat')],meta,main='Mouse-rat',xlab='devAS (5953), Pearson',ylab='devAS (5953), Spearman')

plotSpCor2Cor(orth.seg.cor$lpd,orth.seg.cor$nsd,age.al.i[,c('mouse','rat')],meta,main='Mouse-rat',xlab='devAS (5953), Pearson',ylab='ancient devAS (1517), Spearman')
plotSpCor2Cor(orth.seg.cor$lsd,orth.seg.cor$nsd,age.al.i[,c('mouse','rat')],meta,main='Mouse-rat',xlab='devAS (5953), Spearman',ylab='ancient devAS (1517), Spearman')
dev.off()


#just mouse-rat
p = cbind(orth.seg.ad$mouse$ir,orth.seg.ad$rat$ir)
a = cbind(orth.seg.ad$mouse$seg$type=='ALT',orth.seg.ad$rat$seg$type=='ALT')
d = cbind(apply(orth.per.tissue.age.qv$mouse<0.05,1,sum,na.rm=T)>0,apply(orth.per.tissue.age.qv$rat<0.05,1,sum,na.rm=T)>0)
table(d[,1],d[,2])
table(alt=a[,1] + a[,2],apply(d,1,sum)>0)
#use only both alt
#sum(a[,1] & a[,2] & apply(d,1,sum)>0) 2243
#sum(a[,1] & a[,2] & apply(d,1,sum)==0) 3112
mr.cors = list()
f = a[,1] & a[,2] & apply(is.na(p),1,sum)==0 #& apply(p,1,sd)>0
mr.cors$fpd = cor(p[f & apply(d,1,sum)>0,],use='p',method = 'pear')
mr.cors$fpn = cor(p[f & apply(d,1,sum)==0,],use='p',method = 'pear')
mr.cors$fsd = cor(p[f & apply(d,1,sum)>0,],use='p',method = 'spear')
mr.cors$fsn = cor(p[f & apply(d,1,sum)==0,],use='p',method = 'spear')

mr.cors$npd = cor(p[a[,1] & a[,2] & apply(d,1,sum)>0,],use='p',method = 'pear')
mr.cors$npn = cor(p[a[,1] & a[,2] & apply(d,1,sum)==0,],use='p',method = 'pear')
mr.cors$nsd = cor(p[a[,1] & a[,2] & apply(d,1,sum)>0,],use='p',method = 'spear')
mr.cors$nsn = cor(p[a[,1] & a[,2] & apply(d,1,sum)==0,],use='p',method = 'spear')


m = p[,'MB83M_1']
r = p[,'RB133M_1']
pdf('figures/devAS/mouse-rat.brain.9wpb.cor.dev-vs-not.dev.only.common.pdf',w=9,h=4.5)
par(mfrow=c(1,2),tck=-0.01,mgp=c(1.2,0.2,0),mar=c(2.2,2,1.5,0),oma=c(0,0,1,1))
f = a[,1] & a[,2] & apply(d,1,sum)>0
plot(m[f],r[f],xlab='mouse',ylab='rat',main='Brain, 9wpb, devAS')
legend('topleft',legend=paste0(c('Spearman=','Pearson='),round(c(cor(m[f],r[f],u='p',m='s'),cor(m[f],r[f],u='p',m='p')),3)))
f = a[,1] & a[,2] & apply(d,1,sum)==0
plot(m[f],r[f],xlab='mouse',ylab='rat',main='Brain, 9wpb, non-devAS')
legend('topleft',legend=paste0(c('Spearman=','Pearson='),round(c(cor(m[f],r[f],u='p',m='s'),cor(m[f],r[f],u='p',m='p')),3)))
dev.off()

f = function(p1,p2,aa,thr=0.95){
	ts = unique(meta$tissue)
	r = data.frame(tissue=rep(ts,each=nrow(aa)),stage=rep(aa[,1],times=length(ts)))
	rownames(r) = paste(colnames(aa)[1],r$tissue,r$stage)
	r$pearson=r$spearman=r$seg.cnt=NA
	for(t in ts){
		t1 = paste(colnames(aa)[1],t,aa[,1])
		t2 = paste(colnames(aa)[2],t,aa[,2])
		f = (t1 %in% colnames(p1)) & (t2 %in% colnames(p2))
		t1 = p1[,t1[f]]
		t2 = p2[,t2[f]]
		f = apply(t1 <= thr,1,sum,na.rm=T) > 0 | apply(t2 <= thr,1,sum,na.rm=T) > 0
		t1 = t1[f,]
		t2 = t2[f,]
		for(i in 1:ncol(t1)){
			r[colnames(t1)[i],'pearson']  = cor(t1[,i],t2[,i],u='p',m='p')
			r[colnames(t1)[i],'spearman'] = cor(t1[,i],t2[,i],u='p',m='s')
			r[colnames(t1)[i],'seg.cnt']  = sum(f)
		}
	}
	r
}

#a = cbind(orth.seg.ad$mouse$seg$type=='ALT',orth.seg.ad$rat$seg$type=='ALT')
d = cbind(apply(orth.per.tissue.age.qv$mouse<0.05,1,sum,na.rm=T)>0,apply(orth.per.tissue.age.qv$rat<0.05,1,sum,na.rm=T)>0)

ff = T#a[,1] & a[,2]
dev  = f(orth.seg.ad.tsm$mouse[ff & apply(d,1,sum)> 0,],orth.seg.ad.tsm$rat[ff & apply(d,1,sum)> 0,],age.al.i[,c('mouse','rat')],thr=0.9)
ndev = f(orth.seg.ad.tsm$mouse[ff & apply(d,1,sum)==0,],orth.seg.ad.tsm$rat[ff & apply(d,1,sum)==0,],age.al.i[,c('mouse','rat')],thr=0.9)

pdf('figures/devAS/mouse-rat.cor.dev-vs-not.dev.at.least.one.PSI<0.9.pdf',w=9,h=4.5)
par(mfrow=c(1,2),tck=-0.01,mgp=c(1.2,0.2,0),mar=c(2.2,2,1.5,0),oma=c(0,0,1,1))
plot(dev$pearson,ndev$pearson,pch=19,col=params$tissue.col[dev$tissue],xlab='devAS',ylab='non-devAS',main='Pearson')
plot(dev$spearman,ndev$spearman,pch=19,col=params$tissue.col[dev$tissue],xlab='devAS',ylab='non-devAS',main='Spearman')
dev.off()
# mouse-human
p = cbind(orth.seg.ad$mouse$ir,orth.seg.ad$human$ir)
a = cbind(orth.seg.ad$mouse$seg$type=='ALT',orth.seg.ad$human$seg$type=='ALT')
d = cbind(apply(orth.per.tissue.age.qv$mouse<0.05,1,sum,na.rm=T)>0,apply(orth.per.tissue.age.qv$human<0.05,1,sum,na.rm=T)>0)
table(d[,1],d[,2])
table(alt=a[,1] + a[,2],apply(d,1,sum)>0)
#use only both alt
#sum(a[,1] & a[,2] & apply(d,1,sum)>0) 2039
#sum(a[,1] & a[,2] & apply(d,1,sum)==0) 3132
mh.cors = list()
 f = a[,1] & a[,2] & apply(is.na(p),1,sum)==0 #& apply(p,1,sd)>0
mh.cors$fpd = cor(p[f & apply(d,1,sum)>0,],use='p',method = 'pear')
mh.cors$fpn = cor(p[f & apply(d,1,sum)==0,],use='p',method = 'pear')
mh.cors$fsd = cor(p[f & apply(d,1,sum)>0,],use='p',method = 'spear')
mh.cors$fsn = cor(p[f & apply(d,1,sum)==0,],use='p',method = 'spear')

mh.cors$npd = cor(p[a[,1] & a[,2] & apply(d,1,sum)>0,],use='p',method = 'pear')
mh.cors$npn = cor(p[a[,1] & a[,2] & apply(d,1,sum)==0,],use='p',method = 'pear')
mh.cors$nsd = cor(p[a[,1] & a[,2] & apply(d,1,sum)>0,],use='p',method = 'spear')
mh.cors$nsn = cor(p[a[,1] & a[,2] & apply(d,1,sum)==0,],use='p',method = 'spear')

pdf('figures/devAS/species.cor.dev-vs-not.dev.only.common.pdf',w=7,h=7)
par(mfrow=c(2,2),tck=-0.01,mgp=c(2,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,1,1))
plotSpCor2Cor(mr.cors$npd,mr.cors$npn,age.al.i[,c('mouse','rat')],meta,main='Mouse-rat, ALT in both, Pearson',xlab='devAS (2243)',ylab='not devAS (3112)')
plotSpCor2Cor(mr.cors$nsd,mr.cors$nsn,age.al.i[,c('mouse','rat')],meta,main='Mouse-rat, ALT in both, Spearman',xlab='devAS (2243)',ylab='not devAS (3112)')

#plotSpCor2Cor(mr.cors$fpd,mr.cors$fpn,age.al.i[,c('mouse','rat')],meta,main='Mouse-rat, ALT in both, sd>0, not NA, Pearson',xlab='devAS (287)',ylab='not devAS (610)')
#plotSpCor2Cor(mr.cors$fsd,mr.cors$fsn,age.al.i[,c('mouse','rat')],meta,main='Mouse-rat, ALT in both, sd>0, not NA, Spearman',xlab='devAS (287)',ylab='not devAS (610)')


plotSpCor2Cor(mh.cors$npd,mh.cors$npn,age.al.i[,c('mouse','human')],meta,main='Mouse-human, ALT in both, Pearson',xlab='devAS (2039)',ylab='not devAS (3132)')
plotSpCor2Cor(mh.cors$nsd,mh.cors$nsn,age.al.i[,c('mouse','human')],meta,main='Mouse-human, ALT in both, Spearman',xlab='devAS (2039)',ylab='not devAS (3132)')

#plotSpCor2Cor(mh.cors$fpd,mh.cors$fpn,age.al.i[,c('mouse','human')],meta,main='Mouse-human, ALT in both, sd>0, not NA, Pearson',xlab='devAS (94)',ylab='not devAS (181)')
#plotSpCor2Cor(mh.cors$fsd,mh.cors$fsn,age.al.i[,c('mouse','human')],meta,main='Mouse-human, ALT in both, sd>0, not NA, Spearman',xlab='devAS (94)',ylab='not devAS (181)')
dev.off()

# 
# orth.seg.ad.tsm=readRDS('Rdata/orth.seg.ad.tsm.Rdata')
# 
# cor(orth.seg.ad.tsm$mouse[,'mouse brain 9wpb'],orth.seg.ad.tsm$rat[,'rat brain 16wpb'],u='p',m='s')
# f = orth.seg.ad.tsm$mouse[,'mouse brain 9wpb'] <0.9 & orth.seg.ad.tsm$rat[,'rat brain 16wpb'] <0.9
# f = orth.seg.ad.tsm$mouse[,'mouse brain 9wpb'] <0.9 & orth.seg.ad.tsm$rat[,'rat brain 16wpb'] <0.9 & orth.seg.ad.tsm$mouse[,'mouse kidney 9wpb'] <0.9 & orth.seg.ad.tsm$rat[,'rat kidney 16wpb'] <0.9
# table(f)
# cor(orth.seg.ad.tsm$mouse[f,'mouse brain 9wpb'],orth.seg.ad.tsm$rat[f,'rat brain 16wpb'],u='p',m='p')
# plot(orth.seg.ad.tsm$mouse[f,'mouse brain 9wpb'],orth.seg.ad.tsm$rat[f,'rat brain 16wpb'],pch=19)
# f = !is.na(f) & f
# plot(rank(orth.seg.ad.tsm$mouse[f,'mouse brain 9wpb']),rank(orth.seg.ad.tsm$rat[f,'rat brain 16wpb']),pch=19)
# 
# cor(orth.seg.ad.tsm$mouse[,'mouse kidney 9wpb'],orth.seg.ad.tsm$rat[,'rat kidney 16wpb'],u='p',m='s')
# f = orth.seg.ad.tsm$mouse[,'mouse kidney 9wpb'] <0.9 & orth.seg.ad.tsm$rat[,'rat kidney 16wpb'] <0.9
# table(f)
# cor(orth.seg.ad.tsm$mouse[f,'mouse kidney 9wpb'],orth.seg.ad.tsm$rat[f,'rat kidney 16wpb'],u='p',m='p')
# 
# 
# f = !is.na(f) & f
# plot((orth.seg.ad.tsm$mouse[f,'mouse kidney 9wpb']),(orth.seg.ad.tsm$rat[f,'rat kidney 16wpb']),pch=19)
# plot(rank(orth.seg.ad.tsm$mouse[f,'mouse kidney 9wpb']),rank(orth.seg.ad.tsm$rat[f,'rat kidney 16wpb']),pch=19)
# 
x = cbind(mb=orth.seg.ad.tsm$mouse[,'mouse brain 2wpb'],rb=orth.seg.ad.tsm$rat[,'rat brain 2wpb'],
					ml=orth.seg.ad.tsm$mouse[,'mouse liver 2wpb'],rl=orth.seg.ad.tsm$rat[,'rat liver 2wpb'])

f = orth.dev
f = apply(is.na(x),1,sum)==0
#f = f & apply(x,1,sd)>0

cr = function(x,thr,m){
	f1 = apply(x[,1:2]>=thr & x[,1:2]<=(1-thr),1,sum)>0
	f2 = apply(x[,3:4]>=thr & x[,3:4]<=(1-thr),1,sum)>0
	c(Nb=sum(f1,na.rm=T),Nl=sum(f2,na.rm=T),b=cor(x[f1,1],x[f1,2],m=m,u='p'),l=cor(x[f2,3],x[f2,4],m=m,u='p'))
}
thrs=seq(0,0.3,length.out = 33)
pp = sapply(thrs,function(t)cr(x[ancient,],t,m='p'))
plot(thrs,pp[3,],t='l',col=params$tissue.col['brain'],ylim=range(pp[3:4,]))
lines(thrs,pp[4,],col=params$tissue.col['liver'])

inx=1:2
b = x[apply(is.na(x[,inx]),1,sum)==0 & apply(x[,inx]>0.05 & x[,inx]<0.95,1,sum)>0,inx]
inx=3:4
l = x[apply(is.na(x[,inx]),1,sum)==0 & apply(x[,inx]>0.05 & x[,inx]<0.95,1,sum)>0,inx]


z=plotSpCor2Cor(orth.seg.cor$lpa,orth.seg.cor$lma,age.al.i[,c('mouse','rat')],meta,main='Mouse-human, ALT in both, Spearman',xlab='devAS (2039)',ylab='not devAS (3132)')
x=z[z[,2]=='brain',]
x[order(-x[,1])[1:4],]

c=orth.seg.cor$lsa
c['RK133M_2','MK83M_2',drop=F]
c['RB35M_1','MB34F_1',drop=F]

x = cbind(orth.seg.ad$mouse$ir[,c('MB34F_1','MK83M_2')],orth.seg.ad$rat$ir[,c('RB35M_1','RK133M_2')])
x = x[,c(1,3,2,4)]
cor(x,u='p',m='s')


f = function(r1,r2=0.1,n1=1000,n2=9000){
	x = runif(n1,0,1)
	y = runif(n1,pmax(0,x-r1),pmin(x+r1,1))
	x = c(x,runif(n2,1-r2,1))
	y = c(y,runif(n2,1-r2,1))
	cbind(x,y)
}

r=t(sapply(1:1000,function(i){
	cat('\r',i)
	d1 = f(0.3)
	d2 = f(0.6,n1=1500,n2=8500)
	c(cor(d1[,1],d1[,2],m='p'),cor(d2[,1],d2[,2],m='p'),cor(d1[,1],d1[,2],m='s'),cor(d2[,1],d2[,2],m='s'))
}))

boxplot(r)
stat = round(apply(r,2,quantile,prob=c(0.025,0.5,0.975)),3)

d1 = f(0.3)
d2 = f(0.6,n1=1500,n2=8500)
pdf('figures/pearson.vs.spearman.modelled.pdf',w=9,h=9)
par(mfrow=c(2,2),tck=-0.01,mgp=c(2,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,1,1))

plot(d1,main='Modeled "kidney"',xlab='species 1',ylab='species 2',col=c(rep('red',1000),rep('gray',9000)))
legend('topleft',legend=c(paste0(c('pearson=','spearman='),stat[2,c(1,3)],' [',stat[1,c(1,3)],',',stat[3,c(1,3)],']')))
legend('bottomright',col=c('gray','red'),pch=1,legend=c('Part 1 (close to 1; 9000)','Part 2 (variable; 1000)'))

plot(d2,main='Modeled "brain"',xlab='species 1',ylab='species 2',col=c(rep('red',1500),rep('gray',8500)))
legend('topleft',legend=c(paste0(c('pearson=','spearman='),stat[2,c(2,4)],' [',stat[1,c(2,4)],',',stat[3,c(2,4)],']')))
legend('bottomright',col=c('gray','red'),pch=1,legend=c('Part 1 (close to 1; 8500)','Part 2 (variable; 1500)'))
x = cbind(orth.seg.ad$mouse$ir[,c('MB83M_1','MK83M_1')],orth.seg.ad$rat$ir[,c('RB133M_1','RK133M_2')])
x = x[,c(1,3,2,4)]

f(x[,3:4],main='Kidney 9wpb',xlab='mouse',ylab='rat')
f(x[,1:2],main='Brain 9wpb',xlab='mouse',ylab='rat')
dev.off()

f = function(x,thr=0.9,...){
	x = x[apply(is.na(x),1,sum)==0,]
	ft = x[,1] > thr & x[,2] > thr
	plot(x[!ft,],col='red',...)
	points(x[ft,],col='gray')
	legend('topleft',legend=paste0(c('pearson=','spearman='),round(c(cor(x[,1],x[,2],m='p'),cor(x[,1],x[,2],m='s')),3)))
	legend('bottomright',col=c('gray','red'),pch=1,legend=paste0(c('Part 1 (close to 1; ','Part 2 (variable; '),c(sum(ft),sum(!ft)),')'))
}

cor(x,u='p',m='s')

### proportion of devAS on evolutionary age
# alternification
alt = sapply(orth.seg.ad,function(x)x$seg$type=='ALT')
devAS = sapply(orth.per.tissue.age.qv,function(x)apply(x<0.05,1,sum,na.rm=T)>0)
alt.sp = apply(alt,1,function(x)paste(species$short[x],collapse=''))

sp = c(species$short,'hq','mr','mrb','hqmrb','hqmrbo','hqmrboc')

t = split(apply(devAS,1,sum),alt.sp)[sp]
barplot(sapply(t,function(x)mean(x>0)),las=3)
barplot(sapply(t,mean)/nchar(names(t)),las=3)
