#setwd('~/skoltech/projects/evo.devo/')

source('code/r.functions/load.all.data.F.R')
source('code/r.functions/paper.figures.F.R')
library(SAJR)
library(GenomicRanges)
library(doMC)

species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
anns = readRDS('Rdata/anns.Rdata')
# all.anns = readRDS('Rdata/all.anns.Rdata') 
psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
psi.tsm.ms = readRDS('Rdata/psi.tsm.ms.Rdata')
border.stages = readRDS('Rdata/border.stages.Rdata')
seg2ens = readRDS('Rdata/seg2ens.Rdata')
ens.ge.cod = readRDS('Rdata/ens.ge.cod.Rdata')
# orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
# orth.seg.ad.tsm = readRDS('Rdata/orth.seg.ad.tsm.Rdata')
born.exn.sajr = readRDS('Rdata/born.exn.sajr.Rdata')
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')
age.dpsi = readRDS('Rdata/age.diam.spline4.with.replicates.Rdata')

ens.descr = unique(read.table('input/hs.37.74.gene.descr.txt',sep=',',quote='"',header=TRUE)[,1:3])
ens.descr$Description = sapply(strsplit(ens.descr$Description,' [',TRUE),'[',1)
rownames(ens.descr) = ens.descr$Ensembl.Gene.ID


age.al.i = read.table('input/margarida.age.alignment.final.2017.05.17.improved.tab',sep='\t',header = T)[,c(rownames(species),'to.remove')]
params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)

# devAS ####
# orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
# orth.per.tissue.age.qv = vector('list',nrow(species))
# names(orth.per.tissue.age.qv) = rownames(species)
# 
# registerDoMC(16)
# for(s in rownames(species)){
# 	cat(toupper(s))
# 	tmp = orth.seg.ad[[s]]
# 	if(s == 'human')
# 		tmp = tmp[,!(colnames(tmp$ir) %in% rownames(meta)[meta$stage %in% c('oldermidage','senior')])]
# 	orth.per.tissue.age.qv[[s]] = sapply(unique(meta$tissue),function(t){testASAge(tmp,meta,t,min.cov.sams=0.6)})
# 	colnames(orth.per.tissue.age.qv[[s]]) = unique(meta$tissue)
# 	dimnames(orth.per.tissue.age.qv[[s]]) = setNames(dimnames(orth.per.tissue.age.qv[[s]]),NULL)
# }
# saveRDS(orth.per.tissue.age.qv,'Rdata/orth.per.tissue.age.qv.Rdata')

# born.per.tissue.age.qv = vector('list',nrow(species))
# names(born.per.tissue.age.qv) = rownames(species)
# 
# for(s in rownames(species)){
# 	print(s)
# 	tmp = born.exn.sajr[[s]]
# 	if(s == 'human')
# 		tmp = tmp[,!(colnames(tmp$ir) %in% rownames(meta)[meta$stage %in% c('oldermidage','senior')])]
# 	born.per.tissue.age.qv[[s]] = sapply(unique(meta$tissue),function(t){testASAge(tmp,meta,t,min.cov.sams=0.6)})
# 	colnames(born.per.tissue.age.qv[[s]]) = unique(meta$tissue)
# 	dimnames(born.per.tissue.age.qv[[s]]) = setNames(dimnames(born.per.tissue.age.qv[[s]]),NULL)
# }
# saveRDS(born.per.tissue.age.qv,'Rdata/born.per.tissue.age.qv.Rdata')
# 
# per.tissue.age.qv = vector('list',nrow(species))
# names(per.tissue.age.qv) = rownames(species)
# registerDoMC(16)
# for(s in rownames(species)){
# 	cat(toupper(s))
# 	tmp = readRDS(paste('Rdata/',s,'.as.u.filtered.Rdata',sep=''))
# 	if(s == 'human')
# 		tmp = tmp[,!(colnames(tmp$ir) %in% rownames(meta)[meta$stage %in% c('oldermidage','senior')])] # to remove aging samples
# 	per.tissue.age.qv[[s]] = sapply(unique(meta$tissue),function(t){testASAge(tmp,meta,t,min.cov.sams=0.6)})
# 	colnames(per.tissue.age.qv[[s]]) = unique(meta$tissue)
# 	dimnames(per.tissue.age.qv[[s]]) = setNames(dimnames(per.tissue.age.qv[[s]]),NULL)
# }
# saveRDS(per.tissue.age.qv,'Rdata/per.tissue.age.qv.Rdata')
# 
# per.tissue.age.amp = per.tissue.age.qv
# for(s in rownames(species)){
# 	cat(toupper(s))
# 	tmp = psi.tsm[[s]]
# 	if(s == 'human')
# 		tmp = tmp[,!(colnames(tmp) %in% rownames(meta.tsm)[meta.tsm$stage %in% c('oldermidage','senior')])]
# 	for(t in unique(meta$tissue)){
# 		tt = tmp[,colnames(tmp) %in% rownames(meta.tsm)[meta.tsm$tissue == t]]
# 		per.tissue.age.amp[[s]][,t] = apply(tt,1,max,na.rm=T) - apply(tt,1,min,na.rm=T)
# 	}
# }
# saveRDS(per.tissue.age.amp,'Rdata/per.tissue.age.amp.Rdata')

# devAS before sex maturation ####
# check sex maturation

# oldest "non-mature" sample
my.sex.maturation = setNames(c(1740,530,23,24,44,43,56),rownames(species)) #according to Margarida suggestions

pdf('figures/devAS/befor.sex.maturation/testis.MDS1.on.age.pdf',w=9,h=9)
par(mfrow=c(3,3),tck=-0.01,mgp=c(2,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,1,1))
for(s in rownames(species)){
	m = meta.tsm[meta.tsm$species==s & meta.tsm$tissue=='testis',]
	m = m[order(m$days),]
	p = psi.tsm[[s]][anns[[s]]$sites=='ad',rownames(m)]
	mds = cmdscale(1-cor(p,u='p'),1)
	plot(mds,main=s,xaxt='n',pch=19,col=(m$days>my.sex.maturation[s])+1,ylab='MDS1',xlab='days')
	abline(v=1:nrow(m),col='gray',lty=3)
	axis(1,1:nrow(m),round(m$days),las=3)
}
dev.off()

# per.tissue.age.bm.qv = vector('list',nrow(species))
# names(per.tissue.age.bm.qv) = rownames(species)
# registerDoMC(16)
# for(s in rownames(species)){
# 	cat(toupper(s))
# 	tmp = readRDS(paste('Rdata/',s,'.as.u.filtered.Rdata',sep=''))
# 	if(s == 'human')
# 		tmp = tmp[,!(colnames(tmp$ir) %in% rownames(meta)[meta$stage %in% c('oldermidage','senior')])]
# 	tmp = tmp[,(colnames(tmp$ir) %in% rownames(meta)[meta$days <= my.sex.maturation[s]])]
# 	per.tissue.age.bm.qv[[s]] = sapply(unique(meta$tissue),function(t){testASAge(tmp,meta,t,min.cov.sams=0.6)})
# 	colnames(per.tissue.age.bm.qv[[s]]) = unique(meta$tissue)
# 	dimnames(per.tissue.age.bm.qv[[s]]) = setNames(dimnames(per.tissue.age.bm.qv[[s]]),NULL)
# }
#saveRDS(per.tissue.age.bm.qv,'Rdata/per.tissue.age.bm.qv.Rdata')

# orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
# orth.per.tissue.age.bm.qv = vector('list',nrow(species))
# names(orth.per.tissue.age.bm.qv) = rownames(species)
# 
# registerDoMC(16)
# for(s in rownames(species)){
# 	cat(toupper(s))
# 	tmp = orth.seg.ad[[s]]
# 	if(s == 'human')
# 		tmp = tmp[,!(colnames(tmp$ir) %in% rownames(meta)[meta$stage %in% c('oldermidage','senior')])]
# 	tmp = tmp[,(colnames(tmp$ir) %in% rownames(meta)[meta$days < my.sex.maturation[s]])]
# 	orth.per.tissue.age.bm.qv[[s]] = sapply(unique(meta$tissue),function(t){testASAge(tmp,meta,t,min.cov.sams=0.6)})
# 	colnames(orth.per.tissue.age.bm.qv[[s]]) = unique(meta$tissue)
# 	dimnames(orth.per.tissue.age.bm.qv[[s]]) = setNames(dimnames(orth.per.tissue.age.bm.qv[[s]]),NULL)
# }
# saveRDS(orth.per.tissue.age.bm.qv,'Rdata/orth.per.tissue.age.bm.qv.Rdata')


per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')
per.tissue.age.amp = readRDS('Rdata/per.tissue.age.amp.Rdata')


f = function(a){a$sites=='ad'}
# plot amp distr
pdf('figures/devAS/ad.amp.pdf',w=21,h=21)
par(mfrow=c(7,7),tck=-0.01,mgp=c(2,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,1,1))
brakes = 0:50/50
for(s in names(per.tissue.age.qv)){
	ff = f(anns[[s]])	
	for(t in colnames(per.tissue.age.qv[[1]])){
		fff = ff & !is.na(per.tissue.age.qv[[s]][,t])
		if(sum(fff)==0){
			plot.new()
		}
		else{
			sg = per.tissue.age.amp[[s]][fff & per.tissue.age.qv[[s]][,t]<=0.05,t]
			ns = per.tissue.age.amp[[s]][fff & per.tissue.age.qv[[s]][,t]> 0.05,t]
			ymax = max(hist(sg,brakes,plot=FALSE)$density,hist(ns,brakes,plot=FALSE)$density)
			hist(ns,brakes,freq = FALSE,ylim=c(0,ymax),col='gray',border = NA,main=paste0(s,'-',t,'(',length(sg),'/',length(ns),')'),xlab='max(PSI)-min(PSI)')
			hist(sg,brakes,freq = FALSE,col='#FF000088',border = NA,xlab='',ylab='',main='',add=T)
			abline(v=c(0.2,0.5),lty=2)
		}
	}
}
dev.off()
# plot stat
# orth.per.tissue.age.qv = readRDS('Rdata/orth.per.tissue.age.qv.Rdata')
# 
# orth.tst.cnt   = sapply(rownames(species),function(s){apply(!is.na(orth.per.tissue.age.qv[[s]]),2,sum)})
# orth.sgn.cnt   = sapply(rownames(species),function(s){apply(orth.per.tissue.age.qv[[s]]<0.05,2,sum,na.rm=T)})

tst.cnt   = sapply(rownames(species),function(s){apply(!is.na(per.tissue.age.qv[[s]][f(anns[[s]]),]),2,sum)})
sgn.cnt   = sapply(rownames(species),function(s){apply(per.tissue.age.qv[[s]][f(anns[[s]]),]<0.05,2,sum,na.rm=T)})
sgn.cnt02 = sapply(rownames(species),function(s){apply(per.tissue.age.qv[[s]][f(anns[[s]]),]<0.05 & per.tissue.age.amp[[s]][f(anns[[s]]),]>0.2,2,sum,na.rm=T)})
sgn.cnt05 = sapply(rownames(species),function(s){apply(per.tissue.age.qv[[s]][f(anns[[s]]),]<0.05 & per.tissue.age.amp[[s]][f(anns[[s]]),]>0.5,2,sum,na.rm=T)})
ad2em.05  = sapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,psi.thr = 0.5,border.stages,s))
ad2em.05.sgn.pr = sapply(rownames(species),function(s)sapply(unique(meta$tissue),function(t){
	if(!(t %in% colnames(ad2em.05[[s]])))
		 return(NA)
	a2e = ad2em.05[[s]][,t] %in% c('u','d')
	qv = per.tissue.age.qv[[s]][a2e & f(anns[[s]]),t]
	sum(qv < 0.05,na.rm=T)/length(qv)
	}))


pdf('figures/devAS/ad.sgn.stat.pdf',w=20,h=12)
par(mfrow=c(3,4),tck=-0.01,mgp=c(2,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,1,1))
cs=rep(params$tissue.col,each=7)
b=barplot(t(tst.cnt),beside = T,col=cs,ylab='# of cassette exons',main='Tested')
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.9)
barplot(t(sgn.cnt),beside = T,col=cs,ylab='# of cassette exons',main='FDR < 0.05')
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.9)
barplot(t(sgn.cnt02),beside = T,col=cs,ylab='# of cassette exons',main='FDR < 0.05 & dPSI > 0.2')
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.9)
barplot(t(sgn.cnt05),beside = T,col=cs,ylab='# of cassette exons',main='FDR < 0.05 & dPSI > 0.5')
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.9)
barplot(tst.cnt,beside = T,col=params$tissue.col,ylab='# of cassette exons',main='Tested')

barplot(t(sgn.cnt/tst.cnt*100),beside = T,col=rep(params$tissue.col,each=7),ylab='% sign./tested',main='FDR < 0.05')
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.9)
barplot(t(sgn.cnt02/tst.cnt*100),beside = T,col=rep(params$tissue.col,each=7),ylab='% sign./tested',main='FDR < 0.05 & dPSI > 0.2')
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.9)
barplot(t(sgn.cnt05/tst.cnt*100),beside = T,col=rep(params$tissue.col,each=7),ylab='% sign./tested',main='FDR < 0.05 & dPSI > 0.5')
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.9)

barplot(t(ad2em.05.sgn.pr*100),beside = T,col=rep(params$tissue.col,each=7),ylab='% of significant',main='|adult(PSI)-embryo(PSI)| > 0.5')
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.9)

# barplot(t(orth.sgn.cnt/orth.tst.cnt*100),beside = T,col=rep(params$tissue.col,each=7),ylab='% of significant',main='|adult(PSI)-embryo(PSI)| > 0.5')
# barplot(t(orth.sgn.cnt),beside = T,col=rep(params$tissue.col,each=7),ylab='% of significant',main='|adult(PSI)-embryo(PSI)| > 0.5')
# text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.9)
dev.off()

# sex-maturation sgn cnts ####
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')
per.tissue.age.bm.qv = readRDS('Rdata/per.tissue.age.bm.qv.Rdata')

age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)
age.dpsi$macaque = cbind(age.dpsi$macaque[,1:5],ovary=NaN,age.dpsi$macaque[,6,drop=FALSE])

border.stages.bm = border.stages


for(s in names(border.stages.bm)){
	for(t in rownames(border.stages.bm[[s]])){
		m = meta[meta$species == s & meta$tissue == t & meta$days <= my.sex.maturation[s],]
		if(nrow(m) > 0){
			border.stages.bm[[s]][t,2] = m$stage[order(m$day,decreasing = T)[1]]
		}
	}
}

age.dpsi.bm = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages.bm,s,get.dPSI=T))
names(age.dpsi.bm) = rownames(species)
age.dpsi.bm$macaque = cbind(age.dpsi$macaque[,1:5],ovary=NaN,age.dpsi$macaque[,6,drop=FALSE])

sgn02.stat = sapply(names(anns),function(s)	apply(per.tissue.age.qv[[s]][anns[[s]]$sites == 'ad',] < 0.05 & abs(age.dpsi[[s]][anns[[s]]$sites == 'ad',])>0.2,2,sum,na.rm=T))
sgn02.stat.bm = sapply(names(anns),function(s)	apply(per.tissue.age.bm.qv[[s]][anns[[s]]$sites == 'ad',] < 0.05 & abs(age.dpsi.bm[[s]][anns[[s]]$sites == 'ad',])>0.2,2,sum,na.rm=T))

sgn05.stat = sapply(names(anns),function(s)	apply(per.tissue.age.qv[[s]][anns[[s]]$sites == 'ad',] < 0.05 & abs(age.dpsi[[s]][anns[[s]]$sites == 'ad',])>0.5,2,sum,na.rm=T))
sgn05.stat.bm = sapply(names(anns),function(s)	apply(per.tissue.age.bm.qv[[s]][anns[[s]]$sites == 'ad',] < 0.05 & abs(age.dpsi.bm[[s]][anns[[s]]$sites == 'ad',])>0.5,2,sum,na.rm=T))


pdf('figures/devAS/befor.sex.maturation/ad.sgn.stat.bm.pdf',w=10,h=4)
#jpeg('figures/devAS/befor.sex.maturation/ad.sgn.stat.bm.jpg',units = 'in',res = 600,w=10,h=4)
par(mfrow=c(1,2),tck=-0.01,mgp=c(2,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,1,1))
cs=rep(params$tissue.col,each=7)
b=barplot(t(sgn02.stat),beside = T,col=paste0(cs,'60'),ylab='# of cassette exons',main='FDR < 0.05 & |dPSI| > 0.2',border=NA,cex.names = 0.7)
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.7)
b=barplot(t(sgn02.stat.bm),beside = T,col=cs,ylab='# of cassette exons',main='FDR < 0.05; Before sexual maturation',add=T,border=NA,xaxt='n')

b=barplot(t(sgn05.stat),beside = T,col=paste0(cs,'60'),ylab='# of cassette exons',main='FDR < 0.05 & |dPSI| > 0.5',border=NA,cex.names=0.7)
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.7)
b=barplot(t(sgn05.stat.bm),beside = T,col=cs,ylab='# of cassette exons',main='FDR < 0.05; Before sexual maturation',add=T,border=NA,xaxt='n')

legend('topright',fill=c('#00000060','#000000'),legend=c('Whole lifespan','Before sexual maturation'),bty = 'n')
dev.off()

# paper S4 ######
pdf('figures/paper.figures/5/suppl/S4.pdf',w=8,h=4)
#jpeg('figures/paper.figures/5/suppl/S4.jpg',units = 'in',res = 600,w=10,h=4)
par(mfrow=c(1,1),tck=-0.01,mgp=c(1.5,0.5,0),mar=c(2.5,2.5,.5,0),oma=c(0,0,0,1))
cs=rep(params$tissue.col,each=7)
b=barplot(t(sgn02.stat),beside = T,col=paste0(cs,'60'),ylab='# of cassette exons',main='',border=NA,cex.names = 1)
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=1)
b=barplot(t(sgn02.stat.bm),beside = T,col=cs,ylab='# of cassette exons',main='FDR < 0.05; Before sexual maturation',add=T,border=NA,xaxt='n')
legend('topright',fill=c('#00000060','#000000'),legend=c('Whole lifespan','Before sexual maturation'),bty = 'n')
dev.off()

#compare to Proteomic analysis of postsynaptic proteins in regions of the human neocortex ####
egid2ens = unique(read.table('input/hs.37.74.gene.descr.txt',sep=',',quote='"',header=TRUE)[,c(1,6)])
egid2ens = egid2ens[!is.na(egid2ens$RefSeq.Protein.ID..e.g..NP_001005353.) & egid2ens$RefSeq.Protein.ID..e.g..NP_001005353. != '',]

psp = read.csv('input/41593_2017_25_MOESM4_ESM.PMID: 29203896.csv')
psp$entrez.id = sapply(strsplit(psp$Accessiona,'[|.]'),'[',4)

table(psp$entrez.id %in% egid2ens$RefSeq.Protein.ID..e.g..NP_001005353.)
psp$Accessiona[!(psp$entrez.id %in% egid2ens$RefSeq.Protein.ID..e.g..NP_001005353.)]

#psp = psp[psp$entrez.id %in% egid2ens$RefSeq.Protein.ID..e.g..NP_001005353.,]
egid2ens = egid2ens[egid2ens$RefSeq.Protein.ID..e.g..NP_001005353. %in% psp$entrez.id,]
psp$ens.gid=setNames(egid2ens$Ensembl.Gene.ID,egid2ens$RefSeq.Protein.ID..e.g..NP_001005353.)[psp$entrez.id]

hist(psp$Anova..p.e)
hist(log(psp$X141124_o2_04_so141117_hu_48))
dPSI=0.3
h.age.ens = getAgeASEns(psi.tsm,meta.tsm,dPSI,border.stages,'human')
h.age.ens.me = getAgeASEns(psi.tsm,meta.tsm,dPSI,border.stages,'human',max.len = 27)
fisher.test(table(h.age.ens$brain$u %in% psp$ens.gid,h.age.ens$brain$u %in% h.age.ens.me$brain$u))



as.in.psp=lapply(h.age.ens,function(x){
	tested = unique(unlist(x))
	r = list()
	r$psp2u = fisher.test(table(tested %in% psp$ens.gid,tested %in% x$u))
	r$psp2d = fisher.test(table(tested %in% psp$ens.gid,tested %in% x$d))
	r$psp2b = fisher.test(table(tested %in% psp$ens.gid,tested %in% c(x$d,x$u)))
	
	f = psp$ens.gid %in% tested
	u = wilcox.test(psp$Anova..p.e[psp$ens.gid %in% x$u],psp$Anova..p.e[f & !(psp$ens.gid %in% x$u)],a='l')$p.value
	d = wilcox.test(psp$Anova..p.e[psp$ens.gid %in% x$d],psp$Anova..p.e[f & !(psp$ens.gid %in% x$d)],a='l')$p.value
	b = wilcox.test(psp$Anova..p.e[psp$ens.gid %in% c(x$d,x$u)],psp$Anova..p.e[f & !(psp$ens.gid %in% c(x$u,x$d))],a='l')$p.value
	r$as2pspDE = c(up=u,down=d,both=b)
	r
})

pdf('figures/postsynaptic.proteins.ageAS.dPSI=0.5.pdf',w=12,h=6)
layout(matrix(c(1,2,3,3),ncol=2),widths = c(1,1))
par(tck=-0.01,mgp=c(1.9,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,1,1))
or = log2(sapply(as.in.psp,function(x){c(up=x$psp2u$estimate,down=x$psp2d$estimate,both=x$psp2b$estimate)}))
ci1 = log2(sapply(as.in.psp,function(x){c(up=x$psp2u$conf.int[1],down=x$psp2d$conf.int[1],both=x$psp2b$conf.int[1])}))
ci2 = log2(sapply(as.in.psp,function(x){c(up=x$psp2u$conf.int[2],down=x$psp2d$conf.int[2],both=x$psp2b$conf.int[2])}))
b=barplot(or,beside = T,legend.text = c('up','down','up or down'),ylab='log2(odds ratio)',ylim=range(or,ci1,ci2),main='Assotiation of ageAS with PSP')
segments(b,ci1,b,ci2)

x = sapply(as.in.psp,function(x){x$as2pspDE})
barplot(-log10(x),beside = T,legend.text = c('up','down','up or down'),ylab='-log10(p-value)',main='Assotiation of ageAS with diff. PSP')
abline(h=-log10(c(0.05,max(x[p.adjust(x,m='BH')<0.05]))),col='red')
as = unique(c(h.age.ens$brain$u,h.age.ens$brain$d))
f = psp$ens.gid %in% c(h.age.ens$brain$n,as)
h = t(apply(table(psp$Highestg[f],psp$ens.gid[f] %in% as)[,2:1],1,my.binom.test))
l = t(apply(table(psp$Lowesth[f] ,psp$ens.gid[f] %in% as)[,2:1],1,my.binom.test))
plot(1,t='n',xlim=range(h),ylim=range(l),xlab='ageAS prop in highest regions',ylab='ageAS prop in lowest regions',main='ageAS in different brain regions')
segments(h[,2],l[,1],h[,3],l[,1],col='gray')
segments(h[,1],l[,2],h[,1],l[,3],col='gray')
text(h[,1],l[,1],rownames(h))

dev.off()





a = unique(unlist(h.age.ens$heart))
fisher.test(table(a %in% psp$ens.gid,a %in% h.age.ens$heart$d))

# ohnologs#####
sapply(h.age.ens,function(x)sapply(x,length))
ohnologs = list()
ohnologs$m = read.table('input/ohnologs/MOUSE.Pairs.Intermediate.2R.txt',sep='\t',header = T)
ohnologs$h = read.table('input/ohnologs/HUMAN.Pairs.Intermediate.2R.txt',sep='\t',header = T)
oh = union(mo$Ohnolog.1.Id,mo$Ohnolog.2.Id)



#Haploinsufficiency
hapin1=read.table('input/ohnologs/Characterising and Predicting Haploinsufficiency in the Human Genome/Dataset_S1.txt',sep='\t',skip = 1)
hapin2=read.table('input/ohnologs/Characterising and Predicting Haploinsufficiency in the Human Genome/Dataset_S2.txt',sep='\t',skip = 1)
hapin1$gene.name = sapply(strsplit(hapin1$V4,'|',fixed = T),'[',1)
hapin2$gene.name = sapply(strsplit(hapin2$V4,'|',fixed = T),'[',1)

table(hapin1$gene.name %in% ens.descr$Associated.Gene.Name)
table(hapin2$gene.name %in% ens.descr$Associated.Gene.Name)

hapin1 = hapin1[hapin1$gene.name %in% ens.descr$Associated.Gene.Name,]
hapin2 = hapin2[hapin2$gene.name %in% ens.descr$Associated.Gene.Name,]

n2e = setNames(rownames(ens.descr),ens.descr$Associated.Gene.Name)
rownames(hapin1) = n2e[hapin1$gene.name]
rownames(hapin2) = n2e[hapin2$gene.name]

hi1.by.as = lapply(h.age.ens,function(x)lapply(x,function(z){
	hapin1[intersect(z,rownames(hapin1)),5]
	}))

hi2.by.as = lapply(h.age.ens,function(x)lapply(x,function(z){
	hapin2[intersect(z,rownames(hapin2)),5]
}))

sapply(hi1.by.as,function(x)sapply(x,length))


sapply(m.age.ens,function(x)sapply(x,length))

xxx = lapply(m.age.ens,function(x)lapply(x,function(z)intersect(z,rownames(ens.ge.cod$mouse$gene))))
t = lapply(xxx,function(x)lapply(x,function(z)setdiff(z,c(ohnologs$m$Ohnolog.1.Id,ohnologs$m$Ohnolog.2.Id))))
noto = sapply(t,function(x)sapply(x,length))
t = lapply(xxx,function(x)lapply(x,function(z)intersect(z,c(ohnologs$m$Ohnolog.1.Id,ohnologs$m$Ohnolog.2.Id))))
ino = sapply(t,function(x)sapply(x,length))

plot(apply(noto[-1,,drop=F],2,sum)/apply(noto,2,sum),apply(ino[-1,,drop=F],2,sum)/apply(ino,2,sum),pch=19,col=params$tissue.col,xlab='freq. of ageAS in not ohnologs',ylab='freq. of ageAS in ohnologs')
abline(a=0,b=1)
dev.off()
f = function(dPSI){
	h.age.ens = getAgeASEns(psi.tsm,meta.tsm,dPSI,border.stages,'human')
	m.age.ens = getAgeASEns(psi.tsm,meta.tsm,dPSI,border.stages,'mouse')
	
	plotOhnologsFreq(m.age.ens,c(ohnologs$m$Ohnolog.1.Id,ohnologs$m$Ohnolog.2.Id),main='Mouse',filter = rownames(ens.ge.cod$mouse$gene))
	plotOhnologsFreq(h.age.ens,c(ohnologs$h$Ohnolog.1.Id,ohnologs$h$Ohnolog.2.Id),main='Human',filter = rownames(ens.ge.cod$human$gene))
	
	plotOhnologsFreq(h.age.ens,rownames(hapin1)[hapin1$V5>0.8],main='Human dataset1',filter = rownames(hapin1),ylab='proportion of genes with p(Haploinsufficiency)>0.8')
	plotOhnologsFreq(h.age.ens,rownames(hapin2)[hapin2$V5>0.8],main='Human dataset2',filter = rownames(hapin2),ylab='proportion of genes with p(Haploinsufficiency)>0.8',legend = T)
	mtext(paste0('dPSI>',dPSI),3,outer = T)
}

pdf('figures/ohnologs.in.ageAS.pdf',w=12,h=8)
par(mfrow=c(2,2),tck=-0.01,mgp=c(1.9,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,1,1))
f(0.5)
f(0.3)
dev.off()

m.age.ens = getAgeASEns(psi.tsm,meta.tsm,0.5,border.stages,'mouse')
plotOhnologsFreq(m.age.ens,c(ohnologs$m$Ohnolog.1.Id,ohnologs$m$Ohnolog.2.Id),main='Mouse',filter = ge.info.m$Mouse_ID[ge.info.m$Age==0])
plotOhnologsFreq(m.age.ens,c(ohnologs$m$Ohnolog.1.Id,ohnologs$m$Ohnolog.2.Id),main='Mouse',filter = ge.info.m$Mouse_ID[ge.info.m$Age> 0])
a1=plotOhnologsFreq(m.age.ens,ge.info.m$Mouse_ID[ge.info.m$Age==0],main='Mouse',filter = ge.info.m$Mouse_ID)
a2=plotOhnologsFreq(m.age.ens,ge.info.m$Mouse_ID[ge.info.m$Age==0],main='Mouse',filter = intersect(ge.info.m$Mouse_ID,c(ohnologs$m$Ohnolog.1.Id,ohnologs$m$Ohnolog.2.Id)))
a3=plotOhnologsFreq(m.age.ens,ge.info.m$Mouse_ID[ge.info.m$Age==0],main='Mouse',filter = setdiff(ge.info.m$Mouse_ID,c(ohnologs$m$Ohnolog.1.Id,ohnologs$m$Ohnolog.2.Id)))

ups.ens = table(unlist(lapply(h.age.ens,'[[','u')))
dws.ens = table(unlist(lapply(h.age.ens,'[[','d')))
uds.ens = table(unlist(lapply(h.age.ens,'[[','b')))
x = unique(unlist(h.age.ens))
t = data.frame(gid=x)
rownames(t) = x
t$up = t$dw = t$b = 0
t[names(ups.ens),'up'] = ups.ens
t[names(dws.ens),'dw'] = dws.ens
t[names(uds.ens),'b' ] = uds.ens
table(t$dw,t$up)
t$dir = 'n'
t$dir[ t$up> 0 & t$dw==0  & t$b==0] = paste0('u',t$up[t$up> 0 & t$dw==0 & t$b==0])
t$dir[ t$up==0 & t$dw> 0  & t$b==0] = paste0('d',t$dw[t$up==0 & t$dw> 0 & t$b==0])
t$dir[(t$up> 0 & t$dw> 0) | t$b> 0] = paste0('b',(pmax(t$up,t$dw)+t$b)[(t$up> 0 & t$dw> 0) | t$b> 0])
t$dir.sum = t$dir
t$dir.sum[t$dir.sum %in% c('u2','u3','u4','u5','u6','u7')] = 'un'
t$dir.sum[t$dir.sum %in% c('d2','d3','d4','d5','d6','d7')] = 'dn'
t$dir.sum[t$dir.sum %in% c('b1','b2','b3','b4','b5','b6','b7')] = 'b'

sort(table(t$dir.sum))
t$hapin1 = hapin2[t$gid,5]
boxplotWithSgn(split(t$hapin1 , t$dir.sum)[c('n','u1','un','d1','dn','b')],notch=T)

ge.info.m = read.csv('/home/mazin/skoltech/projects/evo.devo/input/gene.info.from.marg/Mouse.Indexes.csv')
rownames(ge.info.m) = ge.info.m$Mouse_ID
t = table(ge.info.m$Age,ge.info.m$Mouse_ID %in% c(ohnologs$m$Ohnolog.1.Id,ohnologs$m$Ohnolog.2.Id))
plot(sweep(t,1,apply(t,1,sum),'/')[,2])

fisher.test(table(ge.info.m$Age==0,ge.info.m$Mouse_ID %in% c(ohnologs$m$Ohnolog.1.Id,ohnologs$m$Ohnolog.2.Id)))

fisher.test(table(hapin2$V5>0.8,oh=rownames(hapin2) %in% c(ohnologs$h$Ohnolog.1.Id,ohnologs$h$Ohnolog.2.Id)))
# early all-tissues-alt exons ####
s = 'mouse'
DPSI = 0.3

m = getAgeASchanges(psi.tsm,meta.tsm,DPSI,border.stages,s)
mm = apply(m,2,function(x){x %in% c('u','d')})
psi =psi.tsm[[s]][anns[[s]]$sites=='ad' & apply(mm,1,sum)>0,]
psi.=psi.tsm[[s]][anns[[s]]$sites=='ad',]
tr.compl.m     = getAltExonStat(psi ,meta.tsm,0.1,tissues = c('brain','heart','liver','ovary','testis'),na.as.cnst=TRUE)
tr.compl.m.all = getAltExonStat(psi.,meta.tsm,0.1,tissues = c('brain','heart','liver','ovary','testis'),na.as.cnst=TRUE)
tr.compl.m.    = getAltExonStat(psi ,meta.tsm,0.1,tissues = c('brain','heart','liver','ovary','testis'),na.as.cnst=F)

dim(psi)
m = meta.tsm[colnames(psi),]
m = m[!(m$tissue %in% c('cerebellum','kidney')),]
table(m$days,m$tissue)
psi = psi[,rownames(m)]

areaplot(tr.compl.m,col=rainbow(11))

mmeta = meta.tsm[colnames(psi),]
is.alt = psi > 0.1 & psi < 0.9

alt.early = apply(is.alt[,mmeta$days==10.5],1,sum,na.rm=F)
table(alt.early)


plotTissueAgeProile(apply(psi[alt.early==5,],2,mean,na.rm=T),meta.tsm)

f = alt.early>1 & alt.early<5
f = alt.early == 5
table(f)
c = cor(t(psi[f,]),u='p')
c = cleanNADistMatrix(c)
dim(c)
h = hclust(as.dist(1-c))
c = reorderClustersBySize(cutree(h,k=20))
par(mfrow=c(4,5),tck=-0.01,mgp=c(1.1,0.2,0),mar=c(6,3,1.5,0),oma=c(0,0,2,1))
for(i in 1:max(c)){
	p = psi[names(c)[c==i],]
	#	p = t(scale(t(p)))
	plotTissueAgeProile(apply(p,2,mean,na.rm=T),meta.tsm,main=paste0('c',i,' (',sum(c==i),')'))
}

plotTissueAgeProile(psi[names(c)[c==1][6],],meta.tsm,main=paste0('c',i,' (',sum(c==i),')'))
boxplot(apply(is.na(psi),1,mean) ~ alt.early,notch=T,outline=F)
dev.off()

x=lapply(0:5,function(i){
	gids = unique(anns$mouse[names(alt.early)[!is.na(alt.early) & alt.early==i],'gene_id'])
	gids = intersect(gids,rownames(my.ge.cod$mouse$rpkm))
	apply(my.ge.cod$mouse$rpkm[gids,]>10,1,mean)
	#log2(apply(my.ge.cod$mouse$rpkm[gids,],1,mean))
})
boxplot(x,notch=T,outline=F)
dev.off()


####
m = readRDS('Rdata/mouse.as.u.filtered.Rdata')
m = m$ir[m$seg$sites=='ad',]
ens2seg = revList(seg2ens$mouse)
ohno.seg = rownames(m) %in% unlist(ens2seg[unique(c(ohnologs$m$Ohnolog.1.Id,ohnologs$m$Ohnolog.2.Id))])
table(ohno.seg)
all.mds = cmdscale(1-cor(m            ,u='p'),k=2)
ohn.mds = cmdscale(1-cor(m[ ohno.seg,],u='p'),k=2)
noh.mds = cmdscale(1-cor(m[!ohno.seg,],u='p'),k=2)
ohn.mds[,2] = -ohn.mds[,2]

par(mfrow=c(2,2),tck=-0.01,mgp=c(1.9,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,1,1))
mm = meta[colnames(m),]
SAJR::plotMDS(points=all.mds,col=mm$col,cex=mm$cex,pch=19,main=paste0('All mouse cassette exons (',nrow(m),')'))
SAJR::plotMDS(points=ohn.mds,col=mm$col,cex=mm$cex,pch=19,main=paste0('Ohnolog mouse cassette exons (',sum(ohno.seg),')'))
SAJR::plotMDS(points=noh.mds,col=mm$col,cex=mm$cex,pch=19,main=paste0('Not ohnolog mouse cassette exons (',sum(!ohno.seg),')'))

### clustering #####
age.al = age.al.i[age.al.i$to.remove==0,1:7]
apply(age.al,2,function(x)sum(table(x[x!=''])>1))
age.al$chicken[c(6,8)] = ''

psi.tsm.al = vector('list',7)
names(psi.tsm.al) = names(psi.tsm)
ts = unique(meta$tissue)
m = data.frame(mouse.stage=rep(age.al$mouse,times=length(ts)),tissue=rep(ts,each=nrow(age.al)))
rownames(m) = paste(m$tissue,m$mouse.stage)
d = unique(meta[meta$species=='mouse',c('days','stage')])
m$days = setNames(d$days,d$stage)[m$mouse.stage]

for(s in names(psi.tsm)){
	psi.tsm.al[[s]] = matrix(NA,ncol=nrow(m),nrow=nrow(psi.tsm[[s]]),dimnames = list(rownames(psi.tsm[[s]]),rownames(m)))
	for(i in 1:nrow(m)){
		cat('\r',s,i,'   ')
		stage = age.al[age.al$mouse==m$mouse.stage[i],s]
		for(t in ts){
			st = paste(s,t,stage)
			if(st %in% colnames(psi.tsm[[s]]))
				psi.tsm.al[[s]][,paste(t,m$mouse.stage[i])] = psi.tsm[[s]][,st]
		}
	}
}
sapply(psi.tsm.al,dim)

undef = sapply(psi.tsm.al,function(x)apply(!is.na(x),2,sum)==0)
undef = apply(undef[,c(1,3:6)],1,sum) # I'll not use macaque and chicken for clustering
undef=cbind(m,undef=undef)[undef>0 & m$tissue%in% c('brain','heart','testis'),]
undef
f = !(rownames(m) %in% rownames(undef))
table(f)

dPSI = 0.5
tissues = c('brain','heart','testis')
seg2clust = lapply(names(psi.tsm.al),function(s){
	t = psi.tsm.al[[s]][anns[[s]]$sites=='ad',f & m$tissue %in% tissues]
	mm = m[f & m$tissue %in% tissues,]
	ts = lapply(split(mm,mm$tissue),function(x){x=x[order(x$days),];x$mouse.stage[c(1,nrow(x))]})
	fd = apply(sapply(tissues,function(tis)apply(t[,mm$tissue==tis & mm$mouse.stage %in% ts[[tis]]],1,function(x){max(x,na.rm=TRUE)-min(x,na.rm=TRUE)>dPSI})),1,sum)>0
	fa = apply(sapply(tissues,function(tis)apply(t[,mm$tissue==tis],1,function(x){max(x,na.rm=TRUE)-min(x,na.rm=TRUE)>dPSI})),1,sum)>0
	fs = apply(sapply(unique(mm$mouse.stage),function(stage)apply(t[,mm$mouse.stage==stage,drop=F],1,function(x){max(x,na.rm=TRUE)-min(x,na.rm=TRUE)>dPSI})),1,sum)>0
	rownames(t)[apply(!is.na(t),1,sum)>=10 & (fa | fs)]
	})

names(seg2clust) = names(psi.tsm.al)#[c(1,3:6)]
sapply(seg2clust,length)
seg2clust.list = seg2clust

seg2clust = do.call(rbind,lapply(names(seg2clust.list)[c(1,3:6)],function(s)psi.tsm.al[[s]][seg2clust.list[[s]],f & m$tissue %in% c('brain','heart','testis')]))

table(substr(rownames(seg2clust),1,3),apply(is.na(seg2clust),1,sum)==0)
f = apply(is.na(seg2clust),1,sum)==0
table(f)

corr  = cor(t(seg2clust[f,]),u='p')
table(is.na(corr))
library(cluster)
clusters = list()
clusters$cl14 = pam(1-corr,diss = TRUE,k = 14)
clusters$cl16 = pam(1-corr,diss = TRUE,k = 16)
clusters$cl20 = pam(1-corr,diss = TRUE,k = 20)
saveRDS(clusters,'Rdata/hmrbo.pam.clust.Rdata')
# clusters = readRDS('Rdata/hmrbo.pam.clust.Rdata')
# clusters.sum = lapply(clusters,'[',c('id.med','clustering'))
# saveRDS(clusters.sum,'Rdata/hmrbo.pam.clust.sum.Rdata')
clusters.sum = readRDS('Rdata/hmrbo.pam.clust.sum.Rdata')
t=table(substr(names(clusters.sum$cl14$clustering),1,3),clusters.sum$cl14$clustering)
barplot(sweep(t,2,apply(t,2,sum),'/'))

assign2clusts = function(psi,meds){
	c=cor(t(psi),t(psi[meds,]),u='p')
	apply(c,1,which.max)
}


p = do.call(rbind,lapply(names(seg2clust.list),function(s)psi.tsm.al[[s]][seg2clust.list[[s]],colnames(seg2clust)]))
p = seg2clust[f,]
c = cor(t(p),t(p[clusters.sum$cl20$id.med,]),u='p')
cc = apply(c,1,which.max)
table(clusters.sum$cl20$clustering,cc)
c[clusters.sum$cl20$clustering==1 & cc==3,]


clusters.sum.all = lapply(clusters.sum,function(x){(assign2clusts(p,names(x$clustering)[x$id.med]))})



psi = seg2clust#psi.tsm$human
colnames(psi) = paste('mouse',colnames(psi))
psi = t(scale(t(psi)))
m = meta.tsm

cl = (clusters.sum$cl20$clustering)
par(mfrow=c(4,5),tck=-0.01,mgp=c(1.9,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,0,1))

for(i in 1:max(cl)){
	p=psi[rownames(psi) %in% names(cl)[cl ==i],,drop=F]
	
	plotTissueAgeProile(apply(p,2,mean,na.rm=T),m,main=paste0('c',i,' (',nrow(p),')'))
}
plotTissueAgeProile(p[106,],m,main=paste0('c',i,' (',nrow(p),')'))


psi = seg2clust[f,]
m = meta.tsm[paste('mouse',colnames(psi)),]
rownames(m) = colnames(psi)


### mouse only
s = 'human'
dPSI = 0.5
mp = psi.tsm[[s]][anns[[s]]$sites=='ad',meta.tsm[colnames(psi.tsm[[s]]),'tissue'] %in% c('brain','heart','testis')]
m = meta.tsm[colnames(mp),]
#f1 = apply(mp,1,function(x){max(x,na.rm=T)-min(x,na.rm=T) > 0.5})

fd = getAgeASchanges(setNames(list(mp),s),meta.tsm,dPSI,border.stages,s)
fa = sapply(unique(m$tissue),function(tis){t = mp[,m$tissue==tis];apply(t,1,max,na.rm=T)-apply(t,1,min,na.rm=T)>dPSI})
fs = apply(sapply(unique(m$stage) ,function(s)  {t = mp[,m$stage==s,drop=F];apply(t,1,max,na.rm=T)-apply(t,1,min,na.rm=T)>dPSI}),1,sum)>0

f  = apply(is.na(mp),1,sum) == 0

table(f1,f)
table(fs | fa,f1)
table(fs,f)

cl16l = pam(1-cor(t(mp[f & f1,])),k=16)
cl16a = pam(1-cor(t(mp[f & fa,])),k=16)
cl16s = pam(1-cor(t(mp[f & fs,])),k=16)
cl16sa = pam(1-cor(t(mp[f & (fs | fa),])),k=16)

cl = reorderClustersBySize(cl16sa$clustering)
par(mfrow=c(4,4),tck=-0.01,mgp=c(1.9,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,1,1))
for(i in 1:max(cl)){
	t = mp[names(cl)[cl==i],]
	#t=t(scale(t(t)))
	plotTissueAgeProile(apply(t,2,mean,na.rm=T),meta.tsm,main=paste0('c',i,' (',nrow(t),')'))
}
mtext('no-norm, tissue+age',3,outer = T)

#
l = list()
l$bh = psi.tsm$human[anns$human$sites=='ad',meta.tsm[colnames(psi.tsm$human),'tissue']=='brain']
l$th = psi.tsm$human[anns$human$sites=='ad',meta.tsm[colnames(psi.tsm$human),'tissue']=='testis']
l$bm = psi.tsm$mouse[anns$mouse$sites=='ad',meta.tsm[colnames(psi.tsm$mouse),'tissue']=='brain']
l$tm = psi.tsm$mouse[anns$mouse$sites=='ad',meta.tsm[colnames(psi.tsm$mouse),'tissue']=='testis']
l = lapply(l,function(x)x[,order(meta.tsm[colnames(x),'days'])])

n = 'bm'
table(apply(l[[n]],1,function(x){max(x,na.rm=T)-min(x,na.rm=T) > dPSI}),abs(l[[n]][,1]-l[[n]][,ncol(l[[n]])]) > dPSI)

# GE of SFs: MBNL, RBFOX, NOVA and PTBP ####
# inspired by "Precise temporal regulation of alternative splicing during neural development"
sfs.mm = ens.descr.mm[grep('nova|rbfox|ptbp|mbnl',ens.descr.mm$Associated.Gene.Name,ignore.case = TRUE),]
sfs.mm = sfs.mm[order(sfs.mm$Associated.Gene.Name),]
# ENSMUSG00000108585 - Retired
sfs.mm = sfs.mm[sfs.mm$Ensembl.Gene.ID!='ENSMUSG00000108585',]
sfs.mm = sfs.mm[c(1:3,6:11,4:5),]
o = read.csv('input/ens.orths.txt.gz')
colnames(o) = rownames(species)
o = o[o[,3] %in% sfs.mm$Ensembl.Gene.ID,]
rownames(o) = o$mouse
o = o[sfs.mm$Ensembl.Gene.ID,]
pdf('figures/ptbp.nova.rbfox.mbnl.exp.pdf',w=33,h=21)
par(mfcol=c(7,11),tck=-0.01,mgp=c(1.9,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,1,1))
for(i in 1:nrow(o))
	for(s in colnames(o)){
		if(o[i,s] != '')
			plotTissueAgeProile(ens.ge.cod[[s]]$rpkm[o[i,s],],meta,main=paste(sfs.mm$Associated.Gene.Name[i],s))
		else
			plot.new()
}
dev.off()

###
m = psi.tsm$mouse[anns$mouse$sites=='ad',]
par(mfrow=c(2,2),tck=-0.01,mgp=c(1.9,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,1,1))
hist(apply(m,1,median,na.rm=T),0:100/100)
hist(m,0:100/100)
hist(apply(m,1,function(x){sum(x>0.1 & x<0.9,na.rm=T)/length(x)}),0:100/100)


# dPSI on age ######

age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)

pdf('figures/devAS/dPSI.on.age.devAS>0.2.pdf',w=14,h=12)
par(mfrow=c(6,7),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,1,0))
for(s in rownames(species)[c(3:5,1,6:7)])
	plotdPSIOnAge(s,0.2)
mtext('# of exons with dPSI > 0.2 between adj. stages',outer = T)


for(s in rownames(species)[c(3:5,1,6:7)]){
	xlim=unique(meta.tsm[meta.tsm$species==s,c('age.rank','stage')])
	xlim = setNames(xlim$age.rank,xlim$stage)
	for(t in unique(meta$tissue)){
		c=getAdjStageCor(s,t)
		plotArea(xlim[rownames(c)],c,new=T,xaxt='n',xlab='Age',ylab='Correlation',main=paste(s,t),col=params$tissue.col[t],lwd=3,xlim=range(xlim))
		axis(1,xlim[rownames(c)],rownames(c))
		m = meta.tsm[meta.tsm$species==s,]
		b=m$age.rank[order(abs(species[s,'gestation']-m$days))[1]]
		abline(v=b,lty=3)
	}
}
mtext('correlation of adj. stages',outer = T)

for(s in rownames(species)[c(3:5,1,6:7)]){
	xlim=unique(meta.tsm[meta.tsm$species==s,c('age.rank','stage')])
	xlim = setNames(xlim$age.rank,xlim$stage)
	for(t in unique(meta$tissue)){
		mb = getSplineDevdPSI(s,t,df=4,shuffle = F,dpsi.thr=0.2)
		
		plotArea(mb$x[-length(mb$x)],mb$up,new=T,ylim=range(mb$up,mb$dw),lwd=3,xaxt='n',xlab='Age',ylab='mean dPSI/stage',main=paste(s,t),col=params$tissue.col[t],xlim=range(xlim))
		plotArea(mb$x[-length(mb$x)],mb$dw,col=params$tissue.col[t],new=F,area.den = 20,lwd=3,lty=2)
		axis(1,xlim,names(xlim))
		m = meta.tsm[meta.tsm$species==s,]
		b=m$age.rank[order(abs(species[s,'gestation']-m$days))[1]]
		abline(v=b,lty=3)
	}
}
mtext('mean spline dPSI (df=4)',outer = T)
# for(s in rownames(species)[c(3:5,1,6:7)])
# 	plotdPSIOnAge(s,0.3)
# mtext('dPSI>0.3',outer = T)
# for(s in rownames(species)[c(3:5,1,6:7)])
# 	plotdPSIOnAge(s,0.4)
# mtext('dPSI>0.4',outer = T)
dev.off()

# _compare to GE ######
ge.cnt = read.csv('input/number.of.de.gene.on.dev.from.Margarida_UPDATED.csv')
t=strsplit(ge.cnt$TimeComp2,'[_-]')
ge.cnt$t1 = sapply(t,'[',1)
ge.cnt$t2 = sapply(t,'[',2)
ge.cnt$t1[substr(ge.cnt$t2,nchar(ge.cnt$t2)-2,10)=='wpc'] = paste0(ge.cnt$t1[substr(ge.cnt$t2,nchar(ge.cnt$t2)-2,10)=='wpc'],'wpc')

hs=c(newb='newborn',inf='infant',tod='toddler',sch='school',teen='teenager',yteen='youngteenager',oteen='oldTeenager',ya='youngadult',yma='youngmidage',oma='oldermidage',sen='senior')

setdiff(setdiff(unique(ge.cnt$t1[ge.cnt$species=='human']),meta$stage),names(hs))
setdiff(setdiff(unique(ge.cnt$t2[ge.cnt$species=='human']),meta$stage),names(hs))

ge.cnt$t1[ge.cnt$t1 %in% names(hs)] = hs[ge.cnt$t1[ge.cnt$t1 %in% names(hs)]]
ge.cnt$t2[ge.cnt$t2 %in% names(hs)] = hs[ge.cnt$t2[ge.cnt$t2 %in% names(hs)]]

table(tolower(ge.cnt$t1) %in% tolower(meta$marg.stage),ge.cnt$species)
table(paste(ge.cnt$species,tolower(ge.cnt$t2))%in% paste(meta$species,tolower(meta$marg.stage)),ge.cnt$species)

f = ge.cnt$species=='opossum'
setdiff(ge.cnt$t1[f],meta$marg.stage[meta$species=='opossum'])
setdiff(meta$marg.stage[meta$species=='opossum'],ge.cnt$t1[f])

f = ge.cnt$species=='opossum' & substr(ge.cnt$t1,1,1) == 'P'
ge.cnt$t1[f] = as.numeric(substr(ge.cnt$t1[f],2,20))+14
f = ge.cnt$species=='opossum' & substr(ge.cnt$t1,1,1) == 'e'
ge.cnt$t1[f] = substr(ge.cnt$t1[f],2,20)

f = ge.cnt$species=='opossum' & substr(ge.cnt$t2,1,1) == 'P'
ge.cnt$t2[f] = as.numeric(substr(ge.cnt$t2[f],2,20))+14
f = ge.cnt$species=='opossum' & substr(ge.cnt$t2,1,1) == 'e'
ge.cnt$t2[f] = substr(ge.cnt$t2[f],2,20)
ge.cnt$tissue = tolower(ge.cnt$tissue)

mm = unique(paste(meta$species,meta$marg.stage))
mm = setNames(mm,tolower(mm))
table(paste(ge.cnt$species,tolower(ge.cnt$t2)) %in% names(mm))
ge.cnt$t1 = sapply(strsplit(mm[paste(ge.cnt$species,tolower(ge.cnt$t1))],' '),'[',2)
ge.cnt$t2 = sapply(strsplit(mm[paste(ge.cnt$species,tolower(ge.cnt$t2))],' '),'[',2)

table(paste(ge.cnt$species,ge.cnt$tissue,tolower(ge.cnt$t2)) %in% unlist(lapply(psi.tsm.ms,colnames)),ge.cnt$species)

# I'll perform analysis on Margarida stages, so code below is depricated
# my2marg.stage = unique(meta[,c('species','tissue','stage','marg.stage')])
# mystage = paste(my2marg.stage$species,my2marg.stage$tissue, my2marg.stage$stage)
# mrstage = paste(my2marg.stage$species,my2marg.stage$tissue,my2marg.stage$marg.stage)
# x = split(mystage,mrstage)
# y = split(mrstage,mystage)
# x[sapply(x,length)>1]
# y[sapply(y,length)>1]
# 
# # code below is not completely correct since my stages has ambiguous correspondence to Margarida stages
# my2marg.stage = setNames(paste(my2marg.stage$species,my2marg.stage$tissue,my2marg.stage$stage),tolower(paste(my2marg.stage$species,my2marg.stage$tissue,my2marg.stage$marg.stage)))
# id = tolower(paste(ge.cnt$species,ge.cnt$tissue,ge.cnt$t1))
# id[!(id %in% names(my2marg.stage))]
# rownames(ge.cnt) = my2marg.stage[id]


'rabbit heart 9mpb' %in% rownames(ge.cnt);#my2marg.stage

dpsi = 0.2
pdf('figures/devAS/GE.vs.AS.number.of.events.on.dev.pdf',w=14,h=12)
par(mfrow=c(6,7),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,1,0))
for(sp in rownames(species)[c(3:5,1,6:7)]){
	dpsi.age = calcdPSIonAge(psi.tsm[[sp]],meta.tsm,stages2use=NULL)
	for(t in unique(meta$tissue)){
		f = grep(t,colnames(dpsi.age))
		up=apply(abs(dpsi.age[anns[[sp]]$sites=='ad' & per.tissue.age.qv[[sp]][,t]<0.05,f])>  dpsi,2,function(x){x = x[!is.na(x)];if(length(x)==0){NA}else{sum(x)}})
		up = up[!is.na(up) & names(up) %in% rownames(ge.cnt)]
		plot(ge.cnt[names(up),'total']+1,up+1,log='xy',pch=19,col=params$tissue.col[t],cex=meta.tsm[names(up),'cex'],xlab='DE cnt + 1',ylab='DAS cnt + 1',main=paste(sp,t))
		lm = lm(log(up+1)~log(ge.cnt[names(up),'total']+1))
		p = exp(predict.lm(lm,interval='confidence'))
		plotArea(ge.cnt[names(up),'total']+1,p,col=params$tissue.col[t],lwd=3)
		cr = round(cor(ge.cnt[names(up),'total'],up,m='sp'),3)
		c=10^par('usr')
		text(c[1],c[4],paste0('Sp. cor=',cr),adj = c(0,1))
	}
}
mtext('both directions - both directions',outer=T)
dev.off()

# _peak change 190205 ####
# expression pseudocount is 0.01...
#ens.ge.marg.tsm = readRDS('Rdata/ens.ge.marg.tsm.Rdata')
dpsi.age = lapply(rownames(species),function(s)calcdPSIonAge(psi.tsm[[s]][anns[[s]]$sites=='ad',],meta.tsm,stages2use=NULL))
#dgex.age = lapply(rownames(species),function(s)calcdPSIonAge(log2(ens.ge.marg.tsm[[s]]+1e-1),meta.tsm,stages2use=NULL))
#names(dgex.age) = names(dpsi.age) = rownames(species)

plot(apply(abs(dpsi.age$mouse[,grep('brain',colnames(dpsi.age$mouse))])>0.2,2,sum,na.rm=T),t='b')
plot(apply(abs(dgex.age$mouse[,grep('brain',colnames(dgex.age$mouse))])>1,2,sum,na.rm=T),t='b',ylim=c(0,7000))

pv=od=matrix(NA,ncol=7,nrow=7,dimnames = list(rownames(species),unique(meta$tissue)))
peak.ct = list()
for(s in rownames(species)){
	peak.ct[[s]] = list()
	for(t in unique(meta$tissue)){
		g = dgex.age[[s]][,grep(t,colnames(dgex.age[[s]]))]
		a = dpsi.age[[s]][,grep(t,colnames(dpsi.age[[s]]))]
		cnts = c('all as'=nrow(a))
		fc = apply(!is.na(a),2,sum)==0
		cnts = c(cnts,'stages.removed'=sum(fc))
		a = a[,!fc]
		a = a[apply(is.na(a),1,sum)==0,]
		cnts = c(cnts,'filt as'=nrow(a))
		cat(s,t,cnts,'\n')
		s2e = seg2ens[[s]][rownames(a)]
	
		a = a[sapply(s2e,length)==1,]
		s2e = unlist(s2e[rownames(a)])
		f = s2e %in% rownames(dgex.age[[s]])
		s2e = s2e[f]
		a = a[f,]
		
		bkg.gid = unique(s2e)
		as.gid  = unique(s2e[apply(abs(a)>0.2,1,sum,na.rm=T)>0])
		tab=table(as=factor(bkg.gid %in% as.gid,levels=c(FALSE,TRUE)),
							ge=factor(bkg.gid %in% rownames(g)[apply(abs(g) > 1,1,sum,na.rm=T)>0],levels=c(FALSE,TRUE)))
		peak.ct[[s]][[t]] = tab
		ft=fisher.test(tab)
		pv[s,t] = ft$p.value
		od[s,t] = ft$estimate
	}
}
od[od==0] = NA
barplot(log2(od),beside=T,las=3,col=paste0(rep(params$tissue.col,each=7),ifelse(pv<0.001,'FF','22')),border=NA,ylab='log2(odds ratio)')

plot(apply(abs(a)>0.2,2,sum),t='b')
plot(apply(abs(dgex.age[[s]][s2e[rownames(a)],grep('brain',colnames(dgex.age$mouse))])>1,2,sum),t='b')
plot(apply(abs(dgex.age[[s]][s2e[rownames(a)],grep('brain',colnames(dgex.age$mouse))])>1,2,sum),t='b')
plot(apply(abs(dgex.age[[s]][,grep('brain',colnames(dgex.age$mouse))])>1,2,sum),t='b')

as.that.ge = sapply(peak.ct,function(x)sapply(x,function(y){
	y[2,2]/sum(y[2,])
}))

na.that.ge = sapply(peak.ct,function(x)sapply(x,function(y){
	y[1,2]/sum(y[1,])
}))


ge.that.as = sapply(peak.ct,function(x)sapply(x,function(y){
	y[2,2]/sum(y[,2])
}))

ne.that.as = sapply(peak.ct,function(x)sapply(x,function(y){
	y[2,1]/sum(y[,1])
}))

pdf('figures/devAS/overlap.peakAS-peakGE.pdf',w=14,h=5)
par(mfrow=c(1,2),tck=-0.01,mgp=c(1.2,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,0))
barplot(t(as.that.ge),beside=T,col=paste0(rep(params$tissue.col,each=7),ifelse(pv<0.001,'FF','22')),border=NA,ylab='fraq of peakAS genes that are peakGE',main='AS')
barplot(t(na.that.ge),beside=T,den=30,border=NA,add=T,col='black')
legend(8,0.7,fill='black',den=c(-1,30),legend=c('peak','not peak'))

barplot(t(ge.that.as),beside=T,col=paste0(rep(params$tissue.col,each=7),ifelse(pv<0.001,'FF','22')),border=NA,ylab='fraq of peakGE genes that are peakAS',main='GE')
b=barplot(t(ne.that.as),beside=T,den=30,border=NA,add=T,col='black')
dev.off()

tab=table(as=factor(bkg.gid %in% as.gid.o,levels=c(FALSE,TRUE)),
					ge=factor(bkg.gid %in% rownames(g)[apply(g> 1,1,sum,na.rm=T)>0],levels=c(FALSE,TRUE)))

# _compare peak genes at peak stages ####
# m2m = NULL
# ts = unique(meta$tissue)
# for(s in rownames(species)){
# 	m = meta[meta$species==s,]
# 	m = sort(sapply(split(m$days,m$stage),mean))
# 	m = names(m)[-length(m)]
# 	r = data.frame(my=paste(s,rep(ts,each=length(m)),rep(m,times=length(ts))),marg=paste0('T',rep(1:length(m),times=length(ts))))
# 	r = r[r$my %in% rownames(meta.tsm),]
# 	m2m = rbind(m2m,r)
# }
# 
# 
# rownames(m2m) = m2m$my
# table(m2m$my %in% rownames(ge.cnt))
# cmn = intersect(m2m$my , rownames(ge.cnt))
# table(m2m[cmn,2] == ge.cnt[cmn,'timeCom'])
# t = m2m[cmn,2] != ge.cnt[cmn,'timeCom']
# m2m[cmn[t],][1:4,]
# ge.cnt[cmn[t],][1:4,]
# # lets just fix to margarida stages
# m2m[cmn[t],2] = ge.cnt[cmn[t],'timeCom']
# my2marg = m2m


# meta.tsm.ms = readRDS('Rdata/meta.tsm.ms.Rdata')
# dpsi.age = lapply(rownames(species),function(s)calcdPSIonAge(psi.tsm.ms[[s]],meta.tsm.ms,stages2use=NULL))
# names(dpsi.age) = rownames(species)

peak.ge.gids = lapply(paste0(firstToupper(rownames(species)[-2]),'_annotations'),function(s){
	r = lapply(firstToupper(unique(meta$tissue)),function(t){
		d = paste0('processed/GE.from.marg/Annotations_peaks_of_dev_change/',s,'/',t,'/')
		fls = gsub('.txt','',list.files(d,pattern = 'txt'))
		r = lapply(fls,function(f)readLines(paste0(d,f,'.txt')))
		names(r) = gsub('a','',gsub(t,'',fls))
		r
	})
	names(r) = unique(meta$tissue)
	r
})
names(peak.ge.gids) = rownames(species)[-2]
sapply(peak.ge.gids,function(x)sapply(x,length))
# remove rabbit, because genes lists from Margarida have T15 point, but in gene counts for peack changes that I use to find stages compared the maximal point is T14
# peak.ge.gids = peak.ge.gids[-4]

#rename to match csv table
s='rabbit'
for(t in names(peak.ge.gids[[s]])){
	mar.stage = as.numeric(gsub('Down|Up|T','',names(peak.ge.gids[[s]][[t]])[-1]))
	o = c('Bckground',paste0('T',rep(mar.stage,each=2),rep(c('Up','Down'),times=length(mar.stage))))
	peak.ge.gids[[s]][[t]] = peak.ge.gids[[s]][[t]][o]
	f = t != 'ovary' | mar.stage != 14
	if(sum(f)>0)
		mar.stage[f] = mar.stage[f] - 1
	names(peak.ge.gids[[s]][[t]])[-1] =paste0('T',rep(mar.stage,each=2),rep(c('Up','Down'),times=length(mar.stage)))
}		
n = names(peak.ge.gids$human$testis)
names(peak.ge.gids$human$testis)[n=='T19Down'] = 'T18Down'
names(peak.ge.gids$human$testis)[n=='T19Up'] = 'T18Up'
# comp numbers to ge.cnt
# r = NULL 
# for(s in names(peak.ge.gids)){
# 	for(t in names(peak.ge.gids[[s]])){
# 		mar.stage = sort(unique(as.numeric(gsub('Down|Up|T','',names(peak.ge.gids[[s]][[t]])[-1]))))
# 		if(length(mar.stage)>0){
# 			r = rbind(r,data.frame(species=s,tissue=t,stage=paste0('T',mar.stage),up  =sapply(peak.ge.gids[[s]][[t]][paste0('T',mar.stage,'Up')],length),
# 																								down=sapply(peak.ge.gids[[s]][[t]][paste0('T',mar.stage,'Down')],length)))
# 		}
# 	}		
# }
# rownames(r) = paste(r$species,r$tissue,r$stage)
# ge.cnt.rn = paste(ge.cnt$species,ge.cnt$tissue,ge.cnt$timeCom)
# r[setdiff(rownames(r) , ge.cnt.rn),]
# inx = match(rownames(r),ge.cnt.rn)
# table(r$up==ge.cnt$up[inx],r$species)
# table(r$down==ge.cnt$down[inx],r$species)
# 
# f = (r$up!=ge.cnt$up[inx] | r$down!=ge.cnt$down[inx])
# cbind(r[f,],ge.cnt[inx[f],])
# 
# f=r$species!='rabbit'
# plot(r$up[f],ge.cnt$up[inx[f]])
# plot(r$down[f],ge.cnt$down[inx[f]])
# 
# comp genes

rownames(ge.cnt) = paste(ge.cnt$species,ge.cnt$tissue,ge.cnt$timeCom)
peak.as.gids = list()
dPSI = 0.2
for(s in names(peak.ge.gids)){
	print(s)
	peak.as.gids[[s]] = list()
	for(t in names(peak.ge.gids[[s]])){
		peak.as.gids[[s]][[t]] = list()
		#ge
		mar.stage = sort(unique(as.numeric(gsub('Down|Up|T','',names(peak.ge.gids[[s]][[t]])[-1]))))
		mar.stage = mar.stage[paste0(s,' ',t,' T',mar.stage) %in% paste(ge.cnt$species,ge.cnt$tissue,ge.cnt$timeCom)]
		if(length(mar.stage)<1)
			next
		o = c('Bckground',paste0('T',rep(mar.stage,each=2),rep(c('Up','Down'),times=length(mar.stage))))
		peak.ge.gids[[s]][[t]] = peak.ge.gids[[s]][[t]][o]
		names(peak.ge.gids[[s]][[t]])[1] = 'bkg'
		#as
		m2m = ge.cnt[ge.cnt$tissue==t & ge.cnt$species==s,]
		m2m = tolower(setNames(m2m$t1,m2m$timeCom))
		#colnames(dpsi.age[[s]]) = tolower(colnames(dpsi.age[[s]]))

		# only devAS - might be confounded by expression level
		sgn = per.tissue.age.qv[[s]][,t] < 0.05 & age.dpsi[[s]][,t]>0.2
		sgn[is.na(sgn)] = FALSE
		#dp = dpsi.age[[s]][anns[[s]]$sites=='ad' & sgn,paste(s,t,m2m[paste0('T',mar.stage)]),drop=FALSE]
		psi = psi.tsm.ms[[s]][anns[[s]]$sites=='ad' & sgn,]
		as.bkg = unique(unlist(seg2ens[[s]][rownames(per.tissue.age.qv[[s]])[anns[[s]]$sites=='ad' & !is.na(per.tissue.age.qv[[s]][,t]) & !is.na(age.dpsi[[s]][,t])]]))
		peak.as.gids[[s]][[t]]$bkg = as.bkg #intersect(as.bkg,peak.ge.gids[[s]][[t]]$bkg)


		# just use dPSI threshold
		#dp = dpsi.age[[s]][anns[[s]]$sites=='ad',paste(s,t,m2m[paste0('T',mar.stage)]),drop=FALSE]
		#dp = dp[apply(is.na(dp),1,sum)==0,,drop=FALSE] #I use only exons that expressed at all stages studied
		#peak.ge.gids[[s]][[t]]$bkg = peak.as.gids[[s]][[t]]$bkg = intersect(unique(unlist(seg2ens[[s]][rownames(dp)])),peak.ge.gids[[s]][[t]]$bkg)

		for(i in 1:length(mar.stage)){
			cmp = ge.cnt[paste0(s,' ',t,' T',mar.stage[i]),]
			dp = psi[,paste(s,t,tolower(cmp$t2))] - psi[,paste(s,t,tolower(cmp$t1))]
			f = !is.na(dp)
			peak.ge.gids[[s]][[t]][[paste0('T',mar.stage[i],'Both')]] = unique(c(peak.ge.gids[[s]][[t]][[paste0('T',mar.stage[i],'Up')]],peak.ge.gids[[s]][[t]][[paste0('T',mar.stage[i],'Down')]]))
			peak.as.gids[[s]][[t]][[paste0('T',mar.stage[i],'Both')]] = unique(unlist(seg2ens[[s]][names(dp)[f & abs(dp)>dPSI]]))
			peak.as.gids[[s]][[t]][[paste0('T',mar.stage[i],'Down')]] = unique(unlist(seg2ens[[s]][names(dp)[f & dp< -dPSI]]))
			peak.as.gids[[s]][[t]][[paste0('T',mar.stage[i],'Up')]]   = unique(unlist(seg2ens[[s]][names(dp)[f & dp>dPSI]]))
		}
		o = c('bkg',paste0('T',rep(mar.stage,each=3),rep(c('Up','Down','Both'),times=length(mar.stage))))
		peak.ge.gids[[s]][[t]] = peak.ge.gids[[s]][[t]][o]
	}
}



checkSetsOverlap = function(s1,s2,b){
	pv=odds =cnt = matrix(NA,nrow=length(s1),ncol=length(s2),dimnames=list(names(s1),names(s2)))
	for(i in 1:length(s1))
		for(j in 1:length(s2)){
			t = fisher.test(factor(b %in% s1[[i]],levels = c(TRUE,FALSE)),factor(b %in% s2[[j]],levels = c(TRUE,FALSE)))
			pv[i,j] = t$p.value
			odds[i,j] = t$estimate
			cnt[i,j] = sum(b %in% s1[[i]] & b %in% s2[[j]])
		}
	list(pv=pv,odd=odds,cnts=cnt)
}

pdf('figures/devAS/peak.change.overlaps.all-species-tissues20201021.devAS_.pdf',w=7*6,h=5*6)
par(mfrow=c(5,7),tck=-0.01,mgp=c(3.3,0.2,0),mar=c(2.5,4.5,1.5,0),oma=c(0,0,1,0))
for(s in names(peak.ge.gids)[c(2,3,1,4,5)])
	for(t in names(peak.ge.gids[[s]])){
		if(length(peak.as.gids[[s]][[t]])==0){
			plot.new()
			next
		}
		z=checkSetsOverlap(peak.ge.gids[[s]][[t]],peak.as.gids[[s]][[t]],intersect(peak.as.gids[[s]][[t]]$bkg,peak.ge.gids[[s]][[t]]$bkg))
		o = round(z$odd,2)
		o[1,] = ''
		o[,1] = ''
		imageWithText(z$pv*sign(log(z$odd)),t = paste0(z$cnt,'\n',o),col=c('white','gray','cyan','blue','red','yellow','gray'),breaks=c(-1,-1,-0.05,-0.001,0,0.001,0.05,1),las=1,xlab='GE',ylab='AS',main=paste(s,t),names.as.labs = TRUE)
	}
for(s in names(peak.ge.gids)[c(2,3,1,4,5)])
	for(t in names(peak.ge.gids[[s]])){
		if(length(peak.as.gids[[s]][[t]])==0){
			plot.new()
			next
		}
		z = checkSetsOverlap(peak.ge.gids[[s]][[t]],peak.as.gids[[s]][[t]],intersect(peak.as.gids[[s]][[t]]$bkg,peak.ge.gids[[s]][[t]]$bkg))
		z = lapply(z,function(x)x[rownames(x)[grep('Both|bkg',rownames(x))],colnames(x)[grep('Both|bkg',colnames(x))],drop=F])
		o = round(z$odd,2)
		o[1,] = ''
		o[,1] = ''
		imageWithText(z$pv*sign(log(z$odd)),t = paste0(z$cnt,'\n',o),col=c('white','gray','cyan','blue','red','yellow','gray'),breaks=c(-1,-1,-0.05,-0.001,0,0.001,0.05,1),las=1,xlab='GE',ylab='AS',main=paste(s,t),names.as.labs = TRUE)
	}
dev.off()

checkSetsOverlap1 = function(s1,s2,b){
	s1 = intersect(s1,b)
	s2 = intersect(s2,b)
	ft = fisher.test(factor(b %in% s1,levels = c(TRUE,FALSE)),factor(b %in% s2,levels = c(TRUE,FALSE)))
	data.frame(pv=ft$p.value,or=ft$estimate,intersect=length(intersect(s1,s2)),len1=length(s1),len2=length(s2))
}

r = NULL
for(s in names(peak.ge.gids)){
	for(t in names(peak.ge.gids[[s]])){
		stages = unique(gsub('Down|Up|Both','',names(peak.ge.gids[[s]][[t]])[-1]))
		for(stage in stages){
			for(gdir in c('Up','Down','Both')){
				for(adir in c('Up','Down','Both')){
					tmp = checkSetsOverlap1(peak.ge.gids[[s]][[t]][[paste0(stage,gdir)]],peak.as.gids[[s]][[t]][[paste0(stage,adir)]],intersect(peak.as.gids[[s]][[t]]$bkg,peak.ge.gids[[s]][[t]]$bkg))
					r = rbind(r,cbind(species=s,tissue=t,stage=stage,GE.dir=gdir,AS.dir=adir,tmp))
				}
			}
		}
	}
}

dim(r)
table(is.na(r))

r[1:13,]

w = r[r$AS.dir=='Both' & r$GE.dir=='Both',]
dim(w)
w[1:2,]
hist(w$pv)
w$fdr = p.adjust(w$pv,m='BH')
table(w$fdr<0.05,w$GE.dir)
boxplot(log2(w$or) ~ w$GE.dir)
f = w$fdr<0.05
table(w$or[f]>1,w$GE.dir[f])
hist(log2(w$or[f]))


hist((w$intersect/pmin(w$len1,w$len2))[f])

rownames(w) = NULL
w$GE.dir = w$AS.dir = NULL
w = w[,c('species','tissue','stage','len1','len2','intersect','or','pv','fdr')]
w[1:2,]
colnames(w)[c(3,4,5,6,7)] = c('Time interval','GE.cnt','AS.cnt','GE-and-AS.cnt','odds ratio')

# m = unique(meta[,c('species','tissue','stage','paper.stages')])
# m = setNames(m$paper.stages,paste(m$species,m$tissue,m$stage))
# mar2my = setNames(m[rownames(ge.cnt)],paste(ge.cnt$species,ge.cnt$tissue,ge.cnt$timeCom))

mar2my = setNames(ge.cnt$TimeComp2,paste(ge.cnt$species,ge.cnt$tissue,ge.cnt$timeCom))
w = cbind(w[,1:3],stages=mar2my[paste(w$species,w$tissue,w$`Time interval`)],w[,-1:-3])

openxlsx::write.xlsx(w,'figures/paper.figures/6/nature.review/peak.changes.overlap.xlsx')

hist(r$pv[r$GE.dir=='Up'])
hist(r$pv[r$GE.dir=='Down'])
hist(log2(r$or[r$GE.dir=='Up']))
hist(log2(r$or[r$GE.dir=='Down']))
r$qv = p.adjust(r$pv,m='BH')
r[r$qv<0.05 & r$GE.dir=='Down',]
f=r$qv<0.05
boxplot(log2(r$or[f]) ~ r$GE.dir[f])
boxplot(log2(r$or[f]) ~ r$AS.dir[f] + r$GE.dir[f])
table(r$AS.dir[f] , r$GE.dir[f])
table(r$GE.dir[f],r$or[f]>1)
hist(r$intersect/(r$len1+r$len2-r$intersect),100)
f = r$AS.dir=='Both'
hist((r$intersect[f]/(r$len2[f])),100)

sapply(mb.as.gids,length)

sapply(mb.as.gids,function(x)table(x %in% mb$Background))

mb.as.gids = lapply(mb.as.gids,intersect,y=mb$Background)

z2=t(apply(sapply(mb,function(x)table(intersect(x,mb.as.gids$t2.bkg) %in% mb.as.gids$t2.both))[2:1,],2,my.binom.test))
z11=t(apply(sapply(mb,function(x)table(intersect(x,mb.as.gids$t11.bkg) %in% mb.as.gids$t11.both))[2:1,],2,my.binom.test))


par(mfrow=c(2,1))
b=barplot(z2[,1],ylim=range(z2,0))
segments(b,z2[,2],b,z2[,3])

b=barplot(z11[,1],ylim=range(z11,0))
segments(b,z11[,2],b,z11[,3])


## max change
ens.ge.marg.tsm = readRDS('Rdata/ens.ge.marg.tsm.Rdata')
age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)
f=anns[[sp]]$sites=='ad' & per.tissue.age.qv$mouse[,'brain'] < 0.05 & abs(age.dpsi$mouse[,'brain'])>0.2 & anns[[sp]]$cod.gene
f[is.na(f)] = F
x=dpsi.age[f,grep('brain',colnames(dpsi.age))]
dim(x)
x = x[apply(is.na(x),1,sum)==0,]
as.max.stage=apply(abs(x),1,which.max)
plot(table(as.max.stage))

my.ge.cod.tsm = readRDS('Rdata/my.ge.cod.tsm.Rdata')
g=my.ge.cod.tsm$mouse[,grep('brain',colnames(my.ge.cod.tsm$mouse))]
g = g[,order(meta.tsm[colnames(g),'days'])]
g = log2(g+1)
g = -(g[,-ncol(g)] - g[,-1])

ge.max.stage=sapply(apply(abs(g[anns$mouse[rownames(x),'gene_id'],]),1,which.max),function(x){ifelse(length(x)==0,NA,x)})
plot(table(ge.max.stage))
plot(ge.max.stage,(as.max.stage-ge.max.stage))

plot(table((as.max.stage-ge.max.stage)[ge.max.stage==11]))
image(table(ge.max.stage,as.max.stage))

## Are inlcusion/exlusion have 0/1 PSI in other tissues
age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)

pdf('figures/devAS/inclusion-exclusion.vs.median.PSI.pdf',w=14,h=12)
par(mfrow=c(6,7),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,1,0))
for(s in rownames(species)[c(3:5,1,6:7)]){
	for(t in unique(meta$tissue)){
		f = anns[[s]]$sites=='ad'
		qv = per.tissue.age.qv[[s]]
		qv[is.na(age.dpsi[[s]])] = NA
		f = f & !is.na(qv[,t]) & abs(age.dpsi[[s]][,t]) > 0.5 & qv[,t] < 0.05
		#med=apply(psi.tsm[[s]][f,meta.tsm[colnames(psi.tsm[[s]]),'tissue'] %in% c('liver','ovary','kidney')],1,median,na.rm=T)
		med=apply(psi.tsm[[s]][f,!grepl(t,colnames(psi.tsm[[s]]))],1,median,na.rm=T)
		dir = ifelse(age.dpsi[[s]][f,t]>0,'u','d')
		hist(med[dir=='u'],0:20/20,col='#FF000090',border = NA,xlab='median(PSI) in other tissues',main=paste(s,t))
		hist(med[dir=='d'],0:20/20,col='#0000FF90',border = NA,add=T)
		legend(0.5,par('usr')[4],xjust=0.5,yjust=1,fill=c('#FF000090','#0000FF90'),
					 legend=paste0(c('incl','exlc'),' (',c(sum(dir=='u'),sum(dir=='d')),')'))
	}
}
mtext('|dPSI| > 0.5',outer = T)

par(mfrow=c(6,7),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,1,0))
for(s in rownames(species)[c(3:5,1,6:7)]){
	for(t in unique(meta$tissue)){
		f = anns[[s]]$sites=='ad'
		qv = per.tissue.age.qv[[s]]
		qv[is.na(age.dpsi[[s]])] = NA
		f = f & !is.na(qv[,t]) & abs(age.dpsi[[s]][,t]) > 0.2 & qv[,t] < 0.05
		#med=apply(psi.tsm[[s]][f,meta.tsm[colnames(psi.tsm[[s]]),'tissue'] %in% c('liver','ovary','kidney')],1,median,na.rm=T)
		med=apply(psi.tsm[[s]][f,!grepl(t,colnames(psi.tsm[[s]]))],1,median,na.rm=T)
		dir = ifelse(age.dpsi[[s]][f,t]>0,'u','d')
		hist(med[dir=='u'],0:20/20,col='#FF000090',border = NA,xlab='median(PSI) in other tissues',main=paste(s,t))
		hist(med[dir=='d'],0:20/20,col='#0000FF90',border = NA,add=T)
		legend(0.5,par('usr')[4],xjust=0.5,yjust=1,fill=c('#FF000090','#0000FF90'),
					 legend=paste0(c('incl','exlc'),' (',c(sum(dir=='u'),sum(dir=='d')),')'))
	}
}
mtext('|dPSI| > 0.2',outer = T)
dev.off()

## exont
library(ontologyIndex)
o = get_OBO('processed/exon.onthology/exont.obo')
seg2exont = readRDS('Rdata/seg2exont.Rdata')


# Intrinsically_unstructured_polypeptide_region
# EXONT:000074

iupr = names(seg2exont)[sapply(seg2exont,function(x)'EXONT:000074' %in% x)]

cnst = rownames(all.anns$human)[all.anns$human$sites=='ad' & all.anns$human$type=='EXN' & rownames(all.anns$human) %in% names(seg2exont)]
cnst.iupr=mean(cnst %in% iupr)


bins = c(-1,-0.7,-0.5,-0.2,0,0.2,0.5,0.7,1)
names=paste0('[',bins[-length(bins)],',',bins[-1],')')
substr(names[length(names)],nchar(names[length(names)]),100) = ']'

findInterval(bins,bins,all.inside = T)

iupr.stat = lapply(unique(meta$tissue),function(t){
	p = age.dpsi$human[anns$human$sites=='ad' & rownames(anns$human) %in% names(seg2exont),t]
	p = p[!is.na(p)]
	p = setNames(findInterval(p,bins,all.inside = T),names(p))
	r = t(apply(table(p,names(p) %in% iupr)[,c('TRUE','FALSE')],1,my.binom.test))
	rownames(r) = names
	r
	})
names(iupr.stat) = unique(meta$tissue)

pdf('figures/paper.figures/3.2/2018.07.04/05.IUPR.in.devAS.pdf',w=9,h=9)
par(mfrow=c(3,3),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,1,0))
for(t in names(iupr.stat)){
	plotArea(1:8,iupr.stat[[t]],col=params$tissue.col[t],lwd=3,ylim=range(cnst.iupr,unlist(iupr.stat)),new=T,xaxt='n',xlab='dPSI',ylab='IUPR fraction',main=t)
	axis(1,1:8,names)
	abline(h=cnst.iupr,lwd=3,lty=2,col='gray')
}
mtext('Intrinsically unstructured polypeptide regions in devAS',outer = T)
dev.off()

age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)

hdevas02 = abs(age.dpsi$human) > 0.2 & per.tissue.age.qv$human<0.05
hdevas02[is.na(hdevas02)] = FALSE

hdevas05 = abs(age.dpsi$human) > 0.5 & per.tissue.age.qv$human<0.05
hdevas05[is.na(hdevas05)] = FALSE

hdevas07 = abs(age.dpsi$human) > 0.7 & per.tissue.age.qv$human<0.05
hdevas07[is.na(hdevas07)] = FALSE


f = anns$human$sites=='ad' & anns$human$cod!='n' & rownames(anns$human) %in% names(seg2exont)

idr02 = apply(hdevas02,2,function(x)my.binom.test(table(rownames(hdevas02)[f & x] %in% iupr)[c('TRUE','FALSE')]))
idr05 = apply(hdevas05,2,function(x)my.binom.test(table(rownames(hdevas05)[f & x] %in% iupr)[c('TRUE','FALSE')]))
idr07 = apply(hdevas07,2,function(x)my.binom.test(table(rownames(hdevas07)[f & x] %in% iupr)[c('TRUE','FALSE')]))


pdf('figures/devAS/IUPR.in.devAS.pre.tissue.pdf',w=9,h=4)
par(mfrow=c(1,3),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(4.5,2.5,1.5,0),oma=c(0,0,0,1))
b=barplot(idr02[1,],col=params$tissue.col,border=NA,las=3,ylab='fraction of IDR',main='|dPSI| > 0.2',ylim=range(0,.8))
abline(h=cnst.iupr,lty=2)
segments(b,idr02[2,],b,idr02[3,])

b=barplot(idr05[1,],col=params$tissue.col,border=NA,las=3,ylab='fraction of IDR',main='|dPSI| > 0.5',ylim=range(0,.8))
abline(h=cnst.iupr,lty=2)
segments(b,idr05[2,],b,idr05[3,])


b=barplot(idr07[1,],col=params$tissue.col,border=NA,las=3,ylab='fraction of IDR',main='|dPSI| > 0.7',ylim=range(0,.8))
abline(h=cnst.iupr,lty=2)
segments(b,idr07[2,],b,idr07[3,])
dev.off()


pdf('figures/devAS/micro-vs-macro.exons.dPSI.pdf',w=7,h=3*7)
par(mfrow=c(7,2),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,2.5,0),oma=c(0,0,0,1))
for(t in unique(meta$tissue)){
	f = anns$human$sites=='ad' & abs(age.dpsi$human[,t]) > 0.2
	f[is.na(f)] = FALSE
	l = anns$human$length <=27
	
	br = seq(-1,1,by=0.1)
	hist(age.dpsi$human[f & !l,t],br,col='#00000060',border=NA,freq=F,xlab='dPSI',main=c(t,'all CE with |dPSI| > 0.2'))
	hist(age.dpsi$human[f &  l,t],br,col=paste0(params$tissue.col[t],'77'),border=NA,freq=F,add=T)
	
	f = anns$human$sites=='ad' & abs(age.dpsi$human[,t]) > 0.2 & per.tissue.age.qv$human[,t]<0.05
	f[is.na(f)] = FALSE
	hist(age.dpsi$human[f & !l,t],br,col='#00000060',border=NA,freq=F,xlab='dPSI',main=c(t,'all CE with |dPSI| > 0.2 & pv.afj < 0.05'))
	hist(age.dpsi$human[f &  l,t],br,col=paste0(params$tissue.col[t],'99'),border=NA,freq=F,add=T)
}
dev.off()

########### dPSI per tissue/direction ######

age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)

pdf('figures/paper.figures/3.2/2018.07.18/dPSI.distr.pdf',w=9,h=9)
#pdf('figures/devAS/dPSI.distr.pdf',w=9,h=9)
s = 'mouse'
for(s in rownames(species)[-2]){
	tissues = unique(meta$tissue)
	dpsi=lapply(tissues,function(t){
		dpsi=age.dpsi[[s]][anns[[s]]$sites=='ad' & per.tissue.age.qv[[s]][,t]<0.05 & abs(age.dpsi[[s]][,t])>=0,t]
		dpsi=dpsi[!is.na(dpsi)]
		list(up=dpsi[dpsi>0],down=abs(dpsi[dpsi<0]))
		})
	names(dpsi) = tissues
	
	
	br = 0:20/20
	freq=FALSE
	par(mfrow=c(3,3),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,1,0))
	for(t in tissues){
		if(freq)
			ylim = range(0,hist(dpsi[[t]]$up,br,plot=F)$counts,hist(dpsi[[t]]$down,br,plot=F)$counts)
		else
			ylim = range(0,hist(dpsi[[t]]$up,br,plot=F)$density,hist(dpsi[[t]]$down,br,plot=F)$density)
		hist(dpsi[[t]]$up,br,col=paste0(params$tissue.col[t],55),freq=freq,ylim=ylim,border = NA,xlab='dev dPSI',main=t)
		hist(dpsi[[t]]$down,br,add=T,col=params$tissue.col[t],den=20,freq=freq,border = NA,xlab='')
		legend('topright',fill=params$tissue.col[t],den=c(-1,20),legend=paste0(c('inclusion','exclusion'),' (',c(length(dpsi[[t]]$up),length(dpsi[[t]]$down)),')'))
	}
	boxplot(unlist(dpsi,recursive = F),col=rep(params$tissue.col[tissues],each=2),pars=list(den=c(-1,40)),notch=T,las=3)
	mtext(s,outer = T)
}
dev.off()


melt(t(dpsi.age[1:5,1:5]))

# devPeakChange: AS vs GE ####
ens.ge.marg.tsm = readRDS('Rdata/ens.ge.marg.tsm.Rdata')


getAsGeAdjStageChangeMatrix = function(d,dpsi.thr,l2fc.thr,as.max=FALSE,ge.max=FALSE,rm.na=FALSE){
	if(rm.na){
		f = apply(is.na(d$as),1,sum)==0
		d$as = d$as[f,]
		d$ge = d$ge[f,]
	}
	r=matrix(0,ncol=ncol(d$as)+1,nrow=ncol(d$as)+1,dimnames=list(paste0('AS.',c('nochange',colnames(d$as))),paste0('GE.',c('nochange',colnames(d$as)))))
	d$as[is.na(d$as)]
	for(s in 1:nrow(d$as)){
		if(as.max){
			a = order(abs(d$as[s,]),decreasing = T)[1]
			a = a[!is.na(d$as[s,a]) & abs(d$as[s,a])>dpsi.thr]
		}else
			a = which(!is.na(d$as[s,]) & abs(d$as[s,])>dpsi.thr)
		
		if(ge.max){
			g = order(abs(d$ge[s,]),decreasing = T)[1]
			g = g[abs(d$ge[s,g])>l2fc.thr]
		}else
			g = which(abs(d$ge[s,]) > l2fc.thr)
		if(length(a)==0) a=0
		if(length(g)==0) g=0
		r[a+1,g+1] = r[a+1,g+1] + 1
	}
	r
}

matrixFT = function(m){
	pv = m
	pv[,] = NA
	or = tounion = tomin = tomax = pv
	cs = apply(m,2,sum)
	rs = apply(m,1,sum)
	as = sum(m)
	for(r in 1:nrow(m))
		for(c in 1:ncol(m)){
			t = matrix(c(as-cs[c]-rs[r]+m[r,c],rs[r]-m[r,c],cs[c]-m[r,c],m[r,c]),ncol=2)
			ft = fisher.test(t)
			pv[r,c] = ft$p.value
			or[r,c] = ft$estimate
			tounion[r,c] = m[r,c]/as
			tomin[r,c] = m[r,c]/min(rs[r],cs[r])
			tomax[r,c] = m[r,c]/max(rs[r],cs[r])
		}
	
	list(m=m,pv=pv,or=or,tounion=tounion,tomin=tomin,tomax=tomax)
}


pdf('figures/devAS/peak.changes-AS-all.vs.GE-all.pdf',w=3*5,h=6*5)
par(mfrow=c(6,3),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(3,2.5,1.5,3),oma=c(0,0,1,0))
sp = 'mouse'
for(sp in c('mouse','rat','human','rabbit','opossum','chicken')){
	print(sp)
	tissue = 'brain'
	dPSI = 0.3
	l2fc = 1
	stageinx = unique(meta.tsm[meta.tsm$species==sp & meta.tsm$tissue==tissue,c('stage','days')])
	stageinx = stageinx[order(stageinx$days),]
	
	x = getAdjStageChange(psi.tsm[[sp]][anns[[sp]]$sites=='ad' & anns[[sp]]$cod !='n',],log2(ens.ge.marg.tsm[[sp]]+0.1),meta.tsm,tissue,seg2ens[[sp]])
	gcnt = as.numeric(table(x$stage[abs(x$de)>l2fc])[stageinx$stage[-nrow(stageinx)]])
	acnt = as.numeric(table(x$stage[abs(x$dpsi)>dPSI])[stageinx$stage[-nrow(stageinx)]])
	plot(gcnt,t='b',col='blue',lwd=3,xlab='age',yaxt='n',xaxt='n',ylab='',main=paste(sp,tissue))
	mtext(paste0('# exons with |l2fc| > ',l2fc),2,1.5,col='blue')
	mtext(paste0('# exons with |dPSI| > ',dPSI),4,1.5,col='red')
	points(acnt/max(acnt)*max(gcnt),t='b',col='red',lwd=3)
	axis(1,1:length(gcnt),stageinx$stage[-nrow(stageinx)])
	axis(2,col='blue',col.axis='blue')
	at = 1:20*100
	axis(4,at/max(acnt)*max(gcnt),at,col='red',col.axis='red')
	
	
	x = getAdjStageChange(psi.tsm[[sp]][anns[[sp]]$sites=='ad' & anns[[sp]]$cod !='n',],log2(ens.ge.marg.tsm[[sp]]+0.1),meta.tsm,tissue,seg2ens[[sp]],melt=FALSE)
	t02 = getAsGeAdjStageChangeMatrix(x,dPSI,l2fc,as.max=F,ge.max = F,rm.na=TRUE)
	colnames(t02) = gsub('GE.','',colnames(t02))
	rownames(t02) = gsub('AS.','',rownames(t02))
	t02. = matrixFT(t02)
	imageWithText(t02.$pv,t02,names.as.labs = T,xlab='AS',ylab='GE',col=c('red','yellow','white'),breaks=c(0,0.05/length(t02.$pv),0.05,1),main=paste(sp,tissue))
	t02 = t02[-1,-1]
	t02. = matrixFT(t02)
	imageWithText(t02.$pv,t02,names.as.labs = T,xlab='AS',ylab='GE',col=c('red','yellow','white'),breaks=c(0,0.05/length(t02.$pv),0.05,1),main=paste(sp,tissue))
}
dev.off()

b.as.max = dpsi.age[,grep('brain',colnames(dpsi.age))]
b.as.max = b.as.max[apply(is.na(b.as.max),1,sum)==0,]
stages = sapply(strsplit(colnames(b.as.max),' ',T),'[',3)
b.as.max = do.call(rbind,apply(b.as.max,1,function(x){o = order(abs(x),decreasing = T)[1];data.frame(stage=stages[o],dpsi=x[o])}))

s2e = seg2ens$mouse[rownames(b.as.max)]
table(sapply(s2e,length))
b.as.max = b.as.max[sapply(s2e,length)==1,]
b.as.max$gene_id = unlist(seg2ens$mouse[rownames(b.as.max)])

b.ge.max = dexp.age[,grep('brain',colnames(dexp.age))]
stages = sapply(strsplit(colnames(b.ge.max),' ',T),'[',3)
b.ge.max = do.call(rbind,apply(b.ge.max,1,function(x){o = order(abs(x),decreasing = T)[1];data.frame(stage=stages[o],dpsi=x[o])}))
table(b.as.max$gene_id %in% rownames(b.ge.max))
b.as.max = b.as.max[b.as.max$gene_id %in% rownames(b.ge.max),]
b.as.max$gene.stage = b.ge.max[b.as.max$gene_id,'stage']
b.as.max$gene.l2fc = b.ge.max[b.as.max$gene_id,'dpsi']

hist(b.as.max$gene.l2fc)

f = abs(b.as.max$dpsi)>0.3 
as.cnt = as.numeric(table(factor(b.as.max$stage[f],levels=stages))[stages])
de.cnt = as.numeric(table(factor(b.as.max$gene.stage[f],levels=stages))[stages])
plot(as.cnt)
plot(de.cnt)
plot(as.cnt,de.cnt)



s = c('-',stages)
as = factor(ifelse(abs(b.as.max$dpsi     )>0.3    ,b.as.max$stage,'-'),levels=s)
ge = factor(ifelse(abs(b.as.max$gene.l2fc)>1 ,b.as.max$gene.stage,'-'),levels=s)

t=table(as,ge)[s,s]
pdf('figures/paper.figures/2018.12.10/01.peakChange.mouse.brain.pdf',w=7,h=7)
imageWithText(log(t+1),t,las=2,names.as.labs=T,xlab='AS',ylab='GE',main='Mouse brain. dPSI > 0.3, l2FC>1')
dev.off()



z=matrixFT(t)
y=matrixFT(t[-1,-1])
image(log(y$pv))

# N3 fraction for exons in [0.25-0.75] on tissue/age ####
sp = 'mouse'
ff = apply(is.na(psi.tsm[[sp]]),1,sum)==0
table(ff)
z=sapply(1:ncol(psi.tsm[[sp]]),function(i){
	f = ff & anns[[sp]]$sites=='ad' & anns[[sp]]$cod.gene & !is.na(psi.tsm[[sp]][,i]) & psi.tsm[[sp]][,i]>0.2 & psi.tsm[[sp]][,i]<0.8
	table(anns[[sp]]$length[f] %% 3)
})

zz = z[1,]/apply(z,2,sum)
names(zz) = colnames(psi.tsm[[sp]])
plotTissueAgeProile(zz,meta.tsm)



