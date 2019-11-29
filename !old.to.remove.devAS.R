#setwd('~/skoltech/projects/evo.devo/')
options(stringsAsFactors = FALSE)
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
border.stages = readRDS('Rdata/border.stages.Rdata')
seg2ens = readRDS('Rdata/seg2ens.Rdata')
ens.ge.cod = readRDS('Rdata/ens.ge.cod.Rdata')
# orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
# orth.seg.ad.tsm = readRDS('Rdata/orth.seg.ad.tsm.Rdata')
born.exn.sajr = readRDS('Rdata/born.exn.sajr.Rdata')
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')

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
# 	tmp = tmp[,(colnames(tmp$ir) %in% rownames(meta)[meta$days < my.sex.maturation[s]])]
# 	per.tissue.age.bm.qv[[s]] = sapply(unique(meta$tissue),function(t){testASAge(tmp,meta,t,min.cov.sams=0.6)})
# 	colnames(per.tissue.age.bm.qv[[s]]) = unique(meta$tissue)
# 	dimnames(per.tissue.age.bm.qv[[s]]) = setNames(dimnames(per.tissue.age.bm.qv[[s]]),NULL)
# }
# saveRDS(per.tissue.age.bm.qv,'Rdata/per.tissue.age.bm.qv.Rdata')

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


pdf('figures/devAS/ad.sgn.stat',w=20,h=12)
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

# compare to GE
ge.cnt = read.csv('input/number.of.de.gene.on.dev.from.Margarida.csv')
t=strsplit(ge.cnt$TimeComp2,'[_-]')
ge.cnt$t1 = sapply(t,'[',1)
ge.cnt$t2 = sapply(t,'[',2)
ge.cnt$t1[substr(ge.cnt$t2,nchar(ge.cnt$t2)-2,10)=='wpc'] = paste0(ge.cnt$t1[substr(ge.cnt$t2,nchar(ge.cnt$t2)-2,10)=='wpc'],'wpc')
hs=c(newb='newborn',inf='infant',tod='toddler',sch='school',teen='teenager',yteen='youngteenager',oteen='oldTeenager',ya='youngadult',yma='youngmidage',oma='oldermidage')
ge.cnt$t1[ge.cnt$t1 %in% names(hs)] = hs[ge.cnt$t1[ge.cnt$t1 %in% names(hs)]]
f = ge.cnt$species=='opossum' & substr(ge.cnt$t1,1,1) == 'P'
ge.cnt$t1[f] = as.numeric(substr(ge.cnt$t1[f],2,20))+14
f = ge.cnt$species=='opossum' & substr(ge.cnt$t1,1,1) == 'e'
ge.cnt$t1[f] = substr(ge.cnt$t1[f],2,20)
ge.cnt$tissue = tolower(ge.cnt$tissue)


my2marg.stage = unique(meta[,c('species','tissue','stage','marg.stage')])
my2marg.stage = setNames(paste(my2marg.stage$species,my2marg.stage$tissue,my2marg.stage$stage),tolower(paste(my2marg.stage$species,my2marg.stage$tissue,my2marg.stage$marg.stage)))
id = tolower(paste(ge.cnt$species,ge.cnt$tissue,ge.cnt$t1))
id[!(id %in% names(my2marg.stage))]
rownames(ge.cnt) = my2marg.stage[id]

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

sp = 'mouse'
t = 'brain'
dpsi.age = calcdPSIonAge(psi.tsm[[sp]],meta.tsm,stages2use=NULL)
ge.cnt[ge.cnt$species==sp & ge.cnt$tissue=='brain' & ge.cnt$timeCom=='T11',]

fls = gsub('.txt|Brain','',list.files('processed/GE.from.marg/Annotations_peaks_of_dev_change/Mouse_annotations/Brain',pattern = 'txt'))
mb=lapply(fls,function(f)readLines(paste0('processed/GE.from.marg/Annotations_peaks_of_dev_change/Mouse_annotations/Brain/Brain',f,'.txt')))
names(mb) = fls
mb$no.change = setdiff(mb$Background,unlist(mb[-1]))
mb = mb[c(1,6,5,4,3,2)]
mb$both.up = intersect(mb$T2Up,mb$T11Up)
mb$T2Up.o = setdiff(mb$T2Up,mb$both.up)
mb$T11Up.o = setdiff(mb$T11Up,mb$both.up)
sapply(mb,length)

mb.t2 = dpsi.age[anns[[sp]]$sites=='ad','mouse brain 11.5']
mb.t2 = mb.t2[!is.na(mb.t2)]

mb.t11 = dpsi.age[anns[[sp]]$sites=='ad','mouse brain 3dpb']
mb.t11 = mb.t11[!is.na(mb.t11)]

mb.as.gids = list()
mb.as.gids$t2.bkg = unique(unlist(seg2ens[[sp]][names(mb.t2)]))
mb.as.gids$t2.up = unique(unlist(seg2ens[[sp]][names(mb.t2)[mb.t2 >  0.2]]))
mb.as.gids$t2.dw = unique(unlist(seg2ens[[sp]][names(mb.t2)[mb.t2 < -0.2]]))
mb.as.gids$t2.both = unique(unlist(seg2ens[[sp]][names(mb.t2)[abs(mb.t2) > 0.2]]))

mb.as.gids$t11.bkg = unique(unlist(seg2ens[[sp]][names(mb.t11)]))
mb.as.gids$t11.up = unique(unlist(seg2ens[[sp]][names(mb.t11)[mb.t11 >  0.2]]))
mb.as.gids$t11.dw = unique(unlist(seg2ens[[sp]][names(mb.t11)[mb.t11 < -0.2]]))
mb.as.gids$t11.both = unique(unlist(seg2ens[[sp]][names(mb.t11)[abs(mb.t11) > 0.2]]))


checkSetsOverlap = function(s1,s2,b){
	pv=odds = matrix(NA,nrow=length(s1),ncol=length(s2),dimnames=list(names(s1),names(s2)))
	for(i in 1:length(s1))
		for(j in 1:length(s2)){
			t = fisher.test(factor(b %in% s1[[i]],levels = c(TRUE,FALSE)),factor(b %in% s2[[j]],levels = c(TRUE,FALSE)))
			pv[i,j] = t$p.value
			odds[i,j] = t$estimate
		}
	list(pv=pv,odd=odds)
}

b=intersect(intersect(mb.as.gids$t11.bkg,mb.as.gids$t2.bkg),mb$Background)
length(b)
z=checkSetsOverlap(mb,mb.as.gids,b)
imageWithText(log2(z$odd[-1,-c(1,5)]),names.as.labs = T,xlab='',ylab='',col=c('blue','white','red'),breaks=c(-100,-0.2,0.2,100))

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



# devPeakChange: AS vs GE ####
sp = 'mouse'
ens.ge.marg.tsm = readRDS('Rdata/ens.ge.marg.tsm.Rdata')
dpsi.age = calcdPSIonAge(psi.tsm[[sp]][anns[[sp]]$sites=='ad' & anns[[sp]]$cod !='n',],meta.tsm,stages2use=NULL)
dexp.age = calcdPSIonAge(log2(ens.ge.marg.tsm[[sp]]+0.1),meta.tsm,stages2use=NULL)

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
