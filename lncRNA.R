options(stringsAsFactors = FALSE)
# source('code/r.functions/ad.on.ge.F.R')
source('code/r.functions/load.all.data.F.R')
# source('code/r.functions/paper.figures.4.F.R')
# source('~/skoltech/r.code/util.R')

library(SAJR)

anns = readRDS('Rdata/anns.Rdata')
seg2ens = readRDS('Rdata/seg2ens.Rdata')
species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')
border.stages = readRDS('Rdata/border.stages.Rdata')
orth.seg.ad.all.id = readRDS('Rdata/orth.seg.ad.all.id.Rdata')
params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)
age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)
age.dpsi$macaque = cbind(age.dpsi$macaque[,1:5],ovary=NaN,age.dpsi$macaque[,6,drop=FALSE])
# new ####
all.anns = readRDS('Rdata/all.anns.plus.cds.pos.Rdata')
anns = lapply(names(anns),function(s)all.anns[[s]][rownames(anns[[s]]),])
names(anns) = names(all.anns)
lvs = c('5utr','p5utr','cds','nc-in-cds','3utr','p3utr','lncRNA','pc-lncRNA','-')
sapply(anns,function(x)table(factor(x$cds[x$sites=='ad'],levels=lvs)))

# _sign stat ####

f = function(a,q,dp,cds,thr,comp,qv.thr=0.05,sites='ad'){
	sapply(names(a),function(s){
		f = a[[s]]$sites %in% sites & a[[s]]$cds %in% cds
		apply(q[[s]][f,] <= qv.thr & comp(dp[[s]][f,],thr),2,sum,na.rm=T)
		})
}

pdf('figures/lncRNA/devAS.stat.pdf',w=8,h=3*9)
par(mfrow=c(9,2),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,0,1))
for(l in lvs){
	sgn02 = f(anns,per.tissue.age.qv,age.dpsi,l,0.2,`>=`) + f(anns,per.tissue.age.qv,age.dpsi,l,-0.2,`<=`)
	sgn05 = f(anns,per.tissue.age.qv,age.dpsi,l,0.5,`>=`) + f(anns,per.tissue.age.qv,age.dpsi,l,-0.5,`<=`)
	teste = f(anns,per.tissue.age.qv,age.dpsi,l,0,`>=`,1) + f(anns,per.tissue.age.qv,age.dpsi,l,0,`<=`,1)
	
	short.tiss = substr(rownames(sgn02),1,1)
	
	b=barplot(t(sgn02),den=50,beside = T,col=rep(params$tissue.col,each=7),ylab='Number of devAS exons',main=l,names.arg=short.tiss)
	barplot(t(sgn05),beside = T,col=rep(params$tissue.col,each=7),ylab='',main='',add=T,xaxt='n')
	text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.7)
	
	b=barplot(t(sgn02/teste)*100,den=50,beside = T,col=rep(params$tissue.col,each=7),ylab='% of devAS exons',main=l,names.arg=short.tiss,ylim=c(0,30))
	barplot(t(sgn05/teste)*100,beside = T,col=rep(params$tissue.col,each=7),ylab='',main='',add=T,xaxt='n')
	text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.7)
}
dev.off()

# _MDS ####
pdf('figures/lncRNA/MDSs.pdf',w=7*3,h=3*9)
par(mfrow=c(9,7),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,0,1))
for(l in lvs){
	for(s in rownames(species)){
		print(paste(s,l))
		p = psi.tsm[[s]][anns[[s]]$sites=='ad' & anns[[s]]$cds == l,]
		c = cor(p,u='p')
		if(sum(is.na(c))>0){
			plot.new()
		}else{
			mds = cmdscale(1-c,k=2)
			m = meta.tsm[colnames(p),]
			cex = (m$cex-min(m$cex))/(max(m$cex)-min(m$cex))*2 + 0.5 
			plot(mds,col=m$col,pch=19,cex=cex,xlab='Dim 1',ylab='Dim 2',main=paste0(s,'; ', l,' (',nrow(p),')'))
		}
	}
}
dev.off()
# compare to GE (?) get it from

# old ####
mouse.lnc = loadSAData('input/lncRNA/mouse.lncRNA.sajr')
mouse.lnc = setSplSiteTypes(mouse.lnc,'input/lncRNA/mouse.lncRNA.sajr')
mouse.lnc = mouse.lnc$seg


human.lnc = loadSAData('input/lncRNA/human.lncRNA.sajr')
human.lnc = setSplSiteTypes(human.lnc,'input/lncRNA/human.lncRNA.sajr')
human.lnc = human.lnc$seg

m=all.anns$mouse
t = setNames(rownames(m),paste(m$chr_id,m$start,m$stop,m$strand))
mouse.lnc$my.sid = t[paste(mouse.lnc$chr_id,mouse.lnc$start,mouse.lnc$stop,mouse.lnc$strand)]


h=all.anns$human
t = setNames(rownames(h),paste(h$chr_id,h$start,h$stop,h$strand))
human.lnc$my.sid = t[paste(human.lnc$chr_id,human.lnc$start,human.lnc$stop,human.lnc$strand)]


table(mouse.lnc$type,mouse.lnc$position)
table(mouse.lnc$type,mouse.lnc$sites)[,c('ad','aa','dd','da')]

table(human.lnc$type,human.lnc$position)
table(human.lnc$type,human.lnc$sites)[,c('ad','aa','dd','da')]


cmn = intersect(na.omit(mouse.lnc$my.sid),rownames(anns$mouse))
table(anns$mouse[cmn,'sites'])
seg2ens$mouse[cmn[anns$mouse[cmn,'cod']=='c']]


apply(abs(age.dpsi$mouse[cmn,]) > 0.2 & per.tissue.age.qv$mouse[cmn,]<0.05,2,sum,na.rm=T)
apply(abs(age.dpsi$mouse      ) > 0.5 & per.tissue.age.qv$mouse      <0.05,2,sum,na.rm=T)

x=cmn[apply(abs(age.dpsi$mouse[cmn,]) > 0.2 & per.tissue.age.qv$mouse[cmn,]<0.05,1,sum,na.rm=T)>0]
f = anns$mouse$sites=='ad'
plotTissueAgeProile(apply(psi.tsm$mouse[f & rownames(anns$mouse) %in% x,],2,mean,na.rm=T),meta.tsm)
plotTissueAgeProile(apply(psi.tsm$mouse[f,],2,mean,na.rm=T),meta.tsm)

sgn = abs(age.dpsi$mouse[cmn,]) > 0.5 & per.tissue.age.qv$mouse[cmn,]<0.05
sgn[is.na(sgn)] = F
anns$mouse[cmn,][sgn[,7],]
plotTissueAgeProile(psi.tsm$mouse['mou.26911.s3',],meta.tsm)
seg2ens$mouse['mou.10173.s3']
anns$mouse['mou.30236.s12',]

z=apply(abs(age.dpsi$mouse[cmn,]) > 0.5 & per.tissue.age.qv$mouse[cmn,]<0.05,1,sum,na.rm=T)
seg2ens$mouse[names(z[z>0])]
intersect(cmn,orth.seg.ad.all.id[,'mouse'])
all.anns$mouse['mou.26911.s2',]
lnc.orth=orth.seg.ad.all.id[orth.seg.ad.all.id[,3] %in% mouse.lnc$my.sid,]
all.anns$human[lnc.orth[,1],]
all.anns$mouse[lnc.orth[,3],]

plotTissueAgeProile(psi.tsm$mouse['mou.4901.s7',],meta.tsm)


cmn = intersect(na.omit(human.lnc$my.sid),rownames(anns$human))
table(anns$human[cmn,'sites'])
apply(abs(age.dpsi$human[cmn,]) > 0.2 & per.tissue.age.qv$human[cmn,]<0.05,2,sum,na.rm=T)

