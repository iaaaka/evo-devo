options(stringsAsFactors = FALSE)
source('code/r.functions/paper.figures.F.R')
source('code/r.functions/load.all.data.F.R')
library(SAJR)

species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
border.stages = readRDS('Rdata/border.stages.Rdata')
anns = readRDS('Rdata/anns.Rdata')

orth.seg.ad.tsm = readRDS('Rdata/orth.seg.ad.tsm.Rdata')
orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
orth.irs = readRDS('Rdata/not.cassette/orth.irs.Rdata')
orth.aas = readRDS('Rdata/not.cassette/orth.aas.Rdata')
orth.dds = readRDS('Rdata/not.cassette/orth.dds.Rdata')

age.al.i = read.table('input/margarida.age.alignment.final.2017.05.17.improved.tab',sep='\t',header = T)[,c(rownames(species),'to.remove')]

psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')

orth.per.tissue.age.qv = readRDS('Rdata/orth.per.tissue.age.qv.Rdata')
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')

params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)


## all AS type ancient orth MDSs (pearson + spearman) ####
psi=list()
psi$ad = do.call(cbind,sapply(orth.seg.ad,'[[','ir'))
psi$ad = psi$ad[apply(sapply(orth.seg.ad,function(x)x$seg$type)=='ALT',1,sum)==7,]
psi$aa = do.call(cbind,sapply(orth.aas,'[[','ir'))
psi$dd = do.call(cbind,sapply(orth.dds,'[[','ir'))
psi$da = do.call(cbind,sapply(orth.irs,'[[','ir'))

cor.p = lapply(psi,function(x)cor(x,u='p',m='p'))
cor.s = lapply(psi,function(x)cor(x,u='p',m='s'))

mds.p = lapply(cor.p,function(x) cmdscale(1-x,k=2))
mds.s = lapply(cor.s,function(x) cmdscale(1-x,k=2))
mdss = lapply(c(mds.p,mds.s),function(x)x[rownames(meta),])

pdf('figures/paper.figures/3.2/2018.07.04/01.MDS.orth.4types.pdf',w=12,h=6)
par(mfrow=c(2,4),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
names = paste0(rep(c('CE','AA','AD','RI'),times=2),' (',rep(sapply(psi,nrow),times=2),', ',rep(c('Pearson','Spearman'),each=4),')')
for(i in 1:length(mdss))
	plot(mdss[[i]],xlab='Dim 1',ylab='Dim 2',col=meta$col,cex=meta$cex,pch=meta$pch,main=names[i])
dev.off()

## fraction of devAS (per tissue/species in all)
age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)
age.dpsi$macaque = cbind(age.dpsi$macaque[,1:5],ovary=NaN,age.dpsi$macaque[,6,drop=FALSE])

pdf('figures/paper.figures/3.2/2018.07.04/02.sgn.counts.4types.pdf',w=12,h=6)
astypes = c(CE='ad',AA='aa',AD='dd',RI='da')
par(mfcol=c(2,4),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
for(i in 1:length(astypes)){
	stat = sapply(names(anns),function(s)apply(!is.na(per.tissue.age.qv[[s]][anns[[s]]$sites == astypes[i],]),2,sum))
	b=barplot(t(stat),beside = T,col=rep(params$tissue.col,each=7),ylab='Number of events',main=names(astypes)[i],names.arg=substr(rownames(stat),1,1))
	text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.7)
	
	sgn02 = sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][anns[[s]]$sites == astypes[i],] < 0.05 & abs(age.dpsi[[s]][anns[[s]]$sites == astypes[i],])>0.2,2,sum,na.rm=T))
	sgn05 = sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][anns[[s]]$sites == astypes[i],] < 0.05 & abs(age.dpsi[[s]][anns[[s]]$sites == astypes[i],])>0.5,2,sum,na.rm=T))
	sgn02 = sgn02/stat*100
	sgn05 = sgn05/stat*100
	
	b=barplot(t(sgn02),den=50,beside = T,col=rep(params$tissue.col,each=7),ylab='% of events',main='% of devAS (of tested)',names.arg=substr(rownames(stat),1,1))
	barplot(t(sgn05),beside = T,col=rep(params$tissue.col,each=7),ylab='',main='',add=T,xaxt='n')
	text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.7)
}
dev.off()

## microexons, early/late; conservation? #####
# this code is deprocated, use microexons.R
orth.age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(orth.seg.ad.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(orth.age.dpsi) = rownames(species)

ff = function(a,s,l=Inf,st='ad'){a[[s]]$sites %in% st & a[[s]]$length<=l}
dPSI=0.5
sgn02u = sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][ff(anns,s),] < 0.05 & age.dpsi[[s]][ff(anns,s),]>  dPSI,2,sum,na.rm=T))
sgn02d = sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][ff(anns,s),] < 0.05 & age.dpsi[[s]][ff(anns,s),]< -dPSI,2,sum,na.rm=T))

sgn02u.me = sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][ff(anns,s,l=27),] < 0.05 & age.dpsi[[s]][ff(anns,s,l=27),]>  dPSI,2,sum,na.rm=T))
sgn02d.me = sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][ff(anns,s,l=27),] < 0.05 & age.dpsi[[s]][ff(anns,s,l=27),]< -dPSI,2,sum,na.rm=T))


micro.timing = lapply(rownames(species), function(s){
	f = anns[[s]]$sites=='ad' & per.tissue.age.qv[[s]][,'brain']<0.05 & age.dpsi[[s]][,'brain'] > 0.5
	f = !is.na(f) & f
	dpsi.age = calcdPSIonAge(psi.tsm[[s]],meta.tsm)
	dpsi.age = dpsi.age[,grep('brain',colnames(dpsi.age))]
	max.brain.dpsi = apply(dpsi.age[f,],1,function(x)meta.tsm[colnames(dpsi.age)[which.max(x)],'days'])
	t = table(micro=anns[[s]]$length[f]<=27,max.brain.dpsi>=meta.tsm[paste(s,'brain',age.al.i[10,s]),'days'])
})

pdf('figures/paper.figures/3.2/2018.07.04/03.microexons.pdf',w=12,h=12)
par(mfrow=c(3,3),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
b=barplot(t(sgn02u.me/sgn02u)*100,beside = T,col=rep(params$tissue.col,each=7),ylab='% of microexons',main='% of microexons in devAS (dPSI > 0.2)',names.arg=substr(rownames(sgn02u.me),1,1))
text(b,t(sgn02u.me/sgn02u)*100+0.2,t(sgn02u.me),srt=90,cex=0.5,adj=c(0,0.5))
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.7)
plotPanelLetter('A')

b=barplot(t(sgn02d.me/sgn02d)*100,beside = T,col=rep(params$tissue.col,each=7),ylab='% of microexons',main='% of microexons in devAS (dPSI < -0.2)',names.arg=substr(rownames(sgn02u.me),1,1),ylim=range(0,sgn02u.me/sgn02u*100,na.rm=T))
text(b,t(sgn02d.me/sgn02d)*100+0.2,t(sgn02d.me),srt=90,cex=0.5,adj=c(0,0.5))
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.7)
plotPanelLetter('B')

names(micro.timing) = rownames(species)
ft=sapply(micro.timing[-2],function(t)as.numeric(fisher.test(t)[c('estimate','p.value')]))
r=sapply(micro.timing[-2],function(x)x[2,]/(x[1,]+x[2,]))*100
b=barplot(r,beside = T,ylab='% of microexons',legend.text = c('before both','after birth'),ylim=c(0,50))
plotPanelLetter('C')

orth.micro = apply(sapply(orth.seg.ad,function(x)x$seg$length<=27),1,sum)==7
anc = apply(sapply(orth.seg.ad,function(x)x$seg$type=='ALT'),1,sum)==7
tis = unique(meta$tissue)[-7]

f = anc &  orth.per.tissue.age.qv$mouse[,'brain']<0.05 & orth.age.dpsi$mouse[,'brain'] > 0.2
f = !is.na(f) & f

dpsi.mb = orth.seg.ad.tsm$mouse[,'mouse brain 9wpb']-orth.seg.ad.tsm$mouse[,'mouse brain 10.5']
orth.macro.eq = names(dpsi.mb) %in% names(equalyseDistrs(dpsi.mb[f & orth.micro],dpsi.mb[f & !orth.micro]))


hist(dpsi.mb[f & orth.micro],col='#FF000080',freq=F,border = NA,xlab='mouse brain dPSI',main='Mouse devAS dPSI distr',ylim=c(0,2.5))
hist(dpsi.mb[f & !orth.micro],col='#0000FF80',add=T,freq=F,border = NA)
hist(dpsi.mb[f & orth.macro.eq],col='#00FF0080',add=T,freq=F,border = NA,xlab='')
legend('topleft',title='Ancient, brain devAS',fill=c('red','green','blue'),
			 legend=paste0(c('microexons','macroexons, same dPSI','macroexons'),' (',c(sum(f & orth.micro),sum(f & orth.macro.eq),sum(f & !orth.micro)),')'))
plotPanelLetter('D')

lineCors = function(x,y,...){
	sapply(1:nrow(x),function(i)cor(x[i,],y[i,],...))
}

#correlation across adult tissues
sps = c('rat','rabbit','human','opossum','chicken')
ma=sapply(lapply(sps,function(s)lineCors(orth.seg.ad.tsm$mouse[f & !orth.micro  ,paste('mouse',tis,'9wpb')],orth.seg.ad.tsm[[s]][f & !orth.micro  ,paste(s,tis,border.stages[[s]][tis,2])],u='p')),function(x){x=x[!is.na(x)];sd=sd(x)/sqrt(length(x));r=mean(x);c(r,r-2*sd,r+2*sd)})
me=sapply(lapply(sps,function(s)lineCors(orth.seg.ad.tsm$mouse[f & orth.macro.eq,paste('mouse',tis,'9wpb')],orth.seg.ad.tsm[[s]][f & orth.macro.eq,paste(s,tis,border.stages[[s]][tis,2])],u='p')),function(x){x=x[!is.na(x)];sd=sd(x)/sqrt(length(x));r=mean(x);c(r,r-2*sd,r+2*sd)})
mi=sapply(lapply(sps,function(s)lineCors(orth.seg.ad.tsm$mouse[f &  orth.micro  ,paste('mouse',tis,'9wpb')],orth.seg.ad.tsm[[s]][f &  orth.micro  ,paste(s,tis,border.stages[[s]][tis,2])],u='p')),function(x){x=x[!is.na(x)];sd=sd(x)/sqrt(length(x));r=mean(x);c(r,r-2*sd,r+2*sd)})

plotArea(1:ncol(ma),t(ma),col='blue',t='b',ylim=range(mi,ma),lwd=3,new=T,xlab='',ylab='mean Pearson corr. to mouse',main='Corelation across adult tissues (no testis)',xaxt='n')
plotArea(1:ncol(ma),t(mi),col='red',t='b',ylim=range(mi,ma),lwd=3)
plotArea(1:ncol(me),t(me),col='green',t='b',ylim=range(mi,ma),lwd=3)
axis(1,1:5,sps)
plotPanelLetter('E')
#correlation across brain development
age.al.i. = age.al.i[age.al.i$human %in% meta.tsm$stage[meta.tsm$species=='human' & meta.tsm$tissue=='brain'],]

ma=sapply(lapply(sps[-5],function(s)lineCors(orth.seg.ad.tsm$mouse[f & !orth.micro   ,paste('mouse brain',age.al.i.$mouse)],orth.seg.ad.tsm[[s]][f & !orth.micro   ,paste(s,'brain',age.al.i.[,s])],u='p')),function(x){x=x[!is.na(x)];sd=sd(x)/sqrt(length(x));r=mean(x);c(r,r-2*sd,r+2*sd)})
mi=sapply(lapply(sps[-5],function(s)lineCors(orth.seg.ad.tsm$mouse[f &  orth.micro   ,paste('mouse brain',age.al.i.$mouse)],orth.seg.ad.tsm[[s]][f &  orth.micro   ,paste(s,'brain',age.al.i.[,s])],u='p')),function(x){x=x[!is.na(x)];sd=sd(x)/sqrt(length(x));r=mean(x);c(r,r-2*sd,r+2*sd)})
me=sapply(lapply(sps[-5],function(s)lineCors(orth.seg.ad.tsm$mouse[f &  orth.macro.eq,paste('mouse brain',age.al.i.$mouse)],orth.seg.ad.tsm[[s]][f &  orth.macro.eq,paste(s,'brain',age.al.i.[,s])],u='p')),function(x){x=x[!is.na(x)];sd=sd(x)/sqrt(length(x));r=mean(x);c(r,r-2*sd,r+2*sd)})

plotArea(1:ncol(ma),t(ma),col='blue',t='b',ylim=range(mi,ma),lwd=3,new = T,xlab='',ylab='mean Pearson corr. tp mouse',main='Corelation across brain development',xaxt='n')
plotArea(1:ncol(ma),t(mi),col='red',t='b',ylim=range(mi,ma),lwd=3)
plotArea(1:ncol(me),t(me),col='green',t='b',ylim=range(mi,ma),lwd=3)
axis(1,1:5,sps)
plotPanelLetter('F')

ma=apply(sapply(orth.age.dpsi,function(x){x[f & !orth.micro,'brain']})>0.2,2,function(x){x=x[!is.na(x)];my.binom.test(sum(x),sum(!x))})[,sps[-5]]
mi=apply(sapply(orth.age.dpsi,function(x){x[f &  orth.micro,'brain']})>0.2,2,function(x){x=x[!is.na(x)];my.binom.test(sum(x),sum(!x))})[,sps[-5]]
me=apply(sapply(orth.age.dpsi,function(x){x[f &  orth.macro.eq,'brain']})>0.2,2,function(x){x=x[!is.na(x)];my.binom.test(sum(x),sum(!x))})[,sps[-5]]


plotArea(1:ncol(ma),t(ma)*100,col='blue',t='b',ylim=range(mi,ma)*100,lwd=3,new = T,xlab='',ylab='% of exons with dPSI > 0.2',main='Conservation of devAS',xaxt='n')
plotArea(1:ncol(ma),t(mi)*100,col='red',t='b',ylim=range(mi,ma),lwd=3)
plotArea(1:ncol(me),t(me)*100,col='green',t='b',ylim=range(mi,ma),lwd=3)
axis(1,1:5,sps)
plotPanelLetter('G')
dev.off()
