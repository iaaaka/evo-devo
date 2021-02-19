options(stringsAsFactors = FALSE)
source('code/r.functions/paper.figures.F.R')
source('code/r.functions/load.all.data.F.R')
library(SAJR)

species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
border.stages = readRDS('Rdata/border.stages.Rdata')
anns = readRDS('Rdata/anns.Rdata')
all.anns = readRDS('Rdata/all.anns.Rdata')
age.al.i = read.table('input/margarida.age.alignment.final.2017.05.17.improved.tab',sep='\t',header = T)[,c(rownames(species),'to.remove')]
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')

orth.per.tissue.age.qv = readRDS('Rdata/orth.per.tissue.age.qv.Rdata')
orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
orth.seg.ad.tsm = readRDS('Rdata/orth.seg.ad.tsm.Rdata')


psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)

age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)
age.dpsi$macaque = cbind(age.dpsi$macaque[,1:5],ovary=NaN,age.dpsi$macaque[,6,drop=FALSE])

devas02 = lapply(rownames(species),function(s){
	apply(per.tissue.age.qv[[s]]<0.05 & abs(age.dpsi[[s]])>0.2,1,sum,na.rm=T)>0
	})
names(devas02) = rownames(species)
sapply(devas02,table)

me.stat = lapply(all.anns,function(a){
	l = a$length[a$sites=='ad']
	t = table(len=ceiling(l / 3)*3,frame=l %% 3)
	cbind(exon.len = as.numeric(rownames(t)),t)
})


plot(table(all.anns$mouse$length[all.anns$mouse$sites=='ad' & all.anns$mouse$type=='EXN' & all.anns$mouse$cod=='n'])[1:150])
plot(table(all.anns$mouse$length[all.anns$mouse$sites=='ad' & all.anns$mouse$type=='EXN' & all.anns$mouse$cod=='n' & is.na(all.anns$mouse$antisense.dupl.rate)])[1:150])

f = function(l,xlim,col=c('green','gray','lightgray'),relative=TRUE,...){
	l = l[l >= xlim[1] & l <= xlim[2]]
	t = table(l %% 3,ceiling(l / 3)*3)
	if(relative)
		t = sweep(t,2,apply(t,2,sum),'/')
	barplot(t,col=col,border=NA,xlab='exon length (nt)',ylab=ifelse(relative,'fraction of exons','# exons'),space = 0,width=3,xaxt='n',...)
	axis(1)
	legend('topright',fill=col,legend = 0:2,title = 'l %% 3')
	abline(v=27,lty=3)
	abline(h=1/3,lty=3)
}


pdf('figures/paper.figures/2019.04.30/microexons.length.pdf',w=4*4,h=4*7)
par(mfrow=c(7,4),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,2.5),oma=c(0,0,0,1))
for(s in names(anns)){
	print(s)
	a  = all.anns[[s]][all.anns[[s]]$sites=='ad' & all.anns[[s]]$is.ce,]
	af = anns[[s]][anns[[s]]$sites=='ad' & anns[[s]]$is.ce,]
	b = seq(0,max(a$length)+3,by=3)
	xlim=c(0,201)
	h1 = hist(a$length,b,xlim=xlim,m=s,xlab='exon length (nt)',ylab='# exons',col='lightgray',border=NA)
	h2 = hist(a$length[a$type=='ALT'],b,xlim=xlim,m=s,xlab='',ylab='',col='orange',border=NA,add=T)
	h3 = hist(af$length,b,plot=F)
	h4 = hist(a$length[rownames(a) %in% names(devas02[[s]])[devas02[[s]]]],b,plot=F)
	lines(b[-1]-1.5,h2$counts/h1$counts*max(h1$counts),col='red',lwd=3)
	lines(b[-1]-1.5,h4$counts/h3$counts*max(h1$counts),col='magenta',lwd=3)
	abline(v=27,lty=3)
	l=legend('topright',fill=c('lightgray','orange'),legend=c('const.','AS'))
	legend(l$rect$left+l$rect$w,l$rect$top-l$rect$h,col=c('red','magenta'),lwd=3,legend=c('AS/ALL','devAS/AS'),xjust = 1)
	at = seq(0,100,by=20)
	axis(4,at/100*max(h1$counts),at,col='red',col.axis='red')
	mtext('% of exons',4,1.3,col='red')
	
	f(a$length[a$type=='EXN'],xlim=xlim,main='Constitutive exons')
	f(a$length[a$type=='ALT'],xlim=xlim,main='AS exons')
	f(anns[[s]]$length[anns[[s]]$site=='ad'],xlim=xlim,main='AS exons passed filtering')
}
dev.off()

s = 'mouse'
af = anns[[s]][anns[[s]]$sites=='ad' & anns[[s]]$is.ce & anns[[s]]$length < 10,]
table(af$length,rownames(af) %in% names(devas02$mouse)[devas02$mouse])
af[af$length==4,]
plotTissueAgeProile(psi.tsm[[s]]['mou.19170.s11',],meta.tsm)

b = list.files('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/mouse/',pattern = 'Heart',full.names = T)
b = b[grep('.bam$',b)]

r = getReadCoverage(b[1:2],'17',8293415-100,8293418+100,-1,T,min.junc.cov = 0)
abline(v=c(18932498,18932569),col='blue')


age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)
age.dpsi$macaque = cbind(age.dpsi$macaque[,1:5],ovary=NaN,age.dpsi$macaque[,6,drop=FALSE])


ff = function(a,s,l=Inf,st='ad'){a[[s]]$sites %in% st & a[[s]]$length<=l}


pdf('figures/paper.figures/2019.02.21/microexons.adult.tissue-spec.pdf',w=3*4,h=4*7)
par(mfrow=c(7,3),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,2.5),oma=c(0,0,0,1))
for(s in names(anns)){
	tissues=unique(meta$tissue)
	if(s=='macaque')
		tissues = tissues[-6]
	tis = substr(tissues,1,1)
	t=psi.tsm[[s]][ff(anns,s,27),paste(s,tissues,border.stages[[s]][tissues,2])]
	max = factor(apply(t,1,function(x)tis[order(-x)[1]]))
	dim = apply(t,1,function(x){max(x,na.rm=T)-min(x,na.rm=T)})
	barplot(table(max),col=params$tissue.col[tissues],border = NA,main=paste0(s,', all (',nrow(t),')'),ylab='# microexons')
	f=apply(!is.na(t),1,sum)==length(tissues)
	barplot(table(max[f]),col=params$tissue.col[tissues],border = NA,main=paste0('exp. in all tissues (',sum(f),')'),ylab='# microexons')
	f=apply(!is.na(t),1,sum)==length(tissues) & !is.na(dim) & dim > 0.5
	barplot(table(max[f]),col=params$tissue.col[tissues],border = NA,main=paste0('exp. in all tissues and adult-dPSI > 0.5 (',sum(f),')'),ylab='# microexons')
}
dev.off()



# fig 5 ####
phastcons = read.table('/home/mazin/skoltech/projects/evo.devo/processed/ad.phastcons.gz',sep='\t')
phastcons = setNames(strsplit(phastcons[,2],',',TRUE),phastcons[,1])
phastcons = lapply(phastcons,as.numeric)



orth.age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(orth.seg.ad.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(orth.age.dpsi) = rownames(species)
ff = function(a,s,l=Inf,st='ad'){a[[s]]$sites %in% st & a[[s]]$length<=l}

dPSI=0.5
me.cnt=sapply(names(anns),function(s)sum(ff(anns,s,l=27)))
me.cnt=do.call(rbind,lapply(1:7,function(i)me.cnt))
me.cnt.tis = sapply(names(anns),function(s)apply(!is.na(per.tissue.age.qv[[s]][ff(anns,s,l=27),]) & !is.na(age.dpsi[[s]][ff(anns,s,l=27),]),2,sum,na.rm=T))
sgn02u = sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][ff(anns,s),] < 0.05 & age.dpsi[[s]][ff(anns,s),]>  dPSI,2,sum,na.rm=T))
sgn02d = sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][ff(anns,s),] < 0.05 & age.dpsi[[s]][ff(anns,s),]< -dPSI,2,sum,na.rm=T))

sgn02u.me = sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][ff(anns,s,l=27),] < 0.05 & age.dpsi[[s]][ff(anns,s,l=27),]>  dPSI,2,sum,na.rm=T))
sgn02d.me = sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][ff(anns,s,l=27),] < 0.05 & age.dpsi[[s]][ff(anns,s,l=27),]< -dPSI,2,sum,na.rm=T))

# only inclusiond and only brain
# 
micro.timing = lapply(rownames(species)[-2], function(s){
	f = anns[[s]]$sites=='ad' & per.tissue.age.qv[[s]][,'brain']<0.05 & age.dpsi[[s]][,'brain'] > dPSI
	f = !is.na(f) & f
	# dpsi.age = calcdPSIonAge(psi.tsm[[s]][f,grep('brain',colnames(psi.tsm[[s]]))],meta.tsm)
	# max.brain.dpsi = apply(dpsi.age,1,function(x)meta.tsm[colnames(dpsi.age)[which.max(x)],'days'])
	# t = table(micro=anns[[s]]$length[f]<=27,max.brain.dpsi>=meta.tsm[paste(s,'brain',age.al.i[10,s]),'days'])
	t = psi.tsm[[s]][f,paste(s,'brain',c(border.stages[[s]]['brain',1],age.al.i[10,s],border.stages[[s]]['brain',2]))]
	t = table(micro=anns[[s]]$length[f]<=27,before.birth = factor((t[,2]-t[,1])>(t[,3]-t[,2]),levels=c(F,T)))
})
names(micro.timing) = rownames(species)[-2]


pdf(paste0('figures/paper.figures/2019.02.21/5.microexons.devAS.dPSI>',dPSI,'.pdf'),w=12,h=12)
par(mfrow=c(3,3),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
b=barplot(t(sgn02u.me/sgn02u)*100,beside = T,col=rep(params$tissue.col,each=7),ylab='% of microexons',main=paste0('% of microexons in devAS (dPSI > ',dPSI,')'),names.arg=substr(rownames(sgn02u.me),1,1),border=NA)
text(b,t(sgn02u.me/sgn02u)*100+0.2,t(sgn02u.me),srt=90,cex=0.5,adj=c(0,0.5))
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.7)
plotPanelLetter('A')

b=barplot(t(sgn02d.me/sgn02d)*100,beside = T,col=rep(params$tissue.col,each=7),ylab='% of microexons',main=paste0('% of microexons in devAS (dPSI < ',dPSI,')'),names.arg=substr(rownames(sgn02u.me),1,1),ylim=range(0,sgn02u.me/sgn02u*100,na.rm=T),border=NA)
text(b,t(sgn02d.me/sgn02d)*100+0.2,t(sgn02d.me),srt=90,cex=0.5,adj=c(0,0.5))
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.7)
plotPanelLetter('B')


cols = paste0(rep(params$tissue.col,each=7),'22')
barplot(t(me.cnt),beside = T,col=cols,ylab='# of microexons',main='devAS in microexons',names.arg=substr(rownames(sgn02u.me),1,1),border=NA)
cols = paste0(rep(params$tissue.col,each=7),'44')
barplot(t(me.cnt.tis),beside = T,col=cols,ylab='',main='',names.arg=substr(rownames(sgn02u.me),1,1),border=NA,add=T)
cols = rep(params$tissue.col,each=7)
barplot(t(sgn02u.me),beside = T,col=cols,ylab='',main='',names.arg=substr(rownames(sgn02u.me),1,1),border=NA,add=T)
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.7)
legend('topright',fill=c('#00000033','#00000077','#000000'),legend = c('all','expressed',paste0('devAS, dPSI>',dPSI)),border = NA,bty = 'n')
plotPanelLetter('C')


#ft=sapply(micro.timing,function(t)as.numeric(fisher.test(t)[c('estimate','p.value')]))
r=sapply(micro.timing,function(x)x[2,]/(x[1,]+x[2,]))*100
b=barplot(r[2:1,],beside = T,ylab='% of microexons',legend.text = c('before both','after birth'),ylim=c(0,50),main=paste('Inclusion devAS (dPSI>',dPSI,')'))
t = paste0(sapply(micro.timing,function(x)x[2,2:1]),'/',sapply(micro.timing,function(x){x[1,2:1]+x[2,2:1]}))
text(b,r[2:1,],t,srt=90,cex=1,adj=c(-0.1,0.5))
plotPanelLetter('D')

orth.micro = apply(sapply(orth.seg.ad,function(x)x$seg$length<=27),1,sum)==7
anc = apply(sapply(orth.seg.ad,function(x)x$seg$type=='ALT'),1,sum)==7
tis = unique(meta$tissue)[-7]

f = anc &  orth.per.tissue.age.qv$mouse[,'brain']<0.05 & orth.age.dpsi$mouse[,'brain'] > dPSI
f = !is.na(f) & f

dpsi.mb = orth.seg.ad.tsm$mouse[,'mouse brain 9wpb']-orth.seg.ad.tsm$mouse[,'mouse brain 10.5']
orth.macro.eq = names(dpsi.mb) %in% names(equalyseDistrs(dpsi.mb[f & orth.micro],dpsi.mb[f & !orth.micro]))


hist(dpsi.mb[f & orth.micro],col='#FF000080',freq=F,border = NA,xlab='mouse brain dPSI',main='Mouse devAS dPSI distr',ylim=c(0,2.5))
hist(dpsi.mb[f & !orth.micro],col='#0000FF80',add=T,freq=F,border = NA)
hist(dpsi.mb[f & orth.macro.eq],col='#00FF0080',add=T,freq=F,border = NA,xlab='')
legend('topleft',title=paste0('Ancient; mouse brain devAS (dPSI>',dPSI,')'),fill=c('red','green','blue'),
			 legend=paste0(c('microexons','macroexons, same dPSI','macroexons'),' (',c(sum(f & orth.micro),sum(f & orth.macro.eq),sum(f & !orth.micro)),')'))
plotPanelLetter('E')

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
plotPanelLetter('F')
#correlation across brain development
age.al.i. = age.al.i[age.al.i$human %in% meta.tsm$stage[meta.tsm$species=='human' & meta.tsm$tissue=='brain'],]

ma=sapply(lapply(sps[-5],function(s)lineCors(orth.seg.ad.tsm$mouse[f & !orth.micro   ,paste('mouse brain',age.al.i.$mouse)],orth.seg.ad.tsm[[s]][f & !orth.micro   ,paste(s,'brain',age.al.i.[,s])],u='p')),function(x){x=x[!is.na(x)];sd=sd(x)/sqrt(length(x));r=mean(x);c(r,r-2*sd,r+2*sd)})
mi=sapply(lapply(sps[-5],function(s)lineCors(orth.seg.ad.tsm$mouse[f &  orth.micro   ,paste('mouse brain',age.al.i.$mouse)],orth.seg.ad.tsm[[s]][f &  orth.micro   ,paste(s,'brain',age.al.i.[,s])],u='p')),function(x){x=x[!is.na(x)];sd=sd(x)/sqrt(length(x));r=mean(x);c(r,r-2*sd,r+2*sd)})
me=sapply(lapply(sps[-5],function(s)lineCors(orth.seg.ad.tsm$mouse[f &  orth.macro.eq,paste('mouse brain',age.al.i.$mouse)],orth.seg.ad.tsm[[s]][f &  orth.macro.eq,paste(s,'brain',age.al.i.[,s])],u='p')),function(x){x=x[!is.na(x)];sd=sd(x)/sqrt(length(x));r=mean(x);c(r,r-2*sd,r+2*sd)})

plotArea(1:ncol(ma),t(ma),col='blue',t='b',ylim=range(mi,ma),lwd=3,new = T,xlab='',ylab='mean Pearson corr. tp mouse',main='Corelation across brain development',xaxt='n')
plotArea(1:ncol(ma),t(mi),col='red',t='b',ylim=range(mi,ma),lwd=3)
plotArea(1:ncol(me),t(me),col='green',t='b',ylim=range(mi,ma),lwd=3)
axis(1,1:5,sps)
plotPanelLetter('G')

ma=apply(sapply(orth.age.dpsi,function(x){x[f & !orth.micro,'brain']})>0.2,2,function(x){x=x[!is.na(x)];my.binom.test(sum(x),sum(!x))})[,sps[-5]]
mi=apply(sapply(orth.age.dpsi,function(x){x[f &  orth.micro,'brain']})>0.2,2,function(x){x=x[!is.na(x)];my.binom.test(sum(x),sum(!x))})[,sps[-5]]
me=apply(sapply(orth.age.dpsi,function(x){x[f &  orth.macro.eq,'brain']})>0.2,2,function(x){x=x[!is.na(x)];my.binom.test(sum(x),sum(!x))})[,sps[-5]]


plotArea(1:ncol(ma),t(ma)*100,col='blue',t='b',ylim=range(mi,ma)*100,lwd=3,new = T,xlab='',ylab='% of exons with dPSI > 0.2',main='Conservation of devAS',xaxt='n')
plotArea(1:ncol(ma),t(mi)*100,col='red',t='b',ylim=range(mi,ma),lwd=3)
plotArea(1:ncol(me),t(me)*100,col='green',t='b',ylim=range(mi,ma),lwd=3)
axis(1,1:5,sps)
plotPanelLetter('H')

n3 = rbind(mi=table(orth.seg.ad$mouse$seg$length[f &  orth.micro] %% 3 == 0),
					 me=table(orth.seg.ad$mouse$seg$length[f &  orth.macro.eq] %% 3 == 0),
					 ma=table(orth.seg.ad$mouse$seg$length[f & !orth.micro] %% 3 == 0))
n3.s = apply(n3[,2:1],1,my.binom.test)
plot(1:3,n3.s[1,],xaxt='n',xlab='Exon set',ylab='proportion of 3N',ylim=range(n3.s),col=c('red','green','blue'),pch=19,cex=1.4,main='mouse')
segments(1:3,n3.s[2,],1:3,n3.s[3,],col=c('red','green','blue'),lwd=3)
plotPanelLetter('I')
dev.off()

t=table(orth.seg.ad$human$seg$length < 28,apply(sapply(orth.seg.ad,function(x)x$seg$type=='ALT'),1,sum))
f = t[2,]/apply(t,2,sum)
b=barplot(f,ylab='proportion of microexons',xlab='# of AS species',ylim=c(0,max(f)*1.2))
text(b,f,t[2,],srt=90,adj=c(-0.1,0.5))
# phastcons
plotPhast = function(phast,ss,...){
	ph = lapply(ss,function(s)t(sapply(phast[intersect(s,names(phast))],function(x)c(x[1:200],rep(NA,30),x[length(x)-(199:0)]))))
	ph.p = lapply(ph,function(p)cbind(apply(p,2,mean,na.rm=T),apply(p,2,sd,na.rm=T)/sqrt(nrow(p))))
	x = 1:nrow(ph.p$mi)
	
	plotArea(x,ph.p$mi,col='red',xaxt='n',new = T,ylab='human mean phastcons',...)
	plotArea(x,ph.p$me,col='green',xaxt='n',new = F)
	plotArea(x,ph.p$ma,col='blue',xaxt='n',new = F)
	plotArea(x,ph.p$unsign,col='gray',xaxt='n',new = F)
	rect(200,0,230,1,col = 'white',border = NA)
	axis(1,c(1,200,230,430),c('-200','acc','don','+200'),las=2)
	legend('topright',fill=c('red','green','blue','gray'),legend = paste0(c('microexons','macroexons, same dPSI','macroexons','!devAS'),' (',sapply(ss,length),')'))
}

dPSI = 0.2

orth.micro = apply(sapply(orth.seg.ad,function(x)x$seg$length<=27),1,sum)==7
anc = apply(sapply(orth.seg.ad,function(x)x$seg$type=='ALT'),1,sum)==7
tis = unique(meta$tissue)[-7]



f = anc &  orth.per.tissue.age.qv$human[,'brain']<0.05 & orth.age.dpsi$human[,'brain'] > dPSI
f = !is.na(f) & f

dpsi.mb = orth.seg.ad.tsm$human[,'human brain 9wpb']-orth.seg.ad.tsm$human[,'human brain 10.5']
orth.macro.eq = names(dpsi.mb) %in% names(equalyseDistrs(dpsi.mb[f & orth.micro],dpsi.mb[f & !orth.micro]))

sids = list(mi=rownames(orth.seg.ad$human$seg)[f &  orth.micro],
						me=rownames(orth.seg.ad$human$seg)[f &  orth.macro.eq],
						ma=rownames(orth.seg.ad$human$seg)[f & !orth.micro],
						unsign=rownames(orth.seg.ad$human$seg)[orth.per.tissue.age.qv$human[,'brain']>0.05])

f = anns$human$sites=='ad' & per.tissue.age.qv$human[,'brain'] < 0.05 & age.dpsi$human[,'brain'] > dPSI
f[is.na(f)] = F
sids.h = list(mi=rownames(anns$human)[f & anns$human$length < 28],
						ma=rownames(anns$human)[f & anns$human$length > 27],
						unsign = rownames(anns$human)[anns$human$sites=='ad' & per.tissue.age.qv$human[,'brain'] > 0.05 & anns$human$length > 27])

sids.h$me = names(equalyseDistrs(age.dpsi$human[sids.h$mi,'brain'],age.dpsi$human[f,'brain']))
sids.h = sids.h[c('mi','me','ma','unsign')]
sapply(sids.h,length)

pdf(paste0('figures/paper.figures/2019.02.21/microexons.phastcons.dPSI>',dPSI,'.pdf'),w=10,h=5)
par(mfrow=c(1,2),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
plotPhast(phastcons,sids,main=paste0('Ancient exons, human devAS (dPSI > ',dPSI,')'),ylim=c(0.15,1.20))
plotPhast(phastcons,sids.h,main=paste0('Human devAS (dPSI > ',dPSI,')'),ylim=c(0.15,1.20))
dev.off()

hist(age.dpsi$human[sids$mi,'brain'],col='#FF000080',freq=F,border = NA,xlab='mouse brain dPSI',main='Mouse devAS dPSI distr',ylim=c(0,2.5))
hist(age.dpsi$human[sids$ma,'brain'],col='#0000FF80',add=T,freq=F,border = NA)
hist(age.dpsi$human[sids$me,'brain'],col='#00FF0080',add=T,freq=F,border = NA,xlab='')


# evo-ages of microexons
olen = sapply(orth.seg.ad,function(x)x$seg$length)
alt = sapply(orth.seg.ad,function(x)x$seg$type)=='ALT'

#ce = sapply(names(orth.seg.ad),function(s)all.anns[[s]][rownames(orth.seg.ad[[s]]$seg),'is.ce'])
#the are all CE)

born.exn.sajr = readRDS('Rdata/born.exn.sajr.Rdata')
exon.birth.one = readRDS('Rdata/exon.birth.one.Rdata')



pdf('figures/paper.figures/2019.04.30/microexons.evo.age.pdf',w=9,h=6)
par(mfrow=c(2,3),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,0,1))
t=table(apply(alt,1,sum),apply(olen,1,max)<=27)
barplot(t[,2]/apply(t,1,sum),xlab='# of ALT species',ylab='fraction of microexons')
barplot(t[,2],xlab='# of ALT species',ylab='# of exons',main='Microexons')
barplot(t[,1],xlab='# of ALT species',ylab='# of exons',main='Macroexons')

x=t(sapply(exon.birth.one,function(x)c(sp.cnt=sum(!is.na(x$length)),len=max(x$length,na.rm=T))))
t=table(x[,1],x[,2]<=27)


barplot(t[,2]/apply(t,1,sum),xlab='# of species with exon (newborn)',ylab='fraction of microexons')
barplot(t[,2],xlab='# of species with exon (newborn)',ylab='# of exons',main='Microexons')
barplot(t[,1],xlab='# of species with exon (newborn)',ylab='# of exons',main='Macroexons')
dev.off()

# pseudomicroexons ######
adj.segments = list()
for(sp in rownames(species)){
	print(sp)
	h = readRDS(paste0('Rdata/',sp,'.as.u.all.Rdata'))
	s = h$seg[order(h$seg$gene_id,h$seg$start),]
	l = nrow(s)
	fu = s$gene_id[-1] == s$gene_id[-l] & s$start[-1] == (s$stop[-l]+1)
	fd = s$gene_id[-l] == s$gene_id[-1] & s$stop[-l] == (s$start[-1]-1)
	table(c(F,fu),c(fd,F))
	
	s$up.ajd.sid = s$dw.ajd.sid = NA
	s$up.ajd.sid[c(F,fu)] = rownames(s)[c(fu,F)]
	s$dw.ajd.sid[c(fd,F)] = rownames(s)[c(F,fd)]
	f = s$strand == -1
	t = s$dw.ajd.sid[f]
	s$dw.ajd.sid[f] = s$up.ajd.sid[f]
	s$up.ajd.sid[f] = t
	
	s = s[rownames(h$seg),]
	s$dw.ajd.med.psi.rate = s$up.ajd.med.psi.rate = NA
	f = !is.na(s$up.ajd.sid)
	s$up.ajd.med.psi.rate[f] = apply(h$ir[rownames(s)[f],]/h$ir[s$up.ajd.sid[f],],1,median,na.rm=T)
	f = !is.na(s$dw.ajd.sid)
	s$dw.ajd.med.psi.rate[f] = apply(h$ir[rownames(s)[f],]/h$ir[s$dw.ajd.sid[f],],1,median,na.rm=T)
	adj.segments[[sp]] = s[,c('up.ajd.sid','dw.ajd.sid','up.ajd.med.psi.rate','dw.ajd.med.psi.rate')]
}
#saveRDS(adj.segments,'Rdata/adj.segments.sids-n-psi.Rdata')

