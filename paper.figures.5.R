options(stringsAsFactors = FALSE)
source('code/r.functions/load.all.data.F.R')
source('~/skoltech/r.code/util.R')
source('code/r.functions/paper.figures.5.F.R')
library(png)
library(reshape)
library(SAJR)
library(ape)
library("extrafont")
#font_import()
#fonts()

# load data #####
species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
anns = readRDS('Rdata/anns.Rdata')
#all.anns = readRDS('Rdata/all.anns.Rdata')
border.stages = readRDS('Rdata/border.stages.Rdata')
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')
orth.seg.ad.tsm = readRDS('Rdata/orth.seg.ad.tsm.Rdata')
seg2ens = readRDS('Rdata/seg2ens.Rdata')
orth.per.tissue.age.qv = readRDS('Rdata/orth.per.tissue.age.qv.Rdata')
age.al.i = read.table('input/margarida.age.alignment.final.2017.05.17.improved.tab',sep='\t',header = T)[,c(rownames(species),'to.remove')]
gene.descrs = readRDS('Rdata/ens.gene.descr.Rdata')
#ens.ge = readRDS('Rdata/ens.ge.Rdata')
#ens.ge.cod = readRDS('Rdata/ens.ge.cod.Rdata')

params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)
lab.cex=1.5



# precalculate #####
age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)
age.dpsi$macaque = cbind(age.dpsi$macaque[,1:5],ovary=NaN,age.dpsi$macaque[,6,drop=FALSE])

alt.sp = readRDS('Rdata/paper.figures/alt.sp.Rdata')
alt.sp = alt.sp[alt.sp!='']

orth.age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(orth.seg.ad.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(orth.age.dpsi) = rownames(species)
orth.age.dpsi$macaque = cbind(orth.age.dpsi$macaque[,1:5],ovary=NA,orth.age.dpsi$macaque[,6,drop=F])

# look for examples #####
# alt.sp = readRDS('Rdata/paper.figures/alt.sp.Rdata')
# o = readRDS('Rdata/orth.seg.ad.all.Rdata')
# o = o$human$seg$north[alt.sp!='']
# table(o)
# saveRDS(o,'Rdata/paper.figures/orth.seg.ad.all.north.Rdata')
hgmd = read.table('input/hgmd/tosend/2017_1_HGMD_ALL_Variants.csv',sep='\t',row.names=1,header=T,quote='',comment.char = '')
colnames(hgmd)
ens.orth = readRDS('Rdata/paper.figures/orth.seg.ad.all.north.Rdata')
table(ens.orth)
f = ens.orth==7 & alt.sp=='hqmrboc' & apply(is.na(orth.age.dpsi$mouse) | is.na(orth.per.tissue.age.qv$mouse),1,sum)==0 & abs(orth.age.dpsi$mouse[,'brain']) > 0.7 & orth.per.tissue.age.qv$mouse[,'brain'] < 0.05
table(f)
sids = rownames(orth.age.dpsi$mouse)[f]
gd = gene.descrs$mouse[unique(unlist(seg2ens$mouse[sids])),]
gd[toupper(gd$gene.name) %in% toupper(hgmd$gene),]
hgmd[hgmd$gene=='DST',]

# sid = sids[sapply(seg2ens$mouse[sids],function(x)'ENSMUSG00000032479' %in% x)] #MAP4 - OK
# sid = sids[sapply(seg2ens$mouse[sids],function(x)'ENSMUSG00000042605' %in% x)] #Atxn2 - not so nice
# sid = sids[sapply(seg2ens$mouse[sids],function(x)'ENSMUSG00000034826' %in% x)] #Nup54 - looks nice
sid = sids[sapply(seg2ens$mouse[sids],function(x)'ENSMUSG00000000881' %in% x)] #Dlg3 - looks nice + Intellectual disability (HGMD)
# sid = sids[sapply(seg2ens$mouse[sids],function(x)'ENSMUSG00000022377' %in% x)] #Asap1 - looks nice + Schizophrenia? (HGMD)
# sid = sids[sapply(seg2ens$mouse[sids],function(x)'ENSMUSG00000039844' %in% x)] #Rapgef1 - looks nice + Intellectual disability? (HGMD)

sid = which(rownames(orth.age.dpsi$mouse) == sid)
anns$mouse[rownames(orth.age.dpsi$mouse)[sid],]
anns$human[rownames(orth.age.dpsi$human)[sid],]

sid = which(rownames(orth.age.dpsi$mouse) == 'mou.49187.s24')
par(mfrow=c(3,3),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(1.5,2.5,1.5,0),oma=c(0,0,0,0))
for(s in rownames(species))
	plotTissueAgeProile(orth.seg.ad.tsm[[s]][sid,],meta.tsm,age.axis = 'rank',main=s)
s = 'human'
id = rownames(orth.age.dpsi[[s]])[sid]
bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$days>25*365 & meta$days<40*365 & meta$tissue=='brain'],'.bam')
b = getReadCoverage(bams,anns[[s]][id,'chr_id'],anns[[s]][id,'start'] - 1000,anns[[s]][id,'stop']+4000,-anns[[s]][id,'strand'])
bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$days>25*365 & meta$days<40*365 & meta$tissue=='testis'],'.bam')
t = getReadCoverage(bams,anns[[s]][id,'chr_id'],anns[[s]][id,'start'] - 1000,anns[[s]][id,'stop']+4000,-anns[[s]][id,'strand'])
bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$days>25*365 & meta$days<40*365 & meta$tissue=='heart'],'.bam')
h = getReadCoverage(bams,anns[[s]][id,'chr_id'],anns[[s]][id,'start'] - 1000,anns[[s]][id,'stop']+4000,-anns[[s]][id,'strand'])

par(mfrow=c(3,1),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(1.5,2.5,1.5,0),oma=c(0,0,0,0))
plotReadCov(b,reverse = (anns[[s]][id,'strand'] == -1),min.junc.cov = 2,bty='n',xlab=paste0('Chr ',anns[[s]][id,'chr_id']),ylab='Reads')
plotReadCov(t,reverse = (anns[[s]][id,'strand'] == -1),min.junc.cov = 2,bty='n',xlab=paste0('Chr ',anns[[s]][id,'chr_id']),ylab='Reads')
plotReadCov(h,reverse = (anns[[s]][id,'strand'] == -1),min.junc.cov = 2,bty='n',xlab=paste0('Chr ',anns[[s]][id,'chr_id']),ylab='Reads')
plot(orth.seg.ad.tsm[[s]]['hum.68349.s28',],orth.seg.ad.tsm[[s]]['hum.68349.s30',],col=meta.tsm[colnames(orth.seg.ad.tsm[[s]]),'col'],pch=19)



# bams = paste0('processed/mapping/hisat2.s/mouse/',meta$fname[meta$species=='mouse' & meta$stage=='13.5' & meta$tissue=='brain'],'.bam')
# dlg3.mdata = list()
# dlg3.mdata$ann = all.anns$mouse[all.anns$mouse$gene_id=='mou.49187',]
# dlg3.mdata$cov = getReadCoverage(bams,dlg3.mdata$ann$chr_id[1],min(dlg3.mdata$ann$start),max(dlg3.mdata$ann$stop),-dlg3.mdata$ann$strand[1])
# dlg3.mdata$zoom.coor = c(100807180,100812776)
#saveRDS(dlg3.mdata,'Rdata/paper.figures/dlg3.mdata.Rdata')



# 1 #####
# _prepare ####
phyl.tree=read.tree(text = '((((human:29.44,macaque:29.44):60.38,((mouse:20.89,rat:20.89):61.25,rabbit:82.14):7.68):68.77,opossum:158.6):153.4,chicken:312);')
phyl.tree$edge.length = rev(phyl.tree$edge.length)
phyl.tree$edge = phyl.tree$edge[nrow(phyl.tree$edge):1,]

plot(phyl.tree)
at = c(300,180,90,25)
axis(1,312-at,at,las=2)
title(xlab='Million years')
# B
astypes = c(CE='ad',AA='aa',AD='dd',RI='da')
tested.stat = lapply(setNames(astypes,astypes),function(ss)sapply(names(per.tissue.age.qv),function(s)apply(!is.na(per.tissue.age.qv[[s]][anns[[s]]$sites==ss,]),2,sum)))

# C
sgn02.stat = lapply(setNames(astypes,astypes),function(ss) {
	sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][anns[[s]]$sites == ss,] < 0.05 & abs(age.dpsi[[s]][anns[[s]]$sites == ss,])>0.2,2,sum,na.rm=T))/tested.stat[[ss]]*100})

astypes.pchs=c(ad=19,aa=2,dd=6,da=13)

# _no of tissue-specific devAS ######
sgn = list()
for(s in rownames(species)){
	sgn[[s]] = abs(age.dpsi[[s]]) > 0.5 & per.tissue.age.qv[[s]] < 0.05
	sgn[[s]][is.na(sgn[[s]])] = FALSE
}

sort(round(sapply(rownames(species),function(s){
	t = table(apply(sgn[[s]][anns[[s]]$sites=='ad',-2],1,sum))
	t['1']/sum(t[-1])
	})*100,1))



# dPSI > 0.2, ad
# chicken.1     rat.1   mouse.1  rabbit.1 opossum.1   human.1 macaque.1 
# 68.0      69.3      70.0      71.7      76.1      77.4      88.3
# dPSI > 0.5, ad
# rat.1 chicken.1   mouse.1  rabbit.1   human.1 opossum.1 macaque.1
# 85.6      85.7      88.7      89.0      91.4      92.0      95.7

#_plot #####

pdf('figures/paper.figures/5/1/1.pdf',w=7.2,h=7.2/3,family='Arial',)
layout(matrix(1:4,ncol=4,byrow = T),widths = c(2,3.25,0.5,3.25))
par(tck=-0.01,mgp=c(1.3,0.2,0),mar=c(1.5,2.5,1.5,0),oma=c(0,0,0,0))
plot.new()
plotPanelLetter('A',lab.cex)
# order by species/tissue
plotAsEventCount(tested.stat,astypes.pchs,by.tissue = T,ylab='# of detected events',main='Detected AS',bty='n')
plotPanelLetter('B',lab.cex)
par(mar=c(2.5,0,1.5,0))
y=plotASTypes(astypes.pchs)
legend(0,y-5,pch=19,col=params$tissue.col,legend = names(params$tissue.col),bty='n',cex = 0.5,xjust=0,yjust=1)
rect(0,95,100,16,xpd=T,col=NA,border='gray')
par(mar=c(1.5,2.5,1.5,0))
plotAsEventCount(sgn02.stat,astypes.pchs,by.tissue = T,ylab='% of devAS',main='DevAS',bty='n')
plotPanelLetter('C',lab.cex)
dev.off()

# 2 #####
# _prepare #####
# cors = list()
# for(s in rownames(species)){
# 	print(s)
# 	d = readRDS(paste0('Rdata/',s,'.as.u.filtered.Rdata'))
# 	d = d$ir[d$seg$sites=='ad',]
# 	gc()
# 	cors[[paste0(species[s,'short'],'')]] = cor(d,u='p',m='p')
# 	cors[[paste0(species[s,'short'],'.tsm')]] = cor(psi.tsm[[s]][anns[[s]]$sites=='ad',],u='p',m='p')
# }
# saveRDS(cors,'Rdata/ad.all.cor.Rdata')

# just plot all I have
# cors = readRDS('Rdata/ad.all.cor.Rdata')
# pdf('figures/MDSs/one-species.MDS.pdf',w=6,h=3*7)
# pows = c(0.1,0.5,1,2,4)
# par(mfrow=c(7,2),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(1.5,2.5,1.5,0),oma=c(0,0,0,0))
# for(p in pows){
# 	mds = lapply(cors,function(x)cmdscale((1-x)^p,k=2))
# 	for(s in rownames(species)){
# 		id = species[s,'short']
# 		m = meta[rownames(mds[[id]]),]
# 		plot(mds[[id]],pch=19,col=m$col,cex=m$cex,xlab='Dim 1',ylab='Dim 2',main=paste0(s,', all samples; pow=',p))
# 		
# 		id = paste0(species[s,'short'],'.tsm')
# 		m = meta.tsm[rownames(mds[[id]]),]
# 		plot(mds[[id]],pch=19,col=m$col,cex=m$cex,xlab='Dim 1',ylab='Dim 2',main=paste0(s,', stage-tissue means; pow=',p))
# 	}
# }
# dev.off()
cors = readRDS('Rdata/ad.all.cor.Rdata')
sp.ce = 'mouse'
mouse.ge = readRDS('Rdata/ens.ge.marg.tsm.Rdata')[[sp.ce]]
mouse.ge = mouse.ge + min(mouse.ge[mouse.ge!=0],na.rm=T)
mouse.ge = log2(mouse.ge)

use.mean.embryo = FALSE
#sd = apply(psi.tsm[[sp.ce]],1,sd)
#psi.cor2embryo = caclCor2Embryo(psi.tsm[[sp.ce]][anns[[sp.ce]]$sites=='ad' & sd > 0.05,],meta.tsm,cor.m = 'sp',use.mean.embryo=use.mean.embryo)
psi.cor2embryo = caclCor2Embryo(psi.tsm[[sp.ce]][anns[[sp.ce]]$sites=='ad',],meta.tsm,cor.m = 'p',use.mean.embryo=use.mean.embryo)
ge.cor2embryo = caclCor2Embryo(mouse.ge,meta.tsm,cor.m = 'p',use.mean.embryo=use.mean.embryo)

# tau
human.tau = read.csv('input/gene.info.from.marg/Human.Indexes.All.csv')
rownames(human.tau) = human.tau[,1]
table(substr(human.tau$Human_ID,1,1),human.tau$Human_ID %in% rownames(ens.ge$human$gene))
human.tau$biotype = setNames(ens.ge$human$gene$biotype,rownames(ens.ge$human$gene))[human.tau$Human_ID]

psi.tsm.ad = lapply(names(psi.tsm),function(s)psi.tsm[[s]][anns[[s]]$sites=='ad',])
names(psi.tsm.ad) = names(psi.tsm)
ts = unique(meta$tissue)
human.tissue.gene.dpsi=lapply(setNames(ts,ts), function(t){print(t);t(getsPSIbyEnsID(psi.tsm.ad,border.stages,t,seg2ens,'human',use.random = F))})


sids = sgn$human[anns$human$sites=='ad','brain']
z = seg2ens$human[names(sids)[sids]]
be = human.tau$Human_ID[!is.na(human.tau$TissueTau) & human.tau$TissueTau<0.5]
mean(sapply(z,function(x)sum(x %in% be)>0)) #81.4% of brain devAS is in tau < 0.5
n = seg2ens$human[names(sids)[!sids]]
mean(sapply(n,function(x)sum(x %in% be)>0)) #77.1%

# _plot #####

pdf('figures/paper.figures/5/2/2.pdf',w=7.2,h=7.2/3*2,family='Arial')
#jpeg('figures/paper.figures/5/2/2.jpg',units = 'in',w=7.2,h=7.2/3*2,quality = 100,res = 600,family='Arial')
layout(matrix(1:6,ncol=3,byrow = F),widths = c(10,10,10))
par(tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,0))
# mouse MDS
mds = cmdscale(1-cors$m,k=2)
m = meta[rownames(mds),]
plot(mds,pch=19,col=m$col,cex=m$cex,xlab='Dim 1',ylab='Dim 2',main='Mouse',bty='n')
plotPNG("figures/paper.figures/5/icons/mouse.png",0.7,0.85,0.15)

x = rep(-0.16+cumsum(c(0,0.035,0.065,0.035)),times=2)[1:7]
y = rep(c(0.13,0.115),each=4)[1:7]
points(x,y,pch=19,col=params$tissue.col,xpd=T)
text(x+0.005,y,names(params$tissue.col),adj = c(0,0.5),xpd=T,cex=0.7)

y = 0.09
points(c(-0.16,-0.148,-0.136),rep(y,3),pch=19,cex = c(0.4,0.8,1.2))
text(-0.136,y,"Early to late",adj=c(-0.15,0.5),cex=0.7)
plotPanelLetter('A',lab.cex)

# oposuum MDS
mds = cmdscale(1-cors$o,k=2)
m = meta[rownames(mds),]
plot(mds,pch=19,col=m$col,cex=m$cex,xlab='Dim 1',ylab='Dim 2',main='Opossum',bty='n')
plotPNG("figures/paper.figures/5/icons/opossum.png",0.8,0.85,0.2)

# divergense to embrio
tiss = unique(meta$tissue)[c(1,3,5:7)]
#par(mar=c(5.5,2.5,1.5,0),mgp=c(1.3,0.2,0))
plotCor2Embryo(psi.cor2embryo[tiss,,],main='Splicing',ylab='Pearson cor. to embryo, PSI',lwd=3,area.transp=0.1,xlab='',bty='n')
plotPNG("figures/paper.figures/5/icons/mouse.png",0.8,0.85,0.15)
plotPanelLetter('B',lab.cex)
plotCor2Embryo(ge.cor2embryo[tiss,,],main='Gene expression',ylab='Pearson cor. to embryo',lwd=3,area.transp=0.1,xlab='',bty='n')
plotPNG("figures/paper.figures/5/icons/mouse.png",0.8,0.85,0.15)
plotPanelLetter('C',lab.cex)

hc = !is.na(human.tau$biotype) & human.tau$biotype=='protein_coding'
plotAsInTauDistr(human.tau[hc,],human.tissue.gene.dpsi$brain,0.2,main='Brain',border=NA)
plotPNG("figures/paper.figures/5/icons/human.png",0.8,0.85,0.15)
legend(4,3000,fill=c('orange','lightgray','darkgray'),legend=c('devAS','AS','non-AS'),bty='n')
plotPanelLetter('D',lab.cex)
# plotAsInTauDistr(human.tau[hc,],human.tissue.gene.dpsi$heart,0.2,main='Heart')
# plotPanelLetter('E',lab.cex)
plotAsInTauDistr(human.tau[hc,],human.tissue.gene.dpsi$liver,0.2,main='Liver',border=NA)
plotPNG("figures/paper.figures/5/icons/human.png",0.8,0.85,0.15)
plotPanelLetter('E',lab.cex)
dev.off()

# 3 #####
# _prepare ####
library(RColorBrewer)
orth.mds2 = readRDS('Rdata/paper.figures/orth.mds2.Rdata')
dlg3.mdata = readRDS('Rdata/paper.figures/dlg3.mdata.Rdata')
# devAS cons
#pdf('figures/paper.figures/5/3/speciesDevAscorr.pdf',w=21,h=15,family='Arial')
# jpeg('figures/paper.figures/5/3/speciesDevAscor.jpg',units = 'in',w=21,h=15,quality = 100,res = 200)
# par(mfrow=c(5,7),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,0))
# s1 = 'mouse'
# for(s2 in rownames(species)[c(4,5,1,6,7)])
# 	for(t in unique(meta$tissue))
# 		plot2Sign(orth.per.tissue.age.qv[[s1]][,t]<0.05,orth.per.tissue.age.qv[[s2]][,t]<0.05,orth.age.dpsi[[s1]][,t],orth.age.dpsi[[s2]][,t],species[s1,'short'],species[s2,'short'],xlab = s1,ylab = s2,main=t,range=1)
# dev.off()
# 
# 
# dev.as.cons = list()
# for(s1 in rownames(species))
# 	for(s2 in setdiff(rownames(species),s1))
# 		for(t in unique(meta$tissue))
# 			for(n in 1:7){
# 				dev.as.cons[[length(dev.as.cons)+1]] = cbind(getDevASCons(s1,s2,t,orth.per.tissue.age.qv,orth.age.dpsi,nchar(alt.sp)==n),sp.cnt=n)
# 			}
# dev.as.cons = do.call(rbind,dev.as.cons)

dev.as.cons.all = list()
sps = c('rat','rabbit','human','opossum','chicken')
f = alt.sp == 'hqmrboc'
f =T
for(s2 in sps)
	for(t in unique(meta$tissue)){
		dev.as.cons.all[[length(dev.as.cons.all)+1]] = getDevASCons('mouse',s2,t,orth.per.tissue.age.qv,orth.age.dpsi,f)
	}

dev.as.cons.all = do.call(rbind,dev.as.cons.all)

# phastcons
age.ad.ph. = readRDS('Rdata/paper.figures/age.ad.ph.50nt.Rdata')
age.ad.ph = age.ad.ph.[1:2]
for(t in unique(meta$tissue))
	age.ad.ph[[t]] = c(age.ad.ph.[[paste(t,'up')]],age.ad.ph.[[paste(t,'dw')]])

# 3N
s='human'
ad.3n = sapply(unique(meta$tissue),function(t){
	f = anns[[s]]$sites=='ad' & per.tissue.age.qv[[s]][,t]<0.05 & abs(age.dpsi[[s]][,t])>0.2 & anns[[s]]$cod=='c'
	f[is.na(f)] = FALSE
	s = sum(anns[[s]][f,'length'] %% 3 == 0);
	t = sum(f)
	c(len3=s,total=t,freq=s/t,conf=binom.test(s,t)$conf.int)
})

sgn = anns[[s]]$sites=='ad' & apply(per.tissue.age.qv[[s]]<0.05 & abs(age.dpsi[[s]])>0.2,1,sum)>0 & anns[[s]]$cod=='c'


exn.len3.freq = c(0.4001031,0.3903307,0.4099356)#my.binom.test(table(all.anns[[s]][all.anns[[s]]$cod=='c' & all.anns[[s]]$sites=='ad' & all.anns[[s]]$type=='EXN' & all.anns[[s]]$gene_id %in% anns[[s]]$gene_id[sgn],'length'] %% 3 == 0)[c('TRUE','FALSE')])
ndevAS.len3.freq = c(0.3834853,0.3619688,0.4053485)#my.binom.test(table(anns[[s]]$length[anns[[s]]$cod=='c' & !sgn & anns[[s]]$sites=='ad' &  anns[[s]]$gene_id %in% anns[[s]]$gene_id[sgn]]  %% 3 == 0)[c('TRUE','FALSE')])

# exn.len3.freq = my.binom.test(table(all.anns[[s]][all.anns[[s]]$sites=='ad' & all.anns[[s]]$type=='EXN' & all.anns[[s]]$gene_id %in% anns[[s]]$gene_id[sgn],'length'] %% 3 == 0)[c('TRUE','FALSE')])
# ndevAS.len3.freq = my.binom.test(table(anns[[s]]$length[!sgn & anns[[s]]$sites=='ad']  %% 3 == 0)[c('TRUE','FALSE')])
ad.3n = t(cbind(cnst=exn.len3.freq,'non-devAS'=ndevAS.len3.freq,ad.3n[-1:-2,]))
ad.3n = ad.3n[,c(2,1,3)]

# IUPR
library(ontologyIndex)
o = get_OBO('processed/exon.onthology/exont.obo')
seg2exont = readRDS('Rdata/seg2exont.Rdata')

iupr = names(seg2exont)[sapply(seg2exont,function(x)'EXONT:000074' %in% x)]

hdevas02 = abs(age.dpsi$human) > 0.2 & per.tissue.age.qv$human<0.05
hdevas02[is.na(hdevas02)] = FALSE
f = anns$human$sites=='ad' & anns$human$cod!='n' & rownames(anns$human) %in% names(seg2exont)
idr02 = apply(hdevas02,2,function(x)my.binom.test(table(rownames(hdevas02)[f & x] %in% iupr)[c('TRUE','FALSE')]))


sgn = anns[[s]]$sites=='ad' & apply(hdevas02,1,sum)>0 & anns[[s]]$cod!='n'
alt = rownames(anns$human)[anns[[s]]$cod!='n' & anns[[s]]$sites=='ad' & !sgn &  anns[[s]]$gene_id %in% anns[[s]]$gene_id[sgn] & rownames(anns$human) %in% names(seg2exont)]

#cnst = rownames(all.anns$human)[all.anns[[s]]$cod!='n' & all.anns[[s]]$sites=='ad' & all.anns[[s]]$type=='EXN' & all.anns[[s]]$gene_id %in% anns[[s]]$gene_id[sgn] & rownames(all.anns$human) %in% names(seg2exont)]
cnst = c(0.3685637,0.3623576,0.3748027)#my.binom.test(table(cnst %in% iupr)[c('TRUE','FALSE')])

idr02 = t(cbind(cnst=cnst,"non-devAS"=my.binom.test(table(alt %in% iupr)[c('TRUE','FALSE')]),idr02))
idr02 = idr02[,c(2,1,3)]

# _plot #####
pdf('figures/paper.figures/5/3/3.pdf',w=7.2,h=7.2/3*2.5,family='Arial')
#jpeg('figures/paper.figures/5/3/3.jpg',units = 'in',w=7.2,h=7.2/3*2.5,quality = 100,res = 600,family='Arial')
layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5),ncol=6,byrow = T),heights = c(3,2))
par(tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
m = meta[rownames(orth.mds2$sp7),]
plot(orth.mds2$sp7,xlab='Dim 1',ylab='Dim 2',col=params$tissue.col[m$tissue],pch=params$species.pch[m$species],cex=m$cex,bty='n')
# z=legend(-0.6,0.7,pch=19,col=params$tissue.col,legend=names(params$tissue.col),ncol = 4,xjust = 0,yjust = 1,bty = 'n')
# z=legend(-0.6,z$rect$top-z$rect$h,pch=params$species.pch,legend=names(params$species.pch),ncol = 4,xjust = 0,yjust = 1,bty = 'n')
# y = z$rect$top-z$rect$h*1.3
x = c(-.6,-.6,-.34,-.34,-.18,-.18,-.02)
y = rep(c(0.54,0.5),times=4)[1:7]
points(x,y,pch=19,col=params$tissue.col,xpd=T)
text(x+0.02,y,names(params$tissue.col),adj = c(0,0.5),xpd=T)

x = c(-.6,-.6,-.38,-.38,-.2,-.2,-.0)
y = y - 0.11
points(x,y,pch=params$species.pch)
text(x+0.02,y,names(params$species.pch),adj = c(0,0.5))

y = y[1]-0.11
points(c(-0.6,-0.55,-0.5),rep(y,3),pch=19,cex = c(0.4,1.15,1.9))
text(-0.5,y,"Early to late",adj=c(-0.2,0.5))
plotPanelLetter('A',lab.cex)
plot.new()
# x = match(dev.as.cons.all$s2,sps)
# plot(x,dev.as.cons.all$also.sgn,pch=params$species.pch[dev.as.cons.all$s2],col=params$tissue.col[dev.as.cons.all$t],bty='n',xaxt='n',xlab='',main='Conservation of mouse devAS',ylab='proportion of devAS',t='n',ylim=c(0,0.7))
# for(t in unique(dev.as.cons.all$t)){
# 	f = dev.as.cons.all$t==t
# 	plotArea(x[f],dev.as.cons.all[f,c('also.sgn','ci1','ci2')],col=params$tissue.col[t],lwd=2)
# }
# axis(1,1:5,sps)
# lp=legend('topright',pch=19,col=params$tissue.col,legend=names(params$tissue.col),bty='n',ncol=2)
#legend(lp$rect$left,lp$rect$top,pch=params$species.pch[names(d)],legend=names(d),bty='n',xjust = 1)
#plotPanelLetter('B',lab.cex)

stat = t(sapply(age.ad.ph,function(x){m=mean(x);s=sd(x)/sqrt(length(x));c(m-2*s,m,m+2*s)}))
x = 1:nrow(stat)
cols=c('gray','black',params$tissue.col[rownames(stat)[-1:-2]])
par(mar=c(3.5,2.5,1.5,2))
plot(x,stat[,2],col=cols,main="Phastcons",ylab="mean intronic PhastCons (±50nt)",pch=19,bty='n',xaxt='n',xlab='',cex=2,ylim=range(stat))
abline(h=stat[1:2,2],lty=3,col=c('gray','black'))
segments(x,stat[,1],x,stat[,3],col=cols)
text(x,0.236,rownames(stat),adj=c(0,0),srt=-45,xpd=T)
plotPanelLetter('C',lab.cex)


plot(x,ad.3n[,2],col=cols,main="3N",ylab="fraction of 3N exons",pch=19,bty='n',xaxt='n',xlab='',cex=2,ylim=range(ad.3n))
abline(h=ad.3n[1:2,2],lty=3,col=c('gray','black'))
segments(x,ad.3n[,1],x,ad.3n[,3],col=cols)
text(x,0.35,rownames(ad.3n),adj=c(0,0),srt=-45,xpd=T)
plotPanelLetter('D',lab.cex)


plot(x,idr02[,2],col=cols,main="IUPR",ylab="fraction of IUPR",pch=19,bty='n',xaxt='n',xlab='',cex=2,ylim=range(idr02))
abline(h=idr02[1:2,2],lty=3,col=c('gray','black'))
segments(x,idr02[,1],x,idr02[,3],col=cols)
text(x,0.25,rownames(idr02),adj=c(0,0),srt=-45,xpd=T)
plotPanelLetter('E',lab.cex)

par(tck=-0.01,mgp=c(1.3,0.2,0),mar=c(0,0,0,0),oma=c(0,0,0,0))
lx=grconvertX(1,'lines','ndc')
plotExampleDLG3(fig=c(0.5,1-lx,2/5,1),'B')
dev.off()

# versions of B
# f = dev.as.cons$s1=='mouse' & !(dev.as.cons$s2 %in% c('macaque','chicken'))
# ll = dev.as.cons[f,]
# d = setNames(seq(-.3,.3,length.out = 4),rownames(species)[c(4,5,1,6)])
# plot(ll$sp.cnt + d[ll$s2],ll$also.sgn,pch=19,col=params$tissue.col[ll$t],t='n',xlab='# of AS species',ylab='fraction of devAS',bty='n',main='DevAS conservation with mouse')
# for(t in unique(ll$t))
# 	for(n in unique(ll$sp.cnt)){
# 		f = which(ll$t==t & ll$sp.cnt==n)
# 		f = f[order(d[ll$s2[f]])]
# 		lines(ll$sp.cnt[f] + d[ll$s2[f]],ll$also.sgn[f],pch=params$species.pch[ll$s2[f]],col=params$tissue.col[t],t='b')
# 	}
# lp=legend('topleft',pch=19,col=params$tissue.col,legend=names(params$tissue.col),bty='n')
# legend(lp$rect$left+lp$rect$w,lp$rect$top,pch=params$species.pch[names(d)],legend=names(d),bty='n')
# plot(ll$sp.cnt + d[ll$s2],ll$also.sgn,pch=19,col=params$tissue.col[ll$t],t='n',xlab='# of AS species',ylab='fraction of devAS',bty='n',main='DevAS conservation with mouse')
# for(t in unique(ll$t))
# 	for(s2 in names(d)){
# 		f = which(ll$t==t & ll$s2==s2)
# 		f = f[order(d[ll$sp.cnt[f]])]
# 		lines(ll$sp.cnt[f] + d[ll$s2[f]],ll$also.sgn[f],pch=params$species.pch[ll$s2[f]],col=params$tissue.col[t],t='b',lwd=2)
# 	}
# lp=legend('topleft',pch=19,col=params$tissue.col,legend=names(params$tissue.col),bty='n')
# legend(lp$rect$left+lp$rect$w,lp$rect$top,pch=params$species.pch[names(d)],legend=names(d),bty='n')
# plotPanelLetter('B1',lab.cex)
# 
# 
# f = dev.as.cons$s1=='mouse' & (dev.as.cons$s2 %in% c('rabbit'))
# ll = dev.as.cons[f,]
# plot(ll$sp.cnt + d[ll$s2],ll$also.sgn,pch=19,col=params$tissue.col[ll$t],t='n',xlab='# of AS species',ylab='fraction of devAS',bty='n',main='DevAS conservation: mouse-rabbit')
# for(t in unique(ll$t))
# 	for(s2 in names(d)){
# 		f = which(ll$t==t & ll$s2==s2)
# 		f = f[order(d[ll$sp.cnt[f]])]
# 		lines(ll$sp.cnt[f] + d[ll$s2[f]],ll$also.sgn[f],pch=params$species.pch[ll$s2[f]],col=params$tissue.col[t],t='b',lwd=2)
# 	}
# lp=legend('topleft',pch=19,col=params$tissue.col,legend=names(params$tissue.col),bty='n')
# #legend(lp$rect$left+lp$rect$w,lp$rect$top,pch=params$species.pch[names(d)],legend=names(d),bty='n')
# plotPanelLetter('B2',lab.cex)



# _improve 3B #####
dev.as.cons.all = list()
sps = c('rat','rabbit','human','opossum','chicken')

for(fi in 1:4){
	for(s2 in sps){
		if(fi == 1) f = T
		if(fi == 2) f = grepl(species['mouse','short'],alt.sp) & grepl(species[s2,'short'],alt.sp)
		if(fi == 3) f = apply(orth.per.tissue.age.qv$mouse < 0.05 & abs(orth.age.dpsi$mouse) > 0.2 & orth.per.tissue.age.qv[[s2]] < 0.05 & abs(orth.age.dpsi[[s2]]) > 0.2,1,sum,na.rm=T)>0
		if(fi == 4) f = alt.sp == 'hqmrboc'
		for(t in unique(meta$tissue)){
			cat(fi,s2,t,sum(f),"\n")
			r =  getDevASCons('mouse',s2,t,orth.per.tissue.age.qv,orth.age.dpsi,f)
			fc = abs(orth.age.dpsi$mouse[,t])>0.2 & orth.per.tissue.age.qv$mouse[,t]<0.05 & abs(orth.age.dpsi[[s2]][,t])>0.2 & orth.per.tissue.age.qv[[s2]][,t]<0.05
			fc[is.na(fc)] = FALSE
			cor = calcDevASCor(orth.seg.ad.tsm,age.al.i,'mouse',s2,t,fc & f)
			r = cbind(filter=fi,r,cor.both.mean=cor[1],cor.both.ci1=cor[2],cor.both.ci2=cor[3])
			
			fc = abs(orth.age.dpsi$mouse[,t])>0.2 & orth.per.tissue.age.qv$mouse[,t]<0.05# & abs(orth.age.dpsi[[s2]][,t])>0.2 & orth.per.tissue.age.qv[[s2]][,t]<0.05
			fc[is.na(fc)] = FALSE
			cor = calcDevASCor(orth.seg.ad.tsm,age.al.i,'mouse',s2,t,fc & f)
			r = cbind(r,cor.mean=cor[1],cor.ci1=cor[2],cor.ci2=cor[3])
			
			dev.as.cons.all[[length(dev.as.cons.all)+1]] = r
		}
	}
}

dev.as.cons.all = do.call(rbind,dev.as.cons.all)

names = c('All','ALT in both','devAS in both','Ancient')
pdf('figures/paper.figures/5/3/3B.pdf',w=8,h=6,family='Arial')
par(mfcol=c(3,4),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(1.5,2.5,1.5,0),oma=c(0,0,0,0))
for(fi in 1:4){
	data = dev.as.cons.all[dev.as.cons.all[,1]==fi,]
	x = match(data$s2,sps)
	plot(x,data$also.sgn,pch=params$species.pch[data$s2],col=params$tissue.col[data$t],bty='n',xaxt='n',xlab='',main=names[fi],ylab='proportion of devAS',t='n',ylim=c(0,1))
	for(t in unique(data$t)){
		f = data$t==t
		plotArea(x[f],data[f,c('also.sgn','ci1','ci2')],col=params$tissue.col[t],lwd=2)
	}
	axis(1,1:5,species[sps,'short'])
	if(fi==1)
		lp=legend('topright',pch=19,col=params$tissue.col,legend=names(params$tissue.col),bty='n',ncol=2)

	
	plot(x,data$cor.mean,pch=params$species.pch[data$s2],col=params$tissue.col[data$t],bty='n',xaxt='n',xlab='',main='',ylab='mean PCC',t='n',ylim=c(0,1))
	for(t in unique(data$t)){
		f = data$t==t
		plotArea(x[f],data[f,c('cor.mean','cor.ci1','cor.ci2')],col=params$tissue.col[t],lwd=2)
	}
	axis(1,1:5,species[sps,'short'])
	
	plot(x,data$cor.both.mean,pch=params$species.pch[data$s2],col=params$tissue.col[data$t],bty='n',xaxt='n',xlab='',main='',ylab='mean PCC',t='n',ylim=c(0,1))
	for(t in unique(data$t)){
		f = data$t==t
		plotArea(x[f],data[f,c('cor.both.mean','cor.both.ci1','cor.both.ci2')],col=params$tissue.col[t],lwd=2)
	}
	axis(1,1:5,species[sps,'short'])
}
text(grconvertX(0.5,'ndc','user'),grconvertY(2/3,'ndc','user'),labels='devAS in mouse',xpd=NA)	
text(grconvertX(0.5,'ndc','user'),grconvertY(1/3,'ndc','user'),labels='devAS in both species',xpd=NA)	
dev.off()


# 4 ######
# _prepare #####
mouse.liver.cor = calcCor2TisOnDev(psi.tsm$mouse[anns$mouse$sites=='ad',],'liver',meta.tsm,cor.meth = 'pearson')

as.in.ge.patterns.mouse = readRDS('Rdata/paper.figures/as.in.ge.patterns.mouse.Rdata') 
DPSI = 0.5
as.in.ge.patterns.mouse.cnt = lapply(as.in.ge.patterns.mouse,function(x){table(factor(x[,3]),pmax(abs(x[,1]),abs(x[,2]))>DPSI)[,c('TRUE','FALSE')]})
as.in.ge.patterns.mouse.stat = sapply(as.in.ge.patterns.mouse.cnt,function(x){cbind(my.binom.test(x['Increasing',]),my.binom.test(x['Decreasing',]))})

m = 'pearson'
s1 = 'mouse'

as.cor.on.dev = list()
s1 = 'mouse'
for(s2 in c('rat','rabbit','human','opossum')){
	qv = orth.per.tissue.age.qv
	f = grepl(species[s1,'short'],alt.sp) & grepl(species[s2,'short'],alt.sp)
	for(s in names(qv)){
		qv[[s]][is.na(orth.age.dpsi[[s]]) | abs(orth.age.dpsi[[s]]) < 0.2] = 1 # to set dPSI threshould
		qv[[s]] = qv[[s]][f,]
	}
	as.cor.on.dev[[s2]]=calcBootstrapSpeciesDiv(lapply(orth.seg.ad.tsm,function(x)x[f,]),c(s1,s2),function(x){x[s1,s2]},age.al.i,0,qv=qv,cor.meth = m,sign.both = F)[,,1]
}


# peak changes 
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
mystage = paste(my2marg.stage$species,my2marg.stage$tissue, my2marg.stage$stage)
mrstage = paste(my2marg.stage$species,my2marg.stage$tissue,my2marg.stage$marg.stage)
x = split(mystage,mrstage)
y = split(mrstage,mystage)
x[sapply(x,length)>1]
y[sapply(y,length)>1]

# code below is not completely correct since my stages has ambiguous correspondence to Margarida stages
my2marg.stage = setNames(paste(my2marg.stage$species,my2marg.stage$tissue,my2marg.stage$stage),tolower(paste(my2marg.stage$species,my2marg.stage$tissue,my2marg.stage$marg.stage)))
id = tolower(paste(ge.cnt$species,ge.cnt$tissue,ge.cnt$t1))
id[!(id %in% names(my2marg.stage))]
rownames(ge.cnt) = my2marg.stage[id]

'rabbit heart 9mpb' %in% rownames(ge.cnt);#my2marg.stage

dpsi = 0.2
peak.changes = list()
for(sp in rownames(species)[-2]){
	dpsi.age = calcdPSIonAge(psi.tsm[[sp]],meta.tsm,stages2use=NULL)
	for(t in unique(meta$tissue)){
		f = grep(t,colnames(dpsi.age))
		as = apply(abs(dpsi.age[anns[[sp]]$sites=='ad' & per.tissue.age.qv[[sp]][,t]<0.05,f])>  dpsi,2,function(x){x = x[!is.na(x)];if(length(x)==0){NA}else{sum(x)}})
		as = as[!is.na(as) & names(as) %in% rownames(ge.cnt)]
		peak.changes[[sp]][[t]] = list(ge=ge.cnt[names(as),'total'],as=as)
	}
}


# _plot #####

pdf('figures/paper.figures/5/4/4.pdf',w=6,h=6,family='Arial')
#jpeg('figures/paper.figures/5/4/4.jpg',units = 'in',w=6,h=6,quality = 100,res = 600,family='Arial')
layout(matrix(c(1,2,3,3,4,4),ncol=2,byrow = T),heights = c(2,1,1))
par(tck=-0.01,mgp=c(1.2,0.3,0),mar=c(4.5,2.5,1.5,0),oma=c(0,0,0,1),las=3)
plotTissueAgeProile(mouse.liver.cor,meta.tsm,age.axis = 'rank',bty='n',xlab='',ylab='Correlation to liver',ylim=c(0.8,1),pch=19,cex=1)
plotPNG("figures/paper.figures/5/icons/mouse.png",0.63,0.9,0.15)
plotPNG("figures/paper.figures/5/icons/liver.png",0.8,0.9,0.15)
title(xlab='Developmental stage',mgp=c(2,0.3,0))
plotPanelLetter('A',lab.cex)

col = c(paste0(rep(params$tissue.col,each=2),c('30','FF')))
b = barplot(as.in.ge.patterns.mouse.stat[c(4,1),],col=col,border=NA,ylim=c(0,max(as.in.ge.patterns.mouse.stat)),beside = T,ylab='Proportion of genes with |dPSI| > 0.5',xaxt='n')
segments(b,as.in.ge.patterns.mouse.stat[c(5,2),],b,as.in.ge.patterns.mouse.stat[c(6,3),])
text(apply(b,2,mean),rep(0,7),colnames(as.in.ge.patterns.mouse.stat),xpd=NA,adj=c(0,1),srt=-45)

plot4B.legend()
plotPanelLetter('B',lab.cex)

par(mar=c(3,1.2,0,0))
plot4C.DivergenceOnAge('mouse','human',as.cor.on.dev$human,'C')
#peak.changes
par(mar=c(3,1.2,1,0))
plot4C.PeakChange('rabbit',peak.changes$rabbit,'D')
dev.off()

# _5 ######
# _prepare #####
age.segs = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,psi.thr = 0.2,border.stages,s)[anns[[s]]$sites=='ad',])
names(age.segs) = rownames(species)
for(s in names(age.segs)) age.segs[[s]][is.na(per.tissue.age.qv[[s]][rownames(age.segs[[s]]),colnames(age.segs[[s]])])] = '-'
for(s in names(age.segs)) age.segs[[s]][age.segs[[s]] != '-' & per.tissue.age.qv[[s]][rownames(age.segs[[s]]),colnames(age.segs[[s]])]>0.05] = 'n'

# B phastcons
age.ad.ph = readRDS('Rdata/paper.figures/age.ad.ph.50nt.Rdata')

#D
s='human'
up.ad.len = apply(age.segs[[s]],2,function(t){
	f = anns[[s]]$cod=='c' & rownames(anns[[s]]) %in% rownames(age.segs[[s]])[t=='u']
	s = sum(anns[[s]][f,'length'] %% 3 == 0);
	t = sum(f)
	c(len3=s,total=t,freq=s/t,conf=binom.test(s,t)$conf.int)
})

dw.ad.len = apply(age.segs[[s]],2,function(t){
	f = anns[[s]]$cod=='c' & rownames(anns[[s]]) %in% rownames(age.segs[[s]])[t=='d']
	s = sum(anns[[s]][f,'length'] %% 3 == 0);
	t = sum(f)
	c(len3=s,total=t,freq=s/t,conf=binom.test(s,t)$conf.int)
})

up2dw.len.pv = sapply(1:7,function(t){prop.test(c(up.ad.len[1,t],dw.ad.len[1,t]),c(up.ad.len[2,t],dw.ad.len[2,t]))$p.value})
if(s=='mouse')
 	cnst.len3.freq=c(0.3996714,0.3942682,0.4050930)#my.binom.test(table(all.anns[[s]][all.anns[[s]]$sites=='ad' & all.anns[[s]]$type=='EXN' & all.anns[[s]]$cod=='c' & all.anns[[s]]$gene_id %in% all.anns[[s]][rownames(age.segs[[s]])[apply(age.segs[[s]] =='u' | age.segs[[s]] =='d',1,sum)>0],'gene_id'],'length'] %% 3 == 0)[c('TRUE','FALSE')])
if(s=='human')
	cnst.len3.freq=c(0.4052597,0.3997698,0.4107675)#my.binom.test(table(all.anns[[s]][all.anns[[s]]$sites=='ad' & all.anns[[s]]$type=='EXN' & all.anns[[s]]$cod=='c' & all.anns[[s]]$gene_id %in% all.anns[[s]][rownames(age.segs[[s]])[apply(age.segs[[s]] =='u' | age.segs[[s]] =='d',1,sum)>0],'gene_id'],'length'] %% 3 == 0)[c('TRUE','FALSE')])
t = anns[[s]][rownames(age.segs[[s]]),]
ndevAS.len3.freq = my.binom.test(table(t$length[t$cod=='c' & apply(age.segs[[s]],1,function(x){sum(x=='n') > 0 && sum(x %in% c('u','d')==0)})]  %% 3 == 0)[c('TRUE','FALSE')])

age.ad.gnomad = readRDS('Rdata/paper.figures/age.ad.gnomad.50nt.Rdata')
# hexamers
hex.dws.age02sgn = readRDS('Rdata/hex.dws.age02sgn.Rdata')
hex.ups.age02sgn = readRDS('Rdata/hex.ups.age02sgn.Rdata')

hex.ups = hex.ups.age02sgn
hex.dws = hex.dws.age02sgn

hex.stat = rbind(apply(hex.ups$up$ih.qv<0.05,2,sum),
								 apply(hex.dws$up$ih.qv<0.05,2,sum),
								 apply(hex.ups$dw$ih.qv<0.05,2,sum),
								 apply(hex.dws$dw$ih.qv<0.05,2,sum))

f = function(m,t){
	rbind(incl.ups = hex.ups$up$pv[m,t,],
				incl.dws = hex.dws$up$pv[m,t,],
				excl.ups = hex.ups$dw$pv[m,t,],
				excl.dws = hex.dws$dw$pv[m,t,])[,7:1]
}

pv.b= f('actaac','brain')
pv.h= f('actaac','heart')

fa = readRDS('Rdata/ad.alt.fa.Rdata')

# proportion of known
hex2mot = read.table('output/hex2mot2sf.tab.gz')
hex.tis.no = apply(hex.ups$up$ih.qv<0.05 | hex.dws$up$ih.qv<0.05 | hex.ups$dw$ih.qv<0.05 | hex.dws$dw$ih.qv<0.05,1,sum)
known.hex.stat = table(pmin(hex.tis.no,5),known=names(hex.tis.no) %in% hex2mot$V1[hex2mot$V2!=''])
rownames(known.hex.stat)[6] = '>5'

k = rownames(hex.ups$up$ih.qv) %in% rownames(hex2mot)[hex2mot$V2!='']
hex.stat.known = rbind(apply(hex.ups$up$ih.qv[k,]<0.05,2,sum),
											 apply(hex.dws$up$ih.qv[k,]<0.05,2,sum),
											 apply(hex.ups$dw$ih.qv[k,]<0.05,2,sum),
											 apply(hex.dws$dw$ih.qv[k,]<0.05,2,sum))

# _plot #####

pdf('figures/paper.figures/5/5/5.pdf',w=7.2,h=7.2,family='Arial')
#jpeg('figures/paper.figures/5/5/5.jpg',units = 'in',w=7.2,h=7.2,quality = 100,res = 600,family='Arial')
par(mfrow=c(3,3),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(4.5,2.5,1.5,0.3),oma=c(0,0,0,1),las=3)
plot.new()
plotPanelLetter('A',lab.cex)

at = (0:22)[-seq(3,22,by = 3)]
at[1] = -1
cols = c('gray','black',rep(params$tissue.col,each=2))
pch=c(1,1,rep(c(19,1),times=7))
xax = setNames(c(-1,1,seq(3.5,22,by=3)),c('const.','non-devAS',colnames(up.ad.len)))
stat = matrix(NA,ncol=3,nrow=16)
stat[1,] = cnst.len3.freq[c(2,1,3)]
stat[2,] = ndevAS.len3.freq[c(2,1,3)]
stat[seq(3,16,by = 2),] = t(up.ad.len)[,c(4,3,5)]
stat[seq(4,16,by = 2),] = t(dw.ad.len)[,c(4,3,5)]

plotASSegStat(at,stat,up2dw.len.pv,cols,xax,lty=c(2,2,rep(1:2,times=7)),main='Proportion of 3N exons',ylab='proportion',pch=pch,bty='n')
plotPanelLetter('B',lab.cex)

stat = t(sapply(age.ad.ph,function(x){m=mean(x);s=sd(x)/sqrt(length(x));c(m-2*s,m,m+2*s)}))
pv = sapply(1:7,function(i)wilcox.test(age.ad.ph[[2*i]],age.ad.ph[[2*i+1]])$p.value)
plotASSegStat(at,stat,pv,cols,xax,lty=c(2,2,rep(1:2,times=7)),main="50nt intron conservation ",ylab="mean PhastCons",pch=pch,bty='n')
plotPanelLetter('C',lab.cex)

par(mar=c(4.5,2.5,1.5,2))
stat = t(sapply(age.ad.gnomad,function(x){m=mean(x);s=sd(x)/sqrt(length(x));c(m-2*s,m,m+2*s)}))
pv = sapply(1:7,function(i)wilcox.test(age.ad.gnomad[[2*i]],age.ad.gnomad[[2*i+1]])$p.value)
plotASSegStat(at,stat,pv,cols,xax,lty=c(2,2,rep(1:2,times=7)),main="50nt intron human constraints" ,ylab="log10(SNP freq)",pch=pch,bty='n')
pv = sapply(2:15,function(i)wilcox.test(age.ad.gnomad$cnst,age.ad.gnomad[[i]])$p.value)
text(at[length(at)-2]+0.3,-3,c('p-value\n* <0.05\n** <0.01\n*** <0.001'),xpd=NA,adj=c(0,1))
plotPanelLetter('D')

par(mar=c(4.5,2.5,1.5,0.3))
den=50
barplotWithText(hex.stat[c(1,3),],col=paste0(rep(params$tissue.col,each=2),'55'),las=3,main='upstream',beside = T,border=NA,den=c(-1,den),srt = 90,adj=c(-0.1,0.5),ylim=c(0,max(hex.stat[c(1,3),])*1.2),xaxt='n',ylab='# hexamers')
b=barplot(hex.stat.known[c(1,3),],col=rep(params$tissue.col,each=2),add=T,beside = T,xaxt='n',yaxt='n',den=c(-1,den),border=NA)

text(apply(b,2,mean),rep(0,7),colnames(hex.stat),xpd=NA,srt=-45,adj=c(0,1))

plotPanelLetter('E')
barplotWithText(hex.stat[c(2,4),],col=paste0(rep(params$tissue.col,each=2),'55'),las=3,main='downstream',beside = T,border=NA,den=c(-1,den),srt = 90,adj=c(-0.1,0.5),ylim=c(0,max(hex.stat[c(2,4),])*1.2),xaxt='n',ylab='# hexamers')
b=barplot(hex.stat.known[c(2,4),],col=rep(params$tissue.col,each=2),add=T,beside = T,xaxt='n',yaxt='n',den=c(-1,den),border=NA)
text(apply(b,2,mean),rep(0,7),colnames(hex.stat),xpd=NA,srt=-45,adj=c(0,1))
legend('topright',fill=c('black','black','black','#00000055'),den=c(-1,den,-1,-1),legend=c('inclusion','exclusion','known','unknown'),bty = 'n',border = NA)
plotPanelLetter('F')

plotMirroredMotFreq(fa,age.segs,'actaac','brain','heart',main = 'ACTAAC in brain devAS')
plotPanelLetter('G')
plotMirroredMotFreq(fa,age.segs,'actaac','heart','brain',main = 'ACTAAC in heart devAS',plot.leg = F)

dev.off()


# 5J
pdf('figures/paper.figures/5/5/5J.pdf',w=6,h=12,family='Arial')
pat = c('tactaac','actaac','actaa','ctaac','ctaa')
par(mfrow=c(length(pat),2),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(4.5,2.5,1.5,0.3),oma=c(0,0,0,1),las=3)
for(p in pat){
	plotMirroredMotFreq(fa,age.segs,p,'brain','heart',main = paste0(toupper(p),' in brain devAS'))
	plotMirroredMotFreq(fa,age.segs,p,'heart','brain',main = paste0(toupper(p),' in heart devAS'),plot.leg = F)
}
dev.off()

# hexamers per species
dim(hex.dws.age02sgn$up$pv)
hex.dws.age02sgn$up$pv[1:10,,'mouse']
apply(hex.dws.age02sgn$up$pv,2:3,function(x)sum(p.adjust(x,m='BH')<0.05))
apply(hex.dws.age02sgn$dw$pv,2:3,function(x)sum(p.adjust(x,m='BH')<0.05))

# Suppl #####
# _2 #####
cols = list(const='gray',alt='orange')
sp=2
wd=6
pdf('figures/paper.figures/5/suppl/S2.pdf',w=7,h=3)
#jpeg('figures/paper.figures/5/suppl/S2.jpg',units = 'in',w=7,h=3,quality = 100,res = 600)
par(mfrow=c(1,2),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(0,0,0,0),oma=c(0,0,0,0))
plotSegDef('A',cols,sp,wd)
plotSAJR.AS.Q('B',wd)
dev.off()

# _3 #####
pdf('figures/paper.figures/5/suppl/S3.pdf',w=7,h=3,family='Arial')
#jpeg('figures/paper.figures/5/suppl/S3.jpg',units = 'in',w=7,h=3,quality = 100,res = 600,family='Arial')
layout(matrix(1:3,ncol=3,byrow = T),widths = c(3.2,0.6,3.2))
par(tck=-0.01,mgp=c(1.3,0.2,0),mar=c(1.5,2.5,1.5,0),oma=c(0,0,0,0))
# order by species/tissue
plotAsEventCount(tested.stat,astypes.pchs,by.tissue = F,ylab='# of detected events',main='Detected AS',bty='n')
plotPanelLetter('A',lab.cex)
par(mar=c(2.5,0,1.5,0))
y=plotASTypes(astypes.pchs)
legend(0,y-5,pch=19,col=params$tissue.col,legend = names(params$tissue.col),bty='n',cex = 0.7,xjust=0,yjust=1)
par(mar=c(1.5,2.5,1.5,0))
plotAsEventCount(sgn02.stat,astypes.pchs,by.tissue = F,ylab='% of devAS',main='DevAS',bty='n')
plotPanelLetter('B',lab.cex)
dev.off()

# _4 ######
# code is in devAS

# _5 #####
ens.ge.marg.tsm = readRDS('Rdata/ens.ge.marg.tsm.Rdata')
use.mean.embryo = FALSE
c2e = list()
for(s in rownames(species)){
	c2e[[s]] = list()
	ge = ens.ge.marg.tsm[[s]]
	ge = ge + min(ge[ge!=0],na.rm=T)
	ge = log2(ge)
	
	#sd = apply(psi.tsm[[s]],1,sd)
	#c2e[[s]]$psi.cor2embryo = caclCor2Embryo(psi.tsm[[s]][anns[[s]]$sites=='ad' & sd > 0.05,],meta.tsm,cor.m = 'sp',use.mean.embryo=use.mean.embryo)
	c2e[[s]]$psi.cor2embryo = caclCor2Embryo(psi.tsm[[s]][anns[[s]]$sites=='ad',],meta.tsm,cor.m = 'p',use.mean.embryo=use.mean.embryo)
	c2e[[s]]$ge.cor2embryo = caclCor2Embryo(ge,meta.tsm,cor.m = 'p',use.mean.embryo=use.mean.embryo)
}

pdf('figures/paper.figures/5/suppl/S5.pdf',w=7.2,h=7.2/3*7.4,family='Arial')
#jpeg('figures/paper.figures/5/suppl/S5.jpg',units = 'in',w=7.2,h=7.2/3*7.4,quality = 100,res = 300,family='Arial')
cors = readRDS('Rdata/ad.all.cor.Rdata')
par(mfrow= c(6,3),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(6,2.5,1.5,0),oma=c(2,0,0,0))
# MDS
for(s in rownames(species)[-3]){
	mds = cmdscale(1-cors[[species[s,'short']]],k=2)
	m = meta[rownames(mds),]
	plot(mds,pch=19,col=m$col,cex=m$cex,xlab='Dim 1',ylab='Dim 2',main=firstToupper(s),bty='n')
# cor-to-embrio
	tiss = intersect(unique(meta$tissue)[c(1,3,5:7)],dimnames(c2e[[s]]$psi.cor2embryo)[[1]])
	plotCor2Embryo(c2e[[s]]$psi.cor2embryo[tiss,,],main=paste('Splicing'),ylab='Pearson cor. to embryo, PSI',lwd=3,area.transp=0.1,xlab='',bty='n')
	plotCor2Embryo(c2e[[s]]$ge.cor2embryo[tiss,,],main=paste('Gene expression'),ylab='Pearson cor. to embryo, log(RPKM)',lwd=3,area.transp=0.1,xlab='',bty='n')
}
dev.off()

# _6 ######
tissues = setdiff(unique(meta$tissue),'liver')
tissues = setNames(tissues,tissues)
m.pe.cor = lapply(tissues,function(t)calcCor2TisOnDev(psi.tsm$mouse[anns$mouse$sites=='ad',],t,meta.tsm,cor.meth = 'pe'))
h.pe.cor = lapply(tissues,function(t)calcCor2TisOnDev(psi.tsm$human[anns$human$sites=='ad',],t,meta.tsm,cor.meth = 'pe'))

pdf('figures/paper.figures/5/suppl/S6.pdf',w=7.2,h=7.2/3*2,family='Arial')
#jpeg('figures/paper.figures/5/suppl/S6.jpg',units = 'in',w=7.2,h=7.2/3*2,quality = 100,res = 300,family='Arial')
par(mfrow=c(2,3),tck=-0.01,mgp=c(1.5,0.3,0),mar=c(4,2.5,1.5,0),oma=c(0,0,0,1),las=3)

for(i in 1:3){
	plotTissueAgeProile(h.pe.cor[[i]],meta.tsm,age.axis = 'rank',bty='n',xlab='',ylab=paste0('Correlatin to ',tissues[i]))
	plotPNG("figures/paper.figures/5/icons/human.png",0.15,0.15,0.15)
	plotPNG(paste0("figures/paper.figures/5/icons/",tissues[i],".png"),0.30,0.15,0.15)
	plotPanelLetter(LETTERS[i],lab.cex)
}

for(i in 1:3){
	plotTissueAgeProile(m.pe.cor[[3+i]],meta.tsm,age.axis = 'rank',bty='n',xlab='',ylab=paste0('Correlatin to ',tissues[i+3]))
	plotPNG("figures/paper.figures/5/icons/mouse.png",0.15,0.15,0.15)
	plotPNG(paste0("figures/paper.figures/5/icons/",tissues[i+3],".png"),0.30,0.15,0.15)
	plotPanelLetter(LETTERS[3+i],lab.cex)
}
dev.off()

# _7 div-on-age #####

pdf('figures/paper.figures/5/suppl/S7.pdf',w=6,h=4.5,family='Arial')
#jpeg('figures/paper.figures/5/suppl/S7.jpg',units = 'in',w=6,h=4.5,quality = 100,res = 600,family='Arial')
par(tck=-0.01,mgp=c(1.3,0.2,0),mar=c(3,1.2,0,0),oma=c(0,0,0,1),las=3,cex=2/3)
plot.new()
for(i in 1:3){
	s = c('rat','rabbit','opossum')[i]
	plot4C.DivergenceOnAge('mouse',s,as.cor.on.dev[[s]],LETTERS[i],yrange = c((3-i)/3,(4-i)/3),plot.xlab=i==3,ylab.pos=0.3,panel.lab.yadj=1,plot.tissue.lab=i==1)
}
dev.off()


# _8 peak change #####

pdf('figures/paper.figures/5/suppl/S8.pdf',w=6,h=5*1.2,family='Arial')
#jpeg('figures/paper.figures/5/suppl/S8.jpg',units = 'in',w=6,h=5*1.2,quality = 100,res = 600,family='Arial')
par(tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,1.2,0,0),oma=c(0,0,0,1),las=3,cex=2/3)
plot.new()
for(i in 1:5){
	s = rownames(species)[c(-2,-5)][i]
	plot4C.PeakChange(s,peak.changes[[s]],LETTERS[i],yrange = c((5-i)/5,(6-i)/5),plot.tissue.lab = i==1,plot.xlab = i==5)
}
dev.off()

# _table S1 ####
#  pc-lncRNA - possible coding
all.anns = readRDS('Rdata/all.anns.plus.cds.pos.Rdata')
s = 'human'
t = rep('',nrow(all.anns[[s]]))
t[!is.na(all.anns[[s]]$ens.transc.cds.pos) & grepl('cds',all.anns[[s]]$ens.transc.cds.pos)] = 'cds'
t[!is.na(all.anns[[s]]$ens.transc.cds.pos) & t != 'cds' & grepl('utr',all.anns[[s]]$ens.transc.cds.pos)] = 'utr'

table(t,all.anns[[s]]$sites)
table(all.anns[[s]]$cds,all.anns[[s]]$sites)

s2e = list()
for(s in rownames(species)){
	print(s)
	cg = unique(ens.exon.transc.cds.pos[[s]]$gene_id[ens.exon.transc.cds.pos[[s]]$gene_biotype=='protein_coding'])
	a = unique(unlist(seg2ens[[s]]))
	a = setNames(a %in% cg,a)
	s2e[[s]] = lapply(seg2ens[[s]],function(x)x[a[x]])
}

# overlap ens CDS
# overlap ens ncod exons of prot cod genes
# connected to prot cod ens gene
# overlap lnkRNA
# overlap any ens
# absolutely new!
ens.exon.transc.cds.pos = readRDS('Rdata/ens.exon.transc.cds.pos.Rdata')
res = list()
for(s in rownames(species)){
	print(s)
	t = setNames(rep('unannotated',nrow(all.anns[[s]])),rownames(all.anns[[s]]))
	# overlap any ens exon
	e = ens.exon.transc.cds.pos[[s]][,c(1,4,5,7)]
	colnames(e) = c('chr_id','start','stop','strand')
	e = getAnnOverlap(all.anns[[s]],e)
	table(e)
	t[e!='-'] = 'ens non-coding genes'
	#lncRNA
	t[all.anns[[s]]$lncRNA!='-' | all.anns[[s]]$pc.lncRNA!='-'] = 'lncRNA'
	#have link to prot-coding ens
	t[names(t) %in% names(seg2ens[[s]])[sapply(s2e[[s]],length)>0]] = 'ens coding genes: new exons'
	# overlap any ens exon from prot-cod genes
	e = ens.exon.transc.cds.pos[[s]][ens.exon.transc.cds.pos[[s]]$gene_biotype=='protein_coding',c(1,4,5,7)]
	colnames(e) = c('chr_id','start','stop','strand')
	e = getAnnOverlap(all.anns[[s]],e)
	table(e)
	t[e!='-'] = 'ens non-coding exons of coding genes'
	# overlap CDS
	t[all.anns[[s]]$cod != 'n'] = 'ens coding exons'
	types = c('ens coding exons','ens non-coding exons of coding genes','ens coding genes: new exons','lncRNA','ens non-coding genes','unannotated')
	res[[s]] = table(t,all.anns[[s]]$sites)[types,c('ad','aa','dd','da')]
}

sapply(res,function(x)x[,'aa'])
(function(){
	cat("Cassette exons\n")
	write.table(sapply(res,function(x)x[,'ad']),sep='\t',quote = F)
	cat("Alternative acceptor sites\n")
	write.table(sapply(res,function(x)x[,'aa']),sep='\t',quote = F)
	cat("Alternative donor sites\n")
	write.table(sapply(res,function(x)x[,'dd']),sep='\t',quote = F)
	cat("Retained introns\n")
	write.table(sapply(res,function(x)x[,'da']),sep='\t',quote = F)
})()

# additional #####
# _devAS on evolution #####
orth.sgn = list()
for(s in names(orth.age.dpsi)){
	orth.sgn[[s]] = abs(orth.age.dpsi[[s]])>0.2 & orth.per.tissue.age.qv[[s]] < 0.05
	orth.sgn[[s]][is.na(orth.sgn[[s]])] = FALSE
}
orth.devAS.cnt = sapply(orth.sgn,function(x)apply(x,2,sum))
barplot(orth.devAS.cnt,beside=T,col=params$tissue.col)
barplot(t(orth.devAS.cnt),beside=T,col=rep(params$tissue.col,each=7))

orth.devAS.prop.evo = lapply(names(orth.sgn),function(s){
	f = grep(species[s,'short'],alt.sp)
	sgn = orth.sgn[[s]][f,]
	sps = nchar(alt.sp[f])
	r = array(NA,dim=c(7,7,4),dimnames = list(colnames(orth.sgn[[s]],1:7,c('mean','ci1','ci2',cnt))))
	for(c in 1:7){
		sgn. = sgn[sps==c,]
		for(t in colnames(orth.sgn[[s]])){
			r[t,c,] = c(my.binom.test(sum(sgn.[,t]),sum(!sgn.[,t])),sum(sgn.[,t]))
		}
	}
	r
	})
names(orth.devAS.prop.evo) = names(orth.sgn)

orth.devAS.prop.evo$mouse[,,4]
pdf('figures/paper.figures/5/devAS.per.species-tissue.on.number-of-alt-species.pdf',w=8,h=4,family='Arial')
par(mfrow=c(2,4),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1,0),oma=c(0,0,0,1))
for(s in names(orth.devAS.prop.evo)){
	d = orth.devAS.prop.evo[[s]]
	plot(1,t='n',bty='n',xlab='# of species',ylab='proportion of devAS',main=s,xlim=c(1,7),ylim=c(0,0.4))
	for(t in dimnames(d)[[1]])
		plotArea(1:7,d[t,,1:3],col = params$tissue.col[[t]])
}
dev.off()
f = apply(orth.sgn$mouse,1,sum)>0 & apply(orth.sgn$opossum,1,sum)>0
fisher.test(table(orth.sgn$mouse[f,'ovary'],orth.sgn$opossum[f,'ovary']))
