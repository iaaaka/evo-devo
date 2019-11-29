options(stringsAsFactors = FALSE)
source('code/r.functions/ad.on.ge.F.R')
source('code/r.functions/load.all.data.F.R')
source('code/r.functions/paper.figures.4.F.R')
source('~/skoltech/r.code/util.R')

library(rphylopic)
library(SAJR)
library(GO.db)
# library(ape)

# library(GenomicRanges)

# load data #####
species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
anns = readRDS('Rdata/anns.Rdata')
border.stages = readRDS('Rdata/border.stages.Rdata')
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')
params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)
orth.seg.ad.tsm = readRDS('Rdata/orth.seg.ad.tsm.Rdata')
seg2ens = readRDS('Rdata/seg2ens.Rdata')
ens.ge.cod = readRDS('Rdata/ens.ge.cod.Rdata')
orth.per.tissue.age.qv = readRDS('Rdata/orth.per.tissue.age.qv.Rdata')
age.al.i = read.table('input/margarida.age.alignment.final.2017.05.17.improved.tab',sep='\t',header = T)[,c(rownames(species),'to.remove')]

# orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
# orth.seg.ad.all = readRDS('Rdata/orth.seg.ad.all.Rdata')
# orth.seg.ad.all.id = readRDS('Rdata/orth.seg.ad.all.id.Rdata')
# ens.ge.cod.tsm = readRDS('Rdata/ens.ge.cod.tsm.Rdata')
# orth.ens.genes = readRDS('Rdata/orth.ens.genes.Rdata')
# all.anns = readRDS('Rdata/all.anns.Rdata')
# 
# ens.ge.marg = readRDS('Rdata/ens.ge.marg.Rdata')
# ens.ge.marg.tsm = readRDS('Rdata/ens.ge.marg.tsm.Rdata')

# orth.per.tissue.age.qv = readRDS('Rdata/orth.per.tissue.age.qv.Rdata')
#per.tissue.age.amp = readRDS('Rdata/per.tissue.age.amp.Rdata')
#per.tissue.age.sgn = lapply(1:7,function(s){per.tissue.age.qv[[s]]<0.05 & per.tissue.age.amp[[s]]>0.2})
#names(per.tissue.age.sgn) = names(per.tissue.age.amp)
# hex.dws.age03 = readRDS('Rdata/hex.dws.age03.Rdata')
# hex.ups.age03 = readRDS('Rdata/hex.ups.age03.Rdata')

# hex.dws.age02sgn = readRDS('Rdata/hex.dws.age02sgn.Rdata')
# hex.ups.age02sgn = readRDS('Rdata/hex.ups.age02sgn.Rdata')
# 
# 
# 
# params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
# params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)
# DPSI = 0.5
# 
# ens.descr.mm = unique(read.table('input/mm.38.84.gene.descr.txt',sep=',',quote='"',header=TRUE))
# rownames(ens.descr.mm) = ens.descr.mm$Ensembl.Gene.ID
# ens.descr.mm$Description = sapply(strsplit(ens.descr.mm$Description,' [',TRUE),'[',1)

# prepare data #######
# # _MDS ####
# orth.all.psi = do.call(cbind,lapply(orth.seg.ad,function(x)x$ir))
# orth.all.psi = orth.all.psi[,rownames(meta)]
# 
# alt.cnt = apply(sapply(orth.seg.ad,function(x)x$seg$type=='ALT'),1,sum)
# table(alt.cnt)
# orth.cor = list()
# orth.cor$all = cor(orth.all.psi,u='p')
# orth.cor$sp1 = cor(orth.all.psi[alt.cnt==1,],u='p')
# orth.cor$sp7 = cor(orth.all.psi[alt.cnt==7,],u='p')
# orth.cor$sp7tsm = cor(do.call(cbind,orth.seg.ad.tsm)[alt.cnt==7,],u='p')
# saveRDS(orth.cor,'Rdata/paper.figures/orth.cor.Rdata')
# orth.mds2 = lapply(orth.cor,function(x)cmdscale(1-x,k=2))
# saveRDS(orth.mds2,'Rdata/paper.figures/orth.mds2.Rdata')

# _GO #####
GOALLANCESTOR = c(as.list(GOCCANCESTOR),as.list(GOMFANCESTOR),as.list(GOBPANCESTOR))
mouse.go = read.table('input/GO/Mus_musculus.GRCm38.84.GO.csv.gz',sep=',',header = T)
mouse.go = mouse.go[mouse.go[,2] != '',]
mouse.go = split(mouse.go$GO.Term.Accession,mouse.go$Ensembl.Gene.ID)
mouse.go.rev = revList(mouse.go)
 
mouse.go.full = lapply(mouse.go,function(x){r=unique(c(x,unlist(GOALLANCESTOR[x])));r[grep('GO:',r,fixed=T)]})
mouse.go.full.rev = revList(mouse.go.full)



### figure 1 #####
orth.mds2 = readRDS('Rdata/paper.figures/orth.mds2.Rdata')
filtered.seg.cnt = sapply(anns,getNoOfEvents,gene=FALSE)
age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)
age.dpsi$macaque = cbind(age.dpsi$macaque[,1:5],ovary=NaN,age.dpsi$macaque[,6,drop=FALSE])


sign.stat02 = list(all =sapply(names(per.tissue.age.qv),function(s)apply((abs(age.dpsi[[s]]) >  0.2 & per.tissue.age.qv[[s]] < 0.05)[anns[[s]]$sites=='ad',],2,sum,na.rm=T)),
									 up  =sapply(names(per.tissue.age.qv),function(s)apply((    age.dpsi[[s]]  >  0.2 & per.tissue.age.qv[[s]] < 0.05)[anns[[s]]$sites=='ad',],2,sum,na.rm=T)),
									 down=sapply(names(per.tissue.age.qv),function(s)apply((    age.dpsi[[s]]  < -0.2 & per.tissue.age.qv[[s]] < 0.05)[anns[[s]]$sites=='ad',],2,sum,na.rm=T)))

sign.stat05 = list(all =sapply(names(per.tissue.age.qv),function(s)apply((abs(age.dpsi[[s]]) >  0.5 & per.tissue.age.qv[[s]] < 0.05)[anns[[s]]$sites=='ad',],2,sum,na.rm=T)),
									 up  =sapply(names(per.tissue.age.qv),function(s)apply((    age.dpsi[[s]]  >  0.5 & per.tissue.age.qv[[s]] < 0.05)[anns[[s]]$sites=='ad',],2,sum,na.rm=T)),
									 down=sapply(names(per.tissue.age.qv),function(s)apply((    age.dpsi[[s]]  < -0.5 & per.tissue.age.qv[[s]] < 0.05)[anns[[s]]$sites=='ad',],2,sum,na.rm=T)))

tested.stat = sapply(names(per.tissue.age.qv),function(s)apply(!is.na(per.tissue.age.qv[[s]][anns[[s]]$sites=='ad',]),2,sum))

short.tiss = substr(rownames(tested.stat),1,1)

# orth.age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(orth.seg.ad.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
# names(orth.age.dpsi) = rownames(species)

# _plot ####
pdf('figures/paper.figures/4/1.pdf',w=9,h=6)
layout(matrix(c(1,1,2:5),byrow = TRUE,ncol=3))
par(tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
plot.new()
text(0.5,0.5,'study design\nspecies x tissue x time',cex=3)
plotPanelLetter('A')

b=barplot(t(tested.stat),beside = T,col=rep(params$tissue.col,each=7),ylab='Number of exons',main='Cassette exons passed filtering',names.arg=short.tiss)
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.7)
plotPanelLetter('B')

m = meta
SAJR::plotMDS(points=-orth.mds2$all,col=m$col,cex=m$cex*1.5,pch=m$pch,main=paste0('All cassette exons (',nrow(orth.seg.ad.tsm$human),')'))
legend('topleft',col=params$tissue.col,pch=19,legend=names(params$tissue.col),bty='n')
plotPanelLetter('C')

b=barplot(t(sign.stat02$all),den=50,beside = T,col=rep(params$tissue.col,each=7),ylab='Number of exons',main='Number of devAS',names.arg=short.tiss)
barplot(t(sign.stat05$all),beside = T,col=rep(params$tissue.col,each=7),ylab='',main='',add=T,xaxt='n')
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.7)
legend(18,2700,den=c(c(-1,50)),legend=c('q<0.05 & dPSI>0.5','q<0.05 & dPSI>0.2'))
plotPanelLetter('D')

b=barplot(t(sign.stat02$all/tested.stat)*100,den=50,beside = T,col=rep(params$tissue.col,each=7),ylab='% of devAS in tested',main='% of devAS in tested',names.arg=short.tiss)
barplot(t(sign.stat05$all/tested.stat)*100,beside = T,col=rep(params$tissue.col,each=7),ylab='',main='',add=T,xaxt='n')
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.7)
plotPanelLetter('E')
dev.off()


# figure 2 (4) ####
# make MDSs
# mdata = readRDS('Rdata/mouse.as.u.filtered.Rdata')
# mouse.ad.cor = cor(mdata$ir[mdata$seg$sites=='ad',],u='p')
# mouse.ad.cor.mds = cmdscale(1-mouse.ad.cor,k=2)
# saveRDS(mouse.ad.cor.mds,'Rdata/paper.figures/mouse.ad.cor.mds.Rdata')

# make as in ge patterns
# m = makeASinGEpatterns('mouse')
# h = makeASinGEpatterns('human')
# saveRDS(m,'Rdata/paper.figures/as.in.ge.patterns.mouse.Rdata') 
# saveRDS(h,'Rdata/paper.figures/as.in.ge.patterns.human.Rdata')
# m = makeASinGEpatterns('mouse',use.random.seg=TRUE)
# h = makeASinGEpatterns('human',use.random.seg=TRUE)
# saveRDS(m,'Rdata/paper.figures/as.in.ge.patterns.mouse.random.exon.Rdata') 
# saveRDS(h,'Rdata/paper.figures/as.in.ge.patterns.human.random.exon.Rdata')


mouse.cassett.cnt = 20186 #sum(mdata$seg$sites=='ad')
ancient.cnt = 1590 #sum(apply(sapply(orth.seg.ad,function(x)x$seg$type=='ALT'),1,sum)==7)
mouse.ad.cor.mds = readRDS('Rdata/paper.figures/mouse.ad.cor.mds.Rdata')

# corr to embryo
sp.ce = 'mouse'
mouse.ge = readRDS('Rdata/ens.ge.marg.tsm.Rdata')[[sp.ce]]
mouse.ge = mouse.ge + min(mouse.ge[mouse.ge!=0],na.rm=T)
mouse.ge = log2(mouse.ge)

use.mean.embryo = FALSE
#sd = apply(psi.tsm[[sp.ce]],1,sd)
#psi.cor2embryo = caclCor2Embryo(psi.tsm[[sp.ce]][anns[[sp.ce]]$sites=='ad' & sd > 0.05,],meta.tsm,cor.m = 'sp',use.mean.embryo=use.mean.embryo)
psi.cor2embryo = caclCor2Embryo(psi.tsm[[sp.ce]][anns[[sp.ce]]$sites=='ad',],meta.tsm,cor.m = 'p',use.mean.embryo=use.mean.embryo)
ge.cor2embryo = caclCor2Embryo(mouse.ge,meta.tsm,cor.m = 'p',use.mean.embryo=use.mean.embryo)
cisbp = read.table('input/cisbp-db/all/RBP_Information_all_motifs.txt',sep='\t',header = T,quote = '',fill = T)
sf.cor2embryo = caclCor2Embryo(mouse.ge[rownames(mouse.ge) %in% cisbp$DBID[cisbp$RBP_Species==c(mouse='Mus_musculus',human='Homo_sapiens')[sp.ce]],],meta.tsm,cor.m = 'p',use.mean.embryo=use.mean.embryo)

# AS in GE patterns
as.in.ge.patterns.mouse = readRDS('Rdata/paper.figures/as.in.ge.patterns.mouse.Rdata') 
as.in.ge.patterns.human = readRDS('Rdata/paper.figures/as.in.ge.patterns.human.Rdata')

# as.in.ge.patterns.mouse = readRDS('Rdata/paper.figures/as.in.ge.patterns.mouse.random.exon.Rdata') 
# as.in.ge.patterns.human = readRDS('Rdata/paper.figures/as.in.ge.patterns.human.random.exon.Rdata')

DPSI = 0.5

as.in.ge.patterns.mouse.cnt.up = lapply(as.in.ge.patterns.mouse,function(x){table(factor(x[,3]),x[,1]>DPSI)[,c('TRUE','FALSE')]})
as.in.ge.patterns.mouse.stat.up = sapply(as.in.ge.patterns.mouse.cnt.up,function(x){cbind(my.binom.test(x['Increasing',]),my.binom.test(x['Decreasing',]))})
as.in.ge.patterns.mouse.cnt.dw = lapply(as.in.ge.patterns.mouse,function(x){table(factor(x[,3]),x[,2]< -DPSI)[,c('TRUE','FALSE')]})
as.in.ge.patterns.mouse.stat.dw = sapply(as.in.ge.patterns.mouse.cnt.dw,function(x){cbind(my.binom.test(x['Increasing',]),my.binom.test(x['Decreasing',]))})

as.in.ge.patterns.human.cnt.up = lapply(as.in.ge.patterns.human,function(x){table(factor(x[,3]),x[,1]>DPSI)[,c('TRUE','FALSE')]})
as.in.ge.patterns.human.stat.up = sapply(as.in.ge.patterns.human.cnt.up,function(x){cbind(my.binom.test(x['Increasing',]),my.binom.test(x['Decreasing',]))})
as.in.ge.patterns.human.cnt.dw = lapply(as.in.ge.patterns.human,function(x){table(factor(x[,3]),x[,2]< -DPSI)[,c('TRUE','FALSE')]})
as.in.ge.patterns.human.stat.dw = sapply(as.in.ge.patterns.human.cnt.dw,function(x){cbind(my.binom.test(x['Increasing',]),my.binom.test(x['Decreasing',]))})

pdf('figures/paper.figures/2018.12.20/4F.random.exon.pdf',w=4,h=4)
plotASinGEPatterns(as.in.ge.patterns.mouse.stat.up,as.in.ge.patterns.mouse.stat.dw)
dev.off()

# duplication and ohnologs
#mouse
ge.info.m = read.csv('/home/mazin/skoltech/projects/evo.devo/input/gene.info.from.marg/Mouse.Indexes.csv')
rownames(ge.info.m) = ge.info.m$Mouse_ID
m.age.ens = getAgeASEns(psi.tsm,meta.tsm,DPSI,border.stages,'mouse')
ohnologs.m = read.table('input/ohnologs/MOUSE.Pairs.Intermediate.2R.txt',sep='\t',header = T)
#human
ge.info.h = read.csv('/home/mazin/skoltech/projects/evo.devo/input/gene.info.from.marg/Human.Indexes.All.csv')
rownames(ge.info.h) = ge.info.h$Human_ID
h.age.ens = getAgeASEns(psi.tsm,meta.tsm,DPSI,border.stages,'human')
ohnologs.h = read.table('input/ohnologs/HUMAN.Pairs.Intermediate.2R.txt',sep='\t',header = T)

### _plot ####
pdf('figures/paper.figures/4/2.pdf',w=9,h=12)
par(mfrow=c(4,3),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,2.5),oma=c(0,0,0,1))

m = meta.tsm[rownames(orth.mds2$sp7tsm),]
SAJR::plotMDS(points=orth.mds2$sp7tsm,col=m$col,cex=m$cex*1.5,pch=m$pch,main=paste0('All cassette exons (',ancient.cnt,')'))
plotPanelLetter('A')

m = meta[rownames(mouse.ad.cor.mds),]
SAJR::plotMDS(points=mouse.ad.cor.mds,col=m$col,cex=m$cex*2,pch=19,main=paste0('All mouse cassette exons (',mouse.cassett.cnt,')'))
plotPanelLetter('B')

plot.new()
legend('topleft',col=params$tissue.col,pch=19,legend=names(params$tissue.col),bty='n')

tiss = unique(meta$tissue)[c(1,3,5:7)]
par(mar=c(5.5,2.5,1.5,0),mgp=c(1.3,0.2,0))
plotCor2Embryo(psi.cor2embryo[tiss,,],main='Splicing',ylab='Pearson cor. to embryo, PSI',lwd=3,area.transp=0.1,xlab='')
plotPanelLetter('C')
plotCor2Embryo(ge.cor2embryo[tiss,,],main='Gene expression',ylab='Pearson cor. to embryo, log(RPKM)',lwd=3,area.transp=0.1,xlab='')
plotPanelLetter('D')
plotCor2Embryo(sf.cor2embryo[tiss,,],main='SF gene expression',ylab='Pearson cor. to embryo, log(RPKM)',lwd=3,area.transp=0.1,xlab='')
plotPanelLetter('E')


plotASinGEPatterns(as.in.ge.patterns.mouse.stat.up,as.in.ge.patterns.mouse.stat.dw)
plotPanelLetter('F')

plotDevASFreq(m.age.ens,ge.info.m$Mouse_ID[!is.na(ge.info.m$Age) & ge.info.m$Age>3],'mam. dupl.','topleft',main='DevAS are depleted in recent duplications',
							filter = intersect(ge.info.m$Mouse_ID[!is.na(ge.info.m$Age)],rownames(ens.ge.cod$mouse$gene)))
plotPanelLetter('G')

plotDevASFreq(m.age.ens,c(ohnologs.m$Ohnolog.1.Id,ohnologs.m$Ohnolog.2.Id),'ohnologs','topright',main='DevAS are enriched in ohnologs')
plotPanelLetter('H')

plotASinGEPatterns(as.in.ge.patterns.human.stat.up,as.in.ge.patterns.human.stat.dw)
plotPanelLetter('F.1')

plotDevASFreq(h.age.ens,ge.info.h$Human_ID[!is.na(ge.info.h$Age) & ge.info.h$Age>3],'mam. dupl.','topleft',main='DevAS are depleted in recent duplications',
							filter = intersect(ge.info.h$Human_ID[!is.na(ge.info.h$Age)],rownames(ens.ge.cod$human$gene)))

plotPanelLetter('G.1')
plotDevASFreq(h.age.ens,c(ohnologs.h$Ohnolog.1.Id,ohnologs.h$Ohnolog.2.Id),'ohnologs','topright',main='DevAS are enriched in ohnologs')
plotPanelLetter('H.1')

dev.off()


# figure 3 (10) #####
#A
s = 'mouse'
age.segs = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,psi.thr = 0.2,border.stages,s)[anns[[s]]$sites=='ad',])
names(age.segs) = rownames(species)
for(s in names(age.segs)) age.segs[[s]][is.na(per.tissue.age.qv[[s]][rownames(age.segs[[s]]),colnames(age.segs[[s]])])] = '-'
for(s in names(age.segs)) age.segs[[s]][age.segs[[s]] != '-' & per.tissue.age.qv[[s]][rownames(age.segs[[s]]),colnames(age.segs[[s]])]>0.05] = 'n'

age.segs.cod = lapply(rownames(species),function(s)age.segs[[s]][anns[[s]][rownames(age.segs[[s]]),'cod']!='n',])
names(age.segs.cod) = rownames(species)

ts = unique(meta$tissue)
up.all = t(sapply(age.segs,function(x)apply(x=='u',2,sum)[ts]))
up.cod = t(sapply(age.segs.cod,function(x)apply(x=='u',2,sum)[ts]))

dw.all = t(sapply(age.segs,function(x)apply(x=='d',2,sum)[ts]))
dw.cod = t(sapply(age.segs.cod,function(x)apply(x=='d',2,sum)[ts]))

#B
ts = unique(meta$tissue)[-2]
up.ntiss =  t(sapply(age.segs[-2],function(x)table(factor(apply(x[,ts]=='u',1,sum),levels=1:length(ts)))))
dw.ntiss = -t(sapply(age.segs[-2],function(x)table(factor(apply(x[,ts]=='d',1,sum),levels=1:length(ts)))))

#C
age.ad.over = lapply(age.segs,function(x){
	x[x=='-'] = NA
	x = cbind(x=='u',x=='d')
	#x = x[apply(x,1,sum,na.rm=T)>0,]
	colnames(x) = paste(colnames(x),rep(c('up','dw'),each=ncol(x)/2))
	caclSegOverlap(x)
})

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

#E phastcons
# h = readRDS('Rdata/all.anns.Rdata')$human
# gc()
# phastcons = read.table('/home/mazin/skoltech/projects/evo.devo/processed/ad.phastcons.gz',sep='\t')
# phastcons = setNames(strsplit(phastcons[,2],',',TRUE),phastcons[,1])
# phastcons = lapply(phastcons,as.numeric)
# age.ad.ph = list()
# as = age.segs$human[anns$human[rownames(age.segs$human),'cod']=='c',]
# inx = 151:200
# h.age.sids = list()
# h.age.sids$cnst = rownames(h)[h$sites=='ad' & h$cod=='c' & h$type=='EXN' & h$gene_id %in% h[rownames(as),'gene_id']]
# s='human'
# t = anns[[s]][rownames(age.segs[[s]]),]
# h.age.sids$`non-devAS` = rownames(t)[t$cod=='c' & apply(age.segs[[s]],1,function(x){sum(x=='n') > 0 && sum(x %in% c('u','d')==0)})]
# for(t in colnames(as)){
# 	h.age.sids[[paste(t,'up')]] = rownames(as)[as[,t]=='u']
# 	h.age.sids[[paste(t,'dw')]] = rownames(as)[as[,t]=='d']
# }
# 
# age.ad.ph = lapply(h.age.sids,function(sids)sapply(phastcons[intersect(names(phastcons),sids)],function(x)mean(x[c(inx,length(x)-inx)])))
# saveRDS(age.ad.ph,'Rdata/paper.figures/age.ad.ph.50nt.Rdata')
age.ad.ph = readRDS('Rdata/paper.figures/age.ad.ph.50nt.Rdata')

#F gnomad
# gnomad = read.table('processed/gnomad201/human.ad.snp.tab.gz',sep='\t')
# colnames(gnomad) = c('chr_id','pos','id','ref','alt','qual','filter','seg_id','alt_cnt','freq','tot_cnt')
# hann = h[unique(unlist(h.age.sids)),]
# gnomad = gnomad[gnomad$alt_cnt > 1 & gnomad$seg_id %in% rownames(hann),]
# s2g = hann[gnomad$seg_id,]
# gnomad$dist2start = gnomad$pos - s2g$start
# gnomad$dist2stop  = s2g$stop - gnomad$pos
# gnomad[s2g$strand==-1,c('dist2start','dist2stop')] = gnomad[s2g$strand==-1,c('dist2stop','dist2start')]
# gnomad$strand = s2g$strand
# 
# age.ad.gnomad = lapply(h.age.sids,function(sids)log10(gnomad$freq[gnomad$seg_id %in% sids & ((gnomad$dist2start<0 & gnomad$dist2start>= -50) | (gnomad$dist2stop<0 & gnomad$dist2stop>= -50))]))
# saveRDS(age.ad.gnomad,'Rdata/paper.figures/age.ad.gnomad.50nt.Rdata')
age.ad.gnomad = readRDS('Rdata/paper.figures/age.ad.gnomad.50nt.Rdata')

#G hgmd
# hgmd = read.table('input/hgmd/tosend/2017_1_HGMD_ALL_Variants.csv',sep='\t',row.names=1,header=T,quote='',comment.char = '')
# hgmd.gr = GRanges(hgmd$chrom_VCF_hg19,IRanges(hgmd$pos_VCF_hg19,hgmd$pos_VCF_hg19+nchar(hgmd$ref_VCF_hg19)-1))
# 
# seg.gr = GRanges(hann$chr_id,IRanges(hann$start,hann$stop))
# hgmd2hadf = findOverlaps(hgmd.gr,seg.gr,maxgap=200,type='any',select='all',ignore.strand=TRUE)
# hgmd2hadf = data.frame(h=hgmd2hadf@from,s=hgmd2hadf@to)
# hgmd2hadf$seg.id = rownames(hann)[hgmd2hadf$s]
# hgmd2hadf$hgmd.id = rownames(hgmd)[hgmd2hadf$h]
# 
# s = hann[hgmd2hadf$s,]
# h = hgmd[hgmd2hadf$h,]
# 
# hgmd2hadf$position = ''
# hgmd2hadf$position = paste0(hgmd2hadf$position,ifelse((s$strand== 1 & s$start > h$pos_VCF_hg19) | (s$strand==-1 & s$stop < h$pos_VCF_hg19+nchar(h$ref_VCF_hg19)-1),'u',''))
# hgmd2hadf$position = paste0(hgmd2hadf$position,ifelse(s$start<=h$pos_VCF_hg19+nchar(h$ref_VCF_hg19)-1 & s$stop>=h$pos_VCF_hg19,'i',''))
# hgmd2hadf$position = paste0(hgmd2hadf$position,ifelse((s$strand==-1 & s$start > h$pos_VCF_hg19) | (s$strand== 1 & s$stop < h$pos_VCF_hg19+nchar(h$ref_VCF_hg19)-1),'d',''))
# table(hgmd2hadf$position)
# # locaction: intron  (u,w), exon (e), ppt (p), acc.ag (A),acc(a), don(d),don.gt(D)
# # mutation annotation priproti: dinucleotides > other nt of canonical sites > PPT > intron
# hgmd2hadf$loc = ''
# hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse(hgmdOverlapLocaction(h,s,0,Inf,T) & hgmdOverlapLocaction(h,s,0,Inf,F),'e',''))
# hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse(hgmdOverlapLocaction(h,s,-2,-1,T),'A',''))
# hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse(hgmdOverlapLocaction(h,s,-2,-1,F),'D',''))
# hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse((!grepl('A|D',hgmd2hadf$loc)) & (hgmdOverlapLocaction(h,s,-3,-3,T) | hgmdOverlapLocaction(h,s,1,1,T)),'a',''))
# hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse((!grepl('A|D',hgmd2hadf$loc)) & (hgmdOverlapLocaction(h,s,-5,-3,F) | hgmdOverlapLocaction(h,s,1,4,F)),'d',''))
# hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse((!grepl('A|D|a|d',hgmd2hadf$loc)) & hgmdOverlapLocaction(h,s,-23,-4,T),'p',''))
# hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse((!grepl('A|D|a|d|p',hgmd2hadf$loc)) & hgmdOverlapLocaction(h,s,-200,-24,T),'u',''))
# hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse((!grepl('A|D|a|d|p',hgmd2hadf$loc)) & hgmdOverlapLocaction(h,s,-200,-6,F),'w',''))
# sort(table(hgmd2hadf$loc))
# hs = unique(hgmd2hadf$seg.id[grep('A|a|p|D|d|u|d',hgmd2hadf$loc)])
# age.ad.hgmd = t(sapply(h.age.sids,function(sids){r=c(sum(sids %in% hs),sum(!(sids %in% hs)));c(r,my.binom.test(r))}))

# saveRDS(age.ad.hgmd,'Rdata/paper.figures/age.ad.hgmd.200nt.Rdata')
age.ad.hgmd = readRDS('Rdata/paper.figures/age.ad.hgmd.200nt.Rdata')

sps = c('mouse','rat','rabbit','human','opossum','chicken')

orth.age.ad2 = lapply(rownames(species),function(s)getAgeASchanges(orth.seg.ad.tsm,meta.tsm,psi.thr = 0.2,border.stages,s))
names(orth.age.ad2) = rownames(species)
for(sp in names(orth.age.ad2)) orth.age.ad2[[sp]][orth.age.ad2[[sp]] != '-' & (is.na(orth.per.tissue.age.qv[[sp]]) | orth.per.tissue.age.qv[[sp]]>0.05)[,colnames(orth.age.ad2[[sp]])]] = 'n'

u2=getASChangeCons(orth.age.ad2,'u',sps)
d2=getASChangeCons(orth.age.ad2,'d',sps)


# _plot #####
pdf('figures/paper.figures/4/3.pdf',w=10,h=11)
pdf('figures/paper.figures/2018.12.20/10.pdf',w=10,h=11)
s = 'mouse'
l = matrix(c(1,2,4,1,2,5,3,3,6,3,3,7),ncol=4)
l = rbind(l,8:11)
layout(l)
par(tck=-0.02,mgp=c(1.1,0.2,0),mar=c(3,2,1.5,0),oma=c(0,0,0,1))
cols=rep(params$tissue.col[colnames(up.cod)],each=nrow(up.all))
b=barplot( up.all,beside = T,col=NA,border = cols,ylim=c(-max(dw.all,na.rm=T)*1.4,max(up.all,na.rm=T)),yaxt='n',ylab='# of exons',main='DevAS cassette exons')
barplot( up.cod,beside = T,col=cols,border = cols,add=T,density = 50,yaxt='n')
barplot(-dw.all,beside = T,col=NA,border = cols,add=T,yaxt='n')
barplot(-dw.cod,beside = T,col=cols,border = cols,add=T,density = 50,yaxt='n')
at = c(-2:5*200)
axis(2,at,abs(at))
text(b[,1],rep(-1500,nrow(up.all)),substr(rownames(up.all),1,3),srt=90,adj = c(0,0.5),cex=0.6)
text(b[21], 800,'inclusion',adj=c(-0.1,1.1))
text(b[21],-800,'exclusion',adj=c(-0.1,-1.1))
plotPanelLetter('A')


cols = rep(gray.colors(nrow(up.ntiss)),times=ncol(up.ntiss))
b=barplot(up.ntiss,beside = T,xlab='# of tissues',col=cols,border=NA,yaxt='n',ylim=range(up.ntiss,dw.ntiss*1.2),ylab='# of exons',main='Number of shared devAS exons',legend.text=T,args.legend = list(x='bottomright'))
barplot(dw.ntiss,beside = T,col=cols,border=NA,yaxt='n',ylim=range(up.ntiss,dw.ntiss),add=TRUE)
abline(h=0)
at = c(2000,1000,0,1000,2000)
axis(2,c(-1,-1,1,1,1)*at,at)
text(b[21], 2000,'inclusion',adj=c(-0.1,1.1))
text(b[21],-2000,'exclusion',adj=c(-0.1,-1.1))
plotPanelLetter('B')

par(mar=c(1.5,6,3,1))
plotAgeSegOverlap(age.ad.over[[s]],main=paste('Overlap of devAS exons across',s,'tissues'))
par(mar=c(3,2,1.5,0))
plotPanelLetter('C')

par(mar=c(5,2,1.5,0))
at = (0:22)[-seq(3,22,by = 3)]
cols = c('black','gray',rep(params$tissue.col,each=2))
pch=c(1,1,rep(c(19,1),times=7))
xax = setNames(c(0,1,seq(3.5,22,by=3)),c('const.','non-devAS',colnames(up.ad.len)))
stat = matrix(NA,ncol=3,nrow=16)
stat[1,] = cnst.len3.freq[c(2,1,3)]
stat[2,] = ndevAS.len3.freq[c(2,1,3)]
stat[seq(3,16,by = 2),] = t(up.ad.len)[,c(4,3,5)]
stat[seq(4,16,by = 2),] = t(dw.ad.len)[,c(4,3,5)]

plotASSegStat(at,stat,up2dw.len.pv,cols,xax,lty=c(2,2,rep(1:2,times=7)),main='Proportion of 3N exons',ylab='proportion',pch=pch)
plotPanelLetter('D')

stat = t(sapply(age.ad.ph,function(x){m=mean(x);s=sd(x)/sqrt(length(x));c(m-2*s,m,m+2*s)}))
pv = sapply(1:7,function(i)wilcox.test(age.ad.ph[[2*i]],age.ad.ph[[2*i+1]])$p.value)
plotASSegStat(at,stat,pv,cols,xax,lty=c(2,2,rep(1:2,times=7)),main="50nt intron conservation ",ylab="mean PhastCons",pch=pch)
plotPanelLetter('E')

stat = t(sapply(age.ad.gnomad,function(x){m=mean(x);s=sd(x)/sqrt(length(x));c(m-2*s,m,m+2*s)}))
pv = sapply(1:7,function(i)wilcox.test(age.ad.gnomad[[2*i]],age.ad.gnomad[[2*i+1]])$p.value)
plotASSegStat(at,stat,pv,cols,xax,lty=c(2,2,rep(1:2,times=7)),main="50nt intron human constrains" ,ylab="log10(SNP freq)",pch=pch)
pv = sapply(2:15,function(i)wilcox.test(age.ad.gnomad$cnst,age.ad.gnomad[[i]])$p.value)
legend('topright',title = 'p-value',legend=c('* <0.05','** <0.01','*** <0.001'),bty = 'n')
plotPanelLetter('F')

pv = sapply(1:7,function(t){s=c(age.ad.hgmd[t*2,1],age.ad.hgmd[t*2+1,1]);prop.test(s,s+c(age.ad.hgmd[t*2,2],age.ad.hgmd[t*2+1,2]))$p.value})
plotASSegStat(at,age.ad.hgmd[,c(4,3,5)],pv,cols,xax,lty=c(2,2,rep(1:2,times=7)),main='200nt intron HGMD mutations',ylab='fraction of exons with mutations',pch=pch)
plotPanelLetter('G')
legend('topright',pch=c(19,1),legend=c('inclusion','exclusion'))

plotTisUpDownCOns(u2$brain[,-1],d2$brain[,-1],'H',main='Brain ageAS conservation')
plotTisUpDownCOns(u2$heart[,-1],d2$heart[,-1],'I',main='Heart ageAS conservation')
plotTisUpDownCOns(u2$liver[,-1],d2$liver[,-1],'J',main='Liver ageAS conservation')
plotTisUpDownCOns(u2$testis[,-1],d2$testis[,-1],'K',main='Testis ageAS conservation')
legend('topright',fill=c('red','blue'),legend=c('Inclusion','Exclusion'))
dev.off()

#### figure 4 hexamers (11) ####
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

pdf('figures/paper.figures/2019.02.21/11.ACTAAC-QKI.pv.pdf',w=8,h=4)
par(mfrow=c(1,2),tck=-0.01,mgp=c(1.1,0.2,0),mar=c(6,4,1.5,1),oma=c(0,0,0,1))
imageWithText(log10(pv.b),names.as.labs = T,xlab='',ylab='',main='brain')
imageWithText(log10(pv.h),names.as.labs = T,xlab='',ylab='',main='heart')
dev.off()

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

# hexamer overlap
all.hex.qv = cbind(hex.ups$up$ih.qv, hex.dws$up$ih.qv, hex.ups$dw$ih.qv, hex.dws$dw$ih.qv)
colnames(all.hex.qv) = paste0(substr(colnames(all.hex.qv),1,1),rep(c('iu','id','eu','ed'),each=7))
all.hex.qv = all.hex.qv[,rep(c(1,8,15,22),times=7)+rep(0:6,each=4)]
f = apply(all.hex.qv<0.05,1,sum)

hex.overlap=calcAllPairsFT(all.hex.qv<0.05)
hex.overlap=calcAllPairsFT(all.hex.qv[apply(all.hex.qv<0.05,1,sum)>0,]<0.05)
# QKI 
DPSI = 0.2
fa = readRDS('Rdata/ad.alt.fa.Rdata')
age.segs = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,DPSI,border.stages,s)[anns[[s]]$sites=='ad',])
names(age.segs) = rownames(species)
for(s in names(age.segs)) age.segs[[s]][age.segs[[s]] != '-' & (is.na(per.tissue.age.qv[[s]]) | per.tissue.age.qv[[s]]>0.05)[rownames(age.segs[[s]]),colnames(age.segs[[s]])]] = 'n'

qki.coor = c(163983952,163988000)
b=c('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/6048sTS.Human.Brain.28ypb.Male.bam',
		'/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/5531sTS.Human.Brain.29ypb.Male.bam',
		'/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/5517sTS.Human.Brain.50ypb.Male.bam',
		'/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/5524sTS.Human.Brain.58ypb.Male.bam')
qki.brain=getReadCoverage(b,'6',qki.coor[1],qki.coor[2],NA,F,0)
b=c('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/6042sTS.Human.Heart.25ypb.Male.bam',
		'/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/6043sTS.Human.Heart.54ypb.Male.bam')
qki.hearta=getReadCoverage(b,'6',qki.coor[1],qki.coor[2],NA,F,0)
b=c('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/3588sTS.Human.Heart.CS13.Male.bam',
		'/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/3671sTS.Human.Heart.CS16.Male.bam',
		'/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/2088sTS.Human.Heart.CS18.Male.bam')
qki.hearte=getReadCoverage(b,'6',qki.coor[1],qki.coor[2],NA,F,0)

human.qki.gtf = read.table('processed/annotation/all.species/ensambl/Homo_sapiens.GRCh37.73.gtf.gz',sep='\t')
human.qki.gtf = unique(human.qki.gtf[(human.qki.gtf$V4 >= qki.coor[1] & human.qki.gtf$V4 <= qki.coor[2]) & 
																		 	(human.qki.gtf$V5 >= qki.coor[1] & human.qki.gtf$V5 <= qki.coor[2])	&
																		 	human.qki.gtf$V3 == 'CDS' & human.qki.gtf$V7 == '+' & human.qki.gtf$V1=='6',c(1,3,4,5,7,8)])
tmp = read.table('processed/annotation/all.species/merged/human.sajr',sep='\t')[,c(1,3,4,5,7,8,9)]
tmp = tmp[(tmp$V4 >= qki.coor[1] & tmp$V4 <= qki.coor[2]) & 
						(tmp$V5 >= qki.coor[1] & tmp$V5 <= qki.coor[2])	&
						tmp$V3 == 'segment' & tmp$V7 == '+' & tmp$V1=='6' & !(grepl('type=INT',tmp$V9)),-7]
human.qki.gtf = rbind(tmp,human.qki.gtf[human.qki.gtf$V3=='CDS',])
#there is a rarery used retained intron, I merged it with the exon
human.qki.gtf$V5[human.qki.gtf$V4==163986978 & human.qki.gtf$V3=='segment'] = 163987752

# _plot #####
pdf('figures/paper.figures/4/4.pdf',w=10,h=13)
m = matrix(c(1:11,11),ncol=3,byrow = T)
layout(m)
den=50
par(tck=-0.01,mgp=c(1.1,0.2,0),mar=c(5,3,1.5,1),oma=c(0,0,0,1))
barplotWithText(hex.stat[c(1,3),],col=paste0(rep(params$tissue.col,each=2),'55'),las=3,main='upstream',beside = T,border=NA,den=c(-1,den),srt = 90,adj=c(-0.1,0.5),ylim=c(0,max(hex.stat[c(1,3),])*1.2))
barplot(hex.stat.known[c(1,3),],col=rep(params$tissue.col,each=2),add=T,beside = T,xaxt='n',yaxt='n',den=c(-1,den),border=NA)
legend('topright',fill=c('black','black','black','#00000055'),den=c(-1,den,-1,-1),legend=c('inclusion','exclusion','known','unknown'),bty = 'n')

plotPanelLetter('A')
barplotWithText(hex.stat[c(2,4),],col=paste0(rep(params$tissue.col,each=2),'55'),las=3,main='downstream',beside = T,border=NA,den=c(-1,den),srt = 90,adj=c(-0.1,0.5),ylim=c(0,max(hex.stat[c(2,4),])*1.2))
barplot(hex.stat.known[c(2,4),],col=rep(params$tissue.col,each=2),add=T,beside = T,xaxt='n',yaxt='n',den=c(-1,den),border=NA)

plotPanelLetter('B')
par(mgp=c(2.1,0.2,0))
barplotWithText(known.hex.stat[,2]/(known.hex.stat[,1]+known.hex.stat[,2]),t = known.hex.stat[,1]+known.hex.stat[,2],srt = 90,adj = c(-0.1,0.5),ylim=c(0,1.2),xlab='# of tissues\nwhere hexamer is sign.',yaxt='n')
par(mgp=c(1.1,0.2,0))
title(ylab='proportion of known motifs')
axis(2,seq(0,1,0.2))
plotPanelLetter('C')
par(mar=c(3,3,1.5,1))

pv.col=c(gray=1,yellow=1e-2,orange=1e-5,red=1e-20,none=-1)
or.col=c(blue=0,'#0000FFAA'=1/3,'#0000FF55'=1/1.5,gray=1,yellow=5,orange=10,red=Inf)
plotFTMotifSimMatrix(hex.overlap,F,diag.text = rep('',7*4),main='Hexamer overlap',pv.col=pv.col,or.col=or.col)
rect(0:6*4,-0.3,1:7*4,-1.3,col=params$tissue.col,border = NA,xpd=T)
mtext('odds ratio',2)
mtext('p-value',4)
abline(v=0:7*4)
abline(h=0:7*4)
plotPanelLetter('D')

par(mar=c(3,3,1.5,3))
dt=c('incl.\nupstr.','incl.\ndownstr.','excl.\nupstr.','excl.\ndownstr.')
plotFTMotifSimMatrix(hex.overlap[1:4,1:4],T,main='Brain',diag.text = dt,pv.col=pv.col,or.col=or.col)
rect(0,-0.3/7,4,-1.3/7,col=params$tissue.col[1],border = NA,xpd=T)
par(mar=c(3,0,1.5,6))
plotFTMotifSimMatrix(hex.overlap[9:12,9:12],T,main='Heart',pv.col=pv.col,or.col=or.col,diag.text = dt)
rect(0,-0.3/7,4,-1.3/7,col=params$tissue.col[3],border = NA,xpd=T)
legend(4.1,4,fill=names(pv.col)[1:4],legend = paste0('<',pv.col[1:4]),title = 'p-value',xpd=T)
legend(4.1,2,fill=names(or.col)[4:7],legend = paste0('<',or.col[4:7]),title = 'odds ratio',xpd=T)
#QKI
par(mar=c(3,3,1.5,1))
plotMirroredMotFreq(fa,age.segs,'actaac','brain','heart',main = 'QKI motif in brain')
plotPanelLetter('E')
plotMirroredMotFreq(fa,age.segs,'actaac','heart','brain',main = 'QKI motif in heart',plot.leg = F)
plotPanelLetter('F')

#plotTissueAgeProile(ens.ge.cod$mouse$rpkm['ENSMUSG00000062078',],meta,ylab='RPKM',ylim=range(0,ens.ge.cod$mouse$rpkm['ENSMUSG00000062078',]),main='Expression of mouse QKI')
plotTissueAgeProile(ens.ge.cod$human$rpkm['ENSG00000112531',],meta,ylab='RPKM',ylim=range(0,ens.ge.cod$human$rpkm['ENSG00000112531',]),main='Expression of human QKI',age.axis = 'rank')
#plotTissueAgeProile(ens.ge.marg$human['ENSG00000112531',],meta,ylab='RPKM',ylim=range(0,ens.ge.marg$human['ENSG00000112531',]),main='Expression of human QKI',age.axis = 'rank')
plotPanelLetter('G')
plotTissueAgeProile(psi.tsm$human['hum.57513.s15',],meta.tsm,ylab='PSI',ylim=c(0,1),main='Splicing of last QKI exon',age.axis = 'rank')
plotPanelLetter('H')
plotQKICov('I')
dev.off()

# age.segs = lapply(age.segs,function(x){
# 	f = (x[,'brain'] == 'u' & x[,'heart'] == 'd') | (x[,'brain'] == 'd' & x[,'heart'] == 'u')
# 	x[f,c('brain','heart')] = 'n'
# 	x
# 	})
# 
# table(age.segs$human[,'brain'],age.segs$human[,'heart'])
# table(age.segs$mouse[,'brain'],age.segs$mouse[,'heart'])



pdf('figures/paper.figures/4/4EF-per.species.pdf',w=5,h=10)
par(mfrow=c(7,2),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
for(s in rownames(species)){
	plotMirroredMotFreq(fa[s],age.segs[s],'actaac','brain','heart',main = paste0('QKI motif in ',s,' brain'),plot.leg = F)
	plotMirroredMotFreq(fa[s],age.segs[s],'actaac','heart','brain',main = paste0('QKI motif in ',s,' heart'),plot.leg = F)
}
dev.off()
### figure 5 (8) #####
alt.sp = readRDS('Rdata/paper.figures/alt.sp.Rdata')
alt.sp = alt.sp[alt.sp!='']
ens.ge.marg.tsm = readRDS('Rdata/ens.ge.marg.tsm.Rdata')
orth = read.table('input/TFs.from.margarida/E85.MRRHMOC.one2one.txt',header=T)
colnames(orth) = tolower(sapply(strsplit(colnames(orth),'_'),'[',1))

f = sapply(names(ens.ge.marg.tsm),function(s)orth[,s] %in% rownames(ens.ge.marg.tsm[[s]]))
f = apply(f,1,sum)==7

ens.ge.cod.tsm.log = lapply(names(ens.ge.marg.tsm),function(s)log(ens.ge.marg.tsm[[s]][orth[f,s],]+1e-4))
names(ens.ge.cod.tsm.log) = names(ens.ge.marg.tsm)


#m = 'pearson'
m = 'spearman'
sign.both=TRUE
ge.mr.cor=calcBootstrapSpeciesDiv(ens.ge.cod.tsm.log,c('mouse','rat'),function(x){x['mouse','rat']},age.al.i,0,cor.meth = m,sign.both=sign.both)
ge.mh.cor=calcBootstrapSpeciesDiv(ens.ge.cod.tsm.log,c('mouse','human'),function(x){x['mouse','human']},age.al.i[-5,],0,cor.meth = m,sign.both=sign.both)
as.mr.cor=calcBootstrapSpeciesDiv(orth.seg.ad.tsm,c('mouse','rat'),function(x){x['mouse','rat']},age.al.i,0,qv=NULL,cor.meth = m,sign.both=sign.both)
as.mh.cor=calcBootstrapSpeciesDiv(orth.seg.ad.tsm,c('mouse','human'),function(x){x['mouse','human']},age.al.i[-5,],0,qv=NULL,cor.meth = m,sign.both=sign.both)

# ALT in both
f = grepl('mr',alt.sp)
as.mr.cor=calcBootstrapSpeciesDiv(lapply(orth.seg.ad.tsm,function(x)x[f,]),c('mouse','rat'),function(x){x['mouse','rat']},age.al.i,0,qv=NULL,cor.meth = m,sign.both=sign.both)
as.mh.cor=calcBootstrapSpeciesDiv(lapply(orth.seg.ad.tsm,function(x)x[f,]),c('mouse','human'),function(x){x['mouse','human']},age.al.i[-5,],0,qv=NULL,cor.meth = m,sign.both=sign.both)

as.mr.cor.devAS=calcBootstrapSpeciesDiv(orth.seg.ad.tsm,c('mouse','rat'),function(x){x['mouse','rat']},age.al.i,0,qv=orth.per.tissue.age.qv,cor.meth = m,sign.both=sign.both)
as.mh.cor.devAS=calcBootstrapSpeciesDiv(orth.seg.ad.tsm,c('mouse','human'),function(x){x['mouse','human']},age.al.i[-5,],0,qv=orth.per.tissue.age.qv,cor.meth = m,sign.both=sign.both)

# _what is wrong with kidney? 
# orth.seg.ad = readRDS('Rdata/orth.seg.ad.all.Rdata')
# seg.len = readRDS('Rdata/paper.figures/alt.seg.len.Rdata')
# table(orth.seg.ad$human$seg$north)
# 
# t = 'kidney'
# s1 = 'mouse'
# s2 = 'rat'
# sgn = orth.per.tissue.age.qv[[s1]][,t] < 0.05 & orth.per.tissue.age.qv[[s2]][,t] < 0.05
# sgn[is.na(sgn)] = FALSE
# table(sgn)
# c1 = paste(s1,t,age.al.i[,s1])
# c2 = paste(s2,t,age.al.i[,s2])
# f = c1 %in% colnames(orth.seg.ad.tsm[[s1]]) & c2 %in% colnames(orth.seg.ad.tsm[[s2]])
# table(f)
# 
# p1 = orth.seg.ad.tsm[[s1]][sgn,c1[f]]
# p2 = orth.seg.ad.tsm[[s2]][sgn,c2[f]]
# 
# 
# plot(seg.len[rownames(p1)],seg.len[rownames(p2)],log='xy')
# table(orth.seg.ad$human$seg[rownames(orth.seg.ad[[s1]]$seg) %in% rownames(p1),'north'])
# z = seg.len[rownames(p1)] == seg.len[rownames(p2)] & orth.seg.ad$human$seg[rownames(orth.seg.ad[[s1]]$seg) %in% rownames(p1),'north'] == 7
# table(z)
# z=T
# image(cor(cbind(p1,p2),u='p'))

#z = apply(is.na(cbind(p1,p2)),1,sum) == 0
#table(z)
#z = z & apply(p1,1,sd,na.rm=T)>0.02 & apply(p2,1,sd,na.rm=T)>0.02
#table(z)
# 
# plotLine(1:ncol(p1),sapply(1:ncol(p1),function(i)cor(p1[z,i],p2[z,i],u='p',m='p')))
# plotLine(1:ncol(p1),sapply(1:ncol(p1),function(i)cor(p1[z,i],p2[z,i],u='p',m='sp')))
# 
# par(mfrow=c(1,2))
# plotLine(p1[z,1],p2[z,1])
# plotLine(p1[z,11],p2[z,11])
# 
# hist(p1[z,11],20)
# plot(table(round(p1[z,11],1))/table(round(p1[z,1],1)))
# 
# table(p1[z,11]>0.1)
#alternification
# 
# orth.seg.ad.all = readRDS('Rdata/orth.seg.ad.all.Rdata')
# alt.sp = sapply(orth.seg.ad.all,function(x)x$seg$type=='ALT')
# alt.sp = apply(alt.sp,1,function(x)paste(species$short[x],collapse = ''))
# names(alt.sp) = rownames(orth.seg.ad.all$human$seg)
# saveRDS(alt.sp,'Rdata/paper.figures/alt.sp.Rdata')
alt.sp = readRDS('Rdata/paper.figures/alt.sp.Rdata')
sp.groups = c(species$short,'hq','mr','mrb','hqmrb')
mean = table(alt.sp)[sp.groups]

#born/lost
born.exn.sajr = readRDS('Rdata/born.exn.sajr.Rdata')
exon.birth.one = readRDS('Rdata/exon.birth.one.Rdata')
sp.birth = sapply(exon.birth.one,function(x){paste(species$short[!is.na(x$seg_id)],collapse='')})

tsp.birth = table(sp.birth)
brth=tsp.birth[sp.groups]
lost=tsp.birth[sapply(sp.groups,function(x){paste(setdiff(species$short,strsplit(x,'')[[1]]),collapse='')})]
names(lost) = sp.groups

plot.points = straight.line=TRUE


# _ plot ####
pdf('figures/paper.figures/4/5.sign.both.spearman.pdf',w=13,h=8)
par(mfrow=c(2,4),tck=-0.01,mgp=c(2,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,0,1))

x=setNames(rep('',11),sp.groups)
plotDivergenceOnAgeWithConf(as.mr.cor,params$tissue.col,ylab='Mouse-rat Pearson corr.',main='AS divergence',plot.points=plot.points,straight.line=straight.line)
plotPanelLetter('A')
plotDivergenceOnAgeWithConf(as.mr.cor.devAS,params$tissue.col,ylab='Mouse-rat Pearson corr.',main='AS divergence (devAS)',plot.points=plot.points,straight.line=straight.line)
plotPanelLetter('B')
plotDivergenceOnAgeWithConf(ge.mr.cor,params$tissue.col,ylab='Mouse-rat Pearson corr.',main='GE divergence',plot.points=plot.points,straight.line=straight.line)
plotPanelLetter('C')
plotSpeciesTree.3('D',x,mean,brth,lost)

plotDivergenceOnAgeWithConf(as.mh.cor,params$tissue.col,ylab='Mouse-human Pearson corr.',main='AS divergence',plot.points=plot.points,straight.line=straight.line)
plotPanelLetter('A.1')
plotDivergenceOnAgeWithConf(as.mh.cor.devAS,params$tissue.col,ylab='Mouse-human Pearson corr.',main='AS divergence (devAS)',plot.points=plot.points,straight.line=straight.line)
plotPanelLetter('B.1')
plotDivergenceOnAgeWithConf(ge.mh.cor,params$tissue.col,ylab='Mouse-human Pearson corr.',main='GE divergence',plot.points=plot.points,straight.line=straight.line)
plotPanelLetter('C.1')
dev.off()


# figure 6 (14)#####
# new exons
sps=list(`m/r` = c('m','r'),mr='mr',mrb='mrb',hqmrb='hqmrb',hqmrbo='hqmrbo')
sps. = c(sps,hqmrboc='hqmrboc')

exon.birth.one = readRDS('Rdata/exon.birth.one.Rdata')

sp.birth = sapply(exon.birth.one,function(x){paste(species$short[!is.na(x$seg_id)],collapse='')})
born.seg.ids = t(sapply(exon.birth.one,function(x){x$seg_id}))
colnames(born.seg.ids) = rownames(exon.birth.one[[1]])

born.exn.sajr = readRDS('Rdata/born.exn.sajr.Rdata')
born.exn.tsm = lapply(born.exn.sajr,function(x){
	psi = x$ir[,colnames(x$ir) %in% rownames(meta)]
	m = meta[colnames(psi),]
	calcMeanCols(psi,paste(m$species,m$tissue,m$stage))
})
born.per.tissue.age.qv = readRDS('Rdata/born.per.tissue.age.qv.Rdata')

age.al.i = read.table('input/margarida.age.alignment.final.2017.05.17.improved.tab',sep='\t',header = T)[,c(rownames(species),'to.remove')]
age.al = age.al.i[age.al.i$to.remove==0,1:7]
apply(age.al,2,function(x)sum(table(x[x!=''])>1))
age.al$chicken[c(6,8)] = ''

ts = unique(meta$tissue)
m = data.frame(mouse.stage=rep(age.al$mouse,times=length(ts)),tissue=rep(ts,each=nrow(age.al)))
rownames(m) = paste(m$tissue,m$mouse.stage)
d = unique(meta[meta$species=='mouse',c('days','stage')])
m$days = setNames(d$days,d$stage)[m$mouse.stage]
meta.tsm.al = m

for(s in names(born.exn.sajr)){
	tsm = born.exn.tsm[[s]]
	born.exn.sajr[[s]]$psi.tsm.al = matrix(NA,ncol=nrow(m),nrow=nrow(tsm),dimnames = list(rownames(tsm),rownames(m)))
	for(i in 1:nrow(m)){
		stage = age.al[age.al$mouse==m$mouse.stage[i],s]
		st = paste(s,m$tissue[i],stage)
		if(st %in% colnames(tsm))
			born.exn.sajr[[s]]$psi.tsm.al[,paste(m$tissue[i],m$mouse.stage[i])] = tsm[,st]
	}
}
meta.tsm.al$col = paste0(params$tissue.col[meta.tsm.al$tissue],rep(c('44','55','66','77','88','99','AA','BB','CC','DD','EE','FF'),times=7))


## cleaned
#minus:
# macaque and chicken
# kidney, ovary and cerebellum
# stages 10.5, 11.5, 3dpb, 9wpb

missed.cond = lapply(born.exn.sajr,function(x)colnames(x$psi.tsm.al)[apply(!is.na(x$psi.tsm.al),2,sum)==0])
lapply(missed.cond[c(-2,-7)],function(x)x[!grepl('10.5|11.5|3dpb|9wpb',x) & !grepl('cerebellum|kidney|ovary',x)])
mm = meta.tsm.al[!(meta.tsm.al$tissue %in% c('cerebellum','kidney','ovary') | meta.tsm.al$mouse.stage %in% c('10.5','11.5','3dpb','9wpb')),]
#mm = meta.tsm.al[!(meta.tsm.al$tissue %in% c('cerebellum','kidney','ovary','heart','liver') | meta.tsm.al$mouse.stage %in% c('10.5','11.5','3dpb','9wpb')),]

tmp = lapply(setNames(born.exn.sajr,NULL)[c(-2,-7)],function(x)x$psi.tsm.al[,rownames(mm)])
max.stage.cl=unlist(lapply(tmp,function(x)setNames(colnames(x)[apply(x,1,function(x){r = which.max(x);ifelse(length(r)==0 || sum(x==x[r],na.rm=T)>1,NA,r)})],rownames(x))))

born.sids = lapply(split.data.frame(born.seg.ids,sp.birth),function(x){x[!is.na(x)]})
tab.exn = sapply(sps, function(x){table(factor(max.stage.cl[intersect(unlist(born.sids[x]),names(max.stage.cl))],levels=rownames(meta.tsm.al)))})

# alternification
orth.seg.ad.all.tsm = readRDS('Rdata/orth.seg.ad.all.tsm.Rdata')
alt.sp = readRDS('Rdata/paper.figures/alt.sp.Rdata')

orth.tsm.al = lapply(names(orth.seg.ad.all.tsm),function(s){
	r = matrix(NA,ncol=nrow(meta.tsm.al),nrow=nrow(orth.seg.ad.all.tsm[[s]]),dimnames = list(rownames(orth.seg.ad.all.tsm[[s]]),rownames(meta.tsm.al)))
	for(i in 1:nrow(meta.tsm.al)){
		ss = age.al[age.al$mouse==meta.tsm.al$mouse.stage[i],s]
		cn = paste(s,meta.tsm.al$tissue[i],ss)
		if(cn %in% colnames(orth.seg.ad.all.tsm[[s]]))
			r[,i] = orth.seg.ad.all.tsm[[s]][,cn]
	}
	r
})
names(orth.tsm.al) = names(orth.seg.ad.all.tsm)
#min.stage = sapply(orth.tsm.al,function(x)colnames(x)[apply(x,1,function(x){r=which.min(x);ifelse(length(r)==0 || sum(x==x[r],na.rm=T)>1,NA,r)})])
#names(min.stage) = rownames(orth.tsm.al$human)

#for(i in 1:nrow(min.stage))	min.stage[i,rownames(species)[!(species$short %in% strsplit(alt.sp[i],'')[[1]])]] = NA

mm = meta.tsm.al[!(meta.tsm.al$tissue %in% c('cerebellum','kidney','ovary') | meta.tsm.al$mouse.stage %in% c('10.5','11.5','3dpb','9wpb')),]
tmp = lapply(orth.tsm.al[c(-2,-7)],function(x)x[,rownames(mm)])

min.stage.cl = sapply(tmp,function(x)colnames(x)[apply(x,1,function(x){r=which.min(x);ifelse(length(r)==0|| sum(x==x[r],na.rm=T)>1,NA,r)})])
rownames(min.stage.cl) = rownames(tmp$human)
for(i in 1:nrow(min.stage.cl))	min.stage.cl[i,colnames(min.stage.cl) %in% rownames(species)[!(species$short %in% strsplit(alt.sp[i],'')[[1]])]] = NA
#min.stage.cl=unlist(lapply(colnames(min.stage.cl),function(s)setNames(min.stage.cl[,s],rownames(tmp[[s]]))))


# orth.seg.ad.all.dpsi = lapply(orth.seg.ad.all.tsm,function(x){
# 	sapply(unique(meta$tissue),function(t){
# 		p = x[alt.sp!='',grep(t,colnames(x))]
# 		apply(p,1,max,na.rm=T) - apply(p,1,min,na.rm=T)
# 		})
# 	})

orth.seg.ad.all.dpsi = lapply(rownames(species),function(s)getAgeASchanges(orth.seg.ad.all.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(orth.seg.ad.all.dpsi) = rownames(species)

# par(mfrow=c(3,3))
# sapply(1:7,function(s){hist(orth.seg.ad.all.dpsi[[s]][orth.per.tissue.age.qv[[s]]>0.05],0:50/50,border = NA,col='gray',freq = F,main=names(born.tis.dpsi)[i]);hist(orth.seg.ad.all.dpsi[[s]][orth.per.tissue.age.qv[[s]]<0.05],0:50/50,border = NA,col='#FF000080',add=T,freq = F,main=names(born.tis.dpsi)[s])})

alt.devAS = sapply(names(orth.per.tissue.age.qv),function(s){apply(orth.per.tissue.age.qv[[s]][,colnames(orth.seg.ad.all.dpsi[[s]])]<0.05 & abs(orth.seg.ad.all.dpsi[[s]][rownames(orth.per.tissue.age.qv[[s]]),])>0.2,1,sum,na.rm=T)>0})
#sps=list(h='h',`m/r` = c('m','r'),hq='hq',mr='mr',mrb='mrb',hqmrb='hqmrb',hqmrbo='hqmrbo',hqmrboc='hqmrboc')
altern.prop.sgn.dpsi0.2=getSgnFraqWithBootstrap(alt.sp[alt.sp!=''],sps.,alt.devAS[,-c(2,7)],N=500)


## %% 3

#sps=list(`m/r` = c('m','r'),mr='mr',mrb='mrb',hqmrb='hqmrb',hqmrbo='hqmrbo')
seg.len = unlist(lapply(setNames(born.exn.sajr,NULL),function(x)setNames(x$seg$length,rownames(x$seg))))
seg.cod = unlist(lapply(setNames(born.exn.sajr,NULL),function(x)setNames(x$seg$cod=='c',rownames(x$seg))))
born.exn.cod = sapply(sps,function(s){r=seg.cod[unique(unlist(born.sids[s]))];my.binom.test(sum(r),sum(!r))})
born.exn.len3 = sapply(sps,function(s){r=seg.len[unique(unlist(born.sids[s]))] %% 3 == 0;my.binom.test(sum(r),sum(!r))})

#orth.seg.ad.all = readRDS('Rdata/orth.seg.ad.all.Rdata')
seg.cod = readRDS('Rdata/paper.figures/alt.seg.cod.Rdata')#seg.cod = unlist(lapply(setNames(orth.seg.ad.all,NULL),function(x)setNames(x$seg$cod=='c',rownames(x$seg))))
#saveRDS(seg.cod,'Rdata/paper.figures/alt.seg.cod.Rdata')
seg.len = readRDS('Rdata/paper.figures/alt.seg.len.Rdata')#seg.len=unlist(lapply(setNames(orth.seg.ad.all,NULL),function(x)setNames(x$seg$length,rownames(x$seg))))
#saveRDS(seg.len,'Rdata/paper.figures/alt.seg.len.Rdata')
alt.exn.len3 = sapply(sps.,function(s){r=seg.len[names(alt.sp)[alt.sp %in% s]] %% 3 == 0;my.binom.test(sum(r),sum(!r))})
alt.exn.cod = sapply(sps.,function(s){r=seg.cod[names(alt.sp)[alt.sp %in% s]];my.binom.test(sum(r),sum(!r))})


### mean PSI
s2s = setNames(rownames(species),species$short)
orth.mean.psi = sapply(orth.seg.ad.all.tsm,function(x)apply(x,1,mean,na.rm=T))
alt.sp. = strsplit(alt.sp,'')
orth.mean.psi. = sapply(1:nrow(orth.mean.psi),function(i)mean(orth.mean.psi[i,s2s[alt.sp.[[i]]]],na.rm=T))

alt.psi.on.ev=sapply(sps.,function(s){
	t=orth.mean.psi.[alt.sp %in% s]
	t = t[!is.na(t)]
	m = mean(t)
	s = sd(t)/sqrt(length(t))*3
	c(mean=m,ci1=m-s,m+s)
})

born.mean.psi = apply(born.seg.ids,1,function(x){
	mean(sapply(1:length(x),function(i){
		if(is.na(x[i])) return(NA)
		mean(born.exn.sajr[[i]]$psi.tsm.al[x[i],],na.rm=T)
	}),na.rm=T)
})

brn.psi.on.ev=sapply(sps,function(s){
	t=born.mean.psi[sp.birth %in% s]
	t = t[!is.na(t)]
	m = mean(t)
	s = sd(t)/sqrt(length(t))*3
	c(mean=m,ci1=m-s,m+s)
})
spes2use = rownames(species)[c(-2,-7)]#c('mouse','rat')
#sps=list(`m/r` = c('m','r'),mr='mr',mrb='mrb',hqmrb='hqmrb',hqmrbo='hqmrbo')
#max.stage.cl. = max.stage.cl[grep(paste(substr(spes2use,1,3),collapse='|'),names(max.stage.cl))]


tab.alt = sapply(sps.,function(s){table(factor(min.stage.cl[alt.sp %in% s,spes2use],levels=rownames(meta.tsm.al)))})


# devAS
born.tis.dpsi = lapply(rownames(species),function(s)getAgeASchanges(born.exn.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(born.tis.dpsi) =  rownames(species)
born.devAS = matrix(NA,ncol=ncol(born.seg.ids),nrow=nrow(born.seg.ids))
colnames(born.devAS) = colnames(born.seg.ids)
for(s in colnames(born.seg.ids)){
	for(i in 1:nrow(born.seg.ids)){
		if(!is.na(born.seg.ids[i,s])){
			born.devAS[i,s] = sum(born.per.tissue.age.qv[[s]][born.seg.ids[i,s],colnames(born.tis.dpsi[[s]])]<0.05 & abs(born.tis.dpsi[[s]][born.seg.ids[i,s],])>0.2,na.rm=TRUE)>0
		}
	}
}
born.devAS = !is.na(born.devAS) & born.devAS
born.devAS. = born.seg.ids[!is.na(born.devAS) & born.devAS]
born.ex.prop.sgn.dpsi0.2=getSgnFraqWithBootstrap(sp.birth,sps,born.devAS[,-c(2,7)],N=500)


tab.exn.dev  = sapply(sps,function(x){table(factor(max.stage.cl[intersect(intersect(unlist(born.sids[x]),names(max.stage.cl)),born.devAS.)],levels=rownames(m)))})
tab.exn.ndev = sapply(sps,function(x){table(factor(max.stage.cl[  setdiff(intersect(unlist(born.sids[x]),names(max.stage.cl)),born.devAS.)],levels=rownames(m)))})

min.stage.cl. = min.stage.cl[rownames(alt.devAS),]
alt.sp. = alt.sp[rownames(alt.devAS)]
tab.alt.dev  = sapply(sps.,function(s){table(factor(min.stage.cl.[alt.sp. %in% s,spes2use][ alt.devAS[alt.sp. %in% s,spes2use]],levels=rownames(meta.tsm.al)))})
tab.alt.ndev = sapply(sps.,function(s){table(factor(min.stage.cl.[alt.sp. %in% s,spes2use][!alt.devAS[alt.sp. %in% s,spes2use]],levels=rownames(meta.tsm.al)))})

t = min.stage.cl.
for(i in 1:ncol(t)) t[,i] = paste(t[,i],colnames(t)[i])
z=t[alt.sp.=='hqmrbo',][alt.devAS[alt.sp. %in% 'hqmrbo',spes2use]]
z=strsplit(z,' ')
z=do.call(rbind,z[sapply(z,length)==3])
table(z[,1],z[,3])



stg = sapply(strsplit(rownames(tab.alt),' '),'[',2)
mm = meta.tsm.al[!(meta.tsm.al$tissue %in% c('cerebellum','kidney','ovary') | meta.tsm.al$mouse.stage %in% c('10.5','11.5','3dpb','9wpb')),]
s = unique(mm[,c('days','mouse.stage')])
s = s$mouse.stage[order(s$days)]


propOnAge = function(x){sapply(x,function(i)my.binom.test(i,sum(x)-i))}

exn.on.age = t(sapply(split.data.frame(tab.exn,stg),function(x)apply(x,2,sum)))[s,]
alt.on.age = t(sapply(split.data.frame(tab.alt,stg),function(x)apply(x,2,sum)))[s,]
poa = list()
poa$e1 = propOnAge(exn.on.age[,1])
poa$e2 = propOnAge(exn.on.age[,5])
poa$a1 = propOnAge(alt.on.age[,1])
poa$a2 = propOnAge(alt.on.age[,5])


# _plot ####
pdf('figures/paper.figures/4/6.pdf',w=9,h=9)
par(mfrow=c(3,3),tck=-0.02,mgp=c(1.3,0.2,0),mar=c(3.5,2.5,1.5,0),oma=c(0,0,0,1))

barplot(sweep(tab.exn,2,apply(tab.exn,2,sum),'/'),col=meta.tsm.al$col,border = NA,las=3,main='New exons',ylab='proportion of exons mostly included in')
plotPanelLetter('A')

barplot(sweep(tab.alt,2,apply(tab.alt,2,sum),'/'),col=meta.tsm.al$col,border = NA,las=3,main='Alternification',ylab='proportion of exons mostly skipped in')
plotPanelLetter('B')
plot.new()
t=params$tissue.col[c(1,3,5,7)]
legend('topleft',fill=t,legend=names(t))

x = 1:ncol(poa[[1]])
cols = c('#FF7777','#FF0000','#7777FF','#0000FF')
lty = c(2,1,2,1)


plot(1,t='n',xaxt='n',xlab='Age (mouse)',main='New AS exons on development',ylab='proportion of alt. exons',ylim=range(unlist(poa)),xlim=range(x))
for(i in 1:length(poa)){
	lines(x,poa[[i]][1,],col=cols[i],lwd=3,lty=lty[i])
	segments(x,poa[[i]][2,],x,poa[[i]][3,],col=cols[i],lty=lty[i])
}
axis(1,x,rownames(exn.on.age),las=3)
legend('topleft',col=c('#FF7777','#FF0000','#7777FF','#0000FF'),lty=c(2,1,2,1),legend=c('recent new exon','mammal new exon','recent alernification','mammal alernification'))
plotPanelLetter('C')


poab = list()
poab$e1 = propOnAge(tab.exn[paste('brain',s),1])
poab$e2 = propOnAge(tab.exn[paste('brain',s),5])
poab$a1 = propOnAge(tab.alt[paste('brain',s),1])
poab$a2 = propOnAge(tab.alt[paste('brain',s),5])

plot(1,t='n',xaxt='n',xlab='Age (mouse)',main='New AS exons on development in brain',ylab='proportion of alt. exons',ylim=range(unlist(poab)),xlim=range(x))
for(i in 1:length(poab)){
	lines(x,poab[[i]][1,],col=cols[i],lwd=3,lty=lty[i])
	segments(x,poab[[i]][2,],x,poab[[i]][3,],col=cols[i],lty=lty[i])
}
axis(1,x,rownames(exn.on.age),las=3)

plotPanelLetter('D')

x=1:length(sps)
plotArea(x,t(born.exn.len3[,x]),col='red',new = T,type='b',lwd=3,ylim=range(born.exn.len3[,x],alt.exn.len3),xlim=c(1,length(sps.)),xaxt='n',xlab='',main='3N exons',ylab='proportion of 3*n exons')
segments(x,born.exn.len3[2,x],x,born.exn.len3[3,x],col='red')
plotArea(x,t(alt.exn.len3[,x]),col='blue',type='b',lwd=3)
segments(x,alt.exn.len3[2,x],x,alt.exn.len3[3,x],col='blue')
points(length(sps.),alt.exn.len3[1,length(sps.)],pch=19,col='orange')
segments(length(sps.),alt.exn.len3[2,length(sps.)],length(sps.),alt.exn.len3[3,length(sps.)],col='orange')
axis(1,1:ncol(alt.exn.len3),colnames(alt.exn.len3),las=3)
abline(h=1/3,lty=2)
plotPanelLetter('E')
legend('topleft',fill=c('red','blue'),legend=c('new exons','alternification'))


plotArea(x,t(born.exn.cod[,x]),col='red',new = T,type='b',lwd=3,ylim=range(born.exn.cod[,x],alt.exn.cod),xlim=c(1,length(sps.)),xaxt='n',xlab='',main='Coding (in Ensebml) exons',ylab='proportion of "coding" exons')
segments(x,born.exn.cod[2,x],x,born.exn.cod[3,x],col='red')
plotArea(x,t(alt.exn.cod[,x]),col='blue',type='b',lwd=3)
segments(x,alt.exn.cod[2,x],x,alt.exn.cod[3,x],col='blue')
points(length(sps.),alt.exn.cod[1,length(sps.)],pch=19,col='orange')
segments(length(sps.),alt.exn.cod[2,length(sps.)],length(sps.),alt.exn.cod[3,length(sps.)],col='orange')
axis(1,1:ncol(alt.exn.len3),colnames(alt.exn.len3),las=3)
plotPanelLetter('F')

plotArea(1:length(sps),p = born.ex.prop.sgn.dpsi0.2[,-2],col='red',lwd=3,new = T,xlim=c(1,length(sps.)),ylim=range(0,born.ex.prop.sgn.dpsi0.2[,-2],altern.prop.sgn.dpsi0.2[,-2]),xaxt='n',xlab='',ylab='proportion of devAS',main='devAS')
plotArea(1:length(sps),p = altern.prop.sgn.dpsi0.2[-length(sps.),-2],col='blue',lwd=3,new = F)
points(length(sps.),altern.prop.sgn.dpsi0.2[length(sps.),1],pch=19,col='orange')
segments(length(sps.),altern.prop.sgn.dpsi0.2[length(sps.),3],length(sps.),altern.prop.sgn.dpsi0.2[length(sps.),4],col='orange')
axis(1,1:length(sps.),colnames(alt.exn.len3),las=3)
abline(h=0,lty=2)
plotPanelLetter('G')

plotArea(1:length(sps),t(alt.psi.on.ev[,-length(sps.)]),col='blue',new=T,ylim=c(0,1),lwd=3,xaxt='n',xlab='',ylab='mean(PSI)',main='mean(PSI)',xlim=c(1,length(sps.)))
plotArea(1:length(sps),t(brn.psi.on.ev[,-length(sps.)]),col='red',lwd=3,new=F)
points(length(sps.),alt.psi.on.ev[1,length(sps.)],pch=19,col='orange')
segments(length(sps.),alt.psi.on.ev[2,length(sps.)],length(sps.),alt.psi.on.ev[3,length(sps.)],col='orange')
axis(1,1:length(sps.),colnames(alt.exn.len3),las=3)
plotPanelLetter('H')
dev.off()


# check newborn dependence on threshouls ####
load('Rdata/tmp.exon.birth.Rdata')
ier.max = unlist(lapply(exon.birth.one,function(x)(max(x$int.cov/x$exon.cov,na.rm=T))))
len =     unlist(lapply(exon.birth.one,function(x)(min(x$length,na.rm=T))))
p = unlist(setNames(lapply(born.exn.sajr,function(x)apply(x$ir[,colnames(x$ir) %in% rownames(meta)],1,function(z)sum(z>0.2,na.rm=T))),NULL))
high.psi = sapply(exon.birth.one,function(x)sum(p[x$seg_id[!is.na(x$seg_id)]]>3))

obs.sp = sapply(exon.birth.one,function(x)paste(species$short[rownames(species) %in% rownames(x)[!is.na(x$seg_id)]],collapse = ''))
table(names(obs.sp) == names(blcl))
s = species[colnames(nb.stat)[1:7],'short']
s=!t(sapply(strsplit(obs.sp,''),function(x){s %in% x}))
fisher.test(table(s,nb.stat[,1:7]>0))
nb.missed.due.n = s & nb.stat[,1:7]>0

p = unlist(setNames(lapply(born.exn.sajr,function(x)apply(x$ir[,colnames(x$ir) %in% rownames(meta)],1,function(z)sum(z>0.2,na.rm=T))),NULL))
high.psi = sapply(exon.birth.one,function(x)sum(p[x$seg_id[!is.na(x$seg_id)]]>3))

pt = unlist(setNames(lapply(born.exn.sajr,function(x)apply(x$ir[,colnames(x$ir) %in% rownames(meta)[meta$tissue=='testis']],1,function(z)sum(z>0.2,na.rm=T))),NULL))
high.psi.t = sapply(exon.birth.one,function(x)sum(pt[x$seg_id[!is.na(x$seg_id)]]>1))
table(high.psi>0,testis=high.psi.t>0)

filters = list(all = rep(TRUE,length(exon.birth.one)),
							 orth                      = nb.stat$orthn == 7,
							 orth.adj                  = nb.stat$orthn == 7 & nb.stat$adj.exons==0,
							 orth.adj.bl               = nb.stat$orthn == 7 & nb.stat$adj.exons==0 & nb.stat$full.obs==1,
							 orth.adj.bl.int           = nb.stat$orthn == 7 & nb.stat$adj.exons==0 & nb.stat$full.obs==1 & ier.max<0.5,
							 orth.adj.bl.int.len       = nb.stat$orthn == 7 & nb.stat$adj.exons==0 & nb.stat$full.obs==1 & ier.max<0.5 & len < 500,
							 orth.adj.bl.int.len.psi   = nb.stat$orthn == 7 & nb.stat$adj.exons==0 & nb.stat$full.obs==1 & ier.max<0.5 & len < 500 & (high.psi > 0 | high.psi.t > 0),
							 orth.adj.bl.int.len.psi.n = nb.stat$orthn == 7 & nb.stat$adj.exons==0 & nb.stat$full.obs==1 & ier.max<0.5 & len < 500 & (high.psi > 0 | high.psi.t > 0) & apply(nb.missed.due.n,1,sum)==0)
sapply(filters,sum)



sps=list(`m/r` = c('m','r'),mr='mr',mrb='mrb',hqmrb='hqmrb',hqmrbo='hqmrbo')
# sps=list(`h/q` = c('h','q'),hq='hq',hqmrb='hqmrb',hqmrbo='hqmrbo')
# sps=list(one = species$short,two=c('hq','mr'),mrb='mrb',hqmrb='hqmrb',hqmrbo='hqmrbo')
# sps=list(one = species$short,two=c('hq','mr'),hqmrbo='hqmrbo')
# sps=list(one = species$short,two=c('hq','mr'),more=c('mrb','hqmrb','hqmrbo'))

propOnAge = function(x){sapply(x,function(i)my.binom.test(i,sum(x)-i))}
s = unique(mm[,c('days','mouse.stage')])
s = s$mouse.stage[order(s$days)]

seg.len.nb = unlist(lapply(setNames(born.exn.sajr,NULL),function(x)setNames(x$seg$length,rownames(x$seg))))
seg.cod.nb = unlist(lapply(setNames(born.exn.sajr,NULL),function(x)setNames(x$seg$cod=='c',rownames(x$seg))))

pdf('figures/paper.figures/4/6(14).on.newborn.thrs.pdf',w=20,h=20)
par(mfrow=c(8,8),tck=-0.02,mgp=c(1.3,0.2,0),mar=c(3.5,2.5,1.5,0),oma=c(0,0,0,1))
for(n in names(filters)){
	f = filters[[n]]
	#main tissue
	born.sids = lapply(split.data.frame(born.seg.ids[f,,drop=F],sp.birth[f]),function(x){x[!is.na(x)]})
	tab.exn = sapply(sps, function(x){table(factor(max.stage.cl[intersect(unlist(born.sids[x]),names(max.stage.cl))],levels=rownames(meta.tsm.al)))})
	barplot(sweep(tab.exn,2,apply(tab.exn,2,sum),'/'),col=meta.tsm.al$col,border = NA,las=3,main=paste0(n,' (',sum(f),')'),ylab='proportion of exons mostly included in')
	#apply(tab.exn,2,sum)
	# plot only brain and testis
	bt = sapply(split.data.frame(tab.exn,sapply(strsplit(rownames(tab.exn),' '),'[',1)),function(x)apply(x,2,sum))
	bt = rbind(`mr-`=bt[1,]+bt[2,],`mrb+`=bt[3,]+bt[4,]+bt[5,])
	b = apply(bt,1,function(x)my.binom.test(x[1],sum(x[-1])))
	t = apply(bt,1,function(x)my.binom.test(x[7],sum(x[-7])))
	c = barplot(rbind(b[1,],t[1,]),beside = T,col=params$tissue.col[c('brain','testis')],border = NA,las=3,ylab='proportion of exons mostly included in',ylim=range(0,t,b))
	segments(c,rbind(b[2,],t[2,]),c,rbind(b[3,],t[3,]))
	
	# main stage
	stg = sapply(strsplit(rownames(tab.exn),' '),'[',2)
	s = unique(mm[,c('days','mouse.stage')])
	s = s$mouse.stage[order(s$days)]
	exn.on.age = t(sapply(split.data.frame(tab.exn,stg),function(x)apply(x,2,sum)))[s,]
	poa = list()
	poa$e1 = propOnAge(exn.on.age[,1])
	poa$e2 = propOnAge(exn.on.age[,length(sps)])
	x = 1:ncol(poa[[1]])
	cols = c('#FF7777','#FF0000','#7777FF','#0000FF')
	lty = c(2,1,2,1)
	

	plot(1,t='n',xaxt='n',xlab='Age (mouse)',main='New AS exons on development',ylab='proportion of alt. exons',ylim=range(unlist(poa)),xlim=range(x))
	for(i in 1:length(poa)){
		lines(x,poa[[i]][1,],col=cols[i],lwd=3,lty=lty[i])
		segments(x,poa[[i]][2,],x,poa[[i]][3,],col=cols[i],lty=lty[i])
	}
	axis(1,x,rownames(exn.on.age),las=3)
	legend('topleft',col=c('#FF7777','#FF0000'),lty=c(2,1,2,1),legend=names(sps)[c(1,length(sps))])
	#main stage brain only
	poab = list()
	poab$e1 = propOnAge(tab.exn[paste('brain',s),1])
	poab$e2 = propOnAge(tab.exn[paste('brain',s),length(sps)])

	plot(1,t='n',xaxt='n',xlab='Age (mouse)',main='New AS exons on development in brain',ylab='proportion of alt. exons',ylim=range(unlist(poab)),xlim=range(x))
	for(i in 1:length(poab)){
		lines(x,poab[[i]][1,],col=cols[i],lwd=3,lty=lty[i])
		segments(x,poab[[i]][2,],x,poab[[i]][3,],col=cols[i],lty=lty[i])
	}
	axis(1,x,rownames(exn.on.age),las=3)
	# 3N exons
	born.exn.len3 = sapply(sps,function(s){r=seg.len.nb[unique(unlist(born.sids[s]))] %% 3 == 0;my.binom.test(sum(r),sum(!r))})
	x=1:length(sps)
	plotArea(x,t(born.exn.len3[,x]),col='red',new = T,type='b',lwd=3,ylim=range(born.exn.len3[,x],alt.exn.len3[3,ncol(alt.exn.len3)]),xlim=c(1,length(sps)+1),xaxt='n',xlab='',main='3N exons',ylab='proportion of 3*n exons')
	segments(x,born.exn.len3[2,x],x,born.exn.len3[3,x],col='red')
	points(length(sps)+1,alt.exn.len3[1,length(sps.)],pch=19,col='orange')
	segments(length(sps)+1,alt.exn.len3[2,ncol(alt.exn.len3)],length(sps)+1,alt.exn.len3[3,length(sps.)],col='orange')
	axis(1,1:(length(sps)+1),c(names(sps),colnames(alt.exn.len3)[ncol(alt.exn.len3)]),las=3)
	abline(h=1/3,lty=2)
	# coding 
	born.exn.cod = sapply(sps,function(s){r=seg.cod.nb[unique(unlist(born.sids[s]))];my.binom.test(sum(r),sum(!r))})
	plotArea(x,t(born.exn.cod[,x]),col='red',new = T,type='b',lwd=3,ylim=range(born.exn.cod[,x],alt.exn.cod[3,ncol(alt.exn.cod)]),xlim=c(1,length(sps)+1),xaxt='n',xlab='',main='Coding (in Ensebml) exons',ylab='proportion of "coding" exons')
	segments(x,born.exn.cod[2,x],x,born.exn.cod[3,x],col='red')
	points(length(sps)+1,alt.exn.cod[1,ncol(alt.exn.cod)],pch=19,col='orange')
	segments(length(sps)+1,alt.exn.cod[2,ncol(alt.exn.cod)],length(sps)+1,alt.exn.cod[3,ncol(alt.exn.cod)],col='orange')
	axis(1,1:(length(sps)+1),c(names(sps),colnames(alt.exn.len3)[ncol(alt.exn.cod)]),las=3)

	# devAS
	born.ex.prop.sgn.dpsi0.2=getSgnFraqWithBootstrap(sp.birth[f],sps,born.devAS[f,],N=500)
	
	
	plotArea(1:length(sps),p = born.ex.prop.sgn.dpsi0.2[,-2],col='red',lwd=3,new = T,xlim=c(1,length(sps)+1),ylim=range(0,born.ex.prop.sgn.dpsi0.2[,-2],altern.prop.sgn.dpsi0.2[,-2]),xaxt='n',xlab='',ylab='proportion of devAS',main='devAS')
	points(length(sps)+1,altern.prop.sgn.dpsi0.2[nrow(altern.prop.sgn.dpsi0.2),1],pch=19,col='orange')
	segments(length(sps)+1,altern.prop.sgn.dpsi0.2[nrow(altern.prop.sgn.dpsi0.2),3],length(sps)+1,altern.prop.sgn.dpsi0.2[nrow(altern.prop.sgn.dpsi0.2),4],col='orange')
	axis(1,1:(length(sps)+1),c(names(sps),rownames(altern.prop.sgn.dpsi0.2)[nrow(altern.prop.sgn.dpsi0.2)]),las=3)
	abline(h=0,lty=2)

	# mean PSI
	brn.psi.on.ev=sapply(sps,function(s){
		t=born.mean.psi[f & sp.birth %in% s]
		t = t[!is.na(t)]
		m = mean(t)
		s = sd(t)/sqrt(length(t))*3
		c(mean=m,ci1=m-s,m+s)
	})
	
	
	plotArea(1:length(sps),t(brn.psi.on.ev[,-length(sps.)]),col='red',lwd=3,new=T,ylim=c(0,1),xaxt='n',xlab='',ylab='mean(PSI)',main='mean(PSI)',xlim=c(1,length(sps)+1))
	points(length(sps)+1,alt.psi.on.ev[1,ncol(alt.psi.on.ev)],pch=19,col='orange')
	segments(length(sps)+1,alt.psi.on.ev[2,ncol(alt.psi.on.ev)],length(sps)+1,alt.psi.on.ev[3,ncol(alt.psi.on.ev)],col='orange')
	axis(1,1:(length(sps)+1),c(names(sps),colnames(alt.psi.on.ev)[ncol(alt.psi.on.ev)]),las=3)
	
}
dev.off()

# alternif on threshould #####
getAltSp = function(p,thr,n){
	apply(sapply(p,function(x){apply(x<thr,1,sum,na.rm=T)>n}),1,function(x)paste(species$short[x],collapse=''))
}

alternif = list(by.ann = alt.sp,
								psi95.2 = getAltSp(orth.seg.ad.all.tsm,0.95,2),
								psi95.4 = getAltSp(orth.seg.ad.all.tsm,0.95,4),
								psi95.8 = getAltSp(orth.seg.ad.all.tsm,0.95,8),
								psi90.2 = getAltSp(orth.seg.ad.all.tsm,0.9 ,2),
								psi90.4 = getAltSp(orth.seg.ad.all.tsm,0.9 ,4),
								psi90.8 = getAltSp(orth.seg.ad.all.tsm,0.9 ,8),
								psi80.2 = getAltSp(orth.seg.ad.all.tsm,0.8 ,2),
								psi80.4 = getAltSp(orth.seg.ad.all.tsm,0.8 ,4))
sps=list(`m/r` = c('m','r'),mr='mr',mrb='mrb',hqmrb='hqmrb',hqmrbo='hqmrbo')

propOnAge = function(x){sapply(x,function(i)my.binom.test(i,sum(x)-i))}
s = unique(mm[,c('days','mouse.stage')])
s = s$mouse.stage[order(s$days)]

seg.len.nb = unlist(lapply(setNames(born.exn.sajr,NULL),function(x)setNames(x$seg$length,rownames(x$seg))))
seg.cod.nb = unlist(lapply(setNames(born.exn.sajr,NULL),function(x)setNames(x$seg$cod=='c',rownames(x$seg))))




pdf('figures/paper.figures/4/6(14).on.alternif.thrs.pdf',w=20,h=20/8*length(alternif))
par(mfrow=c(length(alternif),8),tck=-0.02,mgp=c(1.3,0.2,0),mar=c(3.5,2.5,1.5,0),oma=c(0,0,0,1))
for(n in names(alternif)){
	#main tissue
	tab.alt = sapply(sps.,function(s){table(factor(min.stage.cl[alternif[[n]] %in% s,],levels=rownames(meta.tsm.al)))})
	barplot(sweep(tab.alt,2,apply(tab.alt,2,sum),'/'),col=meta.tsm.al$col,border = NA,las=3,main=n,ylab='proportion of exons mostly skipped in')
	
	# plot only brain and testis
	bt = sapply(split.data.frame(tab.alt,sapply(strsplit(rownames(tab.alt),' '),'[',1)),function(x)apply(x,2,sum))
	bt = rbind(`mr-`=bt[1,]+bt[2,],`mrb+`=bt[3,]+bt[4,]+bt[5,])
	b = apply(bt,1,function(x)my.binom.test(x[1],sum(x[-1])))
	t = apply(bt,1,function(x)my.binom.test(x[7],sum(x[-7])))
	c = barplot(rbind(b[1,],t[1,]),beside = T,col=params$tissue.col[c('brain','testis')],border = NA,las=3,ylab='proportion of exons mostly skipped in',ylim=range(0,t,b))
	segments(c,rbind(b[2,],t[2,]),c,rbind(b[3,],t[3,]))
	
	# main stage
	stg = sapply(strsplit(rownames(tab.alt),' '),'[',2)
	s = unique(mm[,c('days','mouse.stage')])
	s = s$mouse.stage[order(s$days)]
	alt.on.age = t(sapply(split.data.frame(tab.alt,stg),function(x)apply(x,2,sum)))[s,]
	poa = list()
	poa$a1 = propOnAge(alt.on.age[,1])
	poa$a2 = propOnAge(alt.on.age[,length(sps)])
	x = 1:ncol(poa[[1]])
	cols = c('#7777FF','#0000FF')
	lty = c(2,1,2,1)
	
	
	plot(1,t='n',xaxt='n',xlab='Age (mouse)',main='Alternification on development',ylab='proportion of alt. exons',ylim=range(unlist(poa)),xlim=range(x))
	for(i in 1:length(poa)){
		lines(x,poa[[i]][1,],col=cols[i],lwd=3,lty=lty[i])
		segments(x,poa[[i]][2,],x,poa[[i]][3,],col=cols[i],lty=lty[i])
	}
	axis(1,x,rownames(alt.on.age),las=3)
	legend('topleft',col=c('#7777FF','#FF0000'),lty=c(2,1,2,1),legend=names(sps)[c(1,length(sps))])
	
	#main stage brain only
	poab = list()
	poab$e1 = propOnAge(tab.alt[paste('brain',s),1])
	poab$e2 = propOnAge(tab.alt[paste('brain',s),length(sps)])
	
	plot(1,t='n',xaxt='n',xlab='Age (mouse)',main='Alternification on development in brain',ylab='proportion of alt. exons',ylim=range(unlist(poab)),xlim=range(x))
	for(i in 1:length(poab)){
		lines(x,poab[[i]][1,],col=cols[i],lwd=3,lty=lty[i])
		segments(x,poab[[i]][2,],x,poab[[i]][3,],col=cols[i],lty=lty[i])
	}
	axis(1,x,rownames(alt.on.age),las=3)
	
	# 3N exons
	alt.exn.len3 = sapply(sps.,function(s){r=seg.len[names(alternif[[n]])[alternif[[n]] %in% s]] %% 3 == 0;my.binom.test(sum(r),sum(!r))})
	x=1:length(sps)
	plotArea(x,t(alt.exn.len3[,x]),col='blue',new = T,type='b',lwd=3,ylim=range(alt.exn.len3),xlim=c(1,length(sps.)),xaxt='n',xlab='',main='3N exons',ylab='proportion of 3*n exons')
	segments(x,alt.exn.len3[2,x],x,alt.exn.len3[3,x],col='blue')
	points(length(sps)+1,alt.exn.len3[1,length(sps.)],pch=19,col='orange')
	segments(length(sps)+1,alt.exn.len3[2,ncol(alt.exn.len3)],length(sps)+1,alt.exn.len3[3,length(sps.)],col='orange')
	axis(1,1:(length(sps)+1),c(names(sps),colnames(alt.exn.len3)[ncol(alt.exn.len3)]),las=3)
	abline(h=1/3,lty=2)
	# coding 
	alt.exn.cod = sapply(sps.,function(s){r=seg.cod[names(alternif[[n]])[alternif[[n]] %in% s]];my.binom.test(sum(r),sum(!r))})
	plotArea(x,t(alt.exn.cod[,x]),col='blue',new = T,type='b',lwd=3,ylim=range(alt.exn.cod),xlim=c(1,length(sps)+1),xaxt='n',xlab='',main='Coding (in Ensebml) exons',ylab='proportion of "coding" exons')
	segments(x,alt.exn.cod[2,x],x,alt.exn.cod[3,x],col='blue')
	points(length(sps)+1,alt.exn.cod[1,ncol(alt.exn.cod)],pch=19,col='orange')
	segments(length(sps)+1,alt.exn.cod[2,ncol(alt.exn.cod)],length(sps)+1,alt.exn.cod[3,ncol(alt.exn.cod)],col='orange')
	axis(1,1:(length(sps)+1),c(names(sps),colnames(alt.exn.len3)[ncol(alt.exn.cod)]),las=3)
	 
	# devAS
	ff = intersect(names(alternif[[n]])[alternif[[n]]!=''],rownames(alt.devAS)) # age test was performed only for exonas that are alt in at least one species
	altern.prop.sgn.dpsi0.2 = getSgnFraqWithBootstrap(alternif[[n]][ff],sps.,alt.devAS[ff,-c(2,7)],N=500)
	
	plotArea(1:length(sps),p = altern.prop.sgn.dpsi0.2[,-2],col='blue',lwd=3,new = T,xlim=c(1,length(sps)+1),ylim=range(0,altern.prop.sgn.dpsi0.2[,-2]),xaxt='n',xlab='',ylab='proportion of devAS',main='devAS')
	points(length(sps)+1,altern.prop.sgn.dpsi0.2[nrow(altern.prop.sgn.dpsi0.2),1],pch=19,col='orange')
	segments(length(sps)+1,altern.prop.sgn.dpsi0.2[nrow(altern.prop.sgn.dpsi0.2),3],length(sps)+1,altern.prop.sgn.dpsi0.2[nrow(altern.prop.sgn.dpsi0.2),4],col='orange')
	axis(1,1:(length(sps)+1),c(names(sps),rownames(altern.prop.sgn.dpsi0.2)[nrow(altern.prop.sgn.dpsi0.2)]),las=3)
	abline(h=0,lty=2)

	# mean PSI
	alt.psi.on.ev=sapply(sps.,function(s){
		t=orth.mean.psi.[alternif[[n]] %in% s]
		t = t[!is.na(t)]
		m = mean(t)
		s = sd(t)/sqrt(length(t))*3
		c(mean=m,ci1=m-s,m+s)
	})

	plotArea(1:length(sps),t(alt.psi.on.ev[,-length(sps.)]),col='blue',lwd=3,new=T,ylim=c(0,1),xaxt='n',xlab='',ylab='mean(PSI)',main='mean(PSI)',xlim=c(1,length(sps)+1))
	points(length(sps)+1,alt.psi.on.ev[1,ncol(alt.psi.on.ev)],pch=19,col='orange')
	segments(length(sps)+1,alt.psi.on.ev[2,ncol(alt.psi.on.ev)],length(sps)+1,alt.psi.on.ev[3,ncol(alt.psi.on.ev)],col='orange')
	axis(1,1:(length(sps)+1),c(names(sps),colnames(alt.psi.on.ev)[ncol(alt.psi.on.ev)]),las=3)
}
dev.off()

# look on exon loss ####
gain = c(species$short,'hq','mr','mrb','hqmrb','hqmrbo')
loss = sapply(gain[-12],function(x)gsub(x,'','hqmrboc'))
el = sapply(filters,function(f)table(factor(obs.sp)[f]),xlab='',ylab='')
plotArea(1:ncol(el),t(apply(cbind(el['mrb',],apply(el,2,sum)),1,my.binom.test)),col='red',new = T,ylim=c(0,0.01))
plotArea(1:ncol(el),t(apply(cbind(el['hmr',],apply(el,2,sum)),1,my.binom.test)),col='blue',new = F)
plotArea(1:ncol(el),t(apply(cbind(el['qm',],apply(el,2,sum)),1,my.binom.test)),col='gray',new = F)

plotArea(1:ncol(el),t(apply(cbind(el['hqmboc',],apply(el,2,sum)),1,my.binom.test)),col='red',new = T)
plot(sweep(el,2,apply(el,2,sum),'/')['mr',])
el['mrboc',]

pdf('figures/newborn/species.prop.on.threshoulds.pdf',w=12,h=8)
fnames=c('all','o','oa','oab','oabi','oabil','oabilp','oabilpn')
par(mfrow=c(4,6),tck=-0.02,mgp=c(1.3,0.2,0),mar=c(5.5,2.5,1.5,0),oma=c(0,0,0,1))
for(g in gain){
	plotArea(1:ncol(el),t(apply(cbind(el[g,],apply(el,2,sum)),1,my.binom.test)),col='red',new = T,main=g,xlab='',xaxt='n',ylab='proportion')
	axis(1,1:ncol(el),paste0(fnames,' (',el[g,],')'),las=3)
}

for(l in loss){
	plotArea(1:ncol(el),t(apply(cbind(el[l,],apply(el,2,sum)),1,my.binom.test)),col='blue',new = T,main=l,xlab='',xaxt='n',ylab='proportion')
	axis(1,1:ncol(el),paste0(fnames,' (',el[l,],')'),las=3)
}

par(mfrow=c(4,6),tck=-0.02,mgp=c(1.3,0.2,0),mar=c(5.5,2.5,1.5,0),oma=c(0,0,0,1))
incons = setdiff(rownames(el),c(loss,gain))
incons = incons[order(apply(el[incons,],1,sum),decreasing = T)]
for(l in incons){
	plotArea(1:ncol(el),t(apply(cbind(el[l,],apply(el,2,sum)),1,my.binom.test)),col='gray',new = T,main=l,xlab='',xaxt='n',ylab='proportion')
	axis(1,1:ncol(el),paste0(fnames,' (',el[l,],')'),las=3)
}
dev.off()


pdf('figures/newborn/loss-n-gain.cnt.on.threshoulds.pdf',w=10,h=5)
par(mfrow=c(1,2),tck=-0.02,mgp=c(2.3,0.4,0),mar=c(3.5,3.5,1.5,0),oma=c(0,0,0,1))
t=t(el[rev(gain),])
rownames(t) = fnames
f = sweep(t,1,apply(el,2,sum),'/')
imageWithText(sweep(f,2,apply(f,2,sum),'/'),t,names.as.labs = T,xlab='threshould',ylab='observed species',main='exon gain')

t=t(el[rev(loss),])
rownames(t) = fnames
colnames(t) = rev(gain[-12])
f = sweep(t,1,apply(el,2,sum),'/')
imageWithText(sweep(f,2,apply(f,2,sum),'/'),t,names.as.labs = T,xlab='threshould',ylab='lost in',main='exon loss')
dev.off()




# check microexons in newborn ##### 
born.exn.sajr = readRDS('Rdata/born.exn.sajr.Rdata')
sp.birth = sapply(exon.birth.one,function(x){paste(species$short[!is.na(x$seg_id)],collapse='')})
lens=lapply(exon.birth.one,function(x)x$length[!is.na(x$length)])
hist(unlist(lens),seq(0,100000,by=1),xlim=c(0,50),col=c('gray','darkgray','lightgreen'))
hist(sapply(lens,min),seq(0,100000,by=1),xlim=c(0,50),col=c('gray','darkgray','lightgreen'))
f=sapply(lens,min) < 28
sort(table(sp.birth[f]))

psi.mean=unlist(lapply(born.exn.sajr,function(x){apply(x$ir,1,mean,na.rm=T)}))
psi.mean=t(sapply(exon.birth.one,function(x)psi.mean[paste0(rownames(x),'.',x$seg_id)]))

boxplot(apply(psi.mean,1,mean,na.rm=T) ~ nchar(sp.birth))
boxplot(apply(psi.mean,1,max,na.rm=T) ~ nchar(sp.birth))

psi.mean[apply(is.na(psi.mean),1,sum)==7,]
f = apply(is.na(psi.mean),1,sum)<7 & apply(psi.mean,1,max,na.rm=T) > 0.2
table(f)
hist(sapply(lens[f],min),seq(0,100000,by=1),xlim=c(0,150),col=c('gray','darkgray','lightgreen'))
sort(table(sp.birth[f & sapply(lens,min) < 28]))
exon.birth.one[f & sapply(lens,min) < 28 & sp.birth=='hqmrb']



#t553 is always included (except rabbit)

# SNAP91 ENSG00000065609 synaptosomal-associated protein, 6nt
par(mfrow=c(3,3),tck=-0.02,mgp=c(1.3,0.2,0),mar=c(3.5,2.5,1.5,0),oma=c(0,0,0,1))
for(i in 1:5)
	plotTissueAgeProile(born.exn.sajr[[i]]$ir[exon.birth.one[['t3475']]$seg_id[i],],meta,main=rownames(species)[i],ylim=c(0,1))
hgmd[hgmd$chrom_VCF_hg19=='6' & hgmd$pos_VCF_hg19 >= (84324576 - 150) & hgmd$pos_VCF_hg19 <= (84324581 + 150),]
# NSFL1C ENSG00000088833 NSFL1 (p97) cofactor (p47) 6nt
par(mfrow=c(3,3),tck=-0.02,mgp=c(1.3,0.2,0),mar=c(3.5,2.5,1.5,0),oma=c(0,0,0,1))
for(i in 1:5)
	plotTissueAgeProile(born.exn.sajr[[i]]$ir[exon.birth.one[['t2143']]$seg_id[i],],meta,main=rownames(species)[i],ylim=c(0,1))
hgmd[hgmd$chrom_VCF_hg19=='21' & hgmd$pos_VCF_hg19 >= (1436359 - 150) & hgmd$pos_VCF_hg19 <= (1436364 + 150),]

# APP but it is not microexon
par(mfrow=c(3,3),tck=-0.02,mgp=c(1.3,0.2,0),mar=c(3.5,2.5,1.5,0),oma=c(0,0,0,1))
for(i in 1:5)
	plotTissueAgeProile(born.exn.sajr[[i]]$ir[exon.birth.one[['t2267']]$seg_id[i],],meta,main=rownames(species)[i],ylim=c(0,1))
hgmd[hgmd$chrom_VCF_hg19=='21' & hgmd$pos_VCF_hg19 >= (27369675 - 150) & hgmd$pos_VCF_hg19 <= (27369731 + 150),]

exon.birth.one[['t2267']]

# alternif on thrs fig 6 #####


# 2018.12.10 #####
# _3 four patterns ####
age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)

getDevASpattern = function(psi,sgn,dpsi,m,f){
	psi = psi[f,]
	sgn = sgn[f,]
	dpsi=dpsi[f,]
	sgn[is.na(sgn)] = FALSE
	m = m[colnames(psi),]
	r = data.frame(dir=character(nrow(psi)),mean.psi=NA)
	rownames(r) = rownames(psi)
	for(i in 1:nrow(psi)){
		dir = unique(sign(dpsi[i,sgn[i,]]))
		if(length(dir)==1)
			r$dir[i] = ifelse(dir>0,'u','d')
		else if(length(dir)==0)
			r$dir[i] = 'n'
		else
			r$dir[i] = 'b'
		r$mean.psi[i] = mean(psi[i,!(m$tissue %in% colnames(sgn)[sgn[i,]])],na.rm=T)
	}
	r
}

m03 = getDevASpattern(psi.tsm$mouse,abs(age.dpsi$mouse)>0.3 & per.tissue.age.qv$mouse<0.05,age.dpsi$mouse,meta.tsm,anns$mouse$sites=='ad')
m05 = getDevASpattern(psi.tsm$mouse,abs(age.dpsi$mouse)>0.5 & per.tissue.age.qv$mouse<0.05,age.dpsi$mouse,meta.tsm,anns$mouse$sites=='ad')
table(m03$dir,m03$mean.psi>0.5)
table(m05$dir,m05$mean.psi>0.5)

all03=sapply(rownames(species)[-2],function(s){
	x=getDevASpattern(psi.tsm[[s]],abs(age.dpsi[[s]])>0.3 & per.tissue.age.qv[[s]]<0.05,age.dpsi[[s]],meta.tsm,anns[[s]]$sites=='ad')
	c(u0=sum(x$dir=='u' & x$mean.psi<0.5,na.rm=T),u1=sum(x$dir=='u' & x$mean.psi>0.5,na.rm=T),d0=sum(x$dir=='d' & x$mean.psi<0.5,na.rm=T),d1=sum(x$dir=='d' & x$mean.psi>0.5,na.rm=T))
	})

all05=sapply(rownames(species)[-2],function(s){
	x=getDevASpattern(psi.tsm[[s]],abs(age.dpsi[[s]])>0.5 & per.tissue.age.qv[[s]]<0.05,age.dpsi[[s]],meta.tsm,anns[[s]]$sites=='ad')
	c(u0=sum(x$dir=='u' & x$mean.psi<0.5,na.rm=T),u1=sum(x$dir=='u' & x$mean.psi>0.5,na.rm=T),d0=sum(x$dir=='d' & x$mean.psi<0.5,na.rm=T),d1=sum(x$dir=='d' & x$mean.psi>0.5,na.rm=T))
})



plotDevASpatterns = function(d,...){
	b=barplot(t(d),beside=T,ylab='# of exons',xaxt='n',...)
	text(b,rep(0,length(b)),rep(species$short[-2],times=4),xpd=T,adj=c(0.5,1.1))
	
	y=-c(0.62,0.2) * par('usr')[4]
	segments(b[1,1],y[1],b[6,1],y[1],xpd=T,lwd=4)
	segments(b[1,1],y[1],b[6,1],y[2],xpd=T,lwd=2)
	
	segments(b[1,2],y[2],b[6,2],y[2],xpd=T,lwd=4)
	segments(b[1,2],y[1],b[6,2],y[2],xpd=T,lwd=2)
	
	segments(b[1,3],y[1],b[6,3],y[1],xpd=T,lwd=4)
	segments(b[1,3],y[2],b[6,3],y[1],xpd=T,lwd=2)
	
	segments(b[1,4],y[2],b[6,4],y[2],xpd=T,lwd=4)
	segments(b[1,4],y[2],b[6,4],y[1],xpd=T,lwd=2)
	text(-b[1,1],mean(y),'PSI\npattern',adj=c(0.5,0.5),xpd=T,srt=90)
}

pdf('figures/paper.figures/2018.12.10/03.devAS.4patterns.pdf',w=4,h=4)
par(mfrow=c(2,1),tck=-0.02,mgp=c(1.3,0.2,0),mar=c(3.5,2.5,1.5,0),oma=c(0,0,0,1))
plotDevASpatterns(all03,main='dPSI > 0.3')
plotDevASpatterns(all05,main='dPSI > 0.5')
dev.off()

# check AS in early/late ####
mc = read.csv(paste0('processed/GE.from.marg/MouseClusters.csv'),row.names = 1)
as.in.ge.patterns.mouse = readRDS('Rdata/paper.figures/as.in.ge.patterns.mouse.Rdata') # see paper.figures.R

x=lapply(names(as.in.ge.patterns.mouse),function(t){
	apply(table(mc[[paste0(firstToupper(t),'Pattern')]],rownames(mc) %in% rownames(as.in.ge.patterns.mouse[[t]]))[c('Decreasing','Increasing'),c('TRUE','FALSE')],1,my.binom.test)
})

pdf('figures/paper.figures/2018.12.10/04.proportion.of.AS.in.early-late.pdf',w=4,h=3)
par(tck=-0.02,mgp=c(1.3,0.2,0),mar=c(4.5,2.5,1.5,0),oma=c(0,0,0,1))
b=barplot(sapply(x,'[',1,1:2),beside = T,col=rep(params$tissue.col,each=2),den=c(-1,30),las=3,ylab='proportion of genes with AS',names.arg = names(as.in.ge.patterns.mouse))
segments(b,sapply(x,'[',2,1:2),b,sapply(x,'[',3,1:2))
dev.off()

# devAS cons on evo-age ###
dpsi = 0.2
sps = c('mr','mrb','hmrb','hmrbo','hmrboc')
species2comp = c('mouse','rat')
tis = 'brain'

getDevASConsOnEvolAge = function(dpsi,sps,species2comp,tis){
	orth.age.ad2 = lapply(rownames(species),function(s)getAgeASchanges(orth.seg.ad.tsm,meta.tsm,psi.thr = 0.2,border.stages,s))
	names(orth.age.ad2) = rownames(species)
	for(sp in names(orth.age.ad2)) orth.age.ad2[[sp]][orth.age.ad2[[sp]] != '-' & (is.na(orth.per.tissue.age.qv[[sp]]) | orth.per.tissue.age.qv[[sp]]>0.05)[,colnames(orth.age.ad2[[sp]])]] = 'n'
	
	orth.alt = sapply(orth.seg.ad,function(x)x$seg$type)[,-2]
	orth.alt = apply(orth.alt=='ALT',1,function(x)paste(species$short[-2][x],collapse=''))
	
	
	sgn = split.data.frame(cbind(orth.age.ad2[[species2comp[1]]][,tis],orth.age.ad2[[species2comp[2]]][,tis]),orth.alt)[sps]
	
	t(sapply(sgn,function(x){
		tot = nrow(x)
		tst = sum(x[,1] != '-' & x[,2] != '-')
		x = x[x[,1] %in% c('u','d'),]
		r=c(total=tot,tested=tst,cnt=nrow(x),cons=sum(x[,1]==x[,2]))
		c(r,my.binom.test(r[4],r[3]-r[4]))
	}))
}



cons.stat=lapply(unique(meta$tissue),function(t)getDevASConsOnEvolAge(0.2,c('mr','mrb','hmrb','hmrbo','hmrboc'), c('mouse','rat'),t))
names(cons.stat) = unique(meta$tissue)

pdf('figures/devAS/mouse-rat.devAS.cons.on.evo.age.dpsi=0.2.pdf',w=9,h=9)
par(mfrow=c(3,3),tck=-0.02,mgp=c(1.3,0.2,0),mar=c(4.5,2.5,1.5,0),oma=c(0,0,0,1))
for(t in names(cons.stat)){
	plotArea(1:nrow(cons.stat[[t]]),cons.stat[[t]][,5:7],col=params$tissue.col[t],new = T,main=t,ylab='proportion of mouse devAS conserved in rat',xaxt='n',ylim=c(0,1))
	axis(1,1:nrow(cons.stat[[t]]),rownames(cons.stat[[t]]))
}
dev.off()

par(mfrow=c(3,3),tck=-0.02,mgp=c(1.3,0.2,0),mar=c(4.5,2.5,1.5,0),oma=c(0,0,0,1))
for(t in names(cons.stat)){
	z=t(apply(cons.stat[[t]],1,function(x)my.binom.test(x[3],x[2]-x[3])))
	plotArea(1:nrow(cons.stat[[t]]),z,col=params$tissue.col[t],new = T,main=t,ylab='proportion of mouse devAS in tested',xaxt='n')
	axis(1,1:nrow(cons.stat[[t]]),rownames(cons.stat[[t]]))
}

# brain is more divergent? ####
alt.sp = readRDS('Rdata/paper.figures/alt.sp.Rdata')

s1 = 'mouse'
s2 = 'rabbit'


psi = cbind(orth.seg.ad.tsm[[s1]],orth.seg.ad.tsm[[s2]])
filters = data.frame(alt.sp = alt.sp[alt.sp!=''])
c1 = paste(s1,rep(unique(meta$tissue),times=nrow(age.al.i)),rep(age.al.i[[s1]],each=7))
c2 = paste(s2,rep(unique(meta$tissue),times=nrow(age.al.i)),rep(age.al.i[[s2]],each=7))
table(c1 %in% colnames(psi),c2 %in% colnames(psi))
f = c1 %in% colnames(psi) & c2 %in% colnames(psi)
psi = psi[,c(c1[f],c2[f])]
filters$na.cnt = apply(is.na(psi),1,sum)
#hist(filters$na.cnt)
table(filters$na.cnt==0)
filters$sd1 = apply(psi[,c1[f]],1,sd,na.rm=T)
filters$sd2 = apply(psi[,c2[f]],1,sd,na.rm=T)
plot(filters$sd1,filters$sd2)
# so there are exons expressed in all tissues that are ALT in both species
ff = filters$sd1 > 0.05 & filters$sd2 > 0.05 & filters$na.cnt == 0
ff = grepl('hqmrb',filters$alt.sp) & filters$na.cnt == 0
table(ff)
plotTissueAgeProile(setNames(sapply(which(f),function(i)cor(psi[ff,c1[i]],psi[ff,c2[i]],m='p',u='p')),c1[f]),meta.tsm,age.axis = 'rank')

plotLine(psi[ff,paste(s1,'brain',age.al.i[1,s1])],psi[ff,paste(s2,'brain',age.al.i[1,s2])])
plotLine(psi[ff,paste(s1,'brain',age.al.i[14,s1])],psi[ff,paste(s2,'brain',age.al.i[14,s2])])
plotLine(psi[ff,paste(s1,'liver',age.al.i[14,s1])],psi[ff,paste(s2,'liver',age.al.i[14,s2])])

hist(psi[ff,paste(s1,'heart',age.al.i[13,s1])])

plotLine(psi[,'mouse brain 9wpb'],psi[,'rat brain 16wpb'])

plotLine(psi[,'mouse testis 9wpb'],psi[,'rat testis 16wpb'])

# devAS species direction ####
alt.sp = readRDS('Rdata/paper.figures/alt.sp.Rdata')

orth.age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(orth.seg.ad.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(orth.age.dpsi) = rownames(species)


compDevAS = function(f,q,dp,s1,s2,dpt=0.1){
	q1 = q[[s1]][f,]
	q2 = q[[s2]][f,]			 
	d1 = dp[[s1]][f,]
	d2 = dp[[s2]][f,]
	sg1 = q1 < 0.05 & abs(d1) > dpt
	sg2 = q2 < 0.05 & abs(d2) > dpt
	data.frame(tested1 = apply(!is.na(sg1),2,sum),
				tested2 = apply(!is.na(sg2),2,sum),
				testedB = apply(!is.na(sg1) & !is.na(sg2),2,sum),
				sgn1 = apply(sg1,2,sum,na.rm=T),
				sgn2 = apply(sg2,2,sum,na.rm=T),
				sgnT1 = apply(!is.na(sg1) & !is.na(sg2) & sg1,2,sum,na.rm=T),
				sgnT2 = apply(!is.na(sg1) & !is.na(sg2) & sg2,2,sum,na.rm=T),
				sgnB = apply(sg1 & sg2,2,sum,na.rm=T),
				sgnBD = apply(sg1 & sg2 & (sign(d1) == sign(d2)),2,sum,na.rm=T)
	)
}

#f = grepl('mr',alt.sp[alt.sp!=''])
s1 = 'mouse'
s2s = c('rat','rabbit','human','opossum')

pdf('figures/devAS/devAS.conservation.pdf',w=12,h=4.5)
par(mfrow=c(1,4),tck=-0.01,mgp=c(2,0.2,0),mar=c(5,3,1.5,0),oma=c(0,0,1,1))
f = alt.sp[alt.sp!=''] == 'hqmrboc'
for(s2 in s2s){
	z=compDevAS(f,orth.per.tissue.age.qv,orth.age.dpsi,s1,s2,0.1)
	barplot(z$sgnB/pmax(z$sgnT1,z$sgnT2),col=paste0(params$tissue.col[rownames(z)],'60'),border=NA,names.arg = rownames(z),ylab='exon fraction',main=paste(s1,s2,sep='-'),ylim=c(0,0.9),las=3)
	barplot(z$sgnBD/pmax(z$sgnT1,z$sgnT2),col=paste0(params$tissue.col[rownames(z)]),border=NA,add=T)
}
legend('toprigh',fill=c('#000000','#00000060'),legend=c('sgn12/max(sgn1,sgn2)','+same dir'))
mtext('Ancient exons',3,outer = T)


for(s2 in s2s){
	f = grepl(species[s1,'short'],alt.sp[alt.sp!='']) & grepl(species[s2,'short'],alt.sp[alt.sp!=''])
	z=compDevAS(f,orth.per.tissue.age.qv,orth.age.dpsi,s1,s2,0.1)
	barplot(z$sgnB/pmax(z$sgnT1,z$sgnT2),col=paste0(params$tissue.col[rownames(z)],'60'),border=NA,names.arg = rownames(z),ylab='exon fraction',main=paste(s1,s2,sep='-'),ylim=c(0,0.9),las=3)
	barplot(z$sgnBD/pmax(z$sgnT1,z$sgnT2),col=paste0(params$tissue.col[rownames(z)]),border=NA,add=T)
}
legend('toprigh',fill=c('#000000','#00000060'),legend=c('sgn12/max(sgn1,sgn2)','+same dir'))
mtext('Exons alt in both species',3,outer = T)
dev.off()

# devAS amplitude #####
age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)
age.dpsi$macaque = cbind(age.dpsi$macaque[,1:5],ovary=NaN,age.dpsi$macaque[,6,drop=FALSE])

pdf('figures/devAS/ad.dPSI.pdf',w=21,h=21)
par(mfrow=c(7,8),tck=-0.01,mgp=c(2,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,1,1))
brakes = 0:50/50
for(s in rownames(species)){
	f = anns[[s]]$sites=='ad'
	dp = lapply(colnames(per.tissue.age.qv[[s]]),function(t)abs(age.dpsi[[s]][f & per.tissue.age.qv[[s]][,t]<0.05,t]))
	names(dp) = colnames(per.tissue.age.qv[[s]])
	boxplot(dp,col=params$tissue.col[	names(dp) ],las=3,ylab='dPSI')
	for(t in 	names(dp))
		hist(dp[[t]],brakes,col=params$tissue.col[t],border=NA,main=t,xlab='dPSI')
}
dev.off()


# rphylopic #####
hs = name_search('human')
par(mfrow=c(3,3),tck=-0.01,mgp=c(2,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,1,1))
for(i in hs$canonicalName$uid){
	plot(1, 1, type="n", main="")
	add_phylopic_base(image_data(i,size=128)[[1]], 0.5, 0.5, 1)
}