source('code/r.functions/load.all.data.F.R')
source('~/skoltech/r.code/util.R')
source('code/r.functions/paper.figures.5.F.R')

library(openxlsx)
library(png)
library(reshape)
library(SAJR)
library(ape)
library(ontologyIndex)
library(extrafont)
library(doMC)
#https://askubuntu.com/questions/651441/how-to-install-arial-font-and-other-windows-fonts-in-ubuntu
# font_import()
# fonts()

# load data #####
species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
anns = readRDS('Rdata/anns.Rdata')
all.anns = readRDS('Rdata/all.anns.Rdata')
#border.stages = readRDS('Rdata/border.stages.Rdata')
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')
orth.seg.ad.tsm = readRDS('Rdata/orth.seg.ad.tsm.Rdata')
seg2ens = readRDS('Rdata/seg2ens.Rdata')
orth.per.tissue.age.qv = readRDS('Rdata/orth.per.tissue.age.qv.Rdata')
age.al.i = read.table('input/margarida.age.alignment.final.2017.05.17.improved.tab',sep='\t',header = T)[,c(rownames(species),'to.remove')]
gene.descrs = readRDS('Rdata/ens.gene.descr.Rdata')
orth.seg.ad.all.tsm = readRDS('Rdata/orth.seg.ad.all.tsm.Rdata')
age.segs = readRDS('Rdata/devAS.4patt.Rdata')

orth.ads.all.sp = readRDS('Rdata/orth.segs.all.types.all.phylo.groups.Rdata')
orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
orth.ens.genes = readRDS('Rdata/orth.ens.genes.Rdata')
ens.ge.marg.tsm = readRDS('Rdata/ens.ge.marg.tsm.Rdata')
ens.ge = readRDS('Rdata/ens.ge.Rdata')

# my.ge = readRDS('Rdata/my.ge.Rdata')
#ens.ge.cod = readRDS('Rdata/ens.ge.cod.Rdata')

params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)
lab.cex=1.5

plotPanelLetter = function(l,cex=1.2,adj=c(0,1.1),...){
	l = tolower(l)
	x=grconvertX(0,from='nfc',to='user')
	y=grconvertY(1,from='nfc',to='user')
	text(x=x,y=y,labels=l,adj=adj,font=2,cex=cex,xpd=NA)
}

# precalculate #####
phyl.tree=read.tree(text = '((((human:29.44,macaque:29.44):60.38,((mouse:20.89,rat:20.89):61.25,rabbit:82.14):7.68):68.77,opossum:158.6):153.4,chicken:312);')
phyl.tree$edge.length = rev(phyl.tree$edge.length)
phyl.tree$edge = phyl.tree$edge[nrow(phyl.tree$edge):1,]

# age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
# names(age.dpsi) = rownames(species)
# age.dpsi$macaque = cbind(age.dpsi$macaque[,1:5],ovary=NaN,age.dpsi$macaque[,6,drop=FALSE])
age.dpsi = readRDS('Rdata/age.diam.spline4.with.replicates.Rdata')

# alt.sp = readRDS('Rdata/paper.figures/alt.sp.Rdata')
# alt.sp = alt.sp[alt.sp!='']
alt.sp =getAltSp(orth.seg.ad.all.tsm,0.9 ,4)[rownames(orth.per.tissue.age.qv$human)]

# orth.age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(orth.seg.ad.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
# names(orth.age.dpsi) = rownames(species)
# orth.age.dpsi$macaque = cbind(orth.age.dpsi$macaque[,1:5],ovary=NA,orth.age.dpsi$macaque[,6,drop=F])
orth.age.dpsi = readRDS('Rdata/orth.ad.age.diam.spline4.with.replicates.Rdata')

# make "corresponding age rank"
car = unlist(lapply(colnames(age.al.i)[1:7],function(s){
	f = age.al.i[,s] != ''
	sapply(split((1:14)[f],paste(s,age.al.i[f,s])),mean)
}))
meta.tsm$corr.age.rank = car[paste(meta.tsm$species,meta.tsm$stage)]

# _MDS ####
#orth.cor = readRDS('Rdata/paper.figures/orth.cor5.Rdata')
# orth.all.psi = do.call(cbind,lapply(orth.seg.ad,function(x)x$ir))
# orth.all.psi = orth.all.psi[,rownames(meta)]
# 
# orth.cor = list()
# #orth.cor$all = cor(orth.all.psi,u='p')
# #orth.cor$sp1 = cor(orth.all.psi[nchar(alt.sp)==1,],u='p')
# orth.cor$sp7 = cor(orth.all.psi[nchar(alt.sp)==7,],u='p')
# #orth.cor$sp7tsm = cor(do.call(cbind,orth.seg.ad.tsm)[alt.cnt==7,],u='p')
# orth.seg.ad.ann = readRDS('Rdata/orth.seg.ad.all.ann.Rdata')
# f = rownames(orth.seg.ad.ann$human) %in% rownames(orth.all.psi)
# l = sapply(orth.seg.ad.ann,function(x)x$length[f])
# l = apply(l<=27,1,sum)
# table(l,nchar(alt.sp))
# orth.cor$sp7.macro = cor(orth.all.psi[nchar(alt.sp)==7 & l == 0,],u='p')
# orth.cor$sp7.micro = cor(orth.all.psi[nchar(alt.sp)==7 & l == 7,],u='p')
# saveRDS(orth.cor,'Rdata/paper.figures/orth.cor5.Rdata')
# orth.mds = lapply(orth.cor,function(x)cmdscale(1-x,k=2))
# saveRDS(orth.mds,'Rdata/paper.figures/orth.mds5.Rdata')



# look for examples #####
# orth.seg.ad.all = readRDS('Rdata/orth.seg.ad.all.Rdata')
# o = setNames(orth.seg.ad.all$human$seg$north,rownames(orth.seg.ad.all$human$seg))
# table(o)
# saveRDS(o,'Rdata/paper.figures/orth.seg.ad.all.north.Rdata')

hgmd = read.table('input/hgmd/tosend/2017_1_HGMD_ALL_Variants.csv',sep='\t',row.names=1,header=T,quote='',comment.char = '')
colnames(hgmd)
ens.orth = readRDS('Rdata/paper.figures/orth.seg.ad.all.north.Rdata')
table(ens.orth)
f = alt.sp != ''
f = ens.orth[f]==7 & alt.sp[f]=='hqmrboc' & apply(is.na(orth.age.dpsi$mouse) | is.na(orth.per.tissue.age.qv$mouse),1,sum)==0 & abs(orth.age.dpsi$mouse[,'brain']) > 0.7 & orth.per.tissue.age.qv$mouse[,'brain'] < 0.05
table(f)
sids = rownames(orth.age.dpsi$mouse)[f]
gd = gene.descrs$mouse[unique(unlist(seg2ens$mouse[sids])),]
gd[toupper(gd$gene.name) %in% toupper(hgmd$gene),]
hgmd[hgmd$gene=='DLG3',]

# sid = sids[sapply(seg2ens$mouse[sids],function(x)'ENSMUSG00000032479' %in% x)] #MAP4 - OK
# sid = sids[sapply(seg2ens$mouse[sids],function(x)'ENSMUSG00000042605' %in% x)] #Atxn2 - not so nice
# sid = sids[sapply(seg2ens$mouse[sids],function(x)'ENSMUSG00000034826' %in% x)] #Nup54 - looks nice
sid = sids[sapply(seg2ens$mouse[sids],function(x)'ENSMUSG00000000881' %in% x)] #Dlg3 - looks nice + Intellectual disability (HGMD)
# sid = sids[sapply(seg2ens$mouse[sids],function(x)'ENSMUSG00000022377' %in% x)] #Asap1 - looks nice + Schizophrenia? (HGMD)
# sid = sids[sapply(seg2ens$mouse[sids],function(x)'ENSMUSG00000039844' %in% x)] #Rapgef1 - looks nice + Intellectual disability? (HGMD)

dlg3.inx = c(12502,8396,12506)
dlg3.hsegs = anns$human[rownames(orth.age.dpsi$human)[dlg3.inx],]
lapply(1:nrow(dlg3.hsegs),function(i){
		hgmd[hgmd$chrom_VCF_hg19==dlg3.hsegs$chr_id[i] & hgmd$pos_VCF_hg19 >= (dlg3.hsegs$start[i]-100) & hgmd$pos_VCF_hg19 <= (dlg3.hsegs$stop[i]+100),]
	})
hgmd[hgmd$chrom_VCF_hg19=='X' & hgmd$pos_VCF_hg19 >= 69711304 & hgmd$pos_VCF_hg19 <= 69718995,]
sid = which(rownames(orth.age.dpsi$mouse) %in% sid)
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

# m = meta
# m$old.stage = m$stage
# m$stage[grep('wpc',m$stage)] = paste(gsub('wpc','',m$stage[grep('wpc',m$stage)]),'wpc')
# m$stage[!is.na(as.numeric(m$stage))] = paste0('e',m$stage[!is.na(as.numeric(m$stage))])
# m$stage[grep('teenager',m$stage)] = 'teen'
# m$stage[m$stage=='youngadult']  = '25-35 y'
# m$stage[m$stage=='youngmidage'] = '36-45 y'
# m$stage[m$stage=='oldermidage'] = '46-55 y'
# m$stage[m$stage=='senior']      = '56-65 y'
# 
# m$stage[m$species=='macaque' & m$stage=='p152'] = 'p6m'
# m$stage[m$species=='macaque' & m$stage=='p365'] = 'p1y'
# m$stage[m$species=='macaque' & m$stage=='p1095'] = 'p3y'
# m$stage[m$species=='macaque' & m$stage=='p3285'] = 'p9y'
# m$stage[m$species=='macaque' & m$stage=='p5475'] = 'p15y'
# m$stage[m$species=='macaque' & m$stage=='p8030'] = 'p20-26y'
# 
# f = !(m$species %in% c('human','macaque')) & substr(m$stage,1,1)!='e'
# m$stage[f] = paste0('p',m$days[f] - species[m$species[f],'gestation'])
# m$stage[substr(m$stage,1,1)=='p'] = first2Upper(m$stage[substr(m$stage,1,1)=='p'] )
# 
# meta$paper.stages = m$stage
# saveRDS(meta,'Rdata/main.set.meta.Rdata')
# 1 #####
# _prepare ####

# plot(phyl.tree)
# at = c(300,180,90,25)
# axis(1,312-at,at,las=2)
# title(xlab='Million years')

# B
# change to getSegTestDevAsStat
astypes = c(CE='ad',AA='aa',AD='dd',RI='da')
tested.stat = lapply(setNames(astypes,astypes),function(ss)sapply(names(per.tissue.age.qv),function(s)apply(!is.na(per.tissue.age.qv[[s]][anns[[s]]$sites==ss,]),2,sum)))

# C
sgn02.stat = lapply(setNames(astypes,astypes),function(ss) {
	sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][anns[[s]]$sites == ss,] < 0.05 & abs(age.dpsi[[s]][anns[[s]]$sites == ss,])>0.2,2,sum,na.rm=T))/tested.stat[[ss]]*100})

astypes.pchs=c(ad=19,aa=2,dd=6,da=13)

sgn.cnt=lapply(setNames(astypes,astypes),function(ss) {
	sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][anns[[s]]$sites == ss,] < 0.05 & abs(age.dpsi[[s]][anns[[s]]$sites == ss,])>0.2,2,sum,na.rm=T))})

table(sgn.cnt$ad > sgn.cnt$da)
table(sgn02.stat$ad > sgn02.stat$dd)
# _no of tissue-specific devAS ######
sgn = list()
for(s in rownames(species)){
	sgn[[s]] = abs(age.dpsi[[s]]) > 0.2 & per.tissue.age.qv[[s]] < 0.05
	sgn[[s]][is.na(sgn[[s]])] = FALSE
}

sort(round(sapply(rownames(species),function(s){
	t = table(apply(sgn[[s]][anns[[s]]$sites=='ad',-2],1,sum))
	t['1']/sum(t[-1])
	})*100,1))



# dPSI > 0.2, ad
# chicken.1  rabbit.1     rat.1   mouse.1   human.1 opossum.1 macaque.1 
# 63.9      64.5      65.1      65.2      69.6      71.8      83.9 
# dPSI > 0.5, ad
# rat.1 chicken.1   mouse.1  rabbit.1   human.1 opossum.1 macaque.1 
# 84.7      85.5      86.7      86.8      87.8      90.7      95.6 

#_plot #####
pdf('figures/paper.figures/6/1/1a.pdf',w=7.2,h=7.2/3,family='Arial')
layout(matrix(1:4,ncol=4,byrow = T),widths = c(1.8,3.1,1,3.1))
par(tck=-0.01,mgp=c(1,0.2,0),mar=c(1.5,2,1.5,0),oma=c(0,0,0,0))
#plot1A('A')
plot.new()
plotPanelLetter('A',lab.cex)
# order by species/tissue
plotAsEventCount(tested.stat,astypes.pchs,by.tissue = T,ylab='# of detected events',main='Detected AS',bty='n')
plotPanelLetter('B',lab.cex)
par(mar=c(0,0,1.5,0))
y=plotASTypes(astypes.pchs)
legend(0,y-2,pch=19,col=params$tissue.col,legend = names(params$tissue.col),bty='n',cex = 0.8,xjust=0,yjust=1)
rect(3,2,100,100,xpd=NA,col=NA,border='gray')
par(mar=c(1.5,2,1.5,0))
plotAsEventCount(sgn02.stat,astypes.pchs,by.tissue = T,ylab='% of devAS',main='DevAS',bty='n')
plotPanelLetter('C',lab.cex)
dev.off()

ad.stat = lapply(all.anns, function(h){
	h = split(h,h$gene_id)
	t = sapply(h,function(x)sum(x$sites=='ad'))
	c(mult.exn.genes=sum(t>0),internal.exons=sum(t),sd=sd(t[t>0]))
})

ad.stat = do.call(rbind,ad.stat)
ad.stat[,2]/ad.stat[,1]

pdf('figures/paper.figures/6/1/ann.stat.pdf',w=10,h=3,family='Arial')
par(mfcol=c(1,4),mar=c(6,3,1,0),las=0,oma=c(0,0,1.2,1),mgp=c(1,0.2,0),las=3)
barplot(ad.stat[,1],ylab='# of genes with internal exon[s]',main='Genes with internal exons',border=NA)
barplot(ad.stat[,2],ylab='# of internal exons',border=NA,main='Internal exons')
y = ad.stat[,2]/ad.stat[,1]
sd = ad.stat[,3]/sqrt(ad.stat[,1])
b=barplot(y,ylab='# of internal exons per gene',border=NA,main='Internal exons per gene',ylim=c(0,max(y+2*sd)))
segments(b,y-2*sd,b,y+2*sd)
z=sapply(species$short,function(s)table(grepl(s,alt.sp)))
ci=sapply(z[2,],function(x)binom.test(x,83888)$conf*83888)
b=barplot(z[2,],names.arg = rownames(species),border=NA,ylab='# of orthologous cassette exons',main='CE among 83888 orth. exons',ylim=range(0,ci))
segments(b,ci[1,],b,ci[2,])
dev.off()


pdf('figures/paper.figures/6/1/event.stat.pdf',w=12,h=6,family='Arial')
par(mfcol=c(3,7),mar=c(2,3,1,0),las=0,oma=c(0,0,1.2,1),mgp=c(1,0.2,0))
astypes.pchs=c(ad=19,aa=2,dd=6,da=13)
for(at in names(astypes.pchs)){
	for(t in unique(meta$tissue)){
		barplot(tested.stat[[at]][t,],names.arg = species$short,border=NA,col=params$tissue.col[t],ylab='# of detected events',main='Detected AS')
		barplot(sgn.cnt[[at]][t,]    ,names.arg = species$short,border=NA,col=params$tissue.col[t],ylab='# of devAS events',main='DevAS')
		barplot(sgn02.stat[[at]][t,] ,names.arg = species$short,border=NA,col=params$tissue.col[t],ylab='% of devAS events',main='DevAS')
	}
	mtext(c(ad='Cassette exons',aa='Alternative donors',dd='Alternative donors',da='Retained introns')[at],3,line = 0,outer = T)
}
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
#mouse.ge = readRDS('Rdata/ens.ge.marg.tsm.Rdata')[[sp.ce]]
mouse.ge = readRDS('Rdata/ens.ge.marg.Rdata')[[sp.ce]]
mouse.ge = mouse.ge + min(mouse.ge[mouse.ge!=0],na.rm=T)
mouse.ge = log2(mouse.ge)

# use.mean.embryo = FALSE
#sd = apply(psi.tsm[[sp.ce]],1,sd)
#psi.cor2embryo = caclCor2Embryo(psi.tsm[[sp.ce]][anns[[sp.ce]]$sites=='ad' & sd > 0.05,],meta.tsm,cor.m = 'sp',use.mean.embryo=use.mean.embryo)

# psi.cor2embryo = caclCor2Embryo(psi.tsm[[sp.ce]][anns[[sp.ce]]$sites=='ad',],meta.tsm,cor.m = 'p',use.mean.embryo=use.mean.embryo,use.bootstrap = 1000)
# ge.cor2embryo  = caclCor2Embryo(mouse.ge,meta.tsm,cor.m = 'p',use.mean.embryo=use.mean.embryo,use.bootstrap = 1000)
# saveRDS(psi.cor2embryo,'Rdata/paper.figures/psi.cor2embryo.tsm.boot1000.Rdata')
# saveRDS(ge.cor2embryo,'Rdata/paper.figures/ge.cor2embryo.tsm.boot1000.Rdata')

# registerDoMC(3)
# psi.for.ce = readRDS(paste0('Rdata/',sp.ce,'.as.u.filtered.Rdata'))
# psi.cor2embryo = caclCor2EmbryoAllStat(psi.for.ce$ir[psi.for.ce$seg$sites=='ad',],meta,cor.m = 'p',boot=1000)
# ge.cor2embryo = caclCor2EmbryoAllStat(mouse.ge,meta,cor.m = 'p',boot=1000)
# saveRDS(psi.cor2embryo,'Rdata/paper.figures/psi.cor2embryo.boot1000.Rdata')
# saveRDS(ge.cor2embryo,'Rdata/paper.figures/ge.cor2embryo.boot1000.Rdata')
psi.cor2embryo = readRDS('Rdata/paper.figures/psi.cor2embryo.boot1000.Rdata')
ge.cor2embryo = readRDS('Rdata/paper.figures/ge.cor2embryo.boot1000.Rdata')

# tau
human.tau = read.csv('input/gene.info.from.marg/Human.Indexes.All.csv')
rownames(human.tau) = human.tau[,1]
table(substr(human.tau$Human_ID,1,1),human.tau$Human_ID %in% rownames(ens.ge$human$gene))
human.tau$biotype = setNames(ens.ge$human$gene$biotype,rownames(ens.ge$human$gene))[human.tau$Human_ID]

psi.tsm.ad = lapply(names(psi.tsm),function(s)psi.tsm[[s]][anns[[s]]$sites=='ad',])
names(psi.tsm.ad) = names(psi.tsm)
ts = unique(meta$tissue)
# human.tissue.gene.dpsi=lapply(setNames(ts,ts), function(t){print(t);t(getsPSIbyEnsID(age.dpsi$human[,t],seg2ens$human,use.random = F))})
# saveRDS(human.tissue.gene.dpsi,'Rdata/paper.figures/human.tissue.gene.dpsi.Rdata')
human.tissue.gene.dpsi = readRDS('Rdata/paper.figures/human.tissue.gene.dpsi.Rdata')

# % of devAS in tau < 0.5
hb.devAS = rownames(anns$human)[anns$human$sites=='ad' & !is.na(per.tissue.age.qv$human[,'brain']) & per.tissue.age.qv$human[,'brain'] < 0.05 & age.dpsi$human[,'brain'] > 0.2]
hb.devAS = rownames(anns$human)[anns$human$sites=='ad' & !is.na(per.tissue.age.qv$human[,'brain']) &  age.dpsi$human[,'brain'] > 0.2]
table(is.na(hb.devAS))

hb.devAS.tau = sapply(seg2ens$human[hb.devAS],function(x)min(human.tau[x,'TestisTau'],na.rm=T))
sum(hb.devAS.tau<0.5,na.rm=T)/length(hb.devAS.tau) # 0.6929326 
sum(hb.devAS.tau<0.5,na.rm=T)/sum(!is.infinite(hb.devAS.tau)) # 0.7515419

# _plot #####
pdf('figures/paper.figures/6/2/2.v3.boot1000a.pdf',w=7.2/3*2,h=7.2,family='Arial')
#jpeg('figures/paper.figures/6/2/2.v2.jpg',units = 'in',w=7.2/3*2,h=7.2,quality = 100,res = 600,family='Arial')
layout(matrix(1:6,ncol=2,byrow = T),widths = c(10,10,10))
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
hc = !is.na(human.tau$biotype) & human.tau$biotype=='protein_coding'
sum(hc) #22519 protein coding genes 13308 and 12284 of them have AS in brain and liver respectively
t=plotAsInTauDistr(human.tau[hc,],human.tissue.gene.dpsi$brain,0.2,main='Brain',border=NA)
#fisher.test(table(t$as,t$TissueTau<0.5)[2:3,])
plotPNG("figures/paper.figures/5/icons/human.png",0.8,0.85,0.15)
legend(4,3000,fill=c('orange','lightgray','darkgray'),legend=c('devAS','AS','non-AS'),bty='n')
plotPanelLetter('C',lab.cex)
# plotAsInTauDistr(human.tau[hc,],human.tissue.gene.dpsi$heart,0.2,main='Heart')
# plotPanelLetter('E',lab.cex)
plotAsInTauDistr(human.tau[hc,],human.tissue.gene.dpsi$liver,0.2,main='Liver',border=NA)
plotPNG("figures/paper.figures/5/icons/human.png",0.8,0.85,0.15)
dev.off()

# 3 #####
# _prepare ####
# phastcons = read.table('/home/mazin/skoltech/projects/evo.devo/processed/ad.phastcons.gz',sep='\t')
# phastcons = setNames(strsplit(phastcons[,2],',',TRUE),phastcons[,1])
# phastcons = lapply(phastcons,as.numeric)
# age.ad.ph = list()
# inx = 151:200
# s='human'
# h.age.sids = list()
# h.aged.sids = list()
# for(t in unique(meta$tissue)){
# 	f = anns[[s]]$sites=='ad' & per.tissue.age.qv[[s]][,t]<0.05 & abs(age.dpsi[[s]][,1]) > 0.2
# 	f[is.na(f)] = FALSE
# 	h.age.sids[[t]] = rownames(anns[[s]])[f]
# 	h.aged.sids[[paste(t,'up')]] = rownames(anns[[s]])[f & age.segs[[s]][,t]=='u']
# 	h.aged.sids[[paste(t,'dw')]] = rownames(anns[[s]])[f & age.segs[[s]][,t]=='d']
# }
# all = unique(unlist(h.age.sids))
# h = readRDS("Rdata/human.as.u.all.Rdata")$seg;gc()
# h.aged.sids$cnst = h.age.sids$cnst = rownames(h)[h$sites=='ad' & h$type=='EXN' & h$gene_id %in% h[all,'gene_id']]
# h.aged.sids$`non-devAS` = h.age.sids$`non-devAS` = rownames(anns[[s]])[anns[[s]]$sites=='ad' & !(rownames(anns[[s]]) %in% all)]
# sapply(h.aged.sids,length)
# 
# age.ad.ph = lapply(h.age.sids,function(sids)sapply(phastcons[intersect(names(phastcons),sids)],function(x)mean(x[c(inx,length(x)-inx)])))
# age.ad.ph = age.ad.ph[c('cnst','non-devAS',unique(meta$tissue))]
# 
# aged.ad.ph = lapply(h.aged.sids,function(sids)sapply(phastcons[intersect(names(phastcons),sids)],function(x)mean(x[c(inx,length(x)-inx)])))
# aged.ad.ph = aged.ad.ph[c('cnst','non-devAS',setdiff(names(aged.ad.ph),c('cnst','non-devAS')))]
# 
# saveRDS(age.ad.ph,'Rdata/paper.figures/age.ad.ph.50nt.v2.Rdata')
# saveRDS(aged.ad.ph,'Rdata/paper.figures/age.directed.ad.ph.50nt.v2.Rdata')


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

# dev.as.cons.all = list()
# sps = c('rat','rabbit','human','opossum','chicken')
# f = alt.sp == 'hqmrboc'
# f = T
# for(s2 in sps)
# 	for(t in unique(meta$tissue)){
# 		dev.as.cons.all[[length(dev.as.cons.all)+1]] = getDevASCons('mouse',s2,t,orth.per.tissue.age.qv,orth.age.dpsi,f)
# 	}
# 
# dev.as.cons.all = do.call(rbind,dev.as.cons.all)

# phastcons
library(RColorBrewer)
orth.mds2 = readRDS('Rdata/paper.figures/orth.mds5.Rdata')
dlg3.mdata = readRDS('Rdata/paper.figures/dlg3.mdata.Rdata')
age.ad.ph = readRDS('Rdata/paper.figures/age.ad.ph.50nt.v2.Rdata')
# age.ad.ph = age.ad.ph.[1:2]
# for(t in unique(meta$tissue))
# 	age.ad.ph[[t]] = c(age.ad.ph.[[paste(t,'up')]],age.ad.ph.[[paste(t,'dw')]])

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
ndevAS.len3.freq = c(0.3935389,0.3770803,0.4101806)#my.binom.test(table(anns[[s]]$length[anns[[s]]$cod=='c' & !sgn & anns[[s]]$sites=='ad' &  anns[[s]]$gene_id %in% anns[[s]]$gene_id[sgn]]  %% 3 == 0)[c('TRUE','FALSE')])

ndevas = as.numeric(table(anns[[s]]$length[anns[[s]]$cod=='c' & !sgn & anns[[s]]$sites=='ad' &  anns[[s]]$gene_id %in% anns[[s]]$gene_id[sgn]]  %% 3 == 0)[c('TRUE','FALSE')])

apply(ad.3n,2,function(x){r=prop.test(c(ndevas[1],x[1]),c(sum(ndevas),x[2]));c(r$p.value,r$estimate)})

# exn.len3.freq = my.binom.test(table(all.anns[[s]][all.anns[[s]]$sites=='ad' & all.anns[[s]]$type=='EXN' & all.anns[[s]]$gene_id %in% anns[[s]]$gene_id[sgn],'length'] %% 3 == 0)[c('TRUE','FALSE')])
# ndevAS.len3.freq = my.binom.test(table(anns[[s]]$length[!sgn & anns[[s]]$sites=='ad']  %% 3 == 0)[c('TRUE','FALSE')])
ad.3n = t(cbind(cnst=exn.len3.freq,'non-devAS'=ndevAS.len3.freq,ad.3n[-1:-2,]))
ad.3n = ad.3n[,c(2,1,3)]

# IDR
#o = get_OBO('processed/exon.onthology/exont.obo')
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

ndev.cnt = table(alt %in% iupr)[c('TRUE','FALSE')]
idr02.cnt = apply(hdevas02,2,function(x)table(rownames(hdevas02)[f & x] %in% iupr)[c('TRUE','FALSE')])
apply(idr02.cnt,2,function(x){r=prop.test(c(ndev.cnt[1],x[1]),c(sum(ndev.cnt),sum(x)));c(r$p.value,r$estimate)})

# _plot #####
pdf('figures/paper.figures/6/3/3a.pdf',w=7.2,h=7.2/3*2.5,family='Arial')
#png('figures/paper.figures/6/3/3.png',units = 'in',w=7.2,h=7.2/3*2.5,res = 600,family='Arial')
#layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5),ncol=6,byrow = T),heights = c(3,2))
layout(matrix(c(1,1,1,1,2,2,2,2,3,4,4,7,5,5,6,6),ncol=8,byrow = T),heights = c(3,2))
par(tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,1),oma=c(0,0,0,1))
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
plotPanelLetter('C',lab.cex)
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
plot.new()
plotPanelLetter('B',lab.cex)
stat = t(sapply(age.ad.ph,function(x){m=mean(x);s=sd(x)/sqrt(length(x));c(m-2*s,m,m+2*s)}))
x = 1:nrow(stat)
cols=c('gray','black',params$tissue.col[rownames(stat)[-1:-2]])
par(mar=c(3.5,2.5,1.5,2))
plot(x,stat[,2],col=cols,main="PhastCons",ylab="mean intronic PhastCons (±50nt)",pch=19,bty='n',xaxt='n',xlab='',cex=2,ylim=range(stat),xpd=T)
abline(h=stat[1:2,2],lty=3,col=c('gray','black'))
#segments(x,stat[,1],x,stat[,3],col=cols)
arrows(x,stat[,1],x,stat[,3],col=cols,angle=90,code=3,length=0.03,xpd=T)
text(x,grconvertY(-0.02,'npc','user'),rownames(stat),adj=c(0,0),srt=-45,xpd=T)


plot(x,ad.3n[,2],col=cols,main="3N",ylab="fraction of 3N exons",pch=19,bty='n',xaxt='n',xlab='',cex=2,ylim=range(ad.3n),xpd=T)
abline(h=ad.3n[1:2,2],lty=3,col=c('gray','black'))
#segments(x,ad.3n[,1],x,ad.3n[,3],col=cols)
arrows(x,ad.3n[,1],x,ad.3n[,3],col=cols,angle=90,code=3,length=0.03,xpd=T)
text(x,grconvertY(-0.02,'npc','user'),rownames(ad.3n),adj=c(0,0),srt=-45,xpd=T)
plotPanelLetter('D',lab.cex)


plot(x,idr02[,2],col=cols,main="IDR",ylab="fraction of IDR",pch=19,bty='n',xaxt='n',xlab='',cex=2,ylim=range(idr02),xpd=T)
abline(h=idr02[1:2,2],lty=3,col=c('gray','black'))
#segments(x,idr02[,1],x,idr02[,3],col=cols)
arrows(x,idr02[,1],x,idr02[,3],col=cols,angle=90,code=3,length=0.03,xpd=T)
text(x,grconvertY(-0.02,'npc','user'),rownames(idr02),adj=c(0,0),srt=-45,xpd=T)
plotPanelLetter('E',lab.cex)

par(tck=-0.01,mgp=c(1.3,0.2,0),mar=c(0,0,0,0),oma=c(0,0,0,0))
lx=grconvertX(1,'lines','ndc')
plotExampleDLG3(fig=c(0.5,1-lx,2/5,1))
dev.off()

pdf('figures/paper.figures/5/3/micro-macro.ancient.mds.pdf',w=9,h=3,family='Arial')
par(mfrow=c(1,3),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
m = meta[rownames(orth.mds2$sp7),]
plot(orth.mds2$sp7,xlab='Dim 1',ylab='Dim 2',col=params$tissue.col[m$tissue],pch=params$species.pch[m$species],cex=m$cex,bty='n',main='All ancient (1441)')
plot(orth.mds2$sp7.macro,xlab='Dim 1',ylab='Dim 2',col=params$tissue.col[m$tissue],pch=params$species.pch[m$species],cex=m$cex,bty='n',main='Macro ancient (1197)')
plot(orth.mds2$sp7.micro,xlab='Dim 1',ylab='Dim 2',col=params$tissue.col[m$tissue],pch=params$species.pch[m$species],cex=m$cex,bty='n',main='Micro ancient (223)')
dev.off()


names(age.ad.ph)
sapply(age.ad.ph,function(x)wilcox.test(x,age.ad.ph$`non-devAS`,a='g')$p.value)

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
# as.in.ge.patterns.mouse = list()
# sp = 'mouse'
# mc = read.csv(paste0('processed/GE.from.marg/',firstToupper(sp),'Clusters.csv'),row.names = 1)
# colnames(mc) = tolower(colnames(mc))
# for(tis in unique(meta$tissue)){
# 	cat(tis)
# 	t = getsPSIbyEnsID(age.dpsi[[sp]][anns[[sp]]$sites=='ad',tis],seg2ens[[sp]])
# 	t = t(t)
# 	t = t[intersect(rownames(t),rownames(mc)[!is.na(mc[,paste0(tis,'pattern')])]),]
# 	as.in.ge.patterns.mouse[[tis]] = cbind(data.frame(t),ge.pattern=mc[rownames(t),paste0(tis,'pattern')])
# }
# saveRDS(as.in.ge.patterns.mouse,'Rdata/paper.figures/as.in.ge.patterns.mouse.v2.Rdata')

mouse.liver.cor = calcCor2TisOnDev(psi.tsm$mouse[anns$mouse$sites=='ad',],'liver',meta.tsm,cor.meth = 'pearson')

as.in.ge.patterns.mouse = readRDS('Rdata/paper.figures/as.in.ge.patterns.mouse.v2.Rdata') 
DPSI = 0.5
as.in.ge.patterns.mouse.cnt = lapply(as.in.ge.patterns.mouse,function(x){table(factor(x[,3]),pmax(abs(x[,1]),abs(x[,2]))>DPSI)[,c('TRUE','FALSE')]})
as.in.ge.patterns.mouse.stat = sapply(as.in.ge.patterns.mouse.cnt,function(x){cbind(my.binom.test(x['Increasing',]),my.binom.test(x['Decreasing',]))})

m = 'pearson'

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
ge.cnt = read.csv('input/number.of.de.gene.on.dev.from.Margarida_UPDATED.csv')
table(ge.cnt$down+ge.cnt$up==ge.cnt$total)
ge.cnt$total = ge.cnt$down+ge.cnt$up
# comp with old:
# o = read.csv('input/number.of.de.gene.on.dev.from.Margarida.csv')
# dim(o)
# apply(o != ge.cnt,2,sum)
# table(o$down+o$up==o$total)
# table(o$species,o$up!=ge.cnt$up | o$down!=ge.cnt$down)
# f = o$up!=ge.cnt$up | o$down!=ge.cnt$down
# plot(o$up[f],ge.cnt$up[f],log='xy')
# abline(a=0,b=1)
# plot(o$down[f],ge.cnt$down[f],log='xy')
# abline(a=0,b=1)
# 
####

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
		as = apply(abs(dpsi.age[anns[[sp]]$sites=='ad' & per.tissue.age.qv[[sp]][,t]<0.05 & abs(age.dpsi[[sp]][,t]) > 0.2,f])>  dpsi,2,function(x){x = x[!is.na(x)];if(length(x)==0){NA}else{sum(x)}})
		as = as[!is.na(as) & names(as) %in% rownames(ge.cnt)]
		peak.changes[[sp]][[t]] = list(ge=ge.cnt[names(as),'total'],as=as)
	}
}


# _plot #####
pdf('figures/paper.figures/6/4/4a.pdf',w=6,h=6.2,family='Arial')
#jpeg('figures/paper.figures/6/4/4.jpg',units = 'in',w=6,h=6.2,quality = 100,res = 600,family='Arial')
layout(matrix(c(1,2,3,3,4,4),ncol=2,byrow = T),heights = c(3,1.5,1.7))
par(tck=-0.01,mgp=c(1.2,0.3,0),mar=c(4.5,2.5,1.5,0),oma=c(0,0,0,1),las=3)
plotTissueAgeProile(mouse.liver.cor,meta.tsm,age.axis = 'rank',bty='n',xlab='',ylab='Correlation to liver',ylim=c(0.8,1),pch=19,cex=1)
plotPNG("figures/paper.figures/5/icons/mouse.png",0.63,0.9,0.15)
plotPNG("figures/paper.figures/5/icons/liver.png",0.8,0.9,0.15)
title(xlab='Developmental stage',mgp=c(2,0.3,0))
plotPanelLetter('A',lab.cex)

col = c(paste0(rep(params$tissue.col,each=2),c('30','FF')))
b = barplot(as.in.ge.patterns.mouse.stat[c(4,1),],col=col,border=NA,ylim=c(0,max(as.in.ge.patterns.mouse.stat)),beside = T,ylab='Proportion of genes with |dPSI| > 0.5',xaxt='n')
#segments(b,as.in.ge.patterns.mouse.stat[c(5,2),],b,as.in.ge.patterns.mouse.stat[c(6,3),])
arrows(b,as.in.ge.patterns.mouse.stat[c(5,2),],b,as.in.ge.patterns.mouse.stat[c(6,3),],angle=90,code=3,length=0.02,xpd=T)
text(apply(b,2,mean),rep(0,7),colnames(as.in.ge.patterns.mouse.stat),xpd=NA,adj=c(0,1),srt=-45)

plot4B.legend()
plotPanelLetter('B',lab.cex)

par(mar=c(3,1.2,0,0))
plot4C.DivergenceOnAge('mouse','human',as.cor.on.dev$human,'c',yrange=c(1.7,3.2)/6.2)
#peak.changes
par(mar=c(3,1.2,1,0))
plot4D.PeakChange('rabbit',peak.changes$rabbit,'d',yrange=c(0,1.5)/6.2)
dev.off()

# 5 ######
# _prepare #####
# _four patterns #####
# age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
# names(age.dpsi) = rownames(species)
# 
# all02=sapply(rownames(species)[-2],function(s){
# 	x=getDevASpattern(psi.tsm[[s]],abs(age.dpsi[[s]])>0.2 & per.tissue.age.qv[[s]]<0.05,age.dpsi[[s]],meta.tsm,anns[[s]]$sites=='ad',at.least.two.sgn=F)
# 	c(u0=sum(x$dir=='u' & x$mean.psi<0.5,na.rm=T),
# 		u1=sum(x$dir=='u' & x$mean.psi>0.5,na.rm=T),
# 		d0=sum(x$dir=='d' & x$mean.psi<0.5,na.rm=T),
# 		d1=sum(x$dir=='d' & x$mean.psi>0.5,na.rm=T),
# 		total=sum(x$dir != 'n'))
# })
# 
# sweep(all02,2,all02[5,],'/')*100

# divide exons into: devAS 1 tissue, more than 1 and same dir, different dir
# dir.class.all=t(sapply(rownames(species)[-2],function(s){
# 	x=getDevASpattern(psi.tsm[[s]],abs(age.dpsi[[s]])>0.2 & per.tissue.age.qv[[s]]<0.05,age.dpsi[[s]],meta.tsm,anns[[s]]$sites=='ad',at.least.two.sgn=F)
# 	c(one=sum(x$ntissue==1),same=sum(x$ntissue>1 & x$dir %in% c('u','d')),diff=sum(x$dir=='b'))
# }))
# 
# dir.class.all.wo.cbc=t(sapply(rownames(species)[-2],function(s){
# 	x=getDevASpattern(psi.tsm[[s]],abs(age.dpsi[[s]])>0.2 & per.tissue.age.qv[[s]]<0.05,age.dpsi[[s]],meta.tsm[meta.tsm$tissue!='cerebellum',],anns[[s]]$sites=='ad',at.least.two.sgn=F)
# 	c(one=sum(x$ntissue==1),same=sum(x$ntissue>1 & x$dir %in% c('u','d')),diff=sum(x$dir=='b'))
# }))
# # 
# pdf('figures/paper.figures/5/5/devAS.dir.stat.pdf',w=6,h=6,family='Arial')
# par(mfrow=c(2,2),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(5,2.5,2.5,0),oma=c(0,0,0,1),las=3)
# barplot(t(dir.class.all),ylab='# of exons',args.legend = list(bty='n'),legend.text = c('one tissue','same direction','different direction'),ylim=c(0,10000))
# barplot(t(dir.class.all.wo.cbc),main='Without cerebellum',ylab='# of exons',ylim=c(0,10000))
# barplot(t(sweep(dir.class.all,1,apply(dir.class.all,1,sum),'/')*100),ylab='% of exons')
# barplot(t(sweep(dir.class.all.wo.cbc,1,apply(dir.class.all.wo.cbc,1,sum),'/')*100),main='Without cerebellum',ylab='% of exons')
# dev.off()

# age.segs = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,psi.thr = 0.2,border.stages,s)[anns[[s]]$sites=='ad',])
# names(age.segs) = rownames(species)
# for(s in names(age.segs)) age.segs[[s]][is.na(per.tissue.age.qv[[s]][rownames(age.segs[[s]]),colnames(age.segs[[s]])])] = '-'
# for(s in names(age.segs)) age.segs[[s]][age.segs[[s]] != '-' & per.tissue.age.qv[[s]][rownames(age.segs[[s]]),colnames(age.segs[[s]])]>0.05] = 'n'

# B phastcons
age.ad.ph = readRDS('Rdata/paper.figures/age.directed.ad.ph.50nt.v2.Rdata')

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

# gnomad
# for h.aged.sids see fig 3
# h = readRDS('Rdata/all.anns.Rdata')$human
# gc()
# gnomad = read.table('processed/gnomad201/human.ad.snp.tab.gz',sep='\t')
# colnames(gnomad) = c('chr_id','pos','id','ref','alt','qual','filter','seg_id','alt_cnt','freq','tot_cnt')
# hann = h[unique(unlist(h.aged.sids)),]
# gnomad = gnomad[gnomad$alt_cnt > 1 & gnomad$seg_id %in% rownames(hann),]
# s2g = hann[gnomad$seg_id,]
# gnomad$dist2start = gnomad$pos - s2g$start
# gnomad$dist2stop  = s2g$stop - gnomad$pos
# gnomad[s2g$strand==-1,c('dist2start','dist2stop')] = gnomad[s2g$strand==-1,c('dist2stop','dist2start')]
# gnomad$strand = s2g$strand
# 
# age.ad.gnomad = lapply(h.aged.sids,function(sids)log10(gnomad$freq[gnomad$seg_id %in% sids & ((gnomad$dist2start<0 & gnomad$dist2start>= -50) | (gnomad$dist2stop<0 & gnomad$dist2stop>= -50))]))
# age.ad.gnomad = age.ad.gnomad[c('cnst','non-devAS',setdiff(names(age.ad.gnomad),c('cnst','non-devAS')))]
# saveRDS(age.ad.gnomad,'Rdata/paper.figures/age.ad.gnomad.50nt.v2.Rdata')

age.ad.gnomad = readRDS('Rdata/paper.figures/age.ad.gnomad.50nt.v2.Rdata')

# hexamers
hex.dws.age02sgn = readRDS('Rdata/hex.dws.age02sgn.Rdata')
hex.ups.age02sgn = readRDS('Rdata/hex.ups.age02sgn.Rdata')

fa = readRDS('Rdata/ad.alt.fa.Rdata')

# proportion of known
hex2mot = read.table('output/hex2mot2sf.tab.gz')
# hex.tis.no = apply(hex.ups$up$ih.qv<0.05 | hex.dws$up$ih.qv<0.05 | hex.ups$dw$ih.qv<0.05 | hex.dws$dw$ih.qv<0.05,1,sum)
# known.hex.stat = table(pmin(hex.tis.no,5),known=names(hex.tis.no) %in% hex2mot$V1[hex2mot$V2!=''])
# rownames(known.hex.stat)[6] = '>5'
# k = rownames(hex.dws.age02sgn$up$ih.qv) %in% rownames(hex2mot)[hex2mot$V2!='']

hex.stat = getHexStat(hex.ups.age02sgn,hex.dws.age02sgn,0.05)

# _plot #####
pdf('figures/paper.figures/6/5/5a.v2.pdf',w=7.2,h=7.2,family='Arial')
#png('figures/paper.figures/6/5/5.png',units = 'in',w=7.2,h=7.2,res = 600,family='Arial')
l = rbind(c(1,2,5,6),
			c(3,4,5,6),
			c(7,7,8,9),
			c(10,10,11,12))
layout(l,widths = c(1,1,2,2),heights = c(1,1,2,2))
par(tck=-0.01,mgp=c(1.2,0.3,0),mar=c(0.5,1.3,2.5,0),oma=c(0,0,0,1))
#plot5A('A')

plotDevAS4Pattern('u','Up')
mtext('PSI',2,cex=0.8)
plotPNG("figures/paper.figures/5/icons/mouse.png",0.13,0.9,0.2)
plotPNG("figures/paper.figures/5/icons/brain.png",0.33,0.9,0.2)
plotPanelLetter('A',lab.cex)
plotDevAS4Pattern('d','Down')
par(mar=c(1.5,1.3,1.5,0))
plotDevAS4Pattern('ud','Up-down')
mtext('PSI',2,cex=0.8)
mtext('Development',1,0.2,cex=0.8)
plotDevAS4Pattern('du','Down-up')
mtext('Development',1,0.2,cex=0.8)

par(mar=c(4.5,2.5,1.5,0.3))

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

par(mar=c(4.5,2.5,1.5,2.3))
plotASSegStat(at,stat,up2dw.len.pv,cols,xax,lty=c(2,2,rep(1:2,times=7)),main='3N',ylab='fraction of 3N exons',pch=pch,bty='n')#,xlim=c(-4,max(at)))
plotPanelLetter('B',lab.cex)

par(mar=c(4.5,4.5,1.5,0.3))
stat = t(sapply(age.ad.ph,function(x){m=mean(x);s=sd(x)/sqrt(length(x));c(m-2*s,m,m+2*s)}))
pv = sapply(1:7,function(i)wilcox.test(age.ad.ph[[2*i+1]],age.ad.ph[[2*i+2]])$p.value)
plotASSegStat(at,stat,pv,cols,xax,lty=c(2,2,rep(1:2,times=7)),ylab="mean intronic PhastCons (±50nt)",main="PhastCons",pch=pch,bty='n')

legend(-14,0.60,pch=c(19,1),legend=c('Up','Down'),bt='n',xpd=NA)
#text(18,-2.95,c('p-value\n* <0.05\n** <0.01\n*** <0.001'),xpd=NA,adj=c(0,1),cex=1)
lw=grconvertY(0:1,'line','user')
text(-14,seq(0.50,by=(lw[1]-lw[2]),length.out=4),c(expression(paste(italic('P'),'-value')),'* <0.05','** <0.01','*** <0.001'),xpd=NA,adj=c(0,1),cex=1)


plotPanelLetter('C',lab.cex)#,adj = c(-2.5,1.1))

par(mar=c(4.5,2.5,1.5,2))
stat = t(sapply(age.ad.gnomad,function(x){m=mean(x);s=sd(x)/sqrt(length(x));c(m-2*s,m,m+2*s)}))
pv = sapply(1:7,function(i)wilcox.test(age.ad.gnomad[[2*i+1]],age.ad.gnomad[[2*i+2]])$p.value)
plotASSegStat(at,stat,pv,cols,xax,lty=c(2,2,rep(1:2,times=7)),main="Human intron constraints" ,ylab="log10(MAF)",pch=pch,bty='n')
pv = sapply(2:15,function(i)wilcox.test(age.ad.gnomad$cnst,age.ad.gnomad[[i]])$p.value)
#legend(17,-2.7,pch=c(19,1),legend=c('Up','Down'),bt='n',xpd=T)
#text(18,-2.95,c('p-value\n* <0.05\n** <0.01\n*** <0.001'),xpd=NA,adj=c(0,1),cex=1)
#lw=grconvertY(0:1,'line','user')
#text(18,seq(-2.95,by=(lw[1]-lw[2]),length.out=4),c(expression(paste(italic('P'),'-value')),'* <0.05','** <0.01','*** <0.001'),xpd=NA,adj=c(0,1),cex=1)

plotPanelLetter('D',lab.cex)

par(mar=c(4.5,2.5,1.5,0.3))
plotHexStat(hex.stat$hex.stat,hex.stat$hex.stat.known,c('E','F'))

plotMirroredMotFreq(fa,age.segs,'actaac','brain','heart',main = 'ACTAAC in brain devAS')
plotPanelLetter('G',lab.cex)
plotMirroredMotFreq(fa,age.segs,'actaac','heart','brain',main = 'ACTAAC in heart devAS',plot.leg = F)

par(mar=c(0,0,0,0))
plot5G3()
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

# 6 ######
# _prepare #####
# micrexons in ancient
sapply(rownames(species),function(s){
	f = anns[[s]]$sites=='ad' & anns[[s]]$length<28
	mean(rownames(anns[[s]])[f] %in% rownames(orth.seg.ad[[s]]$seg)[alt.sp=='hqmrboc'])
	})

sapply(rownames(species),function(s){
	f = anns[[s]]$sites=='ad' & anns[[s]]$length>27
	mean(rownames(anns[[s]])[f] %in% rownames(orth.seg.ad[[s]]$seg)[alt.sp=='hqmrboc'])
})

# condition on orth
sapply(rownames(species),function(s){
	a = orth.seg.ad[[s]]$seg
	mean(alt.sp[grepl(species[s,'short'],alt.sp) & a$sites=='ad' & a$length<28]=='hqmrboc')
})


sapply(rownames(species),function(s){
	a = orth.seg.ad[[s]]$seg
	mean(alt.sp[grepl(species[s,'short'],alt.sp) & a$sites=='ad' & a$length>27]=='hqmrboc')
})

# for figure
adj.segments = readRDS('Rdata/adj.segments.sids-n-psi.Rdata')
for(s in names(adj.segments)){
	adj.segments[[s]] = adj.segments[[s]][rownames(anns[[s]]),]
	adj.segments[[s]]$f = (is.na(adj.segments[[s]]$up.ajd.med.psi.rate) | adj.segments[[s]]$up.ajd.med.psi.rate>2) & 
		                    (is.na(adj.segments[[s]]$dw.ajd.med.psi.rate) | adj.segments[[s]]$dw.ajd.med.psi.rate>2)
}
sapply(names(adj.segments),function(s)table(adj.segments[[s]]$f,anns[[s]]$sites))


ff = function(a,s,l=Inf,st='ad',check.fake = FALSE){a[[s]]$sites %in% st & a[[s]]$length<=l & (!check.fake | adj.segments[[s]]$f)}
dPSI=0.2

me.cnt=sapply(names(anns),function(s)sum(ff(anns,s,l=27)))
me.cnt=do.call(rbind,lapply(1:7,function(i)me.cnt))
me.cnt.tis = sapply(names(anns),function(s)apply(!is.na(per.tissue.age.qv[[s]][ff(anns,s,l=27),]) & !is.na(age.dpsi[[s]][ff(anns,s,l=27),]),2,sum,na.rm=T))

sgn02u = sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][ff(anns,s),] < 0.05 & age.segs[[s]][ff(anns,s),] == 'u',2,sum,na.rm=T))
sgn02d = sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][ff(anns,s),] < 0.05 & age.segs[[s]][ff(anns,s),] == 'd',2,sum,na.rm=T))

sgn02u.me = sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][ff(anns,s,l=27),] < 0.05 & age.segs[[s]][ff(anns,s,l=27),] == 'u',2,sum,na.rm=T))
sgn02d.me = sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][ff(anns,s,l=27),] < 0.05 & age.segs[[s]][ff(anns,s,l=27),] == 'd',2,sum,na.rm=T))

# B-C v2
s = 'human'
mi.vs.ma.devAS = lapply(colnames(per.tissue.age.qv[[s]]),function(t){
	f = anns[[s]]$sites=='ad'
	sgn = per.tissue.age.qv[[s]][f,t] < 0.05 & abs(age.dpsi[[s]])[f,t] >  dPSI
	micro = anns[[s]]$length[f] <= 27
	table(devAS=sgn,micro=micro)
	})

in.vs.ex.micro = lapply(colnames(per.tissue.age.qv[[s]]),function(t){
	f = anns[[s]]$sites=='ad' & per.tissue.age.qv[[s]][,t] < 0.05 & abs(age.dpsi[[s]])[,t] >  dPSI
	micro = anns[[s]]$length[f] <= 27
	table(incl=age.segs[[s]][f,t] == 'u',micro=micro)
})
names(mi.vs.ma.devAS) = names(in.vs.ex.micro) = colnames(per.tissue.age.qv[[s]])

sapply(mi.vs.ma.devAS,function(x){ft=fisher.test(x);c(pv=ft$p.value,or=ft$estimate)})
sapply(in.vs.ex.micro,function(x){ft=fisher.test(x);c(pv=ft$p.value,or=ft$estimate)})

# D
micro.timing = lapply(rownames(species)[-c(2,7)], function(s){
	f = anns[[s]]$sites=='ad' & age.segs[[s]][,'brain'] == 'u'
	f = !is.na(f) & f
	t = psi.tsm[[s]][f,paste(s,'brain',c(border.stages[[s]]['brain',1],age.al.i[10,s],border.stages[[s]]['brain',2]))]
	t = table(micro=anns[[s]]$length[f]<=27,before.birth = factor((t[,2]-t[,1])>(t[,3]-t[,2]),levels=c(F,T)))
})
names(micro.timing) = rownames(species)[-c(2,7)]

micro.timing.b=sapply(micro.timing,function(x)my.binom.test(x[2,2:1]))*100
macro.timing.b=sapply(micro.timing,function(x)my.binom.test(x[1,2:1]))*100

sapply(micro.timing,function(x){ft=fisher.test(x);c(pv=ft$p.value,or=ft$estimate)})

phastcons = read.table('/home/mazin/skoltech/projects/evo.devo/processed/ad.phastcons.gz',sep='\t')
phastcons = setNames(strsplit(phastcons[,2],',',TRUE),phastcons[,1])
phastcons = lapply(phastcons,as.numeric)

s = 'human'
t = 'brain'
phast.sids = list('micro devAS'=rownames(anns[[s]])[ff(anns,s,27) & age.segs[[s]][,t] == 'u'],
									'macro devAS'=rownames(anns[[s]])[anns[[s]]$sites=='ad' & anns[[s]]$length>27 & age.segs[[s]][,t] == 'u'],
									'macro non-devAS'=rownames(anns[[s]])[anns[[s]]$sites=='ad' & anns[[s]]$length>27 & age.segs[[s]][,t] == 'n'])
phast.sids = lapply(phast.sids,function(x)x[!is.na(x)])
sapply(phast.sids,length)

# 3N
micro.3n = array(NA,dim=c(7,3,3),dimnames = list(rownames(species),names(phast.sids),c('mean','lCI','uCI')))
for(s in rownames(species)){
	l = list('micro devAS'=anns[[s]]$length[ff(anns,s,27) & age.segs[[s]][,t] == 'u'],
						'macro devAS'=anns[[s]]$length[anns[[s]]$sites=='ad' & anns[[s]]$length>27 & age.segs[[s]][,t] == 'u'],
						'macro non-devAS'=anns[[s]]$length[anns[[s]]$sites=='ad' & anns[[s]]$length>27 & age.segs[[s]][,t] == 'n'])
	l = lapply(l,function(x)x[!is.na(x)])
	micro.3n[s,,] = t(sapply(l,function(x)my.binom.test(sum(x%%3==0),sum(x%%3!=0))))
}

l = sapply(orth.seg.ad,function(x)x$seg$length)
dim(l)
length(alt.sp)

apply(l,2,function(x)sapply(split(x<28,nchar(alt.sp)==7),mean))

# _plot #####
points = FALSE
pdf('figures/paper.figures/6/6/6a.pdf',w=7.2,h=7.2/3*2,family='Arial')
#png('figures/paper.figures/6/6/6.png',units = 'in',w=7.2,h=7.2/3*2,res = 600,family='Arial')
par(mfrow=c(2,3),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(3,2.5,1.5,1),oma=c(0,0,0,1))
cols = paste0(rep(params$tissue.col,each=7),'44')
b=barplot(t(me.cnt.tis),beside = T,col=cols,names.arg=,xaxt='n',border=NA,ylab='# of microexons',main='devAS in microexons',ylim=c(0,770))
cols = rep(params$tissue.col,each=7)
barplot(t(sgn02u.me+sgn02d.me),beside = T,col=cols,ylab='',main='',xaxt='n',border=NA,add=T)
#text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.5)
legend(33,810,fill=c('#00000077','#000000'),legend = c('detected AS','devAS'),border = NA,bty = 'n',xpd=T)

x = b
segments(x,grconvertY(-0.005,'npc','user'),x,grconvertY(-0.02,'npc','user'),xpd=T)
x = x[,1]
xx = x*2 - mean(x)*1.5
# adjust q and r
xx[2] = xx[2] - (max(xx)-min(xx))/60
xx[4] = xx[4] + (max(xx)-min(xx))/60
segments(x[1:7],grconvertY(-0.02,'npc','user'),xx,grconvertY(-0.03,'npc','user'),xpd=T)
text(xx,grconvertY(-0.03,'npc','user'),species[colnames(me.cnt.tis),'short'],adj = c(0.5,1),xpd=T,cex=1)

plotPanelLetter('A',lab.cex)

# sps = c('human','mouse','opossum')
# if(points){
# 	b = seq(1.5,by=1,length.out=7*length(sps)) + rep(0:6,each=3)
# 	plot(b,as.numeric(t(sgn02u.me/sgn02u)[sps,]*100),col=rep(params$tissue.col,each=length(sps)),xlab='',ylab='% of microexons',main='Inclusion microexons',bty='n',xaxt='n',ylim=range(0,22),pch=19)
# 	text(sapply(split(b,rep(1:7,each=3)),'mean'),rep(grconvertY(-0.05,'npc','user'),7),substr(rownames(sgn02u.me),1,1),xpd=TRUE)
# }else{
# 	b=barplot(t(sgn02u.me/sgn02u)[sps,]*100,beside = T,col=rep(params$tissue.col,each=length(sps)),ylab='% of microexons',main='Inclusion microexons',names.arg=substr(rownames(sgn02u.me),1,1),border=NA,ylim=range(0,22))
# }
# text(b,t(sgn02u.me/sgn02u)[sps,]*100+0.2,t(sgn02u.me)[sps,],srt=90,cex=0.8,adj=c(0,0.5))
# text(b,0,species[sps,'short'],adj = c(0.5,1),xpd=T,cex=0.8)

at = (0:19)[-seq(3,19,by = 3)]
cols = c(rep(params$tissue.col,each=2))
pch=c(rep(c(19,1),times=7))
xax = setNames(seq(0.5,by=3,length.out = 7),names(mi.vs.ma.devAS))

stat = matrix(NA,ncol=3,nrow=14)
stat[0:6*2+1,] = t(sapply(mi.vs.ma.devAS,function(x)my.binom.test(x[c('TRUE','FALSE'),'TRUE'])))[,c(2,1,3)]
stat[0:6*2+2,] = t(sapply(mi.vs.ma.devAS,function(x)my.binom.test(x[c('TRUE','FALSE'),'FALSE'])))[,c(2,1,3)]

pv = sapply(mi.vs.ma.devAS,function(x){ft=fisher.test(x)$p.value})
plotASSegStat(at,stat,pv,cols,xax,skip.for.pv=NULL,lty=c(rep(1:2,times=7)),main='',ylab='Proportion of devAS',pch=pch,bty='n',first2horiz=FALSE)
plotPNG("figures/paper.figures/5/icons/human.png",0.93,0.9,0.1)
legend(grconvertX(0.5,'npc','user'),grconvertY(0.95,'npc','user'),pch=c(19,1),col='black',legend=c('microexons','macroexons'),xjust = 0.5,yjust=0,bty='n',xpd=TRUE,ncol=2)
lw=grconvertY(0:1,'line','user')
text(16.5,seq(0.34,by=0.9*(lw[1]-lw[2]),length.out=4),c(expression(paste(italic('P'),'-value')),'* <0.05','** <0.01','*** <0.001'),xpd=NA,adj=c(0,1),cex=0.9)

plotPanelLetter('B',lab.cex)

# if(points){
# 	b = seq(1.5,by=1,length.out=7*length(sps)) + rep(0:6,each=3)
# 	plot(b,as.numeric(t(sgn02d.me/sgn02d)[sps,]*100),col=rep(params$tissue.col,each=length(sps)),xlab='',ylab='% of microexons',main='Inclusion microexons',bty='n',xaxt='n',ylim=range(0,22),pch=19)
# 	text(sapply(split(b,rep(1:7,each=3)),'mean'),rep(grconvertY(-0.05,'npc','user'),7),substr(rownames(sgn02u.me),1,1),xpd=TRUE)
# }else{
# 	b=barplot(t(sgn02d.me/sgn02d)[sps,]*100,beside = T,col=rep(params$tissue.col,each=length(sps)),ylab='% of microexons',main='Exclusion microexons',names.arg=substr(rownames(sgn02u.me),1,1),ylim=range(0,22),border=NA)
# }
# text(b,t(sgn02d.me/sgn02d)[sps,]*100+0.2,t(sgn02d.me)[sps,],srt=90,cex=0.8,adj=c(0,0.5))
# text(b,0,species[sps,'short'],adj = c(0.5,1),xpd=T,cex=0.8)

stat = matrix(NA,ncol=3,nrow=14)
stat[0:6*2+1,] = t(sapply(in.vs.ex.micro,function(x)my.binom.test(x[c('TRUE','FALSE'),'TRUE'])))[,c(2,1,3)]
stat[0:6*2+2,] = t(sapply(in.vs.ex.micro,function(x)my.binom.test(x[c('TRUE','FALSE'),'FALSE'])))[,c(2,1,3)]

pv = sapply(in.vs.ex.micro,function(x){ft=fisher.test(x)$p.value})
plotASSegStat(at,stat,pv,cols,xax,skip.for.pv=NULL,lty=c(rep(1:2,times=7)),main='',ylab='Proportion of up',pch=pch,bty='n',first2horiz=FALSE)
plotPNG("figures/paper.figures/5/icons/human.png",0.73,0.9,0.1)
legend(grconvertX(0.5,'npc','user'),grconvertY(0.95,'npc','user'),pch=c(19,1),col='black',legend=c('microexons','macroexons'),xjust = 0.5,yjust=0,bty='n',xpd=TRUE,ncol=2)
plotPanelLetter('C',lab.cex)

par(mar=c(3,2.5,1.5,1))
x = seq(from=1,by=4,length.out=ncol(micro.timing.b))
plot(x,micro.timing.b[1,],pch=19,ylim=c(0,100),xlim=range(1,x+1),bty='n',xaxt='n',xlab='',ylab='% of exons mostly changed before birth')
points(x+1,macro.timing.b[1,],pch=1,cex=2)
arrows(x,micro.timing.b[2,],x,micro.timing.b[3,],angle=90,code=3,length=0.03,xpd=T)
arrows(x+1,macro.timing.b[2,],x+1,macro.timing.b[3,],angle=90,code=3,length=0.03,xpd=T)
axis(1,x+0.5,NA)
text(x+0.5,rep(-5,ncol(micro.timing.b)),colnames(micro.timing.b),srt=-45,adj=c(0,1),xpd=NA)
plotPNG("figures/paper.figures/5/icons/brain.png",0.25,0.32,0.3)
legend(10,40,bty='n',pch=c(19,1),pt.cex=c(1,2),legend=c('Microexons','Macroexons'))
plotPanelLetter('D',lab.cex)
micro.col=c('orange','olivedrab3','black')
plotPhast(phastcons,phast.sids,main='',ylim=c(0.15,1.20),bty='n',xlab='',col=micro.col)
plotPNG("figures/paper.figures/5/icons/brain.png",0.1,0.45,0.15)
plotPanelLetter('E',lab.cex)

x = seq(1,by=3,length.out=dim(micro.3n)[1])
x = c(x,x+1)
plot(x,micro.3n[,1:2,1],pch=19,col=rep(micro.col[1:2],each= dim(micro.3n)[1]),ylim=range(micro.3n[,1:2,]),ylab='fraction of 3N exons',bty='n',xaxt='n',xlab='',cex=2)
arrows(x,micro.3n[,1:2,2],x,micro.3n[,1:2,3],col=rep(micro.col[1:2],each= dim(micro.3n)[1]),angle=90,code=3,length=0.03,xpd=T)
text(x[1:dim(micro.3n)[1]]+0.5,grconvertY(0.02,'npc','user'),dimnames(micro.3n)[[1]],srt=-45,xpd=NA,adj=c(0,1))
plotPNG("figures/paper.figures/5/icons/brain.png",0.1,0.45,0.15)
plotPanelLetter('F',lab.cex)
dev.off()

# _brain<->micro Fisher ######
ct = lapply(rownames(species),function(s){
	sgn = apply(per.tissue.age.qv[[s]]<0.05 & abs(age.dpsi[[s]]) > 0.2,1,sum,na.rm=T)>0 & anns[[s]]$sites=='ad'
	b = per.tissue.age.qv[[s]][,'brain']<0.05 & abs(age.dpsi[[s]][,'brain']) > 0.2
	b[is.na(b)] = FALSE
	m = anns[[s]]$length<=27
	table(brain=b[sgn],micro=m[sgn])
	})
sapply(ct,function(x){ft=fisher.test(x);c(pv=ft$p.value,or=ft$estimate)})

z = ct[[1]]
for(i in 2:7) z = z + ct[[1]]
fisher.test(z)
# same for heart (exclude brain)
ct = lapply(rownames(species),function(s){
	f = !(colnames(per.tissue.age.qv[[s]]) %in% c('brain','cerebellum'))
	sgn = apply(per.tissue.age.qv[[s]][,f]<0.05 & abs(age.dpsi[[s]][,f]) > 0.2,1,sum,na.rm=T)>0 & anns[[s]]$sites=='ad'
	b = per.tissue.age.qv[[s]][,'heart']<0.05 & abs(age.dpsi[[s]][,'heart']) > 0.2
	b[is.na(b)] = FALSE
	m = anns[[s]]$length<=27
	table(brain=b[sgn],micro=m[sgn])
})
sapply(ct,function(x){ft=fisher.test(x);c(pv=ft$p.value,or=ft$estimate)})

z = ct[[1]]
for(i in 2:7) z = z + ct[[1]]
fisher.test(z)$p.value


# 7 ######
# _prepare #####
exon.birth.one = readRDS('Rdata/exon.birth.one.Rdata')
born.per.tissue.age.qv = readRDS('Rdata/born.per.tissue.age.qv.Rdata')
born.exn.sajr = readRDS('Rdata/born.exn.sajr.Rdata')


load('Rdata/tmp.exon.birth.Rdata')
len =     unlist(lapply(exon.birth.one,function(x)(min(x$length,na.rm=T))))
nb.stat[1:2,]
dim(nb.stat)
nbf = nb.stat$adj.exons== 0 & nb.stat$full.obs==1 & len < 500

sps=list(m = 'm',mr='mr',mrb='mrb',hqmrb='hqmrb',hqmrbo='hqmrbo')
sps. = c(sps,hqmrboc='hqmrboc')

sp.birth = sapply(exon.birth.one,function(x){paste(species$short[!is.na(x$seg_id)],collapse='')})
born.seg.ids = t(sapply(exon.birth.one,function(x){x$seg_id}))
colnames(born.seg.ids) = rownames(exon.birth.one[[1]])

born.exn.tsm = lapply(born.exn.sajr,function(x){
	psi = x$ir[,colnames(x$ir) %in% rownames(meta)]
	m = meta[colnames(psi),]
	calcMeanCols(psi,paste(m$species,m$tissue,m$stage))
})

sapply(born.exn.sajr,function(x)mean(x$seg$type=='ALT'))

sapply(rownames(species),function(s)mean(sapply(exon.birth.one[sp.birth==species[s,'short']],function(x)x[s,'type'])=='ALT'))
a = exon.birth.one[sp.birth=='h']
table(sapply(a,function(x)x['human','type']))

# check which species/tissue/stages are present in DS
sts = matrix('',ncol=7,nrow=nrow(age.al.i),dimnames = list(age.al.i$mouse,rownames(species)))
for(s in rownames(species))
	for(t in unique(meta$tissue))
		for(a in 1:nrow(sts))
			if(paste(s,t,age.al.i[a,s]) %in% rownames(meta.tsm))
				sts[a,s] = paste0(sts[a,s],substr(t,1,1))
		
apply(sts[age.al.i$to.remove==0,-c(2,7)],1:2,function(x)grepl('b',x) & grepl('h',x)& grepl('l',x)& grepl('t',x)& grepl('t',x))

species.sps = c(human='bhlt',mouse='bhlt',rat='bhlt',opossum='bhlt',rabbit='bhlt')
phylo.groups = c(species$short,'hq','mr','mrb','hqmrb','hqmrbo')

max.stages = list()
for(s in names(species.sps))
	max.stages[[s]] = getMaxStage(born.exn.tsm[[s]],species.sps[s],s,age.al.i[age.al.i$to.remove==0,])
max.stages.stat = lapply(max.stages,getMaxStageStat,born.seg.ids[nbf,],na.thr=Inf)

# pdf('figures/paper.figures/5/7/7B-D.pdf',w=9,h=9,family='Arial')
# par(mfcol=c(5,5),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(4,1.2,1.2,1),oma=c(0,0,0,1))
# for(s in names(max.stages.stat)){
# 	sps = phylo.groups[grep(species[s,'short'],phylo.groups)]
# 	sps = sps[order(nchar(sps))]
# 	z = max.stages.stat[[s]]$tissue[sps,]
# 	barplot(t(sweep(z,1,apply(z,1,sum),'/')),col=params$tissue.col[colnames(z)],main=s,border=NA,las=3)
# 	z = max.stages.stat[[s]]$all[sps,]
# 	ts = unique(meta$tissue)
# 	ts = setNames(ts,substr(ts,1,1))[strsplit(species.sps[s],'')[[1]]]
# 	for(t in ts){
# 		col = col2rgb(params$tissue.col[t])
# 		col=rgb(col[1],col[2],col[3],seq(255,100,length.out = nrow(z)),maxColorValue = 255)
# 		barplot(z[,paste(s,t,colnames(max.stages.stat[[s]]$stage))],main=t,border = NA,names.arg = colnames(max.stages.stat[[s]]$stage),las=3,col=col)
# 	}
# 	col=rgb(0,0,0,seq(255,100,length.out = nrow(z)),maxColorValue = 255)
# 	legend('topleft',fill=col,legend = rownames(z),bty='n')
# }
# dev.off()

orth.seg.ad.all.id = readRDS('Rdata/orth.seg.ad.all.id.Rdata')
orth.seg.ad.all.id = orth.seg.ad.all.id[names(alt.sp),]
orth.seg.ad.all.tsm = readRDS('Rdata/orth.seg.ad.all.tsm.Rdata')

min.stages = list()
for(s in names(species.sps))
	min.stages[[s]] = getMaxStage(orth.seg.ad.all.tsm[[s]],species.sps[s],s,age.al.i[age.al.i$to.remove==0,],FUN=min)

# alternif = list(by.ann = alt.sp,
# 								psi95.2 = getAltSp(orth.seg.ad.all.tsm,0.95,2),
# 								psi95.4 = getAltSp(orth.seg.ad.all.tsm,0.95,4),
# 								psi95.8 = getAltSp(orth.seg.ad.all.tsm,0.95,8),
# 								psi90.2 = getAltSp(orth.seg.ad.all.tsm,0.9 ,2),
# 								psi90.4 = getAltSp(orth.seg.ad.all.tsm,0.9 ,4),
# 								psi90.8 = getAltSp(orth.seg.ad.all.tsm,0.9 ,8),
# 								psi80.2 = getAltSp(orth.seg.ad.all.tsm,0.8 ,2),
# 								psi80.4 = getAltSp(orth.seg.ad.all.tsm,0.8 ,4))
# pdf('figures/paper.figures/5/7/7B-D.alternif.na.thr=Inf.pdf',w=9,h=9,family='Arial')
# sss = names(min.stages)
# par(mfcol=c(5,5),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(4,1.2,1.2,1),oma=c(0,0,1,1))
# for(n in names(alternif)){
# 	print(n)
# 	as = alternif[[n]][names(alt.sp)]
# 	min.stages.stat = lapply(setNames(sss,sss),function(s)getMinStageStat(min.stages[[s]],setNames(as,orth.seg.ad.all.id[,s]),na.thr=Inf))
# 	for(s in names(min.stages.stat)){
# 		sps = phylo.groups[grep(species[s,'short'],phylo.groups)]
# 		sps = sps[order(nchar(sps))]
# 	
# 		z = min.stages.stat[[s]]$tissue[sps,]
# 		barplot(t(sweep(z,1,apply(z,1,sum),'/')),col=params$tissue.col[colnames(z)],main=s,border=NA,las=3)
# 		z = min.stages.stat[[s]]$all[sps,]
# 		ts = unique(meta$tissue)
# 		ts = setNames(ts,substr(ts,1,1))[strsplit(species.sps[s],'')[[1]]]
# 		for(t in ts){
# 			col = col2rgb(params$tissue.col[t])
# 			col=rgb(col[1],col[2],col[3],seq(255,100,length.out = nrow(z)),maxColorValue = 255)
# 			barplot(z[,paste(s,t,colnames(min.stages.stat[[s]]$stage))],main=t,border = NA,names.arg = colnames(min.stages.stat[[s]]$stage),las=3,col=col)
# 		}
# 		col=rgb(0,0,0,seq(255,100,length.out = nrow(z)),maxColorValue = 255)
# 		legend('topleft',fill=col,legend = rownames(z),bty='n')
# 	}
# 	mtext(n,3,outer=T)
# }
# dev.off()

# as = alternif$psi90.4[names(alt.sp)]
# m=min.stages$mouse$max.stages
# sort(table(as[!is.na(m$max.ts) & m$max.ts=='mouse testis 0dpb']))
# 
# z=m[!is.na(m$max.ts) & m$max.ts=='mouse testis 0dpb' & as =='m',]
# sum(!is.na(m$max.ts) & m$max.ts=='mouse testis 0dpb' & as =='m')
# dim(z)
# sort(table(as[!is.na(m$max.ts) & m$max.ts=='mouse testis 0dpb']))
# z[1:10,]
# plotTissueAgeProile(orth.seg.ad.all.tsm$mouse['mou.19095.s16',],meta.tsm)
sss = names(min.stages)
min.stages.stat = lapply(setNames(sss,sss),function(s)getMinStageStat(min.stages[[s]],setNames(alt.sp,orth.seg.ad.all.id[,s])))

# z=getAltSp(orth.seg.ad.all.tsm,0.9 ,4)
# plotTissueAgeProile(apply(orth.seg.ad.all.tsm$mouse[z=='m',],2,mean,na.rm=T),meta.tsm,age.axis = 'rank')

#plotTissueAgeProile(apply(born.exn.tsm$mouse,2,mean,na.rm=T),meta.tsm,age.axis = 'rank')

# devAS
sps=list(m = 'm',mr='mr',mrb='mrb',hqmrb='hqmrb',hqmrbo='hqmrbo')
sps. = c(sps,hqmrboc='hqmrboc')

#born.tis.dpsi = lapply(rownames(species),function(s)getAgeASchanges(born.exn.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
born.tis.dpsi = lapply(rownames(species),function(s){
	print(s)
	p = born.exn.sajr[[s]]$ir
	p = p[,colnames(p) %in% rownames(meta)[!(meta$stage %in% c('oldermidage','senior'))]]
	m = meta[colnames(p),]
	r = matrix(NA,nrow=nrow(p),ncol=7,dimnames = list(rownames(p),unique(meta$tissue)))
	for(t in colnames(r)){
		cat('  ',t,'\t')
		if(sum(m$tissue==t)>0){
			z = p[,m$tissue==t]
			m. = m[colnames(z),]
			r[,t] = apply(z,1,function(y)getDiamBySpline(m.$age.use,y,4))
		}
	}
	gc()
	r
})
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

# orth.seg.ad.all.dpsi = lapply(rownames(species),function(s)getAgeASchanges(orth.seg.ad.all.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
# names(orth.seg.ad.all.dpsi) = rownames(species)

alt.devAS = sapply(names(orth.per.tissue.age.qv),function(s){apply(orth.per.tissue.age.qv[[s]][,]<0.05 & abs(orth.age.dpsi[[s]])>0.2,1,sum,na.rm=T)>0})

born.ex.prop.sgn.dpsi0.2 = t(apply(table(sp.birth[nbf],born.devAS[nbf,'mouse'])[names(sps),c('TRUE','FALSE')],1,my.binom.test))
altern.prop.sgn.dpsi0.2  = t(apply(table( alt.sp ,alt.devAS[ ,'mouse'])[names(sps.),c('TRUE','FALSE')],1,my.binom.test))

# _devAS frac in ancient ######
alt.devAS.tissues = lapply(names(orth.per.tissue.age.qv),function(s){orth.per.tissue.age.qv[[s]]<0.05 & abs(orth.age.dpsi[[s]])>0.2})
names(alt.devAS.tissues) = names(orth.per.tissue.age.qv)
alt.devAS.tissues = lapply(alt.devAS.tissues,function(x){x[is.na(x)] = FALSE;x})

t=sapply(rownames(species),function(s)	mean(apply(alt.devAS.tissues[[s]][alt.sp=='hqmrboc',],1,sum)>0))
writeLines(paste(names(t),round(t,2)))
sapply(rownames(species)[-2],function(s)	apply(alt.devAS.tissues[[s]][alt.sp=='hqmrboc',],2,mean))

t=sapply(rownames(species),function(s){
	sgn = per.tissue.age.qv[[s]] < 0.05 & abs(age.dpsi[[s]])>0.2
	#sgn[is.na(sgn)] = FALSE
	c(apply(sgn[anns[[s]]$sites=='ad',],2,mean,na.rm=T),any=mean(apply(sgn[anns[[s]]$sites=='ad',],1,sum,na.rm=T)>0),tested=sum(apply(!is.na(sgn[anns[[s]]$sites=='ad',]),1,sum,na.rm=T)>0))
	})
writeLines(paste(colnames(t),round(t[8,],2),t[9,]))
t

t='liver'
mean(apply(do.call(cbind,lapply(alt.devAS.tissues,function(x)x[alt.sp=='hqmrboc',t])),1,sum)>0)

mean(apply(do.call(cbind,alt.devAS.tissues)[alt.sp=='hqmrboc',],1,sum,na.rm=T)>0)

sapply(rownames(species),function(s)	mean(apply(alt.devAS.tissues[[s]][alt.sp!='hqmrboc',],1,sum)>0))
sapply(rownames(species),function(s)	mean(apply(alt.devAS.tissues[[s]][alt.sp!='hqmrboc' & grepl(species[s,'short'],alt.sp),],1,sum)>0))
mean(apply(do.call(cbind,alt.devAS.tissues)[alt.sp!='hqmrboc',],1,sum,na.rm=T)>0)
sapply(rownames(species),function(s)mean(apply(age.segs[[s]]!='-' & age.segs[[s]]!='n',1,sum)[anns[[s]]$sites=='ad']>0))


### mean PSI
s2s = setNames(rownames(species),species$short)
orth.mean.psi = sapply(orth.seg.ad.all.tsm,function(x)apply(x,1,mean,na.rm=T))
alt.psi.on.ev = split(orth.mean.psi[names(alt.sp),'mouse'],alt.sp)[names(sps.)]

alt.psi.on.ev=t(sapply(alt.psi.on.ev,function(t){
	t = t[!is.na(t)]
	m = mean(t)
	s = sd(t)/sqrt(length(t))*3
	c(mean=m,ci1=m-s,m+s)
}))


brn.psi.on.ev = split(born.seg.ids[nbf,'mouse'],sp.birth[nbf])[names(sps)]

mouse.born.mean.psi = apply(born.exn.sajr$mouse$ir,1,mean,na.rm=T)
 

brn.psi.on.ev=t(sapply(brn.psi.on.ev,function(s){
	t=mouse.born.mean.psi[s]
	t = t[!is.na(t)]
	m = mean(t)
	s = sd(t)/sqrt(length(t))*3
	c(mean=m,ci1=m-s,m+s)
}))

## %% 3 and cod
f = !is.na(born.seg.ids[,'mouse']) & nbf
mouse.sp.birth = setNames(sp.birth[f],born.seg.ids[f,'mouse'])

seg.len = setNames(born.exn.sajr$mouse$seg$length,rownames(born.exn.sajr$mouse$seg))#unlist(lapply(setNames(born.exn.sajr,NULL),function(x)setNames(x$seg$length,rownames(x$seg))))
seg.cod = setNames(born.exn.sajr$mouse$seg$cod!='n',rownames(born.exn.sajr$mouse$seg))#unlist(lapply(setNames(born.exn.sajr,NULL),function(x)setNames(x$seg$cod=='c',rownames(x$seg))))
born.exn.len3 = t(sapply(split(seg.len %% 3 == 0,mouse.sp.birth[names(seg.len)])[names(sps)],function(x)my.binom.test(sum(x),sum(!x)))) #sapply(sps,function(s){r=seg.cod[unique(unlist(born.sids[s]))];my.binom.test(sum(r),sum(!r))})
born.exn.cod = t(sapply(split(seg.cod,mouse.sp.birth[names(seg.cod)])[names(sps)],function(x){my.binom.test(sum(x),sum(!x))}))#sapply(sps,function(s){r=seg.len[unique(unlist(born.sids[s]))] %% 3 == 0;my.binom.test(sum(r),sum(!r))})

#orth.seg.ad.all = readRDS('Rdata/orth.seg.ad.all.Rdata')
seg.cod = readRDS('Rdata/paper.figures/alt.seg.cod.Rdata')#seg.cod = unlist(lapply(setNames(orth.seg.ad.all,NULL),function(x)setNames(x$seg$cod=='c',rownames(x$seg))))
#saveRDS(seg.cod,'Rdata/paper.figures/alt.seg.cod.Rdata')
seg.len = readRDS('Rdata/paper.figures/alt.seg.len.Rdata')#seg.len=unlist(lapply(setNames(orth.seg.ad.all,NULL),function(x)setNames(x$seg$length,rownames(x$seg))))
#saveRDS(seg.len,'Rdata/paper.figures/alt.seg.len.Rdata')
alt.exn.len3 = t(sapply(sps.,function(s){r=seg.len[rownames(orth.per.tissue.age.qv$mouse)[alt.sp %in% s]] %% 3 == 0;my.binom.test(sum(r),sum(!r))}))
alt.exn.cod = t(sapply(sps.,function(s){r=seg.cod[rownames(orth.per.tissue.age.qv$mouse)[alt.sp %in% s]];my.binom.test(sum(r),sum(!r))}))

# attempts to look on dev dynamics of new exons
# getNewExonInclExclStat = function(s){
# 	sps = c(species$short,'hq','mr','mrb','hqmrb')
# 	sids = born.seg.ids[nbf & sp.birth %in% sps[grep(species[s,'short'],sps)],s]
# 	rbind(up=apply(born.per.tissue.age.qv[[s]][sids,]<0.05 & born.tis.dpsi[[s]][sids,] >  0.2,2,sum,na.rm=T),
# 				down=apply(born.per.tissue.age.qv[[s]][sids,]<0.05 & born.tis.dpsi[[s]][sids,] < -0.2,2,sum,na.rm=T),
# 				total = apply(born.per.tissue.age.qv[[s]],2,length))
# }
# getNewExonInclExclStat('human')
# f = function(s,t){
# 	sps = c(species$short,'hq','mr','mrb','hqmrb')
# 	sps = sps[grep(species[s,'short'],sps)]
# 	type = rep('other',nrow(psi.tsm[[s]]))
# 	type[rownames(psi.tsm[[s]]) %in% rownames(orth.age.dpsi[[s]])] = 'ort'
# 	type[rownames(psi.tsm[[s]]) %in% born.seg.ids[grep(species[s,'short'],sp.birth),s]] = 'new.'
# 	type[rownames(psi.tsm[[s]]) %in% born.seg.ids[nbf & sp.birth %in% sps,s]] = 'new'
# 	type[rownames(psi.tsm[[s]]) %in% rownames(orth.age.dpsi[[s]])[grep(species[s,'short'],alt.sp)]] = 'alt.'
# 	type[rownames(psi.tsm[[s]]) %in% rownames(orth.age.dpsi[[s]])[alt.sp %in% sps]] = 'alt'
# 	type[rownames(psi.tsm[[s]]) %in% rownames(orth.age.dpsi[[s]])[alt.sp == 'hqmrboc']] = 'anc'
# 	d = rep('n',nrow(psi.tsm[[s]]))
# 	f = per.tissue.age.qv[[s]][,t] < 0.05 & age.dpsi[[s]][,t] > 0.2
# 	d[!is.na(f) & f] = 'u'
# 	f = per.tissue.age.qv[[s]][,t] < 0.05 & age.dpsi[[s]][,t] < -0.2
# 	d[!is.na(f) & f] = 'd'
# 	table(type,d)
# }

# _plot ######
pdf('figures/paper.figures/6/7/7a.v2.pdf',w=7.2,h=7.2/12*11,family='Arial')
#png('figures/paper.figures/6/7/7.png',units = 'in',w=7.2,h=7.2/12*11,res = 600,family='Arial')
# layout(matrix(c(1,1,2,4,
# 								1,1,3,5,
# 								1,1,6,8,
# 								1,1,7,8,
# 								9:12),ncol=4,byrow = T),heights = c(1,1,1,1,2))

layout(matrix(c(1,2,4,
								1,3,5,
								1,11,6,
								1,11,7,
								8:10),ncol=3,byrow = T),heights = c(2,2,2,2,3))

par(tck=-0.01,mgp=c(1.2,0.3,0),mar=c(5.5,.5,6,0.5),oma=c(0,0,0,1))
#plot7A('A',nb.sp = 'hqmrb',alt.sp = 'hq',substitute(atop('NA',atop("Eutherian-specific",paste("exon in ",italic('APP'))))),substitute(atop("Primate-specific",atop("alternification in",italic('AMPD2')))))
plot7A('A',nb.sp = 'hqmrb',alt.sp = 'hq',c("Eutherian-",'specific',expression(paste("exon in ",italic('APP')))),c("Primate-",'specific',"alternification",expression(paste('in ',italic('AMPD2')))))
#early late
lm = 2.5
rm = 1
um = c(3.5,1)
dm = c(0,2.5)

par(xpd=NA,mar=c(dm[1],lm,um[1],rm))
plotNewExonDevPref('brain',xaxt='n',ylab='% of exons')
title(main='New exons',line=2)
legend(0,110,fill = c('#00000055','#000000'),legend = c('early','late'),border = NA,bty='n',ncol=2,x.intersp=0.5,cex=0.85)
plotPNG('figures/paper.figures/5/icons/brain.png',0.065,1.02,0.10)
plotPanelLetter('B',lab.cex)
par(mar=c(dm[2],lm,um[2],rm))
plotNewExonDevPref('testis',plot.species.icons = T,xaxt='n',ylab='% of exons')
plotPNG('figures/paper.figures/5/icons/testis.png',0.065,1.02,0.10)

# new exon on evo
par(mar=c(dm[1],lm,um[1],rm))
plot7B(max.stages.stat,'brain',main='',ylab='% of exons',xaxt='n')
title(main='New exons',line=2)
legend(0.2,80,fill = c('#00000055','#000000'),legend = c('species-spec.','eutherian-spec.'),ncol=1,border = NA,bty='n',xpd=T,x.intersp=0.5,y.intersp=0.8,cex=0.85)
plotPanelLetter('C',lab.cex)
par(mar=c(dm[2],lm,um[2],rm))
plot7B(max.stages.stat,'testis',ylab='% of exons',xaxt='n',plot.species.icons=T)


# alt exons on evo
par(mar=c(dm[1],lm,um[1],rm))
plot7B(min.stages.stat,'brain',xaxt='n',ylab='% of exons')
title(main='Alternification',line=2)
legend(0.2,80,fill = c('#00000055','#000000'),legend = c('species-spec.','eutherian-spec.'),ncol=1,border = NA,bty='n',xpd=T,x.intersp=0.5,y.intersp=0.8,cex=0.85)
plotPanelLetter('H',lab.cex)
par(mar=c(dm[2],lm,um[2],rm))
plot7B(min.stages.stat,'testis',xaxt='n',plot.species.icons=T,ylab='% of exons')

#plotTissueAgeProile(apply(max.stages.stat$mouse$all[,],2,sum),meta.tsm,age.axis = 'rank',ylab='# of exons',bty='n')
#plotTissueAgeProile(apply(max.stages.stat$mouse$all[c('m','mr','mrb','hqmrb'),],2,sum),meta.tsm,age.axis = 'rank',ylab='# of exons',bty='n')
#plotPNG('figures/paper.figures/5/icons/mouse.png',0.25,0.7,0.25)

fig.y = 3/11
par(fig=c(0,1/4,0,fig.y),new = TRUE,xpd=F,mar=c(4,3,1.5,1))
plotArea(1:length(sps),p = born.ex.prop.sgn.dpsi0.2,col='red',lwd=3,new = T,xlim=c(1,length(sps.)),ylim=range(0,born.ex.prop.sgn.dpsi0.2,altern.prop.sgn.dpsi0.2),xaxt='n',xlab='',ylab='proportion of devAS',main='',bty='n')
plotArea(1:length(sps),p = altern.prop.sgn.dpsi0.2[-length(sps.),],col='blue',lwd=3,new = F)
points(length(sps.),altern.prop.sgn.dpsi0.2[length(sps.),1],pch=19,col='orange')
segments(length(sps.),altern.prop.sgn.dpsi0.2[length(sps.),3],length(sps.),altern.prop.sgn.dpsi0.2[length(sps.),2],col='orange')
axis(1,1:length(sps.),rownames(altern.prop.sgn.dpsi0.2),las=3)
#abline(h=0,lty=2)
plotPanelLetter('D',lab.cex)

par(fig=c(1/4,2/4,0,fig.y),new = TRUE)
plotArea(1:length(sps),alt.psi.on.ev[-length(sps.),],col='blue',new=T,ylim=c(0,1),lwd=3,xaxt='n',xlab='',ylab='mean PSI',main='',xlim=c(1,length(sps.)),bty='n')
plotArea(1:length(sps),brn.psi.on.ev[-length(sps.),],col='red',lwd=3,new=F)
points(length(sps.),alt.psi.on.ev[length(sps.),1],pch=19,col='orange')
segments(length(sps.),alt.psi.on.ev[length(sps.),2],length(sps.),alt.psi.on.ev[length(sps.),3],col='orange')
axis(1,1:length(sps.),rownames(alt.psi.on.ev),las=3)
plotPanelLetter('E',lab.cex)

par(fig=c(2/4,3/4,0,fig.y),new = TRUE)
x=1:length(sps)
plotArea(x,born.exn.len3,col='red',new = T,type='l',lwd=3,ylim=range(born.exn.len3,alt.exn.len3),xlim=c(1,length(sps.)),xaxt='n',xlab='',main='',ylab='proportion of 3N exons',bty='n')
plotArea(x,alt.exn.len3[x,],col='blue',type='l',lwd=3)
points(length(sps.),alt.exn.len3[length(sps.),1],pch=19,col='orange')
segments(length(sps.),alt.exn.len3[length(sps.),2],length(sps.),alt.exn.len3[length(sps.),3],col='orange')
axis(1,1:nrow(alt.exn.len3),rownames(alt.exn.len3),las=3)
abline(h=1/3,lty=2)
plotPanelLetter('F',lab.cex)

par(fig=c(3/4,4/4,0,fig.y),new = TRUE)
plotArea(x,born.exn.cod,col='red',new = T,type='l',lwd=3,ylim=range(born.exn.cod,alt.exn.cod),xlim=c(1,length(sps.)),xaxt='n',xlab='',main='',ylab='proportion of coding exons',bty='n')
plotArea(x,alt.exn.cod[x,],col='blue',type='l',lwd=3)
points(length(sps.),alt.exn.cod[length(sps.),1],pch=19,col='orange')
segments(length(sps.),alt.exn.cod[length(sps.),2],length(sps.),alt.exn.cod[length(sps.),3],col='orange')
axis(1,1:nrow(alt.exn.len3),rownames(alt.exn.len3),las=3)
plotPanelLetter('G',lab.cex)
# plot.new()
# legend('topleft',fill=c('red','blue'),legend=c('new exons','alternification'),bty='n',border=NA)
dev.off()



# Suppl #####
# _2 #####
cols = list(const='gray',alt='orange')
sp=2
wd=6
pdf('figures/paper.figures/6/suppl/new.order/S2a.pdf',w=7,h=3)
#jpeg('figures/paper.figures/6/suppl/S2.jpg',units = 'in',w=7,h=3,quality = 100,res = 600)
par(mfrow=c(1,2),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(0,0,0,0),oma=c(0,0,0,0))
plotSegDef('A',cols,sp,wd)
plotSAJR.AS.Q('B',wd)
dev.off()

# _3 #####
detASboot.sub = readRDS('Rdata/hqmrboc.subsample/detected.as.bootstrap.sub.Rdata')
detASboot.sub.orth = readRDS('Rdata/hqmrboc.subsample/detected.as.bootstrap.sub.orth.Rdata')

all.stat = getSegTestDevAsStat(age.segs,anns)

#_compare with direct detAS definition ######
r = array(NA,dim=c(7,7,4),dimnames=list(rownames(species),unique(meta$tissue),c('ad','da','dd','aa')))
for(s in rownames(species)){
	print(s)
	d = readRDS(paste0('Rdata/',s,'.as.u.all.Rdata'))
	d = d[d$seg$chr_id !='MT' & d$seg$position=='INTERNAL' & d$seg$type!='EXN',colnames(d$ir) %in% rownames(meta)]
	m = meta[colnames(d$ir),]
	na = !is.na(d$ir)
	
	f = rep(FALSE,length(d))
	na = is.na(d$ir)
	for(t in unique(meta$tissue)){
		cinx = m$tissue==t
		f = f | (apply(!na[,cinx],1,mean) > 0.6 & apply(d$ir[,cinx],1,function(x){x=x[!is.na(x)];sum(x>0.1 & x < 0.9)>3}))
	}
	print(sum(f))
	d = d[f,]
	na = na[f,]
	for(t in unique(meta$tissue)){
		for(ss in dimnames(r)[[3]]){
			r[s,t,ss] = sum(apply(!na[d$seg$sites==ss,m$tissue==t],1,mean)>0.6)
		}
	}
}
#saveRDS(r,'Rdata/hqmrboc.subsample/detected.as.all.direct.Rdata')
r = readRDS('Rdata/hqmrboc.subsample/detected.as.all.direct.Rdata')
plotLine(all.stat$tested$ad,t(r[,,'ad']))
all.stat$tested$ad - t(r[,,'ad'])
plotAsEventCount(t(r[,,'ad']),astypes.pchs[1],by.tissue = F,ylab='# of detected events',main='Orth: Detected AS',bty='n')
plotAsEventCount(all.stat$tested$ad,astypes.pchs[1],by.tissue = F,ylab='# of detected events',main='Orth: Detected AS',bty='n')

# what is wrong with liver? #####
library(reshape)
m = readRDS(paste0('Rdata/mouse.as.u.filtered.Rdata'))
h = readRDS(paste0('Rdata/human.as.u.filtered.Rdata'))
d = m
d = d[d$seg$sites=='ad',]
m = meta[colnames(d$ir),]
table(m$tissue)
na.cnt=do.call(rbind,lapply(unique(m$tissue),function(t){data.frame(nna=apply(is.na(d$ir[,m$tissue==t]),1,sum),tissue=t)}))
table(na.cnt$nna,na.cnt$tissue)
m$na.cnt=apply(is.na(d$ir),2,sum)

# cnt=cast(m,tissue ~ days,value='na.cnt',fun.aggregate=mean)
# barplot(t(as.matrix(cnt)),beside = T)

plotTissueAgeProile(setNames(m$na.cnt,rownames(m)),m,age.axis='rank')



pdf('figures/paper.figures/6/suppl/new.order/S3a.v2.pdf',w=9,h=6,family='Arial')
astypes.pchs=c(ad=19,aa=2,dd=6,da=13)
#jpeg('figures/paper.figures/6/suppl/S3.jpg',units = 'in',w=7,h=3,quality = 100,res = 600,family='Arial')
#layout(matrix(c(1:4,6,5),ncol=3,byrow = T),widths = c(3.1,1,3.1,3.1))
par(mfrow=c(2,3),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(1.5,2.5,1.5,0),oma=c(0,0,0,0))
# order by species/tissue
plotAsEventCount(all.stat$tested,astypes.pchs,by.tissue = F,ylab='# of detected events',main='All AS - all samples',bty='n')
plotPanelLetter('A',lab.cex)
plotAsEventCount(all.stat$devasn,astypes.pchs,by.tissue = F,ylab='# of devAS',main='DevAS - all samples',bty='n')
plotPanelLetter('B',lab.cex)
plotAsEventCount(all.stat$devasf,astypes.pchs,by.tissue = F,ylab='% of devAS',main='DevAS - all samples',bty='n')
plotPanelLetter('C',lab.cex)

plotAsEventCount.boot(detASboot.sub,astypes.pchs[1],by.tissue = FALSE,ylab='# of detected events',main='All AS - subsampling',bty='n')
plotPanelLetter('D',lab.cex)
plotAsEventCount.boot(detASboot.sub.orth,astypes.pchs[1],by.tissue = FALSE,ylab='# of detected events',main='All AS - subsampling & orthologous exons',bty='n')
plotPanelLetter('E',lab.cex)

par(mar=c(2.5,2,1.5,12))
y=plotASTypes(astypes.pchs)
legend(0,y-2,pch=19,col=params$tissue.col,legend = names(params$tissue.col),bty='n',cex = 1,xjust=0,yjust=1,xpd=T)
rect(3,-1,100,100,xpd=NA,col=NA,border='gray',xpd=T)
dev.off()

# __in orth segs ######
t = getSegTestDevAsStat(age.segs,anns)
o = getSegTestDevAsStat(age.segs,anns,sapply(orth.seg.ad.tsm,rownames))
o$devasn
by.tissue = F

pdf('figures/paper.figures/6/suppl/new.order/S3a_orth.pdf',w=9,h=6,family='Arial')
par(mfrow=c(2,3),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(1.5,2.5,1.5,0),oma=c(0,0,0,0))
for(by.tissue in c(FALSE,TRUE)){
	plotAsEventCount(t$tested,astypes.pchs[1],by.tissue = by.tissue,ylab='# of detected events',main='All: Detected AS',bty='n')
	plotPanelLetter('A',lab.cex)
	plotAsEventCount(t$devasn,astypes.pchs[1],by.tissue = by.tissue,ylab='# of devAS',main='All: DevAS',bty='n')
	plotPanelLetter('B',lab.cex)
	plotAsEventCount(t$devasf,astypes.pchs[1],by.tissue = by.tissue,ylab='% of devAS',main='All: DevAS',bty='n')
	plotPanelLetter('C',lab.cex)
	
	plotAsEventCount(o$tested,astypes.pchs[1],by.tissue = by.tissue,ylab='# of detected events',main='Orth: Detected AS',bty='n')
	plotPanelLetter('D',lab.cex)
	plotAsEventCount(o$devasn,astypes.pchs[1],by.tissue = by.tissue,ylab='# of devAS',main='Orth: DevAS',bty='n')
	plotPanelLetter('E',lab.cex)
	plotAsEventCount(o$devasf,astypes.pchs[1],by.tissue = by.tissue,ylab='% of devAS',main='Orth: DevAS',bty='n')
	plotPanelLetter('F',lab.cex)
}
dev.off()
# _v3 ######
# devas.boot = readRDS('Rdata/hqmrboc.subsample/devAS.subsampling.summary.Rdata')
# devas.boot. = readRDS('Rdata/hqmrboc.subsample/devAS.subsampling.summary.chicken.12-0dpb.Rdata')
# devas.boot[names(devas.boot.)] = devas.boot.
nrow=4

pdf('figures/paper.figures/6/suppl/new.order/S3a.v3.chicken.12dpc-to-mouse.0dpb.pdf',w=9,h=12,family='Arial')
astypes.pchs=c(ad=19,aa=2,dd=6,da=13)
#jpeg('figures/paper.figures/6/suppl/S3.jpg',units = 'in',w=7,h=3,quality = 100,res = 600,family='Arial')
#layout(matrix(c(1:4,6,5),ncol=3,byrow = T),widths = c(3.1,1,3.1,3.1))
par(fig=c(0,0.4,3/nrow,1),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(1.8,2.5,1.5,0),oma=c(0,0,0,0),cex=2.1/3,bty='n')
# order by species/tissue

plotAsEventCount(all.stat$tested,astypes.pchs,by.tissue = F,ylab='# of detected events',main='All AS - all samples',bty='n',xlab.cex = 1)
plotPanelLetter('A',lab.cex)
par(fig=c(0.4,0.8,3/nrow,1),new=T)
plotAsEventCount(all.stat$devasn,astypes.pchs,by.tissue = F,ylab='# of devAS',main='DevAS - all samples',bty='n',xlab.cex = 1)
plotPanelLetter('B',lab.cex)

par(mar=c(2.5,2,1.5,2),fig=c(0.8,1,3/nrow,1),new=T)
y=plotASTypes(astypes.pchs,cex=0.85)
legend(0,y-2,pch=19,col=params$tissue.col,legend = names(params$tissue.col),bty='n',cex = 1,xjust=0,yjust=1,xpd=T)
rect(3,-5,100,100,xpd=NA,col=NA,border='gray',xpd=T)

par(fig=c(0,0.25,2/nrow,3/nrow),mar=c(1.8,2.5,1.5,0),new=T)
plotAsEventCount(all.stat$tested,astypes.pchs['ad'],by.tissue = F,ylab='# of detected events',main='All AS - all samples',bty='n',short.species.lab=T,xlab.cex=1.1)
plotPanelLetter('C',lab.cex)

par(fig=c(0.25,0.5,2/nrow,3/nrow),new=T)
plotAsEventCount(all.stat$devasn,astypes.pchs['ad'],by.tissue = F,ylab='# of devAS',main='DevAS - all samples',bty='n',short.species.lab=T,xlab.cex=1.1)
#plotPanelLetter('D',lab.cex)

r = devas.boot$hq$all$ad$det
par(fig=c(0.5,0.75,2/nrow,3/nrow),new=T)
plotAsEventCount(r$median,cil=r$cil,cih=r$cih,main='All AS - hq subsampling',ylab='# of detected events',short.species.lab=F,xlab.cex=1.1)
plotPanelLetter('D',lab.cex)

r = devas.boot$hq$all$ad$sgn
par(fig=c(0.75,1,2/nrow,3/nrow),new=T)
plotAsEventCount(r$median,cil=r$cil,cih=r$cih,main='DevAS - hq subsampling',ylab='# of detected events',short.species.lab=F,xlab.cex=1.1)
#plotPanelLetter('F',lab.cex)

r = devas.boot$hc$all$ad$det
par(fig=c(0,0.25,1/nrow,2/nrow),new=T)
plotAsEventCount(r$median,cil=r$cil,cih=r$cih,main='All AS - hc subsampling',ylab='# of detected events',short.species.lab=F,xlab.cex=1.1)
plotPanelLetter('E',lab.cex)

r = devas.boot$hc$all$ad$sgn
par(fig=c(0.25,0.5,1/nrow,2/nrow),new=T)
plotAsEventCount(r$median,cil=r$cil,cih=r$cih,main='DevAS - hc subsampling',ylab='# of detected events',short.species.lab=F,xlab.cex=1.1)
#plotPanelLetter('H',lab.cex)


r = devas.boot$hmrbo$all$ad$det
par(fig=c(0.5,0.75,1/nrow,2/nrow),new=T)
plotAsEventCount(r$median,cil=r$cil,cih=r$cih,main='All AS - hmrbo subsampling',ylab='# of detected events',short.species.lab=T,xlab.cex=1.1)
plotPanelLetter('F',lab.cex)

r = devas.boot$hmrbo$all$ad$sgn
par(fig=c(0.75,1,1/nrow,2/nrow),new=T)
plotAsEventCount(r$median,cil=r$cil,cih=r$cih,main='DevAS - hmrbo subsampling',ylab='# of detected events',short.species.lab=T,xlab.cex=1.1)
#plotPanelLetter('J',lab.cex)



r = devas.boot$hqmr$all$ad$det
par(fig=c(0,0.25,0/nrow,1/nrow),new=T)
plotAsEventCount(r$median,cil=r$cil,cih=r$cih,main='All AS - hqmr subsampling',ylab='# of detected events',short.species.lab=F,xlab.cex=1.1)
plotPanelLetter('G',lab.cex)


r = devas.boot$hqc$all$ad$det
par(fig=c(0.25,0.5,0/nrow,1/nrow),new=T)
plotAsEventCount(r$median,cil=r$cil,cih=r$cih,main='All AS - hqc subsampling',ylab='# of detected events',short.species.lab=F,xlab.cex=1.1)
plotPanelLetter('H',lab.cex)

dev.off()

# _v4 #######
pdf('figures/paper.figures/6/suppl/new.order/S3a.v5.chicken.12dpc-to-mouse.0dpb.pdf',w=9,h=9,family='Arial')
nrow=3
astypes.pchs=c(ad=19,aa=2,dd=6,da=13)
#jpeg('figures/paper.figures/6/suppl/S3.jpg',units = 'in',w=7,h=3,quality = 100,res = 600,family='Arial')
#layout(matrix(c(1:4,6,5),ncol=3,byrow = T),widths = c(3.1,1,3.1,3.1))
par(fig=c(0,0.4,(nrow-1)/nrow,1),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(1.8,2.5,1.5,0),oma=c(0,0,0,0),cex=2.1/3,bty='n')
# order by species/tissue

plotAsEventCount(all.stat$tested,astypes.pchs,by.tissue = F,ylab='# of detected events',main='All AS - all samples',bty='n',xlab.cex = 1)
plotPanelLetter('A',lab.cex)
par(fig=c(0.4,0.8,(nrow-1)/nrow,1),new=T)
plotAsEventCount(all.stat$devasn,astypes.pchs,by.tissue = F,ylab='# of devAS',main='DevAS - all samples',bty='n',xlab.cex = 1)
plotPanelLetter('B',lab.cex)

par(mar=c(2.5,2,1.5,2),fig=c(0.8,1,(nrow-1)/nrow,1),new=T)
y=plotASTypes(astypes.pchs,cex=0.85)
legend(0,y-2,pch=19,col=params$tissue.col,legend = names(params$tissue.col),bty='n',cex = 1,xjust=0,yjust=1,xpd=T)
rect(3,-5,100,100,xpd=NA,col=NA,border='gray',xpd=T)

par(fig=c(0,0.25,(nrow-2)/nrow,(nrow-1)/nrow),mar=c(1.8,2.5,1.5,0),new=T)
plotAsEventCount(all.stat$tested,astypes.pchs['ad'],by.tissue = F,ylab='# of detected events',main='All AS - all samples',bty='n',short.species.lab=T,xlab.cex=1.1)
plotPanelLetter('C',lab.cex)

par(fig=c(0.25,0.5,(nrow-2)/nrow,(nrow-1)/nrow),new=T)
plotAsEventCount(all.stat$devasn,astypes.pchs['ad'],by.tissue = F,ylab='# of devAS',main='DevAS - all samples',bty='n',short.species.lab=T,xlab.cex=1.1)
#plotPanelLetter('D',lab.cex)

r = devas.boot$hqc$all$ad$det
par(fig=c(0.5,0.75,(nrow-2)/nrow,(nrow-1)/nrow),new=T)
plotAsEventCount(r$median,cil=r$cil,cih=r$cih,main='All AS - hqc subsampling',ylab='# of detected events',short.species.lab=F,xlab.cex=1.1)
plotPanelLetter('D',lab.cex)

r = devas.boot$hqc$all$ad$sgn
par(fig=c(0.75,1,(nrow-2)/nrow,(nrow-1)/nrow),new=T)
plotAsEventCount(r$median,cil=r$cil,cih=r$cih,main='DevAS - hqc subsampling',ylab='# of detected events',short.species.lab=F,xlab.cex=1.1)
#plotPanelLetter('F',lab.cex)

r = devas.boot$hmrbo$all$ad$det
par(fig=c(0,0.25,(nrow-3)/nrow,(nrow-2)/nrow),new=T)
plotAsEventCount(r$median,cil=r$cil,cih=r$cih,main='All AS - hmrbo subsampling',ylab='# of detected events',short.species.lab=T,xlab.cex=1.1)
plotPanelLetter('E',lab.cex)

r = devas.boot$hmrbo$all$ad$sgn
par(fig=c(0.25,0.5,(nrow-3)/nrow,(nrow-2)/nrow),new=T)
plotAsEventCount(r$median,cil=r$cil,cih=r$cih,main='DevAS - hmrbo subsampling',ylab='# of detected events',short.species.lab=T,xlab.cex=1.1)
#plotPanelLetter('H',lab.cex)

# r = devas.boot$hqmrboc$all$ad$det
# par(fig=c(0.5,0.75,(nrow-3)/nrow,(nrow-2)/nrow),new=T)
# plotAsEventCount(r$median,cil=r$cil,cih=r$cih,main='All AS - hqmrboc subsampling',ylab='# of detected events',short.species.lab=T,xlab.cex=1.1)
# plotPanelLetter('F',lab.cex)
# 
# r = devas.boot$hqmrboc$all$ad$sgn
# par(fig=c(0.75,1,(nrow-3)/nrow,(nrow-2)/nrow),new=T)
# plotAsEventCount(r$median,cil=r$cil,cih=r$cih,main='DevAS - hqmrboc subsampling',ylab='# of detected events',short.species.lab=T,xlab.cex=1.1)
#plotPanelLetter('J',lab.cex)

dev.off()



# _4 ######
age.dpsi.bm = readRDS('Rdata/age.diam.spline4.with.replicates.before.maturation.Rdata')
per.tissue.age.bm.qv = readRDS('Rdata/per.tissue.age.bm.qv.Rdata')

sgn02.stat = sapply(names(anns),function(s)	apply(anns[[s]]$sites == 'ad' & per.tissue.age.qv[[s]] < 0.05 & abs(age.dpsi[[s]])>0.2,2,sum,na.rm=T))
sgn02.stat.bm = sapply(names(anns),function(s)	apply(anns[[s]]$sites == 'ad' & per.tissue.age.bm.qv[[s]] < 0.05 & abs(age.dpsi.bm[[s]])>0.2,2,sum,na.rm=T))


pdf('figures/paper.figures/6/suppl/S4.pdf',w=8,h=4)
#jpeg('figures/paper.figures/6/suppl/S4.jpg',units = 'in',res = 600,w=10,h=4)
par(mfrow=c(1,1),tck=-0.01,mgp=c(1.5,0.5,0),mar=c(2.5,2.5,.5,0),oma=c(0,0,0,1))
cs=rep(params$tissue.col,each=7)
b=barplot(t(sgn02.stat),beside = T,col=paste0(cs,'60'),ylab='# of cassette exons',main='',border=NA,cex.names = 1,ylim=range(0,sgn02.stat,sgn02.stat.bm))
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=1)
b=barplot(t(sgn02.stat.bm),beside = T,col=cs,ylab='# of cassette exons',main='FDR < 0.05; Before sexual maturation',add=T,border=NA,xaxt='n',den=50)
legend(20,3500,fill=c('#00000060','#000000'),legend=c('Whole lifespan','Before sexual maturation'),bty = 'n',den=c(-1,40),border=NA)
dev.off()

# _5 #####
# ens.ge.marg.tsm = readRDS('Rdata/ens.ge.marg.tsm.Rdata')
# use.mean.embryo = FALSE
# c2e = list()
# for(s in rownames(species)){
# 	c2e[[s]] = list()
# 	ge = ens.ge.marg.tsm[[s]]
# 	ge = ge + min(ge[ge!=0],na.rm=T)
# 	ge = log2(ge)
# 	
# 	#sd = apply(psi.tsm[[s]],1,sd)
# 	#c2e[[s]]$psi.cor2embryo = caclCor2Embryo(psi.tsm[[s]][anns[[s]]$sites=='ad' & sd > 0.05,],meta.tsm,cor.m = 'sp',use.mean.embryo=use.mean.embryo)
# 	c2e[[s]]$psi.cor2embryo = caclCor2Embryo(psi.tsm[[s]][anns[[s]]$sites=='ad',],meta.tsm,cor.m = 'p',use.mean.embryo=use.mean.embryo)
# 	c2e[[s]]$ge.cor2embryo = caclCor2Embryo(ge,meta.tsm,cor.m = 'p',use.mean.embryo=use.mean.embryo)
# }


# v2
# registerDoMC(3)
# ens.ge.marg = readRDS('Rdata/ens.ge.marg.Rdata')
# c2e = list()
# for(s in rownames(species)){
# 	print(s)
# 	c2e[[s]] = list()
# 	ge = ens.ge.marg[[s]]
# 	ge = ge + min(ge[ge!=0],na.rm=T)
# 	ge = log2(ge)
# 	
# 	psi.for.ce = readRDS(paste0('Rdata/',s,'.as.u.filtered.Rdata'))
# 	c2e[[s]]$psi.cor2embryo = caclCor2EmbryoAllStat(psi.for.ce$ir[psi.for.ce$seg$sites=='ad',],meta,cor.m = 'p',boot = 1000)
# 	c2e[[s]]$ge.cor2embryo = caclCor2EmbryoAllStat(ge,meta,cor.m = 'p',boot = 1000)
# 	gc()
# }
#saveRDS(c2e,'Rdata/paper.figures/fig.S5.Cor2Embryo.Rdata')
c2e = readRDS('Rdata/paper.figures/fig.S5.Cor2Embryo.Rdata')

cors = readRDS('Rdata/ad.all.cor.Rdata')

pdf('figures/paper.figures/6/suppl/S5.v3.pdf',w=7.2,h=7.2/3*7.4,family='Arial')
#png('figures/paper.figures/6/suppl/S5.v3.png',units = 'in',w=7.2,h=7.2/3*7.4,res = 600,family='Arial')
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

pdf('figures/paper.figures/6/suppl/S6.pdf',w=7.2,h=7.2/3*2,family='Arial')
#png('figures/paper.figures/6/suppl/S6.png',units = 'in',w=7.2,h=7.2/3*2,res = 300,family='Arial')
par(mfrow=c(2,3),tck=-0.01,mgp=c(1.5,0.3,0),mar=c(4,2.5,1.5,0),oma=c(0,0,0,1),las=3)

for(i in 1:3){
	plotTissueAgeProile(h.pe.cor[[i]],meta.tsm,age.axis = 'rank',bty='n',xlab='',ylab=paste0('Correlatin to ',tissues[i]))
	plotPNG("figures/paper.figures/5/icons/human.png",0.15,0.15,0.15)
	plotPNG(paste0("figures/paper.figures/5/icons/",tissues[i],".png"),0.30,0.15,0.15)
#	plotPanelLetter(LETTERS[i],lab.cex)
}

for(i in 1:3){
	plotTissueAgeProile(m.pe.cor[[3+i]],meta.tsm,age.axis = 'rank',bty='n',xlab='',ylab=paste0('Correlatin to ',tissues[i+3]))
	plotPNG("figures/paper.figures/5/icons/mouse.png",0.15,0.15,0.15)
	plotPNG(paste0("figures/paper.figures/5/icons/",tissues[i+3],".png"),0.30,0.15,0.15)
	#plotPanelLetter(LETTERS[3+i],lab.cex)
}
dev.off()

# _8 div-on-age #####
# prev 7
pdf('figures/paper.figures/6/suppl/new.order/S8a.pdf',w=6,h=4.5,family='Arial')
#png('figures/paper.figures/6/suppl/S7.png',units = 'in',w=6,h=4.5,res = 600,family='Arial')
par(tck=-0.01,mgp=c(1.3,0.2,0),mar=c(3,1.2,0,0),oma=c(0,0,0,1),las=3,cex=2/3)
plot.new()
for(i in 1:3){
	s = c('rat','rabbit','opossum')[i]
	plot4C.DivergenceOnAge('mouse',s,as.cor.on.dev[[s]],
												 letters[i],
												 yrange = c((3-i)/3,(4-i)/3),plot.xlab=i==3,ylab.pos=0.3,panel.lab.yadj=1,plot.tissue.lab=i==1)
}
dev.off()


# _9 peak change #####
#prev 8
pdf('figures/paper.figures/6/suppl/new.order/S9a.pdf',w=6,h=5*1.2,family='Arial')
#jpeg('figures/paper.figures/6/suppl/S8.jpg',units = 'in',w=6,h=5*1.2,quality = 100,res = 600,family='Arial')
par(tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,1.2,0,0),oma=c(0,0,0,1),las=3,cex=2/3)
plot.new()
sp.icon.ys = c(0.5,0.9,0.9,0.4,0.8)
for(i in 1:5){
	s = rownames(species)[c(-2,-5)][i]
	plot4D.PeakChange(s,peak.changes[[s]],letters[i],yrange = c((5-i)/5,(6-i)/5),plot.tissue.lab = i==1,plot.xlab = i==5,plot.inset=F,letter.y.pos=(6-i)/5,sp.icon.y=sp.icon.ys[i])
}
dev.off()

# _9 IDR #####
# inclusion/exclusion
# IDR
o = get_OBO('processed/exon.onthology/exont.obo')
seg2exont = readRDS('Rdata/seg2exont.Rdata')

iupr = names(seg2exont)[sapply(seg2exont,function(x)'EXONT:000074' %in% x)]

hdevas02u = age.segs$human == 'u' #& per.tissue.age.qv$human<0.05
#hdevas02u[is.na(hdevas02u)] = FALSE

hdevas02d = age.segs$human == 'd' #& per.tissue.age.qv$human<0.05
#hdevas02d[is.na(hdevas02d)] = FALSE


f = anns$human$sites=='ad' & anns$human$cod!='n' & rownames(anns$human) %in% names(seg2exont)
idr02u = apply(hdevas02u,2,function(x)my.binom.test(table(rownames(hdevas02u)[f & x] %in% iupr)[c('TRUE','FALSE')]))
idr02d = apply(hdevas02d,2,function(x)my.binom.test(table(rownames(hdevas02d)[f & x] %in% iupr)[c('TRUE','FALSE')]))

s = 'human'
sgn = anns$human$sites=='ad' & apply(age.segs[[s]] !='n' & age.segs[[s]] !='-',1,sum) == 0 & anns$human$cod!='n'
alt = rownames(anns$human)[anns[[s]]$cod!='n' & anns[[s]]$sites=='ad' & !sgn &  anns[[s]]$gene_id %in% anns[[s]]$gene_id[sgn] & rownames(anns$human) %in% names(seg2exont)]

#cnst = rownames(all.anns$human)[all.anns[[s]]$cod!='n' & all.anns[[s]]$sites=='ad' & all.anns[[s]]$type=='EXN' & all.anns[[s]]$gene_id %in% anns[[s]]$gene_id[sgn] & rownames(all.anns$human) %in% names(seg2exont)]
cnst = c(0.3685637,0.3623576,0.3748027)#my.binom.test(table(cnst %in% iupr)[c('TRUE','FALSE')])

idr02 = t(cbind(cnst=cnst,"non-devAS"=my.binom.test(table(alt %in% iupr)[c('TRUE','FALSE')]),idr02u,idr02d))
idr02 = idr02[,c(2,1,3)]

pdf('figures/paper.figures/6/suppl/S9a.pdf',w=4,h=4,family='Arial')
#jpeg('figures/paper.figures/6/suppl/S9.jpg',units = 'in',w=4,h=4,quality = 100,res = 600,family='Arial')
par(tck=-0.01,mgp=c(1.3,0.2,0),mar=c(3,2.3,2,1),oma=c(0,0,0,0))
x = c(1:2,seq(4,length.out=7,by=3),seq(5,length.out=7,by=3))
cols=c('gray','black',rep(params$tissue.col,times=2))

plot(x,idr02[,2],col=cols,main="IDR",ylab="fraction of IDR",pch=c(17,17,rep(c(19,1),each=7)),bty='n',xaxt='n',xlab='',cex=2,ylim=c(0.23,0.6))
abline(h=idr02[1:2,2],lty=3,col=c('gray','black'))
arrows(x,idr02[,1],x,idr02[,3],col=cols,angle=90,code=3,length=0.03)
text(c(x[1:2],x[3:9]+0.5),0.23,rownames(idr02)[1:9],adj=c(0,0),srt=-45,xpd=T)
legend('topleft',col='black',pch=c(19,1),bty='n',legend=c('up','down'))
dev.off()

# _10 ####
# Supp. Fig. 10 -> old Fig. 11D (brain and heart panels only). In this figure Pasha is combining exons across all species, 
# which can inflate the p-values. So Pasha is going to see if this same analysis holds if one looks at individual species.
hex.ups = readRDS('Rdata/hex.dws.age02sgn.Rdata')
hex.dws = readRDS('Rdata/hex.ups.age02sgn.Rdata')

iuq = apply(hex.ups$up$pv,2:3,p.adjust,m='BH')
euq = apply(hex.ups$dw$pv,2:3,p.adjust,m='BH')

idq = apply(hex.dws$up$pv,2:3,p.adjust,m='BH')
edq = apply(hex.dws$dw$pv,2:3,p.adjust,m='BH')



f = function(q,thr=0.05,...){
	b=barplot(t(apply(q<thr,2:3,sum)),beside = T,col=rep(params$tissue.col[dimnames(q)[[2]]],each=7),border = NA,xaxt='n',...)
	text(b,rep(0,49),rep(substr(dimnames(q)[[3]],1,1),times=7),srt=0,adj=c(0.5,1),xpd=T,cex=0.7)
}

# pdf('figures/paper.figures/5/suppl/S10.main5EF.per.species.pdf',w=7.2,h=7.2,family='Arial')
# par(mfrow=c(3,2),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.3,2,1),oma=c(0,0,0,0))
# f(iuq,ylab='# of motifs',main='Inclusion upstream, fdr < 0.05')
# f(idq,ylab='# of motifs',main='Inclusion downstream, fdr < 0.05')
# f(euq,ylab='# of motifs',main='Exclusion upstream, fdr < 0.05')
# f(edq,ylab='# of motifs',main='Exclusion downstream, fdr < 0.05')
# f(euq,thr=0.5,ylab='# of motifs',main='Exclusion upstream, fdr < 0.5')
# f(edq,thr=0.5,ylab='# of motifs',main='Exclusion downstream, fdr < 0.5')
# dev.off()

f1 = function(t,s,thr=0.05,main=paste(s,t),only.sgn=F,filter=NULL,...){
	if(!is.null(s))
		qv = apply(cbind(hex.ups$up$pv[,t,s],hex.dws$up$pv[,t,s],hex.ups$dw$pv[,t,s],hex.dws$dw$pv[,t,s]),2,p.adjust,m='BH')
	else
		qv = cbind(hex.ups$up$ih.qv[,t],hex.dws$up$ih.qv[,t],hex.ups$dw$ih.qv[,t],hex.dws$dw$ih.qv[,t])
	qv = qv<thr
	if(!is.null(filter))
		qv = qv[filter,]
	if(only.sgn)
		qv = qv[apply(qv,1,sum)>0,]
	ov=calcAllPairsFT(qv)
	pv.col=c(gray=1,yellow=1e-2,orange=1e-5,red=1e-20,none=-1)
	or.col=c(blue=0,'#0000FFAA'=1/3,'#0000FF55'=1/1.5,gray=1,yellow=5,orange=10,red=Inf)
	dt=c('up\nupstream','up\ndownstream','down\nupstream','down\ndownstream')
	plotFTMotifSimMatrix(ov,T,pv.col=pv.col,or.col=or.col,diag.text=dt,main=main,...)
	invisible(qv)
}

# all.hex.qv = cbind(hex.ups$up$ih.qv, hex.dws$up$ih.qv, hex.ups$dw$ih.qv, hex.dws$dw$ih.qv)
# colnames(all.hex.qv) = paste0(substr(colnames(all.hex.qv),1,1),rep(c('iu','id','eu','ed'),each=7))
# all.hex.qv = all.hex.qv[,rep(c(1,8,15,22),times=7)+rep(0:6,each=4)]
# f = apply(all.hex.qv<0.05,1,sum)
# 
# hex.overlap=calcAllPairsFT(all.hex.qv<0.05)
# hex.overlap=calcAllPairsFT(all.hex.qv[apply(all.hex.qv<0.05,1,sum)>0,]<0.05)
# hex.overlap=calcAllPairsFT(all.hex.qv<0.05)

# pdf('figures/paper.figures/5/suppl/S10.hex.mirrored.positions.brain.heart.pdf',w=6,h=24,family='Arial')
# par(mfrow=c(8,2),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(1.5,1.3,1,1),oma=c(0,0,1,0))
# for(t in c(0.05,0.5)){
# 	f1('brain',NULL,t,main='All species, brain')
# 	f1('heart',NULL,t,main='All species, heart')
# 	for(s in rownames(species)){
# 		f1('brain',s,t)
# 		f1('heart',s,t)
# 	}
# 	mtext(paste0('FDR < ',t),3,outer=T)
# }
# dev.off()

fisher.test(table(b[,1],b[,4]))

f = apply(cbind(hex.ups$up$ih.qv,hex.dws$up$ih.qv,hex.ups$dw$ih.qv,hex.dws$dw$ih.qv)<0.05,1,sum)>0
pdf('figures/paper.figures/6/suppl/S10.sign.at.least.once.pdf',w=10,h=4,family='Arial')
#jpeg('figures/paper.figures/6/suppl/S10.sign.at.least.once.jpg',units = 'in',w=10,h=4,quality = 100,res = 600,family='Arial')
par(mfrow=c(1,2),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(1.5,1.3,1,1),oma=c(0,0,1,0))
b = f1('brain',NULL,0.05,main='Brain',filter = f)
h = f1('heart',NULL,0.05,main='Heart',filter = f)
dev.off()

# there are 1041 hexamers enriched at least once
# 

t = 'brain'
qv = cbind(hex.ups$up$ih.qv[,t],hex.dws$up$ih.qv[,t],hex.ups$dw$ih.qv[,t],hex.dws$dw$ih.qv[,t]) 
dim(qv)
apply(qv<0.05,2,sum)
f = apply(cbind(hex.ups$up$ih.qv,hex.dws$up$ih.qv,hex.ups$dw$ih.qv,hex.dws$dw$ih.qv)<0.05,1,sum)>0
fisher.test(table(qv[f,2]<0.05,qv[f,3]<0.05))
table(qv[,2]<0.5,qv[,3]<0.5)

t = 'heart'
pv = cbind(hex.ups$up$ih.pv[,t],hex.dws$up$ih.pv[,t],hex.ups$dw$ih.pv[,t],hex.dws$dw$ih.pv[,t]) 
imageWithText(cor(log(pv)))
pairs(qv,log='xy')

# _11 ######
# QKI per species (see fig 5 for data)

pdf('figures/paper.figures/6/suppl/S11.pdf',w=7.2/3*2,h=7.2/3*5,family='Arial')
#jpeg('figures/paper.figures/6/suppl/S11.jpg',units = 'in',w=7.2/3*2,h=7.2/3*5,quality = 100,res = 600,family='Arial')
par(mfrow=c(5,2),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.3,2,1),oma=c(0,0,0,0))
for(s in rownames(species)[-c(2,7)]){
	plotMirroredMotFreq(fa[s],age.segs[s],'actaac','brain','heart',main = paste0('ACTAAC in ',s,' brain devAS'),plot.leg = s=='human')
	plotMirroredMotFreq(fa[s],age.segs[s],'actaac','heart','brain',main = paste0('ACTAAC in heart ',s,' devAS'),plot.leg = F)
}
dev.off()


#  _12 APP ######
# APP ENSG00000142192
# save annotation for APP
hgmd = read.table('input/hgmd/tosend/2017_1_HGMD_ALL_Variants.csv',sep='\t',row.names=1,header=T,quote='',comment.char = '')
all.anns. = readRDS('Rdata/anns.subset.for.newborn.cov.plot.Rdata')
load('Rdata/tmp.exon.birth.Rdata')
exon.birth.one[sapply(exon.birth.one,function(x)x$gene_id[1] == 'hum.40661')] #t860 hum.40661.s11

hgid = 'ENSG00000142192'
# app.ann.h = all.anns$human[all.anns$human$gene_id=='hum.40661' & all.anns$human$sites %in% c('sd','ad','ae'),]
# saveRDS(app.ann.h,'Rdata/paper.figures/app.ann.h.Rdata')
app.ann.h = readRDS('Rdata/paper.figures/app.ann.h.Rdata')

#hgmd
s=anns$human['hum.40661.s11',]
hgmd[hgmd$chrom_VCF_hg19 == s$chr_id & hgmd$pos_VCF_hg19 >= s$start-200 & hgmd$pos_VCF_hg19 <= s$stop+200,]

# born from TE?
h2te[h2te$seg_id=='hum.40661.s11',] #no
# protein domain?
# o = get_OBO('processed/exon.onthology/exont.obo')
# seg2exont = readRDS('Rdata/seg2exont.Rdata')
# seg2exont[['hum.40661.s11']]
# o$name[seg2exont[['hum.40661.s11']]]
#"Protein binding" "Intrinsically_unstructured_polypeptide_region"
# check filters
nb.stat['t860',]

cov = readRDS('Rdata/newborn.exon.cov/hqmrb-hqmrb.t860.Rdata')

bams = paste0('processed/mapping/hisat2.s/human/',meta$fname[meta$species=='human'],'.bam')
# set.seed(1)
# bams = sample(bams)
# app.human.whole.cov=getReadCoverage(bams,'21',27252861,27543446,1)
# saveRDS(app.human.whole.cov,'Rdata/paper.figures/app.human.whole.cov.Rdata')
app.human.whole.cov = readRDS('Rdata/paper.figures/app.human.whole.cov.Rdata')

pdf('figures/paper.figures/6/suppl/S12.app.pdf',w=6,h=8,family='Arial')
#jpeg('figures/paper.figures/6/suppl/S12.app.jpg',units = 'in',w=6,h=8,quality = 100,res = 600,family='Arial')
par(tck=-0.01,mgp=c(1.1,0.2,0),mar=c(2,2.2,1,0),oma=c(1,0,0,1),xpd=F)
layout(matrix(c(1,1,1,1:29),ncol=4,byrow = T),widths = c(0.3,3,1,1))
i = 't860'
maxcov = max(app.human.whole.cov$cov,app.human.whole.cov$juncs$score)
plotReadCov(app.human.whole.cov,reverse = T,min.junc.cov = maxcov*0.05,bty='n',xlab='',ylab='',yaxt='n',xaxt='n',
						main=substitute(paste(italic('APP'), " - amyloid beta (A4) precursor protein")),ylim=c(-0.2*maxcov,maxcov),plot.junc.only.within=T,junc.lwd=2,junc.col='black')

segments(min(app.ann.h$start),-0.12*maxcov,max(app.ann.h$stop),-0.12*maxcov)

rect(app.ann.h$start,-0.2*maxcov,app.ann.h$stop,-0.04*maxcov,col = 'black',border = 'black')
seg = app.ann.h['hum.40661.s11',]
rect(seg$start,-0.2*maxcov,seg$stop,-0.04*maxcov,col = 'red',border = 'red')
xwhole = grconvertX(range(cov$human$x),'user','ndc')
ywhole = grconvertY(-0.23*maxcov,'user','ndc')

par(mar=c(1,2,1,0))
nb = exon.birth.one[[i]]
for(s in rownames(species)){
	plot.new()
	par(xpd=NA,mar=c(0,1,0,0))
	plotPNG(paste0('figures/paper.figures/5/icons/',s,'.png'),0.5,0.5,1.6)
	par(xpd=F,mar=c(1,1,1,0))
	maxcov = max(cov[[s]]$juncs$score[cov[[s]]$juncs$start >= min(cov[[s]]$x) & cov[[s]]$juncs$end <= max(cov[[s]]$x)])#,cov[[s]]$cov)
	plotReadCov(cov[[s]],reverse = (all.anns.[[s]][nb[s,'useg_id'],'strand'] == -1),yaxt='n',xaxt='n',min.junc.cov = maxcov*0.05,bty='n',xlab='',ylab='',
							main='',ylim=c(-0.2*maxcov,maxcov),plot.junc.only.within=T,junc.lwd=2,junc.col='black',xaxt='n')
	abline(h = -0.12*maxcov)
	seg = all.anns.[[s]][nb[s,'useg_id'],]
	rect(seg$start,-0.2*maxcov,seg$stop,-0.04*maxcov,col = 'black',border = 'black')
	seg = all.anns.[[s]][nb[s,'dseg_id'],]
	rect(seg$start,-0.2*maxcov,seg$stop,-0.04*maxcov,col = 'black',border = 'black')
	if(s=='human'){
		xwhole = grconvertX(xwhole,'ndc','user')
		ywhole = grconvertY(ywhole,'ndc','user')
		segments(xwhole[2],ywhole,grconvertX(0.08,'nfc','user'),grconvertY(1,'nfc','user'),xpd=NA,lty=2)
		segments(xwhole[1],ywhole,grconvertX(0.98,'nfc','user'),grconvertY(1,'nfc','user'),xpd=NA,lty=2)
	}
	
	if(!is.na(nb[s,'seg_id'])){
		rect(nb[s,'start'],-0.2*maxcov,nb[s,'stop'],-0.04*maxcov,col = 'red',border =  'red')
	}	
	par(mar=c(1,2.5,1,0))
	if(!is.na(nb[s,'seg_id'])){
		plotTissueAgeProile(psi.tsm[[s]][nb[s,'seg_id'],],meta.tsm,age.axis = 'rank',ylab='PSI',main=ifelse(s=='human','AS',''),bty='n',xlab='',age=meta.tsm$corr.age.rank,xlim=range(meta.tsm$corr.age.rank,na.rm=T),plot.xaxt=F,ylim=c(0,1))
		axis(1,labels = NA)
	}
	else
		plot.new()
	plotTissueAgeProile(ens.ge.marg.tsm[[s]][orth.ens.genes[hgid,s],],meta.tsm,age.axis = 'rank',ylab='RPKM',main=ifelse(s=='human','GE',''),bty='n',xlab='',age=meta.tsm$corr.age.rank,xlim=range(meta.tsm$corr.age.rank,na.rm=T),plot.xaxt=F)
	axis(1,labels = NA)
}
text(grconvertX(0,'nfc','user'),grconvertY(0,'ndc','user'),'Development',adj=c(0.5,-1),xpd=NA)
dev.off()

mou.je.as = readRDS('Rdata/japan.embrio/mou.je.as.filtered.Rdata')
meta.je = readRDS('Rdata/japan.embrio/meta.je.Rdata')
mou.je.as$e['mou.18737.s13',]
plotTissueAgeProile(psi.tsm$mouse['mou.18737.s13',],meta.tsm,age.axis = 'rank',ylab='PSI')

plot(meta.je[colnames(mou.je.as$ir),'days'],mou.je.as$ir['mou.18737.s13',])

# _14 6BC, micorexons ######
at = (0:19)[-seq(3,19,by = 3)]
cols = c(rep(params$tissue.col,each=2))
pch=c(rep(c(19,1),times=7))
xax = setNames(seq(0.5,by=3,length.out = 7),unique(meta$tissue))
dPSI = 0.2

pdf('figures/paper.figures/6/suppl/S14.0.6BC.microexons.pdf',w=5,h=12,family='Arial')
#jpeg('figures/paper.figures/6/suppl/S14.0.6BC.microexons.jpg',units = 'in',w=5,h=12,quality = 100,res = 600,family='Arial')
par(mfrow=c(7,2),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.3,2,1),oma=c(0,0,0,0),xpd=NA)
for(s in rownames(species)){
	mi.vs.ma.devAS = lapply(colnames(per.tissue.age.qv[[s]]),function(t){
		f = anns[[s]]$sites=='ad'
		sgn = per.tissue.age.qv[[s]][f,t] < 0.05 & abs(age.dpsi[[s]])[f,t] >  dPSI
		micro = anns[[s]]$length[f] <= 27
		t = table(devAS=sgn,micro=micro)
		if(nrow(t)==0)
			t = matrix(0,ncol=2,nrow=2,dimnames = list(c('FALSE','TRUE'),c('FALSE','TRUE')))
		t
	})
	
	in.vs.ex.micro = lapply(colnames(per.tissue.age.qv[[s]]),function(t){
		f = anns[[s]]$sites=='ad' & per.tissue.age.qv[[s]][,t] < 0.05 & abs(age.dpsi[[s]])[,t] >  dPSI
		micro = anns[[s]]$length[f] <= 27
		t=table(incl=age.segs[[s]][f,t] == 'u',micro=micro)
		if(nrow(t)==0)
			t = matrix(0,ncol=2,nrow=2,dimnames = list(c('FALSE','TRUE'),c('FALSE','TRUE')))
		t
	})
	names(mi.vs.ma.devAS) = names(in.vs.ex.micro) = colnames(per.tissue.age.qv[[s]])
	
	sapply(mi.vs.ma.devAS,function(x){ft=fisher.test(x);c(pv=ft$p.value,or=ft$estimate)})
	sapply(in.vs.ex.micro,function(x){ft=fisher.test(x);c(pv=ft$p.value,or=ft$estimate)})
	
	
	stat = matrix(NA,ncol=3,nrow=14)
	stat[0:6*2+1,] = t(sapply(mi.vs.ma.devAS,function(x)my.binom.test(x[c('TRUE','FALSE'),'TRUE'])))[,c(2,1,3)]
	stat[0:6*2+2,] = t(sapply(mi.vs.ma.devAS,function(x)my.binom.test(x[c('TRUE','FALSE'),'FALSE'])))[,c(2,1,3)]
	
	pv = sapply(mi.vs.ma.devAS,function(x){ft=fisher.test(x)$p.value})
	plotASSegStat(at,stat,pv,cols,xax,skip.for.pv=NULL,lty=c(rep(1:2,times=7)),main='',ylab='Proportion of devAS',pch=pch,bty='n',first2horiz=FALSE,xlim=c(0,21))
	if(s == 'human'){
		#legend(grconvertX(0.5,'npc','user'),grconvertY(0.95,'npc','user'),pch=c(19,1),col='black',legend=c('microexons','macroexons'),xjust = 0.5,yjust=0,bty='n',xpd=TRUE,ncol=2)
		text(grconvertX(0.85,'npc','user'),grconvertY(0.7,'npc','user'),c('p-value\n* <0.05\n** <0.01\n*** <0.001'),xpd=NA,adj=c(0,1),cex=0.7)
	}
	plotPNG(paste0("figures/paper.figures/5/icons/",s,".png"),1,0.9,0.14)
	if(s == 'human')
		legend(grconvertX(0.5,'npc','user'),grconvertY(0.95,'npc','user'),pch=c(19,1),col='black',legend=c('microexons','macroexons'),xjust = 0.5,yjust=0,bty='n',xpd=TRUE,ncol=2)
	
	stat = matrix(NA,ncol=3,nrow=14)
	stat[0:6*2+1,] = t(sapply(in.vs.ex.micro,function(x)my.binom.test(x[c('TRUE','FALSE'),'TRUE'])))[,c(2,1,3)]
	stat[0:6*2+2,] = t(sapply(in.vs.ex.micro,function(x)my.binom.test(x[c('TRUE','FALSE'),'FALSE'])))[,c(2,1,3)]
	
	pv = sapply(in.vs.ex.micro,function(x){ft=fisher.test(x)$p.value})
	plotASSegStat(at,stat,pv,cols,xax,skip.for.pv=NULL,lty=c(rep(1:2,times=7)),main='',ylab='Proportion of up',pch=pch,bty='n',first2horiz=FALSE)
	if(s == 'human'){
		legend(grconvertX(0.5,'npc','user'),grconvertY(0.95,'npc','user'),pch=c(19,1),col='black',legend=c('microexons','macroexons'),xjust = 0.5,yjust=0,bty='n',xpd=TRUE,ncol=2)
		#	text(grconvertX(0.85,'npc','user'),grconvertY(0.9,'npc','user'),c('p-value\n* <0.05\n** <0.01\n*** <0.001'),xpd=NA,adj=c(0,1),cex=0.7)
	}
}
dev.off()


# _S13 AMPD2 ######
# orth.seg.ad.all = readRDS('Rdata/orth.seg.ad.all.Rdata')
# for(s in names(orth.seg.ad.all))
# 	orth.seg.ad.all[[s]][c('i','e')] = NULL
# gc()

# alt.sp.byann = readRDS('Rdata/paper.figures/alt.sp.Rdata') # I changed alt.sp to new definition
# hgid = 'ENSG00000116337' # hum.3042.s5
# h = revList(seg2ens$human)
# alt.sp[h[[hgid]]]
# plotTissueAgeProile(psi.tsm$human['hum.3042.s5',],meta.tsm)

#hgmd
# s=anns$human['hum.3042.s5',]
# hgmd[hgmd$chrom_VCF_hg19 == s$chr_id & hgmd$pos_VCF_hg19 >= s$start-200 & hgmd$pos_VCF_hg19 <= s$stop+200,]
# # check orths
# orth.seg.ad.all$human$seg['hum.3042.s5',]
# idx=which(rownames(orth.seg.ad.all$human$seg) == 'hum.3042.s5')
# orth.sids = sapply(orth.seg.ad.all,function(x)rownames(x$seg)[idx])
# sapply(names(orth.sids),function(s)seg2ens[[s]][[orth.sids[s]]])
# in chicken my orth gene is on other unassembled contig than Ensembl ortholog
# AADN03018718.1 (my)  vs AADN03013667.1; both contigs are small. looks like the gene was fragmented in chicken assembly (AADN03018718.1 contains just 3 exons)
# so I'll not use chicken

# orth.ens.genes['ENSG00000116337',]
# ens.ge.cod$chicken$gene['ENSGALG00000026860',]
# orth.seg.ad.all$chicken$seg[orth.sids[7],]
# 
# orth.ens.genes['ENSG00000116337',-7] == sapply(names(orth.sids)[-7],function(s)seg2ens[[s]][[orth.sids[s]]])
# 
# 
# # protein domain?
# o = get_OBO('processed/exon.onthology/exont.obo')
# seg2exont = readRDS('Rdata/seg2exont.Rdata')
# seg2exont[['hum.3042.s5']]
# paste(o$name[seg2exont[['hum.3042.s5']]],collapse=', ')
# O-phospho-L-serine, Monomethylated L-arginine, Peroxisome targeting signal, Intrinsically_unstructured_polypeptide_region

# AMPD2.data = list()
# AMPD2.data$seg = all.anns$human[all.anns$human$gene_id=='hum.3042' & all.anns$human$sites %in% c('sd','ad','ae'),]
# bams = paste0('processed/mapping/hisat2.s/human/',meta$fname[meta$species=='human' & meta$tissue=='brain'],'.bam')
# AMPD2.data$human.whole.cov=getReadCoverage(bams,'1',min(AMPD2.data$seg$start),max(AMPD2.data$seg$stop),-1)
# orth.seg.ad.all.id = readRDS('Rdata/orth.seg.ad.all.id.Rdata')
# AMPD2.data$sids = orth.seg.ad.all.id['hum.3042.s5',-7]
# AMPD2.data$adj.segs = lapply(names(all.anns)[1:6],function(s){
# 	sid = AMPD2.data$sids[s]
# 	gid = all.anns[[s]][sid,'gene_id']
# 	en = all.anns[[s]][sid,'exon.number']
# 	all.anns[[s]][all.anns[[s]]$gene_id==gid & all.anns[[s]]$exon.number<=(en+1)& all.anns[[s]]$exon.number>=(en-2) & all.anns[[s]]$sites %in% c('sd','ad','ae'),]
# })
# names(AMPD2.data$adj.segs) = names(all.anns)[1:6]
# lapply(AMPD2.data$adj.segs,function(x)x[,c('sites','strand','length','exon.number')])
# 
# # load PSI and cov and exps
# AMPD2.data$psi = AMPD2.data$cov = AMPD2.data$ge = list()
# for(s in rownames(species)[-7]){
# 	print(s)
# 	psi = orth.seg.ad.all[[s]]$ir[orth.sids[s],colnames(orth.seg.ad.all[[s]]$ir) %in% rownames(meta),drop=F]
# 	m=meta[colnames(psi),]
# 	AMPD2.data$psi[[s]] = calcMeanCols(psi,paste(m$species,m$tissue,m$stage))
# 	bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue=='brain'],'.bam')
# 	AMPD2.data$cov[[s]]=getReadCoverage(bams,AMPD2.data$adj.segs[[s]]$chr_id[1],min(AMPD2.data$adj.segs[[s]]$start),max(AMPD2.data$adj.segs[[s]]$stop),-AMPD2.data$adj.segs[[s]]$strand[1])
# 	AMPD2.data$ge[[s]] = ens.ge.marg.tsm[[s]][orth.ens.genes['ENSG00000116337',s],]
# }
# 
# 
# saveRDS(AMPD2.data,'Rdata/paper.figures/AMPD2.data.Rdata')
AMPD2.data = readRDS('Rdata/paper.figures/AMPD2.data.Rdata')

pdf('figures/paper.figures/6/suppl/S13.ampd2.pdf',w=6,h=7,family='Arial')
#jpeg('figures/paper.figures/6/suppl/S13.ampd2.jpg',units = 'in',w=6,h=7,quality = 100,res = 600,family='Arial')
par(tck=-0.01,mgp=c(1.1,0.2,0),mar=c(2,2.2,1,0),oma=c(1,0,0,1),xpd=F)
layout(matrix(c(1,1,1,1:25),ncol=4,byrow = T),widths = c(0.3,3,1,1))

maxcov = max(AMPD2.data$human.whole.cov$cov,AMPD2.data$human.whole.cov$juncs$score)
plotReadCov(AMPD2.data$human.whole.cov,reverse = F,min.junc.cov = maxcov*0.01,bty='n',xlab='',ylab='',yaxt='n',xaxt='n',
						main=substitute(paste(italic('AMPD2'), " - adenosine monophosphate deaminase 2")),ylim=c(-0.2*maxcov,maxcov),plot.junc.only.within=T,junc.lwd=2,junc.col='black')

segments(min(AMPD2.data$seg$start),-0.12*maxcov,max(AMPD2.data$seg$stop),-0.12*maxcov)

rect(AMPD2.data$seg$start,-0.2*maxcov,AMPD2.data$seg$stop,-0.04*maxcov,col = 'black',border = 'black')
seg = AMPD2.data$seg['hum.3042.s5',]
rect(seg$start,-0.2*maxcov,seg$stop,-0.04*maxcov,col = 'red',border = 'red')
xwhole = grconvertX(range(AMPD2.data$cov$human$x),'user','ndc')
ywhole = grconvertY(-0.23*maxcov,'user','ndc')

par(mar=c(1,2,1,0))

for(s in rownames(species)[-7]){
	plot.new()
	par(xpd=NA,mar=c(0,1,0,0))
	plotPNG(paste0('figures/paper.figures/5/icons/',s,'.png'),0.5,0.5,1.6)
	par(xpd=F,mar=c(1,1,1,0))

	cov = AMPD2.data$cov
	seg. = AMPD2.data$adj.segs[[s]]
	maxcov = max(cov[[s]]$juncs$score[cov[[s]]$juncs$start >= min(cov[[s]]$x) & cov[[s]]$juncs$end <= max(cov[[s]]$x)])#,cov[[s]]$cov)
	plotReadCov(cov[[s]],reverse = seg.$strand[1] == -1,yaxt='n',xaxt='n',min.junc.cov = maxcov*0.05,bty='n',xlab='',ylab='',
							main='',ylim=c(-0.2*maxcov,maxcov),plot.junc.only.within=T,junc.lwd=2,junc.col='black',xaxt='n')
	segments(min(seg.$start),-0.12*maxcov,max(seg.$stop),-0.12*maxcov)
	
	rect(seg.$start,-0.2*maxcov,seg.$stop,-0.04*maxcov,col = 'black',border = 'black')
	seg = seg.[AMPD2.data$sids[s],]
	rect(seg$start,-0.2*maxcov,seg$stop,-0.04*maxcov,col = 'red',border = 'red')
	
	if(s=='human'){
		xwhole = grconvertX(xwhole,'ndc','user')
		ywhole = grconvertY(ywhole,'ndc','user')
		segments(xwhole[1],ywhole,grconvertX(0.08,'nfc','user'),grconvertY(1,'nfc','user'),xpd=NA,lty=2)
		segments(xwhole[2],ywhole,grconvertX(0.98,'nfc','user'),grconvertY(1,'nfc','user'),xpd=NA,lty=2)
	}
	
	par(mar=c(1,2.5,1,0))
	plotTissueAgeProile(AMPD2.data$psi[[s]][1,],meta.tsm,age.axis = 'rank',ylab='PSI',main=ifelse(s=='human','AS',''),bty='n',xlab='',age=meta.tsm$corr.age.rank,xlim=range(meta.tsm$corr.age.rank,na.rm=T),plot.xaxt=F,ylim=c(0,1))
	axis(1,labels = NA)

	plotTissueAgeProile(AMPD2.data$ge[[s]],meta.tsm,age.axis = 'rank',ylab='RPKM',main=ifelse(s=='human','GE',''),bty='n',xlab='',age=meta.tsm$corr.age.rank,xlim=range(meta.tsm$corr.age.rank,na.rm=T),plot.xaxt=F)
	axis(1,labels = NA)
}
text(grconvertX(0,'nfc','user'),grconvertY(0,'ndc','user'),'Development',adj=c(0.5,-1),xpd=NA)
dev.off()

# _S15 GDPD5 microexon brain ####
human.tau = read.csv('input/gene.info.from.marg/Human.Indexes.All.csv')
rownames(human.tau) = human.tau$Human_ID
s = which(anns$human$sites=='ad' & anns$human$length<=27 & per.tissue.age.qv$human[,'brain'] < 0.05 & age.dpsi$human[,'brain'] > 0.5)
s = data.frame(seg.id = rownames(anns$human)[s],length=anns$human$length[s],cod=anns$human$cod[s])
s$ens_id = sapply(s$seg.id,function(s)paste(seg2ens$human[[s]],collapse=','))
s$tau = sapply(s$seg.id,function(s)max(human.tau[seg2ens$human[[s]],'TissueTau'],na.rm=T))
s$orth7 = s$seg.id %in% rownames(orth.seg.ad.tsm$human)
s$descr = sapply(s$seg.id,function(s)paste(gene.descrs$human[seg2ens$human[[s]],'descr'],collapse=', '))
s$gene.name = sapply(s$seg.id,function(s)paste(gene.descrs$human[seg2ens$human[[s]],'gene.name'],collapse=', '))

s[s$orth7 & s$length<16 & s$tau<0.5,]

plotTissueAgeProile(psi.tsm$human['hum.11477.s7',],meta.tsm,age.axis = 'rank')
plotTissueAgeProile(ens.ge.cod$human['hum.11477',],meta.tsm)
plotTissueAgeProile(my.ge$human$rpkm['hum.11477',],meta,age.axis = 'rank')
plotTissueAgeProile( ens.ge.marg.tsm$human['ENSG00000158555',],meta.tsm,age.axis = 'rank')


gene.descrs$human[seg2ens$human[['hum.11477.s7']],]
ens.ge$human$gene['ENSG00000158555',]

anns$human['hum.11477.s7',]

# gdpd5.data = list()
# gdpd5.data$sid = 'hum.11477.s7'
# gdpd5.data$ens_id = 'ENSG00000158555'
# s = 'human'
# gdpd5.data$species = s
# gdpd5.data$segs = all.anns[[s]][all.anns[[s]]$gene_id == all.anns[[s]][gdpd5.data$sid,'gene_id'],]
# 
# bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s],'.bam')
# gdpd5.data$all.cov  = getReadCoverage(bams,gdpd5.data$segs$chr_id[1],min(gdpd5.data$segs$start),max(gdpd5.data$segs$stop),-gdpd5.data$segs$strand[1])
# 
# 
# gdpd5.data$tissue.stage.cov = list()
# f = function(t,d1,d2) getReadCoverage(paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue==t & meta$days >=d1 & meta$days < d2],'.bam'),gdpd5.data$segs$chr_id[1],75147981,75151006,-gdpd5.data$segs$strand[1])
# gdpd5.data$tissue.stage.cov$l1 = f('liver',-1,200)
# gdpd5.data$tissue.stage.cov$l2 = f('liver',280,1700)
# gdpd5.data$tissue.stage.cov$l3 = f('liver',5000,15000)
# 
# gdpd5.data$tissue.stage.cov$b1 = f('brain',-1,200)
# gdpd5.data$tissue.stage.cov$b2 = f('brain',280,1700)
# gdpd5.data$tissue.stage.cov$b3 = f('brain',5000,15000)
# 
# saveRDS(gdpd5.data,'Rdata/paper.figures/gdpd5.data.Rdata')
gdpd5.data = readRDS('Rdata/paper.figures/gdpd5.data.Rdata')

# look for orthologs:
oinx = orth.ads.all.sp[!is.na(orth.ads.all.sp[,'human']) & orth.ads.all.sp[,'human'] == gdpd5.data$sid,]


l = rbind(c(1, 1, 1, 1, 1, 1, 1), # whole cov
					c(2, 3, 3, 3, 4, 4, 4), # young cov
					c(5, 6, 6, 6, 7, 7, 7), # adult cov
					c(8, 9, 9, 9,10,10,10), # old cov
					matrix(11:31,ncol=7)
)
# __plot ######
pdf('figures/paper.figures/6/suppl/S15.GDPD5.brain.microexon.pdf',w=6,h=6,family='Arial')
#png('figures/paper.figures/6/suppl/S15.GDPD5.brain.microexon.png',units = 'in',w=6,h=6,res = 600,family='Arial')
exon.center = c(gdpd5.data$segs[gdpd5.data$sid,'start']/2+gdpd5.data$segs[gdpd5.data$sid,'stop']/2)
layout(l,heights = c(1.2,0.8,0.8,0.8,0.6,1,1))
par(tck=-0.01,mgp=c(1.1,0.2,0),mar=c(0,0,1,0),oma=c(0,1,0,1))
maxcov = max(gdpd5.data$all.cov$cov,gdpd5.data$all.cov$juncs$score)
plotReadCov(gdpd5.data$all.cov,reverse = T,min.junc.cov = maxcov*0.05,bty='n',xlab='',ylab='',yaxt='n',xaxt='n',
						main=substitute(paste(italic('GDPD5'), " - glycerophosphodiester phosphodiesterase domain containing 5")),ylim=c(-0.2*maxcov,maxcov),plot.junc.only.within=T,junc.lwd=2,junc.col='black')
par(xpd=NA)
plotPNG(paste0('figures/paper.figures/5/icons/human.png'),0.06,0.8,0.08)
par(xpd=T)
segments(min(gdpd5.data$segs$start),-0.12*maxcov,max(gdpd5.data$segs$stop),-0.12*maxcov)
f = !(gdpd5.data$segs$sites %in% c('da','dd','aa'))
rect(gdpd5.data$segs$start[f],-0.2*maxcov,gdpd5.data$segs$stop[f],-0.04*maxcov,col = 'black',border = 'black')
seg = gdpd5.data$segs[gdpd5.data$sid,]
rect(seg$start,-0.2*maxcov,seg$stop,-0.04*maxcov,col = 'red',border = 'red')

# reg.coor = grconvertX(c(papss2.data$ann.m['mou.35838.s2','start'],papss2.data$ann.m['mou.35838.s4','stop']),'user','ndc')
reg.coor = grconvertX(exon.center,'user','ndc')
reg.coory = grconvertY(0,'npc','ndc')
# plot age.cov
xcoor = range(gdpd5.data$tissue.stage.cov$l1$x)
age.text = c('prenatal','young','adult')
for(i in 1:3){
	plot.new()
	text(grconvertX(0.5,'npc','user'),grconvertY(0.3,'npc','user'),age.text[i],xpd=NA,adj=c(0.5,0.5),font=2,cex=1.2)
	plotReadCov(gdpd5.data$tissue.stage.cov[[paste0('l',i)]],reverse = T,min.junc.cov = 5,bty='n',xlab='',ylab='',yaxt='n',xaxt='n',
							main='',plot.junc.only.within=T,junc.lwd=2,junc.col='black',xlim=xcoor)
	if(i==1){
		mtext('Liver',3,-1,font=2,cex=1)
		reg.coor.  = grconvertX(reg.coor,'ndc','user')
		reg.coory.  = grconvertY(reg.coory,'ndc','user')
		y = grconvertY(1,'npc','user')
		arrows(reg.coor.,reg.coory.,exon.center,y,xpd=NA,length=0.1,col='red')
		# segments(xcoor,c(y,y),reg.coor.,c(reg.coory.,reg.coory.),xpd=NA,lty=2)
	}
	plotReadCov(gdpd5.data$tissue.stage.cov[[paste0('b',i)]],reverse = F,min.junc.cov = 5,bty='n',xlab='',ylab='',yaxt='n',xaxt='n',
							main='',plot.junc.only.within=T,junc.lwd=2,junc.col='black',xlim=xcoor)
	if(i==1){
		mtext('Brain',3,-1,font=2,cex=1)
		reg.coor.  = grconvertX(reg.coor,'ndc','user')
		reg.coory.  = grconvertY(reg.coory,'ndc','user')
		y = grconvertY(1,'npc','user')
		arrows(reg.coor.,reg.coory.,exon.center,y,xpd=NA,length=0.1,col='red')
		# segments(xcoor,c(y,y),reg.coor.,c(reg.coory.,reg.coory.),xpd=NA,lty=2)
	}
}

# plot PSI/RPKM
ogids = orth.ens.genes[orth.ens.genes[,'human'] == gdpd5.data$ens_id,]
#t(sapply(rownames(species),function(s)gene.descrs[[s]][ogids[1,s],]))
#oinx = which(rownames(orth.seg.ad.tsm$mouse) == 'mou.23510.s9')

par(mar=c(1.5,1.5,1,0),xpd=NA)
for(s in rownames(species)){
	plot.new()
	plotPNG(paste0('figures/paper.figures/5/icons/',s,'.png'),0.5,-0.5,0.7)
	plotTissueAgeProile(ens.ge.marg.tsm[[s]][ogids[1,s],],meta.tsm,age.axis = 'rank',ylab=ifelse(s=='human','RPKM',''),main='',bty='n',xlab='',age=meta.tsm$corr.age.rank,xlim=range(meta.tsm$corr.age.rank,na.rm=T),plot.xaxt=F,cex=1)
	axis(1,labels = NA)
	if(is.na(oinx[s]))
		plot.new()
	else
		plotTissueAgeProile(psi.tsm[[s]][oinx[s],],meta.tsm,age.axis = 'rank',ylab=ifelse(s=='human','PSI',''),main='',bty='n',xlab='',age=meta.tsm$corr.age.rank,xlim=range(meta.tsm$corr.age.rank,na.rm=T),plot.xaxt=F,ylim=c(0,1),cex=1)
	axis(1,labels = NA)
}
dev.off()


# _S16 Tmed2 microexon heart ####
# gene.descrs$mouse[seg2ens$mouse[['mou.35838.s3']],]
# ens.ge$mouse$gene['ENSMUSG00000029390',]
# anns$mouse['mou.35838.s3',]

# tmed2.data = list()
# s = 'mouse'
# bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s],'.bam')
# tmed2.data$all.cov  = getReadCoverage(bams,'5',124540772,124607430,-1)
# 
# tmed2.data$tissue.stage.cov = list()
# f = function(t,d1,d2) getReadCoverage(paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue==t & meta$days >=d1 & meta$days < d2],'.bam'),'5',124540695,124550506,-1)
# tmed2.data$tissue.stage.cov$h1 = f('heart',-1,15)
# tmed2.data$tissue.stage.cov$h2 = f('heart',20,45)
# tmed2.data$tissue.stage.cov$h3 = f('heart',45,1000)
# 
# tmed2.data$tissue.stage.cov$b1 = f('brain',-1,15)
# tmed2.data$tissue.stage.cov$b2 = f('brain',20,45)
# tmed2.data$tissue.stage.cov$b3 = f('brain',45,1000)
# 
# tmed2.data$ann.m = all.anns$mouse[all.anns$mouse$gene_id=='mou.35838' & !(all.anns$mouse$sites %in% c('da','dd','aa')),]
# saveRDS(tmed2.data,'Rdata/paper.figures/tmed2.data.Rdata')
tmed2.data = readRDS('Rdata/paper.figures/tmed2.data.Rdata')
l = rbind(c(1, 1, 1, 1, 1, 1, 1), # whole cov
					c(2, 3, 3, 3, 4, 4, 4), # young cov
					c(5, 6, 6, 6, 7, 7, 7), # adult cov
					c(8, 9, 9, 9,10,10,10), # old cov
					matrix(11:31,ncol=7)
)
# __ plot ######
pdf('figures/paper.figures/6/suppl/S16.Tmed2.heart.microexon.pdf',w=6,h=6,family='Arial')
#png('figures/paper.figures/6/suppl/S16.Tmed2.heart.microexon.png',units = 'in',w=6,h=6,res = 600,family='Arial')
exon.center = c(tmed2.data$ann.m['mou.35838.s3','start']/2+tmed2.data$ann.m['mou.35838.s3','stop']/2)
layout(l,heights = c(1.2,0.8,0.8,0.8,0.6,1,1))
par(tck=-0.01,mgp=c(1.1,0.2,0),mar=c(0,0,1,0),oma=c(0,1,0,1))
maxcov = max(tmed2.data$all.cov$cov,tmed2.data$all.cov$juncs$score)
plotReadCov(tmed2.data$all.cov,reverse = F,min.junc.cov = maxcov*0.05,bty='n',xlab='',ylab='',yaxt='n',xaxt='n',
						main=substitute(paste(italic('Tmed2'), " - transmembrane p24 trafficking protein 2")),ylim=c(-0.2*maxcov,maxcov),plot.junc.only.within=T,junc.lwd=2,junc.col='black',xlim=c(124540772,124550393))
plotPNG(paste0('figures/paper.figures/5/icons/mouse.png'),0.06,0.8,0.08)
segments(min(tmed2.data$ann.m$start),-0.12*maxcov,124550393,-0.12*maxcov)

rect(tmed2.data$ann.m$start,-0.2*maxcov,tmed2.data$ann.m$stop,-0.04*maxcov,col = 'black',border = 'black')
seg = tmed2.data$ann.m['mou.35838.s3',]
rect(seg$start,-0.2*maxcov,seg$stop,-0.04*maxcov,col = 'red',border = 'red')

# reg.coor = grconvertX(c(tmed2.data$ann.m['mou.35838.s2','start'],tmed2.data$ann.m['mou.35838.s4','stop']),'user','ndc')
reg.coor = grconvertX(exon.center,'user','ndc')
reg.coory = grconvertY(0,'npc','ndc')
# plot age.cov
xcoor = c(124543000,124547000)
age.text = c('prenatal','young','adult')
for(i in 1:3){
	plot.new()
	text(grconvertX(0.5,'npc','user'),grconvertY(0.3,'npc','user'),age.text[i],xpd=NA,adj=c(0.5,0.5),font=2,cex=1.2)
	plotReadCov(tmed2.data$tissue.stage.cov[[paste0('h',i)]],reverse = F,min.junc.cov = 5,bty='n',xlab='',ylab='',yaxt='n',xaxt='n',
							main='',plot.junc.only.within=T,junc.lwd=2,junc.col='black',xlim=xcoor)
	if(i==1){
		mtext('Heart',3,-1,font=2,cex=1)
		reg.coor.  = grconvertX(reg.coor,'ndc','user')
		reg.coory.  = grconvertY(reg.coory,'ndc','user')
		y = grconvertY(0.9,'npc','user')
		arrows(reg.coor.,reg.coory.,exon.center,y,xpd=NA,length=0.1,col='red')
		# segments(xcoor,c(y,y),reg.coor.,c(reg.coory.,reg.coory.),xpd=NA,lty=2)
	}
	plotReadCov(tmed2.data$tissue.stage.cov[[paste0('b',i)]],reverse = F,min.junc.cov = 5,bty='n',xlab='',ylab='',yaxt='n',xaxt='n',
							main='',plot.junc.only.within=T,junc.lwd=2,junc.col='black',xlim=xcoor)
	if(i==1){
		mtext('Brain',3,-1,font=2,cex=1)
		reg.coor.  = grconvertX(reg.coor,'ndc','user')
		reg.coory.  = grconvertY(reg.coory,'ndc','user')
		y = grconvertY(0.9,'npc','user')
		arrows(reg.coor.,reg.coory.,exon.center,y,xpd=NA,length=0.1,col='red')
		# segments(xcoor,c(y,y),reg.coor.,c(reg.coory.,reg.coory.),xpd=NA,lty=2)
	}
}

# plot PSI/RPKM
ogids = orth.ens.genes[orth.ens.genes[,'mouse'] == 'ENSMUSG00000029390',]
oinx = which(rownames(orth.seg.ad.tsm$mouse) == 'mou.35838.s3')
par(mar=c(1.5,1.5,1,0),xpd=NA)
for(s in rownames(species)){
	plot.new()
	plotPNG(paste0('figures/paper.figures/5/icons/',s,'.png'),0.5,-0.5,0.7)
	plotTissueAgeProile(ens.ge.marg.tsm[[s]][ogids[1,s],],meta.tsm,age.axis = 'rank',ylab=ifelse(s=='human','RPKM',''),main='',bty='n',xlab='',age=meta.tsm$corr.age.rank,xlim=range(meta.tsm$corr.age.rank,na.rm=T),plot.xaxt=F,cex=1)
	axis(1,labels = NA)
	plotTissueAgeProile(orth.seg.ad.tsm[[s]][oinx,],meta.tsm,age.axis = 'rank',ylab=ifelse(s=='human','PSI',''),main='',bty='n',xlab='',age=meta.tsm$corr.age.rank,xlim=range(meta.tsm$corr.age.rank,na.rm=T),plot.xaxt=F,ylim=c(0,1),cex=1)
	axis(1,labels = NA)
}
dev.off()


# _S17 Papss2 microexon liver ####
# gene.descrs$mouse[seg2ens$mouse[['mou.23510.s9']],]
# ens.ge$mouse$gene['ENSMUSG00000024899',]
# anns$mouse['mou.23510.s9',]

#all.anns$mouse[all.anns$mouse$gene=='mou.23510',]

# papss2.data = list()
# s = 'mouse'
# bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s],'.bam')
# papss2.data$all.cov  = getReadCoverage(bams,'19',32620005,32667187,-1)
# 
# papss2.data$tissue.stage.cov = list()
# f = function(t,d1,d2) getReadCoverage(paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue==t & meta$days >=d1 & meta$days < d2],'.bam'),'19',32620005,32667187,-1)
# papss2.data$tissue.stage.cov$l1 = f('liver',-1,15)
# papss2.data$tissue.stage.cov$l2 = f('liver',20,45)
# papss2.data$tissue.stage.cov$l3 = f('liver',45,1000)
# 
# papss2.data$tissue.stage.cov$b1 = f('brain',-1,15)
# papss2.data$tissue.stage.cov$b2 = f('brain',20,45)
# papss2.data$tissue.stage.cov$b3 = f('brain',45,1000)
# 
# papss2.data$ann.m = all.anns$mouse[all.anns$mouse$gene_id=='mou.23510' & all.anns$mouse$sites %in% c('sd','ad','ae'),]
# saveRDS(papss2.data,'Rdata/paper.figures/papss2.data.Rdata')
papss2.data = readRDS('Rdata/paper.figures/papss2.data.Rdata')

# look for orthologs:
oinx = orth.ads.all.sp[!is.na(orth.ads.all.sp[,'mouse']) & orth.ads.all.sp[,'mouse'] == 'mou.23510.s9',]
b1 = orth.ads.all.sp[!is.na(orth.ads.all.sp[,'mouse']) & orth.ads.all.sp[,'mouse'] == 'mou.23510.s8',]
b2 = orth.ads.all.sp[!is.na(orth.ads.all.sp[,'mouse']) & orth.ads.all.sp[,'mouse'] == 'mou.23510.s10',]

# lapply(rownames(species),function(s)
# 	all.anns[[s]][all.anns[[s]]$sites=='ad' & all.anns[[s]]$gene_id == all.anns[[s]][b1[s],'gene_id'] & all.anns[[s]]$start>=min(all.anns[[s]][c(b1[s],b2[s]),'start']) & all.anns[[s]]$stop<=max(all.anns[[s]][c(b1[s],b2[s]),'stop']),]
# )
# add manually
oinx['human'] = 'hum.8123.s11'
oinx['opossum'] = 'opo.732.s8' 
oinx['chicken'] = 'chi.25557.s8'  # it is 21nt long instead of 15

l = rbind(c(1, 1, 1, 1, 1, 1, 1), # whole cov
					c(2, 3, 3, 3, 4, 4, 4), # young cov
					c(5, 6, 6, 6, 7, 7, 7), # adult cov
					c(8, 9, 9, 9,10,10,10), # old cov
					matrix(11:31,ncol=7)
)
# __plot ######
pdf('figures/paper.figures/6/suppl/S17.Papss2.liver.microexon.pdf',w=6,h=6,family='Arial')
#png('figures/paper.figures/6/suppl/S17.Papss2.liver.microexon.png',units = 'in',w=6,h=6,res = 600,family='Arial')
exon.center = c(papss2.data$ann.m['mou.23510.s9','start']/2+papss2.data$ann.m['mou.23510.s9','stop']/2)
layout(l,heights = c(1.2,0.8,0.8,0.8,0.6,1,1))
par(tck=-0.01,mgp=c(1.1,0.2,0),mar=c(0,0,1,0),oma=c(0,1,0,1))
maxcov = max(papss2.data$all.cov$cov,papss2.data$all.cov$juncs$score)
plotReadCov(papss2.data$all.cov,reverse = F,min.junc.cov = maxcov*0.05,bty='n',xlab='',ylab='',yaxt='n',xaxt='n',
						main=substitute(paste(italic('Papss2'), " - 3'-phosphoadenosine 5'-phosphosulfate synthase 2")),ylim=c(-0.2*maxcov,maxcov),plot.junc.only.within=T,junc.lwd=2,junc.col='black')
plotPNG(paste0('figures/paper.figures/5/icons/mouse.png'),0.06,0.8,0.08)
segments(32619961,-0.12*maxcov,max(papss2.data$ann.m$stop),-0.12*maxcov)

rect(papss2.data$ann.m$start,-0.2*maxcov,papss2.data$ann.m$stop,-0.04*maxcov,col = 'black',border = 'black')
seg = papss2.data$ann.m['mou.23510.s9',]
rect(seg$start,-0.2*maxcov,seg$stop,-0.04*maxcov,col = 'red',border = 'red')

# reg.coor = grconvertX(c(papss2.data$ann.m['mou.35838.s2','start'],papss2.data$ann.m['mou.35838.s4','stop']),'user','ndc')
reg.coor = grconvertX(exon.center,'user','ndc')
reg.coory = grconvertY(0,'npc','ndc')
# plot age.cov
xcoor = c(32641300,32652101)
age.text = c('prenatal','young','adult')
for(i in 1:3){
	plot.new()
	text(grconvertX(0.5,'npc','user'),grconvertY(0.3,'npc','user'),age.text[i],xpd=NA,adj=c(0.5,0.5),font=2,cex=1.2)
	plotReadCov(papss2.data$tissue.stage.cov[[paste0('l',i)]],reverse = F,min.junc.cov = 5,bty='n',xlab='',ylab='',yaxt='n',xaxt='n',
							main='',plot.junc.only.within=T,junc.lwd=2,junc.col='black',xlim=xcoor)
	if(i==1){
		mtext('Liver',3,-1,font=2,cex=1)
		reg.coor.  = grconvertX(reg.coor,'ndc','user')
		reg.coory.  = grconvertY(reg.coory,'ndc','user')
		y = grconvertY(0.6,'npc','user')
		arrows(reg.coor.,reg.coory.,exon.center,y,xpd=NA,length=0.1,col='red')
		# segments(xcoor,c(y,y),reg.coor.,c(reg.coory.,reg.coory.),xpd=NA,lty=2)
	}
	plotReadCov(papss2.data$tissue.stage.cov[[paste0('b',i)]],reverse = F,min.junc.cov = 5,bty='n',xlab='',ylab='',yaxt='n',xaxt='n',
							main='',plot.junc.only.within=T,junc.lwd=2,junc.col='black',xlim=xcoor)
	if(i==1){
		mtext('Brain',3,-1,font=2,cex=1)
		reg.coor.  = grconvertX(reg.coor,'ndc','user')
		reg.coory.  = grconvertY(reg.coory,'ndc','user')
		y = grconvertY(0.6,'npc','user')
		arrows(reg.coor.,reg.coory.,exon.center,y,xpd=NA,length=0.1,col='red')
		# segments(xcoor,c(y,y),reg.coor.,c(reg.coory.,reg.coory.),xpd=NA,lty=2)
	}
}

# plot PSI/RPKM
ogids = orth.ens.genes[orth.ens.genes[,'mouse'] == 'ENSMUSG00000024899',]
#oinx = which(rownames(orth.seg.ad.tsm$mouse) == 'mou.23510.s9')
orth.ads.all.sp[!is.na(orth.ads.all.sp[,'mouse']) & orth.ads.all.sp[,'mouse'] == 'mou.23510.s9',]
orth.ads.all.sp[!is.na(orth.ads.all.sp[,'mouse']) & orth.ads.all.sp[,'mouse'] == 'mou.23510.s8',]
par(mar=c(1.5,1.5,1,0),xpd=NA)
for(s in rownames(species)){
	plot.new()
	plotPNG(paste0('figures/paper.figures/5/icons/',s,'.png'),0.5,-0.4,0.7)
	plotTissueAgeProile(ens.ge.marg.tsm[[s]][ogids[1,s],],meta.tsm,age.axis = 'rank',ylab=ifelse(s=='human','RPKM',''),main='',bty='n',xlab='',age=meta.tsm$corr.age.rank,xlim=range(meta.tsm$corr.age.rank,na.rm=T),plot.xaxt=F,cex=1)
	axis(1,labels = NA)
	if(is.na(oinx[s]))
		plot.new()
	else
		plotTissueAgeProile(psi.tsm[[s]][oinx[s],],meta.tsm,age.axis = 'rank',ylab=ifelse(s=='human','PSI',''),main='',bty='n',xlab='',age=meta.tsm$corr.age.rank,xlim=range(meta.tsm$corr.age.rank,na.rm=T),plot.xaxt=F,ylim=c(0,1),cex=1)
	axis(1,labels = NA)
}
dev.off()



# _S10 dPSI/patterns def ######
# prev 18
which(anns$mouse$sites == 'ad' & age.segs$mouse[,'brain'] == 'ud' & age.dpsi$mouse[,'brain'] > 0.5)
plotTissueAgeProile(psi.tsm$mouse['mou.33599.s33',],meta.tsm,age.axis = 'log') #mou.45495.s3

devAS.change.stat = readRDS('Rdata/devAS.change.stat.Rdata')
s = 'mouse'
t = 'brain'
sid = 'mou.33599.s33'
psi.example = readRDS(paste0("Rdata/",s,".as.u.filtered.Rdata"))$ir[sid,]
meta.example = meta[names(psi.example),]
meta.example = meta.example[meta.example$tissue == t,]
meta.example = meta.example[order(meta.example$days),]
psi.example = psi.example[rownames(meta.example)]
spline.fit = smooth.spline(meta.example$age.use,psi.example,df=4)

age1000 = seq(min(meta.example$age.use),max(meta.example$age.use),length.out=1000)
psi.pred1000 = predict(spline.fit,age1000)$y
psi.pred = predict(spline.fit)

pdf('figures/paper.figures/6/suppl/new.order/S10a.dPSI.definition.pdf',w=5,h=6,family='Arial')
#jpeg('figures/paper.figures/6/suppl/S17.dPSI.definition.jpg',units = 'in',w=6,h=6,quality = 100,res = 600,family='Arial')
par(mfrow=c(3,2),tck=-0.01,mgp=c(1.4,0.4,0),mar=c(2.5,2.5,1,0),oma=c(0,0,0,1))
plotS18A()
plotPanelLetter('A',lab.cex)
plotS18B()
plotPanelLetter('B',lab.cex)
plotS18C()
plotPanelLetter('C',lab.cex)
plotS18D()
plotPanelLetter('D',lab.cex)
plot.new()
plotS18E()
plotPanelLetter('E',lab.cex)
dev.off()

# _s11 4patterns mouse details #####
# prev 19
devAS.change.stat = readRDS('Rdata/devAS.change.stat.Rdata')

sites ='ad'
s = 'mouse'
thr = 0.3
pdf('figures/paper.figures/6/suppl/new.order/S11a.mouse.4patterns.pdf',w=10,h=10,family='Arial')
#png('figures/paper.figures/6/suppl/S19.mouse.4patterns.png',units = 'in',w=10,h=10,res = 400,family='Arial')
par(mfrow=c(7,7),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,1,1))

for(t in unique(meta$tissue)){
	y = devAS.change.stat[[s]][[t]]
	y = y[anns[[s]][rownames(y),'sites']==sites,]
	u = y[,1]/(y[,1]+y[,2])
	tim = y[,3]/y[,1]-y[,4]/y[,2]
	tim[is.na(tim)] = 0
	plot(u,tim,col=params$tissue.col[t],xlab='up/(up+down)',ylab='up_time - down_time',main=t)
	if(t == 'brain') plotPanelLetter('A',lab.cex)
	segments(thr,0,1-thr,0,lty=2)
	abline(v=c(thr,1-thr),lty=2)
	text(c(thr/2,1-thr/2,0.5,0.5),c(0,0,0.7,-0.7),c('up','down','up-down','down-up'),adj=0.5,cex=0.8)
	
	h=hist(u,0:20/20,xlab='up/(up+down)',col=params$tissue.col[t],border=NA,main='')
	abline(v=c(thr,1-thr),lty=2)
	text(c(thr/2,1-thr/2,0.5),rep(max(h$counts)/2,3),c('up','down','up-down\nor\ndown-up'),adj=0.5,cex=0.8)
	if(t == 'brain') plotPanelLetter('B',lab.cex)
	
	f = u > thr & u < (1-thr)
	h=hist(tim[f],20,xlab='up_time - down_time',col=params$tissue.col[t],border=NA,main='')
	abline(v=0,lty=2)
	xx = range(tim[f])
	text(c(xx[1]/2,xx[2]/2),rep(max(h$counts)/2,2),c('up-down','down-up'),adj=0.5,cex=0.8)
	if(t == 'brain') plotPanelLetter('C',lab.cex)
	
	mm = meta.tsm
	col = params$tissue.col
	col[names(col)!=t] = paste0(col[names(col)!=t],'40')
	mm$col = col[mm$tissue]
	for(p in c('u','d','ud','du')){
		sids = rownames(age.segs[[s]])[age.segs[[s]][,t]==p & anns[[s]]$sites==sites]
		plotTissueAgeProile(apply(psi.tsm[[s]][sids,],2,mean,na.rm=T),mm,main=paste0(p,' (',length(sids),')'),bty='n',xaxt='n',xlab='',ylab='PSI')
		if(t == 'brain') plotPanelLetter(c(u='D',d='E',ud='F',du='G')[p],lab.cex)
	}
}
dev.off()
# _s20 4patterns stat #####
age.segs = readRDS('Rdata/devAS.4patt.Rdata')
age.segs = readRDS('Rdata/devAS.4patt.rank.Rdata')
age.segs = readRDS('Rdata/devAS.4patt.rare.Rdata')

up = sapply(rownames(species),function(s){
	pp = age.segs[[s]][anns[[s]]$sites=='ad',]	
	stat = apply(pp,2,function(p)c(up=sum(p=='u'),down=sum(p=='d'),up_down=sum(p=='ud'),down_up=sum(p=='du')))
	stat['up',]/apply(stat,2,sum)
	})

down = sapply(rownames(species),function(s){
	pp = age.segs[[s]][anns[[s]]$sites=='ad',]	
	stat = apply(pp,2,function(p)c(up=sum(p=='u'),down=sum(p=='d'),up_down=sum(p=='ud'),down_up=sum(p=='du')))
	stat['down',]/apply(stat,2,sum)
})

range(up,na.rm=T)
range(down,na.rm=T)
range(down+up,na.rm=T)
sites ='ad'

#pdf('figures/paper.figures/6/suppl/S20.4patterns.counts.bars.rare.pdf',w=9,h=6,family='Arial')
pdf('figures/paper.figures/6/suppl/02NG/S16.4patterns.counts.bars.rare.pdf',w=9,h=6,family='Arial')
#png('figures/paper.figures/6/suppl/S20.4patterns.counts.bars.png',units = 'in',w=9,h=6,res = 600,family='Arial')
par(mfrow=c(3,3),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(4.5,2.5,1.5,0),oma=c(0,0,1,1))
for(s in rownames(species)){
	pp = age.segs[[s]][anns[[s]]$sites==sites,]
	stat = apply(pp,2,function(p)c(up=sum(p=='u'),down=sum(p=='d'),up_down=sum(p=='ud'),down_up=sum(p=='du')))
	b=barplot(stat,beside=T,border=NA,main=s,ylab='# of exons',col=rep(params$tissue.col[colnames(stat)],each=4),names.arg = rep('',7))
	
	segments(b,grconvertY(0.0,'npc','user'),b,grconvertY(-0.02,'npc','user'),xpd=T)
	xx = b[,4]
	xx = xx + (xx - mean(xx))

	segments(b[,4],grconvertY(-0.02,'npc','user'),xx,grconvertY(-0.05,'npc','user'),xpd=T)
	text(xx,grconvertY(-0.05,'npc','user'),c('up','down','up-down','down-up'),srt=90,adj = c(1.1,0.5),xpd=T,cex=1)
}
plot.new()
legend('topleft',border=NA,bty='n',fill=params$tissue.col,legend=names(params$tissue.col))
dev.off()

pdf('figures/paper.figures/6/suppl/02NG/S14.4patterns.counts.bars.pdf',w=9,h=10,family='Arial')
par(mfcol=c(7,3),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(4.5,2.5,1.5,0),oma=c(0,0,1,1))
for(i in 1:3){
	if(i==1) age.segs = readRDS('Rdata/devAS.4patt.Rdata')
	if(i==2) age.segs = readRDS('Rdata/devAS.4patt.rare.Rdata')
	if(i==3) age.segs = readRDS('Rdata/devAS.4patt.rank.Rdata')
	
	up = sapply(rownames(species),function(s){
		pp = age.segs[[s]][anns[[s]]$sites=='ad',]	
		stat = apply(pp,2,function(p)c(up=sum(p=='u'),down=sum(p=='d'),up_down=sum(p=='ud'),down_up=sum(p=='du')))
		stat['up',]/apply(stat,2,sum)
	})
	
	down = sapply(rownames(species),function(s){
		pp = age.segs[[s]][anns[[s]]$sites=='ad',]	
		stat = apply(pp,2,function(p)c(up=sum(p=='u'),down=sum(p=='d'),up_down=sum(p=='ud'),down_up=sum(p=='du')))
		stat['down',]/apply(stat,2,sum)
	})
	
	for(s in rownames(species)){
		pp = age.segs[[s]][anns[[s]]$sites==sites,]
		stat = apply(pp,2,function(p)c(up=sum(p=='u'),down=sum(p=='d'),up_down=sum(p=='ud'),down_up=sum(p=='du')))
		b=barplot(stat,beside=T,border=NA,main=s,ylab='# of exons',col=rep(params$tissue.col[colnames(stat)],each=4),names.arg = rep('',7))
		
		segments(b,grconvertY(0.0,'npc','user'),b,grconvertY(-0.02,'npc','user'),xpd=T)
		xx = b[,4]
		xx = xx + (xx - mean(xx))
		
		segments(b[,4],grconvertY(-0.02,'npc','user'),xx,grconvertY(-0.05,'npc','user'),xpd=T)
		text(xx,grconvertY(-0.05,'npc','user'),c('up','down','up-down','down-up'),srt=90,adj = c(1.1,0.5),xpd=T,cex=1)
		if(s == rownames(species)[1]) plotPanelLetter(LETTERS[i],lab.cex)
	}
	text(grconvertX(c(1,3,5)/6,'ndc','user'),grconvertY(1,'ndc','user'),c('All','Rarefaction','Rank'),adj=c(0.5,1),xpd=NA,cex=1.5)
}
dev.off()


# _comp 2 patt ####
a = readRDS('Rdata/devAS.4patt.Rdata')
b = readRDS('Rdata/devAS.4patt.rare.Rdata')

r = NULL
for(s in rownames(species))
	for(t in unique(meta$tissue)){
		r = rbind(r,data.frame(species=s,tissue=t,a=a[[s]][,t],b=b[[s]][,t]))
	}
o = c('-','n','u','d','ud','du')
table(r$a,r$b)[o,o]

# _table S1 ####
#  pc-lncRNA - possible coding
all.anns = readRDS('Rdata/all.anns.plus.cds.pos.Rdata')
ens.exon.transc.cds.pos = readRDS('Rdata/ens.exon.transc.cds.pos.Rdata')
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
	ss = all.anns[[s]]$sites
	ss[ss == 'ad'] = ifelse(all.anns[[s]]$type[ss == 'ad']=='ALT','ad','cnst')
	res[[s]] = table(t,ss)[types,c('cnst','ad','aa','dd','da')]
}

sapply(res,function(x)x[,'aa'])
(function(){
	cat("Constitutive exons\n")
	write.table(sapply(res,function(x)x[,'cnst']),sep='\t',quote = F)
	cat("Cassette exons\n")
	write.table(sapply(res,function(x)x[,'ad']),sep='\t',quote = F)
	cat("Alternative acceptor sites\n")
	write.table(sapply(res,function(x)x[,'aa']),sep='\t',quote = F)
	cat("Alternative donor sites\n")
	write.table(sapply(res,function(x)x[,'dd']),sep='\t',quote = F)
	cat("Retained introns\n")
	write.table(sapply(res,function(x)x[,'da']),sep='\t',quote = F)
})()

# _% of ens covered #####
h = loadSAData('processed/annotation/all.species/ensambl/Homo_sapiens.GRCh37.73.sajr.gz')
h = setSplSiteTypes(h,'processed/annotation/all.species/ensambl/Homo_sapiens.GRCh37.73.sajr.gz')
h = h$seg
t = table(paste(h$type,h$sites),paste(h$chr_id,h$start,h$stop,h$strand) %in% paste(all.anns$human$chr_id,all.anns$human$start,all.anns$human$stop,all.anns$human$strand))
t
t[c('ALT ad','ALT dd','ALT aa','INT da'),]
r = getAnnOverlap(h[h$position=='INTERNAL' & h$type !='EXN',],all.anns$human,return.matches = F)
table(h[h$position=='INTERNAL' & h$type !='EXN','sites'],r)

hei = loadIntrons('processed/annotation/all.species/ensambl/Homo_sapiens.GRCh37.73.sajr.gz')
hmi = loadIntrons('processed/annotation/all.species/merged/human.sajr.gz')
table(paste(hei$chr_id,hei$start,hei$stop,hei$strand) %in% paste(hmi$chr_id,hmi$start,hmi$stop,hmi$strand))
dim(hei)
245197/382046

# additional #####
# _same dir devAS in tissues ######

z=sapply(rownames(species),function(s){
	x = age.segs[[s]][anns[[s]]$sites=='ad',]
	#x[x!='u' & x!='d'] = '-'
	x[x=='n'] = '-'
	sgnc = apply(x!='-',1,sum)
	dir.cnt = apply(x,1,function(z)length(unique(z[z!='-']))) 
	c(sgn1=sum(sgnc>0),sgn21=sum(sgnc > 1 & dir.cnt==1),sgn22=sum(sgnc > 1 & dir.cnt>1))
	})
z

1-z[3,]/(z[3,]+z[2,])
1-z[3,]/(z[3,]+z[2,] + z[1,])
sum(z[2,])/sum(z[2:3,])


age.segs = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,psi.thr = 0.2,border.stages,s)[anns[[s]]$sites=='ad',])
names(age.segs) = rownames(species)
for(s in names(age.segs)) age.segs[[s]][is.na(per.tissue.age.qv[[s]][rownames(age.segs[[s]]),colnames(age.segs[[s]])])] = '-'
for(s in names(age.segs)) age.segs[[s]][age.segs[[s]] != '-' & per.tissue.age.qv[[s]][rownames(age.segs[[s]]),colnames(age.segs[[s]])]>0.05] = 'n'


z=lapply(age.segs,devASSameDir)
c = t(sapply(z,function(x){
	f = x[upper.tri(x)]
	c = t(x)[upper.tri(x)]
	c(sum(f*c),sum(c))
}))
sum(c[,1])/sum(c[,2])

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

# _per species hex p-value ######
hex.dws = readRDS('Rdata/hex.dws.age02sgn.Rdata')
hex.ups = readRDS('Rdata/hex.ups.age02sgn.Rdata')

f = function(m,t){
	rbind(incl.ups = hex.ups$up$pv[m,t,],
				incl.dws = hex.dws$up$pv[m,t,],
				excl.ups = hex.ups$dw$pv[m,t,],
				excl.dws = hex.dws$dw$pv[m,t,])[,7:1]
}

pv.b= f('actaac','brain')
pv.h= f('actaac','heart')

hex.ups$up$qv = apply(hex.ups$up$pv,2:3,p.adjust,m='BH')
table(hex.ups$up$qv[,1,1] - p.adjust(hex.ups$up$pv[,1,1],m='BH'))

# check ens orth vs mine ######

# orth.seg.ad.all = readRDS('Rdata/orth.seg.ad.all.Rdata')
# orth.all.ann = lapply(orth.seg.ad.all,function(x)x$seg)
# saveRDS(orth.all.ann,'Rdata/orth.seg.ad.all.ann.Rdata')
orth.all.ann = readRDS('Rdata/orth.seg.ad.all.ann.Rdata')
sapply(orth.all.ann,dim)
table(orth.all.ann$human$north)
table(orth.all.ann$human$north %in% c(1:6))
# so, out of 83888 orth exons, 
# 60596 come from one-to-one ensembl orthologs; 
# 19405 - in neither species belong to one-to-one ensembl orthologs
# 3887 - in some specie
set.seed(1)
#inx=sample(which(orth.all.ann$human$north == 6))
inx=sample(which(orth.all.ann$human$north == 0))

i=1
#inx=sample(which(orth.all.ann$human$north == 6))
#inx[1] 34657 ENSG00000173473, another gene in rabbit, ensembl one (ENSOCUG00000025458) looks like pseudogene, in Margarita correct one is just missed
#inx[2] 22471 ENSG00000188385, again rabbit, mine is correct (it is what currently used in ensembl and by margarida)
#inx[3] 77225 ENSG00000142798, my exon is not linked to ens gene in chicken; there are no gene in Margarida. But, according to ensembl mRNA and unipot alignemnts it is correct place, and it fits new assembly!


#inx=sample(which(orth.all.ann$human$north == 1))
#inx[1] 68630 ENSG00000254093, there are no ens orthologs in chicken, actually it seems to be annotation issue, gene is there and my chicken exon seems to be correct; margarida also do not have chicken here


#inx=sample(which(orth.all.ann$human$north == 0))
#inx[1] 75600 ENSG00000184203 retroduplication in Opossum, mine seems to be correct; also missed in Margarida DB
#inx[2] 20620 ENSG00000156970 retroduplication in rabbit, presedt in Margarida (but without rabbit)

z=sapply(rownames(species),function(s){
	seg2ens[[s]][[rownames(orth.all.ann[[s]])[inx[i]]]]
	})
z
z==orth.ens.genes[z[[1]],]
o.[o.[,1]==z[[1]],]


do.call(rbind,lapply(rownames(species),function(s){
	gene.descrs[[s]][z[[s]],]
	}))

do.call(rbind,lapply(rownames(species),function(s){
	gene.descrs[[s]][orth.ens.genes[z[[1]],s],]
}))
par(mfrow=c(1,2))
plotTissueAgeProile(ens.ge.marg.tsm$rabbit['ENSOCUG00000001940',],meta.tsm) #mine
plotTissueAgeProile(ens.ge.marg.tsm$rabbit['ENSOCUG00000025458',],meta.tsm) #mine

ens.ge$chicken$gene['ENSGALG00000016646',]
orth.all.ann$chicken[68630,]


# allow missed
o. = read.csv('input/ens.orths.txt.gz')
o = o.[apply(o.!='',1,sum)>1,]
dim(o)
f = rep(T,nrow(o))
for(i in 1:ncol(o)){
	t = table(o[,i])
	f = f & o[,i] %in% c('',names(t)[t==1])
}
table(f)
o = o[f,]

xorth.ens.genes = o
rownames(xorth.ens.genes) = xorth.ens.genes[,1]
colnames(xorth.ens.genes) = rownames(species)
table(orth.ens.genes[,1] %in% xorth.ens.genes[,1])


# r = matrix(NA,ncol=2,nrow=nrow(orth.all.ann$human))
# eidx = unlist(xorth.ens.genes)
# eidx = setNames(rep(1:nrow(xorth.ens.genes),times=ncol(xorth.ens.genes)),eidx)
# length(eidx)
# table(names(eidx) == '')
# eidx = eidx[names(eidx) != '']
# 
# for(i in 1:nrow(r)){
# 	cat('\r',i)
# 	egids = unlist(lapply(names(orth.all.ann),function(s){seg2ens[[s]][[rownames(orth.all.ann[[s]])[i]]]}))
# 	egids = egids[egids != '']
# 	t = sort(table(eidx[egids]),decreasing=T)
# 	if(length(t) > 0)
# 	r[i,] = c(t[1],sum(xorth.ens.genes[as.numeric(names(t)[1]),]!=''))
# }
#saveRDS(r,'xnorth.Rdata')
r = readRDS('xnorth.Rdata')
r[1:10,]
dim(r)
apply(is.na(r),2,sum)

r[is.na(r)] = 0

table(r[,1],total=r[,2])
table(orth.all.ann$human$north)
table(orth.all.ann$human$north,r[,2])
inx=which(orth.all.ann$human$north==7 & r[,2]==0)
i=1
z=sapply(rownames(species),function(s){
	seg2ens[[s]][[rownames(orth.all.ann[[s]])[inx[i]]]]
})
z
z==orth.ens.genes[z[[1]],]
xorth.ens.genes[z[[1]],]

o.[o.[,6]==z[6],]

do.call(rbind,lapply(rownames(species),function(s){
	gene.descrs[[s]][z[[s]],]
}))

do.call(rbind,lapply(rownames(species),function(s){
	gene.descrs[[s]][orth.ens.genes[z[[1]],s],]
}))


# microexons in heart, kidney #####
# are the any that work in hert only? not in brain?

adj.segments = readRDS('Rdata/adj.segments.sids-n-psi.Rdata')
s = 'mouse'
t = 'kidney'
f = anns[[s]]$sites == 'ad' & anns[[s]]$length<28 & per.tissue.age.qv[[s]][,t] < 0.05 & age.dpsi[[s]][,t] > 0.5
f[is.na(f)] = FALSE
sum(f)

pdf('figures/microexons.not.brain/all.mouse.kidney.devAS.microexons.pdf',w=6,h=6)
par(mfrow=c(3,3),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(2.5,2.5,1.5,1),oma=c(0,0,0,0))
for(i in which(f))
	plotTissueAgeProile(psi.tsm[[s]][i,],meta.tsm,main=rownames(anns[[s]])[i])
dev.off()

anns[[s]][f,]
adj.segments[[s]][rownames(anns[[s]])[f],]

plotTissueAgeProile(psi.tsm[[s]]['mou.35838.s3',],meta.tsm)
anns[[s]]['mou.35838.s3',]

cbind(adj.segments[[s]][rownames(anns[[s]])[f],],rownames(anns[[s]])[f] %in% rownames(orth.seg.ad.tsm$mouse))
# mou.35838.s3: do not have dd and aa and orthologous
oinx = which(rownames(orth.seg.ad.tsm$mouse) == 'mou.35838.s3')

gene.descrs$mouse[seg2ens$mouse[['mou.35838.s3']],]

pdf('figures/microexons.not.brain/Tmed2.mou.35838.s3.orth.heart.microexons.pdf',w=6,h=6)
par(mfrow=c(3,3),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(2.5,2.5,1.5,1),oma=c(0,0,0,0))
for(ss in rownames(species))
	plotTissueAgeProile(orth.seg.ad.tsm[[ss]][oinx,],meta.tsm,main=ss)
dev.off()

s
id = 'mou.35838.s3'
bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$days>40 & meta$tissue=='heart'],'.bam')
hl = getReadCoverage(bams,anns[[s]][id,'chr_id'],anns[[s]][id,'start'] - 1000,anns[[s]][id,'stop']+4000,-anns[[s]][id,'strand'])
bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$days<20 & meta$tissue=='heart'],'.bam')
he = getReadCoverage(bams,anns[[s]][id,'chr_id'],anns[[s]][id,'start'] - 1000,anns[[s]][id,'stop']+4000,-anns[[s]][id,'strand'])
bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$days>40 & meta$tissue=='brain'],'.bam')
bl = getReadCoverage(bams,anns[[s]][id,'chr_id'],anns[[s]][id,'start'] - 1000,anns[[s]][id,'stop']+4000,-anns[[s]][id,'strand'])

pdf('figures/microexons.not.brain/Tmed2.mou.35838.s3.heart.mouse.cov.pdf',w=6,h=8)
par(mfrow=c(3,1),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(2.5,2.5,1.5,1),oma=c(0,0,0,0))
plotReadCov(he,plot.junc.only.within = T,min.junc.cov = 10,xlim=c(124543000,124547000),bty='n',main='Prenatal heart')
plotReadCov(hl,plot.junc.only.within = T,min.junc.cov = 10,xlim=c(124543000,124547000),bty='n',main='Postnatal heart')
plotReadCov(bl,plot.junc.only.within = T,min.junc.cov = 10,xlim=c(124543000,124547000),bty='n',main='Postnatal brain')
dev.off()


# kidney
gene.descrs$mouse[seg2ens$mouse[['mou.23510.s9']],]
id = 'mou.23510.s9'
bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$days>40 & meta$tissue=='kidney'],'.bam')
kl = getReadCoverage(bams,anns[[s]][id,'chr_id'],anns[[s]][id,'start'] - 5000,anns[[s]][id,'stop']+6000,-anns[[s]][id,'strand'])
bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$days<13 & meta$tissue=='kidney'],'.bam')
ke = getReadCoverage(bams,anns[[s]][id,'chr_id'],anns[[s]][id,'start'] - 5000,anns[[s]][id,'stop']+6000,-anns[[s]][id,'strand'])
bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue=='brain'],'.bam')
b = getReadCoverage(bams,anns[[s]][id,'chr_id'],anns[[s]][id,'start'] - 5000,anns[[s]][id,'stop']+6000,-anns[[s]][id,'strand'])
bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue=='heart'],'.bam')
h = getReadCoverage(bams,anns[[s]][id,'chr_id'],anns[[s]][id,'start'] - 5000,anns[[s]][id,'stop']+6000,-anns[[s]][id,'strand'])


orth.ads.all.sp[!is.na(orth.ads.all.sp[,'mouse']) & orth.ads.all.sp[,'mouse']=='mou.23510.s9',]



pdf('figures/microexons.not.brain/Papss2.mou.23510.s9.orth.kidney.microexons.pdf',w=6,h=2)
par(mfrow=c(1,3),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(2.5,2.5,1.5,1),oma=c(0,0,0,0))
plotTissueAgeProile(psi.tsm$mouse['mou.23510.s9',],meta.tsm,main='Mouse')
plotTissueAgeProile(psi.tsm$rat['rat.4535.s8',],meta.tsm,main='Rat')
plotTissueAgeProile(psi.tsm$rabbit['rab.16045.s17',],meta.tsm,main='Rabbit')
dev.off()

pdf('figures/microexons.not.brain/Papss2.mou.23510.s9.kidney.mouse.cov.pdf',w=8,h=11)
par(mfrow=c(4,1),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(2.5,2.5,1.5,1),oma=c(0,0,0,0))
plotReadCov(ke,plot.junc.only.within = T,min.junc.cov = 10,bty='n',main='Prenatal kidney')
plotReadCov(kl,plot.junc.only.within = T,min.junc.cov = 10,bty='n',main='Postnatal kidney')
plotReadCov(b,plot.junc.only.within = F,min.junc.cov = 10,bty='n',main='Brain')
plotReadCov(h,plot.junc.only.within = F,min.junc.cov = 10,bty='n',main='Heart')
dev.off()

gene.descrs$mouse[seg2ens$mouse[['mou.46908.s35']],]
pdf('figures/microexons.not.brain/Myo6.mou.46908.s35.orth.kidney.microexons.pdf',w=6,h=6)
sids = orth.ads.all.sp[!is.na(orth.ads.all.sp[,'mouse']) & orth.ads.all.sp[,'mouse']=='mou.46908.s35',]
par(mfrow=c(3,3),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(2.5,2.5,1.5,1),oma=c(0,0,0,0))
for(s in rownames(species))
	if(!is.na(sids[s]))
		plotTissueAgeProile(psi.tsm[[s]][sids[s],],meta.tsm,main=s)
dev.off()




# age dPSI #####
# age.diamsss4 = lapply(rownames(species),function(s){
# 	print(s)
# 	p = readRDS(paste0("Rdata/",s,".as.u.filtered.Rdata"))$ir
# 	p = p[,colnames(p) %in% rownames(meta)[!(meta$stage %in% c('oldermidage','senior'))]]
# 	m = meta[colnames(p),]
# 	r = matrix(NA,nrow=nrow(p),ncol=7,dimnames = list(rownames(p),unique(meta$tissue)))
# 	for(t in colnames(r)){
# 		cat('  ',t,'\t')
# 		if(sum(m$tissue==t)>0){
# 			z = p[,m$tissue==t]
# 			m. = m[colnames(z),]
# 			r[,t] = apply(z,1,function(y)getDiamBySpline(m.$age.use,y,4))
# 		}
# 	}
# 	gc()
# 	r
# })
# names(age.diamsss4) = rownames(species)
# age.diamsss4 = lapply(age.diamsss4,function(x)apply(x,1:2,min,1))
# saveRDS(age.diamsss4,'Rdata/age.diam.spline4.with.replicates.Rdata')

# my.sex.maturation = setNames(c(1740,530,23,24,44,43,56),rownames(species)) #according to Margarida suggestions
# age.diamsss4.bm = lapply(rownames(species),function(s){
# 	print(s)
# 	p = readRDS(paste0("Rdata/",s,".as.u.filtered.Rdata"))$ir
# 	p = p[,colnames(p) %in% rownames(meta)[!(meta$stage %in% c('oldermidage','senior')) & meta$days <= my.sex.maturation[s]]]
# 	m = meta[colnames(p),]
# 	r = matrix(NA,nrow=nrow(p),ncol=7,dimnames = list(rownames(p),unique(meta$tissue)))
# 	for(t in colnames(r)){
# 		cat('  ',t,'\t')
# 		if(sum(m$tissue==t)>0){
# 			z = p[,m$tissue==t]
# 			m. = m[colnames(z),]
# 			r[,t] = apply(z,1,function(y)getDiamBySpline(m.$age.use,y,4))
# 		}
# 	}
# 	gc()
# 	r
# })
# names(age.diamsss4.bm) = rownames(species)
# age.diamsss4.bm = lapply(age.diamsss4.bm,function(x)apply(x,1:2,min,1))
# saveRDS(age.diamsss4.bm,'Rdata/age.diam.spline4.with.replicates.before.maturation.Rdata')
# 
# orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
# orth.age.dpsi = lapply(rownames(species),function(s){
# 		print(s)
# 		p = orth.seg.ad[[s]]$ir
# 		p = p[,colnames(p) %in% rownames(meta)[!(meta$stage %in% c('oldermidage','senior'))]]
# 		m = meta[colnames(p),]
# 		r = matrix(NA,nrow=nrow(p),ncol=7,dimnames = list(rownames(p),unique(meta$tissue)))
# 		for(t in colnames(r)){
# 			cat('  ',t,'\t')
# 			if(sum(m$tissue==t)>0){
# 				z = p[,m$tissue==t]
# 				m. = m[colnames(z),]
# 				r[,t] = apply(z,1,function(y)getDiamBySpline(m.$age.use,y,4))
# 			}
# 		}
# 		gc()
# 		r
# 	})
# names(orth.age.dpsi) = rownames(species)
# orth.age.dpsi = lapply(orth.age.dpsi,function(x)apply(x,1:2,min,1))
# saveRDS(orth.age.dpsi,'Rdata/orth.ad.age.diam.spline4.with.replicates.Rdata')

compTwoAmps = function(qv,a1,a2,t1,t2,sites='ad'){
	r = list()
	r$lost2=r$lost1=r$s2=r$s1=matrix(NA,ncol=7,nrow=7,dimnames = list(rownames(species),unique(meta$tissue)))
	for(s in names(qv)){
		for(t in colnames(r$s1)){
			f = anns[[s]]$sites %in% sites & !is.na(qv[[s]][,t]) & qv[[s]][,t] <0.05
			a1. = abs(a1[[s]][,t])
			a2. = abs(a2[[s]][,t])
			a1.[is.na(a1.)] = 0
			a2.[is.na(a2.)] = 0
			r$s1[s,t] = sum(f & a1. > t1)
			r$s2[s,t] = sum(f & a2. > t2)
			r$lost1[s,t] = sum(f & a1. < t1 & a2. > t2)
			r$lost2[s,t] = sum(f & a1. > t1 & a2. < t2)
		}
	}
	r
}


age.diamss4 = readRDS('Rdata/age.diamss4.Rdata')
age.diamsss4 = readRDS('Rdata/age.diam.spline4.with.replicates.Rdata')
age.diamsss5 = readRDS('Rdata/age.diam.spline5.with.replicates.Rdata')

age.diam = lapply(psi.tsm,function(x){
	m = meta.tsm[colnames(x),]
	r = matrix(NA,nrow=nrow(x),ncol=7,dimnames = list(rownames(x),unique(meta$tissue)))
	for(t in colnames(r)){
		if(sum(m$tissue==t)>0){
			z = x[,m$tissue==t]
			r[,t] = apply(z,1,max,na.rm=T)-apply(z,1,min,na.rm=T)
		}
	}	
	r
	})

age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)
age.dpsi$macaque = cbind(age.dpsi$macaque[,1:5],ovary=NaN,age.dpsi$macaque[,6,drop=FALSE])

age.diamss4$mouse[1:2,]
f = !is.na(per.tissue.age.qv$mouse[,1]) & per.tissue.age.qv$mouse[,1] < 0.05
ce = anns$mouse$sites == 'ad' & !is.na(per.tissue.age.qv$mouse[,1])
f1 = function(b,s,thrs)t(sapply(thrs,function(t)c(sum(b>t,na.rm=T),sum(s>t,na.rm=T))))
breaks=seq(0,1,length.out = 100)
dpsi.stat = f1(abs(age.dpsi$mouse[ce & !f,1]),abs(age.dpsi$mouse[ce & f,1]),breaks)
diam.stat = f1(abs(age.diam$mouse[ce & !f,1]),abs(age.diam$mouse[ce & f,1]),breaks)
diamss4.stat = f1(abs(age.diamss4$mouse[ce & !f,1]),abs(age.diamss4$mouse[ce & f,1]),breaks)
diamsss4.stat = f1(abs(age.diamsss4$mouse[ce & !f,1]),abs(age.diamsss4$mouse[ce & f,1]),breaks)
age.diamsss4.ar.stat = f1(abs(age.diamsss4.ar$mouse[ce & !f,1]),abs(age.diamsss4.ar$mouse[ce & f,1]),breaks)

#gplots::col2hex('orange')
par(mfrow=c(3,3))
hist(abs(age.dpsi$mouse[ce & !f,1]),seq(0,1,length.out = 100),border = NA,col='#66666666',freq = F,ylim=c(0,12),main='abs(PSIadult - PSIearliest)')
hist(abs(age.dpsi$mouse[ce & f,1]),seq(0,1,length.out = 100),border = NA,col='#0000FF66',add=T,freq = F)

hist(age.diam$mouse[ce & !f,1],seq(0,1,length.out = 100),border = NA,col='#66666666',freq = F,main='max(PSI)-min(PSI)')
hist(age.diam$mouse[ce & f,1],seq(0,1,length.out = 100),border = NA,col='#FF000066',add=T,freq = F)

hist(age.diamss4$mouse[ce & !f,1],seq(0,1,length.out = 100),border = NA,col='#66666666',freq = F,main='smooth max(PSI)-min(PSI)')
hist(age.diamss4$mouse[ce & f,1],seq(0,1,length.out = 100),border = NA,col='#FFA50066',add=T,freq = F)

hist(age.diamsss4$mouse[ce & !f,1],seq(0,1,length.out = 100),border = NA,col='#66666666',freq = F,main='smooth4 max(PSI)-min(PSI)')
hist(age.diamsss4$mouse[ce & f,1],seq(0,1,length.out = 100),border = NA,col='#00FF0066',add=T,freq = F)

hist(age.diamsss5$mouse[ce & !f,1],seq(0,1,length.out = 100),border = NA,col='#66666666',freq = F,main='smooth5 max(PSI)-min(PSI)')
hist(age.diamsss5$mouse[ce & f,1],seq(0,1,length.out = 100),border = NA,col='#00FFFF66',add=T,freq = F)

hist(age.diamsss4.ar$mouse[ce & !f,1],seq(0,1,length.out = 100),border = NA,col='#66666666',freq = F,main='smooth5 max(PSI)-min(PSI)')
hist(age.diamsss4.ar$mouse[ce & f,1],seq(0,1,length.out = 100),border = NA,col='#FF00FF66',add=T,freq = F)


plot(breaks,(dpsi.stat[,1]/dpsi.stat[1,1])/(dpsi.stat[,2]/dpsi.stat[1,2]),t='l',col='blue',log='y',xlab='threshould',ylab='exp/observed',lwd=3)
lines(breaks,(diam.stat[,1]/diam.stat[1,1])/(diam.stat[,2]/diam.stat[1,2]),col='red',lwd=3)
lines(breaks,(diamss4.stat[,1]/diamss4.stat[1,1])/(diamss4.stat[,2]/diamss4.stat[1,2]),col='orange',lwd=3)
lines(breaks,(diamsss4.stat[,1]/diamsss4.stat[1,1])/(diamsss4.stat[,2]/diamsss4.stat[1,2]),col='green',lwd=3)
lines(breaks,(diamsss5.stat[,1]/diamsss5.stat[1,1])/(diamsss5.stat[,2]/diamsss5.stat[1,2]),col='#00FFFF',lwd=3)
lines(breaks,(age.diamsss4.ar.stat[,1]/age.diamsss4.ar.stat[1,1])/(age.diamsss4.ar.stat[,2]/age.diamsss4.ar.stat[1,2]),col='#FF00FF',lwd=3)
abline(h=0.05,lty=3)
# spline 4 look bit better than spline 5, so I'll use it

apply(per.tissue.age.qv$mouse[ce,]<0.05 & abs(age.dpsi$mouse[ce,])>0.2,2,sum,na.rm=T)
apply(per.tissue.age.qv$mouse[ce,]<0.05 & age.diamss4$mouse[ce,]>0.2,2,sum,na.rm=T)
#apply(per.tissue.age.qv$mouse[ce,]<0.05 & abs(age.dpsi$mouse[ce,])<0.2 & age.diamss4$mouse[ce,]>0.2,2,sum,na.rm=T)
apply(per.tissue.age.qv$mouse[ce,]<0.05 & (is.na(age.dpsi$mouse[ce,]) | abs(age.dpsi$mouse[ce,])<0.2) & age.diamss4$mouse[ce,]>0.2,2,sum,na.rm=T)
apply(per.tissue.age.qv$mouse[ce,]<0.05 & abs(age.dpsi$mouse[ce,])>0.2 & age.diamss4$mouse[ce,]<0.2,2,sum,na.rm=T)

compTwoAmps(per.tissue.age.qv,age.diamsss4.ar,age.diamsss4,0.2,0.2)



s = 'mouse'
t = 'brain'
f = anns[[s]]$sites=='ad' & per.tissue.age.qv$mouse[,t]<0.05 & (is.na(age.dpsi$mouse[,t]) | abs(age.dpsi$mouse[,t])<0.2) & age.diamss4$mouse[,t]>0.2
f = anns[[s]]$sites=='ad' & per.tissue.age.qv$mouse[,t]<0.05  & age.diamss4$mouse[,t]>0.2

f[is.na(f)] = FALSE
p = psi.tsm[[s]][f,grep(t,colnames(psi.tsm[[s]]))]
dim(p)
cr = cor(t(p),u='p')
cr = cleanNADistMatrix(cr)
dim(cr)

clust = hclust(as.dist(1 - cr))
plot(clust)
cl = renameClustsBySize(cutree(clust,4))
par(mfrow=c(2,2))
for(i in 1:max(cl)){
	clp = apply(psi.tsm[[s]][names(cl)[cl==i],],2,mean,na.rm=T)
	plotTissueAgeProile(clp,meta.tsm,age.axis = 'rank',main=paste0('c',i,' (',sum(cl==i),')'))
}

table(cl,names(cl) %in% rownames(orth.seg.ad.all.tsm$mouse))
ff = anns[[s]]$sites=='ad' & per.tissue.age.qv$mouse[,t]<0.05 & abs(age.dpsi$mouse[,t])>0.2
ff[is.na(ff)] = FALSE
table(rownames(anns[[s]])[ff] %in% rownames(orth.seg.ad.all.tsm$mouse))

# devAS patterns #####

#pdf('figures/paper.figures/6/nature.review/03.1.rarefication.for.4patt.pdf',w=12,h=6,family='Arial')
pdf('figures/paper.figures/6/suppl/02NG/S15.rarefaction.samples.for.4patt.pdf',w=12,h=6,family='Arial')
par(mfrow=c(2,4),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,0,1))

for(s in rownames(species)){
	plot(1,t='n',xlab='log(DPC)',yaxt='n',ylab='',ylim=c(0.8,7.2),xlim=range(meta$age.use[meta$species==s]),bty='n',main=s)
	y = 7
	for(t in unique(meta$tissue)){
		a = meta$age.use[meta$species==s & meta$tissue==t]
		if(length(a)>0){
			f = rarefy2uniform(a,10,seed=19345)
			points(a[!f],runif(sum(!f),y-0.2,y+0.2),col=params$tissue.col[t],pch=1,cex=0.5)
			points(a[ f],runif(sum( f),y-0.2,y+0.2),col=params$tissue.col[t],pch=1,cex=2,lwd=2)
		}
		y = y -1
	}
}
dev.off()

select.ages = 'no'
# select.ages = 'rank'
# select.ages = 'rare'

devAS.change.stat = list()
for(s in rownames(species)){
	print(s)
	devAS.change.stat[[s]] = list()
	p = readRDS(paste0("Rdata/",s,".as.u.filtered.Rdata"))$ir
	p = p[,colnames(p) %in% rownames(meta)[!(meta$stage %in% c('oldermidage','senior'))]]
	for(t in unique(meta$tissue)){
		cat(' ',t,'\n')
		f = per.tissue.age.qv[[s]][,t] <= 0.05 & age.dpsi[[s]][,t] > 0.2
		f[is.na(f)] = FALSE
		if(sum(f)==0){
			devAS.change.stat[[s]][[t]] = matrix(1,ncol=4,nrow=0)
			colnames(devAS.change.stat[[s]][[t]]) = c('up','down','up.time','down.time')
			next
		}
		m = meta[colnames(p),]
		pp = p[f,m$tissue==t]
		m = meta[colnames(pp),]
		if(select.ages == 'rank'){
			stages = sort(sapply(split(m$age.use,m$stage),mean))
			stages[] = 1:length(stages)
			m$age.use = stages[m$stage]
		}else if(select.ages == 'rare'){
			f = rarefy2uniform(m$age.use,n = 10,seed = 19345)
			m = m[f,]
			pp = pp[,f]
		}
		devAS.change.stat[[s]][[t]] = t(apply(pp,1,function(y)getDevASPattern(m$age.use,y,4)))
	}
	gc()
}


# if(select.ages=='rank'){
# 	saveRDS(devAS.change.stat,'Rdata/devAS.change.stat.rank.Rdata')
# }else if(select.ages=='rare'){
# 	saveRDS(devAS.change.stat,'Rdata/devAS.change.stat.rare.Rdata')
# }else
# 	saveRDS(devAS.change.stat,'Rdata/devAS.change.stat.Rdata')

devAS.change.stat = readRDS('Rdata/devAS.change.stat.Rdata')
#devAS.change.stat = readRDS('Rdata/devAS.change.stat.rank.Rdata')
#devAS.change.stat = readRDS('Rdata/devAS.change.stat.rare.Rdata')

thr = 0.3
devAS.4patt = list()
for(s in rownames(species)){
	devAS.4patt[[s]] = matrix('-',ncol=7,nrow=nrow(per.tissue.age.qv[[s]]))
	dimnames(devAS.4patt[[s]]) = dimnames(per.tissue.age.qv[[s]])
	for(t in colnames(per.tissue.age.qv[[s]])){
		devAS.4patt[[s]][!is.na(per.tissue.age.qv[[s]][,t]),t] = 'n'
		y = devAS.change.stat[[s]][[t]]
		y = y[apply(is.na(y),1,sum)==0 & (y[,1]+y[,2])>0,]
		u = y[,1]/(y[,1]+y[,2])
		tim = y[,3]/y[,1]-y[,4]/y[,2]
		devAS.4patt[[s]][rownames(y),t] = 'du'
		devAS.4patt[[s]][rownames(y)[!is.na(tim) & tim < 0],t] = 'ud'
		devAS.4patt[[s]][rownames(y)[u < thr],t] = 'd'
		devAS.4patt[[s]][rownames(y)[u > (1 - thr)],t] = 'u'
	}
}

apply(devAS.4patt$macaque,2,function(x)table(factor(x,levels=c('u','d','ud','du','n','-'))))
table(devAS.4patt$human[,1],anns$human$sites)[c('u','d','ud','du','n','-'),c('ad','da','dd','aa')]

# if(select.ages=='rank'){
# 	saveRDS(devAS.4patt,'Rdata/devAS.4patt.rank.Rdata')
# }else if(select.ages=='rare'){
# 	saveRDS(devAS.4patt,'Rdata/devAS.4patt.rare.Rdata')
# }else
# 	saveRDS(devAS.4patt,'Rdata/devAS.4patt.Rdata')

sites ='ad'
pdf('figures/paper.figures/6/devAS.change.dir.pdf',w=10,h=13,family='Arial')
par(mfrow=c(8,5),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,1,1))
for(s in rownames(species)){
	patt = devAS.4patt[[s]][anns[[s]]$sites==sites,]
	stat = apply(patt,2,function(p)table(factor(p[p %in% c('u','d','ud','du')],levels=c('u','d','ud','du'))))
	barplot(stat,beside=T,col=rep(params$tissue.col,each=4),border=NA,names.arg =rep('',7))
	plot.new()
	plot.new()
	plot.new()
	plot.new()
	for(t in unique(meta$tissue)){
		y = devAS.change.stat[[s]][[t]]
		if(nrow(y) == 0){
			for(i in 1:5)
				plot.new()
			next
		}
		y = y[anns[[s]][rownames(y),'sites']==sites,]
		u = y[,1]/(y[,1]+y[,2])
		tim = y[,3]/y[,1]-y[,4]/y[,2]
		tim[is.na(tim)] = 0
		plot(u,tim,col=params$tissue.col[t],xlab='up/(up+down)',ylab='up_time - down_time',main='')
		abline(h=0,col='gray',lty=1)
		abline(v=c(thr,1-thr),col='gray',lty=1)
		mm = meta.tsm
		col = params$tissue.col
		col[names(col)!=t] = paste0(col[names(col)!=t],'40')
		mm$col = col[mm$tissue]
		for(p in c('u','d','ud','du')){
			sids = rownames(patt)[patt[,t]==p]
			plotTissueAgeProile(apply(psi.tsm[[s]][sids,],2,mean,na.rm=T),mm,main=paste0(p,' (',length(sids),')'),bty='n',xaxt='n',xlab='',ylab='PSI')
		}
	}
	mtext(s,outer=T)
}
dev.off()

# try supervised clustering ######
get4Patterns = function(x){
	mn = min(x)
	mx = max(x)
	x. = (x-mn)/(mx-mn)
	cbind(up=x.,
				dw=1-x.,
				ud=1-4*(x.-0.5)^2,
				du=4*(x.-0.5)^2)
}

s = 'mouse'
t = 'testis'
f = anns[[s]]$sites=='ad' & per.tissue.age.qv$mouse[,t]<0.05  & age.diamsss4$mouse[,t]>0.2
f[is.na(f)] = FALSE
p = psi.tsm[[s]][f,grep(t,colnames(psi.tsm[[s]]))]
m = meta.tsm[colnames(p),]

p4 = get4Patterns(m$age.rank)
# plot(p4[order(m$age.use),4],t='p')
# plot(m$age.rank,p4[,4],t='p')

cor=cor(t(p),p4,u='p')
cl = colnames(cor)[apply(cor,1,which.max)]
table(cl)
names(cl) = rownames(cor)
par(mfrow=c(2,2))
for(i in colnames(p4)){
	clp = apply(psi.tsm[[s]][names(cl)[cl==i],],2,mean,na.rm=T)
	plotTissueAgeProile(clp,meta.tsm,age.axis = 'rank',main=paste0('c',i,' (',sum(cl==i),')'))
}

cor. = t(apply(cor,1,sort))
cor.[1:10,]
plot(cor.[,4],cor.[,3])
abline(a=0,b=1,col='red',)

hist(cor.[,4]-cor.[,3])


# check MDS for subsets of tissues #####
# s = 'mouse'
# all.cors = all.mdss = list()
# for(s in rownames(species)){
# 	print(s)
# 	d = readRDS(paste0('Rdata/',s,'.as.u.filtered.Rdata'))
# 	d = d$ir[d$seg$sites=='ad',]
# 	gc()
# 	
# 	m = meta[colnames(d),]
# 	tissues = unique(meta$tissue)
# 	
# 	cors = list()
# 	cors$bchlkot = cor(d[,m$tissue %in% tissues],u='p',m='p')
# 	cors$hlkot = cor(d[,m$tissue %in% tissues[3:7]],u='p',m='p')
# 	cors$hlko = cor(d[,m$tissue %in% tissues[3:6]],u='p',m='p')
# 	cors$lko = cor(d[,m$tissue %in% tissues[4:6]],u='p',m='p')
# 	cors$lk = cor(d[,m$tissue %in% tissues[4:5]],u='p',m='p')
# 	mdss = lapply(cors,function(x)cmdscale(1-x,k=2))
# 	all.cors[[s]] = cors
# 	all.mdss[[s]] = mdss
# }
# saveRDS(all.mdss,'Rdata/all.mdss.for.tissues.subsets.Rdata')

pdf('figures/paper.figures/6/MDS.substets.of.tissues.pdf',w=15,h=21,family='Arial')
par(mfrow=c(7,5),tck=-0.01,mgp=c(1,0.2,0),mar=c(2,2,1.5,0),oma=c(0,0,0,0))
for(s in rownames(species)){
	mdss = all.mdss[[s]] 
	for(i in 1:length(mdss)){
		m = meta[rownames(mdss[[i]]),]
		cex = (m$cex-min(m$cex))/(max(m$cex)-min(m$cex))
		plot(mdss[[i]],main=paste(s,names(mdss)[i]),xlab='Dim 1',ylab='Dim 2',col=m$col,pch=19,cex=0.2+cex*3)
	}
}
dev.off()

# S21 TE.human-n-mouse #####
# see new.exons.R


# Table S3: devAS list ######
t = do.call(rbind,lapply(rownames(species),function(s){
	r = cbind(seg.id=rownames(anns[[s]]),anns[[s]][,c('chr_id','start','stop','strand','sites')],age.segs[[s]],age.dpsi[[s]],per.tissue.age.qv[[s]])
	r$sites = c(ad='CE',aa='AA',dd='AD',da='RI')[r$sites]
	r$sites[is.na(r$sites)] = 'complex'
	colnames(r)[6] = 'as.type'
	r = cbind(r[,1,drop=FALSE],ens_id=sapply(seg2ens[[s]][r$seg.id],paste,collapse=','),r[,-1])
	colnames(r)[-1:-7] = paste0(rep(c('pattern.','dpsi.','bh.adj.pv.'),each=7),rep(colnames(age.segs[[s]]),times=3))
	write.csv(r,paste0('output/shiny/downloads/',s,'.devAS'))
	r
	}))

t[1:10,]
dim(t)
write.csv(t,'figures/paper.figures/6/table.S4.all.devAS.csv')

# check filter PSI in [0.1,0.9] #####
s = 'mouse'
t = readRDS(paste0("Rdata/",s,".as.u.all.Rdata"))

st = length(t)
t = t[t$seg$chr_id !='MT' & t$seg$position=='INTERNAL' & t$seg$type!='EXN',colnames(t$ir) %in% meta$name]
st[2] = length(t)
na = is.na(t$ir)
f = rep(FALSE,length(t))
for(tis in unique(meta$tissue)){
	cinx = meta[colnames(na),'tissue']==tis
	f = f | (apply(!na[,cinx],1,mean) > 0.6 & apply(t$ir[,cinx],1,function(x){x=x[!is.na(x)];sum(x>0.1 & x < 0.9)>3}))
}
st[3] = sum(f)

f = rep(FALSE,length(t))
for(tis in unique(meta$tissue)){
	cinx = meta[colnames(na),'tissue']==tis
	f = f | (apply(!na[,cinx],1,mean) > 0.6 & apply(t$ir[,cinx],1,function(x){x=x[!is.na(x)];sum(x>0.1) > 3 & sum(x < 0.9)>3}))
}
st[4] = sum(f)
st

# tissue devAS similarity #####
myCramer = function(t,...){
	ct = chisq.test(t,...)
	list(estimate=sqrt(chisq.test(t)$statistic/sum(t)/min(nrow(t)-1,ncol(t)-1)),p.value=ct$p.value)
}


f = function(p,levs= c('n','u','d','ud','du'),use.patterns=TRUE){
	e = matrix(NA,ncol=7,nrow=7)
	colnames(e) = rownames(e) = colnames(p)
	pv = e
	for(t1 in 1:6)
		for(t2 in (t1+1):7){
			f = (p[,t1] %in% levs) & (p[,t2] %in% levs)
			if(use.patterns)
				t = myCramer(table(p[f,t1],p[f,t2]))
			else
				t = myCramer(table(p[f,t1]=='n',p[f,t2]=='n'))
			e[t1,t2] = e[t2,t1] = t$estimate
			pv[t1,t2] = pv[t2,t1] = t$p.value
		}
	list(estimate=e,p.value=pv)
}

plotTissueDevASoverlap = function(p,...){
	library(RColorBrewer)
	# takes only segs that are a) devAS in at least one tissue b) tested in both compared
	p = p[apply(p !='-' & p !='n',1,sum)>0,]
	r = array(NA,dim=c(7,7,5),dimnames = list(colnames(p),colnames(p),c('or','pv','intersect-to-min','intersect-to-max','intersect-to-union')))
	for(t1 in 1:6)
		for(t2 in (t1+1):7){
			f = p[,t1] != '-' & p[,t2] != '-'
			ft = fisher.test(factor(p[f,t1]!='n',levels = c(FALSE,TRUE)),factor(p[f,t2]!='n',levels = c(FALSE,TRUE)))
			r[t1,t2,1] = r[t2,t1,1] = ft$estimate
			r[t1,t2,2] = r[t2,t1,2] = ft$p.value
			r[t1,t2,3] = r[t2,t1,3] = sum(p[f,t1]!='n' & p[f,t2]!='n')/min(sum(p[f,t1]!='n'),sum(p[f,t2]!='n'))
			r[t1,t2,4] = r[t2,t1,4] = sum(p[f,t1]!='n' & p[f,t2]!='n')/max(sum(p[f,t1]!='n'),sum(p[f,t2]!='n'))
			r[t1,t2,5] = r[t2,t1,5] = sum(p[f,t1]!='n' & p[f,t2]!='n')/sum(p[f,t1]!='n' | p[f,t2]!='n')
		}
	
	t = round(r[,,1],1)
	t[lower.tri(t)] = round(r[,,5][lower.tri(t)],2)
	diag(t) = toupper(substr(colnames(p),1,1))
	pv = r[,,2][lower.tri(t)]
	bh.thr = max(pv[p.adjust(pv,m='BH')<0.05])
	
	col = brewer.pal(9, 'RdYlBu')
	col = c('#EEEEEE',col[3:1],col[9:7],'#EEEEEE')
	pv = r[,,2]

	pv[!is.na(r[,,1]) & r[,,1]>1] = - pv[!is.na(r[,,1]) & r[,,1]>1]
	imageWithText(pv,t,names.as.labs = F,xlab='',ylab='',bty='n',xaxt='n',yaxt='n',breaks=c(-2,-bh.thr,-0.05/7/3,-1e-10,0,1e-10,0.05/7/3,bh.thr,2),col=col,...)	
	invisible(r)
}

pdf('figures/paper.figures/6/tissueDevASoverlap.no.macaque.pdf',w=7.2,h=7.2,family='Arial')
#png('figures/paper.figures/6/tissueDevASoverlap.no.macaque.png',units = 'in',w=7.2,h=7.2,res = 600,family='Arial')
par(mfrow=c(3,3),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(1,1,1.5,0.5),oma=c(0,0,0,1))
r = array(NA,dim=c(7,7,7,5),dimnames = list(rownames(species),unique(meta$tissue),unique(meta$tissue),c('or','pv','intersect-to-min','intersect-to-max','intersect-to-union')))
for(s in rownames(species)[-2]){
	sites = 'ad'
	p = age.segs[[s]][anns[[s]]$sites == 'ad', ]
	r[s,,,] = plotTissueDevASoverlap(p,main=s)
}

col = brewer.pal(9, 'RdYlBu')
plot.new()
legend('topleft' ,fill=c('#EEEEEE',col[3:1]),legend = c('not sign.','< BH 0.05','< Bonf 0.05','< 1e-10'),bty='n',border=NA,title='overrepresentation')
legend('topright',fill=c('#EEEEEE',col[7:9]),legend = c('not sign.','< BH 0.05','< Bonf 0.05','< 1e-10'),bty='n',border=NA,title='underrepresentation')
dev.off()
r[,'brain','cerebellum',]
r[,'kidney','liver',]

pdf('figures/paper.figures/6/brain.dPSI.2.tissue.cor.pdf',w=18,h=21,family='Arial')
par(mfrow=c(7,6),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(2.5,2.5,1.5,0.5),oma=c(0,0,0,1))
for(s in rownames(species)){
	f = anns[[s]]$sites=='ad' & per.tissue.age.qv[[s]][,'brain'] < 0.05 & age.dpsi[[s]][,'brain']>0.2
	for(t in unique(meta$tissue)[-1]){
		by = psi.tsm[[s]][f,paste(s,'brain',border.stages[[s]]['brain',1])]
		if(is.na(border.stages[[s]][t,2]))
			plot.new()
		else
			plotLine(psi.tsm[[s]][f,paste(s,'brain',border.stages[[s]]['brain',2])]-by,psi.tsm[[s]][f,paste(s,t,border.stages[[s]][t,2])] - by,col=params$tissue.col[t],line.col = 'black',
							 xlab='brain_late - brain_early',ylab=paste0(t,'_late - brain_early'),main=s,plot.ci = T,line.lwd = 3)
	}
}
dev.off()

btissue = 'liver'
pdf(paste0('figures/paper.figures/6/',btissue,'.dPSI.2.tissue.cor.pdf'),w=18,h=21,family='Arial')
par(mfrow=c(7,6),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(2.5,2.5,1.5,0.5),oma=c(0,0,0,1))
for(s in rownames(species)){
	f = anns[[s]]$sites=='ad' & per.tissue.age.qv[[s]][,btissue] < 0.05 & age.dpsi[[s]][,btissue]>0.2
	for(t in setdiff(unique(meta$tissue),btissue)){
		by = psi.tsm[[s]][f,paste(s,btissue,border.stages[[s]][btissue,1])]
		if(is.na(border.stages[[s]][t,2]))
			plot.new()
		else
			plotLine(psi.tsm[[s]][f,paste(s,btissue,border.stages[[s]][btissue,2])]-by,psi.tsm[[s]][f,paste(s,t,border.stages[[s]][t,2])] - by,col=params$tissue.col[t],line.col = 'black',
							 xlab=paste0(btissue,'_late - ',btissue,'_early'),ylab=paste0(t,'_late - ',btissue,'_early'),main=s,plot.ci = T,line.lwd = 3)
	}
}
dev.off()


r = array(NA,dim=c(7,7,7,4),dimnames=list(rownames(species),unique(meta$tissue),unique(meta$tissue),c('rho','ci1','ci2','pv')))
for(s in rownames(species)){
	for(bt in unique(meta$tissue)){
		if(is.na(border.stages[[s]][bt,2])) next
		f = anns[[s]]$sites=='ad' & per.tissue.age.qv[[s]][,bt] < 0.05 & age.dpsi[[s]][,bt]>0.2
		for(t in setdiff(unique(meta$tissue),bt)){
			by = psi.tsm[[s]][f,paste(s,bt,border.stages[[s]][bt,1])]
			if(is.na(border.stages[[s]][t,2])) next
			ct = cor.test(psi.tsm[[s]][f,paste(s,bt,border.stages[[s]][bt,2])]-by,psi.tsm[[s]][f,paste(s,t,border.stages[[s]][t,2])] - by,m='pear')
			r[s,bt,t,1] = ct$estimate
			r[s,bt,t,2:3] = ct$conf.int
			r[s,bt,t,4] = ct$p.value
		}
	}
}
range(r[3,,,4],na.rm=T)

library(RColorBrewer)

pdf('figures/paper.figures/6/devAStissue.cor.pdf',w=9,h=9,family='Arial')
par(mfrow=c(3,3),tck=-0.01,mgp=c(.2,0.3,0),mar=c(1.5,1.5,1.5,0.5),oma=c(0,0,0,1))
for(s in dimnames(r)[[1]]){
	col = brewer.pal(11, 'RdYlBu')
	col = c('#EEEEEE',col[5:1])
	pv = r[s,,,4]
	bh.thr = max(pv[p.adjust(pv,m='BH')<0.05],na.rm=T)
	t = round(r[s,,,1],2)
	diag(t) = toupper(substr(colnames(t),1,1))
	cor = r[s,,,1]
	cor[pv > 0.05/49] = 0
	imageWithText(t(cor)[,7:1],t(t)[,7:1],names.as.labs = F,xlab='compare with',ylab='DevAS tissue',bty='n',xaxt='n',yaxt='n',breaks=c(0,0.3,0.45,0.6,0.7,0.8,1),col=col,main=s)	
}
dev.off()

# _suppl orth seg table ######
orth.ads.all.sp[1:2,]
o = orth.ads.all.sp
for(s in colnames(o)){
	t = all.anns[[s]][o[!is.na(o[,s]),s],]
	o[!is.na(o[,s]),s] = paste0(rownames(t),':',t$chr_id,':',ifelse(t$strand==1,'+','-'),':',t$start,':',t$stop)
}
# write.csv(o,'output/paper/NG/supp/table.S8.orth.seg.csv',row.names=F)
# _load hgmd ######
hgmd = read.table('input/hgmd/tosend/2017_1_HGMD_ALL_Variants.csv',sep='\t',row.names=1,header=T,quote='',comment.char = '')
colnames(hgmd)
table(nchar(hgmd$ref_VCF_hg19),hgmd$tag)
hgmd.gr = GRanges(hgmd$chrom_VCF_hg19,IRanges(hgmd$pos_VCF_hg19,hgmd$pos_VCF_hg19+nchar(hgmd$ref_VCF_hg19)-1))
hann = all.anns$human
seg.gr = GRanges(hann$chr_id,IRanges(hann$start,hann$stop))
hgmd2hadf = findOverlaps(hgmd.gr,seg.gr,maxgap=200,type='any',select='all',ignore.strand=TRUE)
hgmd2hadf = data.frame(h=hgmd2hadf@from,s=hgmd2hadf@to)
hgmd2hadf$seg.id = rownames(hann)[hgmd2hadf$s]
hgmd2hadf$hgmd.id = rownames(hgmd)[hgmd2hadf$h]

# s = hann[hgmd2hadf$s,]
# h = hgmd[hgmd2hadf$h,]
# 
# hgmd2hadf$position = ''
# hgmd2hadf$position = paste0(hgmd2hadf$position,ifelse((s$strand== 1 & s$start > h$pos_VCF_hg19) | (s$strand==-1 & s$stop < h$pos_VCF_hg19+nchar(h$ref_VCF_hg19)-1),'u',''))
# hgmd2hadf$position = paste0(hgmd2hadf$position,ifelse(s$start<=h$pos_VCF_hg19+nchar(h$ref_VCF_hg19)-1 & s$stop>=h$pos_VCF_hg19,'i',''))
# hgmd2hadf$position = paste0(hgmd2hadf$position,ifelse((s$strand==-1 & s$start > h$pos_VCF_hg19) | (s$strand== 1 & s$stop < h$pos_VCF_hg19+nchar(h$ref_VCF_hg19)-1),'d',''))
# table(hgmd2hadf$position)
# # locaction: intron  (u,w), exon (e), ppt (p), acc.ag (A),acc(a), don(d),don.gt(D)
# # mutation annotation prioroti: dinucleotides > other nt of canonical sites > PPT > intron
# hgmd2hadf$loc = ''
# hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse(hgmdOverlapLocaction(h,s,0,Inf,T) & hgmdOverlapLocaction(h,s,0,Inf,F),'e',''))
# hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse(hgmdOverlapLocaction(h,s,-2,-1,T),'A',''))
# hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse(hgmdOverlapLocaction(h,s,-2,-1,F),'D',''))
# hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse((!grepl('A|D',hgmd2hadf$loc)) & (hgmdOverlapLocaction(h,s,-3,-3,T) | hgmdOverlapLocaction(h,s,1,1,T)),'a',''))
# hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse((!grepl('A|D',hgmd2hadf$loc)) & (hgmdOverlapLocaction(h,s,-5,-3,F) | hgmdOverlapLocaction(h,s,1,4,F)),'d',''))
# hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse((!grepl('A|D|a|d',hgmd2hadf$loc)) & hgmdOverlapLocaction(h,s,-23,-4,T),'p',''))
# hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse((!grepl('A|D|a|d|p',hgmd2hadf$loc)) & hgmdOverlapLocaction(h,s,-200,-24,T),'u',''))
# hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse((!grepl('A|D|a|d|p',hgmd2hadf$loc)) & hgmdOverlapLocaction(h,s,-200,-6,F),'w',''))

# _suppl alternif table ######
segcoor = c()
for(s in names(all.anns)){
	t = all.anns[[s]]
	segcoor = c(segcoor,setNames(paste0(rownames(t),':',t$chr_id,':',ifelse(t$strand==1,'+','-'),':',t$start,':',t$stop),rownames(t)))
}
segcoor[1:2]

gs = c(species$short,'hq','mr','mrb','hqmrb')
gs=c(gs,sapply(gs,sub,'','hqmrboc'))
length(gs)
length(unique(gs))

alt.sp[1:2]
sort(table(alt.sp))
table(alt.sp %in% gs)
alt.sp. = alt.sp[alt.sp %in% gs]
table(alt.sp.)[gs]


alternif = orth.ads.all.sp[!is.na(orth.ads.all.sp[,1]) & orth.ads.all.sp[,1] %in% names(alt.sp.),]
rownames(alternif) = alternif[,1]
alternif = cbind(alternif,alternative.in = alt.sp.[rownames(alternif)],HGMD.id=NA)
table(alternif[,8])
hg=sapply(split(hgmd2hadf$hgmd.id,hgmd2hadf$seg.id),paste,collapse=',')
cmn = intersect(names(hg),rownames(alternif))
length(cmn)
alternif[cmn,'HGMD.id'] = hg[cmn]
table(sapply(strsplit(hg[cmn],','),length))
table(alternif[,1:7] %in% names(segcoor))
for(s in rownames(species))
	alternif[,s] = segcoor[alternif[,s]]
#write.csv(alternif,'output/paper/NG/supp/table.S13.alternification.csv',row.names=F)

# _suppl newborn table ######
load('Rdata/tmp.exon.birth.Rdata')
exon.birth.one = readRDS('Rdata/exon.birth.one.Rdata')
sp.birth = sapply(exon.birth.one,function(x){paste(species$short[!is.na(x$seg_id)],collapse='')})
sort(table(sp.birth))
table(sp.birth %in% gs)

sort(table(sp.birth[!(sp.birth %in% gs)]))

len =     unlist(lapply(exon.birth.one,function(x)(min(x$length,na.rm=T))))
nb.stat[1:2,]
dim(nb.stat)
nbf = nb.stat$adj.exons== 0 & nb.stat$full.obs==1 & len < 500
table(nbf)

nbe = exon.birth.one[nbf & sp.birth %in% gs]
for(i in 1:length(nbe)){
	nbe[[i]]$orth.id = names(nbe)[i]
	nbe[[i]]$species = rownames(nbe[[i]])
}
nbe = do.call(rbind,nbe)
nbe$HGMD.id = NA
nbe = nbe[,c('species','seg_id','useg_id','dseg_id','orth.id','HGMD.id')]
rownames(nbe) = NULL
cmn = intersect(names(hg),nbe$seg_id[!is.na(nbe$seg_id)])
length(cmn)
nbe[match(cmn,nbe$seg_id),'HGMD.id'] = hg[cmn]
table(sapply(strsplit(hg[cmn],','),length))


segcoor[c(NA,'hum.6238.s1')]
table(nbe$seg_id[!is.na(nbe$seg_id)] %in% names(segcoor))
table(nbe$dseg_id %in% names(segcoor))

nbe$seg_id = segcoor[nbe$seg_id]
nbe$useg_id = segcoor[nbe$useg_id]
nbe$dseg_id = segcoor[nbe$dseg_id]
write.csv(nbe,'output/paper/NG/supp/table.S12.new.exons.csv',row.names=F)

# Nature review ######
# _sample heatmap #####


#pdf('figures/paper.figures/6/nature.review/01.2.sample.table.pdf',w=6,h=8.2,family='Arial')
pdf('figures/paper.figures/6/suppl/02NG/S01.sample.table.pdf',w=6,h=8.2,family='Arial')
par(tck=-0.01,mgp=c(2.1,0.4,0),mar=c(4,1,1.5,3),oma=c(0,0,0,0))
plotSampleTable1(meta,rownames(species)[c(1,2,4:7,3)],sspace = 4,cex=0.5)
dev.off()

# pdf('figures/paper.figures/6/nature.review/01.sample.table.pdf',w=10,h=7.2,family='Arial')
# par(mfrow=c(1,2),oma=c(0,0,0,0),mar=c(4,0,1,1))
# plotSampleTable(T,main='Exactly matched stages')
# plotSampleTable(F,main='All uniq samples in-between')
# dev.off()

# _cerebellum microexons #####
t = 'brain'
t = 'cerebellum'
micro.timing.cbm = lapply(rownames(species)[-c(2,7)], function(s){
	f = anns[[s]]$sites=='ad' & age.segs[[s]][,t] == 'u'
	f = !is.na(f) & f
	clns = paste(s,t,c(border.stages[[s]][t,1],age.al.i[10,s],border.stages[[s]][t,2]))
	if(sum(!(clns %in% colnames(psi.tsm[[s]])))>0) return(NULL)
	t = psi.tsm[[s]][f,clns]
	t = table(micro=anns[[s]]$length[f]<=27,before.birth = factor((t[,2]-t[,1])>(t[,3]-t[,2]),levels=c(F,T)))
})
names(micro.timing.cbm) = rownames(species)[-c(2,7)]
micro.timing.cbm = micro.timing.cbm[!sapply(micro.timing.cbm,is.null)]

(micro.timing.cbm.b=sapply(micro.timing.cbm,function(x)my.binom.test(x[2,2:1]))*100)
(macro.timing.cbm.b=sapply(micro.timing.cbm,function(x)my.binom.test(x[1,2:1]))*100)


pdf('figures/paper.figures/6/nature.review/02.microexons.in.cerebellum.pdf',w=8,h=4,family='Arial')
par(mfrow=c(1,2),mar=c(3,2.5,1.5,1),tck=-0.02,mgp=c(1.1,0.2,0),oma=c(0,0,0,1))
x = seq(from=1,by=4,length.out=ncol(micro.timing.b))
plot(x,micro.timing.b[1,],pch=19,ylim=c(0,100),xlim=range(1,x+1),bty='n',xaxt='n',xlab='',ylab='% of exons mostly changed before birth',main='Microexons in brain')
points(x+1,macro.timing.b[1,],pch=1,cex=2)
arrows(x,micro.timing.b[2,],x,micro.timing.b[3,],angle=90,code=3,length=0.03,xpd=T)
arrows(x+1,macro.timing.b[2,],x+1,macro.timing.b[3,],angle=90,code=3,length=0.03,xpd=T)
axis(1,x+0.5,NA)
text(x+0.5,rep(-5,ncol(micro.timing.b)),colnames(micro.timing.b),srt=-45,adj=c(0,1),xpd=NA)


x = seq(from=1,by=4,length.out=ncol(micro.timing.cbm.b))
plot(x,micro.timing.cbm.b[1,],pch=19,ylim=c(0,100),xlim=range(1,x+1),bty='n',xaxt='n',xlab='',ylab='% of exons mostly changed before birth',main='Micorexons in cerebellum')
points(x+1,macro.timing.cbm.b[1,],pch=1,cex=2)
arrows(x,micro.timing.cbm.b[2,],x,micro.timing.cbm.b[3,],angle=90,code=3,length=0.03,xpd=T)
arrows(x+1,macro.timing.cbm.b[2,],x+1,macro.timing.cbm.b[3,],angle=90,code=3,length=0.03,xpd=T)
axis(1,x+0.5,NA)
text(x+0.5,rep(-5,ncol(micro.timing.cbm.b)),colnames(micro.timing.cbm.b),srt=-45,adj=c(0,1),xpd=NA)
dev.off()

t = 'brain'
micro.timing.ctx.sids = lapply(rownames(species)[-c(2,7)], function(s){
	f = anns[[s]]$sites=='ad' & age.segs[[s]][,t] == 'u'
	f = !is.na(f) & f
	clns = paste(s,t,c(border.stages[[s]][t,1],age.al.i[10,s],border.stages[[s]][t,2]))
	if(sum(!(clns %in% colnames(psi.tsm[[s]])))>0) return(NULL)
	t = psi.tsm[[s]][f,clns]
	data.frame(micro=anns[[s]]$length[f]<=27,before.birth = (t[,2]-t[,1])>(t[,3]-t[,2]))
})
names(micro.timing.ctx.sids) = rownames(species)[-c(2,7)]

pdf('figures/paper.figures/6/nature.review/03.brain.microexons.in.cerebellum.pdf',w=12.5,h=5,family='Arial')
par(mfcol=c(2,5),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,0,1))
for(s in names(micro.timing.ctx.sids)){
	h = micro.timing.ctx.sids[[s]]
	plotTissueAgeProile(apply(psi.tsm[[s]][rownames(h)[h$micro & !is.na(h$before.birth) & h$before.birth],],2,mean,na.rm=T),meta.tsm,main=paste0(s,' (before birth)'),xlab='Stage',age.axis = 'rank',ylab='PSI',bty='n')
	plotTissueAgeProile(apply(psi.tsm[[s]][rownames(h)[h$micro & !is.na(h$before.birth) & !h$before.birth],],2,mean,na.rm=T),meta.tsm,main=paste0(s,' (after birth)'),xlab='Stage',age.axis = 'rank',ylab='PSI',bty='n')
}
dev.off()

# _subsample to macaque ####
# see hqm.subsampling.R
# _macroexons #####

orth.mds2 = readRDS('Rdata/paper.figures/orth.mds5.Rdata')
m = meta[rownames(orth.mds2$sp7),]

astypes = c(CE='ad',AA='aa',AD='dd',RI='da')
tested.stat = lapply(setNames(astypes,astypes),function(ss)sapply(names(per.tissue.age.qv),function(s)apply(!is.na(per.tissue.age.qv[[s]][anns[[s]]$sites==ss,]),2,sum)))

sgn02.stat = lapply(setNames(astypes,astypes),function(ss) {
	sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][anns[[s]]$sites == ss,] < 0.05 & abs(age.dpsi[[s]][anns[[s]]$sites == ss,])>0.2,2,sum,na.rm=T))/tested.stat[[ss]]*100})

sgn.cnt=lapply(setNames(astypes,astypes),function(ss) {
	sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][anns[[s]]$sites == ss,] < 0.05 & abs(age.dpsi[[s]][anns[[s]]$sites == ss,])>0.2,2,sum,na.rm=T))})

tested.stat.macro = lapply(setNames(astypes,astypes),function(ss)sapply(names(per.tissue.age.qv),function(s)apply(!is.na(per.tissue.age.qv[[s]][anns[[s]]$sites==ss & anns[[s]]$length > 27,]),2,sum)))

sgn02.stat.macro = lapply(setNames(astypes,astypes),function(ss) {
	sapply(names(anns),function(s){f = anns[[s]]$sites == ss & anns[[s]]$length > 27
	apply(per.tissue.age.qv[[s]][f,] < 0.05 & abs(age.dpsi[[s]][f,])>0.2,2,sum,na.rm=T)})/tested.stat.macro[[ss]]*100})

sgn.cnt.macro=lapply(setNames(astypes,astypes),function(ss) {
	sapply(names(anns),function(s){
		f = anns[[s]]$sites == ss & anns[[s]]$length > 27
		apply(per.tissue.age.qv[[s]][f,] < 0.05 & abs(age.dpsi[[s]][f,])>0.2,2,sum,na.rm=T)})})


astypes.pchs=c(ad=19,aa=2,dd=6,da=13)


#pdf('figures/paper.figures/6/nature.review/S23.brain.macroexons.pdf',w=7.2,h=7.2/3*2,family='Arial')
pdf('figures/paper.figures/6/suppl/02NG/S21.macroexons.pdf',w=7.2,h=7.2/3*2,family='Arial')
par(mfcol=c(2,3),tck=-0.01,mgp=c(1.4,0.4,0),mar=c(2.5,2.5,1,0),oma=c(0,0,0,1))
plot(orth.mds2$sp7,xlab='Dim 1',ylab='Dim 2',col=params$tissue.col[m$tissue],pch=params$species.pch[m$species],cex=m$cex,bty='n',main='All')
plot(orth.mds2$sp7.macro,xlab='Dim 1',ylab='Dim 2',col=params$tissue.col[m$tissue],pch=params$species.pch[m$species],cex=m$cex,bty='n',main='Macroexons')
plotAsEventCount(sgn.cnt['ad'],astypes.pchs[1],by.tissue = T,ylab='# of devAS',main='DevAS all',bty='n',ylim=c(0,3900))
plotAsEventCount(sgn.cnt.macro['ad'],astypes.pchs[1],by.tissue = T,ylab='# of devAS',main='DevAS macroexons',bty='n',ylim=c(0,3900))
plotAsEventCount(sgn02.stat['ad'],astypes.pchs[1],by.tissue = T,ylab='% of devAS',main='DevAS all',bty='n')
plotAsEventCount(sgn02.stat.macro['ad'],astypes.pchs[1],by.tissue = T,ylab='% of devAS',main='DevAS macroexons',bty='n')
dev.off()

# _tau vs no-of-exons #####
human.tau = read.csv('input/gene.info.from.marg/Human.Indexes.All.csv')
rownames(human.tau) = human.tau[,1]
human.tau$biotype = setNames(ens.ge$human$gene$biotype,rownames(ens.ge$human$gene))[human.tau$Human_ID]

h.as.ens = loadSAData('processed/annotation/all.species/ensambl/Homo_sapiens.GRCh37.73.sajr.gz')
h.as.ens = setSplSiteTypes(h.as.ens,'processed/annotation/all.species/ensambl/Homo_sapiens.GRCh37.73.sajr.gz')
h.int.ex.cnt = sapply(split(h.as.ens$seg$sites,h.as.ens$seg$gene_id),function(x){sum(x %in% c('ad'))})


mouse.tau = read.csv('input/gene.info.from.marg/Mouse.Indexes.csv')
rownames(mouse.tau) = mouse.tau[,1]
mouse.tau$biotype = setNames(ens.ge$mouse$gene$biotype,rownames(ens.ge$mouse$gene))[mouse.tau$Mouse_ID]

m.as.ens = loadSAData('processed/annotation/all.species/ensambl/Mus_musculus.GRCm38.84.sajr.gz')
m.as.ens = setSplSiteTypes(m.as.ens,'processed/annotation/all.species/ensambl/Mus_musculus.GRCm38.84.sajr.gz')
m.int.ex.cnt = sapply(split(m.as.ens$seg$sites,m.as.ens$seg$gene_id),function(x){sum(x %in% c('ad'))})


# psi.tsm.ad = lapply(names(psi.tsm),function(s)psi.tsm[[s]][anns[[s]]$sites=='ad',])
# names(psi.tsm.ad) = names(psi.tsm)
# ts = unique(meta$tissue)
# human.tissue.gene.dpsi=lapply(setNames(ts,ts), function(t){print(t);t(getsPSIbyEnsID(age.dpsi$human[,t],seg2ens$human,use.random = F))})
# mouse.tissue.gene.dpsi=lapply(setNames(ts,ts), function(t){print(t);t(getsPSIbyEnsID(age.dpsi$mouse[,t],seg2ens$mouse,use.random = F))})
# saveRDS(human.tissue.gene.dpsi,'Rdata/paper.figures/human.tissue.gene.dpsi.Rdata')
# saveRDS(mouse.tissue.gene.dpsi,'Rdata/paper.figures/mouse.tissue.gene.dpsi.Rdata')

human.tissue.gene.dpsi = readRDS('Rdata/paper.figures/human.tissue.gene.dpsi.Rdata')
mouse.tissue.gene.dpsi = readRDS('Rdata/paper.figures/mouse.tissue.gene.dpsi.Rdata')


human.tau$exon.no = h.int.ex.cnt[human.tau$Human_ID]
mouse.tau$exon.no = m.int.ex.cnt[mouse.tau$Mouse_ID]


pdf('figures/paper.figures/6/suppl/02NG/S8.tau-no.exon.pdf',w=9,h=10,family='Arial')
#par(mfrow=c(4,6),tck=-0.01,mgp=c(1.4,0.4,0),mar=c(2.5,2.5,1,0),oma=c(0,0,0,1))
par(tck=-0.01,mgp=c(1.4,0.4,0),mar=c(2.5,2.5,1.3,0),oma=c(0,0,0,1))
layout(matrix(1:16,ncol=4,byrow = T),widths=c(1,3,3,3))
cols = c('darkgray','lightgray','orange')
ord = c('non-AS','AS','devAS')
taus = list(human=human.tau,mouse=mouse.tau)
dpsis = list(human=human.tissue.gene.dpsi,mouse=mouse.tissue.gene.dpsi)
lab = T
for(s in c('human','mouse')){
	for(t in c('brain','liver')){
		plot.new()
		par(mar=c(0,0,0,0),xpd=NA)
		plotPNG(paste0("figures/paper.figures/5/icons/",s,".png"),0.5,0.7,0.55)#,0.2,0.5,0.35)
		plotPNG(paste0("figures/paper.figures/5/icons/",t,".png"),0.5,0.3,0.55)#,0.8,0.5,0.35)
		par(mar=c(2.5,2.5,1.3,0),xpd=T)
		tau = taus[[s]][!is.na(taus[[s]]$TissueTau) & !is.na(taus[[s]]$biotype) & taus[[s]]$biotype =='protein_coding',]
		tau$as = ifelse(rownames(tau) %in% rownames(dpsis[[s]][[t]]),'AS','non-AS')
		tau$as[rownames(tau) %in% rownames(dpsis[[s]][[t]])[dpsis[[s]][[t]][,'up']>0.2]] = 'devAS'
		tau = tau[tau$exon.no>0,]
		tau$exon.no.char = tau$exon.no
		tau$exon.no.char[tau$exon.no.char>20] = '>20'
		tau$exon.no.char = factor(tau$exon.no.char,levels = c(1:20,'>20'))
		tau$TissueTau.round = floor(tau$TissueTau*10)/10
		tau$TissueTau.round[tau$TissueTau.round==1] = 0.9
		tau$TissueTau.round = tau$TissueTau.round + 0.05
		boxplot(tau$TissueTau ~ tau$exon.no.char,frame=F,xlab='# of exons',ylab='Tissue Tau')
		if(lab) plotPanelLetter('A',lab.cex)
		z = table(tau$exon.no.char, tau$as)
		z = sweep(z,1,apply(z,1,sum),'/')
		barplot(t(z[,ord]),xlab='# of exons',ylab='gene fraction',col=cols,border=NA,main=paste0('All genes (N=',nrow(tau),')'))
		if(lab) plotPanelLetter('B',lab.cex)
		
		tau$rtau = residuals(lm(tau$TissueTau ~ tau$exon.no.char))
		tau$rtau  = floor(tau$rtau *10)/10 + 0.05
		barplot(table(tau$as,tau$rtau)[ord,],col=cols,border=NA,xlab='normalized TissueTau',ylab='# of genes')
		legend(6,2800,fill=c('orange','lightgray','darkgray'),legend=c('devAS','AS','non-AS'),bty='n')
		if(lab) plotPanelLetter('C',lab.cex)
		
		
		#barplot(table(tau$as,tau$TissueTau.round)[ord,],col=cols,border=NA,xlab='TissueTau',ylab='# of genes')
		# f = tau$exon.no == 4
		# barplot(table(tau$as[f],tau$TissueTau.round[f])[ord,],col=cols,border=NA,xlab='TissueTau',ylab='# of genes',main=paste0('Only genes with 4 exons (N=',sum(f),')'))
		# 
		# f = tau$exon.no == 10
		# barplot(table(tau$as[f],tau$TissueTau.round[f])[ord,],col=cols,border=NA,xlab='TissueTau',ylab='# of genes',main=paste0('Only genes with 10 exons (N=',sum(f),')'))
		lab = F
	}
}
dev.off()



pdf('figures/paper.figures/6/nature.review/reviewer.fig/R2.tau-no.exon.pdf',w=7,h=8,family='Arial')
#par(mfrow=c(4,6),tck=-0.01,mgp=c(1.4,0.4,0),mar=c(2.5,2.5,1,0),oma=c(0,0,0,1))
par(tck=-0.01,mgp=c(1.4,0.4,0),mar=c(2.5,2.5,1.3,0),oma=c(0,0,0,1))
layout(matrix(1:12,ncol=3,byrow = T),widths=c(1,3,3))
cols = c('darkgray','lightgray','orange')
ord = c('non-AS','AS','devAS')
taus = list(human=human.tau,mouse=mouse.tau)
dpsis = list(human=human.tissue.gene.dpsi,mouse=mouse.tissue.gene.dpsi)
leg = T
for(s in c('human','mouse')){
	for(t in c('brain','liver')){
		plot.new()
		par(mar=c(0,0,0,0),xpd=NA)
		plotPNG(paste0("figures/paper.figures/5/icons/",s,".png"),0.5,0.7,0.55)#,0.2,0.5,0.35)
		plotPNG(paste0("figures/paper.figures/5/icons/",t,".png"),0.5,0.3,0.55)#,0.8,0.5,0.35)
		par(mar=c(2.5,2.5,1.3,0),xpd=T)
		tau = taus[[s]][!is.na(taus[[s]]$TissueTau) & !is.na(taus[[s]]$biotype) & taus[[s]]$biotype =='protein_coding',]
		tau$as = ifelse(rownames(tau) %in% rownames(dpsis[[s]][[t]]),'AS','non-AS')
		tau$as[rownames(tau) %in% rownames(dpsis[[s]][[t]])[dpsis[[s]][[t]][,'up']>0.2]] = 'devAS'
		tau = tau[tau$exon.no>0,]
		tau$TissueTau.round = floor(tau$TissueTau*10)/10
		tau$TissueTau.round[tau$TissueTau.round==1] = 0.9
		tau$TissueTau.round = tau$TissueTau.round + 0.05

		f = tau$exon.no == 4
		barplot(table(tau$as[f],tau$TissueTau.round[f])[ord,],col=cols,border=NA,xlab='TissueTau',ylab='# of genes',main=paste0('Only genes with 4 exons (N=',sum(f),')'))

		f = tau$exon.no == 10
		barplot(table(tau$as[f],tau$TissueTau.round[f])[ord,],col=cols,border=NA,xlab='TissueTau',ylab='# of genes',main=paste0('Only genes with 10 exons (N=',sum(f),')'))
		if(leg){
			legend('topright',fill=c('orange','lightgray','darkgray'),legend=c('devAS','AS','non-AS'),bty='n')
			leg = F
		}
		
	}
}
dev.off()



pdf('figures/paper.figures/6/nature.review/tau-no.exon-ad.human.brain.pdf',w=9,h=9,family='Arial')
par(mfrow=c(3,3),tck=-0.01,mgp=c(1.4,0.4,0),mar=c(2.5,2.5,1,0),oma=c(0,0,0,1))
plot(exn.cnt,xlim=c(0,130),bty='n',xlab='# of exons',ylab='# of genes')
t$exon.no[t$exon.no>20] = 21
boxplot(t$TissueTau ~ t$exon.no,xlab='# of exons',ylab='TissueTau',frame=F)
boxplot(t$exon.no ~ round(t$TissueTau,1),ylab='# of exons',xlab='TissueTau',frame=F)

z=table( t$exon.no, t$as)
z = sweep(z,1,apply(z,1,sum),'/')
barplot(t(z[,ord]),xlab='# ofexons',ylab='gene fraction',col=cols,border=NA)

barplot(table(t$as,round(t$TissueTau,1))[ord,],col=cols,border=NA,xlab='TissueTau',ylab='# of genes',legend.text = T)

f = t$exon.no == 4
barplot(table(t$as[f],round(t$TissueTau[f],1))[ord,],col=cols,border=NA,legend.text = T,xlab='TissueTau',ylab='# of genes',main='Only genes with 4 exons')
f = t$exon.no == 8
barplot(table(t$as[f],round(t$TissueTau[f],1))[ord,],col=cols,border=NA,legend.text = T,xlab='TissueTau',ylab='# of genes',main='Only genes with 8 exons')
f = t$exon.no == 10
barplot(table(t$as[f],round(t$TissueTau[f],1))[ord,],col=cols,border=NA,legend.text = T,xlab='TissueTau',ylab='# of genes',main='Only genes with 10 exons')
f = t$exon.no == 21
barplot(table(t$as[f],round(t$TissueTau[f],1))[ord,],col=cols,border=NA,legend.text = T,xlab='TissueTau',ylab='# of genes',main='Only genes with >20 exons')
dev.off()

# _hexamer exon downsampling ######
hex2mot = read.table('output/hex2mot2sf.tab.gz')

hex.dws.age02sgn = readRDS('Rdata/hex.dws.age02sgn.Rdata')
hex.ups.age02sgn = readRDS('Rdata/hex.ups.age02sgn.Rdata')
hex.dws.age02sgn.ds = readRDS('Rdata/hex.dws.age02sgn.downsumpling.Rdata')
hex.ups.age02sgn.ds = readRDS('Rdata/hex.ups.age02sgn.downsumpling.Rdata')
hex.dws.age02sgn.dspt = readRDS('Rdata/hex.dws.age02sgn.downsampling.per.tissue.Rdata')
hex.ups.age02sgn.dspt = readRDS('Rdata/hex.ups.age02sgn.downsampling.per.tissue.Rdata')


hist(log2(hex.ups.age02sgn$up$or[,'testis',]))
hist(log2(hex.ups.age02sgn$up$or[,'brain',]))
wilcox.test(hex.ups.age02sgn$up$or[,'testis',],hex.ups.age02sgn$up$or[,'brain',],alternative = 'l')
apply(hex.ups.age02sgn$up$or[,'testis',]>1.5,2,sum)


hex.all  = getHexStat(hex.ups.age02sgn     ,hex.dws.age02sgn,0.05)
hex.dspt = getHexStat(hex.ups.age02sgn.dspt,hex.dws.age02sgn.dspt,0.05)
hex.ds   = getHexStat(hex.ups.age02sgn.ds  ,hex.dws.age02sgn.ds,0.05)
#hex.ds02 = getHexStat(hex.ups.age02sgn.ds,hex.dws.age02sgn.ds,0.2)
#hex.ds03 = getHexStat(hex.ups.age02sgn.ds,hex.dws.age02sgn.ds,0.3)

#pdf('figures/paper.figures/6/nature.review/hexamer.segment.downsampling.pdf',w=6,h=9,family='Arial')
#par(mfrow=c(3,2),mar=c(4.5,2.5,1.5,1.2),tck=-0.01,mgp=c(1.4,0.4,0))
pdf('figures/paper.figures/6/suppl/02NG/S17.hexamer.segment.downsampling.pdf',w=6,h=3.5,family='Arial')
par(mfrow=c(1,2),mar=c(4,2.5,1.5,1.2),tck=-0.01,mgp=c(1.4,0.4,0))
#plotHexStat(hex.all$hex.stat,hex.all$hex.stat.known,c('a','b'))
#mtext('All data, pv < 0.05',4)
plotHexStat(hex.dspt$hex.stat,hex.dspt$hex.stat.known,c('a','b'),plot.leg=T,leg.h=330)
#mtext('Seg downsample per tissue, pv < 0.05',4)
# plotHexStat(hex.ds$hex.stat,hex.ds$hex.stat.known,c('e','f'))
# mtext('Seg downsample all tissues, pv < 0.05',4)
dev.off()

plotHexStat(hex.ds03$hex.stat,hex.ds03$hex.stat.known,c('c','d'))

cnts = c(sum(apply(cbind(hex.dws.age02sgn$up$ih.qv,hex.ups.age02sgn$up$ih.qv)<0.05,1,sum)>0),sum(apply(cbind(hex.dws.age02sgn$dw$ih.qv,hex.ups.age02sgn$dw$ih.qv)<0.05,1,sum)>0))
prop.test(cnts,c(4^6,4^6))$p.value

cnts = c(sum(apply(cbind(hex.dws.age02sgn.ds$up$ih.qv,hex.ups.age02sgn.ds$up$ih.qv)<0.05,1,sum)>0),sum(apply(cbind(hex.dws.age02sgn.ds$dw$ih.qv,hex.ups.age02sgn.ds$dw$ih.qv)<0.05,1,sum)>0))
prop.test(cnts,c(4^6,4^6))$p.value


# _QKI expr and splicing ####
# gene.descrs$human[gene.descrs$human[,1]=='QKI',]# ENSG00000112531
# orth.ens.genes['ENSG00000112531',] # no one-to-one

o = read.csv('input/ens.orths.txt.gz')
qkio = o[o[,1]=='ENSG00000112531',]
gene.descrs$rat[qkio$Rat.Ensembl.Gene.ID,]

e2s = revList(seg2ens$human)
anns$human[rownames(anns$human) %in% e2s[['ENSG00000112531']],]

orth.qki.seg = orth.ads.all.sp[orth.ads.all.sp[,'human'] %in% e2s[['ENSG00000112531']],]
orth.qki.seg = orth.qki.seg[apply(is.na(orth.qki.seg),1,sum)==0,]
qkio
qkio = sapply(rownames(species),function(s)unique(unlist(seg2ens[[s]][orth.qki.seg[,s]]))) #agree with Ens orths (except rat paralog ENSRNOG00000017174)
par(mfrow=c(3,3),tck=-0.01,mgp=c(1.4,0.4,0),mar=c(2.5,2.5,1,0),oma=c(0,0,0,1))
for(s in rownames(species)){
	#plotTissueAgeProile(ens.ge.marg.tsm[[s]][qkio[s],],meta.tsm,age.axis = 'rank',main=s,ylab='CPM')
	plotTissueAgeProile(ens.ge.cod[[s]]$rpkm[qkio[s],],meta,age.axis = 'rank',main=s,ylab='RPKM',bty='n')
}


par(mfrow=c(7,1),tck=-0.01,mgp=c(1.4,0.4,0),mar=c(1,0,1,0),oma=c(0,0,0,1))
for(s in rownames(species)){
	t = all.anns[[s]][all.anns[[s]]$gene_id== all.anns[[s]][orth.qki.seg[1,s],'gene_id'],]
	plotIntronExonStructure(t,bty='n',main=s)
}



hens = loadSAData('processed/annotation/all.species/ensambl/Homo_sapiens.GRCh37.73.sajr.gz')
hens = setSplSiteTypes(hens,'processed/annotation/all.species/ensambl/Homo_sapiens.GRCh37.73.sajr.gz')

e = hens$seg[hens$seg$gene_id=='ENSG00000112531',]

s='human'
t = all.anns[[s]][all.anns[[s]]$gene_id== all.anns[[s]][orth.qki.seg[1,s],'gene_id'],]
xlim=c(min(e$start,t$start),max(e$stop,t$stop))
xlim=c(163980000,164000000)
plotIntronExonStructure(e,bty='n',xlim=xlim,main='Ensabmle')
plotIntronExonStructure(t,bty='n',xlim=xlim,main='Mine')
t[t$stop>=163980000 & t$start<=163990000 & t$sites=='ad',]
abline(v=c(163984476,163987753))
# so QKI6-7 are between hum.57513.s12 and hum.57513.s17
qki67.segs=lapply(rownames(species),function(s){
	t = orth.qki.seg[c('hum.57513.s12','hum.57513.s17'),s]
	t = all.anns[[s]][t,]
	r = all.anns[[s]][all.anns[[s]]$gene_id == t$gene_id[1] & all.anns[[s]]$start>= min(t$start) & all.anns[[s]]$stop<= max(t$stop),]
	r[order(r$strand*r$start),]
	})
names(qki67.segs) = rownames(species)
lapply(qki67.segs,function(x)x[,c('start','stop','type','length','sites')])
qki7 = c('hum.57513.s14','mac.56888.s11','mou.19182.s23','rat.718.s30','rab.6709.s17','opo.18986.s13','chi.19164.s9') # not sure about opossum (and chicken)
qki67= c('hum.57513.s16','mac.56888.s13','mou.19182.s22','rat.718.s31','rab.6709.s18','opo.18986.s14','chi.19164.s8') #or hum.57513.s15
#qki67.v2=c('hum.57513.s15','mac.56888.s13')
names(qki7) = names(qki67) = rownames(species)

pdf('figures/paper.figures/6/nature.review/QKI.expression-n-splicing.pdf',w=7,h=15,family='Arial')
par(mfrow=c(7,3),tck=-0.01,mgp=c(1.4,0.4,0),mar=c(2.5,2.5,1,0),oma=c(0,0,0,1))
for(s in rownames(species)){
	plotTissueAgeProile(ens.ge.cod[[s]]$rpkm[qkio[s],],meta,age.axis = 'rank',main=s,ylab='RPKM',bty='n')
	plotTissueAgeProile(psi.tsm[[s]][qki7[s],],meta.tsm,age.axis = 'rank',main='QKI7',ylab='PSI',bty='n')
	plotTissueAgeProile(psi.tsm[[s]][qki67[s],],meta.tsm,age.axis = 'rank',main='QKI 6+7',ylab='PSI',bty='n')
}
dev.off()

qki67.segs$mouse
plotTissueAgeProile(psi.tsm$mouse['mou.19182.s21',],meta.tsm,age.axis = 'rank',main='mou.19182.s21')
plotTissueAgeProile(psi.tsm$mouse['mou.19182.s22',],meta.tsm,age.axis = 'rank',main='mou.19182.s22')
plotTissueAgeProile(psi.tsm$mouse['mou.19182.s23',],meta.tsm,age.axis = 'rank',main='mou.19182.s23')
plotTissueAgeProile(psi.tsm$mouse['mou.19182.s24',],meta.tsm,age.axis = 'rank',main='mou.19182.s24')

# _radical changes ######
age.al.i
meta[1:2,]

patt.comp = NULL
sps = c(1,3:7)
for(i1 in 1:(length(sps)-1))
	for(i2 in (i1+1):length(sps)){
		a = age.al.i
		s1 = rownames(species)[sps[i1]]
		s2 = rownames(species)[sps[i2]]
		
		if(s1 == 'opossum' | s2 == 'opossum')
			a = a[a$to.remove==0,]
		if(s1 == 'chicken' | s2 == 'chicken')
			a = a[!(a$mouse %in% c('17.5','18.5')),]
		for(t in unique(meta$tissue)){
			cat('\r',s1,s2,t,'        ')
			r = getSpeciesCorForTissue(s1,s2,t,a)
			r = r[r$p1 != '-' & r$p2 != '-',]
			patt.comp = rbind(patt.comp,r)
		}
	}
table(patt.comp$species1,patt.comp$species2)[rownames(species)[sps[-6]],rownames(species)[sps[-1]]]

s2e = unlist(setNames(seg2ens,NULL),recursive = F)
e1 = s2e[patt.comp$sid1]
e2 = s2e[patt.comp$sid2]
table(sapply(e1,length))

orth.ens.genes[1:2,]
g2h = unlist(lapply(1:ncol(orth.ens.genes),function(i)setNames(orth.ens.genes[,1],orth.ens.genes[,i])))
g2h[1:2]

patt.comp$ens.orth = F
for(i in 1:nrow(patt.comp)){
	g1 = g2h[e1[[i]]]
	g2 = g2h[e2[[i]]]
	patt.comp$ens.orth[i] = sum(g1 %in% g2) > 0
}
table(patt.comp$ens.orth)

#mrb = getSpeciesCor('mouse','opossum','brain',age.al.i[age.al.i$to.remove==0,])
dim(patt.comp)

o = c('n','u','d','ud','du')
table(patt.comp$p1,patt.comp$p2)[o,o]
f = patt.comp$p1 %in% o[-1] & patt.comp$p2 %in% o[-1]
hist(patt.comp$cor[f])
hist(patt.comp$diff[f])
plotLine(patt.comp$diff[f],patt.comp$cor[f],pch='.')
#saveRDS(patt.comp,'Rdata/20201014.pattern.changes.Rdata')
patt.comp = readRDS('Rdata/20201014.pattern.changes.Rdata')

patt.comp.sel = patt.comp[patt.comp$p1 %in% c('u','d') & 
													patt.comp$p2 %in% c('u','d') & 
													patt.comp$p1 != patt.comp$p2 & 
													patt.comp$nnna > 5 &
													patt.comp$cor < -0.5 &
													patt.comp$ens.orth & 
													patt.comp$dpsi1>0.5 & 
													patt.comp$dpsi2>0.5,]
hist(patt.comp.sel$cor)
hist(patt.comp.sel$dpsi1)
hist(patt.comp.sel$dpsi2)

table(patt.comp.sel$tissue)
table(patt.comp.sel$species1,patt.comp.sel$species2)
dim(patt.comp.sel)

pdf('figures/paper.figures/6/nature.review/pattern.changes.pdf',w=21,h=15,family='Arial')
par(mfrow=c(5,7),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(2.5,2.5,3.5,0.5),oma=c(0,0,0,1))
for(i in 1:nrow(patt.comp.sel)){
	osids = match(patt.comp.sel$sid1[i],rownames(orth.seg.ad.tsm[[patt.comp.sel$species1[i]]]))
	osids = sapply(orth.seg.ad.tsm,function(x)rownames(x)[osids])
	for(s in rownames(species)){
		d1 = gene.descrs[[s]][seg2ens[[s]][[osids[s]]],]
		plotTissueAgeProile(orth.seg.ad.tsm[[s]][osids[s],],m=meta.tsm,col=paste0(meta.tsm$col,'50'),ylim = c(0,1),age.axis = 'rank',ylab='PSI',tissues = setdiff(unique(meta$tissue),patt.comp.sel$tissue[i]),bty='n',main=paste0(rownames(d1)[1],"\n",d1[1,1],' ',d1[1,2],'\n',s))
		plotTissueAgeProile(orth.seg.ad.tsm[[s]][osids[s],],m=meta.tsm,lwd=4,ylim = c(0,1),age.axis = 'rank',ylab='PSI',tissues = patt.comp.sel$tissue[i],bty='n',add = TRUE,plot.xaxt = FALSE)
	}
}
dev.off()

par(mfrow=c(2,2),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(2.5,2.5,1,0.5),oma=c(0,0,0,1))
plotTissueAgeProile(psi.tsm$human['hum.36710.s45',],meta.tsm,age.axis = 'rank',ylim=c(0,1),bty='n')
plotTissueAgeProile(psi.tsm$rat['rat.25212.s37',],meta.tsm,age.axis = 'rank',ylim=c(0,1),bty='n')
plotTissueAgeProile(psi.tsm$human['hum.16430.s12',],meta.tsm,age.axis = 'rank',ylim=c(0,1),bty='n')
plotTissueAgeProile(psi.tsm$rat['rat.9251.s10',],meta.tsm,age.axis = 'rank',ylim=c(0,1),bty='n')
plotTissueAgeProile(psi.tsm$human['hum.16430.s12',],meta.tsm,age.axis = 'rank',ylim=c(0,1),bty='n')
plotTissueAgeProile(psi.tsm$rat['rat.9251.s10',],meta.tsm,age.axis = 'rank',ylim=c(0,1),bty='n')






table(is.na(mrb$dpsi1),mrb$p1)
mrb[mrb$p1 %in% o[-1:-2] & mrb$p2 %in% o[-1:-2] & mrb$cor< -0.1 & mrb$nnna>10,]
plotTissueAgeProile(psi.tsm$mouse['mou.39205.s18',],meta.tsm,age.axis = 'rank',ylim=c(0,1),bty='n')
plotTissueAgeProile(psi.tsm$opossum['opo.49365.s42',],meta.tsm,age.axis = 'rank',ylim=c(0,1),bty='n')
anns$opossum['opo.49365.s42',]
anns$mouse['mou.39205.s18',]
gene.descrs$opossum[seg2ens$opossum[['opo.49365.s42']],]
gene.descrs$mouse[seg2ens$mouse[['mou.39205.s18']],]
orth.ens.genes['ENSG00000133958',]



getCommonStages = function(s1,s2,tissues,age.al){
	age.al = age.al[,c(s1,s2)]
	age.al = age.al[apply(age.al=='',1,sum)==0,]
	as1 = paste(s1,rep(tissues,each=nrow(age.al)),rep(age.al[,s1],times=length(tissues)))
	as2 = paste(s2,rep(tissues,each=nrow(age.al)),rep(age.al[,s2],times=length(tissues)))
	f = as1 %in% rownames(meta.tsm) & as2 %in% rownames(meta.tsm)
	as1 = as1[f]
	as2 = as2[f]
	r = data.frame(as1,as2)
	colnames(r) = c(s1,s2)
	r
}

hm = getCommonStages('opossum','mouse',unique(meta$tissue),age.al.i[age.al.i$to.remove==0,])
hmc = sapply(1:nrow(orth.seg.ad.tsm$opossum),function(i)cor(orth.seg.ad.tsm$opossum[i,hm$opossum],orth.seg.ad.tsm$mouse[i,hm$mouse],u='p'))
hmn = apply(!is.na(orth.seg.ad.tsm$opossum[,hm$opossum]) & !is.na(orth.seg.ad.tsm$mouse[,hm$mouse]),1,sum)
hist(hmn)
i=which(hmn > 20 & hmc < -0.5)[4]
i
rownames(orth.seg.ad.all.tsm$mouse)[i]
m = meta.tsm[hm$mouse,]
plotLine(orth.seg.ad.tsm$opossum[i,hm$opossum],orth.seg.ad.tsm$mouse[i,hm$mouse],bty='n',pch=19,col=m$col)
gene.descrs$opossum[seg2ens$opossum[['opo.49874.s28']],]
gene.descrs$mouse[seg2ens$mouse[['mou.42866.s7']],]

orth.ens.genes[orth.ens.genes$mouse=='ENSMUSG00000002949',]
orth.ens.genes[orth.ens.genes$opossum=='ENSMODG00000018544',]
