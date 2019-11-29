#setwd('~/skoltech/projects/evo.devo/')
options(stringsAsFactors = FALSE)
source('code/r.functions/load.all.data.F.R')
source('code/r.functions/paper.figures.F.R')
source('code/r.functions/ad.on.ge.F.R')
library(SAJR)
library(ape)
library(GO.db)
library(GenomicRanges)

species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
orth.seg.ad.tsm = readRDS('Rdata/orth.seg.ad.tsm.Rdata')
orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
#orth.seg.ad.all = readRDS('Rdata/orth.seg.ad.all.Rdata')
orth.seg.ad.all.id = readRDS('Rdata/orth.seg.ad.all.id.Rdata')
ens.ge.cod.tsm = readRDS('Rdata/ens.ge.cod.tsm.Rdata')
orth.ens.genes = readRDS('Rdata/orth.ens.genes.Rdata')
psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
anns = readRDS('Rdata/anns.Rdata')
#all.anns = readRDS('Rdata/all.anns.Rdata')
border.stages = readRDS('Rdata/border.stages.Rdata')
seg2ens = readRDS('Rdata/seg2ens.Rdata')
ens.ge.cod = readRDS('Rdata/ens.ge.cod.Rdata')

orth.per.tissue.age.qv = readRDS('Rdata/orth.per.tissue.age.qv.Rdata')
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')
per.tissue.age.amp = readRDS('Rdata/per.tissue.age.amp.Rdata')
per.tissue.age.sgn = lapply(1:7,function(s){per.tissue.age.qv[[s]]<0.05 & per.tissue.age.amp[[s]]>0.2})
names(per.tissue.age.sgn) = names(per.tissue.age.amp)

hex.dws.age03 = readRDS('Rdata/hex.dws.age03.Rdata')
hex.ups.age03 = readRDS('Rdata/hex.ups.age03.Rdata')

age.al.i = read.table('input/margarida.age.alignment.final.2017.05.17.improved.tab',sep='\t',header = T)[,c(rownames(species),'to.remove')]
params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)
DPSI = 0.5

ens.descr.mm = unique(read.table('input/mm.38.84.gene.descr.txt',sep=',',quote='"',header=TRUE))
rownames(ens.descr.mm) = ens.descr.mm$Ensembl.Gene.ID
ens.descr.mm$Description = sapply(strsplit(ens.descr.mm$Description,' [',TRUE),'[',1)

########################
### prepare data #######
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

mouse.go = read.table('input/GO/Mus_musculus.GRCm38.84.GO.csv.gz',sep=',',header = T)
mouse.go = mouse.go[mouse.go[,2] != '',]
mouse.go = split(mouse.go$GO.Term.Accession,mouse.go$Ensembl.Gene.ID)
# GOALLANCESTOR = c(as.list(GOCCANCESTOR),as.list(GOMFANCESTOR),as.list(GOBPANCESTOR))
# mouse.go.full = lapply(mouse.go,function(x){r=unique(c(x,unlist(GOALLANCESTOR[x])));r[grep('GO:',r,fixed=T)]})
# mouse.go.full.rev = revList(mouse.go.full)
mouse.go.rev = revList(mouse.go)

human.go = read.table('input/GO/Homo_sapiens.GRCh37.74.GO.csv.gz',sep=',',header = T)
human.go = human.go[human.go[,2] != '',]
human.go = split(human.go$GO.Term.Accession,human.go$Ensembl.Gene.ID)
human.go.rev = revList(human.go)


# I DO NOT USE NON-CODING SEGMENTS
as.in.ge.patterns.human = list()
sp = 'human'
mc = read.csv(paste0('processed/GE.from.marg/',firstToupper(sp),'Clusters.csv'),row.names = 1)
colnames(mc) = tolower(colnames(mc))
for(tis in unique(meta$tissue)){
	cat(tis)
	t = getsPSIbyEnsID(list(human=psi.tsm[[sp]][anns[[sp]]$sites=='ad' & anns[[sp]]$cod!='n',]),border.stages,tis,seg2ens,sp)
	t = t(t)
	t = t[intersect(rownames(t),rownames(mc)[!is.na(mc[,paste0(tis,'pattern')])]),]
	as.in.ge.patterns.human[[tis]] = cbind(data.frame(t),ge.pattern=mc[rownames(t),paste0(tis,'pattern')])
}
saveRDS(as.in.ge.patterns.human,'Rdata/paper.figures/as.in.ge.patterns.human.Rdata')



##################
### figure 1 #####
orth.mds2 = readRDS('Rdata/paper.figures/orth.mds2.Rdata')

filtered.seg.cnt = sapply(anns,getNoOfEvents,gene=FALSE)
age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)

sign.stat = list(all=sapply(names(per.tissue.age.sgn),function(s)apply(per.tissue.age.sgn[[s]][anns[[s]]$sites=='ad',],2,sum,na.rm=T)),
								 up=sapply(names(per.tissue.age.sgn)[-2],function(s)apply((age.dpsi[[s]] > 0 & per.tissue.age.sgn[[s]])[anns[[s]]$sites=='ad',],2,sum,na.rm=T)),
								 down=sapply(names(per.tissue.age.sgn)[-2],function(s)apply((age.dpsi[[s]] < 0 & per.tissue.age.sgn[[s]])[anns[[s]]$sites=='ad',],2,sum,na.rm=T)))

p = cbind(orth.seg.ad$mouse$ir,orth.seg.ad$rat$ir)
a = orth.seg.ad$mouse$seg$type=='ALT' & orth.seg.ad$rat$seg$type=='ALT'
d = cbind(apply(orth.per.tissue.age.qv$mouse<0.05,1,sum,na.rm=T)>0,apply(orth.per.tissue.age.qv$rat<0.05,1,sum,na.rm=T)>0)
mr.cors = list(dev =cor(p[a & apply(d,1,sum)>0,],u='p',m='sp'),
							 ndev=cor(p[a & apply(d,1,sum)==0,],u='p',m='sp'))

pdf('figures/paper.figures/3.2/1.pdf',w=9,h=6)
#layout(matrix(c(3,1,3,2),ncol=2))
layout(matrix(c(1,1,2:5),byrow = TRUE,ncol=3))
par(tck=-0.01,mgp=c(1.9,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,0,1))
plot.new()
text(0.5,0.5,'study design\nspecies x tissue x time',cex=3)
plotPanelLetter('A')
#col=rep(c('red','orange','magenta','blue'),each=3)
den=rep(c(-1,60,20))#,times=4)
barplot(filtered.seg.cnt[1:3,],den=den,ylab='# of exons',las=3,main='Cassette exons passed filtering',legend.text=c('cod.','partially cod.','non cod.'))#,col=col)
#legend('topright',fill=c('red','orange','magenta','blue','black','black','black'),den=c(-1,-1,-1,-1,-1,den[-1]),legend=c('cassette','alt. acc.','alt. don.','int. ret.','cod.','partially cod.','non cod.'),bty='n')
plotPanelLetter('B')
m = meta
SAJR::plotMDS(points=-orth.mds2$all,col=m$col,cex=m$cex*1.5,pch=m$pch,main=paste0('All cassette exons (',nrow(orth.seg.ad.tsm$human),')'))
legend('topleft',col=params$tissue.col,pch=19,legend=names(params$tissue.col),bty='n')
plotPanelLetter('C')

b=barplot(t(sign.stat$all),beside = T,col=rep(params$tissue.col,each=7),ylab='Number of devAS',main='FDR < 0.05 & dPSI > 0.2')
text(b,0,species$short,adj = c(0.5,1),xpd=T,cex=0.9)
plotPanelLetter('D')
plotSpCor2Cor(mr.cors$dev,mr.cors$ndev,age.al.i[,c('mouse','rat')],meta,main='Mouse-rat Spearman correlation',xlab='devAS (2243)',ylab='not devAS (3112)')

plotPanelLetter('E')
dev.off()



#########################
### figure 2 prepare ####
# library(ape)
# library(TKF)
# mdata = readRDS('Rdata/mouse.as.u.filtered.Rdata')
# mouse.ad.cor = cor(mdata$ir[mdata$seg$sites=='ad',],u='p')
# mouse.ad.cor.mds = cmdscale(1-mouse.ad.cor,k=2)
# mouse.ad.tsm.mds = cmdscale(1-cor(psi.tsm$mouse[anns$mouse$sites=='ad',],u='p'),k=2)
# saveRDS(mouse.ad.tsm.mds,'Rdata/paper.figures/mouse.ad.tsm.mds.Rdata')
# saveRDS(mouse.ad.cor.mds,'Rdata/paper.figures/mouse.ad.cor.mds.Rdata')
#splicing correlation
# p = psi.tsm$mouse[anns$mouse$sites=='ad',]
# d = sort(unique(meta$days[meta$species=='mouse']))
# m = meta.tsm[colnames(p),]
# tissues = unique(meta$tissue)
# base = apply(p[,m$days==d[1]],1,mean,na.rm=T)
# mouse.stage.cors = lapply(d,function(day){
# 	r = p[,m$days==day]
# 	colnames(r) = substr(m$tissue[m$days==day],1,1)
# 	r = cbind(base=base,r)
# 	cor(r,u='p')})
# 
# #gene expression correlation
# mouse.ge = readRDS('Rdata/ens.ge.cod.tsm.Rdata')$mouse
# gc()
# mouse.ge = mouse.ge[,rownames(m)]
# mouse.ge = mouse.ge + min(mouse.ge[mouse.ge!=0],na.rm=T)
# mouse.ge =log2(mouse.ge)
# base = apply(mouse.ge[,m$days==d[1]],1,mean,na.rm=T)
# mouse.ge.stage.cors = lapply(d,function(day){
# 	r = mouse.ge[,m$days==day]
# 	colnames(r) = substr(m$tissue[m$days==day],1,1)
# 	r = cbind(base=base,r)
# 	cor(r,u='p')})
# 
# sfs = rownames(mouse.ge) %in% mouse.go.rev[['GO:0008380']]
# 
# 
# mouse.sf.stage.cors = lapply(d,function(day){
# 	r = mouse.ge[sfs,m$days==day]
# 	colnames(r) = substr(m$tissue[m$days==day],1,1)
# 	r = cbind(base=base[sfs],r)
# 	cor(r,u='p')})
# names(mouse.stage.cors) = names(mouse.ge.stage.cors) = names(mouse.sf.stage.cors) = d

# mouse.stage.trees=list()
# mouse.stage.trees$as = makeNJTreesByCorLastTopology(mouse.stage.cors)
# mouse.stage.trees$ge = makeNJTreesByCorLastTopology(mouse.ge.stage.cors)
# mouse.stage.trees$sf = makeNJTreesByCorLastTopology(mouse.sf.stage.cors)
# saveRDS(mouse.stage.trees,'Rdata/paper.figures/mouse.stage.trees.Rdata')

# calc divergence to embrio
mouse.ge = readRDS('Rdata/ens.ge.cod.tsm.Rdata')$mouse
gc()

mouse.ge = mouse.ge + min(mouse.ge[mouse.ge!=0],na.rm=T)
mouse.ge =log2(mouse.ge)

psi.cor2embryo = caclCor2Embryo(psi.tsm$mouse[anns$mouse$sites=='ad',],meta.tsm,cor.m = 'p')
ge.cor2embryo = caclCor2Embryo(mouse.ge,meta.tsm,cor.m = 'p')
sf.cor2embryo = caclCor2Embryo(mouse.ge[rownames(mouse.ge) %in% unlist(mouse.go.rev[c('GO:0008380')]),],meta.tsm,cor.m = 'p')

# load MDSs
mouse.cassett.cnt = 20186 #sum(mdata$seg$sites=='ad')
ancient.cnt = 1590 #sum(apply(sapply(orth.seg.ad,function(x)x$seg$type=='ALT'),1,sum)==7)
mouse.ad.cor.mds = readRDS('Rdata/paper.figures/mouse.ad.cor.mds.Rdata')
#mouse.stage.trees = readRDS('Rdata/paper.figures/mouse.stage.trees.Rdata')
m = getAgeASchanges(psi.tsm,meta.tsm,DPSI,border.stages,'mouse')
mm = apply(m,2,function(x){x %in% c('u','d')})
#tr.compl.m =getAltExonStat(psi.tsm$mouse[anns$mouse$sites=='ad' & apply(mm,1,sum)>0,],meta.tsm,0.1,tissues = c('brain','heart','liver','ovary','testis'),na.as.cnst=TRUE)

as.in.ge.patterns.mouse = readRDS('Rdata/paper.figures/as.in.ge.patterns.mouse.Rdata')
as.in.ge.patterns.mouse.cnt.up = lapply(as.in.ge.patterns.mouse,function(x){table(factor(x[,3]),x[,1]>DPSI)[,c('TRUE','FALSE')]})
as.in.ge.patterns.mouse.stat.up = sapply(as.in.ge.patterns.mouse.cnt.up,function(x){cbind(my.binom.test(x['Increasing',]),my.binom.test(x['Decreasing',]))})
as.in.ge.patterns.mouse.cnt.dw = lapply(as.in.ge.patterns.mouse,function(x){table(factor(x[,3]),x[,2]< -DPSI)[,c('TRUE','FALSE')]})
as.in.ge.patterns.mouse.stat.dw = sapply(as.in.ge.patterns.mouse.cnt.dw,function(x){cbind(my.binom.test(x['Increasing',]),my.binom.test(x['Decreasing',]))})


as.in.ge.patterns.human = readRDS('Rdata/paper.figures/as.in.ge.patterns.human.Rdata')
as.in.ge.patterns.human.cnt.up = lapply(as.in.ge.patterns.human,function(x){table(factor(x[,3]),x[,1]>DPSI)[,c('TRUE','FALSE')]})
as.in.ge.patterns.human.stat.up = sapply(as.in.ge.patterns.human.cnt.up,function(x){cbind(my.binom.test(x['Increasing',]),my.binom.test(x['Decreasing',]))})
as.in.ge.patterns.human.cnt.dw = lapply(as.in.ge.patterns.human,function(x){table(factor(x[,3]),x[,2]< -DPSI)[,c('TRUE','FALSE')]})
as.in.ge.patterns.human.stat.dw = sapply(as.in.ge.patterns.human.cnt.dw,function(x){cbind(my.binom.test(x['Increasing',]),my.binom.test(x['Decreasing',]))})


#Evolution of GE
#mouse
ge.info.m = read.csv('/home/mazin/skoltech/projects/evo.devo/input/gene.info.from.marg/Mouse.Indexes.csv')
rownames(ge.info.m) = ge.info.m$Mouse_ID
m.age.ens = getAgeASEns(psi.tsm,meta.tsm,DPSI,border.stages,'mouse')
ohnologs.m = read.table('input/ohnologs/MOUSE.Pairs.Intermediate.2R.txt',sep='\t',header = T)

table(rec.dupl=ge.info.m$Age>0,ohn=ge.info.m$Mouse_ID %in% c(ohnologs.m$Ohnolog.1.Id,ohnologs.m$Ohnolog.2.Id))
#human
ge.info.h = read.csv('/home/mazin/skoltech/projects/evo.devo/input/gene.info.from.marg/Human.Indexes.All.csv')
rownames(ge.info.h) = ge.info.h$Human_ID
h.age.ens = getAgeASEns(psi.tsm,meta.tsm,DPSI,border.stages,'human')
ohnologs.h = read.table('input/ohnologs/HUMAN.Pairs.Intermediate.2R.txt',sep='\t',header = T)


#### panAS/ageAS, TF and SF
# mouse
m = getAgeASchanges(psi.tsm,meta.tsm,DPSI,border.stages,'mouse')
tested.ens.m = unique(unlist(seg2ens$mouse[rownames(anns$mouse)[apply(m != '-',1,sum)>0 & anns$mouse$sites=='ad']]))
tested.ens.m = data.frame(ageAS = tested.ens.m %in% unlist(seg2ens$mouse[rownames(anns$mouse)[apply(m == 'u' | m == 'd',1,sum)>0 & anns$mouse$sites=='ad']]),row.names = tested.ens.m)
#panAS >80% stages in [0.1,0.9]
panas=apply(psi.tsm$mouse,1,function(x)sum(x>0.1 & x<0.9,na.rm=TRUE)/length(x))

tested.ens.m$panAS = rownames(tested.ens.m) %in% unlist(seg2ens$mouse[rownames(anns$mouse)[panas > 0.8 & anns$mouse$sites=='ad']])
tested.ens.m$TF = rownames(tested.ens.m) %in% ge.info.m$Mouse_ID[ge.info.m$TF=='Mouse_TF']
#GO:0008380 = RNA splicing, GO:0003723 - RNA binding
tested.ens.m$SF = rownames(tested.ens.m) %in% mouse.go.rev[['GO:0008380']]
tested.ens.m$`RNA-binding` = rownames(tested.ens.m) %in% mouse.go.rev[['GO:0003723']]

inageAS.m=apply(tested.ens.m[,c(3,4,5,2)],2,function(x){t = fisher.test(table(tested.ens.m$ageAS,x));c(t$estimate,t$p.value,t$conf.int)})
inpanAS.m=apply(tested.ens.m[,c(3,4,5,1)],2,function(x){t = fisher.test(table(tested.ens.m$panAS,x));c(t$estimate,t$p.value,t$conf.int)})

# human
h = getAgeASchanges(psi.tsm,meta.tsm,DPSI,border.stages,'human')
tested.ens.h = unique(unlist(seg2ens$human[rownames(anns$human)[apply(h != '-',1,sum)>0 & anns$human$sites=='ad']]))
tested.ens.h = data.frame(ageAS = tested.ens.h %in% unlist(seg2ens$human[rownames(anns$human)[apply(h == 'u' | h == 'd',1,sum)>0 & anns$human$sites=='ad']]),row.names = tested.ens.h)
#panAS >80% stages in [0.1,0.9]
panas=apply(psi.tsm$human,1,function(x)sum(x>0.1 & x<0.9,na.rm=TRUE)/length(x))

tested.ens.h$panAS = rownames(tested.ens.h) %in% unlist(seg2ens$human[rownames(anns$human)[panas > 0.8 & anns$human$sites=='ad']])
tested.ens.h$TF = rownames(tested.ens.h) %in% read.table('input/Homo_sapiens_transcription_factors_gene_list.AnimalTFDB.txt',sep='\t')[,1]
#GO:0008380 = RNA splicing, GO:0003723 - RNA binding
tested.ens.h$SF = rownames(tested.ens.h) %in% human.go.rev[['GO:0008380']]
tested.ens.h$`RNA-binding` = rownames(tested.ens.h) %in% human.go.rev[['GO:0003723']]

# tested.ens.h$TF = rownames(tested.ens.h) %in% ge.info.m$Human_ID[ge.info.m$TF=='Mouse_TF']
# #GO:0008380 = RNA splicing, GO:0003723 - RNA binding
# tested.ens.h$SF = rownames(tested.ens.h) %in% orth.ens.genes[orth.ens.genes[,3] %in% mouse.go.rev[['GO:0008380']],1]
# tested.ens.h$`RNA-binding` = rownames(tested.ens.h) %in% orth.ens.genes[orth.ens.genes[,3] %in%  mouse.go.rev[['GO:0003723']],1]

inageAS.h=apply(tested.ens.h[,c(3,4,5,2)],2,function(x){t = fisher.test(table(tested.ens.h$ageAS,x));c(t$estimate,t$p.value,t$conf.int)})
inpanAS.h=apply(tested.ens.h[,c(3,4,5,1)],2,function(x){t = fisher.test(table(tested.ens.h$panAS,x));c(t$estimate,t$p.value,t$conf.int)})


######################
### figure 2.2 plot ####

pdf('figures/paper.figures/3.2/2.2.pdf',w=12,h=12)
par(mfrow=c(4,4),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,2.5),oma=c(0,0,0,1))

m = meta.tsm[rownames(orth.mds2$sp7tsm),]
#SAJR::plotMDS(points=orth.mds2$sp7,col=m$col,cex=m$cex*1.5,pch=m$pch,main=paste0('All cassette exons (',ancient.cnt,')'))
SAJR::plotMDS(points=orth.mds2$sp7tsm,col=m$col,cex=m$cex*1.5,pch=m$pch,main=paste0('All cassette exons (',ancient.cnt,')'))
plotPanelLetter('A')
m = meta.tsm[rownames(mouse.ad.tsm.mds),]
#SAJR::plotMDS(points=mouse.ad.cor.mds,col=m$col,cex=m$cex*2,pch=19,main=paste0('All mouse cassette exons (',mouse.cassett.cnt,')'))
SAJR::plotMDS(points=mouse.ad.tsm.mds,col=m$col,cex=m$cex*2,pch=19,main=paste0('All mouse cassette exons (',mouse.cassett.cnt,')'))
plotPanelLetter('B')
plot.new()
legend('topleft',col=params$tissue.col,pch=19,legend=names(params$tissue.col),bty='n')
plot.new()


par(mar=c(5.5,2.5,1.5,0),mgp=c(1.3,0.2,0))
plotCor2Embryo(psi.cor2embryo,main='Splicing',ylab='Pearson cor. to embryo, PSI',lwd=3,area.transp=0.1,xlab='')
plotPanelLetter('C')
plotCor2Embryo(ge.cor2embryo,main='Gene expression',ylab='Pearson cor. to embryo, log(RPKM)',lwd=3,area.transp=0.1,xlab='')
plotPanelLetter('D')
plotCor2Embryo(sf.cor2embryo,main='SF gene expression',ylab='Pearson cor. to embryo, log(RPKM)',lwd=3,area.transp=0.1,xlab='')
plotPanelLetter('E')
plot.new()


c = rep(params$tissue.col[colnames(as.in.ge.patterns.mouse.stat.up)],each=2)
b = barplot(as.in.ge.patterns.mouse.stat.up[c(1,4),],beside = T,col=c,den=c(-1,30),las=3,ylim=c(-max(as.in.ge.patterns.mouse.stat.dw[c(3,6),]),max(as.in.ge.patterns.mouse.stat.up[c(3,6),])),ylab='proportion of genes with devAS',main='DevAS in mouse GE clusters',yaxt='n')
segments(b,as.in.ge.patterns.mouse.stat.up[c(2,5),],b,as.in.ge.patterns.mouse.stat.up[c(3,6),])
b = barplot(-as.in.ge.patterns.mouse.stat.dw[c(1,4),],beside = T,col=c,den=c(-1,30),add=T,xaxt='n',yaxt='n')
segments(b,-as.in.ge.patterns.mouse.stat.dw[c(2,5),],b,-as.in.ge.patterns.mouse.stat.dw[c(3,6),])
abline(h=0)
at=-1:3*0.05
axis(2,at,abs(at))
legend('topright',fill='black',den=c(30,-1),legend=c('Early genes','Late Genes'))
text(b[2],par('usr')[4],'AS up',adj=c(-0.1,1.1))
text(b[2],par('usr')[3],'AS down',adj=c(-0.1,-1.1))
plotPanelLetter('F')

#plotOhnologsFreq(m.age.ens,ge.info.m$Mouse_ID[!is.na(ge.info.m$Age) & ge.info.m$Age>3],main='DevAS are depleted in recent duplications',
#								 filter = intersect(ge.info.m$Mouse_ID[!is.na(ge.info.m$Age)],rownames(ens.ge.cod$mouse$gene)),ylab='Proportion of recent genes',legend=TRUE,merge.up.down=TRUE)
plotDevASFreq(m.age.ens,ge.info.m$Mouse_ID[!is.na(ge.info.m$Age) & ge.info.m$Age>3],'mam. dupl.','topleft',main='DevAS are depleted in recent duplications',
							filter = intersect(ge.info.m$Mouse_ID[!is.na(ge.info.m$Age)],rownames(ens.ge.cod$mouse$gene)))
plotPanelLetter('G')
#plotOhnologsFreq(m.age.ens,c(ohnologs.m$Ohnolog.1.Id,ohnologs.m$Ohnolog.2.Id),main='DevAS are enriched in ohnologs',
##								 filter = intersect(ge.info.m$Mouse_ID[!is.na(ge.info.m$Age)], #remove filter from here
#								 rownames(ens.ge.cod$mouse$gene)),merge.up.down=TRUE)
plotDevASFreq(m.age.ens,c(ohnologs.m$Ohnolog.1.Id,ohnologs.m$Ohnolog.2.Id),'ohnologs','topright',main='DevAS are enriched in ohnologs')
plotPanelLetter('H')

t=log2(rbind(inageAS.m[1,],inpanAS.m[1,]))
ci1 = log2(rbind(inageAS.m[3,],inpanAS.m[3,]))
ci2 = log2(rbind(inageAS.m[4,],inpanAS.m[4,]))
#t[2,4] = ci1[2,4] = ci1[2,4] = NA
b=barplot(t,beside = T,ylab='log2(odds ratio)',ylim=range(t,ci1,ci2,na.rm=T),names.arg = c('TFs','SFs','RNA-binding','Pan/DevAS'),legend.text = c('',''),args.legend = list(x='bottomright',title='Assotiation with:',legend=c('AgeAS','PanAS')))
segments(b,ci1,b,ci2)
plotPanelLetter('I')

c = rep(params$tissue.col[colnames(as.in.ge.patterns.human.stat.up)],each=2)
b = barplot(as.in.ge.patterns.human.stat.up[c(1,4),],beside = T,col=c,den=c(-1,30),las=3,ylim=c(-max(as.in.ge.patterns.human.stat.dw[c(3,6),]),max(as.in.ge.patterns.human.stat.up[c(3,6),])),ylab='proportion of genes with devAS',main='DevAS in human GE clusters',yaxt='n')
segments(b,as.in.ge.patterns.human.stat.up[c(2,5),],b,as.in.ge.patterns.human.stat.up[c(3,6),])
b = barplot(-as.in.ge.patterns.human.stat.dw[c(1,4),],beside = T,col=c,den=c(-1,30),add=T,xaxt='n',yaxt='n')
segments(b,-as.in.ge.patterns.human.stat.dw[c(2,5),],b,-as.in.ge.patterns.human.stat.dw[c(3,6),])
abline(h=0)
at=-1:3*0.05
axis(2,at,abs(at))
legend('topright',fill='black',den=c(30,-1),legend=c('Early genes','Late Genes'))
text(b[2],par('usr')[4],'AS up',adj=c(-0.1,1.1))
text(b[2],par('usr')[3],'AS down',adj=c(-0.1,-1.1))
plotPanelLetter('F.1')

#plotOhnologsFreq(h.age.ens,ge.info.h$Human_ID[!is.na(ge.info.h$Age) & ge.info.h$Age>3],main='DevAS are depleted in recent duplications',
#								 filter = intersect(ge.info.h$Human_ID[!is.na(ge.info.h$Age)],rownames(ens.ge.cod$human$gene)),ylab='Proportion of recent genes',merge.up.down=TRUE)
plotDevASFreq(h.age.ens,ge.info.h$Human_ID[!is.na(ge.info.h$Age) & ge.info.h$Age>3],'mam. dupl.','topleft',main='DevAS are depleted in recent duplications',
							filter = intersect(ge.info.h$Human_ID[!is.na(ge.info.h$Age)],rownames(ens.ge.cod$human$gene)))

plotPanelLetter('G.1')
#plotOhnologsFreq(h.age.ens,c(ohnologs.h$Ohnolog.1.Id,ohnologs.h$Ohnolog.2.Id),main='DevAS are enriched in ohnologs',
##								 filter = intersect(ge.info.h$Human_ID[!is.na(ge.info.h$Age)], #remove filter from here
#								 rownames(ens.ge.cod$human$gene)),merge.up.down=TRUE)
plotDevASFreq(h.age.ens,c(ohnologs.h$Ohnolog.1.Id,ohnologs.h$Ohnolog.2.Id),'ohnologs','topright',main='DevAS are enriched in ohnologs')
plotPanelLetter('H.1')

t=log2(rbind(inageAS.h[1,],inpanAS.h[1,]))
ci1 = log2(rbind(inageAS.h[3,],inpanAS.h[3,]))
ci2 = log2(rbind(inageAS.h[4,],inpanAS.h[4,]))
#t[2,4] = ci1[2,4] = ci1[2,4] = NA
b=barplot(t,beside = T,ylab='log2(odds ratio)',ylim=range(t,ci1,ci2,na.rm=T),names.arg = c('TFs','SFs','RNA-binding','Pan/DevAS'),legend.text = c('',''),args.legend = list(x='bottomright',title='Assotiation with:',legend=c('AgeAS','PanAS')))
segments(b,ci1,b,ci2)
plotPanelLetter('I.1')

dev.off()


### figure 2.1 plot ####

pdf('figures/paper.figures/3/2.pdf',w=14,h=10)
l=rbind(cbind(matrix(rep(1:2,each=9),nrow=3),matrix(3:20,byrow = T,nrow=3)),
				matrix(c(rep(21:22,each=6),rep(23:24,each=9),rep(25,each=6)),nrow=3),
				matrix(c(rep(26:27,each=6),rep(28:29,each=9),rep(30,each=6)),nrow=3))
layout(l)
par(tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,2.5),oma=c(0,0,0,1))
m = meta
SAJR::plotMDS(points=orth.mds2$sp7,col=m$col,cex=m$cex*1.5,pch=m$pch,main=paste0('All cassette exons (',ancient.cnt,')'))
legend('topleft',col=params$tissue.col,pch=19,legend=names(params$tissue.col),bty='n')
plotPanelLetter('A')
m = meta[rownames(mouse.ad.cor.mds),]
SAJR::plotMDS(points=mouse.ad.cor.mds,col=m$col,cex=m$cex*2,pch=19,main=paste0('All mouse cassette exons (',mouse.cassett.cnt,')'))

par(mar=c(0.2,0,1.5,0.2))
plotPanelLetter('B')
plotMouseTrees(mouse.stage.trees$as,'Splicing','C')
plotMouseTrees(mouse.stage.trees$ge,'Gene expression','D')
plotMouseTrees(mouse.stage.trees$sf,'Gene expression of SF','E')
mar=par(mar=c(4,2,1.5,1.5))
col = unique(meta[,c('tissue','col')])
col = setNames(col$col,substr(col$tissue,1,1))
col = c(col,const='black',bh='orange','2ts'='#666666','3ts'='#777777','4ts'='#AAAAAA','5ts'='#DDDDDD')

areaplot(tr.compl.m,col=col[rownames(tr.compl.m)],ylab='# of exons with PSI in [0.1,0.9]',main='Mouse transc. complexity',xaxt='n',xlab='Age (days from conception)')
axis(1,1:ncol(tr.compl.m),colnames(tr.compl.m))
plotPanelLetter('F')

par(mar=c(5.5,3,1.5,0))
c = rep(params$tissue.col[colnames(as.in.ge.patterns.mouse.stat.up)],each=2)
b = barplot(as.in.ge.patterns.mouse.stat.up[c(1,4),],beside = T,col=c,den=c(-1,30),las=3,ylim=c(-max(as.in.ge.patterns.mouse.stat.dw[c(3,6),]),max(as.in.ge.patterns.mouse.stat.up[c(3,6),])),ylab='proportion of genes with age AS',main='AS in mouse GE clusters',yaxt='n')
segments(b,as.in.ge.patterns.mouse.stat.up[c(2,5),],b,as.in.ge.patterns.mouse.stat.up[c(3,6),])
b = barplot(-as.in.ge.patterns.mouse.stat.dw[c(1,4),],beside = T,col=c,den=c(-1,30),add=T,xaxt='n',yaxt='n')
segments(b,-as.in.ge.patterns.mouse.stat.dw[c(2,5),],b,-as.in.ge.patterns.mouse.stat.dw[c(3,6),])
abline(h=0)
at=-1:3*0.05
axis(2,at,abs(at))
legend('topright',fill='black',den=c(30,-1),legend=c('Early genes','Late Genes'))
text(b[2],par('usr')[4],'AS up',adj=c(-0.1,1.1))
text(b[2],par('usr')[3],'AS down',adj=c(-0.1,-1.1))
plotPanelLetter('G')

plotOhnologsFreq(m.age.ens,ge.info.m$Mouse_ID[!is.na(ge.info.m$Age) & ge.info.m$Age>0],main='AgeAS are depleted in recent genes',
								 filter = intersect(ge.info.m$Mouse_ID[!is.na(ge.info.m$Age)],rownames(ens.ge.cod$mouse$gene)),ylab='Proportion of recent genes',legend=TRUE)
plotPanelLetter('H')
plotOhnologsFreq(m.age.ens,c(ohnologs.m$Ohnolog.1.Id,ohnologs.m$Ohnolog.2.Id),main='AgeAS is enriched in ohnologs',
								 filter = intersect(ge.info.m$Mouse_ID[!is.na(ge.info.m$Age)],rownames(ens.ge.cod$mouse$gene)))
plotPanelLetter('I')

t=log2(rbind(inageAS.m[1,],inpanAS.m[1,]))
ci1 = log2(rbind(inageAS.m[3,],inpanAS.m[3,]))
ci2 = log2(rbind(inageAS.m[4,],inpanAS.m[4,]))
#t[2,4] = ci1[2,4] = ci1[2,4] = NA
b=barplot(t,beside = T,ylab='log2(odds ratio)',ylim=range(t,ci1,ci2,na.rm=T),names.arg = c('TFs','SFs','RNA-binding','Pan/AgeAS'),legend.text = c('',''),args.legend = list(x='bottomright',title='Assotiation with:',legend=c('AgeAS','PanAS')))
segments(b,ci1,b,ci2)
plotPanelLetter('J')

plot.new()
text(0.5,0.5,'Human:',cex=3)

c = rep(params$tissue.col[colnames(as.in.ge.patterns.human.stat.up)],each=2)
b = barplot(as.in.ge.patterns.human.stat.up[c(1,4),],beside = T,col=c,den=c(-1,30),las=3,ylim=c(-max(as.in.ge.patterns.human.stat.dw[c(3,6),]),max(as.in.ge.patterns.human.stat.up[c(3,6),])),ylab='proportion of genes with age AS',main='AS in human GE clusters',yaxt='n')
segments(b,as.in.ge.patterns.human.stat.up[c(2,5),],b,as.in.ge.patterns.human.stat.up[c(3,6),])
b = barplot(-as.in.ge.patterns.human.stat.dw[c(1,4),],beside = T,col=c,den=c(-1,30),add=T,xaxt='n',yaxt='n')
segments(b,-as.in.ge.patterns.human.stat.dw[c(2,5),],b,-as.in.ge.patterns.human.stat.dw[c(3,6),])
abline(h=0)
at=-1:3*0.05
axis(2,at,abs(at))
legend('topright',fill='black',den=c(30,-1),legend=c('Early genes','Late Genes'))
text(b[2],par('usr')[4],'AS up',adj=c(-0.1,1.1))
text(b[2],par('usr')[3],'AS down',adj=c(-0.1,-1.1))
plotPanelLetter('G.1')

plotOhnologsFreq(h.age.ens,ge.info.h$Human_ID[!is.na(ge.info.h$Age) & ge.info.h$Age>0],main='AgeAS are depleted in recent genes',
								 filter = intersect(ge.info.h$Human_ID[!is.na(ge.info.h$Age)],rownames(ens.ge.cod$human$gene)),ylab='Proportion of recent genes')
plotPanelLetter('H.1')
plotOhnologsFreq(h.age.ens,c(ohnologs.h$Ohnolog.1.Id,ohnologs.h$Ohnolog.2.Id),main='AgeAS is enriched in ohnologs',
								 filter = intersect(ge.info.h$Human_ID[!is.na(ge.info.h$Age)],rownames(ens.ge.cod$human$gene)))
plotPanelLetter('I.1')

t=log2(rbind(inageAS.h[1,],inpanAS.h[1,]))
ci1 = log2(rbind(inageAS.h[3,],inpanAS.h[3,]))
ci2 = log2(rbind(inageAS.h[4,],inpanAS.h[4,]))
#t[2,4] = ci1[2,4] = ci1[2,4] = NA
b=barplot(t,beside = T,ylab='log2(odds ratio)',ylim=range(t,ci1,ci2,na.rm=T),names.arg = c('TFs','SFs','RNA-binding','Pan/AgeAS'),legend.text = c('',''),args.legend = list(x='bottomright',title='Assotiation with:',legend=c('AgeAS','PanAS')))
segments(b,ci1,b,ci2)
plotPanelLetter('J.1')

dev.off()



##################
### figure 3 #####

f = sapply(names(ens.ge.cod.tsm),function(s)orth.ens.genes[,s] %in% rownames(ens.ge.cod.tsm[[s]]))
f = apply(f,1,sum)==7
ens.ge.cod.tsm.log = lapply(names(ens.ge.cod.tsm),function(s)log2(ens.ge.cod.tsm[[s]][orth.ens.genes[f,s],]+0.1))
names(ens.ge.cod.tsm.log) = names(ens.ge.cod.tsm)

# tissue.stage.cor.as = lapply(unique(meta$tissue),function(t){
# 	r=lapply(1:nrow(age.al.i),function(i)getSpeciesCor(orth.seg.ad.tsm,age.al.i[i,c(1:7)],t))
# 	names(r) = age.al.i[,'mouse']
# 	r
# })
# 
# tissue.stage.cor.ge = lapply(unique(meta$tissue),function(t){
# 	r=lapply(1:nrow(age.al.i),function(i)getSpeciesCor(ens.ge.cod.tsm.log,age.al.i[i,c(1:7)],t))
# 	names(r) = age.al.i[,'mouse']
# 	r
# })
# 
# names(tissue.stage.cor.as) = names(tissue.stage.cor.ge) = unique(meta$tissue)
# 
# sp= rownames(species)[c(3:6)]
# njLength = function(c){
# 	if(sum(is.na(c)) > 0) return(NA)
# 	sum(nj(1-c)$edge.length)
# }
# nj.tree.len.mrbo.as=calcBootstrapSpeciesDiv(orth.seg.ad.tsm,sp,njLength,age.al.i,100)
# nj.tree.len.mrbo.ge=calcBootstrapSpeciesDiv(ens.ge.cod.tsm.log,sp,njLength,age.al.i,100)
# sp= rownames(species)[c(1,3:6)]
# nj.tree.len.hmrbo.as=calcBootstrapSpeciesDiv(orth.seg.ad.tsm,sp,njLength,age.al.i,100)
# nj.tree.len.hmrbo.ge=calcBootstrapSpeciesDiv(ens.ge.cod.tsm.log,sp,njLength,age.al.i,100)
# save(nj.tree.len.mrbo.as,nj.tree.len.mrbo.ge,nj.tree.len.hmrbo.as,nj.tree.len.hmrbo.ge,file='Rdata/paper.figures/nj.trees.bootstrap.Rdata')
# number of events
#alternification
sp.groups = c(species$short,'hq','mr','mrb','hqmrb')
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

mean = table(alt.sp[f])[sp.groups]

#born/lost
born.exn.sajr = readRDS('Rdata/born.exn.sajr.Rdata')
exon.birth.one = readRDS('Rdata/exon.birth.one.Rdata')
sp.birth = sapply(exon.birth.one,function(x){paste(species$short[!is.na(x$seg_id)],collapse='')})

tsp.birth = table(sp.birth)
brth=tsp.birth[sp.groups]
lost=tsp.birth[sapply(sp.groups,function(x){paste(setdiff(species$short,strsplit(x,'')[[1]]),collapse='')})]
names(lost) = sp.groups

pdf('figures/paper.figures/3/3.pdf',w=10,h=8)
par(mfrow=c(2,3),tck=-0.01,mgp=c(2,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,0,1))

x=setNames(rep('',11),sp.groups)
plotDivergenceOnAgeWithConf(nj.tree.len.mrbo.as,params$tissue.col,ylab='Total tree length (1-cor)',main='AS divergence (mrbo)')
plotPanelLetter('A')
plotDivergenceOnAgeWithConf(nj.tree.len.mrbo.ge,params$tissue.col,ylab='Total tree length (1-cor)',main='GE divergence (mrbo)')
plotPanelLetter('B')
plotSpeciesTree.3('C',x,mean,brth,lost)

plotDivergenceOnAgeWithConf(nj.tree.len.hmrbo.as,params$tissue.col,ylab='Total tree length (1-cor)',main='AS divergence (hmrbo)')
plotPanelLetter('A.1')
plotDivergenceOnAgeWithConf(nj.tree.len.hmrbo.ge,params$tissue.col,ylab='Total tree length (1-cor)',main='GE divergence (hmrbo)')
plotPanelLetter('B.1')
dev.off()





##########################
### figure 4 prepare #####
#A
DPSI=0.5
s = 'mouse'
age.segs = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,DPSI,border.stages,s)[anns[[s]]$sites=='ad',])
names(age.segs) = rownames(species)
age.segs.cod = lapply(rownames(species),function(s)age.segs[[s]][anns[[s]][rownames(age.segs[[s]]),'cod']!='n',])
names(age.segs.cod) = rownames(species)

ts = unique(meta$tissue)
up.all = t(sapply(age.segs,function(x)apply(x=='u',2,sum)[ts]))
up.cod = t(sapply(age.segs.cod,function(x)apply(x=='u',2,sum)[ts]))

dw.all = t(sapply(age.segs,function(x)apply(x=='d',2,sum)[ts]))
dw.cod = t(sapply(age.segs.cod,function(x)apply(x=='d',2,sum)[ts]))
#B
ts = unique(meta$tissue)[-2]
up.ntiss = log(t(sapply(age.segs[-2],function(x)table(factor(apply(x[,ts]=='u',1,sum),levels=1:length(ts)))))+1)
dw.ntiss = -log(t(sapply(age.segs[-2],function(x)table(factor(apply(x[,ts]=='d',1,sum),levels=1:length(ts)))))+1)

#C
age.ad.over = lapply(age.segs,function(x){
	x[x=='-'] = NA
	x = cbind(x=='u',x=='d')
	#x = x[apply(x,1,sum,na.rm=T)>0,]
	colnames(x) = paste(colnames(x),rep(c('up','dw'),each=ncol(x)/2))
	caclSegOverlap(x)
})

#D
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
cnst.len3.freq=c(0.4017701,0.3943019,0.4092728)#my.binom.test(table(all.anns[[s]][all.anns[[s]]$sites=='ad' & all.anns[[s]]$type=='EXN' & all.anns[[s]]$cod=='c' & all.anns[[s]]$gene_id %in% all.anns[[s]][rownames(age.segs[[s]])[apply(age.segs[[s]] =='u' | age.segs[[s]] =='d',1,sum)>0],'gene_id'],'length'] %% 3 == 0)[c('TRUE','FALSE')])
#E phastcons
h = all.anns$human
#rm(all.anns);gc()
phastcons = read.table('/home/mazin/skoltech/projects/evo.devo/processed/ad.phastcons.gz',sep='\t')
phastcons = setNames(strsplit(phastcons[,2],',',TRUE),phastcons[,1])
phastcons = lapply(phastcons,as.numeric)
age.ad.ph = list()
as = age.segs$human[anns$human[rownames(age.segs$human),'cod']=='c',]
inx = 151:200
h.age.sids = list()
h.age.sids$cnst = rownames(h)[h$sites=='ad' & h$cod=='c' & h$type=='EXN' & h$gene_id %in% h[rownames(as),'gene_id']]
for(t in colnames(as)){
	h.age.sids[[paste(t,'up')]] = rownames(as)[as[,t]=='u']
	h.age.sids[[paste(t,'dw')]] = rownames(as)[as[,t]=='d']
}

age.ad.ph = lapply(h.age.sids,function(sids)sapply(phastcons[intersect(names(phastcons),sids)],function(x)mean(x[c(inx,length(x)-inx)])))

#F gnomad
gnomad = read.table('processed/gnomad201/human.ad.snp.tab.gz',sep='\t')
colnames(gnomad) = c('chr_id','pos','id','ref','alt','qual','filter','seg_id','alt_cnt','freq','tot_cnt')
hann = all.anns$human[unique(unlist(h.age.sids)),]
gnomad = gnomad[gnomad$alt_cnt > 1 & gnomad$seg_id %in% rownames(hann),]
s2g = hann[gnomad$seg_id,]
gnomad$dist2start = gnomad$pos - s2g$start
gnomad$dist2stop  = s2g$stop - gnomad$pos
gnomad[s2g$strand==-1,c('dist2start','dist2stop')] = gnomad[s2g$strand==-1,c('dist2stop','dist2start')]
gnomad$strand = s2g$strand

age.ad.gnomad = lapply(h.age.sids,function(sids)log10(gnomad$freq[gnomad$seg_id %in% sids & ((gnomad$dist2start<0 & gnomad$dist2start>= -50) | (gnomad$dist2stop<0 & gnomad$dist2stop>= -50))]))

#G hgmd
hgmd = read.table('input/hgmd/tosend/2017_1_HGMD_ALL_Variants.csv',sep='\t',row.names=1,header=T,quote='',comment.char = '')
colnames(hgmd)
table(nchar(hgmd$ref_VCF_hg19),hgmd$tag)
hgmd.gr = GRanges(hgmd$chrom_VCF_hg19,IRanges(hgmd$pos_VCF_hg19,hgmd$pos_VCF_hg19+nchar(hgmd$ref_VCF_hg19)-1))

seg.gr = GRanges(hann$chr_id,IRanges(hann$start,hann$stop))
hgmd2hadf = findOverlaps(hgmd.gr,seg.gr,maxgap=200,type='any',select='all',ignore.strand=TRUE)
hgmd2hadf = data.frame(h=hgmd2hadf@from,s=hgmd2hadf@to)
hgmd2hadf$seg.id = rownames(hann)[hgmd2hadf$s]
hgmd2hadf$hgmd.id = rownames(hgmd)[hgmd2hadf$h]

s = hann[hgmd2hadf$s,]
h = hgmd[hgmd2hadf$h,]

hgmd2hadf$position = ''
hgmd2hadf$position = paste0(hgmd2hadf$position,ifelse((s$strand== 1 & s$start > h$pos_VCF_hg19) | (s$strand==-1 & s$stop < h$pos_VCF_hg19+nchar(h$ref_VCF_hg19)-1),'u',''))
hgmd2hadf$position = paste0(hgmd2hadf$position,ifelse(s$start<=h$pos_VCF_hg19+nchar(h$ref_VCF_hg19)-1 & s$stop>=h$pos_VCF_hg19,'i',''))
hgmd2hadf$position = paste0(hgmd2hadf$position,ifelse((s$strand==-1 & s$start > h$pos_VCF_hg19) | (s$strand== 1 & s$stop < h$pos_VCF_hg19+nchar(h$ref_VCF_hg19)-1),'d',''))
table(hgmd2hadf$position)
# locaction: intron  (u,w), exon (e), ppt (p), acc.ag (A),acc(a), don(d),don.gt(D)
# mutation annotation priproti: dinucleotides > other nt of canonical sites > PPT > intron
hgmd2hadf$loc = ''
hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse(hgmdOverlapLocaction(h,s,0,Inf,T) & hgmdOverlapLocaction(h,s,0,Inf,F),'e',''))
hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse(hgmdOverlapLocaction(h,s,-2,-1,T),'A',''))
hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse(hgmdOverlapLocaction(h,s,-2,-1,F),'D',''))
hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse((!grepl('A|D',hgmd2hadf$loc)) & (hgmdOverlapLocaction(h,s,-3,-3,T) | hgmdOverlapLocaction(h,s,1,1,T)),'a',''))
hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse((!grepl('A|D',hgmd2hadf$loc)) & (hgmdOverlapLocaction(h,s,-5,-3,F) | hgmdOverlapLocaction(h,s,1,4,F)),'d',''))
hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse((!grepl('A|D|a|d',hgmd2hadf$loc)) & hgmdOverlapLocaction(h,s,-23,-4,T),'p',''))
hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse((!grepl('A|D|a|d|p',hgmd2hadf$loc)) & hgmdOverlapLocaction(h,s,-200,-24,T),'u',''))
hgmd2hadf$loc = paste0(hgmd2hadf$loc,ifelse((!grepl('A|D|a|d|p',hgmd2hadf$loc)) & hgmdOverlapLocaction(h,s,-200,-6,F),'w',''))
sort(table(hgmd2hadf$loc))
hs = unique(hgmd2hadf$seg.id[grep('A|a|p|D|d|u|d',hgmd2hadf$loc)])
age.ad.hgmd = t(sapply(h.age.sids,function(sids){r=c(sum(sids %in% hs),sum(!(sids %in% hs)));c(r,my.binom.test(r))}))

sps = c('mouse','rat','rabbit','human','opossum','chicken')
orth.age.ad5 = lapply(rownames(species),function(s)getAgeASchanges(orth.seg.ad.tsm,meta.tsm,psi.thr = 0.5,border.stages,s))
names(orth.age.ad5) = rownames(species)
u5=getASChangeCons(orth.age.ad5,'u',sps)
d5=getASChangeCons(orth.age.ad5,'d',sps)

orth.age.ad3 = lapply(rownames(species),function(s)getAgeASchanges(orth.seg.ad.tsm,meta.tsm,psi.thr = 0.3,border.stages,s))
names(orth.age.ad3) = rownames(species)
u3=getASChangeCons(orth.age.ad3,'u',sps)
d3=getASChangeCons(orth.age.ad3,'d',sps)





#######################
### figure 4 plot #####
pdf('figures/paper.figures/3/4.pdf',w=10,h=11)
s = 'mouse'
l = matrix(c(1,2,4,1,2,5,3,3,6,3,3,7),ncol=4)
l = rbind(l,8:11)
layout(l)
par(tck=-0.02,mgp=c(1.1,0.2,0),mar=c(3,2,1.5,0),oma=c(0,0,2,1))
cols=rep(params$tissue.col[colnames(up.cod)],each=nrow(up.all))
b=barplot( up.all,beside = T,col=NA,border = cols,ylim=c(-max(dw.all,na.rm=T),max(up.all,na.rm=T)),yaxt='n',ylab='# of exons',main='AgeAS cassette exons')
barplot( up.cod,beside = T,col=cols,border = cols,add=T,density = 50,yaxt='n')
barplot(-dw.all,beside = T,col=NA,border = cols,add=T,yaxt='n')
barplot(-dw.cod,beside = T,col=cols,border = cols,add=T,density = 50,yaxt='n')
at = c(-5:5*100)
axis(2,at,abs(at))
text(b[1:nrow(up.all)],rep(-500,nrow(up.all)),rownames(up.all),srt=90,adj = c(0,0.5),cex=0.6)
text(b[21], 500,'AS up',adj=c(-0.1,1.1))
text(b[21],-500,'AS down',adj=c(-0.1,-1.1))
plotPanelLetter('A')


cols = rep(gray.colors(nrow(up.ntiss)),times=ncol(up.ntiss))
b=barplot(up.ntiss,beside = T,xlab='# of tissues',col=cols,border=NA,yaxt='n',ylim=range(up.ntiss,dw.ntiss),ylab='# of exons',main='Number of shared age-related exons',legend.text=T,args.legend = list(x='bottomright'))
barplot(dw.ntiss,beside = T,col=cols,border=NA,yaxt='n',ylim=range(up.ntiss,dw.ntiss),add=TRUE)
abline(h=0)
at = c(1000,100,10,0,10,100,1000)
axis(2,c(-1,-1,-1,1,1,1,1)*log(at+1),at)
text(b[21], log(500),'AS up',adj=c(-0.1,1.1))
text(b[21],-log(500),'AS down',adj=c(-0.1,-1.1))
plotPanelLetter('B')

par(mar=c(1.5,6,3,1))
plotAgeSegOverlap(age.ad.over[[s]],main=paste('Overlap of ageAS exons across',s,'tissues'))
par(mar=c(3,2,1.5,0))
plotPanelLetter('C')

par(mar=c(5,2,1.5,0))
at = (0:21)[-seq(2,21,by = 3)]
cols = c('#000000',rep(params$tissue.col,each=2))
pch=c(1,rep(c(19,1),times=7))
xax = setNames(c(0,seq(2.5,21,by=3)),c('const.',colnames(up.ad.len)))
stat = matrix(NA,ncol=3,nrow=15)
stat[1,] = cnst.len3.freq[c(2,1,3)]
stat[seq(2,15,by = 2),] = t(up.ad.len)[,c(4,3,5)]
stat[seq(3,15,by = 2),] = t(dw.ad.len)[,c(4,3,5)]

plotASSegStat(at,stat,up2dw.len.pv,cols,xax,lty=c(2,rep(1:2,times=7)),main='Proportion of 3N exons',ylab='proportion',pch=pch)
plotPanelLetter('D')

stat = t(sapply(age.ad.ph,function(x){m=mean(x);s=sd(x)/sqrt(length(x));c(m-2*s,m,m+2*s)}))
pv = sapply(1:7,function(i)wilcox.test(age.ad.ph[[2*i]],age.ad.ph[[2*i+1]])$p.value)
plotASSegStat(at,stat,pv,cols,xax,lty=c(2,rep(1:2,times=7)),main='Splice site conservation',ylab="3'-ss PhastCons",pch=pch)
plotPanelLetter('E')

stat = t(sapply(age.ad.gnomad,function(x){m=mean(x);s=sd(x)/sqrt(length(x));c(m-2*s,m,m+2*s)}))
pv = sapply(1:7,function(i)wilcox.test(age.ad.gnomad[[2*i]],age.ad.gnomad[[2*i+1]])$p.value)
plotASSegStat(at,stat,pv,cols,xax,lty=c(2,rep(1:2,times=7)),main='Splice site human constrains',ylab="log10(SNP freq)",pch=pch)
pv = sapply(2:15,function(i)wilcox.test(age.ad.gnomad$cnst,age.ad.gnomad[[i]])$p.value)
legend('topright',title = 'p-value',legend=c('* <0.05','** <0.01','*** <0.001'))
plotPanelLetter('F')

pv = sapply(1:7,function(t){s=c(age.ad.hgmd[t*2,1],age.ad.hgmd[t*2+1,1]);prop.test(s,s+c(age.ad.hgmd[t*2,2],age.ad.hgmd[t*2+1,2]))$p.value})
plotASSegStat(at,age.ad.hgmd[,c(4,3,5)],pv,cols,xax,lty=c(2,rep(1:2,times=7)),main='Human disease mutations',ylab='exons with splice site mutations',pch=pch)
plotPanelLetter('G')

plotTisUpDownCOns(u3$brain,d3$brain,'H',main='Brain ageAS conservation')
plotTisUpDownCOns(u3$heart,d3$heart,'I',main='Heart ageAS conservation')
plotTisUpDownCOns(u3$liver,d3$liver,'J',main='Liver ageAS conservation')
plotTisUpDownCOns(u3$testis,d3$testis,'K',main='Testis ageAS conservation')
legend('topright',fill=c('red','blue'),legend=c('Inclusion','Exclusion'))
dev.off()

#########################
### 5 Direction new/alternification
# new exons
sps=list(h='h',`m/r` = c('m','r'),hq='hq',mr='mr',mrb='mrb',hqmrb='hqmrb',hqmrbo='hqmrbo')
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
	calcMeanCols(psi,paste(m$tissue,m$stage))
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
	tsm = born.exn.sajr[[s]]$ir[,colnames(born.exn.sajr[[s]]$ir) %in% rownames(meta)]
	t = meta[colnames(tsm),]
	tsm = born.exn.tsm[[s]]
	born.exn.sajr[[s]]$psi.tsm.al = matrix(NA,ncol=nrow(m),nrow=nrow(tsm),dimnames = list(rownames(tsm),rownames(m)))
	for(i in 1:nrow(m)){
		stage = age.al[age.al$mouse==m$mouse.stage[i],s]
		for(t in ts){
			st = paste(t,stage)
			if(st %in% colnames(tsm))
				born.exn.sajr[[s]]$psi.tsm.al[,paste(t,m$mouse.stage[i])] = tsm[,st]
		}
	}
}
meta.tsm.al$col = paste0(params$tissue.col[meta.tsm.al$tissue],rep(c('44','55','66','77','88','99','AA','BB','CC','DD','EE','FF'),times=7))



max.stage=unlist(lapply(setNames(born.exn.sajr,NULL),function(x)setNames(colnames(x$psi.tsm.al)[apply(x$psi.tsm.al,1,function(x){r = which.max(x);ifelse(length(r)==0 || sum(x==x[r],na.rm=T)>1,NA,r)})],rownames(x$psi.tsm.al))))
sids = lapply(split.data.frame(born.seg.ids,sp.birth),function(x){x[!is.na(x)]})

born.tis.dpsi = lapply(born.exn.tsm,function(x){
	r=sapply(unique(meta$tissue),function(t){
		y = x[,grep(t,colnames(x))]
		apply(y,1,max,na.rm=T)-apply(y,1,min,na.rm=T)
		})
	colnames(r) = unique(meta$tissue)
	r
})
#par(mfrow=c(3,3))
#sapply(1:7,function(i){hist(born.tis.dpsi[[i]][born.per.tissue.age.qv[[i]]>0.05],0:50/50,border = NA,col='gray',freq = F,main=names(born.tis.dpsi)[i]);hist(born.tis.dpsi[[i]][born.per.tissue.age.qv[[i]]<0.05],0:50/50,border = NA,col='#FF000080',add=T,freq = F,main=names(born.tis.dpsi)[i])})
born.devAS = matrix(NA,ncol=ncol(born.seg.ids),nrow=nrow(born.seg.ids))
colnames(born.devAS) = colnames(born.seg.ids)
for(s in colnames(born.seg.ids)){
	for(i in 1:nrow(born.seg.ids)){
		if(!is.na(born.seg.ids[i,s])){
				born.devAS[i,s] = sum(born.per.tissue.age.qv[[s]][born.seg.ids[i,s],]<0.05 & born.tis.dpsi[[s]][born.seg.ids[i,s],]>0.2,na.rm=TRUE)>0
		}
	}
}


born.ex.prop.sgn.dpsi0.2=getSgnFraqWithBootstrap(sp.birth,sps,born.devAS[,-c(2,7)],N=500)

## cleaned
#minus:
# macaque and chicken
# kidney, ovary and cerebellum
# stages 10.5, 11.5, 3dpb, 9wpb

missed.cond = lapply(born.exn.sajr,function(x)colnames(x$psi.tsm.al)[apply(!is.na(x$psi.tsm.al),2,sum)==0])
lapply(missed.cond[c(-2,-7)],function(x)x[!grepl('10.5|11.5|3dpb|9wpb',x) & !grepl('cerebellum|kidney|ovary',x)])
mm = meta.tsm.al[!(meta.tsm.al$tissue %in% c('cerebellum','kidney','ovary') | meta.tsm.al$mouse.stage %in% c('10.5','11.5','3dpb','9wpb')),]
tmp = lapply(setNames(born.exn.sajr,NULL)[c(-2,-7)],function(x)x$psi.tsm.al[,rownames(mm)])
max.stage.cl=unlist(lapply(tmp,function(x)setNames(colnames(x)[apply(x,1,function(x){r = which.max(x);ifelse(length(r)==0 || sum(x==x[r],na.rm=T)>1,NA,r)})],rownames(x))))

born.sids = lapply(split.data.frame(born.seg.ids,sp.birth),function(x){x[!is.na(x)]})

# alternification
orth.seg.ad.all.tsm = readRDS('Rdata/orth.seg.ad.all.tsm.Rdata')
# alt.thrs = c(0.8,0.8)
# min.psi = sapply(orth.seg.ad.all.tsm,function(x)apply(x,1,min,na.rm=T))
# max.na.fraq = sapply(orth.seg.ad.all.tsm,function(x)apply(is.na(x),1,mean))
# hist(apply(max.na.fraq,1,max))
# f = apply(max.na.fraq,1,max)<0.5
# table(f)
# hist(min.psi[f,1],200)
# alt.sp = apply(min.psi,1,function(x){
# 	r = '-'
# 	if(sum(is.na(x)) == 0 && sum(x>alt.thrs[1] & x<alt.thrs[2])==0){
# 		r=paste(species$short[x<alt.thrs[1]],collapse='')
# 	}
# 	r
# })

# alt.sp = sapply(orth.seg.ad.all,function(x)x$seg$type=='ALT')
# alt.sp = apply(alt.sp,1,function(x)paste(species$short[x],collapse = ''))
# names(alt.sp) = rownames(orth.seg.ad.all$human$seg)
# saveRDS(alt.sp,'Rdata/paper.figures/alt.sp.Rdata')
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
min.stage = sapply(orth.tsm.al,function(x)colnames(x)[apply(x,1,function(x){r=which.min(x);ifelse(length(r)==0 || sum(x==x[r],na.rm=T)>1,NA,r)})])
names(min.stage) = rownames(orth.tsm.al$human)

for(i in 1:nrow(min.stage))	min.stage[i,rownames(species)[!(species$short %in% strsplit(alt.sp[i],'')[[1]])]] = NA

mm = meta.tsm.al[!(meta.tsm.al$tissue %in% c('cerebellum','kidney','ovary') | meta.tsm.al$mouse.stage %in% c('10.5','11.5','3dpb','9wpb')),]
tmp = lapply(orth.tsm.al[c(-2,-7)],function(x)x[,rownames(mm)])

min.stage.cl = sapply(tmp,function(x)colnames(x)[apply(x,1,function(x){r=which.min(x);ifelse(length(r)==0|| sum(x==x[r],na.rm=T)>1,NA,r)})])
rownames(min.stage.cl) = rownames(tmp$human)
for(i in 1:nrow(min.stage.cl))	min.stage.cl[i,colnames(min.stage.cl) %in% rownames(species)[!(species$short %in% strsplit(alt.sp[i],'')[[1]])]] = NA
#min.stage.cl=unlist(lapply(colnames(min.stage.cl),function(s)setNames(min.stage.cl[,s],rownames(tmp[[s]]))))


orth.seg.ad.all.dpsi = lapply(orth.seg.ad.all.tsm,function(x){
	sapply(unique(meta$tissue),function(t){
		p = x[alt.sp!='',grep(t,colnames(x))]
		apply(p,1,max,na.rm=T) - apply(p,1,min,na.rm=T)
		})
	})


# par(mfrow=c(3,3))
# sapply(1:7,function(s){hist(orth.seg.ad.all.dpsi[[s]][orth.per.tissue.age.qv[[s]]>0.05],0:50/50,border = NA,col='gray',freq = F,main=names(born.tis.dpsi)[i]);hist(orth.seg.ad.all.dpsi[[s]][orth.per.tissue.age.qv[[s]]<0.05],0:50/50,border = NA,col='#FF000080',add=T,freq = F,main=names(born.tis.dpsi)[s])})

alt.devAS = sapply(names(orth.per.tissue.age.qv),function(s){apply(orth.per.tissue.age.qv[[s]]<0.05 & orth.seg.ad.all.dpsi[[s]]>0.2,1,sum,na.rm=T)>0})
#sps=list(h='h',`m/r` = c('m','r'),hq='hq',mr='mr',mrb='mrb',hqmrb='hqmrb',hqmrbo='hqmrbo',hqmrboc='hqmrboc')
altern.prop.sgn.dpsi0.2=getSgnFraqWithBootstrap(alt.sp[alt.sp!=''],sps.,alt.devAS[,-c(2,7)],N=500)


## %% 3

#sps=list(`m/r` = c('m','r'),mr='mr',mrb='mrb',hqmrb='hqmrb',hqmrbo='hqmrbo')
seg.len = unlist(lapply(setNames(born.exn.sajr,NULL),function(x)setNames(x$seg$length,rownames(x$seg))))
born.exn.len3 = sapply(sps,function(s){r=seg.len[unique(unlist(born.sids[s]))] %% 3 == 0;my.binom.test(sum(r),sum(!r))})

#orth.seg.ad.all = readRDS('Rdata/orth.seg.ad.all.Rdata')
seg.len = readRDS('Rdata/paper.figures/alt.seg.len.Rdata')#unlist(lapply(setNames(orth.seg.ad.all,NULL),function(x)setNames(x$seg$length,rownames(x$seg))))
#saveRDS(seg.len,'Rdata/paper.figures/alt.seg.len.Rdata')
alt.exn.len3 = sapply(sps.,function(s){r=seg.len[names(alt.sp)[alt.sp %in% s]] %% 3 == 0;my.binom.test(sum(r),sum(!r))})


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



#rm(orth.seg.ad.all);gc()

# exont
# EXONT:000074 - Intrinsically_unstructured_polypeptide_region
o = get_OBO('processed/exon.onthology/exont.obo')
seg2exont = readRDS('Rdata/seg2exont.Rdata')
exont2seg = revList(seg2exont)

iup = exont2seg[['EXONT:000074']]
aex = setdiff(names(seg2exont),iup)
#sps=list(h= 'h',hq='hq',hqmrb='hqmrb',hqmrbo='hqmrbo')

iup.new=sapply(sps,function(x){
	sids = unique(unlist(born.sids[x]))
	sids = sids[grep('hum.',sids)]
	#f = intersect(sids,names(seg2exont)) %in% iup
	f = sids %in% iup
	r = c(sum(f),sum(!f))
	c(r,my.binom.test(r))
	})

iup.alt=sapply(sps,function(x){
	sids = names(alt.sp)[alt.sp %in% x]
	#f = intersect(sids,names(seg2exont)) %in% iup
	f = sids %in% iup
	r = c(sum(f),sum(!f))
	c(r,my.binom.test(r))
})

plotArea(1:4,t(iup.new[3:5,]),col='red',new=T)#,ylim=range(iup.new[3:5,],iup.alt[3:5,]))
plotArea(1:4,t(iup.alt[3:5,]),col='blue')

#############
# plot 5 ####
spes2use = rownames(species)[c(-2,-7)]#c('mouse','rat')
#sps=list(`m/r` = c('m','r'),mr='mr',mrb='mrb',hqmrb='hqmrb',hqmrbo='hqmrbo')
#max.stage.cl. = max.stage.cl[grep(paste(substr(spes2use,1,3),collapse='|'),names(max.stage.cl))]

tab.exn = sapply(sps, function(x){table(factor(max.stage.cl[intersect(unlist(born.sids[x]),names(max.stage.cl))],levels=rownames(meta.tsm.al)))})
tab.alt = sapply(sps.,function(s){table(factor(min.stage.cl[alt.sp %in% s,spes2use],levels=rownames(meta.tsm.al)))})

born.devAS. = born.seg.ids[!is.na(born.devAS) & born.devAS]

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

pdf('figures/paper.figures/3.2/5AB.by.devAS.pdf',w=9,h=6)
par(mfrow=c(2,3),tck=-0.02,mgp=c(2,0.2,0),mar=c(3.5,3,1.5,0),oma=c(0,0,2,1))
barplot(sweep(tab.exn,2,apply(tab.exn,2,sum),'/'),col=meta.tsm.al$col,border = NA,las=3,main='New exons, all',ylab='proportion of exons mostly included in')
barplot(sweep(tab.exn.dev,2,apply(tab.exn.dev,2,sum),'/'),col=meta.tsm.al$col,border = NA,las=3,main='New exons, devAS',ylab='proportion of exons mostly included in')
barplot(sweep(tab.exn.ndev,2,apply(tab.exn.ndev,2,sum),'/'),col=meta.tsm.al$col,border = NA,las=3,main='New exons, non-devAS',ylab='proportion of exons mostly included in')

barplot(sweep(tab.alt,2,apply(tab.alt,2,sum),'/'),col=meta.tsm.al$col,border = NA,las=3,main='Alternification, all',ylab='proportion of exons mostly skipped in')
barplot(sweep(tab.alt.dev,2,apply(tab.alt.dev,2,sum),'/'),col=meta.tsm.al$col,border = NA,las=3,main='Alternification, devAS',ylab='proportion of exons mostly skipped in')
barplot(sweep(tab.alt.ndev,2,apply(tab.alt.ndev,2,sum),'/'),col=meta.tsm.al$col,border = NA,las=3,main='Alternification, non-dev',ylab='proportion of exons mostly skipped in')
dev.off()

pdf('figures/paper.figures/3.2/5.pdf',w=9,h=9)
par(mfrow=c(3,3),tck=-0.02,mgp=c(2,0.2,0),mar=c(3.5,3,1.5,0),oma=c(0,0,2,1))

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
segments(length(sps.),alt.exn.len3[2,length(sps.)],8,alt.exn.len3[3,8],col='orange')
axis(1,1:ncol(alt.exn.len3),colnames(alt.exn.len3),las=3)
abline(h=1/3,lty=2)
plotPanelLetter('E')
legend('topleft',fill=c('red','blue'),legend=c('new exons','alternification'))


plotArea(1:length(sps),p = born.ex.prop.sgn.dpsi0.2[,-2],col='red',lwd=3,new = T,xlim=c(1,length(sps.)),ylim=range(0,born.ex.prop.sgn.dpsi0.2[,-2],altern.prop.sgn.dpsi0.2[,-2]),xaxt='n',xlab='',ylab='weighted proportion of devAS',main='devAS')
plotArea(1:length(sps),p = altern.prop.sgn.dpsi0.2[-length(sps.),-2],col='blue',lwd=3,new = F)
points(length(sps.),altern.prop.sgn.dpsi0.2[length(sps.),1],pch=19,col='orange')
segments(length(sps.),altern.prop.sgn.dpsi0.2[length(sps.),3],length(sps.),altern.prop.sgn.dpsi0.2[length(sps.),4],col='orange')
axis(1,1:length(sps.),colnames(alt.exn.len3),las=3)
abline(h=0,lty=2)
plotPanelLetter('F')

plotArea(1:length(sps),t(alt.psi.on.ev[,-length(sps.)]),col='blue',new=T,ylim=c(0,1),lwd=3,xaxt='n',xlab='',ylab='mean(PSI)',main='mean(PSI)',xlim=c(1,length(sps.)))
plotArea(1:length(sps),t(brn.psi.on.ev[,-length(sps.)]),col='red',lwd=3,new=F)
points(length(sps.),alt.psi.on.ev[1,length(sps.)],pch=19,col='orange')
segments(length(sps.),alt.psi.on.ev[2,length(sps.)],length(sps.),alt.psi.on.ev[3,length(sps.)],col='orange')
axis(1,1:length(sps.),colnames(alt.exn.len3),las=3)
plotPanelLetter('G')

dev.off()


########################
#### fig 6 hexamers ####
hex.stat = rbind(apply(hex.ups.age03$up$ih.qv<0.05,2,sum),
								 apply(hex.dws.age03$up$ih.qv<0.05,2,sum),
								 apply(hex.ups.age03$dw$ih.qv<0.05,2,sum),
								 apply(hex.dws.age03$dw$ih.qv<0.05,2,sum))

# proportion of known
hex2mot = read.table('output/hex2mot2sf.tab.gz')
hex.tis.no = apply(hex.ups.age03$up$ih.qv<0.05 | hex.dws.age03$up$ih.qv<0.05 | hex.ups.age03$dw$ih.qv<0.05 | hex.dws.age03$dw$ih.qv<0.05,1,sum)
known.hex.stat = table(hex.tis.no,known=names(hex.tis.no) %in% hex2mot$V1[hex2mot$V2!=''])

# hexamer overlap
all.hex.qv = cbind(hex.ups.age03$up$ih.qv, hex.dws.age03$up$ih.qv, hex.ups.age03$dw$ih.qv, hex.dws.age03$dw$ih.qv)
colnames(all.hex.qv) = paste0(substr(colnames(all.hex.qv),1,1),rep(c('iu','id','eu','ed'),each=7))
all.hex.qv = all.hex.qv[,rep(c(1,8,15,22),times=7)+rep(0:6,each=4)]
f = apply(all.hex.qv<0.05,1,sum)

hex.overlap=calcAllPairsFT(all.hex.qv<0.05)
# QKI 
fa = readRDS('Rdata/ad.alt.fa.Rdata')
age03 = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.3,border.stages,s))
age01 = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s))
names(age01) = names(age03) = rownames(species)

dim(age03$human)
table(age03$human[,1],per.tissue.age.qv$human[,1]<0.05)
age03sgn = age03
for(s in names(age03sgn)) age03sgn[[s]][age03sgn[[s]] %in% c('u','d') & per.tissue.age.qv[[s]][,colnames(age03sgn[[s]])]>0.05] = 'n'


qki.coor = c(163983952,163988000)
b=c('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/6048sTS.Human.Brain.28ypb.Male.bam',
		'/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/5531sTS.Human.Brain.29ypb.Male.bam',
		'/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/5517sTS.Human.Brain.50ypb.Male.bam',
		'/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/5524sTS.Human.Brain.58ypb.Male.bam')
qki.brain=getReadCoverage(b,'6',qki.coor[1],qki.coor[2],F,0)
b=c('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/6042sTS.Human.Heart.25ypb.Male.bam',
		'/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/6043sTS.Human.Heart.54ypb.Male.bam')
qki.hearta=getReadCoverage(b,'6',qki.coor[1],qki.coor[2],F,0)
b=c('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/3588sTS.Human.Heart.CS13.Male.bam',
		'/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/3671sTS.Human.Heart.CS16.Male.bam',
		'/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/2088sTS.Human.Heart.CS18.Male.bam')
qki.hearte=getReadCoverage(b,'6',qki.coor[1],qki.coor[2],F,0)

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

# plot
pdf('figures/paper.figures/3.2/6.pdf',w=10,h=13)
m = matrix(c(1:11,11),ncol=3,byrow = T)
layout(m)
par(tck=-0.01,mgp=c(1.1,0.2,0),mar=c(5,3,1.5,1),oma=c(0,0,2,1))
barplotWithText(hex.stat[c(1,3),],col=rep(params$tissue.col,each=2),las=3,main='upstream',beside = T,den=c(-1,30),srt = 90,adj=c(-0.1,0.5),ylim=c(0,max(hex.stat[c(1,3),])*1.2))
legend('topright',fill='black',den=c(-1,30),legend=c('inclusion','exclusion'),bty = 'n')
plotPanelLetter('A')
barplotWithText(hex.stat[c(2,4),],col=rep(params$tissue.col,each=2),las=3,main='downstream',beside = T,den=c(-1,30),srt = 90,adj=c(-0.1,0.5),ylim=c(0,max(hex.stat[c(2,4),])*1.2))
plotPanelLetter('B')
par(mgp=c(2.1,0.2,0))
barplotWithText(known.hex.stat[,2]/(known.hex.stat[,1]+known.hex.stat[,2]),t = known.hex.stat[,1]+known.hex.stat[,2],srt = 90,adj = c(-0.1,0.5),ylim=c(0,1.2),xlab='# of tissues\nwhere hexamer is sign.',yaxt='n')
par(mgp=c(1.1,0.2,0))
title(ylab='proportion of known motifs')
axis(2,seq(0,1,0.2))
plotPanelLetter('C')
par(mar=c(3,3,1.5,1))

pv.col=c(gray=1,yellow=1e-2,orange=1e-5,red=1e-20,none=0)
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
plotFTMotifSimMatrix(hex.overlap[5:8,5:8],T,main='Heart',pv.col=pv.col,or.col=or.col,diag.text = dt)
rect(0,-0.3/7,4,-1.3/7,col=params$tissue.col[3],border = NA,xpd=T)
legend(4.1,4,fill=names(pv.col)[1:4],legend = paste0('<',pv.col[1:4]),title = 'p-value',xpd=T)
legend(4.1,2,fill=names(or.col)[4:7],legend = paste0('<',or.col[4:7]),title = 'odds ratio',xpd=T)
#QKI
par(mar=c(3,3,1.5,1))
plotMirroredMotFreq(fa,age03sgn,'actaac','brain','heart',main = 'QKI motif in brain')
plotPanelLetter('E')
plotMirroredMotFreq(fa,age03sgn,'actaac','heart','brain',main = 'QKI motif in heart',plot.leg = F)
plotPanelLetter('F')

#plotTissueAgeProile(ens.ge.cod$mouse$rpkm['ENSMUSG00000062078',],meta,ylab='RPKM',ylim=range(0,ens.ge.cod$mouse$rpkm['ENSMUSG00000062078',]),main='Expression of mouse QKI')
plotTissueAgeProile(ens.ge.cod$human$rpkm['ENSG00000112531',],meta,ylab='RPKM',ylim=range(0,ens.ge.cod$human$rpkm['ENSG00000112531',]),main='Expression of human QKI',age.axis = 'rank')
plotPanelLetter('G')
plotTissueAgeProile(psi.tsm$human['hum.57513.s15',],meta.tsm,ylab='PSI',ylim=c(0,1),main='Splicing of last QKI exon',age.axis = 'rank')
plotPanelLetter('H')
plotQKICov('I')
dev.off()

