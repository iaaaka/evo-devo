#setwd('~/skoltech/projects/evo.devo/')
options(stringsAsFactors = FALSE)
source('code/r.functions/load.all.data.F.R')
source('code/r.functions/paper.figures.F.R')
library(SAJR)
library(ape)


species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
orth.seg.ad.tsm = readRDS('Rdata/orth.seg.ad.tsm.Rdata')
orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
ens.ge.cod.tsm = readRDS('Rdata/ens.ge.cod.tsm.Rdata')
orth.ens.genes = readRDS('Rdata/orth.ens.genes.Rdata')
psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
anns = readRDS('Rdata/anns.Rdata')
border.stages = readRDS('Rdata/border.stages.Rdata')
seg2ens = readRDS('Rdata/seg2ens.Rdata')
ens.ge.cod = readRDS('Rdata/ens.ge.cod.Rdata')

age.al.i = read.table('input/margarida.age.alignment.final.2017.05.17.improved.tab',sep='\t',header = T)[,c(rownames(species),'to.remove')]
params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)
#orth.seg.ad = lapply(orth.seg.ad,function(x)x[,colnames(x$ir) %in% rownames(meta)])

library(GO.db)
mouse.go = read.table('input/GO/Mus_musculus.GRCm38.84.GO.csv.gz',sep=',',header = T)
mouse.go = mouse.go[mouse.go[,2] != '',]
mouse.go = split(mouse.go$GO.Term.Accession,mouse.go$Ensembl.Gene.ID)
# GOALLANCESTOR = c(as.list(GOCCANCESTOR),as.list(GOMFANCESTOR),as.list(GOBPANCESTOR))
# mouse.go.full = lapply(mouse.go,function(x){r=unique(c(x,unlist(GOALLANCESTOR[x])));r[grep('GO:',r,fixed=T)]})
mouse.go.rev = revList(mouse.go)

############################
#### figure 1 orth MDS #####
orth.all.psi = do.call(cbind,lapply(orth.seg.ad,function(x)x$ir))
orth.all.psi = orth.all.psi[,rownames(meta)]

alt.cnt = apply(sapply(orth.seg.ad,function(x)x$seg$type=='ALT'),1,sum)
table(alt.cnt)
orth.cor = list()
orth.cor$all = cor(orth.all.psi,u='p')
orth.cor$sp1 = cor(orth.all.psi[alt.cnt==1,],u='p')
orth.cor$sp7 = cor(orth.all.psi[alt.cnt==7,],u='p')
#saveRDS(orth.cor,'Rdata/paper.figures/orth.cor.Rdata')
orth.mds2 = lapply(orth.cor,function(x)cmdscale(1-x,k=2))

pdf('figures/paper.figures/2/figure1.MDSs.pdf',w=8,h=8)
par(mfrow=c(2,2),tck=-0.01,mgp=c(1.9,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,0,1))
plot.new()
plotPanelLetter('A')
m = meta
SAJR::plotMDS(points=-orth.mds2$all,col=m$col,cex=m$cex,pch=m$pch,main=paste0('All cassette exons (',nrow(orth.seg.ad$human$seg),')'))
legend('topleft',col=params$tissue.col,pch=19,legend=names(params$tissue.col),bty='n')
plotPanelLetter('B')
SAJR::plotMDS(points=orth.mds2$sp1,col=params$species.col[m$species],cex=m$cex,pch=m$pch,main=paste0('Species-specific cassette exons (',sum(alt.cnt==1),')'))
legend('topleft',col=params$species.col,pch=19,legend=names(params$species.col),bty='n')
plotPanelLetter('С')
SAJR::plotMDS(points=orth.mds2$sp7,col=m$col,cex=m$cex,pch=m$pch,main=paste0('Ancient cassette exons (',nrow(orth.seg.ad$human$seg),')'))
plotPanelLetter('D')
dev.off()

# f = m$tissue %in% c('liver','kidney','ovary','testis')
# t = cmdscale(1-orth.cor$sp7[f,f],k=4)
# SAJR::plotMDS(points=t[,1:2],col=m$col[f],cex=m$cex[f],pch=m$pch[f],main=paste0('All cassette exons (',nrow(orth.seg.ad$human$seg),')'))
# SAJR::plotMDS(points=t[,3:4],col=m$col,cex=m$cex,pch=m$pch,main=paste0('All cassette exons (',nrow(orth.seg.ad$human$seg),')'))

##########################################
#### figure 2 mouse Age-regulation #######
library(ape)
library(TKF)
mdata = readRDS('Rdata/mouse.as.u.filtered.Rdata')
mouse.ad.cor = cor(mdata$ir[mdata$seg$sites=='ad',],u='p')
mouse.ad.cor.mds = cmdscale(1-mouse.ad.cor,k=2)
# splicing correlation
p = psi.tsm$mouse[anns$mouse$sites=='ad',]
d = sort(unique(meta$days[meta$species=='mouse']))
m = meta.tsm[colnames(p),]
tissues = unique(meta$tissue)
base = apply(p[,m$days==d[1]],1,mean,na.rm=T)
mouse.stage.cors = lapply(d,function(day){
	r = p[,m$days==day]
	colnames(r) = substr(m$tissue[m$days==day],1,1)
	r = cbind(base=base,r)
	cor(r,u='p')})

#gene expression correlation
mouse.ge = readRDS('Rdata/ens.ge.cod.tsm.Rdata')$mouse
gc()
mouse.ge = mouse.ge[,rownames(m)]
mouse.ge = mouse.ge + min(mouse.ge[mouse.ge!=0],na.rm=T)
mouse.ge =log2(mouse.ge)
base = apply(mouse.ge[,m$days==d[1]],1,mean,na.rm=T)
mouse.ge.stage.cors = lapply(d,function(day){
	r = mouse.ge[,m$days==day]
	colnames(r) = substr(m$tissue[m$days==day],1,1)
	r = cbind(base=base,r)
	cor(r,u='p')})

sfs = rownames(mouse.ge) %in% mouse.go.rev[['GO:0008380']]


mouse.sf.stage.cors = lapply(d,function(day){
	r = mouse.ge[sfs,m$days==day]
	colnames(r) = substr(m$tissue[m$days==day],1,1)
	r = cbind(base=base[sfs],r)
	cor(r,u='p')})


names(mouse.stage.cors) = names(mouse.ge.stage.cors) = names(mouse.sf.stage.cors) = d

mouse.as.stage.trees = makeNJTreesByCorLastTopology(mouse.stage.cors)
mouse.ge.stage.trees = makeNJTreesByCorLastTopology(mouse.ge.stage.cors)
mouse.sf.stage.trees = makeNJTreesByCorLastTopology(mouse.sf.stage.cors)

pdf('figures/paper.figures/2/figure2.mouse.trees.pdf',w=12,h=6)
layout(cbind(c(1,1,1),matrix(2:19,nrow=3,byrow = T)),widths = c(3,1,1,1,1,1,1))
par(tck=-0.02,mgp=c(1.1,0.2,0),mar=c(4,2,1.5,3),oma=c(0,0,1,1))
m = meta[rownames(mouse.ad.cor.mds),]
SAJR::plotMDS(points=mouse.ad.cor.mds,col=m$col,cex=m$cex*2,pch=19,main=paste0('All mouse cassette exons (',nrow(mdata$ir),')'))
par(mar=c(0.2,0,1.5,0.2))
plotPanelLetter('A')
plotMouseTrees(mouse.as.stage.trees,'Splicing','B')
plotMouseTrees(mouse.ge.stage.trees,'Gene expression','C')
plotMouseTrees(mouse.sf.stage.trees,'Gene expression of SF','D')
dev.off()

###########################
### figure 3 AS vs GE #####
# I DO NOT USE NON-CODING SEGMENTS
# as.in.ge.patterns.mouse = list()
# sp = 'mouse'
# mc = read.csv(paste0('processed/GE.from.marg/',firstToupper(sp),'Clusters.csv'),row.names = 1)
# colnames(mc) = tolower(colnames(mc))
# for(tis in unique(meta$tissue)){
# 	cat(tis)
# 	t = getsPSIbyEnsID(list(mouse=psi.tsm$mouse[anns$mouse$sites=='ad' & anns$mouse$cod!='n',]),border.stages,tis,seg2ens,sp)
# 	t = t(t)
# 	t = t[intersect(rownames(t),rownames(mc)[!is.na(mc[,paste0(tis,'pattern')])]),]
# 	as.in.ge.patterns.mouse[[tis]] = cbind(data.frame(t),ge.pattern=mc[rownames(t),paste0(tis,'pattern')])
# }
# saveRDS(as.in.ge.patterns.mouse,'Rdata/paper.figures/as.in.ge.patterns.mouse.Rdata')
#GE patterns


DPSI = 0.5
as.in.ge.patterns.mouse = readRDS('Rdata/paper.figures/as.in.ge.patterns.mouse.Rdata')
as.in.ge.patterns.mouse.cnt.up = lapply(as.in.ge.patterns.mouse,function(x){table(factor(x[,3]),x[,1]>DPSI)[,c('TRUE','FALSE')]})
as.in.ge.patterns.mouse.stat.up = sapply(as.in.ge.patterns.mouse.cnt.up,function(x){cbind(my.binom.test(x['Increasing',]),my.binom.test(x['Decreasing',]))})
as.in.ge.patterns.mouse.cnt.dw = lapply(as.in.ge.patterns.mouse,function(x){table(factor(x[,3]),x[,2]< -DPSI)[,c('TRUE','FALSE')]})
as.in.ge.patterns.mouse.stat.dw = sapply(as.in.ge.patterns.mouse.cnt.dw,function(x){cbind(my.binom.test(x['Increasing',]),my.binom.test(x['Decreasing',]))})

#Evolution of GE
ge.info.m = read.csv('/home/mazin/skoltech/projects/evo.devo/input/gene.info.from.marg/Mouse.Indexes.csv')
rownames(ge.info.m) = ge.info.m$Mouse_ID
m.age.ens = getAgeASEns(psi.tsm,meta.tsm,DPSI,border.stages,'mouse')
ohnologs.m = read.table('input/ohnologs/MOUSE.Pairs.Intermediate.2R.txt',sep='\t',header = T)


pdf(paste0('figures/paper.figures/2/figures3.GE.dPSI=',DPSI,'.pdf'),w=10,h=10/3*4)
par(mfrow=c(4,3),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(4,2,1.5,0),oma=c(0,0,1,1))

c = rep(params$tissue.col[colnames(as.in.ge.patterns.mouse.stat.up)],each=2)
b = barplot(as.in.ge.patterns.mouse.stat.up[c(1,4),],beside = T,col=c,den=c(-1,30),las=3,ylim=c(-max(as.in.ge.patterns.mouse.stat.dw[c(3,6),]),max(as.in.ge.patterns.mouse.stat.up[c(3,6),])),ylab='proportion of genes with age AS',main='AS in mouse GE clusters',yaxt='n')
segments(b,as.in.ge.patterns.mouse.stat.up[c(2,5),],b,as.in.ge.patterns.mouse.stat.up[c(3,6),])
b = barplot(-as.in.ge.patterns.mouse.stat.dw[c(1,4),],beside = T,col=c,den=c(-1,30),add=T,xaxt='n',yaxt='n')
segments(b,-as.in.ge.patterns.mouse.stat.dw[c(2,5),],b,-as.in.ge.patterns.mouse.stat.dw[c(3,6),])
abline(h=0)
at=-1:3*0.05
axis(2,at,abs(at))
legend('topright',fill='black',den=c(-1,30),legend=c('Increasing GE','Decreasing GE'))
text(b[2],par('usr')[4],'AS up',adj=c(-0.1,1.1))
text(b[2],par('usr')[3],'AS down',adj=c(-0.1,-1.1))
plotPanelLetter('A')

plotOhnologsFreq(m.age.ens,ge.info.m$Mouse_ID[!is.na(ge.info.m$Age) & ge.info.m$Age==0],main='AgeAS is enriched in ancient genes',
									 filter = intersect(ge.info.m$Mouse_ID[!is.na(ge.info.m$Age)],rownames(ens.ge.cod$mouse$gene)),ylab='Proportion of ancient genes')
plotPanelLetter('B')
plotOhnologsFreq(m.age.ens,c(ohnologs.m$Ohnolog.1.Id,ohnologs.m$Ohnolog.2.Id),main='AgeAS is enriched in ohnologs',
								 filter = intersect(ge.info.m$Mouse_ID[!is.na(ge.info.m$Age)],rownames(ens.ge.cod$mouse$gene)))
plotPanelLetter('C')
dev.off()


#######################################################
### figure 4 Species divergence on age and tissue #####
f = sapply(names(ens.ge.cod.tsm),function(s)orth.ens.genes[,s] %in% rownames(ens.ge.cod.tsm[[s]]))
f = apply(f,1,sum)==7
ens.ge.cod.tsm.log = lapply(names(ens.ge.cod.tsm),function(s)log2(ens.ge.cod.tsm[[s]][orth.ens.genes[f,s],]+0.1))
names(ens.ge.cod.tsm.log) = names(ens.ge.cod.tsm)

tissue.stage.cor = lapply(unique(meta$tissue),function(t){
	r=lapply(1:nrow(age.al.i),function(i)getSpeciesCor(orth.seg.ad.tsm,age.al.i[i,c(1:7)],t))
	names(r) = age.al.i[,'mouse']
	r
})

f = apply(sapply(orth.seg.ad,function(x)x$seg$type=='ALT'),1,sum)==7
table(f)
tissue.stage.cor.anc = lapply(unique(meta$tissue),function(t){
	r=lapply(1:nrow(age.al.i),function(i)getSpeciesCor(orth.seg.ad.tsm,age.al.i[i,c(1:7)],t,f=f))
	names(r) = age.al.i[,'mouse']
	r
})

tissue.stage.cor.anc.sp = lapply(unique(meta$tissue),function(t){
	r=lapply(1:nrow(age.al.i),function(i)getSpeciesCor(orth.seg.ad.tsm,age.al.i[i,c(1:7)],t,f=f,method = 'sp'))
	names(r) = age.al.i[,'mouse']
	r
})

tissue.stage.cor.ge = lapply(unique(meta$tissue),function(t){
	r=lapply(1:nrow(age.al.i),function(i)getSpeciesCor(ens.ge.cod.tsm.log,age.al.i[i,c(1:7)],t))
	names(r) = age.al.i[,'mouse']
	r
})


names(tissue.stage.cor) = names(tissue.stage.cor.anc) = names(tissue.stage.cor.anc.sp) = names(tissue.stage.cor.ge) = unique(meta$tissue)


pdf('figures/paper.figures/2/figure4.divergence.pdf',w=7,h=3.5)
layout(matrix(c(1,1,1,1,2,3,4,5),ncol=4))
par(tck=-0.01,mgp=c(1.9,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,2,1))
plotSpeciesDiv.fig3(lapply(tissue.stage.cor,function(x)x[-5]),c(3:6),c(-0.0,0.2),c(-0.0,0.2),0.9,main='AS divergence')
mtext('no e14.5, no human, 0.9-cor',3,outer = TRUE)
plotSpeciesDiv.fig3(lapply(tissue.stage.cor,function(x)x[-5]),c(1,3:6),c(-0.0,0.25),c(-0.0,0.25),0.9,main='AS divergence')
mtext('no e14.5, 0.9-cor',3,outer = TRUE)
plotSpeciesDiv.fig3(tissue.stage.cor,c(1,3:6),c(-0.0,0.35),c(-0.0,0.35),1,main='AS divergence')
mtext('with e14.5, 1-cor',3,outer = TRUE)
plotSpeciesDiv.fig3(tissue.stage.cor.anc,c(1,3:6),c(-0.0,0.25),c(-0.0,0.25),1,main='ancient AS divergence')
mtext('with e14.5, 1-cor',3,outer = TRUE)
plotSpeciesDiv.fig3(tissue.stage.cor.ge,c(1,3:6),c(-0.0,0.3),c(-0.0,0.3),1,main='GE divergence')
mtext('cor(log(RPKM + 0.1)), 1-cor',3,outer = TRUE)
dev.off()




