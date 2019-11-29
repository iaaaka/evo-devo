options(stringsAsFactors = FALSE)
source('code/r.functions/load.all.data.F.R')
source('code/r.functions/ad.on.ge.F.R')
source('code/r.functions/paper.figures.4.F.R')
library(SAJR)

species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
orth.seg.ad.all.id = readRDS('Rdata/orth.seg.ad.all.id.Rdata')

anns = readRDS('Rdata/anns.Rdata')
psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
my.ge.cod.tsm = readRDS('Rdata/my.ge.cod.tsm.Rdata')
my.ge.cod = readRDS('Rdata/my.ge.cod.Rdata')
ens.ge.marg.tsm = readRDS('Rdata/ens.ge.marg.tsm.Rdata')
ens.ge.marg = readRDS('Rdata/ens.ge.marg.Rdata')
ens.ge.cod = readRDS('Rdata/ens.ge.cod.Rdata')
ens.ge = readRDS('Rdata/ens.ge.Rdata')
seg2ens = readRDS('Rdata/seg2ens.Rdata')
MEL = 500
MAS = 5

anns.ad.cod = lapply(anns,function(x)x[x$length<=MEL & x$alt.size<=MAS & x$sites=='ad' & x$cod.gene, ])
psi.tsm.ad.cod = lapply(rownames(species),function(s)psi.tsm[[s]][rownames(anns.ad.cod[[s]]),])
ge.tsm.ad.cod= lapply(rownames(species),function(s)my.ge.cod.tsm[[s]][anns.ad.cod[[s]]$gene_id,])
names(ge.tsm.ad.cod) = names(psi.tsm.ad.cod) = rownames(species)
border.stages = readRDS('Rdata/border.stages.Rdata')

params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)
exon.cnt = readRDS('Rdata/exon.cnt.per.transc.Rdata')
all.anns = readRDS('Rdata/all.anns.Rdata')

# border.stages = list()
# for(s in rownames(species)){
# 	border.stages[[s]] = list()
# 	for(t in unique(meta$tissue)){
# 		m = meta.tsm[meta.tsm$species==s & meta.tsm$tissue==t,]
# 		border.stages[[s]][[t]] = m$stage[order(m$days)[c(1,nrow(m))]]
# 	}
# 	border.stages[[s]] = do.call(rbind,border.stages[[s]])
# }
# border.stages$human[c(1,2,5,7),2] = 'youngmidage'
# border.stages$human[c(3),2] = 'youngadult'
# border.stages$rabbit[,2] = '6mpb'
# border.stages$opossum[,2] = '120dpb'
# saveRDS(border.stages,'Rdata/border.stages.Rdata')

orth.seg.ad.ids = sapply(orth.seg.ad,function(x)rownames(x$seg))
f = apply(sapply(orth.seg.ad,function(x)x$seg$type)=='ALT',1,sum)==7

pdf('figures/as.on.ge/dPSI.on.exp.andFC.pearson.dPSI-exp.ancient.exons.pdf',w=13,h=28)
cor.stat = list()
for(s in rownames(species)[-2]){
	print(s)
	par(mfrow=c(8,3),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(2,2,1.5,0),oma=c(0,0,1.5,1))
	psi  = psi.tsm.ad.cod[[s]]#[rownames(psi.tsm.ad.cod[[s]]) %in% orth.seg.ad.ids[f,s],]
	rpkm = ge.tsm.ad.cod[[s]]#[rownames(psi.tsm.ad.cod[[s]]) %in% orth.seg.ad.ids[f,s],]
	
	cor.stat[[s]] = plotAllProportionOFSegOnExp(psi,rpkm,by.gene = TRUE,bins=20,stages = border.stages[[s]],cor.by.dPSI=TRUE)
	mtext(s,side = 3,outer = T)
}
range = max(unlist(sapply(cor.stat,function(x){x[,,1]})))
for(s in rownames(species)[-2]){
	plotCorStat(cor.stat[[s]],s,pv.thr = 0.01,ylim=c(-range,range))
}
mtext('pv < 0.001',outer=T) #in 
dev.off()



# as.in.ge.patterns = list()
# 
# for(sp in rownames(species)[c(1,3)]){
# 	print(sp)
# 	as.in.ge.patterns[[sp]] = list()
# 	mc = read.csv(paste0('processed/GE.from.marg/',firstToupper(sp),'Clusters.csv'),row.names = 1)
# 	colnames(mc) = tolower(colnames(mc))
# 	for(tis in unique(meta$tissue)){
# 		cat(tis)
# 		t = getsPSIbyEnsID(psi.tsm.ad.cod,border.stages,tis,seg2ens,sp)
# 		t = t(t)
# 		t = t[intersect(rownames(t),rownames(mc)[!is.na(mc[,paste0(tis,'pattern')])]),]
# 		as.in.ge.patterns[[sp]][[tis]] = cbind(data.frame(t),ge.pattern=mc[rownames(t),paste0(tis,'pattern')])
# 	}
# }
# saveRDS(as.in.ge.patterns,'Rdata/as.in.ge.patterns.Rdata')
as.in.ge.patterns = readRDS('Rdata/as.in.ge.patterns.Rdata')
pdf('figures/as.on.ge/AS.onGE.cluster.types.pdf',w=16,h=16)
cols = unique(meta[,c('tissue','col')])
cols = setNames(cols$col,cols$tissue)
par(mfrow=c(6,7),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(4,2,1.5,0),oma=c(0,0,1,1))
for(thr in c(0.5,0.4,0.3))
	for(sp in rownames(species)[c(1,3)]){
		for(tis in unique(meta$tissue)){
			t = as.in.ge.patterns[[sp]][[tis]] 
			tab = t(sapply(split(factor(pmax(t$up,-t$down) > thr,levels = c(TRUE,FALSE)), t$ge.pattern),table))
			plotBinomWithConf(tab,ylab=paste0('proportion of |dPSI|>',thr),main=paste0(sp,' ',tis,' abs (',nrow(t),')'),col=cols[tis])
		}
	}
for(thr in c(0.5,0.4,0.3))
	for(sp in rownames(species)[c(1,3)]){
		for(tis in unique(meta$tissue)){
			t = as.in.ge.patterns[[sp]][[tis]] 
			tab = t(sapply(split(factor(t$up > thr,levels = c(TRUE,FALSE)), t$ge.pattern),table))
			plotBinomWithConf(tab,ylab=paste0('proportion of dPSI>',thr),main=paste0(sp,' ',tis,' up (',nrow(t),')'),col=cols[tis])
		}
	}
for(thr in c(0.5,0.4,0.3))
	for(sp in rownames(species)[c(1,3)]){
		for(tis in unique(meta$tissue)){
			t = as.in.ge.patterns[[sp]][[tis]] 
			tab = t(sapply(split(factor(t$down < -thr,levels = c(TRUE,FALSE)), t$ge.pattern),table))
			plotBinomWithConf(tab,ylab=paste0('proportion of dPSI<-',thr),main=paste0(sp,' ',tis,' down (',nrow(t),')'),col=cols[tis])		}
	}
dev.off()



mouse.cl = read.csv(paste0('processed/GE.from.marg/MouseClusters.csv'),row.names = 1)
colnames(mouse.cl) = tolower(colnames(mouse.cl))

ge.info.m = read.csv('/home/mazin/skoltech/projects/evo.devo/input/gene.info.from.marg/Mouse.Indexes.csv')
rownames(ge.info.m) = ge.info.m$Mouse_ID
t = t(getsPSIbyEnsID(psi.tsm.ad.cod,border.stages,'brain',seg2ens,'mouse'))
t = t[rownames(t) %in% rownames(mc),]
mouse.cl$dPSI = NA
mouse.cl[rownames(t),'dPSI'] = pmax(t[,'up'],-t[,'down'])
mouse.cl$pat.spl = paste0(ifelse(is.na(mouse.cl$dPSI),'noAS',ifelse(mouse.cl$dPSI>0.3,'psi>0.3','no_change')),'\n',mouse.cl$brainpattern)
mouse.cl$age = ge.info.m[rownames(mouse.cl),'Age']
table(mouse.cl$pat.spl)
table(is.na(mouse.cl$pat.spl))
f = !is.na(mouse.cl$brainpattern)
table(f)

pdf('figures/as.on.ge/mouse.brain.ge-as.evol.age.pdf',w=7,h=6)
par(tck=-0.02,mgp=c(1.7,0.6,0),mar=c(6,3,1.5,0),oma=c(0,0,0,1))
plotBinomWithConf(t(sapply(split(mouse.cl$age[f]==0,mouse.cl$pat.spl[f]),table))[,2:1],col=rep(c('black','gray','red'),each=3),ylab='fraction of age==0',add.totals.to.lab = T)
legend('topleft',col=c('black','gray','red'),pch=19,legend=c('no cassette exons','has cassette, but |dPSI|<0.3','|dPSI| > 0.3'))
dev.off()


me = readRDS('Rdata/ens.ge.Rdata')$mouse$gene
mouse.cl$cod = me[rownames(mouse.cl),'biotype'] == 'protein_coding'
f = mouse.cl$cod 
table(f)
mouse.cl$tf = ge.info.m[rownames(mouse.cl),'TF']=='Mouse_TF'
fisher.test(table(mouse.cl$tf,mouse.cl$age==0))
fisher.test(table(mouse.cl$tf[f],mouse.cl$age[f]==0))

table(mouse.cl$tf,has.as=!is.na(mouse.cl$dPSI),cod=mouse.cl$cod)

fisher.test(table(mouse.cl$tf,!is.na(mouse.cl$dPSI)))
fisher.test(table(mouse.cl$tf[f],!is.na(mouse.cl$dPSI[f])))
fisher.test(table(mouse.cl$tf[f],mouse.cl$dPSI[f]>0.6))

plotBinomWithConf(t(sapply(split(factor(mouse.cl$tf[f]),floor(mouse.cl$dPSI[f]*10)/10),table))[,2:1])

library(GO.db)
mouse.go = read.table('input/GO/Mus_musculus.GRCm38.84.GO.csv.gz',sep=',',header = T)
mouse.go = mouse.go[mouse.go[,2] != '',]
mouse.go = split(mouse.go$GO.Term.Accession,mouse.go$Ensembl.Gene.ID)
GOALLANCESTOR = c(as.list(GOCCANCESTOR),as.list(GOMFANCESTOR),as.list(GOBPANCESTOR))
mouse.go.full = lapply(mouse.go,function(x){r=unique(c(x,unlist(GOALLANCESTOR[x])));r[grep('GO:',r,fixed=T)]})
plot(sapply(mouse.go,length),sapply(mouse.go.full,length),log='xy')
hist(sapply(mouse.go.full,length)/sapply(mouse.go,length),breaks = 0:100,xlim=c(0,25))
abline(a=0,b=1,col='red')
mouse.go.full.rev = revList(mouse.go.full)
mouse.go.rev = revList(mouse.go)

mouse.change.05=getAgeASchanges(psi.tsm,meta.tsm,0.5,border.stages,'mouse',get.dPSI=TRUE)
mouse.change.05.ad = mouse.change.05[anns$mouse$sites=='ad',]

getGeneFreqOndPSI = function(dpsi,s2e,gids,bins=10){
	dpsi = abs(dpsi)
	dpsi = dpsi[!is.na(dpsi)]
	e2s = revList(s2e[names(dpsi)])
	e2dpsi = sapply(e2s,function(sids)max(dpsi[sids]))
	r = split(names(e2dpsi) %in% gids,floor(e2dpsi*bins))
	r = t(sapply(r,function(x)table(factor(x,levels=c(TRUE,FALSE)))))
	t(apply(r,1,function(x){r = prop.test(x[1],sum(x[1:2]));c(x,freq=r$estimate,r$conf.int)}))
}

spl.freq = lapply(unique(meta$tissue),function(t)getGeneFreqOndPSI(mouse.change.05.ad[,t],seg2ens$mouse,mouse.go.full.rev[['GO:0008380']]))
rna.freq = lapply(unique(meta$tissue),function(t)getGeneFreqOndPSI(mouse.change.05.ad[,t],seg2ens$mouse,mouse.go.full.rev[['GO:0003723']]))
tfs.freq = lapply(unique(meta$tissue),function(t)getGeneFreqOndPSI(mouse.change.05.ad[,t],seg2ens$mouse,rownames(ge.info.m)[ge.info.m$TF=='Mouse_TF']))


names(tfs.freq) = names(rna.freq) = names(spl.freq) = unique(meta$tissue)

pdf('figures/freq.of.SF.TF.pdf',w=10,h=3)
par(mfrow=c(1,3),tck=-0.02,mgp=c(1.7,0.6,0),mar=c(3,3,1.5,0),oma=c(0,0,0,1))
maxy = max(sapply(spl.freq,function(x)max(x[1:10,3])))
x = 0:9/10
plot(1,t='n',xlim=c(0,1),ylim=c(0,maxy),xlab='dPSI',ylab='freq of RNA splicing genes',main='Splicing')
for(t in unique(meta$tissue)){
	lines(x,spl.freq[[t]][1:10,3],lwd=3,col=params$tissue.col[t])
#	segments(x,spl.freq[[t]][1:10,4],x,spl.freq[[t]][1:10,5],col=params$tissue.col[t])
}

maxy = max(sapply(rna.freq,function(x)max(x[1:10,3])))
x = 0:9/10
plot(1,t='n',xlim=c(0,1),ylim=c(0,maxy),xlab='dPSI',ylab='freq of RNA binding genes',main='RNA binding')
for(t in unique(meta$tissue)){
	lines(x,rna.freq[[t]][1:10,3],lwd=3,col=params$tissue.col[t])
	#	segments(x,rna.freq[[t]][1:10,4],x,rna.freq[[t]][1:10,5],col=params$tissue.col[t])
}

maxy = max(sapply(tfs.freq,function(x)max(x[1:10,3])))
x = 0:9/10
plot(1,t='n',xlim=c(0,1),ylim=c(0,maxy),xlab='dPSI',ylab='freq of TFs',main='TFs')
for(t in unique(meta$tissue)){
	lines(x,tfs.freq[[t]][1:10,3],lwd=3,col=params$tissue.col[t])
	#	segments(x,tfs.freq[[t]][1:10,4],x,tfs.freq[[t]][1:10,5],col=params$tissue.col[t])
}
dev.off()



#GO:0008380 = RNA splicing, GO:0003723 - RNA binding
mouse.cl$SF = rownames(mouse.cl) %in% mouse.go.full.rev[['GO:0008380']]
mouse.cl$rna.bind = rownames(mouse.cl) %in% mouse.go.full.rev[['GO:0003723']]

plotBinomWithConf(t(sapply(split(factor(mouse.cl$rna.bind[f]),floor(mouse.cl$dPSI[f]*10)/10),table))[,2:1])
plotBinomWithConf(t(sapply(split(factor(mouse.cl$SF[f]),floor(mouse.cl$dPSI[f]*10)/10),table))[,2:1])
plotBinomWithConf(t(sapply(split(factor(mouse.cl$tf[f]),floor(mouse.cl$dPSI[f]*10)/10),table))[,2:1])

mrpkm = readRDS('Rdata/ens.ge.Rdata')$mouse$rpkm
mrpkm = mrpkm[,colnames(mrpkm) %in% rownames(meta)]
mpsi = readRDS('Rdata/mouse.as.u.filtered.Rdata')
m = meta[colnames(mrpkm),]
mpsi = mpsi[,rownames(m)]
sf.cor = cor(log(mrpkm[intersect(rownames(mrpkm),mouse.go.full.rev[['GO:0008380']]),m$days==83]+0.1),m='p',u='p')
al.cor = cor(log(mrpkm[,m$days==83]+0.1),m='p',u='p')
as.cor = cor(log(mpsi$ir[mpsi$seg$sites=='ad',m$days==83]+0.1),m='p',u='p')
library(ape)
par(mfrow=c(2,2))
plot.phylo(nj(as.dist(1-as.cor)),tip.color=m$col[m$days==83],type = 'unrooted',main='AS')
plot.phylo(nj(as.dist(1-sf.cor)),tip.color=m$col[m$days==83],type = 'unrooted',main='SF')
plot.phylo(nj(as.dist(1-al.cor)),tip.color=m$col[m$days==83],type = 'unrooted',main='GE')



ag. = ag[!is.na(ag$dPSI),]
split(ag.$tf,number2bin,)


s = 'mouse'
anc.sp = cor(log(my.ge.cod.tsm[[s]][intersect(rownames(my.ge.cod.tsm[[s]]), anns[[s]][orth.seg.ad.ids[f,s],'gene_id']),]+0.1),u='p',m='sp')
anc.pe = cor(log(my.ge.cod.tsm[[s]][intersect(rownames(my.ge.cod.tsm[[s]]), anns[[s]][orth.seg.ad.ids[f,s],'gene_id']),]+0.1),u='p',m='p')
ort.sp = cor(log(my.ge.cod.tsm[[s]][intersect(rownames(my.ge.cod.tsm[[s]]), all.anns[[s]][orth.seg.ad.all.id[,s],'gene_id']),]+0.1),u='p',m='sp')
ort.pe = cor(log(my.ge.cod.tsm[[s]][intersect(rownames(my.ge.cod.tsm[[s]]), all.anns[[s]][orth.seg.ad.all.id[,s],'gene_id']),]+0.1),u='p',m='p')
all.sp = cor(log(my.ge.cod.tsm[[s]]+0.1),u='p',m='sp')
all.pe = cor(log(my.ge.cod.tsm[[s]]+0.1),u='p',m='p')
anc.sp.mds = cmdscale(1-anc.sp,k = 2)
anc.pe.mds = cmdscale(1-anc.pe,k = 2)
ort.pe.mds = cmdscale(1-ort.pe,k = 2)
ort.sp.mds = cmdscale(1-ort.sp,k = 2)
all.pe.mds = cmdscale(1-all.pe,k = 2)
all.sp.mds = cmdscale(1-all.sp,k = 2)

pdf('figures/as.on.ge/mouse.ge.MDS.pdf',w=6,h=9)
m = meta.tsm[rownames(anc.sp),]
par(mfrow=c(3,2),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(2,2,2.2,0),oma=c(0,0,1.5,1))
plotMDS(points=anc.sp.mds,pch=19,col=m$col,main='Gene with ancient cassettes (915)\nSpearman')
plotMDS(points=anc.pe.mds,pch=19,col=m$col,main='Gene with ancient cassettes (915)\nPearson (log)')
plotMDS(points=ort.sp.mds,pch=19,col=m$col,main='Gene with orth exons (9154)\nSpearman')
plotMDS(points=ort.pe.mds,pch=19,col=m$col,main='Gene with orth exons (9154)\nPearson (log)')
plotMDS(points=all.sp.mds,pch=19,col=m$col,main='All genes (17177)\nSpearman')
plotMDS(points=all.pe.mds,pch=19,col=m$col,main='All genes (17177)\nPearson (log)')
dev.off()




pdf('figures/as.on.ge/as2ge.div.on.age.Spearman.dpsi=0.5.pdf',w=10,h=10)
age.dpsi.thr = 0.5
cor.method = 'sp'
par(mfrow=c(7,7),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(2,1,1.2,0),oma=c(0,0,1.5,1))
ts = unique(meta$tissue)

for(s in rownames(species)[-2]){
	xlim=range(meta.tsm$age.rank[meta.tsm$species==s])
	for(t1 in ts){
		for(t2 in ts){
			if(t1 == t2){
					plot(1.1,t='n',bty='n',xaxt='n',yaxt='n',xlab='',ylab='')
					text(1,1,t1,adj=c(0.5,0.5),cex=1.4)				
				}
			else{
				z=getdPSI2lfcCorOnAge(psi.tsm.ad.cod,ge.tsm.ad.cod,t1,t2,meta.tsm,s,cor.method = cor.method,age.dpsi.thr = age.dpsi.thr)
				plot(z$age.rank,z$rho,pch=ifelse(z$p.value<5e-2,19,1),t='b',ylim=c(-0.15,0.35),ylab='',xlim=xlim,xaxt='n',col=ifelse(z$p.value<0.005,'red','black'),xlab='')
				abline(h=0,lty=2)
				axis(1,z$age.rank,z$stage)
			}
		}
	}
	mtext(s,outer=T)
}
dev.off()


ge.info.m = read.csv('/home/mazin/skoltech/projects/evo.devo/input/gene.info.from.marg/Mouse.Indexes.csv')
rownames(ge.info.m) = ge.info.m$Mouse_ID
d = psi.tsm.ad.cod$mouse[,'mouse brain 9wpb'] - psi.tsm.ad.cod$mouse[,'mouse brain 10.5']
d = d[!is.na(d)]
length(d)
all.ens = unique(unlist(seg2ens$mouse[names(d)]))
up.ens = unique(unlist(seg2ens$mouse[names(d)[d > 0.5]]))
dw.ens = unique(unlist(seg2ens$mouse[names(d)[d < -0.5]]))
ne.ens = setdiff(all.ens,c(up.ens,dw.ens))
boxplotWithSgn(list(ge.info.m[ne.ens,'BrainTau'],ge.info.m[up.ens,'BrainTau'],ge.info.m[dw.ens,'BrainTau']),notch = T)
boxplotWithSgn(list(ge.info.m[ne.ens,'TissueTau'],ge.info.m[up.ens,'TissueTau'],ge.info.m[dw.ens,'TissueTau']),notch = T)
boxplotWithSgn(list(ge.info.m[ne.ens,'Age'],ge.info.m[up.ens,'Age'],ge.info.m[dw.ens,'Age']),notch = T)
boxplotWithSgn(list(ge.info.m[ne.ens,'medianTau'],ge.info.m[up.ens,'medianTau'],ge.info.m[dw.ens,'medianTau']),notch = T)


s='mouse'
t1='brain'
t2='cerebellum'
ds = sort(unique(meta.tsm[meta.tsm$species==s,'days']))[-1:-3]
f = abs(psi.tsm.ad.cod[[s]][,paste(s,t1,'10.5')]-psi.tsm.ad.cod[[s]][,paste(s,t1,'9wpb')])>0.5 | abs(psi.tsm.ad.cod[[s]][,paste(s,t2,'10.5')]-psi.tsm.ad.cod[[s]][,paste(s,t2,'9wpb')])>0.5
f = !is.na(psi.tsm.ad.cod[[s]][,paste(s,t1,'10.5')]-psi.tsm.ad.cod[[s]][,paste(s,t1,'9wpb')]+psi.tsm.ad.cod[[s]][,paste(s,t2,'10.5')]-psi.tsm.ad.cod[[s]][,paste(s,t2,'9wpb')])
c = sapply(ds,function(d){
	d = unique(meta.tsm[meta.tsm$species==s & meta.tsm$days==d,'stage'])
	dpsi = psi.tsm.ad.cod[[s]][f,paste(s,t1,d)]-psi.tsm.ad.cod[[s]][f,paste(s,t2,d)]
	lfc = log2(ge.tsm.ad.cod[[s]][f,paste(s,t1,d)]+0.1)-log2(ge.tsm.ad.cod[[s]][f,paste(s,t2,d)]+0.1)
	cor(abs(dpsi),abs(lfc),u='p',m='sp')
	})
plot(c)
d=10.5
d=83
plotLine(abs(dpsi),abs(lfc),pch=19,col='#00000020')
abline(h=0,col='red');abline(v=0,col='red')


s = 'mouse'
m = meta.tsm[meta.tsm$species==s & meta.tsm$tissue=='brain',]
m = m[order(m$days),]
smpl = rownames(m)
t = nrow(m)

f = which((psi.tsm.ad.cod[[s]][,rownames(m)[nrow(m)]] - psi.tsm.ad.cod[[s]][,rownames(m)[1]]) > 0.5)
length(f)
l = list()
l$c.4 = sapply(f,function(i){cor(psi.tsm.ad.cod[[s]][i,smpl[(-t):(-t+3)]],ge.tsm.ad.cod[[s]][i,smpl[-1:-4]],u='p',m='sp')})
l$c.3 = sapply(f,function(i){cor(psi.tsm.ad.cod[[s]][i,smpl[(-t):(-t+2)]],ge.tsm.ad.cod[[s]][i,smpl[-1:-3]],u='p',m='sp')})
l$c.2 = sapply(f,function(i){cor(psi.tsm.ad.cod[[s]][i,smpl[(-t):(-t+1)]],ge.tsm.ad.cod[[s]][i,smpl[-1:-2]],u='p',m='sp')})
l$c.1 = sapply(f,function(i){cor(psi.tsm.ad.cod[[s]][i,smpl[-t]],ge.tsm.ad.cod[[s]][i,smpl[-1]],u='p',m='sp')})
l$c0 = sapply(f,function(i){cor(psi.tsm.ad.cod[[s]][i,smpl],ge.tsm.ad.cod[[s]][i,smpl],u='p',m='sp')})
l$c1 = sapply(f,function(i){cor(psi.tsm.ad.cod[[s]][i,smpl[-1]],ge.tsm.ad.cod[[s]][i,smpl[-t]],u='p',m='sp')})
l$c2 = sapply(f,function(i){cor(psi.tsm.ad.cod[[s]][i,smpl[-1:-2]],ge.tsm.ad.cod[[s]][i,smpl[(-t):(-t+1)]],u='p',m='sp')})
l$c3 = sapply(f,function(i){cor(psi.tsm.ad.cod[[s]][i,smpl[-1:-3]],ge.tsm.ad.cod[[s]][i,smpl[(-t):(-t+2)]],u='p',m='sp')})
l$c4 = sapply(f,function(i){cor(psi.tsm.ad.cod[[s]][i,smpl[-1:-4]],ge.tsm.ad.cod[[s]][i,smpl[(-t):(-t+3)]],u='p',m='sp')})
l$c5 = sapply(f,function(i){cor(psi.tsm.ad.cod[[s]][i,smpl[-1:-5]],ge.tsm.ad.cod[[s]][i,smpl[(-t):(-t+4)]],u='p',m='sp')})
l$c6 = sapply(f,function(i){cor(psi.tsm.ad.cod[[s]][i,smpl[-1:-6]],ge.tsm.ad.cod[[s]][i,smpl[(-t):(-t+5)]],u='p',m='sp')})
l$c7 = sapply(f,function(i){cor(psi.tsm.ad.cod[[s]][i,smpl[-1:-7]],ge.tsm.ad.cod[[s]][i,smpl[(-t):(-t+6)]],u='p',m='sp')})


boxplot(l)
hist(l$c.4)
plot(sapply(l,function(x)mean(x> 0.9,na.rm=T)))


ge = ge.tsm.ad.cod[[s]]
rownames(ge) = rownames(psi.tsm.ad.cod[[s]])
plotTissueAgeProile(apply(psi.tsm.ad.cod[[s]][names(l$c0)[l$c0>0.9],smpl],2,mean,na.rm=T),meta.tsm)
plotTissueAgeProile(apply(log2(ge[names(l$c0)[l$c0>0.9],smpl]),2,mean,na.rm=T),meta.tsm)

v=apply(psi.tsm.ad.cod[[s]][names(l$c0)[l$c0<0.9],smpl],2,mean,na.rm=T)
plotTissueAgeProile((v-min(v))/(max(v)-min(v)),meta.tsm,ylim=c(0,1),age.axis = 'rank')
v=apply(log2(ge[names(l$c0)[l$c0<0.9],smpl]),2,mean,na.rm=T)
plotTissueAgeProile((v-min(v))/(max(v)-min(v)),meta.tsm,add=T,col='red',age.axis = 'rank')


v=apply(psi.tsm.ad.cod[[s]][f,smpl],2,mean,na.rm=T)
plotTissueAgeProile((v-min(v))/(max(v)-min(v)),meta.tsm,ylim=c(0,1),age.axis = 'rank')
v=apply(log2(ge[f,smpl]+0.1),2,mean,na.rm=T)
v=apply(normRows(ge[f,smpl]+0.1),2,mean,na.rm=T)
plotTissueAgeProile((v-min(v))/(max(v)-min(v)),meta.tsm,add=T,col='red',age.axis = 'rank')

f = which((psi.tsm.ad.cod[[s]][,rownames(m)[nrow(m)]] - psi.tsm.ad.cod[[s]][,rownames(m)[1]]) > 0.5 & (ge.tsm.ad.cod[[s]][,rownames(m)[nrow(m)]] > ge.tsm.ad.cod[[s]][,rownames(m)[1]]))
length(f)

z=t(sapply(f,function(i){
	p = psi.tsm.ad.cod[[s]][i,smpl]
	e =  (ge.tsm.ad.cod[[s]][i,smpl]+0.1)
	(e-min(e))/(max(e)-min(e)) - (p-min(p))/(max(p)-min(p))
	}))
boxplot(z)

s='mouse'
f = which((psi.tsm.ad.cod[[s]][,rownames(m)[nrow(m)]] - psi.tsm.ad.cod[[s]][,rownames(m)[1]]) < -0.5)
psi = psi.tsm.ad.cod[[s]][f,smpl]
exp = my.ge.cod.tsm[[s]][unique(anns[[s]][names(f),'gene_id']),smpl]
psi.cl = reorderClustersBySize(cutree(hclust(as.dist(1-cor(t(psi),m='sp',u='p'))),k=4))
exp.cl = reorderClustersBySize(cutree(hclust(as.dist(1-cor(t(exp),m='sp',u='p'))),k=4))
table(exp.cl)
exp.cl = exp.cl[anns[[s]][names(f),'gene_id']]
table(psi.cl)
table(psi.cl,exp.cl)
par(mfrow=c(4,4),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(2,2,2.2,0),oma=c(0,0,1.5,1))
for(i in 1:4){
	for(j in 1:4){
		sid = names(psi.cl)[psi.cl==i & exp.cl==j]
		gid = unique(anns[[s]][sid,'gene_id'])
		if(length(sid)==0)
			plot.new()
		else{
			v=-apply(psi.tsm.ad.cod[[s]][sid,smpl,drop=F],2,mean,na.rm=T)
			plotTissueAgeProile((v-min(v))/(max(v)-min(v)),meta.tsm,ylim=c(0,1),age.axis = 'rank',main=paste0(i,':',j,' (',length(sid),':',length(gid),')'))
			v=apply(log2(my.ge.cod.tsm[[s]][gid,smpl,drop=F]+0.1),2,mean,na.rm=T)
			plotTissueAgeProile((v-min(v))/(max(v)-min(v)),meta.tsm,add=T,col='red',age.axis = 'rank')
		}
	}
}


# there were many attempts to make this analysis, the newest are on top
# devAS vs TissueTau ######
human.tau = read.csv('input/gene.info.from.marg/Human.Indexes.All.csv')
mouse.tau = read.csv('input/gene.info.from.marg/Mouse.Indexes.csv')
rownames(human.tau) = human.tau[,1]
rownames(mouse.tau) = mouse.tau[,1]
table(substr(mouse.tau$Mouse_ID,1,1),mouse.tau$Mouse_ID %in% rownames(ens.ge$mouse$gene))
table(substr(human.tau$Human_ID,1,1),human.tau$Human_ID %in% rownames(ens.ge$human$gene))
mouse.tau$biotype = setNames(ens.ge$mouse$gene$biotype,rownames(ens.ge$mouse$gene))[mouse.tau$Mouse_ID]
human.tau$biotype = setNames(ens.ge$human$gene$biotype,rownames(ens.ge$human$gene))[human.tau$Human_ID]

psi.tsm.ad = lapply(names(psi.tsm),function(s)psi.tsm[[s]][anns[[s]]$sites=='ad',])
names(psi.tsm.ad) = names(psi.tsm)

mouse.br.gene.dpsi=t(getsPSIbyEnsID(psi.tsm.ad,border.stages,'brain',seg2ens,'mouse',use.random = F))
human.br.gene.dpsi=t(getsPSIbyEnsID(psi.tsm.ad,border.stages,'brain',seg2ens,'human',use.random = F))

mouse.te.gene.dpsi=t(getsPSIbyEnsID(psi.tsm.ad,border.stages,'testis',seg2ens,'mouse',use.random = F))
human.te.gene.dpsi=t(getsPSIbyEnsID(psi.tsm.ad,border.stages,'testis',seg2ens,'human',use.random = F))


f1 = function(tt,tt.thr=0.7,...){
	h=hist(tt,seq(0,1,0.05),xlab='TissueTau',...)
	abline(v=tt.thr,col='red',lty=2)
	text(1,max(h$counts),paste0(round(mean(tt>tt.thr,na.rm=T)*100,1),'%'),col='red',adj=c(1,0))
}

f2 = function(tau,dpsi,psi.thr,sp){
	f1(tau[!(rownames(tau) %in% rownames(dpsi)),'TissueTau'],main=paste0(sp,'; no AS'))
	f1(tau[rownames(dpsi)[apply(abs(dpsi)< psi.thr,1,sum)>0],'TissueTau'],main=paste0(sp,'; has AS, devAS<',psi.thr))
	f1(tau[rownames(dpsi)[apply(abs(dpsi)>=psi.thr,1,sum)>0],'TissueTau'],main=paste0(sp,'; devAS>',psi.thr))
}

pdf('figures/as.on.ge/tissue.tau/tissue.tau.distr.brain.human-mouse.pdf',w=9,h=12)
par(mfrow=c(4,3),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(2,2,2.2,0),oma=c(0,0,1.5,1))
f2(human.tau,human.br.gene.dpsi,0.2,'human')
f2(mouse.tau,mouse.br.gene.dpsi,0.2,'mouse')
f2(human.tau,human.br.gene.dpsi,0.5,'human')
f2(mouse.tau,mouse.br.gene.dpsi,0.5,'mouse')
dev.off()

f3 = function(ge.tau,dpsi,dpsi.thr=0.2,cols=c('darkgray','lightgray','orange'),norm.by.col=F,...){
	ge.tau$as = ifelse(rownames(ge.tau) %in% rownames(dpsi),'y','!AS')
	ge.tau$as[ge.tau$as=='y'] = ifelse(apply(abs(dpsi[rownames(ge.tau)[ge.tau$as=='y'],]),1,max)>dpsi.thr,'devAS','!devAS')
	ge.tau$TissueTau[!is.na(ge.tau$TissueTau) & ge.tau$TissueTau==1] = max(ge.tau$TissueTau[!is.na(ge.tau$TissueTau) & ge.tau$TissueTau<1])
	r=table(ge.tau$as,floor(ge.tau$TissueTau*10)/10)[c('!AS','!devAS','devAS'),]
	if(norm.by.col)
		r = sweep(r,2,apply(r,2,sum),'/')
	barplot(r,names.arg = seq(0.05,0.95,by=0.1)[1:ncol(r)],col=cols,ylab='# genes',xlab='TissueTau',...)
}

pdf('figures/as.on.ge/tissue.tau/tissue.tau.distr.on.AS.brain-testis.human-mouse_.pdf',w=12,h=12)
par(mfrow=c(4,4),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(2,2,2.2,0),oma=c(0,0,1.5,1))
hc = !is.na(human.tau$biotype) & human.tau$biotype=='protein_coding'
mc = !is.na(mouse.tau$biotype) & mouse.tau$biotype=='protein_coding'
f3(human.tau,human.br.gene.dpsi,0.2,main='human brain, all',legend.text=T,args.legend=list(x='topleft'))
f3(human.tau[hc,],human.br.gene.dpsi,0.2,main='human  brain, protein-coding')
f3(human.tau[hc & human.tau$TissueMax!='TestisMax',],human.br.gene.dpsi,0.2,main='human brain, not TestisMax')
f3(human.tau[hc & human.tau$TissueMax=='BrainMax',] ,human.br.gene.dpsi,0.2,main='human brain, brainMax')

f3(human.tau,human.te.gene.dpsi,0.2,main='human testis, all')
f3(human.tau[hc,],human.te.gene.dpsi,0.2,main='human  testis, protein-coding')
f3(human.tau[hc & human.tau$TissueMax!='TestisMax',],human.te.gene.dpsi,0.2,main='human testis, not TestisMax')
f3(human.tau[hc & human.tau$TissueMax=='TestisMax',] ,human.te.gene.dpsi,0.2,main='human testis, TestisMax')

f3(mouse.tau,mouse.br.gene.dpsi,0.2,main='mouse brain, all')
f3(mouse.tau[mc,],mouse.br.gene.dpsi,0.2,main='mouse  brain, protein-coding')
f3(mouse.tau[mc & mouse.tau$TissueMax!='TestisMax',],mouse.br.gene.dpsi,0.2,main='mouse brain, not TestisMax')
f3(mouse.tau[mc & mouse.tau$TissueMax=='BrainMax',] ,mouse.br.gene.dpsi,0.2,main='mouse brain, brainMax')

f3(mouse.tau,mouse.te.gene.dpsi,0.2,main='mouse testis, all')
f3(mouse.tau[mc,],mouse.te.gene.dpsi,0.2,main='mouse  testis, protein-coding')
f3(mouse.tau[mc & mouse.tau$TissueMax!='TestisMax',],mouse.te.gene.dpsi,0.2,main='mouse testis, not TestisMax')
f3(mouse.tau[mc & mouse.tau$TissueMax=='TestisMax',] ,mouse.te.gene.dpsi,0.2,main='mouse testis, TestisMax')
mtext('|dPSI| > 0.2',3,outer=T)
dev.off()

# look on tau-expression-no-of-exon interaction
h.as.ens = loadSAData('processed/annotation/all.species/ensambl/Homo_sapiens.GRCh37.73.sajr')
h.as.ens = setSplSiteTypes(h.as.ens,'processed/annotation/all.species/ensambl/Homo_sapiens.GRCh37.73.sajr')
h.int.ex.cnt = sapply(split(h.as.ens$seg$sites,h.as.ens$seg$gene_id),function(x){sum(x %in% c('.d','a.','ad'))})
f = rownames(human.tau) %in% names(h.int.ex.cnt)
table(f)
human.tau$internal.exon.count = NA
human.tau[f,'internal.exon.count'] = h.int.ex.cnt[rownames(human.tau)[f]]
table(rownames(ens.ge.marg$human)==rownames(human.tau))
ht = calcMeanCols(log2(ens.ge.marg$human+1e-3),meta[colnames(ens.ge.marg$human),'tissue'])
ht[1:10,]
plot(ht[,1],human.tau$BrainMax)

f = !is.na(human.tau$biotype) & human.tau$biotype=='protein_coding' & human.tau$TissueMax=='BrainMax' & ht[,1] > 3 & !is.na(human.tau$internal.exon.count) & human.tau$internal.exon.count > 2
tt = floor(ifelse(human.tau$TissueTau==1,0.99,human.tau$TissueTau)*10)
barplot(table(tt[f]))
boxplot(human.tau$internal.exon.count[f] ~ tt[f],outline=F)
boxplot(human.tau$BrainMax[f] ~ tt[f],outline=F)
boxplot(ht[f,1] ~ tt[f],outline=F)
f3(human.tau[f,] ,human.br.gene.dpsi,0.2,main='human brain, brainMax',norm.by.col = T)

human.tau$has.as = (human.tau$Human_ID %in% rownames(human.br.gene.dpsi)) + 0
human.tau$BrainMean = ht[,1]

f = !is.na(human.tau$biotype) & human.tau$biotype=='protein_coding'
m = glm(has.as ~ internal.exon.count+I(internal.exon.count^2)+BrainMean+I(BrainMean^2) + TissueMax + TissueTau,data = human.tau[f,],family='quasibinomial')
m
anova(m,test='Chisq')

# exon.cnt. = list.files('processed/annotation/all.species/ensambl/','gtf.gz',full.names = T)
# exon.cnt = lapply(exon.cnt.,getExonNumberByGene)
# names(exon.cnt) = c('chicken','human','macaque','opossum','mouse','rabbit','rat')
# exon.cnt = exon.cnt[rownames(species)]
# saveRDS(exon.cnt,'Rdata/exon.cnt.per.transc.Rdata')


table(mouse.tau$Mouse_ID %in% names(mouse.exon.cnt),substr(mouse.tau$Mouse_ID,1,3))


mouse.exon.cnt = exon.cnt$mouse

image(log(1+table(sapply(mouse.exon.cnt, length),sapply(mouse.exon.cnt, max)))[1:30,1:50])
f = intersect(rownames(mouse.tau),names(mouse.exon.cnt)[sapply(mouse.exon.cnt, max)==15])
length(f)
f3(mouse.tau[f,],mouse.br.gene.dpsi,0.2,main='mouse testis 20 exons, all')
boxplot(mouse.tau[f,'ExpressionMax'] ~ round(mouse.tau[f,'TissueTau'],1))

plot(mouse.tau$TissueTau,mouse.tau$ExpressionMax)
f = !is.na(mouse.tau$CDS_Begin)
boxplot(mouse.tau$ExpressionMax[f] ~ round(mouse.tau$TissueTau[f],1))
boxplot(sapply(mouse.exon.cnt[mouse.tau$Mouse_ID[f]],max) ~ round(mouse.tau$TissueTau[f],1),outline=F)
barplot(table(sapply(mouse.exon.cnt[mouse.tau$Mouse_ID[f]],max)>15 , floor(mouse.tau$TissueTau[f]*10)/10))

hist(sapply(mouse.exon.cnt[mouse.tau$Mouse_ID[f]],max)[mouse.tau$TissueTau[f]<0.7],0:200,xlim=c(0,50),col='#FF000060',border=NA,freq=F)
hist(sapply(mouse.exon.cnt[mouse.tau$Mouse_ID[f]],max)[mouse.tau$TissueTau[f]>=0.7],0:200,xlim=c(0,50),col='#0000FF60',border=NA,add=T,freq=F)


# dPSI > 0.2 is indeed bit enriched with high tissueTau
#plot(sapply(split(apply(abs(mouse.br.gene.dpsi)>0.5,1,sum)>0,number2bin(mouse.tau[rownames(mouse.br.gene.dpsi),'TissueTau'],10)),mean))



hb.max.tis = human.tau[rownames(human.br.gene.dpsi),'TissueMax']
hb.max.tis[human.tau[rownames(human.br.gene.dpsi),'TissueTau']<0.5] = 'abig'

par(mfrow=c(2,2),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(2,2,2.2,0),oma=c(0,0,1.5,1))
t = table(hb.max.tis,apply(abs(human.br.gene.dpsi)>=0.2,1,sum)>0)
barplot(sweep(t,2,apply(t,2,sum),'/'),col=c('gray',params$tissue.col),border = NA,main='Human devAS>0.2, in tissue-spec')
barplot(sweep(t,2,apply(t,2,sum),'/')[-1,],col=c(params$tissue.col),border = NA,main='Human devAS>0.2, in tissue-spec')


t = table(mb.max.tis,apply(abs(mouse.br.gene.dpsi)>=0.2,1,sum)>0)
barplot(sweep(t,2,apply(t,2,sum),'/'),col=c('gray',params$tissue.col),border = NA,main='Mouse devAS>0.2, in tissue-spec')
barplot(sweep(t,2,apply(t,2,sum),'/')[-1,],col=c(params$tissue.col),border = NA,main='Mouse devAS>0.2, in tissue-spec')



age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)

mouse.tau$hasCE = mouse.tau$Mouse_ID    %in% unique(unlist(seg2ens$mouse[rownames(all.anns$mouse)[all.anns$mouse$sites=='ad' & all.anns$mouse$type=='ALT']]))
mouse.tau$hasDevAS = mouse.tau$Mouse_ID %in% unique(unlist(seg2ens$mouse[rownames(anns$mouse    )[anns$mouse$sites=='ad' & apply(abs(age.dpsi$mouse)>0.2,1,sum)>0]]))
mouse.tau$maxExonCnt = sapply(exon.cnt$mouse,max)[mouse.tau$Mouse_ID]
mouse.tau$as = ifelse(mouse.tau$hasCE,ifelse(mouse.tau$hasDevAS,'devAS','!devAS'),'!AS')
mt = mouse.tau[!is.na(mouse.tau$maxExonCnt) & !is.na(mouse.tau$biotype) & mouse.tau$biotype=='protein_coding',] 
mt$trans.cnt = sapply(exon.cnt$mouse,length)[mt$Mouse_ID]
dim(mt)
table(mt$as)

f = mt$maxExonCnt==15
barplot(table(mt$as[f],floor(mt$TissueTau[f]*10)/10),col=c('darkgray','lightgray','orange'))

mt$brain.bin = number2bin(mt$BrainMax,3)
mt$heart.bin = number2bin(mt$HeartMax,3)
t = cast(mt[mt$as!='!AS' & mt$maxExonCnt==15,],brain.bin ~ heart.bin,value = 'as',fun.aggregate = function(x)mean(x=='devAS'))
t = cast(mt[mt$maxExonCnt==15,],brain.bin ~ heart.bin,value = 'as',fun.aggregate = function(x)mean(x!='!AS'))
image(as.matrix(t))

mouse.br.gene.dpsi
mbg = as.data.frame(mouse.br.gene.dpsi[rownames(mouse.br.gene.dpsi) %in% rownames(mouse.tau)[mouse.tau$biotype=='protein_coding'],]) 
mbg$maxExonCnt = sapply(exon.cnt$mouse,max)[rownames(mbg)]
mbg$brain.bin = number2bin(mouse.tau[rownames(mbg),'BrainMax'],3)
mbg$heart.bin = number2bin(mouse.tau[rownames(mbg),'HeartMax'],3)
mbg$liver.bin = number2bin(mouse.tau[rownames(mbg),'LiverMax'],3)
mbg$adPSI = pmax(mbg$up,-mbg$down)

t = cast(mbg[mbg$maxExonCnt==9,],brain.bin ~ liver.bin,value = 'adPSI',fun.aggregate = function(x)mean(x>0.2))
image(as.matrix(t))

cor.test(mt$trans.cnt,mt$maxExonCnt,m='sp')
image(table(trans.cnt=mt$trans.cnt,mt$maxExonCnt)[1:30,1:30])

t = cast(mt[mt$maxExonCnt==9,],brain.bin ~ heart.bin,value = 'trans.cnt',fun.aggregate = function(x)mean(x>1))
t = cast(mt[mt$maxExonCnt==9,],brain.bin ~ heart.bin,value = 'trans.cnt',fun.aggregate = function(x)sum(x>1))
image(as.matrix(t))
# check is brain devAS depends on heart expression ####
# look on exon level 
library(reshape)
mouse.br.gene.dpsi
m = readRDS('Rdata/mouse.as.u.all.Rdata')
m = m[m$seg$sites=='ad',]
t = seg2ens$mouse[rownames(m$seg)]
table(sapply(t,length))
m = m[sapply(t,length)==1,]
mm = meta[colnames(m$ir),]

all.ad.br.dpsi = data.frame(sid = rownames(m$seg),ens.id = unlist(seg2ens$mouse[rownames(m$seg)]),dPSI.br = apply(m$ir[,mm$tissue=='brain' & mm$stage==border.stages$mouse['brain',2]],1,mean,na.rm=T) - 
					 																																								 									apply(m$ir[,mm$tissue=='brain' & mm$stage==border.stages$mouse['brain',1]],1,mean,na.rm=T))

all.ad.br.dpsi = all.ad.br.dpsi[!is.na(all.ad.br.dpsi$dPSI.br),]
dim(all.ad.br.dpsi)
all.ad.br.dpsi = all.ad.br.dpsi[all.ad.br.dpsi$ens.id %in% rownames(ens.ge.marg.tsm$mouse),]
m.ge.tissue.mean = calcMeanCols(ens.ge.marg.tsm$mouse,meta.tsm[colnames(ens.ge.marg.tsm$mouse),'tissue'])
all.ad.br.dpsi = cbind(all.ad.br.dpsi,m.ge.tissue.mean[all.ad.br.dpsi$ens.id,])
all.ad.br.dpsi$brain.bin = number2bin(all.ad.br.dpsi$brain,10)
all.ad.br.dpsi$heart.bin = number2bin(all.ad.br.dpsi$heart,11)
all.ad.br.dpsi$liver.bin = number2bin(all.ad.br.dpsi$liver,12)
plot(sapply(split(abs(all.ad.br.dpsi$dPSI.br)>0.5, all.ad.br.dpsi$brain.bin),mean))
plot(sapply(split(abs(all.ad.br.dpsi$dPSI.br)>0.5, all.ad.br.dpsi$heart.bin),mean))
t02 = cast(all.ad.br.dpsi,brain.bin ~ liver.bin,value = 'dPSI.br',fun.aggregate = function(x)sum(abs(x)>0.2))
image(as.matrix(t02))

plot(sapply(split(abs(all.ad.br.dpsi$dPSI.br)>0.5 ,  number2bin(all.ad.br.dpsi$brain/all.ad.br.dpsi$heart,20)),mean))
# look in gene level


## devAS in tissue-spec GE ####
sp = 'human'
tissue='brain'
t=list(psi.tsm[[sp]][anns[[sp]]$sites=='ad',])
names(t) = sp
mb=t(getsPSIbyEnsID(t,border.stages,tissue,seg2ens,sp,use.random = T))
if(sp=='human'){
	tis.spec.ge = read.csv('input/gene.info.from.marg/Human.Indexes.All.csv')
}else
	tis.spec.ge = read.csv('input/gene.info.from.marg/Mouse.Indexes.csv')
rownames(tis.spec.ge) = tis.spec.ge[,1]
table(rownames(mb) %in% rownames(tis.spec.ge))
mb = mb[rownames(mb) %in% rownames(tis.spec.ge),]
ts = tis.spec.ge[rownames(mb),]
if(sp=='human')
	ts = ts[,-2]
ts$TissueMaxNoCrb = apply(ts[,c(2,4:8)],1,function(x){o=order(x,decreasing = T);if(x[o[1]]==x[o[2]]){return(NA)};colnames(ts)[c(2,4:8)][o[1]]})


pdf('figures/as.on.ge/human.brain-sp.AS.vs.brain-sp.GE.random.exon.pdf',w=9,h=6)
par(mfrow=c(2,3),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(2,2,2.2,0),oma=c(0,0,1.5,1))
f = function(dPSI){
	devAS = pmax(mb[,'up'],-mb[,'down']) > dPSI
	#devAS = mb[,'up'] > dPSI
	b = number2bin(ts$TissueTau,10)
	z=t(apply(table(b,devAS),1,function(x)c(x,my.binom.test(x[2:1]))))
	plotArea(1:10,z[,3:5],col=params$tissue.col[tissue],new=T,xlab='TissueTau bin',ylab=paste0('proportion of genes with dPSI > ',dPSI),main=paste(sp, tissue))
	
	b = number2bin(ts[[paste0(firstToupper(tissue),'Max')]],10)
	z=t(apply(table(b,devAS),1,function(x)c(x,my.binom.test(x[2:1]))))
	plotArea(1:10,z[,3:5],col=params$tissue.col[tissue],new=T,xlab=paste0(firstToupper(tissue),'Max bin'),ylab=paste0('proportion of genes with dPSI > ',dPSI),main=paste(sp, tissue))
	
	z = table(devAS,ts$TissueMaxNoCrb)
	z = t(apply(t(z[,order(-z[2,]/apply(z,2,sum))]),1,function(x)c(x,my.binom.test(x[2:1]))))
	rownames(z) = substr(rownames(z),1,1)
	
	b=barplot(z[,3],ylim=range(0,z[,3:5]),xlab='TissueMax (without cerebellum)',ylab=paste0('proportion of genes with dPSI > ',dPSI),main=paste(sp, tissue))
	segments(b,z[,4],b,z[,5])
}
f(0.2)
f(0.5)
dev.off()

b = number2bin(ts$BrainMax-ts$ExpressionMax,10)


plotLine(ts$TissueTau,pmax(-mb[,'down'],mb[,'up']))
plotLine(ts$TissueTau,mb[,'up'])
z=cbind(number2bin(ts$BrainMax,10),number2bin(ts$HeartMax,10))
f = mb[,'up'] > 0.3
imageWithText(table(br=z[f,1],hr=z[f,2])/table(br=z[,1],hr=z[,2]),names.as.labs = T)

plot(ts$BrainMax,ts$HeartMax)
f=mb[,'up']>0.5
points(ts$BrainMax[f],ts$HeartMax[f],col='red',pch=19)

table(max2=ts$Tissue2MaxNoCrb,ts$TissueMaxNoCrb)
table(ts$TissueMax,ts$TissueMaxNoCrb)



# double check with devAS enrichment in late genes
ge.patt = read.csv('processed/GE.from.marg/MouseClusters.csv',row.names = 1)
ge.patt = ge.patt[rownames(mb),]
z=table(mb[,'up']>0.8 , ge.patt$BrainPattern)
barplot(z[2,]/apply(z,2,sum))

boxplotWithSgn(split(ts$BrainMax , ge.patt$BrainPattern))


# brain-heart diff AS: is it in genes expressed in both tissues? ####
border.stages$mouse
bh = psi.tsm$mouse[,'mouse brain 9wpb'] - psi.tsm$mouse[,'mouse heart 9wpb']
bh = bh[anns$mouse$sites=='ad']
hist(bh)
bh.ge = log2(my.ge.cod.tsm$mouse[,'mouse brain 9wpb'] / my.ge.cod.tsm$mouse[,'mouse heart 9wpb'])
bh.ge = bh.ge[anns$mouse[names(bh),'gene_id']]
cor.test(abs(bh.ge),abs(bh),m='sp')
plotArea(1:10,t(apply(table(number2bin(abs(bh.ge),10),abs(bh)>0.5)[,2:1],1,my.binom.test)),col='red',new=T)

# brain devAS in broadly expressed genes  #####
sp = 'mouse'
t=list(psi.tsm[[sp]][anns[[sp]]$sites=='ad',])
names(t) = sp
g.dpsi = lapply(unique(meta$tissue),function(tissue)t(getsPSIbyEnsID(t,border.stages,tissue,seg2ens,sp,use.random = F)))
names(g.dpsi) = unique(meta$tissue)
#g.dpsi.random = g.dpsi

if(sp=='human'){
	tis.spec.ge = read.csv('input/gene.info.from.marg/Human.Indexes.All.csv')
}else
	tis.spec.ge = read.csv('input/gene.info.from.marg/Mouse.Indexes.csv')
rownames(tis.spec.ge) = tis.spec.ge[,1]
dsi.gids = table(unlist(lapply(g.dpsi,rownames)))
table(dsi.gids)

cmn = intersect(names(dsi.gids),rownames(tis.spec.ge))
length(cmn)

f = function(dp,gids,dpsi.thr,...){
	z=sapply(g.dpsi,function(x)x[gids,'up'])
	z=t(apply(z>dpsi.thr,2,function(x){s=sum(x);f=sum(!x);c(s,f,my.binom.test(s,f))}))
	b=barplot(z[,3],ylim=range(0,z[,3:5]),xlab='tissue',ylab='# genes with devAS',names.arg = substr(rownames(z),1,1),...)
	segments(b,z[,4],b,z[,5])
}

pdf('figures/as.on.ge/mouse.devAS.in.broadly-expressed-genes.pdf',w=9,h=6)
par(mfrow=c(2,3),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(2,2,2.2,0),oma=c(0,0,1.5,1))
boxplot(tis.spec.ge[cmn,'TissueTau'] ~ dsi.gids[cmn],xlab='# of tissues with detectable splicing',ylab='TissueTau')
f(g.dpsi,names(dsi.gids)[dsi.gids==7],0.2,main='AS detected in 7 tissues; dPSI>0.2')
f(g.dpsi,names(dsi.gids)[dsi.gids==7],0.5,main='AS detected in 7 tissues; dPSI>0.5')
plot.new()
f(g.dpsi,names(dsi.gids)[dsi.gids==7 & names(dsi.gids) %in% rownames(tis.spec.ge)[tis.spec.ge$TissueTau<0.2]],0.2,main='AS detected in 7 tissues,\nTissueTau<0.2; dPSI>0.2')
f(g.dpsi,names(dsi.gids)[dsi.gids==7 & names(dsi.gids) %in% rownames(tis.spec.ge)[tis.spec.ge$TissueTau<0.2]],0.5,main='AS detected in 7 tissues,\nTissueTau<0.2; dPSI>0.5')
dev.off()