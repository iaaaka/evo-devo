#setwd('~/skoltech/projects/evo.devo/')
options(stringsAsFactors = FALSE)
source('code/r.functions/load.all.data.F.R')
source('code/r.functions/as.on.age.F.R')
source('code/r.functions/ad.on.ge.F.R')
source('code/r.functions/paper.figures.F.R')
library(SAJR)
library(GenomicRanges)

species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
anns = readRDS('Rdata/anns.Rdata')
all.anns = readRDS('Rdata/all.anns.Rdata') 
psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
border.stages = readRDS('Rdata/border.stages.Rdata')
seg2ens = readRDS('Rdata/seg2ens.Rdata')
ens.ge.cod = readRDS('Rdata/ens.ge.cod.Rdata')
orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
orth.seg.ad.tsm = readRDS('Rdata/orth.seg.ad.tsm.Rdata')

params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)

# microexons
# figure 1: numbers of events ####


getNoOfEnsGenes = function(a,s2e,cod=c('c','n','p')){
	sapply(names(a),function(s)length(intersect(rownames(ens.ge.cod[[s]]$gene),unique(unlist(s2e[[s]][rownames(a[[s]])[a[[s]]$type != 'EXN' & a[[s]]$cod %in% cod]])))))
}

all.seg.cnt = sapply(all.anns,getNoOfEvents,gene=FALSE)
filtered.seg.cnt = sapply(anns,getNoOfEvents,gene=FALSE)

all.ens.cnt = rbind(cod=getNoOfEnsGenes(all.anns,seg2ens,'c'),all=getNoOfEnsGenes(all.anns,seg2ens))
all.ens.cnt[2,] = all.ens.cnt[2,] - all.ens.cnt[1,]

filtered.ens.cnt = rbind(cod=getNoOfEnsGenes(anns,seg2ens,'c'),all=getNoOfEnsGenes(anns,seg2ens))
filtered.ens.cnt[2,] = filtered.ens.cnt[2,] - filtered.ens.cnt[1,]

#sgn03 = lapply(names(psi.tsm),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.3,border.stages,s))
sgn05m = getAgeASchanges(psi.tsm,meta.tsm,0.3,border.stages,'mouse')
up.seg.cnt.m = apply(sgn05m,2,function(x)getNoOfEvents(anns$mouse[x=='u',]))
dw.seg.cnt.m = apply(sgn05m,2,function(x)getNoOfEvents(anns$mouse[x=='d',]))

up.gen.cnt.m = rbind(cod=apply(sgn05m,2,function(x)length(unique(unlist(seg2ens$mouse[rownames(anns$mouse)[anns$mouse$cod == 'c' & x=='u']])))),
										 all=apply(sgn05m,2,function(x)length(unique(unlist(seg2ens$mouse[rownames(anns$mouse)[x=='u']])))))
up.gen.cnt.m[2,] = up.gen.cnt.m[2,] - up.gen.cnt.m[1,]
dw.gen.cnt.m = rbind(cod=apply(sgn05m,2,function(x)length(unique(unlist(seg2ens$mouse[rownames(anns$mouse)[anns$mouse$cod == 'c' & x=='d']])))),
										 all=apply(sgn05m,2,function(x)length(unique(unlist(seg2ens$mouse[rownames(anns$mouse)[x=='d']])))))
dw.gen.cnt.m[2,] = dw.gen.cnt.m[2,] - dw.gen.cnt.m[1,]

up.gen.cnt.m = apply(sgn05m,1,function(x)getNoOfEvents(anns$mouse[x=='u',]))
col=rep(c('red','orange','magenta','blue'),each=3)
den=rep(c(-1,40,15),times=4)

pdf('figures/paper.figures/1.pdf',w=12,h=8)
par(mfcol=c(2,3),tck=-0.01,mgp=c(1.1,0.2,0),mar=c(6,3,1.5,0),oma=c(0,0,2,1))
barplot(all.seg.cnt,col=col,den=den,ylab='# events',las=3,main='All detected events')
legend('topright',fill=c('red','orange','magenta','blue','black','black','black'),den=c(-1,-1,-1,-1,-1,40,15),legend=c('cassette','alt. acc.','alt. don.','int. ret.','cod.','partially cod.','non cod.'))
plotPanelLetter('A')
barplot(all.ens.cnt,col='#333333',den=c(-1,20),las=3,ylab='# genes',main='Protein-coding Ensembl genes')
legend('topright',fill='black',den=c(-1,20),legend=c('coding AS','non-coding AS'))
plotPanelLetter('B')

barplot(filtered.seg.cnt,col=col,den=den,ylab='# events',las=3,main='Events passed filtering')
#legend('topright',fill=c('red','orange','magenta','blue','black','black','black'),den=c(-1,-1,-1,-1,-1,40,15),legend=c('cassette','alt. acc.','alt. don.','int. ret.','cod.','partially cod.','non cod.'))
plotPanelLetter('C')
barplot(filtered.ens.cnt,col='#333333',den=c(-1,20),las=3,ylab='# genes',main='Protein-coding Ensembl genes')
#legend('topright',fill='black',den=c(-1,20),legend=c('coding AS','non-coding AS'))
plotPanelLetter('D')


barplot(up.seg.cnt.m,col=col,den=den,ylab='# events',las=3,main='Events with dPSI > 0.5 (mouse)',ylim=c(-max(apply(dw.seg.cnt.m,2,sum)),max(apply(up.seg.cnt.m,2,sum))))
barplot(-dw.seg.cnt.m,col=col,den=den,add=T,xaxt='n')
legend('topright',legend = 'up',bty='n')
legend('bottomleft',legend = 'down',bty='n')
plotPanelLetter('E')

barplot(up.gen.cnt.m,col='#333333',den=c(-1,20),ylab='# Ensamble protein-coding genes',las=3,main='Events with dPSI > 0.5 (mouse)',ylim=c(-max(apply(dw.gen.cnt.m,2,sum)),max(apply(up.gen.cnt.m,2,sum))))
barplot(-dw.gen.cnt.m,col='#333333',den=c(-1,20),add=T,xaxt='n')
plotPanelLetter('F')
legend('topright',legend = 'up',bty='n')
legend('bottomleft',legend = 'down',bty='n')
dev.off()

# figure 1.1 - overlap of age-AS abd inclusion vs exclusion exons ####
DPSI = 0.5
age.ad      = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,DPSI,border.stages,s)[anns[[s]]$sites=='ad',])
orth.age.ad = lapply(rownames(species),function(s)getAgeASchanges(orth.seg.ad.tsm,meta.tsm,DPSI,border.stages,s))
names(orth.age.ad) = names(age.ad) = rownames(species)
##### exons not present in Ensembl#
hdir = apply(age.ad$human,1,function(x){paste(sort(unique(x[x %in% c('u','d')])),collapse='')})
table(hdir)
table(hdir,anns$human$ens.exon.overlap[anns$human$sites=='ad'])
table(age.ad$human[,'brain'],anns$human$ens.exon.overlap[anns$human$sites=='ad'])
table(age.ad$human[,'testis'],anns$human$ens.exon.overlap[anns$human$sites=='ad'])

u=apply(age.ad$human,2,function(x)table(anns$human[rownames(age.ad$human),'ens.exon.overlap'][x=='u']))
d=apply(age.ad$human,2,function(x)table(anns$human[rownames(age.ad$human),'ens.exon.overlap'][x=='d']))
barplot(sweep(u,2,apply(u,2,sum),'/'),legend.text = T)
barplot(sweep(d,2,apply(d,2,sum),'/'),legend.text = T)


barplot(rbind(sweep(u,2,apply(u,2,sum),'/')['e',],sweep(d,2,apply(d,2,sum),'/')['e',]),col=c('red','blue'),beside = T)
t=anns$human[rownames(age.ad$human)[age.ad$human[,1]=='u'],]
t[t$ens.exon.overlap=='-'& t$length<27,]
plot(phastcons[['hum.7824.s14']])
plotTissueAgeProile(psi.tsm$human['hum.7824.s14',],meta.tsm)
seg2ens$human[['hum.7824.s14']]
plotTissueAgeProile(ens.ge.cod$human$rpkm['ENSG00000197746',],meta)


s1 = 100
s2 = 240
a = getPhastconsProf(phastcons[rownames(t[t$ens.exon.overlap=='-' & t$cds.pos!='cds' & t$length<500,])],s1,s2)
x = 1:((s2-s1+1)*2)
plotArea(x,a,col=cols[i],new = T,ylim=c(0.1,1),xaxt='n',ylab='mean phastcons',sd.mult=2,xlim=c(0,350),xlab='')
#####


apply(age.ad$mouse,2,table)
age.ad.over = lapply(age.ad,function(x){
	x[x=='-'] = NA
	x = cbind(x=='u',x=='d')
	colnames(x) = paste(colnames(x),rep(c('up','dw'),each=ncol(x)/2))
	caclSegOverlap(x)
})
cols = unique(meta[,c('tissue','col')])
cols = setNames(cols$col,cols$tissue)

# exon length
s = 'mouse'
up.ad.len = apply(age.ad[[s]],2,function(t){
	f = anns[[s]]$cod=='c' & rownames(anns[[s]]) %in% rownames(age.ad[[s]])[t=='u']
	s = sum(anns[[s]][f,'length'] %% 3 == 0);
	t = sum(f)
	c(len3=s,total=t,freq=s/t,conf=binom.test(s,t)$conf.int)
})

dw.ad.len = apply(age.ad[[s]],2,function(t){
	f = anns[[s]]$cod=='c' & rownames(anns[[s]]) %in% rownames(age.ad[[s]])[t=='d']
	s = sum(anns[[s]][f,'length'] %% 3 == 0);
	t = sum(f)
	c(len3=s,total=t,freq=s/t,conf=binom.test(s,t)$conf.int)
})

up2dw.len.pv = sapply(1:7,function(t){prop.test(c(up.ad.len[1,t],dw.ad.len[1,t]),c(up.ad.len[2,t],dw.ad.len[2,t]))$p.value})
cnst.len3.freq=mean(all.anns[[s]][all.anns[[s]]$sites=='ad' & all.anns[[s]]$type=='EXN' & all.anns[[s]]$cod=='c' & all.anns[[s]]$gene_id %in% all.anns[[s]][rownames(age.ad[[s]])[apply(age.ad[[s]] =='u' | age.ad[[s]] =='d',1,sum)>0],'gene_id'],'length'] %% 3 == 0)

# GO 
library(GO.db)
library(goseq)
mouse.go = read.table('input/GO/Mus_musculus.GRCm38.84.GO.csv.gz',sep=',',header = T)
mouse.go = mouse.go[mouse.go[,2] != '',]
mouse.go = split(mouse.go$GO.Term.Accession,mouse.go$Ensembl.Gene.ID)
GOALLANCESTOR = c(as.list(GOCCANCESTOR),as.list(GOMFANCESTOR),as.list(GOBPANCESTOR))
mouse.go.full = lapply(mouse.go,function(x){r=unique(c(x,unlist(GOALLANCESTOR[x])));r[grep('GO:',r,fixed=T)]})
mouse.go.full.rev = revList(mouse.go.full)
mouse.go.rev = revList(mouse.go)
min.go.size = 3
age.go.up = lapply(colnames(age.ad[[s]]),function(t){
	tested = unique(unlist(seg2ens[[s]][rownames(age.ad[[s]])[age.ad[[s]][,t]!='']]))
	sign   = unique(unlist(seg2ens[[s]][rownames(age.ad[[s]])[age.ad[[s]][,t]=='u']]))
	getGO(sign,tested,min.go.size=min.go.size,'Hypergeometric',gene2cat=mouse.go.full)
	})

age.go.dw = lapply(colnames(age.ad[[s]]),function(t){
	tested = unique(unlist(seg2ens[[s]][rownames(age.ad[[s]])[age.ad[[s]][,t]!='']]))
	sign   = unique(unlist(seg2ens[[s]][rownames(age.ad[[s]])[age.ad[[s]][,t]=='d']]))
	getGO(sign,tested,min.go.size=min.go.size,'Hypergeometric',gene2cat=mouse.go.full)
})

names(age.go.up) = names(age.go.dw) = colnames(age.ad[[s]])
sapply(age.go.up,function(x)sum(x$qv.over<0.1))
sapply(age.go.dw,function(x)sum(x$qv.over<0.1))

go.up.cnt=sapply(age.go.up,function(x)table(factor(x$ontology[x$qv.over<0.1],levels = c('BP','CC','MF'))))
go.dw.cnt=sapply(age.go.dw,function(x)table(factor(x$ontology[x$qv.over<0.1],levels = c('BP','CC','MF'))))

go.up.descr=sapply(age.go.up,function(x){sapply(c('BP','CC','MF'),function(o){x$term[x$ontology==o][1]})})
go.dw.descr=sapply(age.go.dw,function(x){sapply(c('BP','CC','MF'),function(o){x$term[x$ontology==o][1]})})

go.up.descr = apply(go.up.descr,1:2,wrapText,maxlen=27,max.chunks=2)
go.dw.descr = apply(go.dw.descr,1:2,wrapText,maxlen=27,max.chunks=2)
go.up.descr[go.up.cnt==0] = ''
go.dw.descr[go.dw.cnt==0] = ''
# go.up.cnt = log(go.up.cnt+1)
# go.dw.cnt = log(go.dw.cnt+1)

## consrvation of age-patterns
sps = c('mouse','rat','rabbit','human','opossum','chicken')
#sps = c('human','mouse','opossum','chicken')

cor.age = lapply(unique(meta$tissue),function(t){
	agesb = setNames(sapply(border.stages[sps],function(x)x[t,2]),sps)
	caclOrthPsiCor(orth.seg.ad.tsm,t,agesb,f=T,cor.method = 'pe')
}) 
names(cor.age) = unique(meta$tissue)



orth.age.ad5 = lapply(rownames(species),function(s)getAgeASchanges(orth.seg.ad.tsm,meta.tsm,psi.thr = 0.5,border.stages,s))
names(orth.age.ad5) = rownames(species)
u5=getASChangeCons(orth.age.ad5,'u',sps)
d5=getASChangeCons(orth.age.ad5,'d',sps)

orth.age.ad3 = lapply(rownames(species),function(s)getAgeASchanges(orth.seg.ad.tsm,meta.tsm,psi.thr = 0.3,border.stages,s))
names(orth.age.ad3) = rownames(species)
u3=getASChangeCons(orth.age.ad3,'u',sps)
d3=getASChangeCons(orth.age.ad3,'d',sps)

pdf('figures/paper.figures/1.2.conservation.pdf',w=7,h=11)
par(mfrow=c(3,2),mar=c(5,2.1,3,1),tck=-0.01,mgp=c(1.1,0.2,0))
plotAgeChangeCons(u5,ylim=c(0,1),ylab='fraction of shared',main='Inclusion, dPSI=0.5',xlab='')
plotAgeChangeCons(d5,ylim=c(0,1),ylab='fraction of shared',main='Exclusion, dPSI=0.5',xlab='')

plotAgeChangeCons(u3,ylim=c(0,1),ylab='fraction of shared',main='Inclusion, dPSI=0.3',xlab='')
plotAgeChangeCons(d3,ylim=c(0,1),ylab='fraction of shared',main='Exclusion, dPSI=0.3',xlab='')
plotAgeChangeCons(cor.age,ylim=c(0.6,1),ylab='Pearson correlation',main='Conservation of adult AS',xlab='',inxs=1:3)
dev.off()

pdf('figures/paper.figures/1.1.all.species.dpsi=0.5.pdf',w=15,h=15)
par(mfrow=c(3,3),mar=c(1.5,1,3,6))
for(s in names(age.ad.over))
	plotAgeSegOverlap(age.ad.over[[s]],main=paste('Overlap of ageAS exons across',s,'tissues'))
dev.off()


pdf('figures/paper.figures/1.1.mouse.dpsi=0.3.pdf',w=10,h=9)
s = 'mouse'
par(mar=c(1.5,1,3,6),oma=c(0,0,0,0))
layout(matrix(c(1,1,4,2,3,4),nrow=3),heights = c(1,1,1.5))
plotAgeSegOverlap(age.ad.over[[s]],main=paste('Overlap of ageAS exons across',s,'tissues'))
plotPanelLetter('A')

age.ad.tis.cnt. = rbind(up=table(factor(apply(age.ad[[s]][,-2]=='u',1,sum),levels = 0:6))[-1],
											 dw=table(factor(apply(age.ad[[s]][,-2]=='d',1,sum),levels = 0:6))[-1])
par(mar=c(2.5,2.5,3,1),mgp=c(1.6,0.6,0))
age.ad.tis.cnt = log(age.ad.tis.cnt.+1)
b=barplot(age.ad.tis.cnt[1,],ylim=c(-max(age.ad.tis.cnt[2,]),max(age.ad.tis.cnt[1,]))*1.1,xlab='# of tissues',ylab='# of exons',main='Number of shared age-relaed exons',yaxt='n')
text(b,age.ad.tis.cnt[1,],age.ad.tis.cnt.[1,],adj = c(0.5,-0.05))
b=barplot(-age.ad.tis.cnt[2,],add=T,yaxt='n')
text(b,-age.ad.tis.cnt[2,],age.ad.tis.cnt.[2,],adj = c(0.5, 1.05))
lab = c(0,50,1000)
at = c(-log(rev(lab)+1),log(lab[-1]+1))
lab = c(rev(lab),lab[-1])
axis(2,at,lab,las=2)
plotPanelLetter('B')

x=1:7
par(mar=c(5,2.5,3,1))
ylim=range(up.ad.len[3:5,],dw.ad.len[3:5,],cnst.len3.freq)
ylim[2] = ylim[2]+(ylim[2]-ylim[1])*.1
plot(x,up.ad.len[3,],pch=19,col=cols[colnames(up.ad.len)],xlab='',ylab='N*3 exon freq.',ylim=ylim,xaxt='n',xlim=c(1,7.3),main='Frequency of exons with length dividible by 3')
segments(x,up.ad.len[4,],x,up.ad.len[5,],col=cols[colnames(up.ad.len)])
points(x+0.3,dw.ad.len[3,],pch=1,col=cols[colnames(up.ad.len)])
segments(x+0.3,dw.ad.len[4,],x+0.3,dw.ad.len[5,],col=cols[colnames(up.ad.len)],lty=2)
abline(h=cnst.len3.freq)
axis(1,x+0.15,colnames(up.ad.len),las=1)
par(mar=c(3,3,0.5,1))
plotPanelLetter('C')
for(i in 1:7)
	if(up2dw.len.pv[i]<0.05){
		segments(x[i],ylim[2],x[i],up.ad.len[5,i]+0.04)
		segments(x[i],ylim[2],x[i]+0.3,ylim[2])
		segments(x[i]+0.3,ylim[2],x[i]+0.3,dw.ad.len[5,i]+0.04)
	}

ylim=max(go.dw.cnt,go.up.cnt)
par(mar=c(3,3,0.5,1))
b=barplot( go.up.cnt,beside = T,yaxt='n',ylim=c(-1000,1000),ylab='# of GO terms',legend.text = T,xlim=c(2,29))
text(b,go.up.cnt+10,go.up.descr,srt=90,adj=c(0,0.5),cex=1)
b=barplot(-go.dw.cnt,beside = T,yaxt='n',add=TRUE)
text(b,-go.dw.cnt-10,go.dw.descr,srt=90,adj=c(1,0.5),cex=1)
#lab = c(0,5,20,50,200)
lab = c(0,200)
#at = c(-log(rev(lab)+1),log(lab[-1]+1))
at = c(-rev(lab),lab[-1])
lab = c(rev(lab),lab[-1])
axis(2,at,lab,las=2)
plotPanelLetter('D')
dev.off()

#+domain enrichment
################################################
## figure 2: AS in up genes, AS is old, and depleted in SF and TF 
DPSI = 0.3
# (what is about "always alts"?). Depletion in mutations, conservations (Phastcons + SNP)
# ? transcriptom complexity on age

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
library(GO.db)
mouse.go = read.table('input/GO/Mus_musculus.GRCm38.84.GO.csv.gz',sep=',',header = T)
mouse.go = mouse.go[mouse.go[,2] != '',]
mouse.go = split(mouse.go$GO.Term.Accession,mouse.go$Ensembl.Gene.ID)
GOALLANCESTOR = c(as.list(GOCCANCESTOR),as.list(GOMFANCESTOR),as.list(GOBPANCESTOR))
mouse.go.full = lapply(mouse.go,function(x){r=unique(c(x,unlist(GOALLANCESTOR[x])));r[grep('GO:',r,fixed=T)]})
mouse.go.full.rev = revList(mouse.go.full)

as.in.ge.patterns.mouse = readRDS('Rdata/paper.figures/as.in.ge.patterns.mouse.Rdata')
#as.in.ge.patterns.mouse.cnt = lapply(as.in.ge.patterns.mouse,function(x){table(factor(x[,3]),pmax(x[,1],-x[,2])>DPSI)[,c('TRUE','FALSE')]})
#as.in.ge.patterns.mouse.stat = sapply(as.in.ge.patterns.mouse.cnt,function(x){cbind(my.binom.test(x['Increasing',]),my.binom.test(x['Decreasing',]))})


as.in.ge.patterns.mouse.cnt.up = lapply(as.in.ge.patterns.mouse,function(x){table(factor(x[,3]),x[,1]>DPSI)[,c('TRUE','FALSE')]})
as.in.ge.patterns.mouse.stat.up = sapply(as.in.ge.patterns.mouse.cnt.up,function(x){cbind(my.binom.test(x['Increasing',]),my.binom.test(x['Decreasing',]))})

as.in.ge.patterns.mouse.cnt.dw = lapply(as.in.ge.patterns.mouse,function(x){table(factor(x[,3]),x[,2]< -DPSI)[,c('TRUE','FALSE')]})
as.in.ge.patterns.mouse.stat.dw = sapply(as.in.ge.patterns.mouse.cnt.dw,function(x){cbind(my.binom.test(x['Increasing',]),my.binom.test(x['Decreasing',]))})


ge.info.m = read.csv('/home/mazin/skoltech/projects/evo.devo/input/gene.info.from.marg/Mouse.Indexes.csv')
rownames(ge.info.m) = ge.info.m$Mouse_ID
up.ens.ids = lapply(as.in.ge.patterns.mouse,function(x)rownames(x)[x[,1] >  DPSI])
dw.ens.ids = lapply(as.in.ge.patterns.mouse,function(x)rownames(x)[x[,2] < -DPSI])
ens.tested = unique(unlist(seg2ens$mouse[rownames(anns$mouse)[anns$mouse$sites=='ad' & anns$mouse$cod !='n']]))
f = rownames(ge.info.m) %in% rownames(ens.ge.cod$mouse$gene)
gene.age = my.binom.test(table(ge.info.m$Age[!(rownames(ge.info.m) %in% ens.tested) & f]==0)[c('TRUE','FALSE')])
gene.age = rbind(no.alt = gene.age,no.ch = my.binom.test(table(ge.info.m$Age[f& (rownames(ge.info.m) %in% ens.tested) & !(rownames(ge.info.m) %in% unlist(c(up.ens.ids,dw.ens.ids)))]==0)[c('TRUE','FALSE')]))

#it is wrong because uses only genes with defined patterns
for(t in names(as.in.ge.patterns.mouse)){
	z=as.in.ge.patterns.mouse[[t]]
	gene.age = rbind(gene.age,
									 my.binom.test(table(ge.info.m$Age[f & (rownames(ge.info.m) %in% rownames(z)[z[,1] >  DPSI])]==0)[c('TRUE','FALSE')]),
									 my.binom.test(table(ge.info.m$Age[f & (rownames(ge.info.m) %in% rownames(z)[z[,2] < -DPSI])]==0)[c('TRUE','FALSE')]))
	#gene.age = rbind(gene.age,c(sum(f & (rownames(ge.info.m) %in% rownames(z)[z[,1] >  DPSI])),sum(f & (rownames(ge.info.m) %in% rownames(z)[z[,2] <  -DPSI]))))
	rownames(gene.age)[-(1:(nrow(gene.age)-2))] = paste(t,c('up','dw'))
}


cols = unique(meta[,c('tissue','col')])
cols = setNames(cols$col,cols$tissue)

# transcriptome complexity
m = getAgeASchanges(psi.tsm,meta.tsm,DPSI,border.stages,'mouse')
mm = apply(m,2,function(x){x %in% c('u','d')})
tr.compl.m =getAltExonStat(psi.tsm$mouse[anns$mouse$sites=='ad' & apply(mm,1,sum)>0,],meta.tsm,0.1,tissues = c('brain','heart','liver','ovary','testis'),na.as.cnst=TRUE)
# TFs and SFs
tf.sf.up = do.call(rbind,lapply(as.in.ge.patterns.mouse,function(x)data.frame(rownames(x),x[,1])))
tf.sf.dw = do.call(rbind,lapply(as.in.ge.patterns.mouse,function(x)data.frame(rownames(x),x[,2])))

tf.sf. = data.frame(dpsi.up=sapply(split(tf.sf.up[,2],tf.sf.up[,1]),max),dpsi.dw=sapply(split(tf.sf.dw[,2],tf.sf.dw[,1]),min))
tf.sf = cbind(dpsi=pmax(tf.sf.$dpsi.up,-tf.sf.$dpsi.dw),tf = rownames(tf.sf.) %in% rownames(ge.info.m)[ge.info.m$TF=='Mouse_TF'],tf.sf.)

#GO:0008380 = RNA splicing, GO:0003723 - RNA binding
tf.sf$sf = rownames(tf.sf) %in% mouse.go.full.rev[['GO:0008380']]
tf.sf$rna = rownames(tf.sf) %in% mouse.go.full.rev[['GO:0003723']]
tf.sf = cbind(tf.sf,tf.sf.)
rm(tf.sf.)
# look on PanAS
alt.prop = t(apply(psi.tsm$mouse[anns$mouse$sites=='ad' & anns$mouse$cod!='n',],1,function(x){t = x>0.1 & x<0.9;c(in.exp=mean(t,na.rm=TRUE),in.all=mean(!is.na(t) & t))}))
e2s = revList(seg2ens$mouse[rownames(alt.prop)])[rownames(tf.sf)]
tf.sf$alt.prop.in.exp = sapply(e2s,function(sids){max(alt.prop[sids,1])})
tf.sf$alt.prop.in.all = sapply(e2s,function(sids){max(alt.prop[sids,2])})
f = function(x){rbind(x,pv=sapply(1:ncol(x),function(i)prop.test(x[1,c(1,i)],x[1,c(1,i)]+x[2,c(1,i)])$p.value))}

# rna.on.dpsi = f(sapply(split(factor(tf.sf$rna | tf.sf$sf),floor(tf.sf$dpsi*5)/5),function(x){z=table(x)[c('TRUE','FALSE')];c(z,my.binom.test(z))})[,1:5])
# tf.on.dpsi = f(sapply(split(factor(tf.sf$tf),floor(tf.sf$dpsi*5)/5),function(x){z=table(x)[c('TRUE','FALSE')];c(z,my.binom.test(z))})[,1:5])
# rna.on.alt = f(sapply(split(factor(tf.sf$rna | tf.sf$sf),floor(tf.sf$alt.prop.in.exp*5)/5),function(x){z=table(x)[c('TRUE','FALSE')];c(z,my.binom.test(z))}))
# tf.on.alt = f(sapply(split(factor(tf.sf$tf),floor(tf.sf$alt.prop.in.exp*5)/5),function(x){z=table(x)[c('TRUE','FALSE')];c(z,my.binom.test(z))}))
# 
# f=anns$mouse$sites=='ad' & anns$mouse$cod=='c'
# plot(apply(m,2,function(t){mean(anns$mouse$length[f & t=='u']%%3==0)}),col='red',pch=19,ylim=c(0.35,0.87))
# points(apply(m,2,function(t){mean(anns$mouse$length[f & t=='d']%%3==0)}),col='blue',pch=19)
# points(apply(m,2,function(t){mean(anns$mouse$length[f & t=='n']%%3==0)}),col='black',pch=19)
# abline(h=mean(all.anns$mouse$length[all.anns$mouse$type=='EXN' & all.anns$mouse$cod=='c' & all.anns$mouse$sites=='ad'] %% 3==0),lty=2)

# phastacons
orth.alt = sapply(orth.seg.ad,function(x)x$seg$type)=='ALT'
orth.alt = cbind(orth.alt[,1:2],mrb=apply(orth.alt[,3:5],1,sum)>0,orth.alt[,6:7])
rownames(orth.alt) = rownames(orth.seg.ad$human$seg)
phastcons = read.table('/home/mazin/skoltech/projects/evo.devo/processed/ad.phastcons.gz',sep='\t')
phastcons = setNames(strsplit(phastcons[,2],',',TRUE),phastcons[,1])
phastcons = lapply(phastcons,as.numeric)
age.sids = lapply(rownames(species),function(s){
	r = getAgeASchanges(psi.tsm,meta.tsm,DPSI,border.stages,s)
	rownames(r)[apply(r=='u' | r=='d',1,sum)>0 & anns[[s]]$sites=='ad' & anns[[s]]$cod=='c']
})
names(age.sids) = rownames(species)



#human.alt = rownames(anns$human)[anns$human$sites=='ad' & anns$human$cod=='c']
orth.seg.ad.all.id=readRDS('Rdata/orth.seg.ad.all.id.Rdata')
human.age = getAgeASchanges(psi.tsm,meta.tsm,DPSI,border.stages,'human')
sids = list()
sids$alt =     rownames(human.age)[apply(human.age =='u',1,sum)==0 & apply(human.age =='d',1,sum)==0 & anns$human$sites=='ad' & anns$human$cod=='c']
sids$`human up` =      rownames(human.age)[apply(human.age =='u',1,sum) >0 & apply(human.age =='d',1,sum)==0 & anns$human$sites=='ad' & anns$human$cod=='c']
sids$`human down` =    rownames(human.age)[apply(human.age =='u',1,sum)==0 & apply(human.age =='d',1,sum) >0 & anns$human$sites=='ad' & anns$human$cod=='c']
sids$`ancient up` = intersect(sids$`human up`,rownames(orth.alt)[apply(orth.alt,1,sum)==5])
sids$`ancient down` = intersect(sids$`human down`,rownames(orth.alt)[apply(orth.alt,1,sum)==5])
sids = c(const = list(rownames(all.anns$human)[all.anns$human$sites=='ad' & all.anns$human$cod=='c' & all.anns$human$type=='EXN' & all.anns$human$gene_id %in% all.anns$human[unlist(sids),'gene_id']]),sids)
#sids$const     = setdiff(orth.seg.ad.all.id[,1],rownames(orth.seg.ad$human$seg))
#sids$h.alt     = rownames(orth.seg.ad$human$seg)[ orth.alt[,1] & !orth.alt[,2] & !orth.alt[,3] & !orth.alt[,4] & !orth.alt[,5]]
#only human age
#sids$h.age     = intersect(sids$h.alt,age.sids$human)
#sids$hq.alt    = rownames(orth.seg.ad$human$seg)[ orth.alt[,1] &  orth.alt[,2] & !orth.alt[,3] & !orth.alt[,4] & !orth.alt[,5]]
#sids$hqm.alt   = rownames(orth.seg.ad$human$seg)[ orth.alt[,1] &  orth.alt[,2] &  orth.alt[,3] & !orth.alt[,4] & !orth.alt[,5]]
#sids$hqmo.alt  = rownames(orth.seg.ad$human$seg)[ orth.alt[,1] &  orth.alt[,2] &  orth.alt[,3] &  orth.alt[,4] & !orth.alt[,5]]
#sids$hqmoc.alt = rownames(orth.seg.ad$human$seg)[ orth.alt[,1] &  orth.alt[,2] &  orth.alt[,3] &  orth.alt[,4] &  orth.alt[,5]]
# human+mouse+opossum age, any alt
#sids$hmo.age  = orth.seg.ad.all.id[orth.seg.ad.all.id[,'human'] %in% age.sids$human & orth.seg.ad.all.id[,'mouse'] %in% age.sids$mouse & orth.seg.ad.all.id[,'opossum'] %in% age.sids$opossum,1]
# sids$const = rownames(all.anns$human)[all.anns$human$sites=='ad' & all.anns$human$cod=='c' & all.anns$human$type=='EXN' & all.anns$human$gene_id %in% anns$human[human.alt,'gene_id']]
# sids$alt = setdiff(human.alt,human.age)
# sids$h.age     = setdiff(human.age,rownames(orth.seg.ad$human$seg)[apply(orth.alt,1,sum)>1])
# sids$hq.age    = intersect(human.age,rownames(orth.seg.ad$human$seg)[orth.alt[,1] & !orth.alt[,2] & !orth.alt[,3] & !orth.alt[,4]])
# sids$hqm.age   = intersect(human.age,rownames(orth.seg.ad$human$seg)[orth.alt[,1] &  orth.alt[,2] & !orth.alt[,3] & !orth.alt[,4]])
# sids$hqmo.age  = intersect(human.age,rownames(orth.seg.ad$human$seg)[orth.alt[,1] &  orth.alt[,2] &  orth.alt[,3] & !orth.alt[,4]])
# sids$hqmoc.age = intersect(human.age,rownames(orth.seg.ad$human$seg)[orth.alt[,1] &  orth.alt[,2] &  orth.alt[,3] &  orth.alt[,4]])
# sids$hq.alt    = rownames(orth.seg.ad$human$seg)[orth.alt[,1] &  orth.alt[,2] & !orth.alt[,3] & !orth.alt[,4] & !orth.alt[,5]]
# sids$hqm.alt   = rownames(orth.seg.ad$human$seg)[orth.alt[,1] &  orth.alt[,2] &  orth.alt[,3] & !orth.alt[,4] & !orth.alt[,5]]
# sids$hqmo.alt  = rownames(orth.seg.ad$human$seg)[orth.alt[,1] &  orth.alt[,2] &  orth.alt[,3] &  orth.alt[,4] & !orth.alt[,5]]
# sids$hqmoc.alt = rownames(orth.seg.ad$human$seg)[orth.alt[,1] &  orth.alt[,2] &  orth.alt[,3] &  orth.alt[,4] &  orth.alt[,5]]

sapply(sids,length)

s1=100
s2=240
phast.areas = lapply(sids,function(s)getPhastconsProf(phastcons[s],s1,s2))
rm(phastcons);gc()

## gnomad
gnomad = read.table('processed/gnomad201/human.ad.snp.tab.gz',sep='\t')
colnames(gnomad) = c('chr_id','pos','id','ref','alt','qual','filter','seg_id','alt_cnt','freq','tot_cnt')

hann = all.anns$human[unique(unlist(sids)),]
gnomad = gnomad[gnomad$seg_id %in% rownames(hann),]
gc()

s2g = hann[gnomad$seg_id,]
gnomad$dist2start = gnomad$pos - s2g$start
gnomad$dist2stop  = s2g$stop - gnomad$pos
gnomad[s2g$strand==-1,c('dist2start','dist2stop')] = gnomad[s2g$strand==-1,c('dist2stop','dist2start')]
gnomad$strand = s2g$strand


snp.distr = list()
snp.distr$ppt=getSNPDistr(g=gnomad[gnomad$dist2start<0 & gnomad$dist2start>= -25,],sids)
snp.distr$exon=getSNPDistr(g=gnomad[gnomad$dist2start>=0 & gnomad$dist2stop>=0,],sids)
snp.distr$down=getSNPDistr(g=gnomad[gnomad$dist2stop>=-25 & gnomad$dist2stop<0,],sids)
#plot(snp.distr$exon[1,])
## hgmd
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

h2a = read.table('input/hgmd/tosend/HighLevelConcepts_ALL_DM_2017_1_mut2AgeOfOnset.csv',sep=',',header=T,quote='',comment.char = '')
colnames(h2a) = c('c','p','a','l') #Congenital disorders,Post-natal disorder,Adult onset,All ages
h2a = unlist(lapply(1:ncol(h2a),function(i){x=h2a[h2a[,i]!='',i];setNames(rep(colnames(h2a[i]),length(x)),x)}))
h2a = split(h2a,names(h2a))
sort(table(sapply(h2a,function(x){paste(sort(x),collapse='')})))

h2c = c2h= read.table('input/hgmd/tosend/HighLevelConcepts_ALL_DM_2017_1_mut2concept.csv',sep=',',header=T,quote='',comment.char = '')
clevels = colnames(h2c)
h2c = unlist(lapply(1:ncol(h2c),function(i){x=h2c[h2c[,i]!='',i];setNames(rep(colnames(h2c[i]),length(x)),x)}))
h2c = split(h2c,names(h2c))
table(sapply(h2c,length))
h2c[1:10]



getHGMDexonStat = function(sids,h,locs){
	s = sum(sids %in% h[grep(locs,h$loc),'seg.id'])
	f = length(sids) - s
	c(has.mut = s,no.mut=f,my.binom.test(s,f))
}

hgmd.stat  = lapply(c(PPT='A|a|p',Exon='e',`Donor site`='D|d'),function(loc)sapply(sids,function(x){getHGMDexonStat(x,hgmd2hadf,loc)}))
hgmd.stat = lapply(hgmd.stat,function(t){
	rbind(t,pv=apply(t,2,function(c){fisher.test(cbind(t[1:2,1],c[1:2]))$p.value}))
})




pdf(paste0('figures/paper.figures/2dPSI=',DPSI,'.pdf'),w=10,h=10/3*4)
par(mfrow=c(4,3),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(4,2,1.5,0),oma=c(0,0,1,1))
cols = unique(meta[,c('tissue','col')])
cols = setNames(cols$col,cols$tissue)
c = rep(cols[colnames(as.in.ge.patterns.mouse.stat.up)],each=2)
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
u = c(1:2,1:7*2+1)
b=barplot(gene.age[u,1],col=c('#555555','gray',cols),yaxt='n',ylim=range(-gene.age,gene.age),xaxt='n',ylab='proportion of ancent genes',main='Evolutionary age of genes')
segments(b,gene.age[u,2],b,gene.age[u,3])
u = c(1:2,1:7*2+2)
b=barplot(-gene.age[u,1],col=c('#555555','gray',cols),add=T,xaxt='n',yaxt='n')
segments(b,-gene.age[u,2],b,-gene.age[u,3])
axis(1,b,c('no AS','no change',names(as.in.ge.patterns.mouse)),las=3)
at = -2:2*0.4
axis(2,at,abs(at))
abline(h=0)
abline(h=c(gene.age[1:2,1],-gene.age[1:2,1]),lty=3)
text(b[1],par('usr')[4],'AS up',adj=c(0,1.1))
text(b[1],par('usr')[3],'AS down',adj=c(0,-1.1))

plotPanelLetter('B')
mar=par(mar=c(4,2,1.5,5))
col = unique(meta[,c('tissue','col')])
col = setNames(col$col,substr(col$tissue,1,1))
col = c(col,const='black',bh='orange','2ts'='#666666','3ts'='#777777','4ts'='#AAAAAA','5ts'='#DDDDDD')

areaplot(tr.compl.m,col=col[rownames(tr.compl.m)],ylab='# of exons with PSI in [0.1,0.9]',main='Mouse transcriptom complexity',xaxt='n',xlab='Age (days from conception)')
axis(1,1:ncol(tr.compl.m),colnames(tr.compl.m))
par(mar=mar)
plotPanelLetter('C')

# x = as.numeric(colnames(rna.on.dpsi))
# ylim=range(rna.on.dpsi[3:5,],tf.on.dpsi[3:5,])
# plot(x,rna.on.dpsi[3,],t='b',pch=ifelse(rna.on.dpsi[6,]<0.01,19,1),col='red',ylim=ylim,xlab='dPSI',ylab='fraqtion of TFs or SFs')
# segments(x,rna.on.dpsi[4,],x,rna.on.dpsi[5,],col='red')
# lines(x+0.02,tf.on.dpsi[3,],t='b',pch=ifelse(tf.on.dpsi[6,]<0.01,19,1),col='blue')
# segments(x+0.02,tf.on.dpsi[4,],x+0.02,tf.on.dpsi[5,],col='blue')
# legend('bottomleft',bty='n',col=c('red','blue'),lwd=1,legend=c('RNA binding','TFs'))
# plotPanelLetter('D')
# 
# x = as.numeric(colnames(rna.on.alt))
# ylim=range(rna.on.alt[3:5,],tf.on.alt[3:5,])
# plot(x,rna.on.alt[3,],t='b',pch=ifelse(rna.on.alt[6,]<0.01,19,1),col='red',ylim=ylim,xlab='prop. of samples PSI in [0.1,0.9]',ylab='fraqtion of TFs or SFs')
# segments(x,rna.on.alt[4,],x,rna.on.alt[5,],col='red')
# lines(x+0.02,tf.on.alt[3,],t='b',pch=ifelse(tf.on.alt[6,]<0.01,19,1),col='blue')
# segments(x+0.02,tf.on.alt[4,],x+0.02,tf.on.alt[5,],col='blue')
# plotPanelLetter('E')


plotOddsRatioBars(tf.sf$dpsi > DPSI,cbind(TFs=tf.sf$tf,SFs=tf.sf$sf,RNBPs=tf.sf$rna,panAS=tf.sf$alt.prop.in.exp>0.8),'D',main='Assotiation with age AS')
plotOddsRatioBars(tf.sf$alt.prop.in.exp>0.8,cbind(TFs=tf.sf$tf,SFs=tf.sf$sf,RNBPs=tf.sf$rna,ageAS=tf.sf$dpsi>DPSI),'E',main='Assotiation with panAS')

#cols = c('black','gray','green','violet','brown','orange','red')
cols = c('black','gray','orange','cyan','red','blue')
x = 1:((s2-s1+1)*2)
for(i in 1:length(phast.areas))
	plotArea(x,phast.areas[[i]],col=cols[i],new = i==1,ylim=c(0.1,1),xaxt='n',ylab='mean phastcons',sd.mult=2,xlim=c(0,350),xlab='')
axis(1,c(1,100,180,280),c('-100','acc.','don.',100))
abline(v=c(100,180),lty=2)
#legend('topright',bty='n',col=rev(cols),legend=rev(paste0(c('const.','human','primates','placentals','mammals','ancient','ancien age'),'(',sapply(sids,length),')')),lwd=2)
legend('topright',bty='n',col=rev(cols),legend=rev(paste0(names(sids),'(',sapply(sids,length),')')),lwd=2)
plotPanelLetter('F')
par(mar=c(6,2,1.5,0))
ylim = range(lapply(snp.distr,function(x)range(x[1:3,])))
plotSNPFreq(snp.distr$ppt,cols,ylim,main='PPT')
plotPanelLetter('G')
plotSNPFreq(snp.distr$exon,cols,ylim,main='Exon')
plotPanelLetter('H')
plotSNPFreq(snp.distr$down,cols,ylim,main='Downstream intron [1,25]')
plotPanelLetter('I')

letts = setNames(c('G','K','L'),names(hgmd.stat))
for(r in names(hgmd.stat)){
	b=barplot(hgmd.stat[[r]][3,],beside = T,legend.text = F,main=r,ylim=range(0,hgmd.stat[[r]][4:5,])*1.1,las=3,ylab='proportion of exons with disease mutations')
	segments(b,hgmd.stat[[r]][4,],b,hgmd.stat[[r]][5,])
	text(b[hgmd.stat[[r]][6,]<0.01],par('usr')[4],'*',adj=c(0.5,2),cex=2)
	plotPanelLetter(letts[[r]])
}

dev.off()

##########################
### figure 3 trees ######
library(ape)
library(TKF)
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
cisbp = read.table('input/cisbp-db/all/RBP_Information_all_motifs.txt',sep='\t',quote = '',fill=T,header = 1)
cisbp.mouse.sfs = unique(cisbp$DBID[cisbp$RBP_Species=='Mus_musculus'])

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
#sfs = rownames(mouse.ge) %in% cisbp.mouse.sfs
#sfs = rownames(mouse.ge) %in% rownames(ens.descr.mm)[grepl('splicing',ens.descr.mm$Description,ignore.case = T)]
#rna = rownames(mouse.ge) %in% mouse.go.full.rev[['GO:0003723']]

mouse.sf.stage.cors = lapply(d,function(day){
	r = mouse.ge[sfs,m$days==day]
	colnames(r) = substr(m$tissue[m$days==day],1,1)
	r = cbind(base=base[sfs],r)
	cor(r,u='p')})

tf = rownames(mouse.ge) %in% rownames(ge.info.m)[ge.info.m$TF=='Mouse_TF'] 

mouse.tf.stage.cors = lapply(d,function(day){
	r = mouse.ge[tf,m$days==day]
	colnames(r) = substr(m$tissue[m$days==day],1,1)
	r = cbind(base=base[tf],r)
	cor(r,u='p')})

names(mouse.stage.cors) = names(mouse.tf.stage.cors) = names(mouse.ge.stage.cors) = names(mouse.sf.stage.cors) = d
# saveRDS(mouse.ge.stage.cors,'Rdata/paper.figures/mouse.ge.stage.cors.Rdata')
# saveRDS(mouse.sf.stage.cors,'Rdata/paper.figures/mouse.sf.stage.cors.Rdata')
#mouse.stage.cors=readRDS('Rdata/paper.figures/mouse.stage.cors.Rdata')
#mouse.ge.stage.cors=readRDS('Rdata/paper.figures/mouse.ge.stage.cors.Rdata')
#mouse.sf.stage.cors=readRDS('Rdata/paper.figures/mouse.sf.stage.cors.Rdata')

# PSI trees
mouse.as.stage.trees = makeNJTreesByCorLastTopology(mouse.stage.cors)
mouse.ge.stage.trees = makeNJTreesByCorLastTopology(mouse.ge.stage.cors)
mouse.sf.stage.trees = makeNJTreesByCorLastTopology(mouse.sf.stage.cors)
mouse.tf.stage.trees = makeNJTreesByCorLastTopology(mouse.tf.stage.cors)

pdf('figures/paper.figures/3.pdf',w=10,h=10)
par(mfrow=c(4,6),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(4,2,1.5,0),oma=c(0,0,1,1))
plotMouseTrees(mouse.as.stage.trees,'Splicing','A')
plotMouseTrees(mouse.ge.stage.trees,'Gene expression','B')
plotMouseTrees(mouse.sf.stage.trees,'Gene expression of SF','C')
plotMouseTrees(mouse.tf.stage.trees,'Gene expression of TF','D')
dev.off()

##########################
### figure 4: orth conservations + MDS
orth.irs = readRDS('Rdata/not.cassette/orth.irs.Rdata')
orth.aas = readRDS('Rdata/not.cassette/orth.aas.Rdata')
orth.dds = readRDS('Rdata/not.cassette/orth.dds.Rdata')
orth.aas.all = readRDS('Rdata/not.cassette/orth.aas.all.Rdata')
orth.dds.all = readRDS('Rdata/not.cassette/orth.dds.all.Rdata')
orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
irs.all=readRDS('Rdata/not.cassette/irs.all.Rdata')
ir.orth.inxs = unlist(lapply(irs.all,function(x){x=x[apply(is.na(x),1,sum)==0,];paste(x$orth.inx1,x$orth.inx2)}))
table(table(ir.orth.inxs))
mds = readRDS('Rdata/orth.mds.Rdata')


barplotWithText. = function(x,t=x,...){
	b=barplot(x,ylim=range(x/2,x*1.5),...)
	text(b,x,t,adj=c(0.5,-.1))
}

pdf('figures/paper.figures/4.pdf',w=12,h=9)
par(mfrow=c(3,4),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(3,2,1.5,0),oma=c(0,0,2,1))
barplotWithText.(table(apply(sapply(orth.seg.ad,function(x)x$seg$type)=='ALT',1,sum)),log='y',main='Orthologous cassette exons',xlab='number of species',ylab='number of exons')
barplotWithText.(table(apply(!is.na(orth.dds.all$dd),1,sum)),log='y',main='Orthologous alternative donors',xlab='number of species',ylab='number of exons')
barplotWithText.(table(apply(!is.na(orth.aas.all$aa),1,sum)),log='y',main='Orthologous alternative acceptors',xlab='number of species',ylab='number of exons')
barplotWithText.(table(table(ir.orth.inxs)),log='y',main='Orthologous retained introns',xlab='number of species',ylab='number of exons')
m = meta
SAJR::plotMDS(points=-mds$ad.all,col=m$col,cex=m$cex,pch=m$pch,main=paste0('Cassette exons (',nrow(orth.seg.ad$human$seg),')'))
SAJR::plotMDS(points=mds$dd.all,col=m$col,cex=m$cex,pch=m$pch,main=paste0('Donors (',nrow(orth.dds$human$seg),')'))
SAJR::plotMDS(points=mds$aa.all,col=m$col,cex=m$cex,pch=m$pch,main=paste0('Acceptors (',nrow(orth.aas$human$seg),')'))
SAJR::plotMDS(points=mds$ir.all,col=m$col,cex=m$cex,pch=m$pch,main=paste0('Retained introns (',nrow(orth.irs$human$seg),')'))
f = apply(sapply(orth.seg.ad,function(x)x$seg$type=='ALT'),1,sum)
SAJR::plotMDS(points=mds$ad.anc,col=m$col,cex=m$cex,pch=m$pch,main=paste0('Ancient cassette exons (',sum(f==7),')'))

plot.new()
sp = unique(meta[,c('species','pch')])
ti = unique(meta[,c('tissue','col')])
legend('topleft',pch=c(rep(19,7),sp$pch),col=c(ti$col,rep('black',7)),legend=c(ti$tissue,sp$species),ncol=2)
dev.off()
#+ кратность 3 от кол-ва видов. 
#?? положение новорожденных aa и dd - внутри/вне экзона
###########################
##fig 5 alternification ###
prop.exp = sapply(orth.seg.ad.tsm,function(x)apply(!is.na(x[,meta.tsm[colnames(x),'tissue']=='brain']),1,mean))
min.psi = sapply(orth.seg.ad.tsm,function(x)apply(x[,meta.tsm[colnames(x),'tissue']=='brain'],1,min,na.rm=T))
table(apply(prop.exp,1,min)>0.9)
h=hist(min.psi[,1],0:1000/1000,xlim=c(0.6,1))
barplot(log(h$counts))

spspb=apply(min.psi<0.8,1,function(x){paste(species$short[x],collapse='')})
sort(table(spspb[apply(prop.exp,1,min)>0.5]))

# fig 5 and 6: alternification and exonification: brain and testis 