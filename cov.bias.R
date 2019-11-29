library(SAJR)
library(xlsx)
source('code/r.functions/load.all.data.F.R')
options(stringsAsFactors = FALSE)
species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
ens.ge.cod = readRDS('Rdata/ens.ge.cod.Rdata')
age.al.i = read.table('input/margarida.age.alignment.final.2017.05.17.improved.tab',sep='\t',header = T)[,c(rownames(species),'to.remove')]

# load RQN
# lib.fields = c('library.ID','tube.with.tissue.ID..source.')
# rna.fields = c('Tube.with.tissue.ID..source.','RQN')
# fls = list.files('../evo.devo.pilot/input/sample.info/','HKDB.+.xlsx')
# fls = fls[c(-1,-4,-8)]
# fls = c('/home/mazin/skoltech/projects/evo.devo.pilot/input/HKDB MACAQUE uploaded 30.04.15.xlsx',paste0('../evo.devo.pilot/input/sample.info/',fls))
# 
# lib = do.call(rbind,lapply(fls,function(f)read.xlsx(f,'libraries')[,lib.fields]))
# rna = do.call(rbind,lapply(fls,function(f)read.xlsx(f,'DNA-RNA')[,rna.fields]))
# lib$library.ID = tolower(lib$library.ID)
# lib$tube.with.tissue.ID..source. = tolower(lib$tube.with.tissue.ID..source.)
# rna$Tube.with.tissue.ID..source. = tolower(rna$Tube.with.tissue.ID..source.)
# rna$RQN = as.numeric(rna$RQN)
# lib = lib[lib$library.ID %in% tolower(meta$lib.id) & lib$tube.with.tissue.ID..source. %in% rna$Tube.with.tissue.ID..source.,]
# rna = rna[rna$Tube.with.tissue.ID..source. %in% lib$tube.with.tissue.ID..source.,]
# dim(lib)
# dim(rna)
# t = table(rna$Tube.with.tissue.ID..source.)
# rna[rna$Tube.with.tissue.ID..source. == names(t[t>1])[1],]
# 
# rna=sapply(split(rna$RQN,rna$Tube.with.tissue.ID..source.),mean)
# length(rna)
# rna = rna[lib$tube.with.tissue.ID..source.]
# names(rna) = lib$library.ID
# meta$RQN = rna[tolower(meta$lib.id)]
# table(meta$species,is.na(meta$RQN))
# saveRDS(meta,'Rdata/main.set.meta.Rdata')
h = meta[meta$species=='human',]
m = meta[meta$species=='mouse',]
# load bias
m.cov.bias = read.table('processed/cov.bias/mouse.1000.0.1.false.false.true.-1.outs',sep='\t')
rownames(m.cov.bias) = sapply(strsplit(m.cov.bias$V1,'/',fixed = T),function(x){strsplit(x[11],'.',T)[[1]][1]})
m.cov.bias = as.matrix(m.cov.bias[m$lib.id,-1])
m.cov.bias = t(apply(m.cov.bias,1,function(x){(x-min(x))/(max(x)-min(x))}))
mcb1 = apply(m.cov.bias,1,mean)[m$lib.id]
mcb2 = apply(m.cov.bias,1,function(x){mean((x-min(x))/(max(x)-min(x)))})[m$lib.id] #now it is equal to mcb1

h.cov.bias = read.table('processed/cov.bias/human.1000.0.1.false.false.true.-1.outs',sep='\t')
rownames(h.cov.bias) = sapply(strsplit(h.cov.bias$V1,'/',fixed = T),function(x){strsplit(x[11],'.',T)[[1]][1]})
h.cov.bias = as.matrix(h.cov.bias[h$lib.id,-1])
h.cov.bias = t(apply(h.cov.bias,1,function(x){(x-min(x))/(max(x)-min(x))}))

hcb1 = apply(h.cov.bias,1,mean)[h$lib.id]
hcb2 = apply(h.cov.bias,1,function(x){mean((x-min(x))/(max(x)-min(x)))})[h$lib.id] #now it is equal to hcb1

# load marg expr
ge.h.marg=read.table('processed/GE.from.marg/HumanRpkmMajorTissuesCor90.Norm.txt',header = 1,row.names = 1)
hgt = read.csv('input/gene.info.from.marg/HumanDDGs.Biotype.csv')
ge.h.marg = as.matrix(ge.h.marg[,h$marg.name])
ge.h.marg.cor = cor(log(ge.h.marg[rownames(ge.h.marg) %in% hgt$Human_ID[hgt$Biotype=='protein_coding'],]+1),m='pear',u='p')
#ge.h.my.cor = cor(log(ens.ge.cod$human$rpkm+1),m='pear',u='p')

ge.m.marg=read.table('processed/GE.from.marg/MouseRpkmMajorTissuesCor90.Norm.txt',header = 1,row.names = 1)
mgt = read.csv('input/gene.info.from.marg/MouseDDGs.Biotype.csv')
ge.m.marg = as.matrix(ge.m.marg[,m$marg.name])
ge.m.marg.cor = cor(log(ge.m.marg[rownames(ge.m.marg) %in% mgt$Mouse_ID[mgt$Biotype=='protein_coding'],]+1),m='pear',u='p')

h.best.marg.cor = getBestCor(ge.h.marg.cor,h$tissue)
m.best.marg.cor = getBestCor(ge.m.marg.cor,m$tissue)

#my GE
ge.h.my.cor = cor(log(ens.ge.cod$human$rpkm[,rownames(h)]+1),m='pear',u='p')
h.best.my.cor = getBestCor(ge.h.my.cor,h$tissue)

ge.m.my.cor = cor(log(ens.ge.cod$mouse$rpkm[,rownames(m)]+1),m='pear',u='p')
m.best.my.cor = getBestCor(ge.m.my.cor,m$tissue)
# AS
#human.ad.cor = cor(hdata$ir[hdata$seg$sites=='ad',],u='p')
human.ad.cor = readRDS('Rdata/human.ad.cor.Rdata')
h.best.cor.as = getBestCor(human.ad.cor[rownames(h),rownames(h)],h$tissue)
mouse.ad.cor = readRDS('Rdata/mouse.ad.cor.Rdata')
m.best.cor.as = getBestCor(mouse.ad.cor[rownames(m),rownames(m)],m$tissue)

hcb2.tsm=calcMeanCols(t(as.matrix(hcb2)),paste(h$species,h$tissue,h$stage),max)
h.best.marg.cor.tsm=calcMeanCols(t(as.matrix(h.best.marg.cor)),paste(h$species,h$tissue,h$stage),max)
h.best.my.cor.tsm=calcMeanCols(t(as.matrix(h.best.my.cor)),paste(h$species,h$tissue,h$stage),max)
h.best.cor.as.tsm = calcMeanCols(t(as.matrix(h.best.cor.as)),paste(h$species,h$tissue,h$stage),max)

mcb2.tsm=calcMeanCols(t(as.matrix(mcb2)),paste(m$species,m$tissue,m$stage),max)
m.best.marg.cor.tsm=calcMeanCols(t(as.matrix(m.best.marg.cor)),paste(m$species,m$tissue,m$stage),max)
m.best.my.cor.tsm=calcMeanCols(t(as.matrix(m.best.my.cor)),paste(m$species,m$tissue,m$stage),max)
m.best.cor.as.tsm = calcMeanCols(t(as.matrix(m.best.cor.as)),paste(m$species,m$tissue,m$stage),max)



### plot coverage bias
pdf('figures/cov.bias.pdf',h=6,w=12)
par(mfcol=c(2,4),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
plot(1,t='n',xlim=c(00,100),ylim=c(0,1),xlab="5'->3' (%)",ylab='Mean relative coverage',main='Mouse')
cols=getPal(c('red','gray','green'),nrow(m.cov.bias))[rank(mcb2)]
for(i in 1:nrow(m.cov.bias))
	lines(m.cov.bias[i,],col=cols[i])

plot(1,t='n',xlim=c(00,100),ylim=c(0,1),xlab="5'->3' (%)",ylab='Mean relative coverage',main='Human')
cols=getPal(c('red','gray','green'),nrow(h.cov.bias))[rank(hcb2)]
for(i in 1:nrow(h.cov.bias))
	lines(h.cov.bias[i,],col=cols[i])
plotLine(m$RQN,mcb2,xlab='RQN',ylab='Coverage bias',col=m$col,cex=m$cex*2,pch=19,leg.pos='bottomright',main='RQN vs Coverage bias')
plotLine(h$RQN,hcb2,xlab='RQN',ylab='Coverage bias',col=h$col,cex=h$cex*2,pch=19,leg.pos='bottomright',main='RQN vs Coverage bias')
# 
# plotLine(m$RQN,m.best.my.cor,xlab='RQN',ylab='max(corr. to same tissue)',col=m$col,cex=m$cex*2,pch=19,leg.pos='bottomright',main='Coverage bias correlates with RQN')
# plotLine(h$RQN,h.best.my.cor,xlab='RQN',ylab='max(corr. to same tissue)',col=m$col,cex=m$cex*2,pch=19,leg.pos='bottomright',main='Coverage bias correlates with RQN')
# plotLine(mcb2 ,m.best.my.cor,xlab='Coverage bias',ylab='max(corr. to same tissue)',col=m$col,cex=m$cex*2,pch=19,leg.pos='bottomright',main='Coverage bias correlates with RQN')
# plotLine(hcb2 ,h.best.my.cor,xlab='Coverage bias',ylab='max(corr. to same tissue)',col=m$col,cex=m$cex*2,pch=19,leg.pos='bottomright',main='Coverage bias correlates with RQN')

plotLine(m$RQN,m.best.cor.as,xlab='RQN',ylab='max(corr. to same tissue)',col=m$col,cex=m$cex*2,pch=19,leg.pos='bottomright',main='RQN vs best-AS-corr.')
plotLine(h$RQN,h.best.cor.as,xlab='RQN',ylab='max(corr. to same tissue)',col=m$col,cex=m$cex*2,pch=19,leg.pos='bottomright',main='RQN vs best-AS-corr.')
plotLine(mcb2 ,m.best.cor.as,xlab='Coverage bias',ylab='max(corr. to same tissue)',col=m$col,cex=m$cex*2,pch=19,leg.pos='bottomright',main='Coverage bias vs best-AS-corr.')
plotLine(hcb2 ,h.best.cor.as,xlab='Coverage bias',ylab='max(corr. to same tissue)',col=m$col,cex=m$cex*2,pch=19,leg.pos='bottomright',main='Coverage bias vs best-AS-corr.')


plotTissueAgeProile(mcb2.tsm[1,],meta.tsm[meta.tsm$stage %in% age.al.i$mouse,],age.axis = 'rank',df = -1,pch=19,cex=1,ylab='max coverage bias',main='Max coverage bias per stage',ylim=c(0.4,0.9))
plotTissueAgeProile(hcb2.tsm[1,],meta.tsm[meta.tsm$stage %in% age.al.i$human,],age.axis = 'rank',df = -1,pch=19,cex=1,ylab='max coverage bias',main='Max coverage bias per stage',ylim=c(0.4,0.9))

plotTissueAgeProile(m.best.cor.as.tsm[1,],meta.tsm[meta.tsm$stage %in% age.al.i$mouse,],age.axis = 'rank',df = -1,pch=19,cex=1,ylab='best-AS-corr.',main='Best-AS-corr. per stage',ylim=c(0.925,0.98))
plotTissueAgeProile(h.best.cor.as.tsm[1,],meta.tsm[meta.tsm$stage %in% age.al.i$human,],age.axis = 'rank',df = -1,pch=19,cex=1,ylab='best-AS-corr.',main='Best-AS-corr. per stage',ylim=c(0.925,0.98))

plotTissueAgeProile(m.best.marg.cor.tsm[1,],meta.tsm[meta.tsm$stage %in% age.al.i$mouse,],age.axis = 'rank',df = -1,pch=19,cex=1,ylab='best-GE-corr.',main='Best-GE-corr. per stage',ylim=c(0.945,1))
plotTissueAgeProile(h.best.marg.cor.tsm[1,],meta.tsm[meta.tsm$stage %in% age.al.i$human,],age.axis = 'rank',df = -1,pch=19,cex=1,ylab='best-GE-corr.',main='Best-GE-corr. per stage',ylim=c(0.945,1))

# plotTissueAgeProile(m.best.my.cor.tsm[1,],meta.tsm[meta.tsm$stage %in% age.al.i$mouse,],age.axis = 'rank',df = -1,pch=19,cex=1,ylab='best-GE-corr.',main='Best-GE-corr. per stage',ylim=c(0.945,1))
# plotTissueAgeProile(h.best.my.cor.tsm[1,],meta.tsm[meta.tsm$stage %in% age.al.i$human,],age.axis = 'rank',df = -1,pch=19,cex=1,ylab='best-GE-corr.',main='Best-GE-corr. per stage',ylim=c(0.945,1))
plot.new()
plot.new()
plotLine(m$RQN,m.best.marg.cor,xlab='RQN',ylab='max(corr. to same tissue)',col=m$col,cex=m$cex*2,pch=19,leg.pos='bottomright',main='RQN vs best-GE-corr.')
plotLine(h$RQN,h.best.marg.cor,xlab='RQN',ylab='max(corr. to same tissue)',col=m$col,cex=m$cex*2,pch=19,leg.pos='bottomright',main='RQN vs best-GE-corr.')
plotLine(mcb2 ,m.best.marg.cor,xlab='Coverage bias',ylab='max(corr. to same tissue)',col=m$col,cex=m$cex*2,pch=19,leg.pos='bottomright',main='Coverage bias vs best-GE-corr.')
plotLine(hcb2 ,h.best.marg.cor,xlab='Coverage bias',ylab='max(corr. to same tissue)',col=m$col,cex=m$cex*2,pch=19,leg.pos='bottomright',main='Coverage bias vs best-GE-corr.')
dev.off()

# compare mcb1 and mcb2
plot(mcb1,mcb2)
abline(h=quantile(mcb2,0.2))
abline(v=quantile(mcb1,0.2))

f = paste(ifelse(mcb1<quantile(mcb1,0.25),'1',''),ifelse(mcb2<quantile(mcb2,0.25),'2',''))
table(f)
boxplotWithSgn(split(m$RQN,f))

m.cov.bias. = t(apply(m.cov.bias,1,function(x){(x-min(x))/(max(x)-min(x))}))
plot(1,t='n',xlim=c(00,100),ylim=range(m.cov.bias.))
for(i in 1:nrow(m.cov.bias))
	lines(m.cov.bias.[i,],col='gray')
for(i in which(f=='1 '))
	lines(m.cov.bias.[m$lib.id[i],],col='red')
for(i in which(f=='1 2'))
	lines(m.cov.bias.[m$lib.id[i],],col='green')




mt=cbind(mcb1,mcb2,RQN=m$RQN)#,m.best.cor,m.best.cor.as)
cor(mt[mcb1 > 0.4,],u='p',m='sp')
pairs(mt[mcb1 > 0.4,],col=m$col,pch=19,cex=m$cex)






