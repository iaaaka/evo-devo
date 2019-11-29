#setwd('~/skoltech/projects/evo.devo/')
options(stringsAsFactors = FALSE)
source('code/r.functions/paper.figures.F.R')
source('code/r.functions/load.all.data.F.R')
library(SAJR)


species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
orth.seg.ad.tsm = readRDS('Rdata/orth.seg.ad.tsm.Rdata')
orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
ens.ge.cod.tsm = readRDS('Rdata/ens.ge.cod.tsm.Rdata')

age.al.i = read.table('input/margarida.age.alignment.final.2017.05.17.improved.tab',sep='\t',header = T)[,c(rownames(species),'to.remove')]
params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)
orth.seg.ad = lapply(orth.seg.ad,function(x)x[,colnames(x$ir) %in% rownames(meta)])
orth.ens.genes = readRDS('Rdata/orth.ens.genes.Rdata')


f = apply(sapply(orth.seg.ad,function(x)x$seg$type=='ALT'),1,sum)==7
f = apply(sapply(orth.seg.ad,function(x)x$seg$length<300),1,sum)==7
table(f)

tissue.stage.cor = lapply(unique(meta$tissue),function(t){
	r=lapply(1:nrow(age.al.i),function(i)getSpeciesCor(orth.seg.ad.tsm,age.al.i[i,c(1:7)],t,f=f))
	names(r) = age.al.i[,'mouse']
	r
	})


f = sapply(names(ens.ge.cod.tsm),function(s)orth.ens.genes[,s] %in% rownames(ens.ge.cod.tsm[[s]]))
f = apply(f,1,sum)==7
ens.ge.cod.tsm.log = lapply(names(ens.ge.cod.tsm),function(s)log2(ens.ge.cod.tsm[[s]][orth.ens.genes[f,s],]+0.1))
names(ens.ge.cod.tsm.log) = names(ens.ge.cod.tsm)

tissue.stage.cor.ge = lapply(unique(meta$tissue),function(t){
	r=lapply(1:nrow(age.al.i),function(i)getSpeciesCor(ens.ge.cod.tsm.log,age.al.i[i,c(1:7)],t))
	names(r) = age.al.i[,'mouse']
	r
})

tissue.stage.cor.sp = lapply(unique(meta$tissue),function(t){
	r=lapply(1:nrow(age.al.i),function(i)getSpeciesCor(orth.seg.ad.tsm,age.al.i[i,c(1:7)],t,method = 'sp'))
	names(r) = age.al.i[,'mouse']
	r
})
names(tissue.stage.cor.ge) = names(tissue.stage.cor.sp) = names(tissue.stage.cor) = unique(meta$tissue)


plotDivergenceOnAge(tissue.stage.cor,species=c(1,3:6))
image(tissue.stage.cor$liver$`15.5`[-2,-2])
tissue.stage.cor$testis$`15.5`

samids = rownames(meta[meta$species=='human' & meta$tissue=='liver' & meta$stage=='8wpc',])

t = 'testis'
s2 = 'rabbit'
m = meta[colnames(orth.seg.ad$mouse$ir),]
m = orth.seg.ad$mouse$ir[,m$tissue==t & m$stage==14.5]
hm = meta[colnames(orth.seg.ad[[s2]]$ir),]
hm = hm[hm$tissue==t,]
hm = hm[order(hm$days),]
h = orth.seg.ad[[s2]]$ir[,rownames(hm)]
cc = cor(m,h,u='p')
c = factor(meta[colnames(h),'stage'])
plot(cc[1,],col=c,pch=ifelse(hm$stage=='8wpc',19,1))
points(cc[2,],col=c,pch=ifelse(hm$stage=='8wpc',19,1))


good = bad = NULL
good = c(good,colnames(cc)[cc[1,]>0.7])
bad = c(bad,colnames(cc)[cc[1,]<0.65])

library(xlsx)
hs = read.xlsx('../evo.devo.pilot/input/sample.info/HKDB HUMAN FINAL VERSION.xlsx','libraries')
hs[hs$library.ID %in% meta[bad,'lib.id'],]

t = unique(hs[,c(1,3,7)])
t = t[!is.na(t$library.ID),]
rownames(t) =t$library.ID

mh = meta[meta$lib.id %in% t$library.ID,]
t = t[mh$lib.id,]
f = apply(sapply(orth.seg.ad,function(x)x$seg$type=='ALT'),1,sum)==7

p = orth.seg.ad$human$ir[orth.seg.ad$human$seg$type=='ALT',rownames(mh)]
p = orth.seg.ad$human$ir[f,rownames(mh)]

p. = readRDS('Rdata/human.as.u.filtered.Rdata')
p = p.$ir[p.$seg$sites=='ad' & p.$seg$cod=='c',rownames(mh)]

p = log2(readRDS('Rdata/ens.ge.cod.Rdata')$human$rpkm[,rownames(mh)]+0.1)
col=setNames(c('red','green','blue','cyan'),unique(t$performed.by))

par(mfcol=c(3,3),tck=-0.01,mgp=c(1.1,0.2,0),mar=c(6,3,1.5,0),oma=c(0,0,2,1))
for(ts in unique(meta$tissue)){
	f = mh$tissue==ts & mh$days<3000
	mds = cmdscale(1-cor(p[,f],u='p'),k=2)
	plot(mds,cex=mh$age.rank[f]/max(mh$age.rank[f])*3,pch=19,col=col[t$performed.by[f]],main=ts)
}

mds = cmdscale(1-cor(p,u='p'),k=4)
plot(mds[,1:2],cex=mh$age.rank/max(mh$age.rank)*3,,col=col[t$performed.by],pch=19)
plot(mds[,3:4],cex=mh$age.rank/max(mh$age.rank)*3,,col=col[t$performed.by],pch=19)
plot(mds[,1:2],cex=mh$age.rank/max(mh$age.rank)*3,,col=mh$col,pch=19)
plot(mds[,3:4],cex=mh$age.rank/max(mh$age.rank)*3,,col=mh$col,pch=19)

table(t$performed.by)
boxplot(log(mh$days)~ t$performed.by,las=3)

apply(orth.seg.ad$human$ir[,samids],2,cor,y=orth.seg.ad.tsm$mouse[,'mouse liver 14.5'],u='p')
table(mh$tissue,t$performed.by)



orth.seg.ad. = readRDS('Rdata/old/orth.seg.ad.Rdata')
meta. = readRDS('Rdata/old/main.set.meta.Rdata')

s1 = 'mouse'
s2 = 'mouse'



f = apply(sapply(orth.seg.ad,function(x)x$seg$type=='ALT'),1,sum)==7
f=T

m = meta
m1 = m[m$species==s1 & m$tissue==t,]
m1 = rownames(m1)[order(m1$days)]
m2 = m[m$species==s2 & m$tissue==t,]
m2 = rownames(m2)[order(m2$days)]
c=cor(orth.seg.ad[[s1]]$ir[f,m1],orth.seg.ad[[s2]]$ir[f,m2],u='p')
image(1:length(m1),1:length(m2),c,yaxt='n',xaxt='n',xlab='',ylab='')
axis(1,1:length(m1),m1,las=2)
axis(2,1:length(m2),m2,las=2)

par(mfrow=c(2,2))
plot(orth.seg.ad.tsm$human[,'human brain toddler'],orth.seg.ad.tsm$mouse[,'mouse brain 2wpb'])#,u='p',m='sp')
plotLine(orth.seg.ad.tsm$human[,'human liver 7wpc'],orth.seg.ad.tsm$opossum[,'opossum liver 2dpb'])#,u='p',m='sp')
plotLine(orth.seg.ad.tsm$human[,'human liver 8wpc'],orth.seg.ad.tsm$opossum[,'opossum liver 4dpb'])#,u='p',m='sp')
plotLine(orth.seg.ad.tsm$human[,'human liver 7wpc'],orth.seg.ad.tsm$human[,'human liver 11wpc'])#,u='p',m='sp')
pairs(orth.seg.ad$human$ir[,c('HL53F_1','HL53M_1','HL56F_1','HL56M_1','HL51M_1','HL49F_1')],pch='.')

c=cor(orth.seg.ad.tsm$human[,paste('human liver',age.al.i$human)],u='p')
image(c)
plot(c[5,])
meta[meta$species=='human' & meta$stage=='7wpc',]

c = cor(orth.seg.ad$human$ir,u='p')
m = meta[meta$species=='human',]
m = m[order(m$tissue,m$days),]
c = c[rownames(m),rownames(m)]
plot(c[m$tissue=='testis' & m$stage=='7wpc',][1,],col=meta[colnames(c),'col'],ylim=c(0.75,1))
