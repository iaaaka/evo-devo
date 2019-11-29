library(SAJR)
options(stringsAsFactors = FALSE)
source('code/r.functions/load.all.data.F.R')
orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
ens.ge.cod.tsm = readRDS('Rdata/ens.ge.cod.tsm.Rdata')
ens.ge.cod = readRDS('Rdata/ens.ge.cod.Rdata')
orth.ens.genes = readRDS('Rdata/orth.ens.genes.Rdata')
ens.ge.cod = readRDS('Rdata/ens.ge.cod.Rdata')
orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
h.as = readRDS('Rdata/human.as.u.filtered.Rdata')
alts.filt = readRDS('Rdata/alts.filt.Rdata')
f = sapply(names(ens.ge.cod.tsm),function(s)orth.ens.genes[,s] %in% rownames(ens.ge.cod.tsm[[s]]))
f = apply(f,1,sum)==7
ens.ge.cod.tsm.log = lapply(names(ens.ge.cod.tsm),function(s)log2(ens.ge.cod.tsm[[s]][orth.ens.genes[f,s],]+0.1))
names(ens.ge.cod.tsm.log) = names(ens.ge.cod.tsm)
all.anns = readRDS('Rdata/all.anns.Rdata') 

# for hcb2 see cov.bias.R
sids = rownames(meta)[meta$tissue=='testis' & meta$species=='human' & meta$stage %in% c('7wpc','8wpc','11wpc')]
hcb2[meta[sids,'lib.id']]
sids = sids[-c(2,7)]

ir = h.as$ir
ir[h.as$i+h.as$e < 50 ] = NA

psi = ir[h.as$seg$sites=='ad' & h.as$seg$alt.size<6,sids]
image(cor(psi,u='p',m='p'))
pairs(psi,pch='.')
cor(psi,u='p',m='sp')
plot(psi[,1],psi[,3])
which(!is.na(psi[,1]+psi[,3]) & psi[,1]>0.6 & psi[,3]<0.2)
h.as['hum.12777.s7',sids]
all.anns$human[all.anns$human$gene_id=='hum.37841',]

pdf('figures/cov.bias.example.hum.12777.s7.pdf',h=10,w=6)
par(mfrow=c(6,1),tck=-0.01,mgp=c(1.1,0.2,0),mar=c(1,1,1,1),oma=c(0,0,0,0))
for(i in sids){
	fs=paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/',meta[i,'fname'],'.bam')
	getReadCoverage(fs,'12',4384039,4414631,NA,plot = T,min.junc.cov = 5,main=paste(meta[i,'stage'],'(',i,')'),ylim=c(0,350))
}
dev.off()
#fs=paste0('/home/mazin/skoltech/projects/evo.devo/raw/bams.from.lausanne/',ht[ht$stage==st,'fname'],'.sorted.bam')

ss=alts.filt$human[h.as$seg$alt.id[h.as$seg$sites=='ad'],'sites']
sort(table(ss),decreasing = T)[1:20]


geneExNo = sapply(split(all.anns$human$exon.number,all.anns$human$gene_id),max)
f = h.as$seg$sites=='ad' & geneExNo[h.as$seg$gene_id] > 10
rpos = h.as$seg$exon.number/geneExNo[h.as$seg$gene_id]
plot(h.as$ir[f,sids[1]],h.as$ir[f,sids[3]])
hist(h.as$ir[f,sids[1]]-h.as$ir[f,sids[3]])
plot(rpos[f],(h.as$ir[f,sids[1]]-h.as$ir[f,sids[3]]),col='#00000020')



h=rownames(meta)[meta$tissue=='testis' & meta$species=='human' & meta$stage=='8wpc']
m=rownames(meta)[meta$tissue=='testis' & meta$species=='mouse' & meta$stage=='14.5']

h1=rownames(meta)[meta$tissue=='testis' & meta$species=='human' & meta$stage=='7wpc']
m1=rownames(meta)[meta$tissue=='testis' & meta$species=='mouse' & meta$stage=='13.5']

h2=rownames(meta)[meta$tissue=='testis' & meta$species=='human' & meta$stage=='11wpc']



cor(cbind(orth.seg.ad$human$ir[,h1],orth.seg.ad$mouse$ir[,m1]),u='p',m='sp')
cor(cbind(orth.seg.ad$human$ir[,h],orth.seg.ad$mouse$ir[,m]),u='p',m='sp')

image(cor(orth.seg.ad$human$ir[,c(h1,h,h2)],u='p',m='p'))

pairs(cbind(orth.seg.ad$human$ir[,h],orth.seg.ad$mouse$ir[,m]))
pairs(cbind(orth.seg.ad$human$ir[,h1],orth.seg.ad$human$ir[,h]))

boxplot(orth.seg.ad$human$ir[,h],outline=F)
apply(orth.seg.ad$human$ir[,h1]<0.9,2,sum,na.rm=T)

hist(orth.seg.ad$human$ir[,h],0:100/100,ylim=c(0,1000))


ens.ge.cod.tsm.log

ens.ge.cod.tsm.log$human[1:2,1:6]

htsm = c('human testis 7wpc','human testis 8wpc','human testis 11wpc')
plot(ens.ge.cod.tsm.log$human[,'human testis 8wpc'],ens.ge.cod.tsm.log$mouse[,'mouse testis 14.5'])
plot(ens.ge.cod.tsm.log$human[,'human testis 7wpc'],ens.ge.cod.tsm.log$mouse[,'mouse testis 13.5'])
plot(ens.ge.cod.tsm.log$human[,'human testis 11wpc'],ens.ge.cod.tsm.log$mouse[,'mouse testis 15.5'])

pairs(ens.ge.cod.tsm.log$human[,htsm])
plot(ens.ge.cod.tsm.log$human[,'human testis 11wpc'],ens.ge.cod.tsm.log$human[,'human testis 8wpc']-ens.ge.cod.tsm.log$human[,'human testis 11wpc'])

par(mfrow=c(2,2))
hist(ens.ge.cod.tsm.log$human[,'human testis 7wpc']-ens.ge.cod.tsm.log$human[,'human testis 11wpc'],-30:30/3,xlim=c(-5,5))

meta[meta$stage=='8wpc' & meta$tissue=='testis' ,]



which(abs(ens.ge.cod.tsm.log$human[,'human testis 8wpc']-ens.ge.cod.tsm.log$human[,'human testis 7wpc'])>4)
plotTissueAgeProile(ens.ge.cod.tsm$human['ENSG00000160710',],meta.tsm,tissues = 'testis',age.axis = 'rank',cex=2,pch=19)
plotTissueAgeProile(ens.ge.cod$human$rpkm['ENSG00000160710',],meta,tissues = 'testis',age.axis = 'rank',cex=2,pch=19)


marg.exp1=read.table('processed/GE.from.marg/HumanRpkmMajorTissuesCor90.Norm.txt',header = 1,row.names = 1)
marg.expc=as.matrix(read.table('processed/GE.from.marg/HumanCountsMajorTissuesCor90.Norm.txt',header = 1,row.names = 1))
marg.exp=read.table('processed/GE.from.marg/HumanCPM.txt',header = 1,row.names = 1)
marg.exp = as.matrix(marg.exp)
marg.exp1 = as.matrix(marg.exp1)

plot(marg.exp1['ENSG00000160710',grep('Testis',colnames(marg.exp1))])


plot(marg.exp['ENSG00000160710',115:134])

h=meta.tsm[meta.tsm$species=='human' & meta.tsm$tissue=='testis', ]

mt = colnames(marg.exp)[grep('Testis',colnames(marg.exp))]
pt = rownames(h)[order(h$days)]
cbind(mt,pt)

cmn = intersect(rownames(marg.exp),rownames(ens.ge.cod.tsm$human))
plot(marg.exp[cmn,mt[4]],ens.ge.cod.tsm$human[cmn,pt[4]]/marg.exp[cmn,mt[4]],log='xy',pch='.',col=factor(ens.ge.cod$human$gene[cmn,'strand']))

t=sapply(1:length(mt),function(i)cor(marg.exp[cmn,mt[i]],ens.ge.cod.tsm$human[cmn,pt[i]],u='p',m='sp'))
plot(t)

table(ens.ge.cod.tsm$human[cmn,pt[5]]/marg.exp[cmn,mt[5]] < 0.01)

mr=cor(marg.exp[cmn,mt],u='p',m='sp')
pac=cor(ens.ge.cod.tsm$human[cmn,pt],u='p',m='sp')

image(mr)
hist(mr/pac)

z = cor(marg.exp[cmn,mt],ens.ge.cod.tsm$human[cmn,pt],u='p',m='sp')
z1 = cor(marg.exp1[cmn,mt1],ens.ge.cod$human$rpkm[cmn,pt1],u='p',m='sp')
image(1:ncol(z),z = z)

plot(1:39,t='n',ylim=c(0.5,1))
for(i in setdiff(colnames(z1)[1:20],h))
	lines(z1[,i])
for(i in h) lines(z1[,i],col='red',lwd=3)
	

plot(marg.exp['ENSG00000160710',mt],ens.ge.cod.tsm$human['ENSG00000160710',pt],col=(pt=='human testis 8wpc')+1,pch=19,log='xy')
plot(marg.exp['ENSG00000160710',mt],pch=19,col='red',ylim=c(1,400),log='y')
points(ens.ge.cod.tsm$human['ENSG00000160710',pt],pch=19,col='blue',ylim=c(0,400))


h=meta[meta$species=='human' & meta$tissue=='testis', ]
mt1 = colnames(marg.exp1)[grep('Testis',colnames(marg.exp1))]
pt1 = rownames(h)[order(h$days)]
cbind(mt1,pt1)

plot(marg.exp1['ENSG00000160710',mt1],ens.ge.cod$human$rpkm['ENSG00000160710',pt1],col=(pt=='human testis 8wpc')+1,pch=19,t='l')

plot(marg.exp1['ENSG00000160710',mt1],pch=19,col='red',ylim=c(1,100),log='')
points(ens.ge.cod$human$rpkm['ENSG00000160710',pt1],pch=19,col='blue')
abline(v=0:40,lty=3)

which(cmn=='ENSG00000160710')
points(marg.exp['ENSG00000160710',mt[5]],ens.ge.cod.tsm$human['ENSG00000160710',pt[5]]/marg.exp['ENSG00000160710',mt[5]],pch=19,col='blue')

b=c('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/3588sTS.Human.Heart.CS13.Male.bam',
		'/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/3671sTS.Human.Heart.CS16.Male.bam',
		'/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/2088sTS.Human.Heart.CS18.Male.bam')

meta[h,]

b = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/',meta[h,'fname'],'.bam')
b1 = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/',meta[h1,'fname'],'.bam')
b2 = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/',meta[h2,'fname'],'.bam')
HT46M_1=getReadCoverage(b[2],'1',154554538,154560000,T,min.junc.cov = 5)
HT46M_1=getReadCoverage(b1,'1',154544538,154600475,T,min.junc.cov = 5)
HT46M_1=getReadCoverage(b2,'1',154544538,154600475,T,min.junc.cov = 5)

HT53M_1=getReadCoverage('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/1990sTS.Human.Testis.CS22.Male.bam','1',154544538,154600475,T)
HT56M_1=getReadCoverage('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/1991sTS.Human.Testis.CS23.Male.bam','1',154544538,154600475,T)
HT56M_2=getReadCoverage('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/2020sTSm.Human.Testis.8w.Male.bam','1',154544538,154600475,T)

ens.h.segs = loadSAData('processed/annotation/all.species/ensambl/Homo_sapiens.GRCh37.73.sajr')

e = ens.h.segs$seg[ens.h.segs$seg$gene_id=='ENSG00000160710' & ens.h.segs$seg$type=='EXN',]
rect(e$start,-20,e$stop,0)



ht = meta[meta$species=='human' & meta$tissue=='testis',]
ht = ht[order(ht$days),]

s = unique(ht$stage)
ens.ge.cod$human$gene['ENSG00000160710',]
par(mfrow=c(7,1),tck=-0.01,mgp=c(1.1,0.2,0),mar=c(1,1,1,1),oma=c(0,0,0,0))
for(st in s[1:7]){
	#fs=paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/',ht[ht$stage==st,'fname'],'.bam')
	fs=paste0('/home/mazin/skoltech/projects/evo.devo/raw/bams.from.lausanne/',ht[ht$stage==st,'fname'],'.sorted.bam')
	getReadCoverage(fs,'chr1',154554538,154600475,1,plot = T,min.junc.cov = 5,main=st)
}



e = ens.h.segs$seg[ens.h.segs$seg$type=='EXN',]


# cov.bias
gff = '/uge_mnt/home/mazin/projects/evo.devo/processed/annotation/all.species/ensambl/Mus_musculus.GRCm38.84.sajr'
bams = paste0('/uge_mnt/home/mazin/projects/evo.devo/processed/mapping/hisat2.s/mouse/',meta$fname[meta$species=='mouse'],'.bam')
write.table(cbind(gff,bams,100,1000,0.1,100,'false','false','true',-1),file = 'processed/cov.bias/mouse.in',sep='\t',quote = F,col.names = F,row.names = F)

gff = '/uge_mnt/home/mazin/projects/evo.devo/processed/annotation/all.species/ensambl/Homo_sapiens.GRCh37.73.sajr'
bams = paste0('/uge_mnt/home/mazin/projects/evo.devo/processed/mapping/hisat2.s/human/',meta$fname[meta$species=='human'],'.bam')
write.table(cbind(gff,bams,100,1000,0.1,100,'false','false','true',-1),file = 'processed/cov.bias/human.in',sep='\t',quote = F,col.names = F,row.names = F)


## load RQN
# load human
library(xlsx)
lib = read.xlsx('../evo.devo.pilot/input/sample.info/HKDB HUMAN FINAL VERSION.xlsx','libraries')
dna.rna = read.xlsx('../evo.devo.pilot/input/sample.info/HKDB HUMAN FINAL VERSION.xlsx','DNA-RNA')

table(lib$tube.with.tissue.ID..source. %in% dna.rna$Tube.with.tissue.ID..source.)

lib = lib[!is.na(lib$library.ID) & tolower(lib$library.ID) %in% tolower(meta$lib.id) & lib$tube.with.tissue.ID..source. %in% dna.rna$Tube.with.tissue.ID..source.,]
dna.rna = dna.rna[dna.rna$Tube.with.tissue.ID..source. %in% lib$tube.with.tissue.ID..source.,]
rownames(dna.rna) = dna.rna$Tube.with.tissue.ID..source.

dim(lib)
dim(dna.rna)

dna.rna = dna.rna[lib$tube.with.tissue.ID..source.,]
table(is.na(dna.rna$RIN),is.na(dna.rna$RQN))
rownames(dna.rna) = lib$library.ID

hmeta = meta[meta$species=='human',]

hmeta = cbind(hmeta,dna.rna[hmeta$lib.id,c('RIN','RQN')])
hmeta$RQN = as.numeric(hmeta$RQN)
hmeta$RIN = as.numeric(hmeta$RIN)
hmeta$has.info = hmeta$lib.id %in% lib$library.ID
table(hmeta$has.info,hmeta$stage)
hmeta[!hmeta$has.info,]
hmeta[h,]

# load mouse
lib = read.xlsx('../evo.devo.pilot/input/sample.info/HKDB_Mouse_17.03.2015.xlsx','libraries')
dna.rna = read.xlsx('../evo.devo.pilot/input/sample.info/HKDB_Mouse_17.03.2015.xlsx','DNA-RNA')
mmeta = meta[meta$species=='mouse',]
lib = lib[!is.na(lib$library.ID) & tolower(lib$library.ID) %in% tolower(meta$lib.id) & lib$tube.with.tissue.ID..source. %in% dna.rna$Tube.with.tissue.ID..source.,]
dna.rna = dna.rna[dna.rna$Tube.with.tissue.ID..source. %in% lib$tube.with.tissue.ID..source.,]
t=table(dna.rna$Tube.with.tissue.ID..source.)
dna.rna[dna.rna$Tube.with.tissue.ID..source. == names(t)[t>1],]
dna.rna = dna.rna[rownames(dna.rna)!='808',]
rownames(dna.rna) = dna.rna$Tube.with.tissue.ID..source.
dna.rna = dna.rna[lib$tube.with.tissue.ID..source.,]
rownames(dna.rna) = lib$library.ID
mmeta$RQN = as.numeric(dna.rna[mmeta$lib.id,'RQN'])


m.cov.bias = read.table('processed/cov.bias/mouse.1000.0.1.false.false.true.-1.outs',sep='\t')
rownames(m.cov.bias) = sapply(strsplit(m.cov.bias$V1,'/',fixed = T),function(x){strsplit(x[11],'.',T)[[1]][1]})

h.cov.bias = read.table('processed/cov.bias/human.1000.0.1.false.false.true.-1.outs',sep='\t')
rownames(h.cov.bias) = sapply(strsplit(h.cov.bias$V1,'/',fixed = T),function(x){strsplit(x[11],'.',T)[[1]][1]})

m.cov.bias = as.matrix(m.cov.bias[,-1])[mmeta$lib.id,]
h.cov.bias = as.matrix(h.cov.bias[,-1])[hmeta$lib.id,]

hmeta$cov.bias = apply(h.cov.bias,1,mean)[hmeta$lib.id]
mmeta$cov.bias = apply(m.cov.bias,1,mean)[mmeta$lib.id]

hmeta$cov.bias = apply(h.cov.bias,1,function(x){mean((x-min(x))/(max(x)-min(x)))})[hmeta$lib.id]
mmeta$cov.bias = apply(m.cov.bias,1,function(x){mean((x-min(x))/(max(x)-min(x)))})[mmeta$lib.id]


plot(hmeta$RQN,hmeta$cov.bias,col=hmeta$col,pch=19,cex=(hmeta$cex-0.3)*10)
plot(mmeta$RQN,mmeta$cov.bias,col=mmeta$col,pch=19,cex=(mmeta$cex-0.3)*3)

cor(hmeta$RQN,hmeta$cov.bias,u='p')
cor(mmeta$RQN,mmeta$cov.bias,u='p')

boxplot(hmeta$RQN~ hmeta$stage)
hist(hmeta$RQN)
abline(v=hmeta$RQN[hmeta$stage=='8wpc' & hmeta$tissue %in% c('testis','liver','ovary','kidney') ],col='red')


cb = h.cov.bias
cb = t(apply(cb,1,function(x)(x-min(x))/(max(x)-min(x))))

plot(1,t='n',xlim=c(00,100),ylim=range(0,1))
for(i in 1:nrow(cb))
	lines(cb[i,],col='gray')
#for(i in hmeta[hmeta$stage=='8wpc' & hmeta$tissue %in% c('testis','liver','kidney','ovary'),'lib.id'])
for(i in mmeta[!is.na(mmeta$RQN) & mmeta$RQN<6.5,'lib.id'])
	lines(cb[i,],col='red',lwd=3)

hist(apply(cov.bias1,2,mean),5:50/100,border=NA,col='gray')
hist(apply(cov.bias2,2,mean),5:50/100,border=NA,col='#FF000080',add=T)


hmeta$cov.bias2 = apply(cov.bias2,1,mean)[hmeta$lib.id]

cor.test(hmeta$RQN,hmeta$cov.bias2,u='p',m='p')

x = apply(cov.bias2,1,function(x){x[60]-x[80]})

plot(hmeta$RQN,hmeta$cov.bias2,col=hmeta$col,pch=19,cex=(hmeta$cex-0.3)*10)
plot(hmeta$RQN,x,col=hmeta$col,pch=19,cex=(hmeta$cex-0.3)*10)
cor.test(hmeta$RQN,hmeta$cov.bias2,m='sp')
cor.test(hmeta$RQN,x,m='sp')

points(hmeta[hmeta$stage=='8wpc','RQN'],hmeta[hmeta$stage=='8wpc','cov.bias2'],pch=19,col='red')
hmeta[h,]


hmeta[!is.na(hmeta$cov.bias) & hmeta$cov.bias<0.3,]

tt = hmeta$tissue=='testis'
plot(hmeta$age.rank[tt],hmeta$cov.bias[tt])

barplot(sapply(split(hmeta$cov.bias2[tt],hmeta$stage[tt]),max),las=2)

z=sapply(split(hmeta$cov.bias2,paste(hmeta$species,hmeta$tissue,hmeta$stage)),max)
plotTissueAgeProile(z,meta.tsm[meta.tsm$stage %in% age.al.i$human,],age.axis = 'rank',pch=19,cex=2,df = -1)


plotTissueAgeProile(setNames(hmeta$cov.bias2,rownames(hmeta)),meta,age.axis = 'rank',tissues = 'testis',pch=19,cex=2)
abline(v=1:30,lty=3)

plot(hmeta$age.rank[tt],hmeta$cov.bias2[tt],pch=1,cex=2,col='orange')


has = readRDS('Rdata/human.as.u.filtered.Rdata')
hmeta['HT51M_1',]

c=(has$i+has$e)[has$seg$sites=='ad',]
hist(log10(c),100,xlim=c(0,4))
abline(v=log10(50),col='red')
plotTissueAgeProile(apply(c>10,2,mean),meta,age.axis = 'rank',pch=19,cex=2,tissues = 'testis')

zz=orth.seg.ad$human$ir
zz[orth.seg.ad$human$seg$alt.size>10,] = NA
zz[orth.seg.ad$human$i+orth.seg.ad$human$e<50] = NA
pairs(cbind(zz[,c(h1,h)]))

plot(zz[,h1[1]],zz[,h[1]])
points(zz[ids,h1[1]],zz[ids,h[1]],col='red',pch=19)
ids=rownames(zz)[!is.na(zz[,h1[1]]+zz[,h[1]]) & zz[,h[1]]<0.2 & zz[,h1[1]]>0.6]
points(zz[orth.seg.ad$human$seg$alt.size>20,h1[1]],zz[orth.seg.ad$human$seg$alt.size>20,h[1]],col='red',pch=19)

orth.seg.ad$human$seg[ids,]

par(mfcol=c(2,2),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
hist(ens.ge.cod$human$cnts[,1],xlab='gene count',main='Gene counts')
hist(log10(ens.ge.cod$human$cnts[,1]),xlab='log10(gene count)',main='Gene counts')
plot(ens.ge.cod$human$cnts[,1],ens.ge.cod$human$cnts[,2],xlab='sample 1',ylab='sample 2')
plot(ens.ge.cod$human$cnts[,1],ens.ge.cod$human$cnts[,2],xlab='sample 1',ylab='sample 2',log='xy')
