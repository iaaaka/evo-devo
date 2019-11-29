options(stringsAsFactors = FALSE)
source('code/r.functions/load.all.data.F.R')
library(SAJR)

species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
orth.seg.ad.tsm = readRDS('Rdata/orth.seg.ad.tsm.Rdata')
seg2ens = readRDS('Rdata/seg2ens.Rdata')
orth.ens.genes = readRDS('Rdata/orth.ens.genes.Rdata')
psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
anns = readRDS('Rdata/anns.Rdata')
all.anns = readRDS('Rdata/all.anns.Rdata')
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')
ens.ge.marg.tsm = readRDS('Rdata/ens.ge.marg.tsm.Rdata')

plotTissueAgeProile(ens.ge.marg.tsm$human['ENSG00000155657',],meta.tsm)

orth.ens.genes['ENSG00000155657',]
titin.sids = names(seg2ens$human)[sapply(seg2ens$human,function(x)'ENSG00000155657' %in% x)]
titin.sids.alt = intersect(titin.sids,rownames(anns$human))
titin.sids.orth = intersect(titin.sids,rownames(orth.seg.ad.tsm$human))
anns$human[titin.sids.alt,]
anns$human[titin.sids.alt,'alt.id']


par(mfrow=c(5,6),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,2.5),oma=c(0,0,0,1))
for(id in titin.sids.alt[per.tissue.age.qv$human[titin.sids.alt,'heart']<0.05])
	plotTissueAgeProile(psi.tsm$human[id,],meta.tsm,main=anns$human[id,'alt.id'],age.axis = 'rank')

anns$human['hum.37457.s408',]
intersect(titin.sids.alt[per.tissue.age.qv$human[titin.sids.alt,'heart']<0.05],rownames(orth.seg.ad.tsm$human))
par(mfrow=c(3,3),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,2.5),oma=c(0,0,0,1))
for(s in names(orth.seg.ad.tsm))
	plotTissueAgeProile(orth.seg.ad.tsm[[s]][10409,],meta.tsm,age.axis = 'rank',main=s)


hh1 = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/',meta$fname[meta$species=='human' & meta$tissue=='heart' & meta$stage %in% c('4wpc')],'.bam')
hh2 = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/',meta$fname[meta$species=='human' & meta$tissue=='heart' & meta$stage %in% c('newborn')],'.bam')
hh3 = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/',meta$fname[meta$species=='human' & meta$tissue=='heart' & meta$stage %in% c('youngadult','youngteenager','infant')],'.bam')
id='hum.37457.s408'
id='hum.37457.s315'
getReadCoverage(hh1,anns$human[id,'chr_id'],anns$human[id,'start']-3000,anns$human[id,'stop']+10000,1,T,min.junc.cov = 5,reverse=anns$human[id,'strand']==-1,ylab='coverage',main='Human adult brain',xlab='position (nt)')
getReadCoverage(hh3,anns$human[id,'chr_id'],179390716,179695529,1,T,min.junc.cov = 5,reverse=anns$human[id,'strand']==-1,ylab='coverage',main='Human adult brain',xlab='position (nt)',ylim=c(0,600))

#
h = readRDS('Rdata/human.as.u.all.Rdata')
c38 = read.table('input/TTN/ttn_exons_new_human38.bed')
c37 = read.table('input/TTN/ttn_exons_new_human37.bed')
c37 = paste0(c37[,1],':',c37[,2],'-',c37[,3])
c38 = paste0(c38[,1],':',c38[,2],'-',c38[,3])

coor2sid = setNames(rownames(h$seg),paste0(h$seg$chr_id,':',h$seg$start,'-',h$seg$stop))
f = c37 %in% names(coor2sid)
table(f)
ids = rownames(meta)[meta$species=='human' & meta$tissue=='heart']
psi= h$ir[coor2sid[c37[f]],ids]
colnames(psi) = meta[colnames(psi),'marg.name']
psi = psi[,order(as.numeric(sapply(strsplit(colnames(psi),'.',T),'[',3)))]

count= h$i[coor2sid[c37[f]],ids]
colnames(count) = meta[colnames(count),'marg.name']
count = count[,order(as.numeric(sapply(strsplit(colnames(count),'.',T),'[',3)))]


r = data.frame(coor=c38[f],
							 my.id=rownames(h$seg[coor2sid[c37[f]],]),
							 type = h$seg[coor2sid[c37[f]],'type'])
r$padj = NA
f1 = r$my.id %in% rownames(per.tissue.age.qv$human)
r$padj[f1] =	per.tissue.age.qv$human[r$my.id[f1],'heart']

r = cbind(r,psi)
write.table(r,file='input/TTN/psi.tab',sep='\t',quote = F,row.names = F)
write.table(count,file='input/TTN/counts.tab',sep='\t',quote = F,row.names = F)

all.coor2sid = setNames(rownames(all.anns$human),paste0(all.anns$human$chr_id,':',all.anns$human$start,'-',all.anns$human$stop))
table(c37 %in% names(all.coor2sid))
cmn = intersect(c37,names(all.coor2sid))
table(all.anns$human[all.coor2sid[cmn],'type'])

rc = read.table('processed/mapping/hisat2.f/human/mapping.stat',row.names = 1)
colnames(rc) = c('read.count','mapped.uniq','mapped.mult')
rc = rc[meta$lib.id[meta$species=='human' & meta$tissue=='heart'],]
rownames(rc) = meta$marg.name[meta$species=='human' & meta$tissue=='heart']
rc = rc[colnames(psi),]
write.table(rc,sep='\t',quote = F,file='input/TTN/lib.sizes.tab')
