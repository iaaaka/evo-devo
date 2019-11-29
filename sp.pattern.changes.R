options(stringsAsFactors = FALSE)
source('code/r.functions/load.all.data.F.R')
library(SAJR)

species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
#orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
orth.seg.ad.tsm = readRDS('Rdata/orth.seg.ad.tsm.Rdata')
orth.per.tissue.age.qv = readRDS('Rdata/orth.per.tissue.age.qv.Rdata')
age.al.i = read.table('input/margarida.age.alignment.final.2017.05.17.improved.tab',sep='\t',header = T)[,c(rownames(species),'to.remove')]
seg2ens = readRDS('Rdata/seg2ens.Rdata')
orth.ens.genes = readRDS('Rdata/orth.ens.genes.Rdata')
ens.ge.marg.tsm = readRDS('Rdata/ens.ge.marg.tsm.Rdata')
exon.birth.one = readRDS('Rdata/exon.birth.one.Rdata')
born.exn.sajr = readRDS('Rdata/born.exn.sajr.Rdata')
all.anns = readRDS('Rdata/all.anns.Rdata')

#hmo
hmo.age.al = age.al.i[-c(3,5),c(1,3,6)]
apply(hmo.age.al,2,function(x)length(unique(x)))

getAgeAlPSI = function(p,s,t,as){
	p = p[[s]]
	cn = paste(s,t,as)
	r = matrix(NA,ncol=length(cn),nrow=nrow(p),dimnames = list(rownames(p),cn))
	cn = intersect(cn,colnames(p))
	r[,cn] = p[,cn]
	r
}

hb = getAgeAlPSI(orth.seg.ad.tsm,'human','brain',hmo.age.al$human)
mb = getAgeAlPSI(orth.seg.ad.tsm,'mouse','brain',hmo.age.al$mouse)
ob = getAgeAlPSI(orth.seg.ad.tsm,'opossum','brain',hmo.age.al$opossum)

c = sapply(1:nrow(hb),function(i)cor(hb[i,],mb[i,]))
c = sapply(1:nrow(hb),function(i){mean(abs(hb[i,]-mb[i,]-mean(hb[i,])+mean(mb[i,])))})

table(!is.na(c),sgn=orth.per.tissue.age.qv$human[,'brain']<0.05 & orth.per.tissue.age.qv$mouse[,'brain']<0.05)
sgn = (orth.per.tissue.age.qv$human[,'brain']<0.05) + (orth.per.tissue.age.qv$mouse[,'brain']<0.05)
boxplot(c ~ sgn)
hist(c[orth.per.tissue.age.qv$human[,'brain']<0.05 & orth.per.tissue.age.qv$mouse[,'brain']<0.05])
table(c< -0.5,sgn)
which(c< -0.5 & (orth.per.tissue.age.qv$human[,'brain']<0.05 & orth.per.tissue.age.qv$mouse[,'brain']<0.05))

which(sgn==2 & c > 0.2)

plotTissueAgeProile(orth.seg.ad.tsm$human[11261,],meta.tsm,age.axis = 'rank',tissues = 'brain',ylim=c(0,1))
#plotTissueAgeProile(orth.seg.ad.tsm$macaque[11261,],meta.tsm,age.axis = 'rank',tissues = 'brain',ylim=c(0,1))
plotTissueAgeProile(orth.seg.ad.tsm$mouse[11261,],meta.tsm,age.axis = 'rank',tissues = 'brain',add=T,col='red')
plotTissueAgeProile(orth.seg.ad.tsm$opossum[11261,],meta.tsm,tissues = 'brain',add=T,col='orange',age.axis = 'rank')

which(orth.seg.ad$human$seg$alt.size<5 & abs(mb[,1]-mb[,12])<0.1 & abs(ob[,1]-ob[,12])<0.1 & abs(hb[,1]-hb[,12]) > 0.5 & orth.per.tissue.age.qv$human[,'brain']<0.05 & orth.per.tissue.age.qv$mouse[,'brain']>0.05)
meta$lib.id[meta$species=='human' & meta$tissue=='brain' & meta$stage=='4wpc']

#11261 - coverage plot is bad, seems to be gene merge
#7582 - coverage plot is bad, seems to be gene merge
h1 = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/',meta$fname[meta$species=='human' & meta$tissue=='brain' & meta$stage=='4wpc'],'.bam')
h2 = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/',meta$fname[meta$species=='human' & meta$tissue=='brain' & meta$stage=='youngmidage'],'.bam')
m1 = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/mouse/',meta$fname[meta$species=='mouse' & meta$tissue=='brain' & meta$stage=='10.5'],'.bam')
m2 = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/mouse/',meta$fname[meta$species=='mouse' & meta$tissue=='brain' & meta$stage=='9wpb'],'.bam')
o1 = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/opossum/',meta$fname[meta$species=='opossum' & meta$tissue=='brain' & meta$stage=='13.5'],'.bam')
o2 = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/opossum/',meta$fname[meta$species=='opossum' & meta$tissue=='brain' & meta$stage=='120dpb'],'.bam')


hh2 = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/',meta$fname[meta$species=='human' & meta$tissue=='heart' & meta$stage %in% c('toddler','youngmidage','youngadult')],'.bam')
mh2 = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/mouse/',meta$fname[meta$species=='mouse' & meta$tissue=='heart' & meta$stage %in% c('2wpb','4wpb','9wpb')],'.bam')
oh2 = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/opossum/',meta$fname[meta$species=='opossum' & meta$tissue=='heart' & meta$stage %in% c('60dpb','90dpb','120dpb')],'.bam')

# AMPD2 adenosine monophosphate deaminase 2 ENSG00000116337 inx = 40186
pdf('figures/paper.figures/3.2/examples/AMPD2.adenosine.monophosphate.deaminase.ENSG00000116337.prim-spec.pdf',w=12,h=9)
inx = 40186
layout(matrix(1:9,ncol=3),widths = c(1,1,2))
par(tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,2.5),oma=c(0,0,0,1))
plotTissueAgeProile(orth.seg.ad.tsm$human[inx,],meta.tsm,age.axis = 'rank',main='human',ylim=c(0,1),ylab='PSI')
#plotTissueAgeProile(orth.seg.ad.tsm$macaque[inx,],meta.tsm,age.axis = 'rank',ylim=c(0,1))
plotTissueAgeProile(orth.seg.ad.tsm$mouse[inx,],meta.tsm,age.axis = 'rank',main='mouse',ylim=c(0,1),ylab='PSI')
plotTissueAgeProile(orth.seg.ad.tsm$opossum[inx,],meta.tsm,age.axis = 'rank',main='opossum',ylim=c(0,1),ylab='PSI')

eids = sapply(colnames(hmo.age.al),function(s)seg2ens[[s]][[rownames(orth.seg.ad[[s]]$seg)[inx]]])
plotTissueAgeProile(ens.ge.marg.tsm$human[eids[1],],meta.tsm,age.axis = 'rank',main='human',ylab='RPKM')
plotTissueAgeProile(ens.ge.marg.tsm$mouse[eids[2],],meta.tsm,age.axis = 'rank',main='mouse',ylab='RPKM')
plotTissueAgeProile(ens.ge.marg.tsm$opossum[eids[3],],meta.tsm,age.axis = 'rank',main='opossum',ylab='RPKM')


getReadCoverage(h2,orth.seg.ad$human$seg$chr_id[inx],orth.seg.ad$human$seg$start[inx]-4500,orth.seg.ad$human$seg$stop[inx]+500,NA,T,min.junc.cov = 5,reverse=orth.seg.ad$human$seg$strand[inx]==-1,ylab='coverage',main='Human adult brain',xlab='position (nt)')
abline(v=orth.seg.ad$human$seg$start[inx]/2+orth.seg.ad$human$seg$stop[inx]/2,col='blue',lty=3)
getReadCoverage(m2,orth.seg.ad$mouse$seg$chr_id[inx],orth.seg.ad$mouse$seg$start[inx]-500,orth.seg.ad$mouse$seg$stop[inx]+4500,NA,T,min.junc.cov = 5,reverse=orth.seg.ad$mouse$seg$strand[inx]==-1,ylab='coverage',main='Mouse adult brain',xlab='position (nt)')
abline(v=orth.seg.ad$mouse$seg$start[inx]/2+orth.seg.ad$mouse$seg$stop[inx]/2,col='blue',lty=3)
getReadCoverage(o2,orth.seg.ad$opossum$seg$chr_id[inx],orth.seg.ad$opossum$seg$start[inx]-6000,orth.seg.ad$opossum$seg$stop[inx]+500,NA,T,min.junc.cov = 5,reverse=orth.seg.ad$opossum$seg$strand[inx]==-1,ylab='coverage',main='Opossum adult brain',xlab='position (nt)')
abline(v=orth.seg.ad$opossum$seg$start[inx]/2+orth.seg.ad$opossum$seg$stop[inx]/2,col='blue',lty=3)
dev.off()

# R3HDM1 R3H domain containing 1 ENSG00000048991 inx = 43683
pdf('figures/paper.figures/3.2/examples/R3HDM1.R3H.domain.containing.ENSG00000048991prim-spec.pdf',w=12,h=9)
inx = 43683
layout(matrix(1:9,ncol=3),widths = c(1,1,2))
par(tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,2.5),oma=c(0,0,0,1))
plotTissueAgeProile(orth.seg.ad.tsm$human[inx,],meta.tsm,age.axis = 'rank',main='human',ylim=c(0,1),ylab='PSI')
#plotTissueAgeProile(orth.seg.ad.tsm$macaque[inx,],meta.tsm,age.axis = 'rank',ylim=c(0,1))
plotTissueAgeProile(orth.seg.ad.tsm$mouse[inx,],meta.tsm,age.axis = 'rank',main='mouse',ylim=c(0,1),ylab='PSI')
plotTissueAgeProile(orth.seg.ad.tsm$opossum[inx,],meta.tsm,age.axis = 'rank',main='opossum',ylim=c(0,1),ylab='PSI')

eids = sapply(colnames(hmo.age.al),function(s)seg2ens[[s]][[rownames(orth.seg.ad[[s]]$seg)[inx]]])
plotTissueAgeProile(ens.ge.marg.tsm$human[eids[1],],meta.tsm,age.axis = 'rank',main='human',ylab='RPKM')
plotTissueAgeProile(ens.ge.marg.tsm$mouse[eids[2],],meta.tsm,age.axis = 'rank',main='mouse',ylab='RPKM')
plotTissueAgeProile(ens.ge.marg.tsm$opossum[eids[3],],meta.tsm,age.axis = 'rank',main='opossum',ylab='RPKM')


getReadCoverage(h1,orth.seg.ad$human$seg$chr_id[inx],orth.seg.ad$human$seg$start[inx]-2600,orth.seg.ad$human$seg$stop[inx]+4500,NA,T,min.junc.cov = 15,reverse=orth.seg.ad$human$seg$strand[inx]==-1,ylab='coverage',main='Human embryo brain',xlab='position (nt)')
abline(v=100+orth.seg.ad$human$seg$start[inx]/2+orth.seg.ad$human$seg$stop[inx]/2,col='blue',lty=3)
getReadCoverage(m1,orth.seg.ad$mouse$seg$chr_id[inx],orth.seg.ad$mouse$seg$start[inx]-2300,orth.seg.ad$mouse$seg$stop[inx]+2500,NA,T,min.junc.cov = 15,reverse=orth.seg.ad$mouse$seg$strand[inx]==-1,ylab='coverage',main='Mouse mbryo brain',xlab='position (nt)')
abline(v=orth.seg.ad$mouse$seg$start[inx]/2+orth.seg.ad$mouse$seg$stop[inx]/2,col='blue',lty=3)
getReadCoverage(o1,orth.seg.ad$opossum$seg$chr_id[inx],orth.seg.ad$opossum$seg$start[inx]-2600,orth.seg.ad$opossum$seg$stop[inx]+2500,NA,T,min.junc.cov = 15,reverse=orth.seg.ad$opossum$seg$strand[inx]==-1,ylab='coverage',main='Opossum mbryo brain heart',xlab='position (nt)')
abline(v=orth.seg.ad$opossum$seg$start[inx]/2+orth.seg.ad$opossum$seg$stop[inx]/2,col='blue',lty=3)
dev.off()


sapply(colnames(hmo.age.al),function(s)seg2ens[[s]][[rownames(orth.seg.ad[[s]]$seg)[inx]]])
s='human'
orth.ens.genes[seg2ens[[s]][[rownames(orth.seg.ad[[s]]$seg)[inx]]],colnames(hmo.age.al)]
orth.seg.ad$human$seg[inx,]
orth.seg.ad$mouse$seg[inx,]
orth.seg.ad$opossum$seg[inx,]


orth.seg.ad$opossum$seg[7582,]
seg2ens$human[['hum.4811.s94']]
seg2ens$mouse[['mou.29343.s61']]
seg2ens$opossum[['opo.19719.s16']]
orth.ens.genes['ENSG00000178104',]

# APP ENSG00000142192

app.sids = names(seg2ens$human)[sapply(seg2ens$human,function(x)'ENSG00000142192' %in% x)]
which(sapply(exon.birth.one,function(x)x$seg_id[1] %in% app.sids))
exon.birth.one[[inx]]

pdf('figures/paper.figures/3.2/examples/APP.ENSG00000142192.placental-spec.7th.new.exon.pdf',w=12,h=9)
inx = 860
layout(matrix(1:9,ncol=3),widths = c(1,1,2))
par(tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,2.5),oma=c(0,0,0,1))
plotTissueAgeProile(born.exn.sajr$human$ir[exon.birth.one[[inx]]$seg_id[1],],meta,age.axis = 'rank',main='human',ylim=c(0,1),ylab='PSI')
#plotTissueAgeProile(orth.seg.ad.tsm$macaque[inx,],meta.tsm,age.axis = 'rank',ylim=c(0,1))
plotTissueAgeProile(born.exn.sajr$mouse$ir[exon.birth.one[[inx]]$seg_id[3],],meta,age.axis = 'rank',main='mouse',ylim=c(0,1),ylab='PSI')
plot.new();title(main='opossum')

eids = sapply(colnames(hmo.age.al),function(s)seg2ens[[s]][[exon.birth.one[[inx]][s,'useg_id']]])
orth.ens.genes[eids[1],]
plotTissueAgeProile(ens.ge.marg.tsm$human[eids[1],],meta.tsm,age.axis = 'rank',main='human',ylab='RPKM')
plotTissueAgeProile(ens.ge.marg.tsm$mouse[eids[2],],meta.tsm,age.axis = 'rank',main='mouse',ylab='RPKM')
plotTissueAgeProile(ens.ge.marg.tsm$opossum[eids[3],],meta.tsm,age.axis = 'rank',main='opossum',ylab='RPKM')

sids=exon.birth.one[[inx]][colnames(hmo.age.al),c('useg_id','seg_id','dseg_id')]

getReadCoverage(hh2,all.anns$human[sids[1,1],'chr_id'],all.anns$human[sids[1,1],'stop'],all.anns$human[sids[1,3],'start'],NA,T,min.junc.cov = 15,reverse=all.anns$human[sids[1,1],'strand']==-1,ylab='coverage',main='Human adult heart',xlab='position (nt)')
abline(v=all.anns$human[sids[1,2],'start']/2+all.anns$human[sids[1,2],'stop']/2,col='blue',lty=3)
getReadCoverage(mh2,all.anns$mouse[sids[2,1],'chr_id'],all.anns$mouse[sids[2,1],'stop'],all.anns$mouse[sids[2,3],'start'],NA,T,min.junc.cov = 15,reverse=all.anns$mouse[sids[2,1],'strand']==-1,ylab='coverage',main='Mouse adult heart',xlab='position (nt)')
abline(v=all.anns$mouse[sids[2,2],'start']/2+all.anns$mouse[sids[2,2],'stop']/2,col='blue',lty=3)
getReadCoverage(oh2,all.anns$opossum[sids[3,1],'chr_id'],all.anns$opossum[sids[3,1],'start'],all.anns$opossum[sids[3,3],'stop'],NA,T,min.junc.cov = 15,reverse=all.anns$opossum[sids[3,1],'strand']==-1,ylab='coverage',main='Opossum adult heart',xlab='position (nt)')
abline(v=all.anns$opossum[sids[3,2],'start']/2+all.anns$opossum[sids[3,2],'stop']/2,col='blue',lty=3)
dev.off()

all.anns$human[as.character(sids[1,]),]
all.anns$mouse[as.character(sids[2,]),]
all.anns$opossum[as.character(sids[3,]),]
