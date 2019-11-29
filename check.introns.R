options(stringsAsFactors = FALSE)
library(SAJR)
setwd('~/skoltech/projects/evo.devo/')
source('code/r.functions/check.introns.F.R')

#load('Rdata/intron.stat.Rdata')
species = c(m='mouse',r='rat',b='rabbit',h='human',q='macaque',o='opossum',c='chicken')

#########################
#### load intron stat ###
for(s in species) {
	print(s)
	assign(substr(s,1,3),loadJunctionStat(paste('processed/mapping/hisat2.f/',s,'/',sep='') ,
																				paste('processed/mapping/junctions/',s,'.merged.gff',sep='')))
}
#save(mou,rat,rab,hum,mac,opo,chi,file = 'Rdata/intron.stat.Rdata')


pdf('figures/intron.filtering/hisat.f.intron.stat.pdf',w=12,h=21)
par(mfrow=c(7,2),mar=c(4,4,3,4),mgp=c(2.3,0.8,0))
plotintronStatForSpecies(mou,'Mouse')
plotintronStatForSpecies(rat,'Rat')
plotintronStatForSpecies(rab,'Rabbit')
plotintronStatForSpecies(hum,'Human')
plotintronStatForSpecies(mac,'Macaque')
plotintronStatForSpecies(opo,'Opossum')
plotintronStatForSpecies(chi,'Chicken')
dev.off()

boxplot(apply(hum$seq.stat,1,sum),apply(mou$seq.stat,1,sum),apply(rat$seq.stat,1,sum),apply(rab$seq.stat,1,sum),apply(chi$seq.stat,1,sum),names = c('hum','mou','rat','rab','chi'))


all.int.stat = vector('list',length(species))
names(all.int.stat) =  substr(species,1,3)
for(s in species) {
	print(s)
	short = substr(s,1,3)
	all.int.stat[[short]] = loadMergedLiftovered(paste('processed/mapping/junctions/',s,'.all.sp.merged.gz',sep=''),species)
	s = all.int.stat[[short]]$introns$seq
	r = all.int.stat[[short]]$introns$strand
	all.int.stat[[short]]$introns$canonical = (s %in% c('GTAG', 'GCAG', 'ATAC') & r != '-') | (s %in% c('CTAC', 'CTGC', 'GTAT') & r != '+')
	#saveRDS(all.int.stat[[short]],paste('Rdata/',short,'.intron.spMerged.stat.Rdata',sep=''))
}



pdf('figures/intron.filtering/hisat.f.spMerged.intron.stat.pdf',w=12,h=21)
for(s in species){
	x = substr(s,1,3)
	print(x)
	self = all.int.stat[[x]]$sam.cnts[,s] > 0
	sp.cnt = apply(all.int.stat[[x]]$sam.cnts>0,1,sum)
	sam.cnt = apply(all.int.stat[[x]]$sam.cnts,1,sum)
	seqs = factor(all.int.stat[[x]]$introns$seq)
	
	par(mfrow=c(7,2),mar=c(4,4,3,4),mgp=c(2.3,0.8,0),oma=c(0,0,1,0))
	for(i in 1:7){
		f = self & sp.cnt==i
		t = table(sam.cnt[f],seqs[f])
		
		plotSeqStat(t,xlab='# samples',main=paste('Present in ',x,' and ',(i-1),' other species.' ,sep=''),ylab='freq',xlog=TRUE,yrange = c(0.5,1.5),ylim=c(0,1.5),plot.total = T)
		plotSeqStat(t,ord=c('ATAC','GCAG','CTGC','GTAT'),xlog=TRUE,add=T,yrange = c(0,0.45))
		
		f = !self & sp.cnt==i
		t = table(sam.cnt[f],seqs[f])
		
		if(i<7){
			plotSeqStat(t,xlab='# samples',main=paste('Absent in ',x,' but observed in ',i,' other species.' ,sep=''),ylab='freq',xlog=TRUE,yrange = c(0.5,1.5),ylim=c(0,1.5),plot.total = T)
			plotSeqStat(t,ord=c('ATAC','GCAG','CTGC','GTAT'),xlog=TRUE,add=T,yrange = c(0,0.45))
		}
	}
	mtext(x,3,outer = TRUE)
}
dev.off()


pdf('figures/intron.filtering/intron.sample.cnt.by.species.violin.pdf',w=12,h=30)
par(mfrow=c(7,1),mar=c(5,3,2,3),mgp=c(2.3,0.8,0),oma=c(0,0,1,0))
for(x in substr(species,1,3)){
	print(x)
	plotIntronSampleCntBySpecies(all.int.stat[[x]],xlab='Species',ylab='# of samples',main=x)
}
dev.off()


pdf('figures/intron.filtering/intron.canonical.prop.by.species.pdf',w=8,h=28)
par(mfrow=c(7,2),mar=c(5,3,2,3),mgp=c(2.3,0.8,0),oma=c(0,0,1,0))
for(x in species){
	print(x)
	plotProportionOfCannonicalSites(all.int.stat[[substr(x,1,3)]],x,4,'==',main=x)
	plotProportionOfCannonicalSites(all.int.stat[[substr(x,1,3)]],x,4,'>=',main=x)
}
dev.off()

r = NULL
for(x in species){
	print(x)
	z = all.int.stat[[substr(x,1,3)]]
	self = z$sam.cnts[,x] > 3
	sam.cnt = apply(z$sam.cnts>3,1,sum)
	
	r = rbind(r,
				c(self.canonical = sum(self & z$introns$canonical),
				self.noncanonical = sum(self & !z$introns$canonical),
				nonself.oneSp.canonical = sum(!self & sam.cnt == 1 & z$introns$canonical),
				nonself.moreSp.canonical = sum(!self & sam.cnt > 1 & z$introns$canonical),
				nonself.moreSp.noncanonical = sum(!self & sam.cnt > 1 & !z$introns$canonical) ))
	write.table(z$introns[self | (sam.cnt == 1 & z$introns$canonical) | sam.cnt > 1 ,1:4],file=paste('processed/mapping/junctions/',x,'.spMerged.filtered.splicesites',sep=''),sep='\t',quote=F,col.names=F,row.names=F)
}
rr = r
rownames(r) = species

pdf('figures/intron.filtering/intron.counts.pdf',w=17,h=8)
par(mfrow=c(1,2),mar=c(5,5,2,3),mgp=c(2.3,0.8,0),oma=c(0,0,1,0))
barplot(t(r),las=2,col=c('red','red','blue','green','green'),density = c(-1,10,-1,-1,10),legend.text = T,xlim=c(0,9),ylim=c(0,1.3e6),main='# introns')
rr = sweep(r,1,apply(r,1,sum),'/')
barplot(t(rr),las=2,col=c('red','red','blue','green','green'),density = c(-1,10,-1,-1,10),legend.text = T,xlim=c(0,9),ylim=c(0,1.3),main='# introns')
dev.off()
