library(png)
library(SAJR)
source('code/r.functions/paper.figures.5.F.R')
source('~/skoltech/r.code/util.R')
source('code/r.functions/load.all.data.F.R')

library(extrafont)

species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
age.al.i = read.table('input/margarida.age.alignment.final.2017.05.17.improved.tab',sep='\t',header = T)[,c(rownames(species),'to.remove')]
orth.seg.ad.all.tsm = readRDS('Rdata/orth.seg.ad.all.tsm.Rdata')
orth.seg.hqm.ad.all.tsm = readRDS('Rdata/orth.seg.hqm.ad.tsm.Rdata')
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')
anns = readRDS('Rdata/anns.Rdata')
age.dpsi = readRDS('Rdata/age.diam.spline4.with.replicates.Rdata')
#orth.per.tissue.age.qv = readRDS('Rdata/orth.per.tissue.age.qv.Rdata')
#all.anns = readRDS('Rdata/all.anns.Rdata')


# load subsampling 
submeta = readRDS('Rdata/human-mouse2macaque.subsample.meta.Rdata')
orth.seg.ad.tsm.s = readRDS('Rdata/hqr.subsample/orth.seg.ad.tsm.Rdata')
per.tissue.age.qv.s = readRDS('Rdata/hqr.subsample/per.tissue.age.qv.Rdata')
age.dpsi.s = readRDS('Rdata/hqr.subsample/age.diam.spline4.with.replicates.Rdata')
anns.s = readRDS('Rdata/hqr.subsample/anns.Rdata')
#all.anns.s = readRDS('Rdata/hqr.subsample/all.anns.Rdata')
#orth.seg.ad.s = readRDS('Rdata/hqr.subsample/orth.seg.ad.Rdata')


params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)

# subsample ######
# sstages = lapply(rownames(species),function(s){
# 	t = meta[meta$species==s,]
# 	sort(rank(sapply(split(t$days,t$stage),mean)))
# })
# names(sstages) = rownames(species)
# 
# r = NULL
# set.seed(24892)
# sps = c('human','mouse','macaque')
# for(s in which(age.al.i$macaque!=''))
# 	for(t in unique(meta$tissue)){
# 		sams = list()
# 		for(sp in sps){
# 			if(age.al.i[s,sp]=='') next
# 			f = meta$species==sp & meta$tissue==t & sstages[[sp]][meta$stage] <= sstages[[sp]][age.al.i[s,sp]]
# 			if(age.al.i[s-1,sp]!='')
# 				f = f & sstages[[sp]][meta$stage] > sstages[[sp]][age.al.i[s-1,sp]]
# 			sams[[sp]] = rownames(meta)[f]
# 		}
# 		sam.no = min(sapply(sams[sps],length))
# 		if(sam.no > 0)
# 			for(sp in sps){
# 				sp.sams = sample(sams[[sp]],sam.no)
# 				r = rbind(r,data.frame(mouse.stage=age.al.i$mouse[s],species=sp,tissue=t,sid=sp.sams))
# 			}
# 	}
# r = cbind(meta[r$sid,],mouse.stage=r$mouse.stage)
# table(r$species,r$tissue,r$mouse.stage)
# saveRDS(r,'Rdata/human-mouse2macaque.subsample.meta.Rdata')


#take m and age.al.i_ in sample heatmap section
sps = c('human','macaque','mouse')
pdf('figures/paper.figures/6/nature.review/04.subsampling/samples.pdf',w=5,h=6,family='Arial')
plotSampleTable1(submeta,sps,cex=0.8)
dev.off()


# make commands for stringtie
# t = split(submeta,paste(submeta$species,submeta$tissue,submeta$mouse.stage,sep='_'))
# fnames = sapply(t,function(x)paste0('../../mapping/hisat2.s/',x$species,'/',x$fname,'.bam',collapse = ' '))
# # sp = sapply(strsplit(names(t),'_'),'[',1) should be done in this way
# # but actualy I used sp='macaque'
# sam.cnt = sapply(t,nrow)
# tt = paste(ifelse(sam.cnt==1,'cat ','samtools merge - '),fnames," | stringtie - --bam -f 0.1 -p 12 -j 3 -g 10 --rf -l ",substr(names(t),1,3),".sub -o ",sp,"/",names(t),'.gtf',sep='')
# write.table(tt,quote = F,col.names = F,row.names = F,'processed/annotation/hqm.subsample/stringtie.commands')
# #write.table(tt,quote = F,col.names = F,row.names = F,'~/skoltech.tmp/projects/evo.devo/processed/annotation/hqm.subsample/stringtie.commands')
# 
# r = c()
# for(s in c('human','mouse','macaque')){
# 	r = c(r,paste0("java -Xmx10g -jar $SGE_O_HOME/bin/sajr.jar count_reads sajr.config -batch_in=",
# 								 paste0('../mapping/hisat2.s/',s,'/',submeta$fname[submeta$species==s],'.bam',collapse = ','),' -ann_in=',
# 								 '../annotation/hqm.subsample/merged/',s,'.sajr -batch_out=',
# 								 paste0("hqm.subsample.uniq/",submeta$fname[submeta$species==s],collapse = ',')," -use_mult=false >> hqm.subsample.uniq.log"))
# }
# 
# writeLines(r,'~/skoltech.tmp/projects/evo.devo/processed/sajr/hqm.subsample.uniq.run.sh')
# load SAJR data
# f = function(a,fs,nm,ens){
# 	r = loadSAData(a,fs,nm)
# 	r = setSplSiteTypes(r,a)
# 	r = addIsCogingByEnsGTF(ens,r)
# 	r$ir[r$i + r$e < 10] = NA
# 	r$seg$length = r$seg$stop - r$seg$start + 1
# 	r$seg = getDistanceToClosestSegs(r,0.5)
# 	r$seg = addExonNumber(r$seg)
# 	r$seg = addCDSPos2Ann(r$seg)
# 	r
# }
# 
# ens = c(mouse='Mus_musculus.GRCm38.84',
# 				human='Homo_sapiens.GRCh37.73',
# 				macaque='Macaca_mulatta.MMUL_1.84')
# 
# 
# tmp = lapply(names(ens),function(s){
# 	print(s)
# 	m = submeta[submeta$species==s,]
# 	t = f(paste('processed/annotation/hqm.subsample/merged/',s,'.sajr',sep=''),paste('processed/sajr/hqm.subsample.uniq//',m$fname,sep=''),m$name,paste('processed/annotation/all.species/ensambl/',ens[s],".gtf.gz",sep=''))
# 	saveRDS(t,paste('Rdata/hqr.subsample/',s,'.as.u.all.Rdata',sep=''))
# 	print('NA freq:')
# 	print(table(is.na(t$seg$gene_id)))
# 	gc()
# 	t$seg
# })
# 
# all.anns = tmp
# names(all.anns) = names(ens)
# all.anns = all.anns[c('human','macaque','mouse')]
#saveRDS(all.anns,'Rdata/hqr.subsample/all.anns.Rdata')

# tmp = do.call(rbind,lapply(all.anns,function(x)cbind(rownames(x),x$sites)))
# write.table(tmp,'processed/annotation/hqm.subsample/merged/seg2sites.tab',quote = F,row.names = F,col.names = F,sep='\t')


# orth.seg.ad = loadAltOrthSegs(c('processed/orth.segs/only.ad/hqm.subsample.hqm.0.6.orth.ad.segs','processed/orth.segs/only.ad/hqm.subsample.hqm.0.0.orth.ad.segs'),only.filtered = FALSE)
# o.by.neig = findOrthExonsByNeighbors(all.anns.s,orth.seg.ad,'human')
# t = table(o.by.neig$same.frame,floor(o.by.neig$length.min/o.by.neig$length.max*20)/20)
#  
# par(mfrow=c(1,2),tck=-0.02,mgp=c(1.3,0.4,0),mar=c(3,3,1.5,0),oma=c(0,0,2,1))
# barplot(t,xlab='min(length)/max(length)',ylab='# of orth exons',legend.text=c('different frame','same frame'),args.legend=list(x='topleft'))
# barplot(sweep(t,2,apply(t,2,sum),'/'),xlab='min(length)/max(length)',ylab='proportion of orth exons',)
# abline(h=2/3,col='red')
# 
# 
# orth.seg.ad = rbind(orth.seg.ad,o.by.neig[o.by.neig$length.min == o.by.neig$length.max,colnames(orth.seg.ad)])
# dim(orth.seg.ad)
# 
# loadInfoForOrths. = function(o){
# 	info = vector('list',ncol(o))
# 	names(info) = colnames(o)
# 	for(s in colnames(o)){
# 		print(s)
# 		info[[s]] = readRDS(paste('Rdata/hqr.subsample/',s,'.as.u.all.Rdata',sep=''))[o[!is.na(o[,s]),s],]
# 		gc(verbose = FALSE)
# 	}
# 	info
# }
# 
# orth.seg.ad.s = loadInfoForOrths.(orth.seg.ad)
#saveRDS(orth.seg.ad.s,'Rdata/hqr.subsample/orth.seg.ad.Rdata')
orth.seg.ad.ids.s = sapply(orth.seg.ad.s,function(x)rownames(x$seg))
#saveRDS(orth.seg.ad.ids.s,'Rdata/hqr.subsample/orth.seg.ad.ids.Rdata')

# load all-data hqr
# orth.seg.ad = loadAltOrthSegs(c('processed/orth.segs/only.ad/hqm.0.6.orth.ad.segs','processed/orth.segs/only.ad/hqm.0.0.orth.ad.segs'),only.filtered = FALSE)
# o.by.neig = findOrthExonsByNeighbors(all.anns[colnames(orth.seg.ad)],orth.seg.ad,'human')
# t = table(o.by.neig$same.frame,floor(o.by.neig$length.min/o.by.neig$length.max*20)/20)
# 
# par(mfrow=c(1,2),tck=-0.02,mgp=c(1.3,0.4,0),mar=c(3,3,1.5,0),oma=c(0,0,2,1))
# barplot(t,xlab='min(length)/max(length)',ylab='# of orth exons',legend.text=c('different frame','same frame'),args.legend=list(x='topleft'))
# barplot(sweep(t,2,apply(t,2,sum),'/'),xlab='min(length)/max(length)',ylab='proportion of orth exons',)
# abline(h=2/3,col='red')
# 
# 
# orth.seg.ad = rbind(orth.seg.ad,o.by.neig[o.by.neig$length.min == o.by.neig$length.max,colnames(orth.seg.ad)])
# dim(orth.seg.ad)
# 
# orth.seg.ad = loadInfoForOrths(orth.seg.ad)
# saveRDS(orth.seg.ad,'Rdata/orth.seg.hqm.ad.Rdata')
# 
# orth.seg.hqm.ad.ids = sapply(orth.seg.ad,function(x)rownames(x$seg))
# saveRDS(orth.seg.hqm.ad.ids,'Rdata/orth.seg.hqm.ad.ids.Rdata')

# orth.seg.hqm.ad.tsm = lapply(orth.seg.ad,function(x){
# 		m = meta[colnames(x$ir),]
# 		calcMeanCols(x$ir,paste(m$species,m$tissue,m$stage))
# 	})
# saveRDS(orth.seg.hqm.ad.tsm,'Rdata/orth.seg.hqm.ad.tsm.Rdata')

alt.sp.s =getAltSp(orth.seg.ad.tsm.s,0.9 ,4)
alt.sp = getAltSp(orth.seg.ad.all.tsm,0.9 ,4)
alt.sp.hqm = getAltSp(orth.seg.hqm.ad.all.tsm,0.9 ,4)

table(alt.sp.s)

alts.s = sapply(orth.seg.ad.s,function(x)x$seg$type)
f.s = apply(alts.s=='ALT',1,sum)>0

alts = sapply(orth.seg.ad,function(x)x$seg$type)
f = apply(alts=='ALT',1,sum)>0
table(f)



f = function(alt,ylim=NULL,frac=FALSE,...){
	x=sapply(c('h','q','m'),function(x)sum(grepl(x,alt)))
	ci=sapply(x,function(n)binom.test(n,length(alt))$conf)*100
	if(frac)
		x = x/length(alt)*100
	else
		ci = ci*length(alt)/100
	if(is.null(ylim))
		ylim = range(0,ci)
	b=barplot(x,ylab=paste0(ifelse(frac,'%','#'),' of AS'),ylim=ylim,...)
	arrows(b,ci[1,],b,ci[2,],angle = 90,code=3,length=0.03)
}

pdf('figures/paper.figures/6/nature.review/04.subsampling/numbers.of.as.in.orth.pdf',w=6,h=7,family='Arial')
par(mfrow=c(2,3),tck=-0.01,mgp=c(2,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,0,1))
f(alt.sp    ,ylim=c(9,24000),main='All data')
f(alt.sp.hqm,ylim=c(9,24000),main='All data, HQM only')
f(alt.sp.s  ,ylim=c(9,24000),main='Subsampling')
f(alt.sp    ,ylim=c(0,18),frac=T,main='All data')
f(alt.sp.hqm,ylim=c(0,18),frac=T,main='All data, HQM only')
f(alt.sp.s  ,ylim=c(0,18),frac=T,main='Subsampling')
dev.off()

par(mfcol=c(2,3))
barplot(sapply(c('h','q','m'),function(x)sum(grepl(x,alt.sp.s)))/length(alt.sp.s))
barplot(sapply(c('h','q','m'),function(x)sum(grepl(x,alt.sp)))/length(alt.sp))


barplot(sapply(c('h','q','m'),function(x)sum(grepl(x,alt.sp.s[f.s])))/sum(f.s))
barplot(sapply(c('h','q','m'),function(x)sum(grepl(x,alt.sp[f])))/sum(f))

f = apply(alts[,c('human','macaque','mouse')]=='ALT',1,sum)>0
table(f)
barplot(sapply(c('h','q','m'),function(x)sum(grepl(x,alt.sp.s[f.s])))/sum(f.s))
barplot(sapply(c('h','q','m'),function(x)sum(grepl(x,alt.sp[f])))/sum(f))

# filter
# anns.s = list()
# for(s in names(orth.seg.ad.s)){
# 	t = readRDS(paste('Rdata/hqr.subsample/',s,'.as.u.all.Rdata',sep=''))
# 	st = length(t)
# 	t = t[t$seg$chr_id !='MT' & t$seg$position=='INTERNAL' & t$seg$type!='EXN',]
# 	st[2] = length(t)
# 	f = rep(FALSE,length(t))
# 	na = is.na(t$ir)
# 	for(tis in unique(submeta$tissue)){
# 		cinx = meta[colnames(na),'tissue']==tis
# 		f = f | (apply(!na[,cinx],1,mean) > 0.6 & apply(t$ir[,cinx],1,function(x){x=x[!is.na(x)];sum(x>0.1 & x < 0.9)>3}))
# 	}
# 	t = t[f,]
# 	st[3] = length(t)
# 	cat(s,st,"\t")
# 	saveRDS(t,paste('Rdata/hqr.subsample/',s,'.as.u.filtered.Rdata',sep=''))
# 	anns.s[[s]] = t$seg
# }
# saveRDS(anns.s,'Rdata/hqr.subsample/anns.Rdata')

# run devAS (in mouse scale)
# per.tissue.age.qv.s = list()
# registerDoMC(3)
# smeta = submeta
# z = unique(smeta[smeta$species=='mouse',c('age.use','mouse.stage')])
# z = setNames(z$age.use,z$mouse.stage)
# smeta$age.use = z[smeta$mouse.stage]
# table(is.na(z[smeta$mouse.stage]))
# for(s in names(orth.seg.ad.s)){
# 	cat(toupper(s))
# 	tmp = readRDS(paste('Rdata/hqr.subsample/',s,'.as.u.filtered.Rdata',sep=''))
# 	per.tissue.age.qv.s[[s]] = sapply(unique(smeta$tissue),function(t){testASAge(tmp,smeta,t,min.cov.sams=0.6)})
# 	colnames(per.tissue.age.qv.s[[s]]) = unique(smeta$tissue)
# 	dimnames(per.tissue.age.qv.s[[s]]) = setNames(dimnames(per.tissue.age.qv.s[[s]]),NULL)
# }
# saveRDS(per.tissue.age.qv.s,'Rdata/hqr.subsample/per.tissue.age.qv.Rdata')

# age.diamsss4.s = lapply( names(orth.seg.ad.s),function(s){
# 	print(s)
# 	p = readRDS(paste0("Rdata/hqr.subsample/",s,".as.u.filtered.Rdata"))$ir
# 	m = smeta[colnames(p),]
# 	r = matrix(NA,nrow=nrow(p),ncol=length(unique(smeta$tissue)),dimnames = list(rownames(p),unique(smeta$tissue)))
# 	for(t in colnames(r)){
# 		cat('  ',t,'\t')
# 		if(sum(m$tissue==t)>0){
# 			z = p[,m$tissue==t]
# 			m. = m[colnames(z),]
# 			r[,t] = apply(z,1,function(y)getDiamBySpline(m.$age.use,y,4))
# 		}
# 	}
# 	gc()
# 	r
# })
# names(age.diamsss4.s) = names(orth.seg.ad.s)
# age.diamsss4.s = lapply(age.diamsss4.s,function(x)apply(x,1:2,min,1))
# saveRDS(age.diamsss4.s,'Rdata/hqr.subsample/age.diam.spline4.with.replicates.Rdata')

sapply(names(per.tissue.age.qv.s),function(s)apply(per.tissue.age.qv.s[[s]]<0.05 & age.dpsi.s[[s]]>0.2,2,sum,na.rm=T))

t = sapply(names(per.tissue.age.qv.s),function(s)apply(per.tissue.age.qv.s[[s]]<0.05,2,sum,na.rm=T))
barplot(t,beside = T,col=params$tissue.col[rownames(t)])


astypes = c(CE='ad',AA='aa',AD='dd',RI='da')

tested.stat = lapply(setNames(astypes,astypes),function(ss)sapply(names(per.tissue.age.qv),function(s)apply(!is.na(per.tissue.age.qv[[s]][anns[[s]]$sites==ss,]),2,sum)))

sgn02.stat = lapply(setNames(astypes,astypes),function(ss) {
	sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][anns[[s]]$sites == ss,] < 0.05 & abs(age.dpsi[[s]][anns[[s]]$sites == ss,])>0.2,2,sum,na.rm=T))/tested.stat[[ss]]*100})

sgn02.cnt=lapply(setNames(astypes,astypes),function(ss) {
	sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][anns[[s]]$sites == ss,] < 0.05 & abs(age.dpsi[[s]][anns[[s]]$sites == ss,])>0.2,2,sum,na.rm=T))})

sgn.stat = lapply(setNames(astypes,astypes),function(ss) {
	sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][anns[[s]]$sites == ss,] < 0.05,2,sum,na.rm=T))/tested.stat[[ss]]*100})

sgn.cnt=lapply(setNames(astypes,astypes),function(ss) {
	sapply(names(anns),function(s)apply(per.tissue.age.qv[[s]][anns[[s]]$sites == ss,] < 0.05,2,sum,na.rm=T))})


tested.stat.s = lapply(setNames(astypes,astypes),function(ss)sapply(names(per.tissue.age.qv.s),function(s)apply(!is.na(per.tissue.age.qv.s[[s]][anns.s[[s]]$sites==ss,]),2,sum)))
sgn02.stat.s = lapply(setNames(astypes,astypes),function(ss) {
	sapply(names(anns.s),function(s)apply(per.tissue.age.qv.s[[s]][anns.s[[s]]$sites == ss,] < 0.05 & abs(age.dpsi.s[[s]][anns.s[[s]]$sites == ss,])>0.2,2,sum,na.rm=T))/tested.stat.s[[ss]]*100})

sgn02.cnt.s=lapply(setNames(astypes,astypes),function(ss) {
	sapply(names(anns.s),function(s)apply(per.tissue.age.qv.s[[s]][anns.s[[s]]$sites == ss,] < 0.05 & abs(age.dpsi.s[[s]][anns.s[[s]]$sites == ss,])>0.2,2,sum,na.rm=T))})

sgn.stat.s = lapply(setNames(astypes,astypes),function(ss) {
	sapply(names(anns.s),function(s)apply(per.tissue.age.qv.s[[s]][anns.s[[s]]$sites == ss,] < 0.05,2,sum,na.rm=T))/tested.stat.s[[ss]]*100})

sgn.cnt.s=lapply(setNames(astypes,astypes),function(ss) {
	sapply(names(anns.s),function(s)apply(per.tissue.age.qv.s[[s]][anns.s[[s]]$sites == ss,] < 0.05,2,sum,na.rm=T))})


pdf('figures/paper.figures/6/nature.review/04.subsampling/devAS.stat.pdf',w=12,h=5,family='Arial')
par(mfrow=c(2,3),tck=-0.01,mgp=c(2,0.2,0),mar=c(4,3,1.5,0),oma=c(0,0,0,1))
astypes.pchs=c(ad=19,aa=2,dd=6,da=13)
par(mfrow=c(2,5),tck=-0.01,mgp=c(1.2,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
sps = c('human','macaque','mouse')
tiss = c('brain','kidney','liver','testis')
plotAsEventCount(tested.stat,astypes.pchs,by.tissue = F,sps=sps,tiss = tiss,ylab='# of detected events',main='Detected AS',bty='n')
plotAsEventCount(sgn02.cnt  ,astypes.pchs,by.tissue = F,sps=sps,tiss = tiss,ylab='# of devAS',main='DevAS (dPSI > 0.2)',bty='n')
plotAsEventCount(sgn02.stat ,astypes.pchs,by.tissue = F,sps=sps,tiss = tiss,ylab='% of devAS',main='DevAS (dPSI > 0.2)',bty='n')
plotAsEventCount(sgn.cnt    ,astypes.pchs,by.tissue = F,sps=sps,tiss = tiss,ylab='# of devAS',main='DevAS',bty='n')
plotAsEventCount(sgn.stat   ,astypes.pchs,by.tissue = F,sps=sps,tiss = tiss,ylab='% of devAS',main='DevAS',bty='n')

plotAsEventCount(tested.stat.s,astypes.pchs,by.tissue = F,sps=sps,tiss = tiss,ylab='# of detected events',main='Detected AS',bty='n')
plotAsEventCount(sgn02.cnt.s  ,astypes.pchs,by.tissue = F,sps=sps,tiss = tiss,ylab='# of devAS',main='DevAS (dPSI > 0.2)',bty='n')
plotAsEventCount(sgn02.stat.s ,astypes.pchs,by.tissue = F,sps=sps,tiss = tiss,ylab='% of devAS',main='DevAS (dPSI > 0.2)',bty='n')
plotAsEventCount(sgn.cnt.s    ,astypes.pchs,by.tissue = F,sps=sps,tiss = tiss,ylab='# of devAS',main='DevAS',bty='n')
plotAsEventCount(sgn.stat.s   ,astypes.pchs,by.tissue = F,sps=sps,tiss = tiss,ylab='% of devAS',main='DevAS',bty='n')


plotAsEventCount(tested.stat,astypes.pchs,by.tissue = T,sps=sps,tiss = tiss,ylab='# of detected events',main='Detected AS',bty='n')
plotAsEventCount(sgn02.cnt  ,astypes.pchs,by.tissue = T,sps=sps,tiss = tiss,ylab='# of devAS',main='DevAS (dPSI > 0.2)',bty='n')
plotAsEventCount(sgn02.stat ,astypes.pchs,by.tissue = T,sps=sps,tiss = tiss,ylab='% of devAS',main='DevAS (dPSI > 0.2)',bty='n')
plotAsEventCount(sgn.cnt    ,astypes.pchs,by.tissue = T,sps=sps,tiss = tiss,ylab='# of devAS',main='DevAS',bty='n')
plotAsEventCount(sgn.stat   ,astypes.pchs,by.tissue = T,sps=sps,tiss = tiss,ylab='% of devAS',main='DevAS',bty='n')

plotAsEventCount(tested.stat.s,astypes.pchs,by.tissue = T,sps=sps,tiss = tiss,ylab='# of detected events',main='Detected AS',bty='n')
plotAsEventCount(sgn02.cnt.s  ,astypes.pchs,by.tissue = T,sps=sps,tiss = tiss,ylab='# of devAS',main='DevAS (dPSI > 0.2)',bty='n')
plotAsEventCount(sgn02.stat.s ,astypes.pchs,by.tissue = T,sps=sps,tiss = tiss,ylab='% of devAS',main='DevAS (dPSI > 0.2)',bty='n')
plotAsEventCount(sgn.cnt.s    ,astypes.pchs,by.tissue = T,sps=sps,tiss = tiss,ylab='# of devAS',main='DevAS',bty='n')
plotAsEventCount(sgn.stat.s   ,astypes.pchs,by.tissue = T,sps=sps,tiss = tiss,ylab='% of devAS',main='DevAS',bty='n')
dev.off()


astypes.pchs=c(ad=19,aa=2,dd=6,da=13)
sps = c('human','macaque','mouse')
tiss = c('brain','kidney','liver','testis')
plotPanelLetter = function(l,cex=1.2,adj=c(0,1.1),...){
	l = tolower(l)
	x=grconvertX(0,from='nfc',to='user')
	y=grconvertY(1,from='nfc',to='user')
	text(x=x,y=y,labels=l,adj=adj,font=2,cex=cex,xpd=NA)
}
lab.cex = 1.5
params$tissue.col
pdf('figures/paper.figures/6/nature.review/04.subsampling/S22.subsampling.pdf',w=7.2,h=3,family='Arial')
par(mfrow=c(1,3),tck=-0.01,mgp=c(1.2,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
plotAsEventCount(tested.stat,astypes.pchs,by.tissue = F,sps=sps,tiss = tiss,ylab='# of detected events',main='Detected AS - whole dataset',bty='n',xlab.cex=1)
plotPanelLetter('A',lab.cex)
plotAsEventCount(tested.stat.s,astypes.pchs,by.tissue = F,sps=sps,tiss = tiss,ylab='# of detected events',main='Detected AS - subsampling',bty='n',xlab.cex=1)
plotPanelLetter('B',lab.cex)
t=rbind(tested.all  = tested.stat$ad['brain',sps],
			devAS.all   = sgn02.cnt$ad['brain',sps],
			tested.subs = tested.stat.s$ad['brain',sps],
			devASsubs   = sgn02.cnt.s$ad['brain',sps])

b=barplot(t(t),col=c('#3399CCFF','#3399CCAA','#3399CC55'),beside = T,names.arg = c('detected','devAS','detected','devAS'),legend.text = T,border=NA,args.legend = list(bty='n',border=NA),main='Brain',ylab='# of events')
text(c(mean(b[2,1:2]),mean(b[2,3:4])),grconvertY(0,'nfc','user'),c('All data','Subsampling'),xpd=TRUE,adj=c(0.5,-0.5))
plotPanelLetter('C',lab.cex)
dev.off()