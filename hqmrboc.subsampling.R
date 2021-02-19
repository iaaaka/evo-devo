library(png)
library(SAJR)
library(doMC)
library(plyr)
source('code/r.functions/paper.figures.5.F.R')
source('code/r.functions/load.all.data.F.R')
source('../../r.code/util.R')


library(extrafont)

species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
age.al.i = read.table('input/margarida.age.alignment.final.2017.05.17.improved.tab',sep='\t',header = T)[,c(rownames(species),'to.remove')]
all.anns = readRDS('Rdata/all.anns.Rdata')
orth.seg.ad.all.id = readRDS('Rdata/orth.seg.ad.all.id.Rdata')

orth.seg.ad.all.tsm = readRDS('Rdata/orth.seg.ad.all.tsm.Rdata')
orth.seg.hqm.ad.all.tsm = readRDS('Rdata/orth.seg.hqm.ad.tsm.Rdata')
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')
anns = readRDS('Rdata/anns.Rdata')
age.dpsi = readRDS('Rdata/age.diam.spline4.with.replicates.Rdata')
age.segs = readRDS('Rdata/devAS.4patt.Rdata')
#orth.per.tissue.age.qv = readRDS('Rdata/orth.per.tissue.age.qv.Rdata')
#all.anns = readRDS('Rdata/all.anns.Rdata')


# load subsampling 
submeta = readRDS('Rdata/hqmrboc.subsample.meta.Rdata')
# orth.seg.ad.tsm.s = readRDS('Rdata/hqmrboc.subsample/orth.seg.ad.tsm.Rdata')
per.tissue.age.qv.s = readRDS('Rdata/hqmrboc.subsample/per.tissue.age.qv.Rdata')
age.dpsi.s = readRDS('Rdata/hqmrboc.subsample/age.diam.spline4.with.replicates.Rdata')
anns.s = readRDS('Rdata/hqmrboc.subsample/anns.Rdata')
all.anns.s = readRDS('Rdata/hqmrboc.subsample/all.anns.Rdata')
orth.seg.ad.s = readRDS('Rdata/hqmrboc.subsample/orth.seg.ad.Rdata')

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
# sps = rownames(species)
# for(s in which(age.al.i$macaque!=''))
# 	for(t in unique(meta$tissue)){
# 		sams = list()
# 		for(sp in sps){
# 			if(age.al.i[s,sp]=='') next
# 			f = meta$species==sp & meta$tissue==t & sstages[[sp]][meta$stage] <= sstages[[sp]][age.al.i[s,sp]]
# 			if(age.al.i[s-1,sp]!='' & age.al.i[s-1,sp] != age.al.i[s,sp])
# 				f = f & sstages[[sp]][meta$stage] > sstages[[sp]][age.al.i[s-1,sp]]
# 			sams[[sp]] = setdiff(rownames(meta)[f],r$sid)
# 		}
# 		sam.no = min(sapply(sams[sps],length))
# 		if(sam.no > 0)
# 			for(sp in sps){
# 				sp.sams = sample(sams[[sp]],sam.no)
# 				r = rbind(r,data.frame(mouse.stage=age.al.i$mouse[s],species=sp,tissue=t,sid=sp.sams))
# 			}
# 	}
# r = cbind(meta[r$sid,],mouse.stage=r$mouse.stage)
# z=table(r$tissue,r$mouse.stage,r$species)
# saveRDS(r,'Rdata/hqmrboc.subsample.meta.Rdata')


#take m and age.al.i_ in sample heatmap section
# sps = rownames(species)[c(1,2,4:7,3)]
# pdf('figures/paper.figures/6/nature.review/04.subsampling/hqrmboc.samples.pdf',w=5,h=12,family='Arial')
# plotSampleTable1(submeta,sps,cex=0.8)
# dev.off()


# stringtie-SAJR#####
# t = split(submeta,paste(submeta$species,submeta$tissue,submeta$mouse.stage,sep='_'))
# fnames = sapply(t,function(x)paste0('../../mapping/hisat2.s/',x$species,'/',x$fname,'.bam',collapse = ' '))
# sp = sapply(strsplit(names(t),'_'),'[',1)
# sam.cnt = sapply(t,nrow)
# tt = paste(ifelse(sam.cnt==1,'cat ','samtools merge - '),fnames," | stringtie - --bam -f 0.1 -p 12 -j 3 -g 10 --rf -l ",substr(names(t),1,3),".sub -o ",sp,"/",names(t),'.gtf',sep='')
# write.table(tt,quote = F,col.names = F,row.names = F,'processed/annotation/hqmrboc.subsample/stringtie.commands')

# r = c()
# for(s in rownames(species)){
# 	r = c(r,paste0("java -Xmx10g -jar $SGE_O_HOME/bin/sajr.jar count_reads sajr.config -batch_in=",
# 								 paste0('../mapping/hisat2.s/',s,'/',submeta$fname[submeta$species==s],'.bam',collapse = ','),' -ann_in=',
# 								 '../annotation/hqmrboc.subsample/merged/',s,'.sajr -batch_out=',
# 								 paste0("hqmrboc.subsample.uniq/",submeta$fname[submeta$species==s],collapse = ',')," -use_mult=false >> hqmrboc.subsample.uniq.log"))
# }
# 
# writeLines(r,'~/skoltech/projects/evo.devo/processed/sajr/hqmrboc.subsample.uniq.run.sh')

# load SAJR data ######
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
# 				rat='Rattus_norvegicus.Rnor_5.0.79',
# 				rabbit='Oryctolagus_cuniculus.OryCun2.0.84',
# 				human='Homo_sapiens.GRCh37.73',
# 				macaque='Macaca_mulatta.MMUL_1.84',
# 				opossum='Monodelphis_domestica.BROADO5.84',
# 				chicken='Gallus_gallus.Galgal4.84')
# 
# tmp = lapply(names(ens),function(s){
# 	print(s)
# 	m = submeta[submeta$species==s,]
# 	t = f(paste('processed/annotation/hqmrboc.subsample/merged/',s,'.sajr',sep=''),paste('processed/sajr/hqmrboc.subsample.uniq/',m$fname,sep=''),m$name,paste('processed/annotation/all.species/ensambl/',ens[s],".gtf.gz",sep=''))
# 	saveRDS(t,paste('Rdata/hqmrboc.subsample/',s,'.as.u.all.Rdata',sep=''))
# 	#t = readRDS(paste('Rdata/hqmrboc.subsample/',s,'.as.u.all.Rdata',sep=''))
# 	print('NA freq:')
# 	print(table(is.na(t$seg$gene_id)))
# 	gc()
# 	t$seg
# })
# 
# all.anns = tmp
# names(all.anns) = names(ens)
# all.anns = all.anns[rownames(species)]
# saveRDS(all.anns,'Rdata/hqmrboc.subsample/all.anns.Rdata')

 # tmp = do.call(rbind,lapply(all.anns,function(x)cbind(rownames(x),x$sites)))
 # write.table(tmp,'processed/annotation/hqmrboc.subsample/merged/seg2sites.tab',quote = F,row.names = F,col.names = F,sep='\t')
# 
# load orth ######
# orth.seg.ad = loadAltOrthSegs(c('processed/annotation/hqmrboc.subsample/merged/orth.seg/hqmrboc.0.6.orth.ad.segs','processed/annotation/hqmrboc.subsample/merged/orth.seg/hqmrboc.0.0.orth.ad.segs'),only.filtered = FALSE)
# o.by.neig = findOrthExonsByNeighbors(all.anns.s,orth.seg.ad,'human')
# t = table(o.by.neig$same.frame,floor(o.by.neig$length.min/o.by.neig$length.max*20)/20)
# 
# par(mfrow=c(1,2),tck=-0.02,mgp=c(1.3,0.4,0),mar=c(3,3,1.5,0),oma=c(0,0,2,1))
# barplot(t,xlab='min(length)/max(length)',ylab='# of orth exons',legend.text=c('different frame','same frame'),args.legend=list(x='topleft'))
# barplot(sweep(t,2,apply(t,2,sum),'/'),xlab='min(length)/max(length)',ylab='proportion of orth exons',)
# abline(h=2/3,col='red')
# 
# orth.seg.ad = rbind(orth.seg.ad,o.by.neig[o.by.neig$length.min == o.by.neig$length.max,colnames(orth.seg.ad)])
# dim(orth.seg.ad)
# 
# loadInfoForOrths. = function(o){
# 	info = vector('list',ncol(o))
# 	names(info) = colnames(o)
# 	for(s in colnames(o)){
# 		print(s)
# 		info[[s]] = readRDS(paste('Rdata/hqmrboc.subsample/',s,'.as.u.all.Rdata',sep=''))[o[!is.na(o[,s]),s],]
# 		gc(verbose = FALSE)
# 	}
# 	info
# }

# orth.seg.ad.s = loadInfoForOrths.(orth.seg.ad)
# saveRDS(orth.seg.ad.s,'Rdata/hqmrboc.subsample//orth.seg.ad.Rdata')
# orth.seg.ad.ids.s = sapply(orth.seg.ad.s,function(x)rownames(x$seg))
# saveRDS(orth.seg.ad.ids.s,'Rdata/hqmrboc.subsample/orth.seg.ad.ids.Rdata')


# orth.seg.hqm.ad.tsm = lapply(orth.seg.ad.s,function(x){
# 		m = meta[colnames(x$ir),]
# 		calcMeanCols(x$ir,paste(m$species,m$tissue,m$stage))
# 	})
# saveRDS(orth.seg.hqm.ad.tsm,'Rdata/hqmrboc.subsample/orth.seg.ad.tsm.Rdata')

# filter #####
# anns.s = list()
# for(s in names(orth.seg.ad.s)){
# 	t = readRDS(paste('Rdata/hqmrboc.subsample/',s,'.as.u.all.Rdata',sep=''))
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
# 	saveRDS(t,paste('Rdata/hqmrboc.subsample/',s,'.as.u.filtered.Rdata',sep=''))
# 	anns.s[[s]] = t$seg
# }
# saveRDS(anns.s,'Rdata/hqmrboc.subsample/anns.Rdata')

# run devAS ######
# per.tissue.age.qv.s = list()
# registerDoMC(3)
# smeta = submeta
# z = unique(smeta[smeta$species=='mouse',c('age.use','mouse.stage')])
# z = setNames(z$age.use,z$mouse.stage)
# smeta$age.use = z[smeta$mouse.stage]
# table(is.na(z[smeta$mouse.stage]))
# for(s in names(orth.seg.ad.s)){
# 	cat(toupper(s))
# 	tmp = readRDS(paste('Rdata/hqmrboc.subsample/',s,'.as.u.filtered.Rdata',sep=''))
# 	per.tissue.age.qv.s[[s]] = sapply(unique(smeta$tissue),function(t){testASAge(tmp,smeta,t,min.cov.sams=0.6)})
# 	colnames(per.tissue.age.qv.s[[s]]) = unique(smeta$tissue)
# 	dimnames(per.tissue.age.qv.s[[s]]) = setNames(dimnames(per.tissue.age.qv.s[[s]]),NULL)
# }
# saveRDS(per.tissue.age.qv.s,'Rdata/hqmrboc.subsample/per.tissue.age.qv.Rdata')
# 
# age.diamsss4.s = lapply(rownames(species),function(s){
# 	print(s)
# 	p = readRDS(paste0("Rdata/hqmrboc.subsample/",s,".as.u.filtered.Rdata"))$ir
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
# saveRDS(age.diamsss4.s,'Rdata/hqmrboc.subsample/age.diam.spline4.with.replicates.Rdata')


sapply(names(per.tissue.age.qv.s),function(s)apply(per.tissue.age.qv.s[[s]]<0.05 & age.dpsi.s[[s]]>0.2,2,sum,na.rm=T))
sapply(names(per.tissue.age.qv.s),function(s)apply(per.tissue.age.qv.s[[s]]<0.05,2,sum,na.rm=T))
sapply(names(per.tissue.age.qv.s),function(s)apply(!is.na(per.tissue.age.qv.s[[s]]),2,sum,na.rm=T))


# detected AS ######
detAS = lapply(per.tissue.age.qv.s,function(x){x[!is.na(x)]='n';x[is.na(x)]='-';x})
ts = getSegTestDevAsStat(detAS,anns.s)
os = getSegTestDevAsStat(detAS,anns.s,sapply(orth.seg.ad.s,function(x)rownames(x$seg)))

t = getSegTestDevAsStat(age.segs,anns)
o = getSegTestDevAsStat(age.segs,anns,sapply(orth.seg.ad.all.tsm,rownames))

astypes.pchs=c(ad=19,aa=2,dd=6,da=13)
tiss = c('brain','kidney','liver','testis')

pdf('figures/paper.figures/6/nature.review/04.subsampling/hqmrboc.subsample.detAS.stat.pdf',w=12,h=6,family='Arial')
par(mfrow=c(2,4),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(1.5,2.5,1.5,0),oma=c(0,0,0,0))
plotAsEventCount(t$tested,astypes.pchs[1],tiss = tiss,by.tissue = FALSE,ylab='# of detected events',main='All: Detected AS',bty='n',xlab.cex = 0.6)
plotPanelLetter('A',lab.cex)
plotAsEventCount(o$tested,astypes.pchs[1],tiss = tiss,by.tissue = FALSE,ylab='# of detected events',main='All orth: Detected AS',bty='n',xlab.cex = 0.6)
plotPanelLetter('B',lab.cex)
plotAsEventCount(t$tested,astypes.pchs[1],tiss = tiss,by.tissue = TRUE,ylab='# of detected events',main='All: Detected AS',bty='n')
plotPanelLetter('A',lab.cex)
plotAsEventCount(o$tested,astypes.pchs[1],tiss = tiss,by.tissue = TRUE,ylab='# of detected events',main='All orth: Detected AS',bty='n')
plotPanelLetter('B',lab.cex)


plotAsEventCount(ts$tested,astypes.pchs[1],tiss = tiss,by.tissue = FALSE,ylab='# of detected events',main='Subsample: Detected AS',bty='n',xlab.cex = 0.6)
plotPanelLetter('A',lab.cex)
plotAsEventCount(os$tested,astypes.pchs[1],tiss = tiss,by.tissue = FALSE,ylab='# of detected events',main='Subsample orth: Detected AS',bty='n',xlab.cex = 0.6)
plotPanelLetter('B',lab.cex)
plotAsEventCount(ts$tested,astypes.pchs[1],tiss = tiss,by.tissue = TRUE,ylab='# of detected events',main='Subsample: Detected AS',bty='n')
plotPanelLetter('A',lab.cex)
plotAsEventCount(os$tested,astypes.pchs[1],tiss = tiss,by.tissue = TRUE,ylab='# of detected events',main='Subsample orth: Detected AS',bty='n')
plotPanelLetter('B',lab.cex)
dev.off()


# whait is wrong with Macauqe? ####
submeta. = readRDS('Rdata/human-mouse2macaque.subsample.meta.Rdata')
m = readRDS(paste('Rdata/macaque.as.u.filtered.Rdata',sep=''))
mm = meta[colnames(m$ir),]

bmds = cmdscale(1-cor(m$ir[,mm$tissue=='brain'],u='p'),2)
f2 = rownames(bmds) %in% rownames(submeta)
f1 = rownames(bmds) %in% rownames(submeta.)

plot(bmds,pch=19,cex=mm$cex*2,bty='n')
points(bmds[f1,],pch=19,cex=mm$cex[f1]*2,bty='n',col='red')
points(bmds[f2,],pch=19,cex=mm$cex[f2]*2,bty='n',col='blue')


par(mfrow=c(2,3),mar=c(5,5,1,0))
gem = readRDS('Rdata/ens.ge.marg.Rdata')$macaque
mm1 = meta[colnames(m$ir),]
mm1 = mm1[rownames(mm1) %in% rownames(submeta.) & mm1$tissue=='brain',]
bmds1 = cmdscale(1-cor(m$ir[,rownames(mm1)],u='p'),2)
plot(bmds1,pch=19,cex=mm1$cex*2,bty='n')
imageWithText(cor(m$ir[,rownames(mm1)],u='p'),breaks=seq(0.8,1,length.out = 100),col=rev(heat.colors(99)),names.as.labs = T)
imageWithText(cor(gem[,rownames(mm1)],u='p',m='sp'),breaks=seq(0.8,1,length.out = 100),col=rev(heat.colors(99)),names.as.labs = T)
imageWithText(cor(log(gem[,rownames(mm1)]+0.1),u='p',m='p'),breaks=seq(0.9,1,length.out = 100),col=rev(heat.colors(99)),names.as.labs = T)


mm2 = meta[colnames(m$ir),]
mm2 = mm2[rownames(mm2) %in% rownames(submeta) & mm2$tissue=='brain',]
bmds2 = cmdscale(1-cor(m$ir[,rownames(mm2)],u='p'),2)
plot(bmds2,pch=19,cex=mm2$cex*2,bty='n')
imageWithText(cor(m$ir[,rownames(mm2)],u='p'),breaks=seq(0.8,1,length.out = 100),col=rev(heat.colors(99)),names.as.labs = T)
imageWithText(cor(gem[,rownames(mm2)],u='p',m='sp'),breaks=seq(0.8,1,length.out = 100),col=rev(heat.colors(99)),names.as.labs = T)
imageWithText(cor(log(gem[,rownames(mm2)]+0.1),u='p',m='p'),breaks=seq(0.9,1,length.out = 100),col=rev(heat.colors(99)),names.as.labs = T)


smeta = submeta
z = unique(smeta[smeta$species=='mouse',c('age.use','mouse.stage')])
z = setNames(z$age.use,z$mouse.stage)
smeta$age.use = z[smeta$mouse.stage]
table(is.na(smeta$age.use))

smeta. = submeta.
z = unique(smeta.[smeta.$species=='mouse',c('age.use','mouse.stage')])
z = setNames(z$age.use,z$mouse.stage)
smeta.$age.use = z[smeta.$mouse.stage]
table(is.na(smeta.$age.use))

f = smeta$tissue=='brain' & smeta$species=='macaque'
plot(log(smeta$days[f]),smeta$age.use[f])

f = smeta.$tissue=='brain' & smeta.$species=='macaque'
plot(log(smeta.$days[f]),smeta.$age.use[f])
smeta.[f,]

qall = readRDS(paste('Rdata/macaque.as.u.filtered.Rdata',sep=''))
qsu1 = readRDS(paste('Rdata/hqr.subsample/macaque.as.u.filtered.Rdata',sep=''))
qsu2 = readRDS(paste('Rdata/hqmrboc.subsample/macaque.as.u.filtered.Rdata',sep=''))

qall$seg$coor = paste(qall$seg$chr_id,qall$seg$strand,qall$seg$start,qall$seg$stop)
qsu1$seg$coor = paste(qsu1$seg$chr_id,qsu1$seg$strand,qsu1$seg$start,qsu1$seg$stop)
qsu2$seg$coor = paste(qsu2$seg$chr_id,qsu2$seg$strand,qsu2$seg$start,qsu2$seg$stop)

registerDoMC(3)
f = qall$seg$coor %in% qsu1$seg$coor
qbpv.all  = testASAge(qall[f,],meta  ,'brain',min.cov.sams=0.6,return.pv = TRUE)
qbpv.su1  = testASAge(qall[f,],smeta.,'brain',min.cov.sams=0.6,return.pv = TRUE)
qbpv.su2  = testASAge(qall[f,],smeta ,'brain',min.cov.sams=0.6,return.pv = TRUE)
qbpv.su1.oa=testASAge(qall[f,], meta[rownames(smeta.),] ,'brain',min.cov.sams=0.6,return.pv = TRUE)
qbpv.su2.oa=testASAge(qall[f,], meta[rownames(smeta),] ,'brain',min.cov.sams=0.6,return.pv = TRUE)



qv = cbind(apply(apply(qbpv.all[,-1],2,p.adjust,m='BH'),1,min),
					 apply(apply(qbpv.su1[,-1],2,p.adjust,m='BH'),1,min),
					 apply(apply(qbpv.su2[,-1],2,p.adjust,m='BH'),1,min))
hist(apply(apply(qbpv.su2.oa[,-1],2,p.adjust,m='BH'),1,min),)

hist(qv[,3])
hist(qbpv.sub2.)

qbpv.sub3 = testASAge(tmp,smeta[rownames(smeta) !='QB530M_3',],'brain',min.cov.sams=0.6)
hist(qbpv.sub1[c0 %in% c1])
hist(qbpv.sub1.)

qbpv.sub1. = testASAge(tmp1,smeta.,'brain',min.cov.sams=0.9)
qbpv.sub2. = testASAge(tmp2,smeta,'brain',min.cov.sams=0.9)
qbpv.sub3. = testASAge(tmp2,smeta[rownames(smeta) !='QB530M_3',],'brain',min.cov.sams=0.6)


qbpv.sub1.pv = testASAge(tmp1,smeta.,'brain',min.cov.sams=0.6,return.pv = T)
qbpv.sub2.pv = testASAge(tmp2,smeta ,'brain',min.cov.sams=0.6,return.pv = T)

qbpv.sub2.pv[1:2,]

hist(qbpv.sub2.pv[,2],ylim=c(0,7000))
hist(log(qbpv.sub1.pv[,1]),1000,xlim=c(-5,5))
wilcox.test(qbpv.sub1.pv[,1],qbpv.sub2.pv[,1],alternative = 'g')
# replace QB530M_3 with QB530M_2
meta['QB530M_2',]
z = loadSAData('processed/annotation/hqmrboc.subsample/merged/macaque.sajr',fnames = 'processed/sajr/hqmrboc.subsample.uniq/add/6068sTS.Macaque.Brain.1ypb.Male',lib_names = 'QB530M_2')
z$ir[z$i+z$e < 10] = NA
table(rownames(tmp2$seg) %in% rownames(z$seg))
z = z[rownames(tmp2$seg),]
tmp3 = tmp2[,colnames(tmp2$ir) != 'QB530M_3']
tmp3$ir = cbind(tmp3$ir,QB530M_2=z$ir)
tmp3$i = cbind(tmp3$i,QB530M_2=z$i)
tmp3$e = cbind(tmp3$e,QB530M_2=z$e)
smeta3 = smeta
smeta3['QB530M_3',]
rownames(smeta3)[rownames(smeta3)=='QB530M_3'] = 'QB530M_2'
table(colnames(tmp3$e) %in% rownames(smeta3))
qbpv.sub3. = testASAge(tmp3,smeta3,'brain',min.cov.sams=0.6)

hist(qbpv.sub3.)
m3 = smeta3[colnames(tmp3$ir),]
imageWithText(cor(tmp3$ir[,substr(colnames(tmp3$ir),2,2)=='B'],u='p'),breaks=seq(0.8,1,length.out = 100),col=rev(heat.colors(99)),names.as.labs = T)
imageWithText(cor(tmp1$ir[,substr(colnames(tmp1$ir),2,2)=='B'],u='p'),breaks=seq(0.8,1,length.out = 100),col=rev(heat.colors(99)),names.as.labs = T)



hist(qbpv.all)
hist(qbpv.sub1.)
hist(qbpv.sub2.)
hist(qbpv.sub3.)
meta['QB530M_3',]
submeta[submeta$species=='macaque' & submeta$tissue=='brain' & submeta$mouse.stage=='9wpb',]
submeta.[submeta.$species=='macaque' & submeta.$tissue=='brain' & submeta.$mouse.stage=='9wpb',]


plot(qbpv.sub1,qbpv.sub2,pch='.')
abline(a=0,b=1,col='red')

table(c2 %in% c1)
cmn = intersect(c1,c2)
inx1 = match(cmn,c1)
inx2 = match(cmn,c2)



imageWithText(cor(tmp1$ir[inx1,substr(colnames(tmp1$ir),2,2)=='B'],
									tmp3$ir[inx2,substr(colnames(tmp3$ir),2,2)=='B'],u='p'),breaks=seq(0.8,1,length.out = 100),col=rev(heat.colors(99)),names.as.labs = T)
imageWithText(cor(tmp1$ir[inx1,substr(colnames(tmp1$ir),2,2)=='B'],u='p'),breaks=seq(0.8,1,length.out = 100),col=rev(heat.colors(99)),names.as.labs = T)
imageWithText(cor(tmp3$ir[inx2,substr(colnames(tmp3$ir),2,2)=='B'],u='p'),breaks=seq(0.8,1,length.out = 100),col=rev(heat.colors(99)),names.as.labs = T)

z = intersect(colnames(tmp1$ir),colnames(tmp3$ir))
z = z[substr(z,2,2)=='B']

imageWithText(cor(cbind(tmp1$ir[inx1,z],tmp3$ir[inx2,z]),u='p'),breaks=seq(0.75,1,length.out = 100),col=rev(heat.colors(99)),names.as.labs = T)
ppp = cbind(tmp1$ir[inx1,c('QB93F_1','QB530M_2')],QB530M_2_=tmp3$ir[inx2,'QB530M_2'])
pppi = cbind(tmp1$i[inx1,c('QB93F_1','QB530M_2')],QB530M_2_=tmp3$i[inx2,'QB530M_2'])
pppe = cbind(tmp1$e[inx1,c('QB93F_1','QB530M_2')],QB530M_2_=tmp3$e[inx2,'QB530M_2'])

f = apply(is.na(ppp),1,sum)==0
cor(ppp,u='p')
apply(is.na(ppp),2,sum)
pairs(ppp[f,],pch='.')
pairs(ppp,pch='.',log='',panel = function(x,y,...){points(x,y,...);abline(a=0,b=1,col='red')})

table(is.na(ppp[,3]),pppi[,3]+pppe[,3]>9)

plot(qbpv.sub1.[inx1],qbpv.sub2.[inx2],pch='.')#,log='xy',xlim=c(1e-10,1),ylim=c(1e-10,1))
abline(a=0,b=1,col='red')


table(qbpv.sub1.[inx1]<1e-5,qbpv.sub2.[inx2]>0.5)


match(cmn[which(qbpv.sub1.[inx1]<1e-5 & qbpv.sub2.[inx2]>0.5)],c0)
col=rep('black',nrow(meta))
col[rownames(meta) %in% rownames(submeta.)] = 'blue'
col[rownames(meta) %in% rownames(submeta)] = 'red'
col[rownames(meta) %in% rownames(submeta) & rownames(meta) %in% rownames(submeta.)] = 'green'
plotTissueAgeProile(tmp$ir[64317,],meta,col=col,pch=19,cex=1,age.axis = 'rank',tissues = 'brain',bty='n')
abline(v=1:30,lty=3)


m = meta[colnames(tmp$ir),]
m = m[m$tissue=='brain',]
mds = cmdscale(1-cor(tmp$ir[tmp$seg$sites=='ad' &tmp$seg$cod=='c',rownames(m)],u='p'),k=2)
col=rep('black',nrow(m))
col[rownames(m) %in% rownames(submeta.)] = 'blue'
col[rownames(m) %in% rownames(submeta)] = 'red'
col[rownames(m) %in% rownames(submeta) & rownames(m) %in% rownames(submeta.)] = 'green'


plot(mds,pch=19,cex=(m$cex-0.2)*2,col=col)
mds[mds[,1]>0.2,]


# bootstrap subsampling ######
# check that pv is NA mean low coverage
q = readRDS('Rdata/hqmrboc.subsample/macaque.as.u.all.Rdata')
t = q$ir[q$seg$chr_id !='MT' & q$seg$position=='INTERNAL' & q$seg$type!='EXN',substr(colnames(q$ir),2,2)=='B']
f = (apply(!is.na(t),1,mean) > 0.6 & apply(t,1,function(x){x=x[!is.na(x)];sum(x>0.1 & x < 0.9)>3}))
table(f,rownames(t) %in% rownames(anns.s$macaque))
table(f,rownames(t) %in% rownames(per.tissue.age.qv.s$macaque)[!is.na(per.tissue.age.qv.s$macaque[,'brain'])])


sstages = lapply(rownames(species),function(s){
	t = meta[meta$species==s,]
	sort(rank(sapply(split(t$days,t$stage),mean)))
})
names(sstages) = rownames(species)


r = NULL
sps = rownames(species)
for(s in which(age.al.i$macaque!=''))
	for(t in unique(meta$tissue)){
		sams = list()
		for(sp in sps){
			if(age.al.i[s,sp]=='') next
			f = meta$species==sp & meta$tissue==t & sstages[[sp]][meta$stage] <= sstages[[sp]][age.al.i[s,sp]]
			if(age.al.i[s-1,sp]!='' & age.al.i[s-1,sp] != age.al.i[s,sp])
				f = f & sstages[[sp]][meta$stage] > sstages[[sp]][age.al.i[s-1,sp]]
			sams[[sp]] = setdiff(rownames(meta)[f],r$sid)
		}
		for(sp in sps){
			if(length(sams[[sp]])>0)
				r = rbind(r,data.frame(mouse.stage=age.al.i$mouse[s],species=sp,tissue=t,sid=sams[[sp]]))
		}
	}
r = cbind(meta[r$sid,],mouse.stage=r$mouse.stage)

sps = rownames(species)[c(1,2,4:7,3)]
pdf('figures/paper.figures/6/nature.review/04.subsampling/hqrmboc.all.possible.samples.pdf',w=5,h=12,family='Arial')
plotSampleTable1(r,sps,cex=0.8)
dev.off()

bootstrapDetectedAS = function(m,d,t,replace=TRUE,N=100,ss=c('ad','aa','dd','da')){
	stop('DEPRECATED: there should be at least one tissue that pass thrs, and this tissue should have 60% samples not NA, that is bit different from what is done here. see bootstrapDevAS')
	sps = unique(m$species)
	m = m[m$tissue==t, ]
	stat = apply(table(factor(m$species,levels = sps),m$mouse.stage),2,min)
	d = d[d$seg$chr_id !='MT' & d$seg$position=='INTERNAL' & d$seg$type!='EXN',colnames(d$ir) %in% rownames(m)]
	m = m[colnames(d$ir),]
	if(nrow(m)<4 | sum(stat)<4) return(NULL)
	res = NULL
	res = llply(1:N,function(i){
		#cat('\r',i,'   ')
		sids = c()
		for(i in 1:length(stat)){
			sids = c(sids,sample(rownames(m)[m$mouse.stage==names(stat)[i]],stat[i],replace = replace))
		}
		t = d$ir[,sids,drop=F]
		f = (apply(!is.na(t),1,mean) > 0.6 & apply(t,1,function(x){x=x[!is.na(x)];sum(x>0.1 & x < 0.9)>3}))
		table(d$seg$sites[f])[ss]
	},.parallel = TRUE)
	do.call(rbind,res)
}
# 
# rev1.3: The authors should assess/discuss how incomplete coverage of developmental stages in different tissues/species affects detection of devAS and classification of temporal patterns, as well as cross-species comparison.
# rev3.1: For instance, absolute numbers of cassette exons in macaque (ext Fig 4, and dev AS in general) show a pattern strickingly different to the other species, with much less AS events except in testis.  Is it due to a different developmental sampling or genome quality or is it a biological effect
#         In the external Fig 3, we can see that there are twice more alternative cassette exons in human than in oppossum/rabbit : is it explained by genome quality?
# 
# registerDoMC(3)
# detASboot.sub.orth = detASboot.sub = detASboot.orth =detASboot = list()
# for(s in rownames(species)){
# 	gc()
# 	print(s)
# 	d = readRDS(paste0('Rdata/',s,'.as.u.all.Rdata'))
# 	detASboot[[s]] = list()
# 	for(t in unique(meta$tissue)){
# 		print(paste0('   ',t))
# 		of = rownames(d$seg) %in% rownames(orth.seg.ad.all.tsm[[s]])
# 		sf = paste(d$seg$chr_id,d$seg$strand,d$seg$start,d$seg$stop) %in% paste(all.anns.s[[s]]$chr_id,all.anns.s[[s]]$strand,all.anns.s[[s]]$start,all.anns.s[[s]]$stop)
# 		detASboot[[s]][[t]]          = bootstrapDetectedAS(r,d,t,N=100)
# 		detASboot.orth[[s]][[t]]     = bootstrapDetectedAS(r,d[of,],t,N=100)
# 		detASboot.sub[[s]][[t]]      = bootstrapDetectedAS(r,d[sf,],t,N=100)
# 		detASboot.sub.orth[[s]][[t]] = bootstrapDetectedAS(r,d[of & sf,],t,N=100)
# 	}
# }
# 
# saveRDS(detASboot,'Rdata/hqmrboc.subsample/detected.as.bootstrap.Rdata')
# saveRDS(detASboot.orth,'Rdata/hqmrboc.subsample/detected.as.bootstrap.orth.Rdata')
# saveRDS(detASboot.sub,'Rdata/hqmrboc.subsample/detected.as.bootstrap.sub.Rdata')
# saveRDS(detASboot.sub.orth,'Rdata/hqmrboc.subsample/detected.as.bootstrap.sub.orth.Rdata')

detASboot = readRDS('Rdata/hqmrboc.subsample/detected.as.bootstrap.Rdata')
detASboot.orth = readRDS('Rdata/hqmrboc.subsample/detected.as.bootstrap.orth.Rdata')
detASboot.sub = readRDS('Rdata/hqmrboc.subsample/detected.as.bootstrap.sub.Rdata')
detASboot.sub.orth = readRDS('Rdata/hqmrboc.subsample/detected.as.bootstrap.sub.orth.Rdata')
# f = function(x)lapply(x,function(y){y$ovary=NULL;y})
# detASboot = f(detASboot)
# detASboot.orth = f(detASboot.orth)
# detASboot.sub = f(detASboot.sub)
# detASboot.sub.orth = f(detASboot.sub.orth)


astypes.pchs=c(ad=19,aa=2,dd=6,da=13)




pdf('figures/paper.figures/6/nature.review/04.subsampling/hqrmboc.subsample.bootstrap100.on.whole.ann.pdf',w=6,h=6,family='Arial')
par(mfrow=c(2,2),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(1.5,2.5,1.5,0),oma=c(0,0,0,0))
plotAsEventCount.boot(detASboot,astypes.pchs[1],by.tissue = FALSE,ylab='# of detected events',main='All: Detected AS',bty='n',xlab.cex = 0.6)
plotPanelLetter('A',lab.cex)
plotAsEventCount.boot(detASboot.orth,astypes.pchs[1],by.tissue = FALSE,ylab='# of detected events',main='Orth: Detected AS',bty='n',xlab.cex = 0.6)
plotPanelLetter('B',lab.cex)
plotAsEventCount.boot(detASboot.sub,astypes.pchs[1],by.tissue = FALSE,ylab='# of detected events',main='Sub-ann: Detected AS',bty='n',xlab.cex = 0.6)
plotPanelLetter('C',lab.cex)
plotAsEventCount.boot(detASboot.sub.orth,astypes.pchs[1],by.tissue = FALSE,ylab='# of detected events',main='Orth & Sub-ann: Detected AS',bty='n',xlab.cex = 0.6)
plotPanelLetter('D',lab.cex)
dev.off()

# 20201203 ############
# run hqmrboc.subsampling.cluster.R on cluster
# age.al.i. = age.al.i
# age.al.i.$opossum[c(3,6)] = ''
# age.al.i.$chicken[c(8,10)] = ''
# sapply(age.al.i.[,1:7],function(x){t = table(x[x!='']);t[t>1]})
# 
# sstages = lapply(rownames(species),function(s){
# 	t = meta[meta$species==s,]
# 	sort(rank(sapply(split(t$days,t$stage),mean)))
# })
# names(sstages) = rownames(species)
# 
# 
# meta$mouse.stage = NA
# for(s in rownames(species)){
# 	for(t in names(sort(-sstages[[s]]))){
# 		if(t %in% age.al.i.[[s]]){
# 			ms = age.al.i.[age.al.i.[[s]] == t,'mouse']
# 			if(length(ms)!=1) stop('sww')
# 			meta$mouse.stage[meta$species==s & sstages[[s]][meta$stage] <= sstages[[s]][t]] = ms
# 		}
# 	}
# }
# 
# m = meta[meta$species=='human',]
# plotSampleTable1(meta[!is.na(meta$mouse.stage),],rownames(species),cex=0.8)
# 
# m = meta[!is.na(meta$mouse.stage),]
# registerDoMC(3)
# dh = readRDS(paste0('Rdata/human.as.u.all.Rdata'))
# dh = dh[dh$seg$chr_id !='MT' & dh$seg$position=='INTERNAL' & dh$seg$type!='EXN',colnames(dh$ir) %in% rownames(m)]
# system.time({t = bootstrapDevAS(m,dh,c('human','macaque'),N=100)})
# 
# 
# dh = readRDS(paste0('Rdata/macaque.as.u.all.Rdata'))
# dh = dh[dh$seg$chr_id !='MT' & dh$seg$position=='INTERNAL' & dh$seg$type!='EXN',colnames(dh$ir) %in% rownames(m)]
# system.time({t = bootstrapDevAS(m,dh[1:1000,],c('human','macaque'),N=10)})

#
loadDevAsBootstrap = function(data,oids=NULL,types=list(det=c('c','n','u','d','ud','du'),sgn=c('u','d','ud','du'),u='u',d='d',ud='ud',du='du'),sites=c('ad','aa','dd','da')){
	names(data) = sps
	if(!is.null(oids)){
		for(s in names(data)){
			data[[s]] = lapply(data[[s]],function(x)x[rownames(x) %in% oids[,s],])
		}
	}

	r = lapply(sites,function(site){
				stat = lapply(types,function(t){
					tmp = lapply(names(data),function(species){
						sapply(data[[species]],function(x)apply(x[all.anns[[species]][rownames(x),'sites']==site,],2,function(z)sum(z %in% t)))
					})
					names(tmp) = names(data)
					tmp
				})
				stat$sgn.frac = lapply(names(stat$sgn),function(species){d=stat$sgn[[species]]/stat$det[[species]]})
				names(stat$sgn.frac) = names(stat$sgn)
				lapply(stat,function(x)lapply(x,function(z)apply(z,1,quantile,prob=c(0.025,0.5,0.975),na.rm=TRUE)))
		})
	names(r) = sites
	lapply(r,function(z){
		lapply(z,function(x){
			list(median=sapply(x,function(y)y[2,]),
				 cil=sapply(x,function(y)y[1,]),
				 cih=sapply(x,function(y)y[3,]))
			})
	})
}


sp2test = c(list(c('human','macaque','mouse','rat'),c('human','macaque'),c('human','macaque','chicken'),c('human','chicken'),c('human','mouse','rat','rabbit','opossum'),rownames(species)),as.list(rownames(species)))
devas.boot = list()
for(sps in sp2test){
	n = paste0(species[sps,'short'],collapse = '')
	print(n)

	data = lapply(paste0('Rdata/hqmrboc.subsample/chicken.12-0dpb.20201225/devAS.subsampling-',paste0(sps,collapse = '.'),'-',sps,'.Rdata'),readRDS)
	names(data) = sps

	devas.boot[[n]] = list()

	devas.boot[[n]]$all = loadDevAsBootstrap(data)
	devas.boot[[n]]$orth = loadDevAsBootstrap(data,orth.seg.ad.all.id)
	gc()
}

devas.boot = devas.boot[order(nchar(names(devas.boot)),decreasing = T)]
#saveRDS(devas.boot,'Rdata/hqmrboc.subsample/devAS.subsampling.summary.chicken.12-0dpb.20201225.Rdata')
#saveRDS(devas.boot,'Rdata/hqmrboc.subsample/devAS.subsampling.summary.chicken.12-0dpb.Rdata')
#saveRDS(devas.boot,'Rdata/hqmrboc.subsample/devAS.subsampling.summary.Rdata')
t = readRDS('Rdata/hqmrboc.subsample/detected.as.all.direct.Rdata')
t['mouse',,'ad']
devas.boot$m$all$ad$det$median


r = devas.boot$c$all
plotAsEventCount(r$ad$det$median,cil=r$ad$det$cil,cih=r$ad$det$cih)
plotAsEventCount(r$ad$sgn$median,cil=r$ad$sgn$cil,cih=r$ad$sgn$cih)

# these two are wrong!
# devas.boot. = readRDS('Rdata/hqmrboc.subsample/devAS.subsampling.summary.chicken.12-0dpb.Rdata')
# devas.boot = readRDS('Rdata/hqmrboc.subsample/devAS.subsampling.summary.Rdata')
# devas.boot[names(devas.boot.)] = devas.boot.
devas.boot = readRDS('Rdata/hqmrboc.subsample/devAS.subsampling.summary.chicken.12-0dpb.20201225.Rdata') #it is the only correct one !

pdf('figures/paper.figures/6/nature.review/04.subsampling/devAS.bootstrap.04.pdf',w=6*3,h=6*3,family='Arial')
par(mfrow=c(6,6),tck=-0.01,mgp=c(1.2,0.3,0),mar=c(3,2.5,1.5,1),oma=c(0,0,0,1),bty='n')
for(n in names(devas.boot)){
	for(m in names(devas.boot[[n]])){
		r = devas.boot[[n]][[m]]
		plotAsEventCount(r$ad$det$median,cil=r$ad$det$cil,cih=r$ad$det$cih,main=paste0(first2Upper(m),': detAS'),ylab='# of detected CE')
		plotAsEventCount(r$ad$sgn$median,cil=r$ad$sgn$cil,cih=r$ad$sgn$cih,main=paste0(first2Upper(m),': devAS'),ylab='# of devAS')
		plotAsEventCount(r$ad$sgn.frac$median*100,cil=r$ad$sgn.frac$cil*100,cih=r$ad$sgn.frac$cih*100,main=paste0(first2Upper(m),': devAS'),ylab='% of devAS')
	}
}
dev.off()

t = readRDS('Rdata/hqmrboc.subsample/devAS.subsampling-human.macaque.mouse.rat.rabbit.opossum.chicken-mouse.Rdata')
table(t[[1]][,'heart'])


devas.boot. = readRDS('Rdata/hqmrboc.subsample/devAS.subsampling.summary.Rdata')
plot(devas.boot$hq$all$ad$det$median[,1],devas.boot.$hq$all$ad$det$median[,1])
abline(a=0,b=1)

r = devas.boot.$hqc$all$ad$det
#r = lapply(r,function(x)x[,c('human','macaque')])
plotAsEventCount(r$median,cil=r$cil,cih=r$cih)


# #####################
m = meta[!is.na(meta$mouse.stage) & meta$species=='mouse',]
m$sam.id = rownames(m)
s=ss='mouse'
dh = readRDS(paste0('Rdata/',s,'.as.u.all.Rdata'))
dh = dh[dh$seg$chr_id !='MT' & dh$seg$position=='INTERNAL' & dh$seg$type!='EXN',colnames(dh$ir) %in% rownames(m)]
r = bootstrapDevAS(m,dh,ss,N=1,only.samples = T)
r = r[[1]]
table(r$tissue,r$stage)-table(m$tissue,m$stage)

registerDoMC(3)
m = meta[!is.na(meta$mouse.stage),]
b = bootstrapDevAS(m,dh[which(dh$seg$sites=='ad')[1:1000],],rownames(species)[c(1,2,3)],N=3,only.samples = F)
bb = b
b = bb[[2]]
b = b[all.anns$mouse[rownames(b),'sites']=='ad',]
l = unique(as.vector(b))
apply(b,2,function(x)table(factor(x,levels=l)))

getDetAS = function(m,psi){
	f = rep(F,nrow(psi))
	na = is.na(psi)
	for(t in unique(m$tissue)){
		cinx = m$tissue==t
		f = f | (apply(!na[,cinx],1,mean) > 0.6 & apply(psi[,cinx],1,function(x){x=x[!is.na(x)];sum(x>0.1 & x < 0.9)>3}))
	}
	na = na[f,]
	sapply(unique(m$tissue),function(t){
		sum(apply(!na[,m$tissue==t],1,mean)>0.6)
		})
}
o = getDetAS(m,dh$ir[dh$seg$sites=='ad',m$sam.id])
b = getDetAS(r,dh$ir[dh$seg$sites=='ad',r$sam.id])

dd = readRDS('Rdata/hqmrboc.subsample/detected.as.all.direct.Rdata')
dd['mouse',,'ad']
plot(o,b)
abline(a=0,b=1)



zz = 	readRDS('Rdata/hqmrboc.subsample/chicken.12-0dpb/devAS.subsampling-mouse-mouse.Rdata')
zz[[1]][1:10,]
zzz = zz[[1]][all.anns$mouse[rownames(zz[[1]]),'sites']=='ad',]
apply(zzz,2,function(x)sum(x %in% c('c','n','u','d','ud','du')))
table(zzz[,1])
