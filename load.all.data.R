#setwd('~/skoltech/projects/evo.devo/')
options(stringsAsFactors = FALSE)
source('code/r.functions/load.all.data.F.R')
library(SAJR)
library(GenomicRanges)
library(plyr)
library(doMC)
library(cluster)


species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
# orth.seg.ad.all.id = readRDS('Rdata/orth.seg.ad.all.id.Rdata')
# orth.seg.ad.all = readRDS('Rdata/orth.seg.ad.all.Rdata')
# orth.seg.ad.tsm = readRDS('Rdata/orth.seg.ad.tsm.Rdata')
# orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
# hmo.seg.ad.all = readRDS('Rdata/hmo.seg.ad.all.Rdata')
# meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
# psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
# ens.ge = readRDS('Rdata/ens.ge.Rdata')
# ens.ge.cod = readRDS('Rdata/ens.ge.cod.Rdata')
# my.ge = readRDS('Rdata/my.ge.Rdata')
# my.ge.cod.tsm = readRDS('Rdata/my.ge.cod.tsm.Rdata')
# seg2ens = readRDS('Rdata/seg2ens.Rdata')
# alts = readRDS('Rdata/alts.Rdata')
# alts.filt = readRDS('Rdata/alts.filt.Rdata')
# all.anns = readRDS('Rdata/all.anns.Rdata')
# anns = readRDS('Rdata/anns.Rdata') 
# orth.ens.genes = readRDS('Rdata/orth.ens.genes.Rdata')
# params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
# params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)
write.table(all.anns$human[all.anns$human$sites=='ad',c('chr_id','start','stop','strand')],sep='\t',quote = F,col.names = F,file='processed/gnomad201/human.ad.tab')

ens.descr = unique(read.table('input/hs.37.73.gene.descr.txt',sep=',',quote='"',header=TRUE))
rownames(ens.descr) = ens.descr$Ensembl.Gene.ID
ens.descr$Description = sapply(strsplit(ens.descr$Description,' [',TRUE),'[',1)

ens.descr.mm = unique(read.table('input/mm.38.84.gene.descr.txt',sep=',',quote='"',header=TRUE))
rownames(ens.descr.mm) = ens.descr.mm$Ensembl.Gene.ID
ens.descr.mm$Description = sapply(strsplit(ens.descr.mm$Description,' [',TRUE),'[',1)

#### prepare data #######
# species = data.frame(species=c('human','macaque','mouse','rat','rabbit','opossum','chicken'),
# 										 short = c('h','q','m','r','b','o','c'),
# 										 gestation = c(280,165,20,21,30,15,21),
# 										 weaning = c(639,292,22,25,26,53,NA),
# 										 sex.maturity.m = c(14*365,5.5*365,42,90,NA,122,6*30),
# 										 sex.maturity.f = c(13*365,3.4*365,42,70,180,122,6*30),
# 										 max.longevity = c(90,40,4,3.8,9,5.1,30), #http://genomics.senescence.info/species/
# 										 name = c('Homo sapiens','Macaca mulatta','Mus musculus','Rattus norvegicus','Oryctolagus cuniculus','Monodelphis domestica','Gallus gallus'),
# 										 row.names = 1) #changed rabbit to 180 instead of 730. Margarida comment: So is a trait with a lot of variation between the many different races of rabbits. Our rabbits are surely adults after 6 months but not yet fully adults at 12wpc.
# saveRDS(species,'Rdata/species.Rdata')
#see 'comp.read.counts.R' for source code for this file. It is simply extracted from filenames
# meta = read.csv('input/sample.info.csv',row.names = 1)
# meta = meta[,c(1:5,9,7)]
# meta[,-c(1,7)] = apply(meta[,-c(1,7)],2,tolower)
# meta[1:10,]
# meta$fname = gsub('.sorted.bam','',meta$fname)
# # set age from conseption
# meta$days = as.numeric(meta$age)
# meta = tranfromAges(meta,'dpb',1,TRUE)
# meta = tranfromAges(meta,'wpb',7,TRUE)
# meta = tranfromAges(meta,'mpb',30,TRUE)
# meta = tranfromAges(meta,'ypb',365,TRUE)
# meta = tranfromAges(meta,'dph',1,TRUE)
# meta = tranfromAges(meta,'wph',7,TRUE)
# meta = tranfromAges(meta,'w',7,FALSE)
# meta = tranfromAges(meta,'d',1,FALSE)
# meta$days[meta$species=='chicken' & meta$age == 'adult'] = 155+21
# meta$days[meta$age=='cs13'] = 32
# meta$days[meta$age=='cs14'] = 33
# meta$days[meta$age=='cs16'] = 39
# meta$days[meta$age=='cs17'] = 41
# meta$days[meta$age=='cs18'] = 44
# meta$days[meta$age=='cs19'] = 46
# meta$days[meta$age=='cs20'] = 49
# meta$days[meta$age=='cs21'] = 51
# meta$days[meta$age=='cs22'] = 53
# meta$days[meta$age=='cs23'] = 56
# #
# table(meta$age[is.na(meta$days)])
# table(meta$species,meta$tissue)
# short.tissue = toupper(substr(meta$tissue,1,1))
# short.tissue[meta$tissue=='cam'] = 'Cam'
# short.tissue[meta$tissue=='hindbrain'] = 'Hb'
# short.tissue[meta$tissue=='kidneytestis'] = 'Kt'
# table(short.tissue,meta$tissue)
# meta$sex = substr(meta$sex,1,1)
# 
# short.sam.name = paste(toupper(species[meta$species,'short']),short.tissue,meta$days,toupper(meta$sex),sep='')
# short.sam.name = do.call(rbind,lapply(split(data.frame(1:length(short.sam.name),short.sam.name),short.sam.name),function(x){x[,2] = paste(x[,2],1:nrow(x),sep='_');x}))
# meta$name = short.sam.name[order(short.sam.name[,1]),2]
# meta = meta[order(meta$species,meta$tissue,meta$days,meta$sex),]
# meta[1:10,]
# rownames(meta) = meta$name
# meta$pch = as.numeric(factor(meta$species))
# meta$age.use = log(meta$days)
# 
# meta$col = params$tissue.col[meta$tissue]
# 
# min.age = sapply(split(meta$age.use,meta$species),min,na.rm=T)[meta$species]
# max.age = sapply(split(meta$age.use,meta$species),max,na.rm=T)[meta$species]
# 
# meta$cex = (meta$age.use-min.age)
# meta$cex = 0.4 + meta$cex/max.age*1.5
# meta$age.rank = NA
# for(s in rownames(species)){
# 	f = meta$species==s
# 	s2a = sapply(split(meta$days[f],meta$stage[f]),mean)
# 	s2a = sort(s2a)
# 	s2a[] = 1:length(s2a)
# 	meta$age.rank[f] = s2a[meta$stage[f]]
# }
# saveRDS(meta,'Rdata/whole.meta.Rdata')
# meta = readRDS('Rdata/whole.meta.Rdata')
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
# 
# registerDoMC(7)
# tmp = laply(rownames(species),function(s){
# 	print(s)
# 	m = meta[meta$species==s,]
# 	t = f(paste('processed/annotation/all.species/merged/',s,'.sajr',sep=''),paste('processed/sajr/uniq/',m$fname,sep=''),m$name,paste('processed/annotation/all.species/ensambl/',ens[s],".gtf.gz",sep=''))
# 	alts = makeAlts(t$seg,paste('processed/annotation/all.species/merged/',s,'.sajr',sep=''),TRUE)
# 	t$seg = assignSegs2Alts(alts,t$seg)
# 	saveRDS(t,paste('Rdata/',s,'.as.u.all.Rdata',sep=''))
# 	print('NA freq:')
# 	print(table(is.na(t$seg$gene_id)))
# 	saveRDS(f(paste('processed/annotation/all.species/merged/',s,'.sajr',sep=''),paste('processed/sajr/mult/',m$fname,sep=''),m$name,paste('processed/annotation/all.species/ensambl/',ens[s],".gtf.gz",sep='')),paste('Rdata/',s,'.as.m.all.Rdata',sep=''))
# 	gc()
# 	list(seg=t$seg,alts=alts,alts.filt = filterAlts(alts,filter.introns = TRUE))
# },.parallel = TRUE)

# all.anns = tmp[,1]
# alts = tmp[,2]
# alts.filt = tmp[,3]
# names(all.anns) = names(alts) = names(alts.filt) = rownames(species)
# 
# 
# saveRDS(alts,'Rdata/alts.Rdata')
# saveRDS(alts.filt,'Rdata/alts.filt.Rdata')
# saveRDS(all.anns,'Rdata/all.anns.Rdata')

#load GE
# my.ge = vector('list',7)
# names(my.ge) = rownames(species)
# for(sp in rownames(species)){
# 	print(sp)
# 	m = meta[meta$species==sp,]
# 	my.ge[[sp]] = loadGData(paste('processed/annotation/all.species/merged/',sp,'.sajr',sep=''),paste('processed/sajr/uniq/',m$fname,sep=''),m$name)
# 	s = loadSAData(paste('processed/annotation/all.species/merged/',sp,'.sajr',sep=''))
# 	my.ge[[sp]] = my.ge[[sp]][my.ge[[sp]]$gene$chr_id != 'MT',]
# 	my.ge[[sp]] = calcRPKM(s,my.ge[[sp]],TRUE)
# 	my.ge[[sp]]$gene$cod = rownames(my.ge[[sp]]$gene) %in% all.anns[[sp]]$gene_id[all.anns[[sp]]$cod!='n']
# }

# saveRDS(my.ge,'Rdata/my.ge.Rdata')
# 
# ens.ge = vector('list',7)
# names(ens.ge) = rownames(species)
# for(sp in rownames(species)){
# 	print(sp)
# 	m = meta[meta$species==sp,]
# 	ens.ge[[sp]] = loadGData(paste('processed/annotation/all.species/ensambl/',ens[sp],'.sajr',sep=''),paste('processed/sajr/ens.uniq/',m$fname,sep=''),m$name)
# 	s = loadSAData(paste('processed/annotation/all.species/ensambl/',ens[sp],'.sajr',sep=''))
# 	ens.ge[[sp]] = ens.ge[[sp]][ens.ge[[sp]]$gene$chr_id != 'MT',]
# 	ens.ge[[sp]] = calcRPKM(s,ens.ge[[sp]],TRUE)
# 	bt = getGeneBiotype(paste('processed/annotation/all.species/ensambl/',ens[sp],'.gtf.gz',sep=''))
# 	ens.ge[[sp]]$gene$biotype = bt[rownames(ens.ge[[sp]]$gene)]
# }

# saveRDS(ens.ge ,'Rdata/ens.ge.Rdata')

# filter: remove additional tissues, remove Margarida's defined outliers, and filter exons by types
# marg.names = read.csv('processed/GE.from.marg/marg.sample.names.csv')
# marg.names = marg.names[marg.names$X != 'files' & marg.names$X != 'group',]
# marg.names$lib.id = sapply(strsplit(marg.names$files,'.',TRUE),'[',1)
# marg.names$species = sapply(strsplit(marg.names$files,'.',TRUE),'[',2)
# table(marg.names$lib.id %in% meta$lib.id)
# table(meta$lib.id %in% marg.names$lib.id)#,meta$tissue)
#
# meta$tissue.orig = meta$tissue
# meta$tissue[meta$tissue=='hindbrain'] = 'cerebellum'
# meta$tissue[meta$tissue=='forebrain'] = 'brain'
# meta$tissue[meta$tissue=='wholebrain'] = 'brain'
# meta$tissue[meta$tissue=='kidneytestis'] = 'kidney'
# meta = meta[meta$tissue %in% c('brain','cerebellum','heart','kidney','liver','ovary','testis'),]
# meta = meta[meta$lib.id!='3020sTS',] # it seems to be an outlier
# meta = meta[meta$lib.id %in% marg.names$lib.id,]
# meta$marg.name = setNames(marg.names$X,marg.names$lib.id)[meta$lib.id]
# meta$marg.stage = setNames(marg.names$group,marg.names$lib.id)[meta$lib.id]
# meta$col = params$tissue.col[meta$tissue]
# meta$marg.name[meta$marg.name=='Kidney.P152.3'] = 'Kidney.P183.2' #seem  there is a mistake
# saveRDS(meta,'Rdata/main.set.meta.Rdata')

# make filtered sets
# only internal alternative segments detected in more than 60% of samples of at least one tissue and with at leats 4 samples within [0.1,0.9] (in one tissue)
# anns = vector('list',7)
# names(anns) = rownames(species)
# for(s in rownames(species)){
# 	t = readRDS(paste('Rdata/',s,'.as.u.all.Rdata',sep=''))
# 	st = length(t)
# 	t = t[t$seg$chr_id !='MT' & t$seg$position=='INTERNAL' & t$seg$type!='EXN',colnames(t$ir) %in% meta$name]
# 	st[2] = length(t)
# 	f = rep(FALSE,length(t))
# 	na = is.na(t$ir)
# 	for(tis in unique(meta$tissue)){
# 		cinx = meta[colnames(na),'tissue']==tis
# 		f = f | (apply(!na[,cinx],1,mean) > 0.6 & apply(t$ir[,cinx],1,function(x){x=x[!is.na(x)];sum(x>0.1 & x < 0.9)>3}))
# 	}
# 	t = t[f,]
# 	st[3] = length(t)
# 	cat(s,st,"\t")
# 	saveRDS(t,paste('Rdata/',s,'.as.u.filtered.Rdata',sep=''))
# 	anns[[s]] = t$seg
# }
# saveRDS(anns,'Rdata/anns.Rdata')
# 
# my.ge = readRDS('Rdata/my.ge.Rdata')
# my.ge$human$gene[1:2,]
# my.ge.cod = lapply(my.ge,function(x){x[x$gene$cod & x$gene$chr_id != 'MT',]})
# for(sp in rownames(species)){
# 	print(sp)
# 	s = loadSAData(paste('processed/annotation/all.species/merged/',sp,'.sajr',sep=''))
# 	my.ge.cod[[sp]] = calcRPKM(s,my.ge.cod[[sp]],TRUE)
# }
# saveRDS(my.ge.cod,'Rdata/my.ge.cod.Rdata')

# ens.ge.cod = lapply(ens.ge,function(x){x[x$gene$biotype=='protein_coding' & x$gene$chr_id != 'MT',]})
# for(sp in rownames(species)){
# 	print(sp)
# 	m = meta[meta$species==sp,]
# 	s = loadSAData(paste('processed/annotation/all.species/ensambl/',ens[sp],'.sajr',sep=''))
# 	ens.ge.cod[[sp]] = calcRPKM(s,ens.ge.cod[[sp]],TRUE)
# }
# saveRDS(ens.ge.cod,'Rdata/ens.ge.cod.Rdata')

# ## link segs to ENS #########
# seg2ens = list()
# for(s in rownames(species)){
# 	print(s)
# 	seg2ens[[s]] = makeSeg2Ens(s)
# }
# saveRDS(seg2ens,'Rdata/seg2ens.Rdata')
# 
# ##### alternatives to ens ######
# alt2ens = lapply(rownames(species),function(s){
# 	print(s)
# 	
# 	lapply(split(seg2ens[[s]],all.anns[[s]]$alt.id),function(x){
# 		unique(unlist(x))
# 	})
# })
# names(alt2ens) = rownames(species)
# saveRDS(alt2ens,'Rdata/alt2ens.Rdata')
# 
# #### gene to ens ###########
# gene2ens = lapply(rownames(species),function(s){
# 	print(s)
# 	lapply(split(seg2ens[[s]],all.anns[[s]]$gene_id),function(x){
# 		unique(unlist(x))
# 	})
# })
# names(gene2ens) = rownames(species)
# saveRDS(gene2ens,'Rdata/gene2ens.Rdata')

### species-tissue-stage means
# meta.tsm = split(meta,paste(meta$species,meta$tissue,meta$stage))
# meta.tsm = lapply(meta.tsm,function(x){
# 	data.frame(species=x$species[1],tissue=x$tissue[1],stage=x$stage[1],col=x$col[1],cex=mean(x$cex),pch=x$pch[1])
# 	})
# meta.tsm = do.call(rbind,meta.tsm)
# meta.tsm$days = NA
# for(s in unique(meta.tsm$species))
# 	for(g in unique(meta.tsm$stage))
# 		meta.tsm$days[meta.tsm$stage == g & meta.tsm$species == s] = mean(meta$days[meta$stage == g & meta$species == s])
# meta.tsm$age.use = log(meta.tsm$days)
# rownames(meta.tsm) = paste(meta.tsm$species,meta.tsm$tissue,meta.tsm$stage)
# 
# meta.tsm$age.rank = NA
# for(s in rownames(species)){
# 	f = meta.tsm$species==s
# 	s2a = sapply(split(meta.tsm$days[f],meta.tsm$stage[f]),mean)
# 	s2a = sort(s2a)
# 	s2a[] = 1:length(s2a)
# 	meta.tsm$age.rank[f] = s2a[meta.tsm$stage[f]]
# }

# 
# psi.tsm = vector('list',nrow(species))
# names(psi.tsm) = rownames(species)
# for(s in rownames(species)){
# 	cat(s)
# 	t = readRDS(paste('Rdata/',s,'.as.u.filtered.Rdata',sep=''))
# 	m = meta[colnames(t$ir),]
# 	psi.tsm[[s]] = calcMeanCols(t$ir,paste(m$species,m$tissue,m$stage))
# }

# psi.tsm.ms = vector('list',nrow(species))
# names(psi.tsm.ms) = rownames(species)
# for(s in rownames(species)){
# 	cat(s)
# 	t = readRDS(paste('Rdata/',s,'.as.u.filtered.Rdata',sep=''))
# 	m = meta[colnames(t$ir),]
# 	psi.tsm.ms[[s]] = calcMeanCols(t$ir,paste(m$species,m$tissue,m$marg.stage))
# }
# for(s in names(psi.tsm.ms))colnames(psi.tsm.ms[[s]]) = tolower(colnames(psi.tsm.ms[[s]]))
# saveRDS(psi.tsm.ms,'Rdata/psi.tsm.ms.Rdata')

# meta.tsm.ms = split(meta,paste(meta$species,meta$tissue,meta$marg.stage))
# meta.tsm.ms = lapply(meta.tsm.ms,function(x){
# 	data.frame(species=x$species[1],tissue=x$tissue[1],stage=x$marg.stage[1],col=x$col[1],cex=mean(x$cex),pch=x$pch[1])
# 	})
# meta.tsm.ms = do.call(rbind,meta.tsm.ms)
# meta.tsm.ms$stage = tolower(meta.tsm.ms$stage)
# meta.tsm.ms$days = NA
# for(s in unique(meta.tsm.ms$species))
# 	for(g in unique(meta.tsm.ms$stage))
# 		meta.tsm.ms$days[meta.tsm.ms$stage == g & meta.tsm.ms$species == s] = mean(meta$days[tolower(meta$marg.stage) == g & meta$species == s])
# 
# meta.tsm.ms$age.use = log(meta.tsm.ms$days)
# rownames(meta.tsm.ms) = paste(meta.tsm.ms$species,meta.tsm.ms$tissue,meta.tsm.ms$stage)
# 
# meta.tsm.ms$age.rank = NA
# for(s in rownames(species)){
# 	f = meta.tsm.ms$species==s
# 	s2a = sapply(split(meta.tsm.ms$days[f],meta.tsm.ms$stage[f]),mean)
# 	s2a = sort(s2a)
# 	s2a[] = 1:length(s2a)
# 	meta.tsm.ms$age.rank[f] = s2a[meta.tsm.ms$stage[f]]
# }
# 
# saveRDS(meta.tsm.ms,'Rdata/meta.tsm.ms.Rdata')
# 
# my.ge.tsm = calcGEtsm(my.ge,meta)
# my.ge.cod.tsm = calcGEtsm(my.ge.cod,meta)
# ens.ge.tsm = calcGEtsm(ens.ge,meta)
# ens.ge.cod.tsm = calcGEtsm(ens.ge.cod,meta)
# 
# saveRDS(meta.tsm,'Rdata/meta.tsm.Rdata')
# saveRDS(psi.tsm,'Rdata/psi.tsm.Rdata')
# saveRDS(my.ge.tsm,'Rdata/my.ge.tsm.Rdata')
# saveRDS(my.ge.cod.tsm,'Rdata/my.ge.cod.tsm.Rdata')
# saveRDS(ens.ge.tsm,'Rdata/ens.ge.tsm.Rdata')
# saveRDS(ens.ge.cod.tsm,'Rdata/ens.ge.cod.tsm.Rdata')

# tmp = do.call(rbind,lapply(all.anns,function(x)cbind(rownames(x),x$sites)))
# write.table(tmp,'processed/annotation/all.species/merged/seg2sites.tab',quote = F,row.names = F,col.names = F,sep='\t')
# o = read.csv('input/ens.orths.txt.gz')
# o = o[apply(o=='',1,sum)==0,]
# dim(o)
# f = rep(T,nrow(o))
# for(i in 1:ncol(o)){
# 	t = table(o[,i])
# 	f = f & o[,i] %in% names(t)[t==1]
# }
# table(f)
# o = o[f,]
# 
# orth.ens.genes = o
# rownames(orth.ens.genes) = orth.ens.genes[,1]
# colnames(orth.ens.genes) = rownames(species)
# saveRDS(orth.ens.genes,'Rdata/orth.ens.genes.Rdata')


### orth segments
# orth.seg.ad = loadAltOrthSegs(c('processed/orth.segs/only.ad/hqmrboc.0.6.orth.ad.segs','processed/orth.segs/only.ad/hqmrboc.0.0.orth.ad.segs'),only.filtered = FALSE)
# o.by.neig = findOrthExonsByNeighbors(all.anns,orth.seg.ad,'human')
# t = table(o.by.neig$same.frame,floor(o.by.neig$length.min/o.by.neig$length.max*20)/20)
# 
# pdf('figures/hqmrboc/orth.seg.by.neighbors.hqmrboc.pdf',w=8,h=4)
# par(mfrow=c(1,2),tck=-0.02,mgp=c(1.3,0.4,0),mar=c(3,3,1.5,0),oma=c(0,0,2,1))
# barplot(t,xlab='min(length)/max(length)',ylab='# of orth exons',legend.text=c('different frame','same frame'),args.legend=list(x='topleft'))
# barplot(sweep(t,2,apply(t,2,sum),'/'),xlab='min(length)/max(length)',ylab='proportion of orth exons',)
# abline(h=2/3,col='red')
# dev.off()

# orth.seg.ad = rbind(orth.seg.ad,o.by.neig[o.by.neig$length.min/o.by.neig$length.max >= 0.8,1:7])
# orth.seg.ad = loadInfoForOrths(orth.seg.ad)
# are orth segs come from oth ens genes?
# 
#-2 more than 1 gene (at least one seg assigned to more than 1 gene),
#-1 no orth genes (at least one seg do not have ens gene; -2 is more important than -1)
# 0 contradict (orth segs contradict orth genes),
# 1 genes are not from one-to-one orthologs
# 2 everythings is OK
# f = rep(NA,length(orth.seg.ad$human))
# eidx = unlist(orth.ens.genes)
# eidx = setNames(rep(1:nrow(orth.ens.genes),times=ncol(orth.ens.genes)),eidx)
# for(i in 1:length(f)){
# 	cat('\r',i)
# 	egids = lapply(names(orth.seg.ad),function(s){seg2ens[[s]][[rownames(orth.seg.ad[[s]]$seg)[i]]]})
# 	
# 	# if(sum(sapply(egids,length) > 1)>0){
# 	# 	f[i] = -2
# 	# }else if(sum(sapply(egids,length) == 0)>0){
# 	# 	f[i] = -1
# 	# }else{
# 	# 	egids = eidx[as.character(egids)]
# 	# 	if(sum(is.na(egids))>0)
# 	# 		f[i] = 1
# 	# 	else if(length(unique(egids)) == 1)
# 	# 		f[i] = 2
# 	# 	else
# 	# 		f[i] = 0
# 	# }
# 	f[i] = max(table(eidx[unlist(egids)]))
# }
# 
# f[is.infinite(f)] = 0
# table(f)
# orth.seg.ad$human$seg$north = f
# orth.seg.ad$human$seg$ens.orth.filter = f
# table(orth.seg.ad$human$seg$ens.orth.filter)
# saveRDS(orth.seg.ad,'Rdata/orth.seg.ad.all.Rdata')
# orth.seg.ad.all.id = sapply(orth.seg.ad,function(x){rownames(x$seg)})
# rownames(orth.seg.ad.all.id) = orth.seg.ad.all.id[,1]
# saveRDS(orth.seg.ad.all.id,'Rdata/orth.seg.ad.all.id.Rdata')

# orth.seg.ad.all = readRDS('Rdata/orth.seg.ad.all.Rdata')
# orth.seg.ad.all = lapply(orth.seg.ad.all,function(x){x$seg})
# lens = sapply(orth.seg.ad.all,function(x){x$length})
# len.rate = apply(lens,1,function(x){min(x)/max(x)})
# len.frames.no = apply(lens,1,function(x){length(unique(x%%3))})
# 
# 
# pdf('figures/hqmrboc/orth.seg.length.comparison.pdf',w=8,h=8)
# par(mfrow=c(2,2),tck=-0.02,mgp=c(1.3,0.4,0),mar=c(5,3,1.5,0),oma=c(0,0,2,1))
# hist(len.rate,xlab='max/man orth exon length',main='')
# z=table(len.rate==1,orth.seg.ad.all$human$ens.orth.filter)
# names = paste(c('mult. ens','no ens','contradict','not 1-to-1','agree'),'\n(',apply(z,2,sum),')',sep='')
# barplot(sweep(z,2,apply(z,2,sum),'/')['FALSE',],names.arg = names,las=3,ylab='Proportion of exons with different length')
# 
# 
# z=table(len.frames.no[len.rate!=1]==1,orth.seg.ad.all$human$ens.orth.filter[len.rate!=1])
# names = paste(c('mult. ens','no ens','contradict','not 1-to-1','agree'),'\n(',apply(z,2,sum),')',sep='')
# barplot(sweep(z,2,apply(z,2,sum),'/')['FALSE',],names.arg = names,las=3,ylab='Proportion of exons with different frames',main='Only exons with variable size')
# dev.off()
# a = sapply(orth.seg.ad.all,function(x)x$seg$type)
# table(apply(a=='ALT',1,sum))
# f = apply(a=='ALT',1,sum) > 0
# table(f)
# orth.seg.ad = lapply(orth.seg.ad.all,function(x)x[f,])
# saveRDS(orth.seg.ad,'Rdata/orth.seg.ad.Rdata')
# orth.seg.ad.tsm = lapply(orth.seg.ad,function(x){x$ir = x$ir[,colnames(x$ir) %in% rownames(meta)];m = meta[colnames(x$ir),];calcMeanCols(x$ir,paste(m$species,m$tissue,m$stage))})
# saveRDS(orth.seg.ad.tsm,'Rdata/orth.seg.ad.tsm.Rdata')

# for(s in 1:7){
# 	cat(s)
# 	r=countNumberOfCassettesInAlt(alts.filt[[s]])
# 	names(r) = rownames(alts.filt[[s]])
# 	anns[[s]]$alt.size = r[anns[[s]]$alt.id]
# 	anns[[s]]$alt.size[is.na(anns[[s]]$alt.size)] = 0
# 	all.anns[[s]]$alt.size = r[all.anns[[s]]$alt.id]
# 	all.anns[[s]]$alt.size[is.na(all.anns[[s]]$alt.size)] = 0
# }
# saveRDS(anns,'Rdata/anns.Rdata')
# saveRDS(all.anns,'Rdata/all.anns.Rdata')

# orth.seg.ad.all = readRDS('Rdata/orth.seg.ad.all.Rdata')
# orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
# all.anns = readRDS('Rdata/all.anns.Rdata')
# anns = readRDS('Rdata/anns.Rdata')
# for(s in rownames(species)){
# 	print(s)
# 	a = readRDS(paste0('Rdata/',s,'.as.u.all.Rdata'))
# 	f = readRDS(paste0('Rdata/',s,'.as.u.filtered.Rdata'))
# 	if(sum(rownames(a$seg)!=rownames(all.anns[[s]]))>0)
# 		stop(1)
# 	a$seg = all.anns[[s]]
# 	if(sum(rownames(f$seg)!=rownames(anns[[s]]))>0)
# 		stop(1)
# 	f$seg = anns[[s]]
# 	orth.seg.ad.all[[s]]$seg = all.anns[[s]][rownames(orth.seg.ad.all[[s]]$seg),]
# 	orth.seg.ad[[s]]$seg = all.anns[[s]][rownames(orth.seg.ad[[s]]$seg),]
# 	saveRDS(a,paste0('Rdata/',s,'.as.u.all.Rdata'))
# 	saveRDS(f,paste0('Rdata/',s,'.as.u.filtered.Rdata'))
# }
# saveRDS(orth.seg.ad.all,'Rdata/orth.seg.ad.all.Rdata')
# saveRDS(orth.seg.ad,'Rdata/orth.seg.ad.Rdata')


#### calc fraction of exons that can be liftovered #####
tc = calcProportionOfLiftoveredExons(all.anns,orth.seg.ad.all.ids,'c','human')

par(mfrow=c(1,2),tck=-0.02,mgp=c(1.3,0.4,0),mar=c(3,6,1.5,0),oma=c(0,0,2,1))
hist(tc$orth.freq,0:20/20)
hist(tc$orth.freq[tc$exn.cnt>9],0:20/20)

### hmo
# hmo.seg.ad = loadAltOrthSegs(c('processed/orth.segs/only.ad/hmo.0.6.orth.ad.segs','processed/orth.segs/only.ad/hmo.0.0.orth.ad.segs'),only.filtered = FALSE)
# o.by.neig = findOrthExonsByNeighbors(all.anns[colnames(hmo.seg.ad)],hmo.seg.ad,'human')
# t = table(o.by.neig$same.frame,floor(o.by.neig$length.min/o.by.neig$length.max*20)/20)
# 
# 
# pdf('figures/hmo/orth.seg.by.neighbors.hqmo.pdf',w=8,h=4)
# par(mfrow=c(1,2),tck=-0.02,mgp=c(1.3,0.4,0),mar=c(3,3,1.5,0),oma=c(0,0,2,1))
# barplot(t,xlab='min(length)/max(length)',ylab='# of orth exons',legend.text=c('different frame','same frame'),args.legend=list(x='topleft'))
# barplot(sweep(t,2,apply(t,2,sum),'/'),xlab='min(length)/max(length)',ylab='proportion of orth exons',)
# abline(h=2/3,col='red')
# dev.off()
# 
# hmo.seg.ad = rbind(hmo.seg.ad,o.by.neig[o.by.neig$length.min/o.by.neig$length.max >= 0.8,1:3])
# hmo.seg.ad = loadInfoForOrths(hmo.seg.ad)
# 
# f = rep(NA,length(hmo.seg.ad$human))
# eidx = unlist(orth.ens.genes)
# eidx = setNames(rep(1:nrow(orth.ens.genes),times=ncol(orth.ens.genes)),eidx)
# for(i in 1:length(f)){
# 	cat('\r',i)
# 	egids = lapply(names(hmo.seg.ad),function(s){seg2ens[[s]][[rownames(hmo.seg.ad[[s]]$seg)[i]]]})
# 	if(sum(sapply(egids,length) > 1)>0){
# 		f[i] = -2
# 	}else if(sum(sapply(egids,length) == 0)>0){
# 		f[i] = -1
# 	}else{
# 		egids = eidx[as.character(egids)]
# 		if(sum(is.na(egids))>0)
# 			f[i] = 1
# 		else if(length(unique(egids)) == 1)
# 			f[i] = 2
# 		else
# 			f[i] = 0
# 		
# 	}
# }
# hmo.seg.ad$human$seg$ens.orth.filter = f
# table(hmo.seg.ad$human$seg$ens.orth.filter)
# saveRDS(hmo.seg.ad,'Rdata/hmo.seg.ad.all.Rdata')
# a = sapply(hmo.seg.ad.all,function(x)x$seg$type)
# f = apply(a=='ALT',1,sum) > 0
# table(f)
# hmo.seg.ad = lapply(hmo.seg.ad.all,function(x)x[f,])
# saveRDS(hmo.seg.ad,'Rdata/hmo.seg.ad.Rdata')

# look on CDS position #####
ens = c(mouse='Mus_musculus.GRCm38.84',
				rat='Rattus_norvegicus.Rnor_5.0.79',
				rabbit='Oryctolagus_cuniculus.OryCun2.0.84',
				human='Homo_sapiens.GRCh37.73',
				macaque='Macaca_mulatta.MMUL_1.84',
				opossum='Monodelphis_domestica.BROADO5.84',
				chicken='Gallus_gallus.Galgal4.84')[rownames(species)]
# ens.as.ann = lapply(ens,function(f){
# 	print(f)
# 	r = loadSAData(paste0('processed/annotation/all.species/ensambl/',f,'.sajr'))
# 	r = setSplSiteTypes(r,paste0('processed/annotation/all.species/ensambl/',f,'.sajr'))
# 	r = addIsCogingByEnsGTF(paste0('processed/annotation/all.species/ensambl/',f,'.gtf.gz'),r)
# 	r$seg$length = r$seg$stop - r$seg$start + 1
# 	r$seg = addExonNumber(r$seg)
# 	r$seg = addCDSPos2Ann(r$seg)
# 	r$seg
# })
# ens.as.ann = ens.as.ann[names(all.anns)]
# saveRDS(ens.as.ann,'Rdata/ens.as.ann.Rdata')

# loadExonCDSPositionPerTranscript = function(f){
# 	print(f)
# 	gc()
# 	m = read.table(f,sep='\t')
# 	m = m[m$V3 %in% c('exon','CDS'),]
# 	m = cbind(m[,-9],t(sapply(strsplit(m$V9,';\\s?'),function(x){x = strsplit(x,' ');setNames(sapply(x,'[',2),sapply(x,'[',1))[c('transcript_id','gene_id','gene_biotype','exon_number')]})))
# 	m = split(m,m$transcript_id)
# 
# 	r = lapply(m,function(x){
# 		e = x[x$V3 == 'exon',]
# 		c = x[x$V3 == 'CDS',]
# 		e = e[order(e$V4),]
# 		e$internal = e$exon_number != 1 & e$exon_number != max(e$exon_number)
# 		c = c[order(c$V4),]
# 		e$cds.pos = 'non-cod'
# 		if(nrow(c)==0)
# 			return(e)
# 		j = 1
# 		e$cds.pos = '-'
# 		for(i in 1:nrow(e)){
# 			if(e$V5[i] < c$V4[j]){
# 				e$cds.pos[i] = ifelse(e$V7[i]=='+','5utr','3utr')
# 			}else if(e$V5[i] == c$V5[j] & e$V4[i] < c$V4[j]){
# 				e$cds.pos[i] = ifelse(e$V7[i]=='+','first','last')
# 				j = j + 1
# 			}else if(e$V5[i] == c$V5[j] & e$V4[i] == c$V4[j]){
# 				e$cds.pos[i] = 'cds'
# 				j = j + 1
# 			}else if(e$V5[i] > c$V5[j] & e$V4[i] == c$V4[j]){
# 				e$cds.pos[i] = ifelse(e$V7[i]=='+','last','first')
# 				j = j + 1
# 			}else if(e$V5[i] > c$V5[j] & e$V4[i] < c$V4[j]){
# 				e$cds.pos[i] = 'first-last'
# 				j = j + 1
# 			}else if(e$V4[i] > c$V5[j]){
# 				e$cds.pos[i] = ifelse(e$V7[i]=='+','3utr','5utr')
# 			}
# 			j = min(j,nrow(c))
# 		}
# 		e
# 	})
# 
# 	r = do.call(rbind,r)
# 	r$coor = paste(r$V1,r$V7,r$V4,r$V5)
# 	r
# }
# 
# ens.exon.transc.cds.pos = lapply(ens,function(f)loadExonCDSPositionPerTranscript(paste0('processed/annotation/all.species/ensambl/',f,'.gtf.gz')))
# saveRDS(ens.exon.transc.cds.pos,'Rdata/ens.exon.transc.cds.pos.Rdata')
ens.exon.transc.cds.pos = readRDS('Rdata/ens.exon.transc.cds.pos.Rdata')

# so, for human 2nd column gives information about transcript biotype, but it doesn't work in other species
lapply(ens.exon.transc.cds.pos,function(x)table(x$V2))
lapply(ens.exon.transc.cds.pos,function(x)table(x$gene_biotype))
sapply(ens.exon.transc.cds.pos,function(x)table(x$cds.pos))


# registerDoMC(4)
# ens.exon.cds.pos.cod.gene = llply(ens.exon.transc.cds.pos,function(r){
# 	gc()
# 	print(1)
# 	pc = r[r$gene_biotype=='protein_coding',]
# 	pc = split(pc,pc$coor)
# 	pc= do.call(rbind,llply(pc,function(x){
# 		cbind(x[1,c(1,4,5,7)],cds.poss=paste(sort(unique(x$cds.pos)),collapse = '|'))
# 	},.parallel=F))
# 	setNames(pc$cds.poss,rownames(pc))
# },.parallel=T)
# saveRDS(ens.exon.cds.pos.cod.gene,'Rdata/ens.exon.cds.pos.cod.gene.Rdata')


# idea with NMD-transcripts works only for human
# ens.exon.cds.pos.cod.gene = lapply(ens.exon.transc.cds.pos,function(r){
# 	gc()
# 	print(1)
# 	pc = r[r$gene_biotype=='protein_coding' & r$V2=='protein_coding',]
# 	pc = split(pc,pc$coor)
# 	pc=do.call(rbind,lapply(pc,function(x){
# 		cbind(x[1,c(1,4,5,7)],cds.poss=paste(sort(unique(x$cds.pos)),collapse = '|'))
# 	}))
# 
# 	nmd = r[r$gene_biotype=='protein_coding' & r$V2=='nonsense_mediated_decay',]
# 	table(nmd$coor %in% rownames(pc))
# 	nmd = nmd[!(nmd$coor %in% rownames(pc)),]
# 	nmd = split(nmd,nmd$coor)
# 	nmd=do.call(rbind,lapply(nmd,function(x){
# 		cbind(x[1,c(1,4,5,7)],cds.poss=paste0('nmd:',sort(unique(x$cds.pos)),collapse = '|'))
# 	}))
# 
# 	pc.nmd = rbind(pc,nmd)
# 	setNames(pc.nmd$cds.poss,rownames(pc.nmd))
# }) 
ens.exon.cds.pos.cod.gene = readRDS('Rdata/ens.exon.cds.pos.cod.gene.Rdata')

u = unique(unlist(ens.exon.cds.pos.cod.gene))
t = sapply(ens.exon.cds.pos.cod.gene,function(x)table(factor(x,levels = u)))
t[order(apply(t,1,sum),decreasing = T),]

anns     = lapply(    anns,function(x)cbind(x,coor=paste(x$chr_id,ifelse(x$strand==1,'+','-'),x$start,x$stop)))
all.anns = lapply(all.anns,function(x)cbind(x,coor=paste(x$chr_id,ifelse(x$strand==1,'+','-'),x$start,x$stop)))

t(sapply(names(anns),function(s){table(anns[[s]]$coor[anns[[s]]$sites=='ad'] %in% names(ens.exon.cds.pos.cod.gene[[s]]))}))
t(sapply(names(all.anns),function(s){table(all.anns[[s]]$coor[all.anns[[s]]$sites=='ad'] %in% names(ens.exon.cds.pos.cod.gene[[s]]))}))

smpl = function(x){
	x[grep('cds|first|last',x)] = 'cds'
	x[ grepl('5utr',x) & !grepl('3utr',x)] = '5utr'
	x[!grepl('5utr',x) &  grepl('3utr',x)] = '3utr'
	x
}

for(s in names(ens.exon.cds.pos.cod.gene)){
	anns[[s]]$ens.transc.cds.pos = NA
	f = anns[[s]]$sites=='ad' & anns[[s]]$coor %in% names(ens.exon.cds.pos.cod.gene[[s]])
	anns[[s]]$ens.transc.cds.pos[f] = ens.exon.cds.pos.cod.gene[[s]][anns[[s]]$coor[f]]
	anns[[s]]$ens.transc.cds.pos.smpl[f] = smpl(anns[[s]]$ens.transc.cds.pos[f])
	
	all.anns[[s]]$ens.transc.cds.pos = NA
	f = all.anns[[s]]$sites=='ad' & all.anns[[s]]$coor %in% names(ens.exon.cds.pos.cod.gene[[s]])
	all.anns[[s]]$ens.transc.cds.pos[f] = ens.exon.cds.pos.cod.gene[[s]][all.anns[[s]]$coor[f]]
	all.anns[[s]]$ens.transc.cds.pos.smpl[f] = smpl(all.anns[[s]]$ens.transc.cds.pos[f])
}

t = table(anns$human$ens.transc.cds.pos.smpl,anns$human$cds.pos)
t = table(all.anns$human$ens.transc.cds.pos,all.anns$human$cds.pos)
t[order(apply(t,1,sum),decreasing = T),]


f = anns$human$sites=='ad'
table(is.na(anns$human$ens.transc.cds.pos.smpl[f]),anns$human$cds.pos[f])

u = unique(unlist(lapply(anns,function(x)x$cds.pos)))
t = t(sapply(anns,function(x)table(factor(x$cds.pos,levels = u))))

u = unique(unlist(lapply(all.anns,function(x)x$ens.transc.cds.pos.smpl)))
t = t(sapply(all.anns,function(x)table(factor(x$ens.transc.cds.pos.smpl,levels = u))))

t[order(apply(t,1,sum),decreasing = T),order(apply(t,2,sum),decreasing = T)]

# annotate by lncRNA #####
for(s in rownames(species)){
	print(s)
	h = read.table(paste0('input/lncRNA/20190411/',s,'.xlocs.gtf'),sep='\t')
	h = h[h$V3=='exon',]
	colnames(h)[c(1,4,5,7)] = c('chr_id','start','stop','strand')
	h$gene_biotype = sapply(strsplit(h$V9,';\\s?'),function(x){x = strsplit(x,' ');setNames(sapply(x,'[',2),sapply(x,'[',1))['gene_biotype']})
	print(table(h$gene_biotype))
	all.anns[[s]]$lncRNA = getAnnOverlap(all.anns[[s]],h[h$gene_biotype=='lncRNA',])
	all.anns[[s]]$pc.lncRNA = getAnnOverlap(all.anns[[s]],h[h$gene_biotype=='putative_coding',])
	all.anns[[s]]$cds = '-'
	all.anns[[s]]$cds[all.anns[[s]]$cds.pos=='cds'] = 'nc-in-cds'
	all.anns[[s]]$cds[all.anns[[s]]$cds.pos=='3utr'] = 'p3utr'
	all.anns[[s]]$cds[all.anns[[s]]$cds.pos=='5utr'] = 'p5utr'
	all.anns[[s]]$cds[all.anns[[s]]$pc.lncRNA != '-'] = 'pc-lncRNA'
	all.anns[[s]]$cds[all.anns[[s]]$lncRNA != '-'] = 'lncRNA'
	all.anns[[s]]$cds[!is.na(all.anns[[s]]$ens.transc.cds.pos.smpl) & all.anns[[s]]$ens.transc.cds.pos.smpl=='3utr'] = '3utr'
	all.anns[[s]]$cds[!is.na(all.anns[[s]]$ens.transc.cds.pos.smpl) & all.anns[[s]]$ens.transc.cds.pos.smpl=='5utr'] = '5utr'
	all.anns[[s]]$cds[all.anns[[s]]$cod != 'n'] = 'cds'
	gc()
}

table(is.na(all.anns$human$antisense.dupl.rate),all.anns$human$cds)

sapply(all.anns,function(x)table(x$cds[x$position=='INTERNAL']))
sapply(all.anns,function(x)table(x$cds[x$position=='INTERNAL' & x$sites=='ad']))
sapply(all.anns,function(x)table(x$cds[x$position=='INTERNAL' & is.na(x$antisense.dupl.rate)]))

sapply(names(anns),function(s)table(all.anns[[s]][rownames(anns[[s]])[anns[[s]]$sites=='ad'],'cds']))

# look on previous seg (to check that cassette exons in 3utrs are frequently follow highly retained intron)
for(s in names(all.anns)){
	print(s)
	ha = all.anns[[s]]
	ha = ha[order(ha$gene_id,ha$start),]
	ha$prev.sid = ha$prev.sites = NA
	i = which(ha$strand==1)
	i = i[i>1]
	i = i[ha$stop[i-1]+1 == ha$start[i]]
	ha$prev.sid[i] = rownames(ha)[i-1]
	
	i = which(ha$strand== -1)
	i = i[i<nrow(ha)]
	i = i[ha$start[i+1]-1 == ha$stop[i]]
	ha$prev.sid[i] = rownames(ha)[i+1]
	
	ha$prev.sites[!is.na(ha$prev.sid)] = ha[ha$prev.sid[!is.na(ha$prev.sid)],'sites']
	
	h = readRDS(paste0('Rdata/',s,'.as.u.all.Rdata'))
	mir = apply(h$ir,1,mean,na.rm=T)
	ha$prev.psi[!is.na(ha$prev.sid)] = mir[ha$prev.sid[!is.na(ha$prev.sid)]]
	all.anns[[s]] = ha[rownames(all.anns[[s]]),]
}
# add is.ce
for(s in names(all.anns)){
	print(s)
	gtf = read.table(paste0('processed/annotation/all.species/merged/',s,'.gtf'),sep='\t')
	gtf = gtf[gtf$V3=='exon',]
	gtf = unique(paste(gtf$V1,gtf$V4,gtf$V5))
	all.anns[[s]]$is.ce = paste(all.anns[[s]]$chr_id,all.anns[[s]]$start,all.anns[[s]]$stop) %in% gtf
}
sapply(all.anns,function(x)table(x$is.ce[x$sites=='ad']))
#saveRDS(all.anns,'Rdata/all.anns.Rdata')

# anns = lapply(names(anns),function(s)all.anns[[s]][rownames(anns[[s]]),])
# names(anns) = names(all.anns)
# saveRDS(anns,'Rdata/anns.Rdata') 

lvs = c('5utr','p5utr','cds','nc-in-cds','3utr','p3utr','lncRNA','pc-lncRNA','-')
sapply(all.anns,function(x)table(factor(x$cds[x$sites=='ad'],levels=lvs)))
b=sapply(all.anns,function(x)table(factor(x$cds[x$sites=='ad' & (!is.na(x$prev.psi) & x$prev.sites == 'da' & x$prev.psi > 0.5)],levels=lvs)))/sapply(all.anns,function(x)table(factor(x$cds[x$sites=='ad'],levels=lvs)))
b=sapply(all.anns,function(x)table(factor(x$cds[x$sites=='ad' & (!is.na(x$prev.psi) & x$prev.psi > 0.5)],levels=lvs)))/sapply(all.anns,function(x)table(factor(x$cds[x$sites=='ad'],levels=lvs)))
barplot(t(b),bes=T,las=3,ylab='% with high prev seg')

anns = lapply(names(anns),function(s)all.anns[[s]][rownames(anns[[s]]),])
names(anns) = names(all.anns)

sapply(anns,function(x)table(factor(x$cds[x$sites=='ad'],levels=lvs)))
b=sapply(anns,function(x)table(factor(x$cds[x$sites=='ad' & (!is.na(x$prev.psi) & x$prev.sites == 'da' & x$prev.psi > 0.5)],levels=lvs)))/sapply(anns,function(x)table(factor(x$cds[x$sites=='ad'],levels=lvs)))
b=sapply(anns,function(x)table(factor(x$cds[x$sites=='ad' & (!is.na(x$prev.psi) & x$prev.psi > 0.5)],levels=lvs)))/sapply(anns,function(x)table(factor(x$cds[x$sites=='ad'],levels=lvs)))
barplot(t(b),bes=T,las=3,ylab='% with high prev seg')


anns$human[anns$human$cds=='nc-in-cds' & anns$human$sites=='ad',][1:10,]
f =  anns$human$sites=='ad'
table(anns$human$cds[f],anns$human$ens.exon.overlap[f])
plotTissueAgeProile(psi.tsm$human['hum.6786.s50',],meta.tsm)

anns$human['hum.6786.s50',]

b = list.files('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/',pattern = 'Testis',full.names = T)
b = b[grep('.bam$',b)]


r = getReadCoverage(b,'10',18932498-2000,18932569+7000,1,T,min.junc.cov = 10)
abline(v=c(18932498,18932569),col='blue')

## exont #####
library(GenomicRanges)
library(ontologyIndex)
h = readRDS('Rdata/all.anns.Rdata')$human
gc()
exont = read.table('processed/exon.onthology/exont_annotations_v1.5.0.tsv',sep='\t',header=T)
colnames(exont) = c('exont.id','chr.id','strand','start','end')
exont.gr = GRanges(exont$chr.id,IRanges(exont$start,exont$end),ifelse(exont$strand==1,'+','-'))
seg.gr = GRanges(h$chr_id,IRanges(h$start,h$stop),ifelse(is.na(h$strand),'*',ifelse(h$strand== 1,'+','-')))
s2e = findOverlaps(exont.gr,seg.gr,maxgap=0,type='any',select='all',ignore.strand=FALSE)
s2e = data.frame(exont.inx=s2e@from,exont.id = exont$exont.id[s2e@from],seg.id = rownames(h)[s2e@to])
table(h$cod,rownames(h) %in% s2e$seg.id)
table(table(s2e$exont.inx))
plot(table(table(s2e$seg.id)),xlim=c(0,50))
seg2exont = sapply(split(s2e$exont.id,s2e$seg.id),unique)
plot(table(sapply(seg2exont,length)))
o = get_OBO('processed/exon.onthology/exont.obo')
o$name[get_ancestors(o,c('EXONT:000064','EXONT:000060'))]
#saveRDS(seg2exont,'Rdata/seg2exont.Rdata')

# orth other types #####

aas = lapply(names(all.anns),function(s){
	print(s)
	#alt.ids = rownames(alts.filt[[s]])[alts.filt[[s]]$sites=='daa' & alts.filt[[s]]$segs=='1-2' & alts.filt[[s]]$ints=='0-1;0-2']
	#sids = rownames(all.anns[[s]])[!is.na(all.anns[[s]]$alt.id) & all.anns[[s]]$alt.id %in% alt.ids & all.anns[[s]]$sites=='aa']
	sids = rownames(all.anns[[s]])[all.anns[[s]]$sites=='aa']
	t = all.anns[[s]][all.anns[[s]]$sites %in% c('ad','ae') & all.anns[[s]]$gene_id %in% all.anns[[s]][sids,'gene_id'],]
	t$seg.id = rownames(t)
	rownames(t) = paste(t$gene_id,ifelse(t$strand==1,t$start-1,t$stop+1))
	aaseg = all.anns[[s]][sids,]
	host = t[paste(all.anns[[s]][sids,'gene_id'],ifelse(aaseg$strand==1,aaseg$stop,aaseg$start)),'seg.id']
	print(table(is.na(host)))
	cbind(aaseg.id=sids[!is.na(host)],host.id=host[!is.na(host)])
})
names(aas) = names(all.anns)
sapply(aas,dim)
sapply(aas,function(x){table(x[,2] %in% orth.seg.ad.all.id)})

aa.orth.host = sapply(1:ncol(orth.seg.ad.all.id),function(s){
	orth.seg.ad.all.id[,s] %in% aas[[s]][,2]
})

table(apply(aa.orth.host,1,sum))
t=sort(table(apply(aa.orth.host,1,function(x){paste(species$short[x],collapse='')})))
t[nchar(names(t))==2]

orth.aas.all = list(host = orth.seg.ad.all.id[apply(aa.orth.host,1,sum)>0,])
orth.aas.all$aa = orth.aas.all$host
orth.aas.all$aa[,] = NA
orth.aas.all$host.len = orth.aas.all$aa.len = array(NA,dim=dim(orth.aas.all$aa),dimnames=dimnames(orth.aas.all$aa))
for(s in colnames(orth.aas.all$host)){
	print(s)
	orth.aas.all$host.len[,s] = all.anns[[s]][orth.aas.all$host[,s],'length']
	rownames(aas[[s]]) = aas[[s]][,2]
	cmn = intersect(aas[[s]][,2],orth.aas.all$host[,s])
	cmn.inx = match(cmn,orth.aas.all$host[,s])
	orth.aas.all$aa[cmn.inx,s] = aas[[s]][cmn,1]
	orth.aas.all$aa.len[cmn.inx,s] = all.anns[[s]][orth.aas.all$aa[cmn.inx,s],'length']
}
class(orth.aas.all) = c('sajr','list')
orth.aas.all$host = cbind(orth.aas.all$host,species = apply(!is.na(orth.aas.all$aa),1,function(x){paste(species$short[x],collapse='')}))
orth.aas.all[orth.aas.all$host[,8] == 'hqmrb',]
orth.aas.all$host = as.data.frame(orth.aas.all$host)
orth.aas.all$host$host.same = apply(orth.aas.all$host.len,1,function(x)length(unique(x))) == 1
orth.aas.all$host$aa.same = apply(orth.aas.all$aa.len,1,function(x)length(unique(x[!is.na(x)]))) == 1
z = orth.aas.all$aa.len
z[is.na(z)] = 0
orth.aas.all$host$host.plus.aa.same = apply(orth.aas.all$host.len+z,1,function(x)length(unique(x))) == 1

sort(table(orth.aas.all$host$species[nchar(orth.aas.all$host$species)==6 & orth.aas.all$host$aa.same & orth.aas.all$host$host.plus.aa.same]))

f = apply(aa.orth.host,1,sum) == 7
orth.aas = sapply(1:7,function(s){
	rownames(aas[[s]]) = aas[[s]][,2]
	aas[[s]][orth.seg.ad.all.id[f,s],1]
})
lens = sapply(1:7,function(s){all.anns[[s]][orth.aas[,s],'length']})
table(apply(lens,1,function(x){length(unique(x))}))
f = apply(lens,1,function(x){length(unique(x))}) == 1
table(lens[f,1]%%3)

colnames(orth.aas) = names(aas)
orth.aas = loadInfoForOrths(orth.aas)
orth.aas$human$seg$same.len = apply(sapply(orth.aas,function(x)x$seg$length),1,function(x)length(unique(x))) == 1
simple.alt = sapply(1:7,function(s){a=alts.filt[[s]][orth.aas[[s]]$seg$alt.id,];a$sites=='daa' & a$segs=='1-2' & a$ints=='0-1;0-2'})
orth.aas$human$seg$n.species.simple.alt = apply(simple.alt,1,sum)
hist(orth.aas$human$seg$length[orth.aas$human$seg$same.len],0:10000,col=c('gray','gray','red'),xlim=c(0,50))
# saveRDS(orth.aas,'Rdata/not.cassette/orth.aas.Rdata')
# saveRDS(orth.aas.all,'Rdata/not.cassette/orth.aas.all.Rdata')
# saveRDS(aas,'Rdata/not.cassette/aas.Rdata')

# alt donors
dds = lapply(names(all.anns),function(s){
	print(s)
	sids = rownames(all.anns[[s]])[all.anns[[s]]$sites=='dd']
	t = all.anns[[s]][all.anns[[s]]$sites %in% c('ad','sd') & all.anns[[s]]$gene_id %in% all.anns[[s]][sids,'gene_id'],]
	t$seg.id = rownames(t)
	rownames(t) = paste(t$gene_id,ifelse(t$strand==1,t$stop+1,t$start-1))
	ddseg = all.anns[[s]][sids,]
	host = t[paste(ddseg$gene_id,ifelse(ddseg$strand==1,ddseg$start,ddseg$stop)),'seg.id']
	print(table(is.na(host)))
	#browser()
	cbind(ddseg.id=sids[!is.na(host)],host.id=host[!is.na(host)])
})

names(dds) = names(all.anns)
sapply(dds,dim)
sapply(dds,function(x){table(x[,2] %in% orth.seg.ad.all.id)})

dd.orth.host = sapply(1:ncol(orth.seg.ad.all.id),function(s){
	orth.seg.ad.all.id[,s] %in% dds[[s]][,2]
})

table(apply(dd.orth.host,1,sum))
t=sort(table(apply(dd.orth.host,1,function(x){paste(species$short[x],collapse='')})))
t[nchar(names(t))==6]


orth.dds.all = list(host = orth.seg.ad.all.id[apply(dd.orth.host,1,sum)>0,])
orth.dds.all$dd = orth.dds.all$host
orth.dds.all$dd[,] = NA
orth.dds.all$host.len = orth.dds.all$dd.len = array(NA,dim=dim(orth.dds.all$dd),dimnames=dimnames(orth.dds.all$dd))
for(s in colnames(orth.dds.all$host)){
	print(s)
	orth.dds.all$host.len[,s] = all.anns[[s]][orth.dds.all$host[,s],'length']
	rownames(dds[[s]]) = dds[[s]][,2]
	cmn = intersect(dds[[s]][,2],orth.dds.all$host[,s])
	cmn.inx = match(cmn,orth.dds.all$host[,s])
	orth.dds.all$dd[cmn.inx,s] = dds[[s]][cmn,1]
	orth.dds.all$dd.len[cmn.inx,s] = all.anns[[s]][orth.dds.all$dd[cmn.inx,s],'length']
}
class(orth.dds.all) = c('sajr','list')
orth.dds.all$host = cbind(orth.dds.all$host,species = apply(!is.na(orth.dds.all$dd),1,function(x){paste(species$short[x],collapse='')}))
orth.dds.all[orth.dds.all$host[,8] == 'hqmrb',]
orth.dds.all$host = as.data.frame(orth.dds.all$host)
orth.dds.all$host$host.same = apply(orth.dds.all$host.len,1,function(x)length(unique(x))) == 1
orth.dds.all$host$dd.same = apply(orth.dds.all$dd.len,1,function(x)length(unique(x[!is.na(x)]))) == 1
z = orth.dds.all$dd.len
z[is.na(z)] = 0
orth.dds.all$host$host.plus.dd.same = apply(orth.dds.all$host.len+z,1,function(x)length(unique(x))) == 1

f = apply(dd.orth.host,1,sum) == 7
orth.dds = sapply(1:7,function(s){
	rownames(dds[[s]]) = dds[[s]][,2]
	dds[[s]][orth.seg.ad.all.id[f,s],1]
})

dd.lens = sapply(1:7,function(s){all.anns[[s]][orth.dds[,s],'length']})
table(apply(dd.lens,1,function(x){length(unique(x))}))
f = apply(dd.lens,1,function(x){length(unique(x))}) == 1
table(dd.lens[f,1])
hist(dd.lens[f,1],0:330,xlim=c(0,50),col=c('gray','gray','red'))

colnames(orth.dds) = names(dds)
orth.dds = loadInfoForOrths(orth.dds)
orth.dds$human$seg$same.len = apply(sapply(orth.dds,function(x)x$seg$length),1,function(x)length(unique(x))) == 1
simple.alt = sapply(1:7,function(s){a=alts.filt[[s]][orth.dds[[s]]$seg$alt.id,];a$sites=='dda' & a$segs=='0-1' & a$ints=='0-2;1-2'})
orth.dds$human$seg$n.species.simple.alt = apply(simple.alt,1,sum)
# saveRDS(orth.dds,'Rdata/not.cassette/orth.dds.Rdata')
# saveRDS(orth.dds.all,'Rdata/not.cassette/orth.dds.all.Rdata')
# saveRDS(dds,'Rdata/not.cassette/dds.Rdata')

# retained introns
irs = lapply(names(all.anns),function(s){
	print(s)
	#alt.ids = rownames(alts.filt[[s]])[alts.filt[[s]]$sites=='daa' & alts.filt[[s]]$segs=='1-2' & alts.filt[[s]]$ints=='0-1;0-2']
	#sids = rownames(all.anns[[s]])[!is.na(all.anns[[s]]$alt.id) & all.anns[[s]]$alt.id %in% alt.ids & all.anns[[s]]$sites=='aa']
	sids = rownames(all.anns[[s]])[all.anns[[s]]$sites=='da']
	e1 = all.anns[[s]][all.anns[[s]]$sites %in% c('ad','sd','ae') & all.anns[[s]]$gene_id %in% all.anns[[s]][sids,'gene_id'],]
	e1$seg.id = rownames(e1)
	e2 = e1
	rownames(e1) = paste(e1$gene_id,ifelse(e1$strand== 1,e1$stop+1,e1$start-1))
	rownames(e2) = paste(e2$gene_id,ifelse(e2$strand==-1,e2$stop+1,e2$start-1))
	irseg = all.anns[[s]][sids,]
	host1 = e1[paste(irseg$gene_id,ifelse(irseg$strand== 1,irseg$start,irseg$stop)),'seg.id']
	host2 = e2[paste(irseg$gene_id,ifelse(irseg$strand==-1,irseg$start,irseg$stop)),'seg.id']
	print(table(is.na(host1),is.na(host2)))
	f = !is.na(host1) & !is.na(host2)
	#browser()
	cbind(irseg.id=sids[f],host1.id=host1[f],host2.id=host2[f])
})

names(irs) = names(all.anns)
sapply(irs,dim)
sapply(irs,function(x){table(x[,2] %in% orth.seg.ad.all.id & x[,3] %in% orth.seg.ad.all.id)})

for(i in 1:length(irs))
	irs[[i]] = cbind(as.data.frame(irs[[i]][,1:3]),orth.inx1 = match(irs[[i]][,2],orth.seg.ad.all.id[,i]),orth.inx2 = match(irs[[i]][,3],orth.seg.ad.all.id[,i]))

ir.orth.inxs = unlist(lapply(irs,function(x){x=x[apply(is.na(x),1,sum)==0,];paste(x$orth.inx1,x$orth.inx2)}))
table(table(ir.orth.inxs))

ir.orth.sp = sapply(1:7,function(s){unique(ir.orth.inxs) %in% paste(irs[[s]]$orth.inx1,irs[[s]]$orth.inx2)})
t=sort(table(apply(ir.orth.sp,1,function(x){paste(species$short[x],collapse='')})))
t[nchar(names(t))==2]


orth.irs = table(ir.orth.inxs)
orth.irs = names(orth.irs)[orth.irs==7]

f = apply(dd.orth.host,1,sum) == 7
orth.irs = sapply(1:7,function(s){
	x = irs[[s]]
	x = x[!is.na(x$orth.inx1) & !is.na(x$orth.inx2),]
	rownames(x) = paste(x$orth.inx1,x$orth.inx2)
	x[orth.irs,1]
})
colnames(orth.irs) = names(irs)

ir.lens = sapply(1:7,function(s){all.anns[[s]][orth.irs[,s],'length']})
pairs(ir.lens,log='xy')

# 
# orth.irs = loadInfoForOrths(orth.irs)
# saveRDS(orth.irs,'Rdata/not.cassette/orth.irs.Rdata')
# saveRDS(irs,'Rdata/not.cassette/irs.all.Rdata')
irs.psi = do.call(cbind,lapply(orth.irs,function(x)x$ir))[,rownames(meta)]
irs.cor = cor(irs.psi,u='p')
dim(irs.cor)
irs.mds = cmdscale(1-irs.cor,k=2)


barplotWithText = function(x,t=x,...){
	b=barplot(x,...)
	text(b,x,t,adj=c(0.5,1.1))
}

### load marg GE #######
Sp = paste0(toupper(substr(rownames(species),1,1)),substr(rownames(species),2,200))
ens.ge.marg = lapply(Sp,function(s)read.table(paste0('processed/GE.from.marg/',s,'RpkmMajorTissuesCor90.Norm.txt'),header = 1,row.names = 1))
names(ens.ge.marg) = rownames(species)
sapply(ens.ge.marg,function(x)table(colnames(x) %in% meta$marg.name))

for(s in names(ens.ge.marg)){
	x = meta[meta$species==s,]
	x = setNames(rownames(x),x$marg.name)
	#print(table(colnames(ens.ge.marg[[s]]) %in% names(x)))
	colnames(ens.ge.marg[[s]]) = x[colnames(ens.ge.marg[[s]])]
}


saveRDS(ens.ge.marg,'Rdata/ens.ge.marg.Rdata')
ens.ge.marg.tsm = lapply(ens.ge.marg,function(g){
	m = meta[colnames(g),]
	calcMeanCols(g,paste(m$species,m$tissue,m$stage),FUN=base::mean)
	})
saveRDS(ens.ge.marg.tsm,'Rdata/ens.ge.marg.tsm.Rdata')



#make all MDS
mds.per.species.tsm = cnts = list()
for(s in rownames(species)){
	print(s)
	#t = readRDS(paste0('Rdata/',s,'.as.u.filtered.Rdata'))
	cnts[[s]] = mds.per.species.tsm[[s]] = list()
	for(sites in c('ad','aa','dd','da')){
		cat('\t',sites)
		#d = t$ir[t$seg$sites==sites,]
		d = psi.tsm[[s]][anns[[s]]$sites == sites,]
		cnts[[s]][[sites]] = nrow(d)
		cat(' ',nrow(d),'\n')
		d = cor(d,u='p',m='pears')
		mds.per.species.tsm[[s]][[sites]] = cmdscale(1-d,k=2)
	}
}
mds.per.species.tsm$counts = cnts
#saveRDS(mds.per.species.tsm,'Rdata/mds.per.species.tsm.pearson.Rdata')

mds.per.species.tsm = readRDS('Rdata/mds.per.species.tsm.pearson.Rdata')
pdf('figures/MDS.be.species-sites.pearson.pdf',w=12,h=21)
n = c(ad='cassette',aa='acceptor',dd='donor',da='int. ret.')
par(mfrow=c(7,4),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(3,2,1.5,0),oma=c(0,0,2,1))
for(s in rownames(species)){
	for(sites in names(mds.per.species.tsm[[s]])){
		m = meta.tsm[rownames(mds.per.species.tsm[[s]][[sites]]),]
		plotMDS(points=mds.per.species.tsm[[s]][[sites]],col=m$col,cex=m$cex*1.5,pch=19,main=paste0(s,', ',n[sites],' (',mds.per.species.tsm$counts[[s]][[sites]][1],')'))
	}
}
dev.off()

orth.irs = readRDS('Rdata/not.cassette/orth.irs.Rdata')
orth.aas = readRDS('Rdata/not.cassette/orth.aas.Rdata')
orth.dds = readRDS('Rdata/not.cassette/orth.dds.Rdata')
orth.aas.all = readRDS('Rdata/not.cassette/orth.aas.all.Rdata')
orth.dds.all = readRDS('Rdata/not.cassette/orth.dds.all.Rdata')

orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
#orth.seg.ad.old = readRDS('old.Rdata/orth.seg.ad.Rdata')

pdf('figures/orth.N-species.stat.pdf',w=14,h=5)
par(mfrow=c(1,4),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(3,2,1.5,0),oma=c(0,0,2,1))
barplotWithText(table(apply(sapply(orth.seg.ad,function(x)x$seg$type)=='ALT',1,sum)),log='y',main='Orthologous cassette exons',xlab='number of species',ylab='number of exons')
barplotWithText(table(apply(!is.na(orth.dds.all$dd),1,sum)),log='y',main='Orthologous alternative donors',xlab='number of species',ylab='number of exons')
barplotWithText(table(apply(!is.na(orth.aas.all$aa),1,sum)),log='y',main='Orthologous alternative acceptors',xlab='number of species',ylab='number of exons')
barplotWithText(table(table(ir.orth.inxs)),log='y',main='Orthologous retained introns',xlab='number of species',ylab='number of exons')
dev.off()


irs.psi = do.call(cbind,lapply(orth.irs,function(x)x$ir))[,rownames(meta)]
aas.psi = do.call(cbind,lapply(orth.aas,function(x)x$ir))[,rownames(meta)]
dds.psi = do.call(cbind,lapply(orth.dds,function(x)x$ir))[,rownames(meta)]
ads.psi = do.call(cbind,lapply(orth.seg.ad,function(x)x$ir))[,rownames(meta)]
# ads.psi.old = do.call(cbind,lapply(orth.seg.ad.old,function(x)x$ir))
# 
# colnames(ads.psi.old) = meta.old[colnames(ads.psi.old),'lib.id']
# ads.psi.old = ads.psi.old[,meta$lib.id]

mds = list()
a = apply(sapply(orth.seg.ad,function(x){apply(x$ir,1,function(y){sum(y>0.2 & y<0.8,na.rm=T)})})>3,1,sum) > 0
f = apply(sapply(orth.seg.ad,function(x)x$seg$type)=='ALT',1,sum)
table(a) # 25602
mds$aa.all = cmdscale(1-cor(aas.psi,u='p',m='pear'),k=2)
mds$dd.all = cmdscale(1-cor(dds.psi,u='p',m='pear'),k=2)
mds$ir.all = cmdscale(1-cor(irs.psi,u='p',m='pear'),k=2)
mds$ad.all = cmdscale(1-cor(ads.psi,u='p',m='pear'),k=2)
mds$ad.alt = cmdscale(1-cor(ads.psi[a,],u='p',m='pear'),k=2)
mds$ad.anc = cmdscale(1-cor(ads.psi[f==7,],u='p',m='pear'),k=2)
#saveRDS(mds,'Rdata/orth.mds.Rdata')

mds = readRDS('Rdata/orth.mds.Rdata')

pdf('figures/orth.MDS.pearson.pdf',w=10,h=10)
par(mfrow=c(3,3),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(2,2,1.5,0),oma=c(0,0,0,1))
m = meta
plotMDS(points=mds$aa.all,col=m$col,cex=m$cex,pch=m$pch,main=paste0('Acceptors (',nrow(aas.psi),')'))
plotMDS(points=mds$dd.all,col=m$col,cex=m$cex,pch=m$pch,main=paste0('Donors (',nrow(dds.psi),')'))
plotMDS(points=mds$ir.all,col=m$col,cex=m$cex,pch=m$pch,main=paste0('Retained introns (',nrow(irs.psi),')'))
plotMDS(points=-mds$ad.all,col=m$col,cex=m$cex,pch=m$pch,main=paste0('Cassette exons (',nrow(ads.psi),')'))
plotMDS(points=-mds$ad.alt,col=m$col,cex=m$cex,pch=m$pch,main=paste0('Cassette ([0.2,0.8] > 3) exons (',sum(a),')'))
plotMDS(points=mds$ad.anc,col=m$col,cex=m$cex,pch=m$pch,main=paste0('Ancient cassette exons (',sum(f==7),')'))
plot.new()
sp = unique(meta[,c('species','pch')])
ti = unique(meta[,c('tissue','col')])
legend('topleft',pch=c(rep(19,7),sp$pch),col=c(ti$col,rep('black',7)),legend=c(ti$tissue,sp$species),ncol=2)
dev.off()

### position of alt aa and dd
sp = c('qmrboc','mrboc','oc','c')
sp = c('m','mr','mrb','hqmrb','hqmrbo')
sp = c('h','hq','hqmrb','hqmrbo')
sp = c('h','q','m','r','b','o','c','hq','mr','mrb','hqmrb','hqmrbo')

dd.stat=sapply(sp,function(ss){
	t = orth.dds.all[orth.dds.all$host$species %in% ss & orth.dds.all$host$dd.same ,]
	ti = apply(t$dd.len[t$host$host.same,,drop=F],1,function(x)x[!is.na(x)][1])
	te = apply(t$dd.len[t$host$host.plus.dd.same,,drop=F],1,function(x)x[!is.na(x)][1])
	c(exn3 = sum(te==3),exn3n=sum(te[te!=3]%%3==0),exn=sum(te%%3!=0),int3 = sum(ti==3),int3n=sum(ti[ti!=3]%%3==0),int=sum(ti%%3!=0))
})

aa.stat=sapply(sp,function(ss){
	t = orth.aas.all[orth.aas.all$host$species %in% ss & orth.aas.all$host$aa.same ,]
	ti = apply(t$aa.len[t$host$host.same,,drop=F],1,function(x)x[!is.na(x)][1])
	te = apply(t$aa.len[t$host$host.plus.aa.same,,drop=F],1,function(x)x[!is.na(x)][1])
	c(exn3 = sum(te==3),exn3n=sum(te[te!=3]%%3==0),exn=sum(te%%3!=0),int3 = sum(ti==3),int3n=sum(ti[ti!=3]%%3==0),int=sum(ti%%3!=0))
})
pdf('figures/birth.of.alt.aa.dd.position.and.length',w=5,h=7)
par(mfrow=c(3,1),tck=-0.02,mgp=c(2.1,0.2,0),mar=c(6,3,1.5,0),oma=c(0,0,0,1))
cols = c('#FF0000FF','#FF0000CC','#FF000066','#0000FFFF','#0000FFCC','#0000FF66')
barplot(sweep(aa.stat,2,apply(aa.stat,2,sum),'/'),col=cols,names.arg = paste0(sp,'\n(',apply(aa.stat,2,sum),')'),las=3,ylab='% of segments',main='Alternative acceptors',xlab='species with site')
barplot(sweep(dd.stat,2,apply(dd.stat,2,sum),'/'),col=cols,names.arg = paste0(sp,'\n(',apply(dd.stat,2,sum),')'),las=3,ylab='% of segments',main='Alternative donors',xlab='species with site')
plot.new()
legend('topleft',fill=cols,legend=c('exon 3 nt','exon 3n nt','exon !3n nt','intron 3 nt','intron 3n nt','intron !3n nt'))
dev.off()
table(orth.dds.all$dd.len[orth.dds.all$host$dd.same,1] %% 3,nchar(orth.dds.all$host$species[orth.dds.all$host$dd.same]))
table(orth.aas.all$aa.len[orth.aas.all$host$aa.same,1] %% 3,nchar(orth.aas.all$host$species[orth.aas.all$host$aa.same]))

plot(table(orth.aas.all$aa.len[orth.aas.all$host$aa.same & nchar(orth.aas.all$host$species)>1,1]),xlim=c(0,50),col=c(1,1,2))
plot(table(orth.dds.all$dd.len[orth.dds.all$host$dd.same & nchar(orth.dds.all$host$species)>1,1]),xlim=c(0,50),col=c(1,1,2))

##### plot alts #############
pdf('figures/alts.pdf',w=10,h=10)
par(mfrow=c(7,7),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(3,2,1.5,0),oma=c(0,0,2,1))
for(i in 1:7){
	plotAllAlts(alts[[i]],to.plot=49)
	mtext(rownames(species)[i],outer=TRUE)
}
dev.off()

pdf('figures/alts.no.ints.pdf',w=10,h=10)
par(mfrow=c(7,7),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(3,2,1.5,0),oma=c(0,0,2,1))
for(i in 1:7){
	plotAllAlts(alts.filt[[i]],to.plot=49)
	mtext(rownames(species)[i],outer=TRUE)
}
dev.off()


### compare old data and new ####

o = readRDS('old.Rdata/mouse.as.u.all.Rdata')
o$i = o$e = NULL
gc()
n = readRDS('Rdata/mouse.as.u.all.Rdata')
n$i = n$e = NULL
gc()
on = paste(o$seg$chr_id,o$seg$strand,o$seg$start,o$seg$stop)
nn = paste(n$seg$chr_id,n$seg$strand,n$seg$start,n$seg$stop)
cmn = intersect(on,nn)
length(cmn)
t = table(o$seg$sites,on %in% nn)
t = table(n$seg$sites,nn %in% on)
t = cbind(t,t[,2]/(t[,1]+t[,2]))
t[order(t[,3]),]

cmn.sam = intersect(colnames(o$ir),colnames(n$ir))
o = o[setNames(rownames(o$seg),on)[cmn],cmn.sam]
n = n[setNames(rownames(n$seg),nn)[cmn],cmn.sam]

f = n$seg$sites=='ad' & n$seg$type=='ALT' & o$seg$type == 'ALT'
table(f)
plotLine(o$ir[f,1],n$ir[f,1],pch='.')

cor = cor(cbind(o$ir[f,],n$ir[f,]),u='p')
image(cor)
z=sapply(1:nrow(cor),function(i)order(-cor[i,359:716])[2])
plot(z[1:ncol(o$ir)])
table(z[1:ncol(o$ir)]-1:ncol(o$ir)-ncol(o$ir)==0)


# MEX #####
checkMEXCor = function(a,p){
	a = a[a$sites=='ad',]
	p = p[rownames(a),]
	o = order(a$gene_id,a$exon.number)
	a = a[o,]
	p = p[o,]
	r = list()
	for(i in 2:nrow(a)){
		cat('\r',i,nrow(a))
		if(a$exon.number[i-1]+1 == a$exon.number[i]){
			f = !is.na(p[i-1,]) & !is.na(p[i,])
			r[[length(r)+1]] = data.frame(sid1 = rownames(a)[i-1],sid2 = rownames(a)[i],n=sum(f),mean.sum=mean(p[i-1,]+p[i,],na.rm=T),sd.sum=sd(p[i-1,]+p[i,],na.rm=T),cor=cor(p[i-1,],p[i,],u='p'))
		}
	}
	r = do.call(rbind,r)
	r$l1 = a[r$sid1,'length']
	r$l2 = a[r$sid2,'length']
	r
}



h = checkMEXCor(anns$mouse,psi.tsm$mouse)
hist(h$cor[(h$l1==h$l2)])

hist(h$cor[(h$l1-h$l2)%%3==0])
hist(h$n)

f = h$n>9
table(f)
plot(h$cor[f],h$mean.sum[f],col=((h$l1-h$l2)[f]%%3 == 0) +1)

plot(h$cor[f],h$sd.sum[f],col=((h$l1-h$l2)[f]%%3 == 0) +1)
table(h$cor< -0.8,abs(h$mean.sum-1)<0.1,len=(h$l1-h$l2)%%3 == 0)


mex = h[!is.na(h$cor) & h$cor< -0.7 & abs(h$mean.sum-1)<0.2 ,]
dim(mex)
table(b=age.segs$mouse[mex$sid1,'brain'],h=age.segs$mouse[mex$sid2,'heart'])
table(b=age.segs$mouse[mex$sid2,'brain'],h=age.segs$mouse[mex$sid1,'heart'])
age.segs$mouse[mex$sid2,]


## look on PSI change conservation##
border.stages = readRDS('Rdata/border.stages.Rdata')



orth.ch = lapply(rownames(species)[-2],function(s){
  getAgeASchanges(orth.seg.ad.tsm,meta.tsm,0.3,border.stages,s)
	})
names(orth.ch) = rownames(species)[-2]

getASChangeCons = function(och,tissues,main.sp,filter){
	r = sapply(och,function(x){apply(x[filter,tissues,drop=F],1,paste,collapse='')})
	r = r[apply(r,1,function(x)sum(grepl('-',x)))==0,]
	st = r[,main.sp]
	f = !grepl('n',st)
	r = t(apply(r,2,function(x){c(table(factor(st[f] == x[f],levels = c(F,T))))}))
	cbind(r,freq=r[,'TRUE']/apply(r,1,sum))
}



pdf('figures/hqmrboc/cons.of.brain.AS.human.dPSI>0.3.pdf',w=10,h=4)
f = function(dirs,f=TRUE,...){
	bu = getASChangeCons(orth.ch,'brain','human',orth.ch$human[,'brain'] %in% dirs & f)
	x = 1:nrow(bu)
	plot(bu[,3],t='l',col=params$tissue.col['brain'],ylim=c(0,1),lwd=3,xaxt='n',ylab='Proportion of conserved',xlab='',...)
	axis(1,x,rownames(bu),las=3)
	for(t in setdiff(unique(meta$tissue),'brain')){
		z = getASChangeCons(orth.ch,'brain','human',orth.ch$human[,'brain'] %in% dirs & orth.ch$human[,t]==orth.ch$human[,'brain']& f)
		
		lines(z[,3],t='l',col=params$tissue.col[t],lwd=1)
		pv = sapply(1:nrow(z),function(i)prop.test(c(bu[i,2],z[i,2]),c(bu[i,1]+bu[i,2],z[i,1]+z[i,2]))$p.value)<0.05
		points(x[pv],z[pv,3],pch=19,col=params$tissue.col[t])
	}
}
par(mfrow=c(1,3),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(4,2,1.5,0),oma=c(0,0,2,1))
f('u',main='human, brain, up')
f('d',main='human, brain, down')
f(c('u','d'),main='human, brain, both')
dev.off()

plotUpDownBar = function(d,...){
	up = apply(d=='u',2,sum)
	dw = apply(d=='d',2,sum)
	barplot(up,ylim=c(-max(dw),max(up)),...)
	barplot(-dw,add=T,...)
}

# microexons
pdf('figures/age.dPSI>0.5.ad.and.microexons.stat.pdf',w=21,h=6)
par(mfcol=c(2,7),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(4,2,1.5,0),oma=c(0,0,2,1))
for(s in rownames(species)){
	r=getAgeASchanges(psi.tsm,meta.tsm,0.5,border.stages,s)
	plotUpDownBar(r[anns[[s]]$sites=='ad',],col=params$tissue.col[colnames(r)],las=3,main=paste0(s,', cassette exons'))
	plotUpDownBar(r[anns[[s]]$sites=='ad' & anns[[s]]$length<=27,],col=params$tissue.col[colnames(r)],las=3,main=paste0(s,', microexons'))
}
dev.off()



f = function(tis,main='',...){
	microexons = getASChangeCons(orth.ch,tis,'human',orth.ch$human[,tis] %in% 'u' & orth.seg.ad$human$seg$length<=27)
	macroexons = getASChangeCons(orth.ch,tis,'human',orth.ch$human[,tis] %in% 'u' & orth.seg.ad$human$seg$length> 27)
	cnts = paste0('(',microexons[1,2],'/',macroexons[1,2],')')
	microexons = apply(microexons,1,function(x){r = prop.test(x[2],sum(x[1:2]));c(r$estimate,r$conf.int)})
	macroexons = apply(macroexons,1,function(x){r = prop.test(x[2],sum(x[1:2]));c(r$estimate,r$conf.int)})
	x = 1:ncol(macroexons)
	plot(x,microexons[1,],t='l',lwd=3,xaxt='n',xlab='',ylab='Proportion of conserved',col='red',ylim=c(0,1),main=paste(main,cnts),...)
	lines(x,macroexons[1,],lwd=3,col='blue')
	segments(x,microexons[2,],x,microexons[3,],col='red')
	segments(x,macroexons[2,],x,macroexons[3,],col='blue')
	axis(1,x,colnames(microexons),las=3)
}

pdf('figures/human.microexons.cons.dPSI>0.3.pdf',w=10,h=10)
par(mfrow=c(3,3),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(4,2,1.5,0),oma=c(0,0,2,1))
for(t in unique(meta$tissue))
	f(t,main=t)
dev.off()


# exons not present in Ensembl
he = read.table('processed/annotation/all.species/ensambl/Homo_sapiens.GRCh37.73.gtf.gz',sep='\t',quote = '')
he = unique(he[he$V3=='exon',1:8])

cmn.chr = union(all.anns$human$chr_id,he$V1)
seg.gr = GRanges(all.anns$human$chr_id,IRanges(all.anns$human$start,all.anns$human$stop),ifelse(is.na(all.anns$human$strand),'*',ifelse(all.anns$human$strand== 1,'+','-')),seqinfo = Seqinfo(cmn.chr ))
ens.gr = GRanges(he$V1,IRanges(he$V4,he$V5),he$V7,seqinfo = Seqinfo(cmn.chr ))

s2e.any = findOverlaps(ens.gr,seg.gr,maxgap=0,type='any',select='all',ignore.strand=FALSE)
s2e.any = cbind(s2e.any@from,s2e.any@to)

s2e.win = findOverlaps(ens.gr,seg.gr,maxgap=0,type='within',select='all',ignore.strand=FALSE)
s2e.win = cbind(s2e.win@from,s2e.win@to)

s2e.eql = findOverlaps(ens.gr,seg.gr,maxgap=0,type='equal',select='all',ignore.strand=FALSE)
s2e.eql = cbind(s2e.eql@from,s2e.eql@to)


dim(s2e.eql)

table(1:nrow(all.anns$human) %in% s2e.any[,2],all.anns$human$cod)
table(1:nrow(all.anns$human) %in% s2e.win[,2],all.anns$human$cod)
table(1:nrow(all.anns$human) %in% s2e.eql[,2],all.anns$human$cod)

all.anns$human$ens.exon.overlap = '-'
all.anns$human$ens.exon.overlap[1:nrow(all.anns$human) %in% s2e.any[,2]] = 'o'
all.anns$human$ens.exon.overlap[1:nrow(all.anns$human) %in% s2e.win[,2]] = 'w'
all.anns$human$ens.exon.overlap[1:nrow(all.anns$human) %in% s2e.eql[,2]] = 'e'

# anns$human = all.anns$human[rownames(anns$human),]
# orth.seg.ad$human$seg$ens.exon.overlap = all.anns$human[rownames(orth.seg.ad$human$seg),'ens.exon.overlap']
# orth.seg.ad.all$human$seg$ens.exon.overlap = all.anns$human[rownames(orth.seg.ad.all$human$seg),'ens.exon.overlap']
# saveRDS(all.anns,'Rdata/all.anns.Rdata')
# saveRDS(anns,'Rdata/anns.Rdata')
# saveRDS(orth.seg.ad,'Rdata/orth.seg.ad.Rdata')
# saveRDS(orth.seg.ad.all,'Rdata/orth.seg.ad.all.Rdata')

# load gene descr ####
library(biomaRt)
descrs = list()
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",version = 75)
descrs$human=getBM(attributes=c('ensembl_gene_id','external_gene_id','description'), mart =ensembl)

ensembl = useEnsembl(biomart="ensembl", dataset="mmulatta_gene_ensembl",version = 84)
descrs$macaque=getBM(attributes=c('ensembl_gene_id','external_gene_name','description'), mart =ensembl)

ensembl = useEnsembl(biomart="ensembl", dataset="rnorvegicus_gene_ensembl",version = 79)
descrs$rat=getBM(attributes=c('ensembl_gene_id','external_gene_name','description'), mart =ensembl)

ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl",version = 84)
descrs$mouse=getBM(attributes=c('ensembl_gene_id','external_gene_name','description'), mart =ensembl)

ensembl = useEnsembl(biomart="ensembl", dataset="ocuniculus_gene_ensembl",version = 84)
descrs$rabbit=getBM(attributes=c('ensembl_gene_id','external_gene_name','description'), mart =ensembl)

ensembl = useEnsembl(biomart="ensembl", dataset="mdomestica_gene_ensembl",version = 84)
descrs$opossum=getBM(attributes=c('ensembl_gene_id','external_gene_name','description'), mart =ensembl)

ensembl = useEnsembl(biomart="ensembl", dataset="ggallus_gene_ensembl",version = 84)
descrs$chicken=getBM(attributes=c('ensembl_gene_id','external_gene_name','description'), mart =ensembl)



lapply(descrs,head)
sapply(rownames(species),function(s)table(rownames(ens.ge[[s]]$gene) %in% descrs[[s]][,1])) # wrong for human because I used ens73

descrs$human = unique(read.table('input/hs.37.73.gene.descr.txt',sep=',',quote='"',header=TRUE))

for(s in names(descrs)){
	rownames(descrs[[s]]) = descrs[[s]][,1]
	descrs[[s]] = descrs[[s]][,-1]
	colnames(descrs[[s]]) = c('gene.name','descr')
	descrs[[s]]$descr =  sapply(strsplit(descrs[[s]]$descr,' [',TRUE),'[',1)
}
#saveRDS(descrs,'Rdata/ens.gene.descr.Rdata',version=2)

# load map stat #####
dim(meta)
map.stat = sapply(paste0('processed/mapping/hisat2.s/',meta$species,'/',meta$fname,'.log'),function(f){print(f);parseHisat2LogS(f)})
map.stat = as.data.frame(t(map.stat))
f = factor(meta$species,levels=rownames(species))
boxplot(map.stat$total ~ f)
boxplot(map.stat$once/map.stat$total ~ f)
