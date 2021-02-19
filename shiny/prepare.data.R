library(SAJR)
library(png)
options(stringsAsFactors = FALSE)
source('code/r.functions/load.all.data.F.R')
source('~/skoltech/r.code/util.R')
source('code/r.functions/load.all.data.F.R')
meta = readRDS('Rdata/main.set.meta.Rdata')
species = readRDS('Rdata/species.Rdata')
all.anns = readRDS('Rdata/all.anns.Rdata')
anns = readRDS('Rdata/anns.Rdata')
orth.ads.all.sp = readRDS('output/shiny/Rdata/2020/orth.segs.all.types-species.Rdata')
psi.tsm.shiny = readRDS('output/shiny/Rdata/2020/psi.tsm.types-species.Rdata')
seg2ens = readRDS('Rdata/seg2ens.Rdata')
age.dpsi = readRDS('Rdata/age.diam.spline4.with.replicates.Rdata')
age.segs = readRDS('Rdata/devAS.4patt.Rdata')
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')
# 2020 #####
# _orth ######
# for 7 species I added 3376 exons by "synteny", for other phylo groups I'll add only liftovered
# orth.alt.ids = readRDS('Rdata/orth.seg.ad.all.id.Rdata')
# 
# phylogr = c('hq','mr','mrb','hqmrb','hqmrbo','hmo')
# orth.ads = lapply(phylogr,function(s){
# 	print(s)
# 	loadAltOrthSegs(c(paste0('processed/orth.segs/only.ad/',s,'.0.6.orth.ad.segs'),paste0('processed/orth.segs/only.ad/',s,'.0.0.orth.ad.segs')),only.filtered = FALSE)
# })
# names(orth.ads) = phylogr
# orth.ads$hqmrboc=orth.alt.ids
# orth.ads = orth.ads[c('hqmrboc','hqmrbo','hqmrb','mrb','mr','hq','hmo')]
# orth.adsf = orth.ads
# 
# total = orth.adsf$hqmrboc
# for(i in 2:length(orth.adsf)){
# 	t = matrix(NA,ncol=7,nrow=nrow(orth.adsf[[i]]))
# 	colnames(t) = rownames(species)
# 	t[,colnames(orth.adsf[[i]])] = orth.adsf[[i]]
# 	f = is.na(t[,1]) | !(t[,1] %in% total[,1])
# 	for(j in 2:ncol(total)){
# 		f = f & (is.na(t[,j]) | !(t[,j] %in% total[,j]))
# 	}
# 	total = rbind(total,t[f,])
# 	orth.adsf[[i]] = t[f,]
# }
# sapply(orth.adsf,dim)
# orth.ads.all.sp = do.call(rbind,orth.adsf)
# dim(orth.ads.all.sp)
# apply(!is.na(orth.ads.all.sp), 2, sum)
# table(table(orth.ads.all.sp))
# # add orth not ad
# orth.irs = readRDS('Rdata/not.cassette/orth.irs.Rdata')
# orth.aas = readRDS('Rdata/not.cassette/orth.aas.Rdata')
# orth.dds = readRDS('Rdata/not.cassette/orth.dds.Rdata')
# sapply(orth.dds,length)
# 
# orth.ads.all.sp = rbind(orth.ads.all.sp,
# 			sapply(orth.irs,function(x)rownames(x$seg)),
# 			sapply(orth.aas,function(x)rownames(x$seg)),
# 			sapply(orth.dds,function(x)rownames(x$seg)))
# #saveRDS(orth.ads.all.sp,'Rdata/orth.segs.all.types.all.phylo.groups.Rdata')
# # add remaining and remove these, that didn't pass filters in all species
# sapply(names(anns),function(s){table(rownames(anns[[s]]) %in% orth.ads.all.sp[,s])})
# sapply(names(anns),function(s){table(rownames(anns[[s]])[anns[[s]]$sites=='ad'] %in% orth.ads.all.sp[,s])})
# sapply(names(anns),function(s){table(orth.ads.all.sp[,s] %in% rownames(anns[[s]]))})
# 
# f = rep(F,nrow(orth.ads.all.sp))
# for(s in colnames(orth.ads.all.sp))
# 	f = f | orth.ads.all.sp[,s] %in% rownames(anns[[s]])
# orth.ads.all.sp = orth.ads.all.sp[f,]
# 
# for(s in colnames(orth.ads.all.sp)){
# 	n = setdiff(rownames(anns[[s]]),orth.ads.all.sp[,s])
# 	nm = matrix(NA,ncol=ncol(orth.ads.all.sp),nrow=length(n))
# 	colnames(nm) = colnames(orth.ads.all.sp)
# 	nm[,s] = n
# 	orth.ads.all.sp = rbind(orth.ads.all.sp,nm)
# }
# dim(orth.ads.all.sp)
#saveRDS(orth.ads.all.sp,'output/shiny/Rdata/2020/orth.segs.all.types-species.Rdata')

# prepare PSI,i,e and PSI.tsm for them
orth.ads.all.sp = readRDS('output/shiny/Rdata/2020/orth.segs.all.types-species.Rdata')
psi.tsm.shiny = list()
psi.all.shiny = list()
for(s in rownames(species)){
	print(s)
	d = readRDS(paste0('Rdata/',s,'.as.u.all.Rdata'))
	sids = orth.ads.all.sp[,s]
	sids = sids[!is.na(sids)]
	d = d[sids,colnames(d$ir) %in% rownames(meta)]
	psi.all.shiny[[s]] = d$ir
	m = meta[colnames(d$ir),]
	#psi.tsm.shiny[[s]] = calcMeanCols(d$ir,paste(m$species,m$tissue,m$stage))
	#write.csv(psi.tsm.shiny[[s]], paste0('~/tmp/shiny/downloads/',s,'.tissue-stage.psi'))
	
	colnames(d$ir) = colnames(d$i) = colnames(d$e) = paste0(first2Upper(m$species),'.',m$marg.name)
	write.csv(d$ir,paste0('~/tmp/shiny/downloads/',s,'.psi'))
	write.csv(d$i, paste0('~/tmp/shiny/downloads/',s,'.e'))
	write.csv(d$e, paste0('~/tmp/shiny/downloads/',s,'.i'))
}
# saveRDS(psi.tsm.shiny,'output/shiny/Rdata/2020/psi.tsm.types-species.Rdata')
# saveRDS(psi.all.shiny,'output/shiny/Rdata/2020/psi.all.types-species.Rdata')

# prepare
# all.anns = readRDS('~/skoltech/projects/evo.devo/Rdata/all.anns.Rdata')
# for(s in names(all.anns)){
# 	sgn = per.tissue.age.qv[[s]] <= 0.05 & age.dpsi[[s]] > 0.2
# 	sgn[is.na(sgn)] = FALSE
# 	dpsi = age.dpsi[[s]]
# 	#dpsi[!sgn] = 0
# 	dpsi[is.na(dpsi)] = 0
# 	
# 	all.anns[[s]]$alt = rownames(all.anns[[s]]) %in% orth.ids[,s]
# 	all.anns[[s]]$filtered = rownames(all.anns[[s]]) %in% rownames(sgn)
# 	all.anns[[s]]$devAS = rownames(all.anns[[s]]) %in% rownames(sgn)[apply(sgn,1,sum)>0]
# 	all.anns[[s]]$max.dPSI = 0
# 	all.anns[[s]][rownames(dpsi),'max.dPSI'] = apply(dpsi,1,max)
# 	
# 	o = orth.ids[!is.na(orth.ids[,s]),]
# 	rownames(o) = o[,s]
# 	o = apply(!is.na(o),1,function(x)paste(species[colnames(o)[x],'short'],collapse=''))
# 	all.anns[[s]]$orth = species[s,'short']
# 	all.anns[[s]][names(o),'orth'] = o
# }
# saveRDS(all.anns,'output/shiny/Rdata/all.anns.shiny.Rdata')

# 
# gene2seg = list()
# for(s in names(seg2ens)){
# 	print(s)
# 	t = seg2ens[[s]]
# 	t = t[names(t) %in% orth.ads.all.sp[,s]]
# 	t = revList(t)
# 	names(t) = tolower(names(t))
# 	gene2seg[[s]] = t
# }
# 
# saveRDS(gene2seg,'~/skoltech/projects/evo.devo/output/shiny/Rdata/gene2seg.Rdata')

# load introns
# introns = lapply(rownames(species),function(s){loadIntrons(paste0('processed/annotation/all.species/merged/',s,'.sajr'))})
# names(introns) = rownames(species)
# saveRDS(introns,'~/skoltech/projects/evo.devo/output/shiny/Rdata/introns.Rdata')
# saveRDS(introns,'~/tmp/shiny/Rdata/introns.Rdata')

# 2019 #####
# load data #####
# original data
meta = readRDS('Rdata/main.set.meta.Rdata')
species = readRDS('Rdata/species.Rdata')
all.anns = readRDS('Rdata/all.anns.Rdata')
anns = readRDS('Rdata/anns.Rdata')
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')
border.stages = readRDS('Rdata/border.stages.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')

# shiny data
alt.anns = readRDS('output/shiny/Rdata/alt.anns.Rdata')
alt.tsm = readRDS('output/shiny/Rdata/alt.tsm.Rdata')
orth.alt.ids = readRDS('output/shiny/Rdata/orth.alt.ids.Rdata')

# make orth ids #####
orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
orth.irs = readRDS('Rdata/not.cassette/orth.irs.Rdata')
orth.aas = readRDS('Rdata/not.cassette/orth.aas.Rdata')
orth.dds = readRDS('Rdata/not.cassette/orth.dds.Rdata')

orth.alt.ids = rbind(cbind(sapply(orth.seg.ad,function(x)rownames(x$seg)),sites='ad'),
			cbind(sapply(orth.irs,function(x)rownames(x$seg)),sites='da'),
			cbind(sapply(orth.aas,function(x)rownames(x$seg)),sites='aa'),
			cbind(sapply(orth.dds,function(x)rownames(x$seg)),sites='dd'))
dim(orth.alt.ids)
table(orth.alt.ids[,8])
# saveRDS(orth.alt.ids,'output/shiny/Rdata/orth.alt.ids.Rdata')

# make alt anns ####
alt.anns = lapply(all.anns,function(x)x[x$type != 'EXN' & x$position == 'INTERNAL',])
sapply(all.anns,nrow)
sapply(alt.anns,nrow)

age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(alt.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)
age.dpsi$macaque = cbind(age.dpsi$macaque[,1:5],ovary=NaN,age.dpsi$macaque[,6,drop=FALSE])

cols = c('gene_id','chr_id','start','stop','strand','sites','cod')
alt.anns = lapply(alt.anns,function(x)x[,cols])
for(s in names(alt.anns)){
	print(s)
	alt.anns[[s]]$used = rownames(alt.anns[[s]]) %in% rownames(anns[[s]])
	orth = setNames(1:nrow(orth.alt.ids),orth.alt.ids[,s])
	alt.anns[[s]]$orth.inx = NA
	cmn = intersect(rownames(alt.anns[[s]]),orth.alt.ids[,s])
	alt.anns[[s]][cmn,'orth.inx'] = orth[cmn]
	sgn = per.tissue.age.qv[[s]] < 0.05
	sgn[is.na(sgn)] = FALSE
	dpsi = age.dpsi[[s]][rownames(sgn),]
	dpsi[!sgn] = 0
	sgn = sgn & abs(dpsi) > 0.2
	sgn[is.na(sgn)] = FALSE
	age.dpsi[[s]][,] = 0
	age.dpsi[[s]][rownames(sgn),] = dpsi
	
	t = substr(colnames(sgn),1,1)
	sgn = apply(sgn,1,function(x)paste(t[x],collapse=''))
	alt.anns[[s]]$sgn.tissues = ''
	alt.anns[[s]][names(sgn),'sgn.tissues'] = sgn
	dpsi = apply(abs(age.dpsi[[s]]),1,max,na.rm=T)
	alt.anns[[s]]$max.dpsi = dpsi
	alt.anns[[s]]$dir = apply(age.dpsi[[s]],1,function(x)sign(x[order(abs(x),decreasing = T)[1]]))
}

# saveRDS(alt.anns,'output/shiny/Rdata/alt.anns.Rdata')

# make tsm ####
alt.tsm = list()
for(s in rownames(species)){
	print(s)
	t = readRDS(paste0('Rdata/',s,'.as.u.all.Rdata'))
	t = t$ir[union(rownames(alt.anns[[s]]),orth.alt.ids[,s]),colnames(t$ir) %in% rownames(meta)]
	m = meta[colnames(t),]
	alt.tsm[[s]] = calcMeanCols(t,paste(m$species,m$tissue,m$stage))
	gc()
}
#saveRDS(alt.tsm,'output/shiny/Rdata/alt.tsm.Rdata')

# take only human 10nt chr
osid10 = orth.alt.ids[orth.alt.ids[,'human'] %in% rownames(alt.anns$human)[alt.anns$human$chr_id=='10'],]
alt.tsm.10 = alt.tsm
alt.tsm.10$human = alt.tsm.10$human[rownames(alt.anns$human)[alt.anns$human$chr_id=='10'],]
for(s in rownames(species)[-1])
	alt.tsm.10[[s]] = alt.tsm.10[[s]][osid10[,s],]
#saveRDS(alt.tsm.10,'output/shiny/Rdata/alt.tsm.10.Rdata')

# load gene description ####
seg2ens = readRDS('Rdata/seg2ens.Rdata')
gtfs = list.files('processed/annotation/all.species/ensambl/','gtf.gz',full.names = T)
gtfs=lapply(gtfs,function(f){
	print(f)
	t = read.table(f,sep='\t',comment.char = '#')
	f = t[,3] == 'gene'
	if(sum(f)>0)
		t = t[f,]
	t = strsplit(gsub('"','',t[,9]),'; ',fixed = TRUE)
	t = unique(as.data.frame(do.call(rbind,lapply(t,function(x){
		x=strsplit(x,' ',TRUE)
		setNames(sapply(x,'[',2),sapply(x,'[',1))[c('gene_id','gene_name')]
	}))))
	})
sapply(gtfs,dim)
names(gtfs) = c('chicken','human','macaque','opossum','mouse','rabbit','rat')

# gene2seg = list()
# for(s in names(seg2ens)){
# 	print(s)
# 	t = seg2ens[[s]]
# 	t = t[names(t) %in% orth.ads.all.sp[,s]]
# 	t = revList(t)
# 	names(t) = tolower(names(t))
# 	gene2seg[[s]] = t
# }
# 
# saveRDS(gene2seg,'output/shiny/Rdata/gene2seg.Rdata')
# saveRDS(gene2seg,'~/tmp/shiny/Rdata/gene2seg.Rdata')

# deploy ####
# install.packages('rsconnect')
# rsconnect::setAccountInfo(name='iaaaka',
# 													token='192C88AA096E4BE18FE19F8BF1D89482',
# 													secret='WScwdrDRxrcgPiXEJxypV033b0QNN+R8Kpy868Zv')
library(rsconnect)
rsconnect::deployApp('output/shiny/')
