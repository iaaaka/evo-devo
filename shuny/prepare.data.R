library(SAJR)
options(stringsAsFactors = FALSE)
source('code/r.functions/load.all.data.F.R')
source('~/skoltech/r.code/util.R')
source('code/r.functions/load.all.data.F.R')

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

gene2seg = list()
for(s in names(seg2ens)){
	print(s)
	t = seg2ens[[s]]
	t = t[names(t) %in% rownames(alt.anns[[s]])]
	t = revList(t)
	esn2name = setNames(gtfs[[s]][,2],gtfs[[s]][,1])
	t = c(t,setNames(t,esn2name[names(t)]))
	names(t) = tolower(names(t))
	gene2seg[[s]] = t
}

saveRDS(gene2seg,'output/shiny/Rdata/gene2seg.Rdata')


# deploy ####
# install.packages('rsconnect')
# rsconnect::setAccountInfo(name='iaaaka',
# 													token='192C88AA096E4BE18FE19F8BF1D89482',
# 													secret='WScwdrDRxrcgPiXEJxypV033b0QNN+R8Kpy868Zv')
library(rsconnect)
rsconnect::deployApp('output/shiny/')
