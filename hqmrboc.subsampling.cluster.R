library(png)
library(SAJR)
library(doMC)
library(plyr)
source('code/r.functions/paper.figures.5.F.R')
source('code/r.functions/load.all.data.F.R')
source('../../r.code/util.R')

species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
age.al.i = read.table('input/margarida.age.alignment.final.2017.05.17.improved.tab',sep='\t',header = T)[,c(rownames(species),'to.remove')]

# age.al.i. = age.al.i
# age.al.i.$opossum[c(3,6)] = ''
# age.al.i.$chicken[c(8,9)] = ''
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
# saveRDS(meta,'Rdata/main.set.meta.Rdata')

sp2test. = c(list(c('human','macaque','mouse','rat'),c('human','macaque'),c('human','macaque','chicken'),c('human','chicken'),c('human','mouse','rat','rabbit','opossum'),rownames(species)),as.list(rownames(species)))

species. = c()
sp2test = list()
for(i in 1:length(sp2test.)){
	species. = c(species.,sp2test.[[i]])
	sp2test = c(sp2test,rep(sp2test.[i],length(sp2test.[[i]])))
}


m = meta[!is.na(meta$mouse.stage),]

i = as.numeric(Sys.getenv('PBS_ARRAYID'))

ss = sp2test[[i]]
s = species.[i]

registerDoMC(16)
dh = readRDS(paste0('Rdata/',s,'.as.u.all.Rdata'))
dh = dh[dh$seg$chr_id !='MT' & dh$seg$position=='INTERNAL' & dh$seg$type!='EXN',colnames(dh$ir) %in% rownames(m)]
r = bootstrapDevAS(m,dh,ss,N=100)
saveRDS(r,paste0('Rdata/hqmrboc.subsample/chicken.12-0dpb.20201225/devAS.subsampling-',paste(ss,collapse = '.'),'-',s,'.Rdata'))



#######
# i = 7
# (ss = sp2test[[i]])
# (s = species.[i])
# dh = readRDS(paste0('Rdata/',s,'.as.u.all.Rdata'))
# dh = dh[dh$seg$chr_id !='MT' & dh$seg$position=='INTERNAL' & dh$seg$type!='EXN',colnames(dh$ir) %in% rownames(m)]
# r = bootstrapDevAS(m,dh,ss,N=1,only.samples = T)
# hq = r[[1]]
# hqc = r[[1]]
# table(hq$tissue,hq$mouse.stage,hq$species)-table(hqc$tissue,hqc$mouse.stage,hqc$species)

# check orths #############

