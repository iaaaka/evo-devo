options(stringsAsFactors = FALSE)
source('code/r.functions/load.all.data.F.R')
source('~/skoltech/r.code/util.R')
library(SAJR)
library(seqinr)
library(reshape)


species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
#all.anns = readRDS('Rdata/all.anns.Rdata')
my.ge.cod = readRDS('Rdata/my.ge.cod.Rdata')
exon.birth.one = readRDS('Rdata/exon.birth.one.Rdata')
#psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
born.exn.sajr = readRDS('Rdata/born.exn.sajr.Rdata')
all.anns. = readRDS('Rdata/anns.subset.for.newborn.cov.plot.Rdata') # see newborn.exons.plot.cov.R
load('Rdata/tmp.exon.birth.Rdata')

# new exon birth #######
orth.seg.ad.all.ids = readRDS('Rdata/orth.seg.ad.all.Rdata')
orth.ad.all.hum.ann = orth.seg.ad.all.ids$human$seg
orth.seg.ad.all.ids = sapply(orth.seg.ad.all.ids,function(x){rownames(x$seg)})
gc()
#exon.birthl = lapply(rownames(species),function(s){findOrthExonsByNeighbors(all.anns,orth.seg.ad.all.ids,s,returnBad = T)})
#exon.birth = unlist(exon.birthl,recursive = F)
#saveRDS(exon.birth,'Rdata/exon.birth.Rdata')
exon.birth = readRDS('Rdata/exon.birth.Rdata')
max.exn.cnt = sapply(exon.birth,function(x){max(sapply(x,nrow))})
plot(table(max.exn.cnt),xlab='Maximal numer of not-orthologous exons',ylab='Frequency')
min.exn.cnt = sapply(exon.birth,function(x){min(sapply(x,nrow))})
plot(min.exn.cnt,max.exn.cnt)
table(min.exn.cnt,max.exn.cnt)[1:6,1:6]

exon.birth.12 = exon.birth[max.exn.cnt==2 & min.exn.cnt==1]
one.eq = sapply(exon.birth.12,function(x){
	lens = lapply(x,'[[','length')
	sum(sapply(unique(unlist(lens)),function(l){
		sum(!sapply(lens,function(ls){l %in% ls})) == 0
	}))
})
table(one.eq)

exon.birth.one = exon.birth[max.exn.cnt==1]
cnames = c(colnames(exon.birth.one[[1]]$macaque),'seg_id')
exon.birth.one = lapply(exon.birth.one,function(x){x=lapply(x,function(x){x$seg_id = rownames(x);x[,cnames]});x=do.call(rbind,x)[rownames(species),];rownames(x)=rownames(species);x})
t1 = sapply(exon.birth.one,function(x){sum(!is.na(x$seg_id))})
t2 = sapply(exon.birth.one,function(x){paste(x$seg_id,collapse=',')})
t1 = data.frame(n=t2,c=t1)
t1 = unique(t1)
t1 = setNames(t1$c,t1$n)
t2 = table(t2)
table(t2[names(t1)] == t1)
length(t2)
for(i in 1:length(exon.birth.one)){
	id = paste(exon.birth.one[[i]]$seg_id,collapse=',')
	if(t1[id] != -1){
		t1[id] = -1
	}else
		exon.birth.one[i] = list(NULL)
}
table(sapply(exon.birth.one,is.null))
exon.birth.one = exon.birth.one[!sapply(exon.birth.one,is.null)]

for(i in 1:length(exon.birth.one)){
	cat('\r',i)
	sids = exon.birth.one[[i]]$seg_id
	sp = which(!is.na(sids))[1]
	seg = all.anns[[sp]][sids[sp],]
	g = all.anns[[sp]][all.anns[[sp]]$gene_id == seg$gene_id,]
	g = g[g$sites == 'ad',]
	g = g[order(g$exon.number),]
	inx = which(rownames(g) == sids[sp])
	sid1 = rownames(g)[inx-1]
	sid2 = rownames(g)[inx+1]
	exon.birth.one[[i]]$useg_id = orth.seg.ad.all.ids[orth.seg.ad.all.ids[,sp] == sid1,]
	exon.birth.one[[i]]$dseg_id = orth.seg.ad.all.ids[orth.seg.ad.all.ids[,sp] == sid2,]
	orth.sids = exon.birth.one[[i]]$useg_id
	orth.sids = sapply(strsplit(orth.sids,'.s',TRUE),'[',1)
	if(sum(!is.na(exon.birth.one[[i]]$gene_id) & exon.birth.one[[i]]$gene_id != orth.sids)>0) stop("ERROR")
	exon.birth.one[[i]]$gene_id = orth.sids
}
names(exon.birth.one) = paste0('t',1:length(exon.birth.one))
exon.num.diff = sapply(rownames(species),function(s){
	t = do.call(rbind,lapply(exon.birth.one,function(x)x[s,c('useg_id','dseg_id')]))
	all.anns[[s]][t$dseg_id,'exon.number']-all.anns[[s]][t$useg_id,'exon.number']
})
table(apply(exon.num.diff<1,1,sum)) # so, there are two cases when something is wrong
exon.birth.one = exon.birth.one[apply(exon.num.diff<1,1,sum)==0]


#born.exn.sajr = loadInfoForOrths(born.seg.ids)
#saveRDS(born.exn.sajr,'Rdata/born.exn.sajr.Rdata')
#saveRDS(exon.birth.one,'Rdata/exon.birth.one.Rdata')
exon.birth.one = readRDS('Rdata/exon.birth.one.Rdata')
born.exn.sajr = readRDS('Rdata/born.exn.sajr.Rdata')
#be.tissue.mean = lapply(born.exn.sajr,function(x){t = x$ir[,colnames(x$ir) %in% rownames(meta)];calcMeanCols(t,meta[colnames(t),'tissue'])})
sp.birth = sapply(exon.birth.one,function(x){paste(species$short[!is.na(x$seg_id)],collapse='')})
born.seg.ids = t(sapply(exon.birth.one,function(x){x$seg_id}))
colnames(born.seg.ids) = rownames(species)
rownames(born.seg.ids) = NULL
born.seg.ids=cbind(as.data.frame(born.seg.ids),species=sp.birth)

# look on patterns in human
pdf('figures/born.exon.PSI.on.evol.ages.pdf',w=12,h=5)
par(mfcol=c(2,5),tck=-0.01,mgp=c(1.1,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,2,1))
m=plotExpAndPsiForNewExons('mouse',c('m','mr','mrb','hqmrb','hqmrbo'),born.exn.sajr,born.seg.ids,my.ge.cod,center = T,scale = T,max.na.prop = 1,ylab.fun.name = 'z-score')
mtext('mouse',3,outer = TRUE)
m=plotExpAndPsiForNewExons('human',c('h','hq','hqmrb','hqmrbo'),born.exn.sajr,born.seg.ids,my.ge.cod,center = T,scale = T,max.na.prop = 1,ylab.fun.name = 'z-score')
mtext('human',3,outer = TRUE)

par(mfcol=c(2,5),tck=-0.01,mgp=c(1.1,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,2,1))
m=plotExpAndPsiForNewExons('mouse',c('m','mr','mrb','hqmrb','hqmrbo'),born.exn.sajr,born.seg.ids,my.ge.cod,center = F,scale = F,max.na.prop = 1,ylab.fun.name = '')
mtext('mouse',3,outer = TRUE)
m=plotExpAndPsiForNewExons('human',c('h','hq','hqmrb','hqmrbo'),born.exn.sajr,born.seg.ids,my.ge.cod,center = F,scale = F,max.na.prop = 1,ylab.fun.name = '')
mtext('human',3,outer = TRUE)
dev.off()


# exon duplication ######
getConseqExonLength = function(a){
	a = a[a$sites=='ad',]
	a = split(a,a$gene_id)
	do.call(rbind,lapply(a,function(g){
		g = g[order(g$exon.number),]
		l = nrow(g)
		data.frame(types = paste0(substr(g$type[-l],1,1),substr(g$type[-1],1,1)),e1=g$exon.number[-l],e2=g$exon.number[-1],l1=g$length[-l],l2=g$length[-1])
		}))
}

m = getConseqExonLength(all.anns$mouse)
f = (m$e2-m$e1) == 1
ff = f & m$l1==m$l2
z=table(m$types[f],ff[f])
z[,2]/(z[,1]+z[,2])

hist(m$l1[f]-m$l2[f],-25000:25000,xlim=c(-200,200))
hist(m$l1[f]-m$l2[f],-25000:25000,xlim=c(-20,20))
hist(m$l1[f]-sample(m$l2[f]),-25000:25000,xlim=c(-20,20),add=T,col='gray')
table(m$l1[f]-m$l2[f] == 0)
table(m$l1[f]-sample(m$l2[f]) == 0)

exon.birth.one[[1]]

sp.birth = sapply(exon.birth.one,function(x){paste(species$short[!is.na(x$seg_id)],collapse='')})

getNewBornLen = function(nb,a,sp){
	a = a[[sp]]
	t(sapply(nb,function(x){
		a[as.character(x[sp,c('useg_id','seg_id','dseg_id')]),'length']
		}))
}


t = getNewBornLen(exon.birth.one[sp.birth=='hqmrbo'],all.anns,'mouse')
t
boxplot(t,outline=F)
pairs(t,log='')
cor(t,m='sp')
cor.test(t[,3],t[,1],m='sp')
table((t[,1]-t[,2])%%3)


# compare nt seqs #####
fas = c(chicken='../index/Gallus_gallus.Galgal4.dna.toplevel.cleanNames.fa',
				human='../index/Homo_sapiens.GRCh37.73.dna.primary_assembly.cleanNames.fa',
				macaque='../index/Macaca_mulatta.MMUL_1.dna.toplevel.cleanNames.fa',
				opossum='../index/Monodelphis_domestica.BROADO5.dna.toplevel.cleanNames.fa',
				mouse='../index/Mus_musculus.GRCm38.dna.primary_assembly.cleanNames.fa',
				rabbit='../index/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.cleanNames.fa',
				rat='../index/Rattus_norvegicus.Rnor_5.0.dna.toplevel.cleanNames.fa')



length(exon.birth.one)

dir.create('processed/new.born.exons')
mar.len=4 # I'll add 4 nt intron margins for exons
for(s in rownames(species)){
	print(s)
	# make gtfs for whole region
	sids = do.call(rbind,lapply(exon.birth.one,function(x)x[s,c('useg_id','dseg_id')]))
	u = all.anns.[[s]][sids$useg_id,]
	d = all.anns.[[s]][sids$dseg_id,]
	t = u$strand==1
	#I'll include adjasted exons into database becasue in some cases these "newborn" exons are actually part of adjasted ones
	write.table(cbind(u$chr_id,'.','seq',ifelse(t,u$start,d$start),ifelse(t,d$stop,u$stop),0,ifelse(t,'+','-'),'.',paste0(s,'.iei.',rownames(sids))),sep='\t',row.names = F,col.names = F,quote = F,file=paste0('processed/new.born.exons/',s,'.iei.gtf'))
	# make gtfs for exons
	sids = sapply(exon.birth.one,function(x)x[s,c('seg_id')])
	sids = sids[!is.na(sids)]
	i = all.anns.[[s]][sids,]
	write.table(cbind(i$chr_id,'.','seq',i$start-mar.len,i$stop+mar.len,0,ifelse(i$strand==1,'+','-'),'.',paste0(s,'.e.',names(sids))),sep='\t',row.names = F,col.names = F,quote = F,file=paste0('processed/new.born.exons/',s,'.e.gtf'))
}
writeLines(c(paste0('python ~/projects/evo.devo/code/not.r/extractSeqFromFasta.py < ',names(fas),'.iei.gtf  ',fas,' fa > ',names(fas),'.iei.fa'),
						 paste0('python ~/projects/evo.devo/code/not.r/extractSeqFromFasta.py < ',names(fas),'.e.gtf  '  ,fas,' fa > ',names(fas),'.e.fa')),con = 'processed/new.born.exons/runExtractSeqFromFasta.sh')
# run on cluster

# _load-n-blast ####
# cat *.iei.fa > all.iei.fa
# cat *.e.fa > all.e.fa
# makeblastdb -in all.iei.fa -dbtype nucl
# blastn -db all.iei.fa -word_size 8 -evalue 10000 -num_threads 16 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qcovs qlen' -query all.e.fa  -max_target_seqs 100000 -dust no -soft_masking false | perl -e 'while($l=<>){@t = split("\t",$l);@a=split("[.]",$t[0]);@b=split("[.]",$t[1]);if($a[2] eq $b[2]){print($l)}}' > same.exon.out.final.tab
e = read.fasta('processed/new.born.exons/all.e.fa',as.string = T)
iei = read.fasta('processed/new.born.exons/all.iei.fa',as.string = T)

bl = read.table('processed/new.born.exons/old/same.exon.out.final.tab',sep='\t') # DB didn't include adjasted exons
colnames(bl) = c('qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','qcovhsp','qcovs','qlen')
bl$map.frac = (bl$qend-bl$qstart+1)/bl$qlen
plot(bl$qcovs,bl$map.frac,pch='.')
bl[bl$qcovs > 80 & bl$map.frac < 0.1,][1:10,] # so, probably qcovs is a total coverage of query by hits 

hist(bl$evalue)
hist(bl$map.frac)
plot(bl$qcovhsp,bl$evalue)
length(unique(bl$qseqid))
e[setdiff(names(e),bl$qseqid)] #all were found

bl = bl[bl$map.frac> 0.5,]
e[setdiff(names(e),bl$qseqid)] #all were found

bl$qid = sapply(strsplit(bl$qseqid,'.',T),'[',3)
bl = split(bl,bl$qid)
length(bl) 
bl[[1]]

blast.sp = sapply(bl,function(x){paste(species$short[rownames(species) %in% sapply(strsplit(x$sseqid,'.',T),'[',1)],collapse = '')})
#blast.sp = sapply(bl,function(x){x=x[substr(x$qseqid,1,3)!=substr(x$sseqid,1,3),];paste(species$short[rownames(species) %in% sapply(strsplit(x$sseqid,'.',T),'[',1)],collapse = '')})

sort(table(blast.sp))

obs.sp = sapply(exon.birth.one,function(x)paste(species$short[rownames(species) %in% rownames(x)[!is.na(x$seg_id)]],collapse = ''))
blast.sp = blast.sp[names(obs.sp)]

sps = c(species$short,'hq','mr','mrb','hqmrb','hqmrbo','hqmrboc')
spsa = unique(c(sps,blast.sp,obs.sp))

# look on Ns
nstat = as.data.frame(do.call(rbind,strsplit(names(iei),'.',T)))[,-2]
colnames(nstat) = c('species','nb.exon')
nstat$ncount = nchar(iei) - nchar(gsub('n','',iei))
nstat = as.matrix(cast(nstat,nb.exon ~ species,value = 'ncount',fun.aggregate = identity))

nstat[1:10,]
table(apply(nstat>0,1,sum))
good.nbe = rownames(nstat)[apply(nstat>0,1,sum)==0]

tn = t(table(factor(blast.sp[good.nbe],levels = spsa),factor(obs.sp[good.nbe]       ,levels = spsa))[sps,sps[-13]])
t  = t(table(factor(blast.sp          ,levels = spsa),factor(obs.sp[names(blast.sp)],levels = spsa))[sps,sps[-13]])

tn = t(table(factor(blast.sp[good.nbe],levels = spsa),factor(obs.sp[good.nbe]       ,levels = spsa))[   ,sps[-13]])
t  = t(table(factor(blast.sp          ,levels = spsa),factor(obs.sp[names(blast.sp)],levels = spsa))[   ,sps[-13]])
t = t[,apply(t,2,sum)>4]
tn = tn[,apply(tn,2,sum)>4]
pdf('figures/newborn.exons.RNAseq-vs-blast.pdf',w=12,h=6)
par(mfcol=c(1,2),tck=-0.01,mgp=c(1.1,0.2,0),mar=c(5,5,1.5,0),oma=c(0,0,2,1))
imageWithText(sweep(t ,1,apply(t ,1,max),'/'),t ,names.as.labs = T,xlab='RNA-Seq',ylab='Blast',main='All newborn exon')
imageWithText(sweep(tn,1,apply(tn,1,max),'/'),tn,names.as.labs = T,xlab='RNA-Seq',ylab='Blast',main='No Ns in all species')
dev.off()

good.nbe[blast.sp[good.nbe] == 'hq' & obs.sp[good.nbe] =='h' ]
e['human.e.t116']
exon.birth.one[['t138']]



# check that orth-newborn are blasted to the same position #####
# _use blast in adj exons - to check exon split
bl = read.table('processed/new.born.exons/same.exon.out.final.tab',sep='\t') # DB didn't include adjasted exons
colnames(bl) = c('qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','qcovhsp','qcovs','qlen')
bl$map.frac = (bl$qend-bl$qstart+1)/bl$qlen
bl$qid = sapply(strsplit(bl$qseqid,'.',T),'[',3)
dim(bl)
bl = bl[bl$map.frac> 0.5,]
e[setdiff(names(e),bl$qseqid)] #all were found
bl = split(bl,bl$qid)
bl = bl[names(exon.birth.one)]
length(bl) 
bl[[1]]

# add relative coordinates
exon.birth.one = lapply(exon.birth.one,function(x){
	x$rel.stop = x$rel.start = NA
	for(i in 1:nrow(x)){
		if(is.na(x$strand[i])) next
		if(x$strand[i] == 1){
			# for adj exons not included 
#			x$rel.start[i] = x$start[i] - all.anns.[[rownames(x)[i]]][x$useg_id[i],'stop'] + 1 
#			x$rel.stop[i]  = x$stop[i]  - all.anns.[[rownames(x)[i]]][x$useg_id[i],'stop'] + 1
			x$rel.start[i] = x$start[i] - all.anns.[[rownames(x)[i]]][x$useg_id[i],'start'] + 1
			x$rel.stop[i]  = x$stop[i]  - all.anns.[[rownames(x)[i]]][x$useg_id[i],'start'] + 1
		}else{
			# for adj exons not included 
#			x$rel.start[i] = all.anns.[[rownames(x)[i]]][x$useg_id[i],'start'] + 1 - x$stop[i]
#			x$rel.stop[i]  = all.anns.[[rownames(x)[i]]][x$useg_id[i],'start'] + 1 - x$start[i]
			x$rel.start[i] = all.anns.[[rownames(x)[i]]][x$useg_id[i],'stop'] + 1 - x$stop[i]
			x$rel.stop[i]  = all.anns.[[rownames(x)[i]]][x$useg_id[i],'stop'] + 1 - x$start[i]
		}
	}
	x
	})

for(r in names(bl)){
	bl[[r]]$orth.pos =  NA
	bl[[r]]$overlap.adj.exon = FALSE
	ss = sapply(strsplit(bl[[r]]$sseqid,'.',T),'[',1)
	nb = exon.birth.one[[r]]
	for(i in 1:nrow(bl[[r]])){
		bl[[r]]$orth.pos[i] = bl[[r]]$send[i] >= nb[ss[i],'rel.start'] & bl[[r]]$sstart[i] <= nb[ss[i],'rel.stop']
		bl[[r]]$overlap.adj.exon[i] = bl[[r]]$sstart[i] <= all.anns.[[ss[i]]][nb[ss[i],'useg_id'],'length'] | bl[[r]]$send[i] >= (nchar(iei[[bl[[r]]$sseqid[i]]])-all.anns.[[ss[i]]][nb[ss[i],'dseg_id'],'length'])
	}
}

nb.stat  = data.frame(nstat[names(exon.birth.one),])
nb.stat$adj.exons = sapply(bl,function(x)mean(x$overlap.adj.exon))

hist(nb.stat$adj.exons)
table(nb.stat$adj.exons>0)
which(nb.stat$adj.exons>0 & obs.sp == 'h')
# clean bl: retain only longest from correct matches

blcl = lapply(bl,function(x){
	x = x[!x$overlap.adj.exon,]
	x$from = species[sapply(strsplit(x$qseqid,'.',T),'[',1),'short']
	x$to = species[sapply(strsplit(x$sseqid,'.',T),'[',1),'short']
	x$species=paste0(x$from,'->',x$to)
	xna = x[is.na(x$orth.pos) ,]
	x = x[!is.na(x$orth.pos) & x$orth.pos,]
	x = do.call(rbind,lapply(split(x,x$species),function(z){
		z[order(z$qend-z$qstart,decreasing = T)[1],]
		}))
	x=rbind(x,xna)
	x = x[order(x$to,x$from),]
	x$to = x$from = NULL
	x
	})

nb.stat$full.obs = sapply(1:length(blcl),function(i){
	if(nchar(obs.sp[i])==1) return(1)
	exp = combn(strsplit(obs.sp[i],'')[[1]],2)
	exp = c(paste0(exp[1,],'->',exp[2,]),paste0(exp[2,],'->',exp[1,]))
	mean(exp %in% blcl[[i]]$species)
	})

nb.stat$full.obs1 = sapply(1:length(blcl),function(i){
	if(nchar(obs.sp[i])==1) return(1)
	exp = strsplit(obs.sp[i],'')[[1]]
	o = do.call(rbind,strsplit(blcl[[i]]$species,'->',T))
	mean(exp %in% o[o[,1] != o[,2],2])
})
hist(nb.stat$full.obs1)




# _check gene orthology #####
orth.ens.genes = readRDS('Rdata/orth.ens.genes.Rdata')
# seg2ens = readRDS('Rdata/seg2ens.Rdata')
# egids = lapply(exon.birth.one,function(x){
# 	r = list()
# 	for(s in 1:nrow(x)){
# 		r[[s]] = unique(unlist(seg2ens[[s]][c(x$useg_id[s],x$dseg_id[s])]))
# 	}
# 	r
# 	})
# saveRDS(egids,'Rdata/newborn.exons.ens.ids.Rdata')
egids = readRDS('Rdata/newborn.exons.ens.ids.Rdata')
oinx = unlist(lapply(1:ncol(orth.ens.genes),function(i)setNames(1:nrow(orth.ens.genes),orth.ens.genes[,i])))
nb.stat$orthn = sapply(egids,function(x){
	x = unlist(x)
	x = x[x %in% names(oinx)]
	max(table(oinx[x]))
	})

# check Ns and missed exons #####
# mark cases when exon is not observed in RNA-Seq in given species and there are Ns
table(names(obs.sp) == names(blcl))
s = species[colnames(nb.stat)[1:7],'short']
s=!t(sapply(strsplit(obs.sp,''),function(x){s %in% x}))
fisher.test(table(s,nb.stat[,1:7]>0))
nb.missed.due.n = s & nb.stat[,1:7]>0
table(apply(nb.missed.due.n,1,sum))
table(apply(nb.stat[,1:7]>0,1,sum))



# check that if exon is present in s1 but not in s2 could it be explained by Ns in s2 #####

pv = or = matrix(NA,ncol=7,nrow=7,dimnames = list(rownames(species),rownames(species)))
nthr = 00
for(s1 in rownames(species)){
	for(s2 in rownames(species)){
		if(s1 == s2) next
		f1 = grepl(species[s1,'short'],obs.sp)
		f2 = grepl(species[s2,'short'],obs.sp)
		ft = fisher.test(table(not.obs=factor(!f2[f1],levels=c(F,T)),has.n=factor(nb.stat[[s2]][f1] > nthr,levels=c(F,T))))
		pv[s1,s2] = ft$p.value
		or[s1,s2] = ft$estimate
	}
}

pv1 = or1 = matrix(NA,ncol=7,nrow=6,dimnames = list(1:6,rownames(species)))
for(i in 1:6)
	for(s1 in rownames(species)){
		f1 = nchar(obs.sp) == i
		f2 = grepl(species[s1,'short'],obs.sp)
		ft = fisher.test(table(not.obs=factor(!f2[f1],levels=c(F,T)),has.n=factor(nb.stat[[s1]][f1] > nthr,levels=c(F,T))))
		pv1[i,s1] = ft$p.value
		or1[i,s1] = ft$estimate
	}


pv2 = or2 = matrix(NA,ncol=7,nrow=7,dimnames = list(rownames(species),rownames(species)))
nthr = 00
for(s1 in rownames(species)){
	for(s2 in rownames(species)){
		if(s1 == s2) next
		f1 = grepl(species[s1,'short'],obs.sp) & nchar(obs.sp) < 6
		f2 = grepl(species[s2,'short'],obs.sp)
		ft = fisher.test(table(not.obs=factor(!f2[f1],levels=c(F,T)),has.n=factor(nb.stat[[s2]][f1] > nthr,levels=c(F,T))))
		pv2[s1,s2] = ft$p.value
		or2[s1,s2] = ft$estimate
	}
}

pdf('figures/newborn.exon.loss-vs-Ns.pdf',w=21,h=7)
par(mfrow=c(1,3),tck=-0.01,mgp=c(3,0.5,0),mar=c(4,4,3,6))
or. = round(or,2)
diag(or.) = paste0('n(N)=\n',apply(nb.stat[,rownames(species)]>0,2,sum))

or2. = round(or2,2)
diag(or2.) = paste0('n(N)=\n',apply(nb.stat[,rownames(species)]>0,2,sum))
cols=c('red','yellow','gray')
breaks=c(0,0.001,0.05,1)
imageWithText(pv,or.,names.as.labs = T,xlab='Observed',ylab='Missed',breaks=breaks,col=cols,main='At-least-one-N vs exon loss (Odds ratio)')
legend(7.5,7.5,xpd=T,fill=rev(cols),legend = paste('<',rev(breaks[-1])),title='p-value')
imageWithText(pv1,round(or1,2),names.as.labs = T,xlab='# of obs species',ylab='Missed in species',breaks=breaks,col=cols,main='At-least-one-N vs exon loss (Odds ratio)')
legend(6.5,7.5,xpd=T,fill=rev(cols),legend = paste('<',rev(breaks[-1])),title='p-value')

imageWithText(pv2,or2.,names.as.labs = T,xlab='Observed',ylab='Missed',breaks=breaks,col=cols,main='At-least-one-N vs exon loss (Odds ratio)\nonly observed in five and less species')
legend(7.5,7.5,xpd=T,fill=rev(cols),legend = paste('<',rev(breaks[-1])),title='p-value')
dev.off()

f1 = nchar(obs.sp) > 1
f2 = grepl('r',obs.sp)
fisher.test(table(not.obs=factor(!f2[f1],levels=c(F,T)),has.n=factor(nb.stat[['rat']][f1] > nthr,levels=c(F,T))))

t = data.frame(species=rep(rownames(species),each=nrow(nb.stat)),
					 obs = rep(obs.sp,times=nrow(species)),
					 n=as.numeric(as.matrix(nb.stat[,rownames(species)])))
t$short = species[t$species,'short']
t$here = sapply(1:nrow(t),function(i)grepl(t$short[i],t$obs[i]))
t$obs.minus.short = sapply(1:nrow(t),function(i)gsub(t$short[i],'',t$obs[i]))
table(t$here,nchar(t$obs.minus.short))
nthr = 00
cms = lapply(1:6,function(i){z=t[nchar(t$obs)==i & t$species %in% c('rat','macaque','rabbit','opossum'),];table(obs=!z$here,n=z$n>nthr)})
sapply(cms,function(x)c(fisher.test(x)$p.value,fisher.test(x)$estimate))

t = table(obs.sp,nb.stat[,'rat']>0)
t = t[order(t[,2]/t[,1]),]
t[nchar(rownames(t))==5 & apply(t,1,sum)>0,]

# t = table(obs.sp)
# t[order(nchar(names(t)),t)]
nrow(nb.stat) # 6536
table(nb.stat$orthn)

table(nb.stat$orthn == 7) # 4903
table(nb.stat$orthn == 7 & nb.stat$full.obs==1) # 4215
table(nb.stat$orthn == 7 & nb.stat$full.obs==1 & nb.stat$adj.exons==0) # 4073
table(nb.stat$orthn == 7 & nb.stat$full.obs==1 & nb.stat$adj.exons==0 & apply(nb.missed.due.n,1,sum)==0) # 2406
f = nb.stat$orthn == 7 & nb.stat$full.obs==1 & nb.stat$adj.exons==0 & apply(nb.missed.due.n,1,sum)==0
table(blast.sp,f)[sps,]
table(obs.sp,f)[sps[-13],]

which(nb.stat$adj.exons > 0)
exon.birth.one[116]
table(nb.stat$adj.exons > 0,nchar(obs.sp))
# compare coverage in newborn exon and in adjacent introns
for(t in names(exon.birth.one)){
	c = readRDS(paste0('Rdata/newborn.exon.cov/',obs.sp[t],'-',blast.sp[t],'.',t,'.Rdata'))
	nb = exon.birth.one[[t]]
	exon.birth.one[[t]]$exon.cov = exon.birth.one[[t]]$int.cov = NA
	for(s in rownames(nb)[!is.na(nb$seg_id)]){
		x = c[[s]]$x >= nb[s,'start'] & c[[s]]$x <= nb[s,'stop']
		
		if(nb[s,'strand']==1){
			u = all.anns.[[s]][nb[s,'useg_id'],'stop']+1
			d = all.anns.[[s]][nb[s,'dseg_id'],'start']-1
		}else{
			u = all.anns.[[s]][nb[s,'dseg_id'],'stop']+1
			d = all.anns.[[s]][nb[s,'useg_id'],'start']-1
		}
		xx = !x & c[[s]]$x >= u &  c[[s]]$x <= d
		exon.birth.one[[t]][s,'exon.cov'] = mean(c[[s]]$cov[x])
		exon.birth.one[[t]][s,'int.cov'] = mean(c[[s]]$cov[xx])
		# j = c[[s]]$juncs
		# j[j$start==u & j$end == d,]
		# j[j$start==u & j$end == nb[s,'start']-1,]
		# j[j$start == nb[s,'stop']+1,]# & j$end==d,]
	}
}
# save(exon.birth.one,nb.stat,bl,blcl,file = 'Rdata/tmp.exon.birth.Rdata')
# 

t = do.call(rbind,exon.birth.one)
t = t[!is.na(t$int.cov),]
plotLine(t$length,t$int.cov/t$exon.cov,ylim=c(0,2),xlim=c(0,1000))
fisher.test(table(len=t$length>100,t$int.cov/t$exon.cov>0.5))
boxplot(split(t$int.cov/t$exon.cov,number2bin(t$length,20)),outline=F)
hist(t$length[t$int.cov/t$exon.cov< 0.5],0:1000*10,xlim=c(0,1000),freq = F,border = NA,col='#00FF0088')
hist(t$length[t$int.cov/t$exon.cov>=0.5],0:1000*10,xlim=c(0,1000),add=T,freq = F,border = NA,col='#FF000088')
hist(t$int.cov/t$exon.cov,5000,xlim=c(0,2),xlab='mean(intron cov)/mean(newborn exon cov)',main='')

ier.max = unlist(lapply(exon.birth.one,function(x)(max(x$int.cov/x$exon.cov,na.rm=T))))
ier.min = unlist(lapply(exon.birth.one,function(x)(min(x$int.cov/x$exon.cov,na.rm=T))))
len =     unlist(lapply(exon.birth.one,function(x)(min(x$length,na.rm=T))))
plot(ier.max,ier.min,xlim=c(0,1),ylim=c(0,1))
table(ier.max<0.5,f)
names(bl)[ier.max>0.5 & len < 100 &  f]# & blast.sp == 'h' & obs.sp == 'h']
names(bl)[ier.max>0.5 & len < 100 &  f]# & blast.sp == 'h' & obs.sp == 'h']

# take only these that have at least 4 samples (or at least 2 testis samples) in at least one species with PSI>0.2 
p = unlist(setNames(lapply(born.exn.sajr,function(x)apply(x$ir[,colnames(x$ir) %in% rownames(meta)],1,function(z)sum(z>0.2,na.rm=T))),NULL))
high.psi = sapply(exon.birth.one,function(x)sum(p[x$seg_id[!is.na(x$seg_id)]]>3))

pt = unlist(setNames(lapply(born.exn.sajr,function(x)apply(x$ir[,colnames(x$ir) %in% rownames(meta)[meta$tissue=='testis']],1,function(z)sum(z>0.2,na.rm=T))),NULL))
high.psi.t = sapply(exon.birth.one,function(x)sum(pt[x$seg_id[!is.na(x$seg_id)]]>1))
table(high.psi>0,testis=high.psi.t>0)

sum(ier.max<0.5  &  f & len < 200)
sum(ier.max<0.5  &  f & len < 200 & high.psi > 0)

names(bl)[ier.max<0.5 & len < 100 &  f & high.psi > 0 & obs.sp == 'hq']

hist(ier.min,15000,xlim=c(0,2))
abline()
exon.birth.one[which(is.infinite(ier))[1]]

maxcov = max(c$human$cov,c$human$juncs$score)
plotReadCov(c$human,reverse = (all.anns.[[s]][nb[s,'useg_id'],'strand'] == -1),min.junc.cov = 100,bty='n',xlab=paste0('Chr ',all.anns.[[s]][nb[s,'useg_id'],'chr_id']),ylab='Reads',ylim=c(-0.2*maxcov,maxcov),plot.junc.only.within = T)


t  = t(table(factor(blast.sp[f]          ,levels = spsa),factor(obs.sp[names(blast.sp)[f]],levels = spsa))[sps,sps[-13]])
t = t[,apply(t,2,sum)>4]
imageWithText(sweep(t ,1,apply(t ,1,max),'/'),t ,names.as.labs = T,xlab='RNA-Seq',ylab='Blast',main='All newborn exon')

names(bl)[f & blast.sp == 'h' & obs.sp == 'h']

s='human'
s='macaque'
bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue=='brain' & meta$days<300],'.bam')
nb = exon.birth.one[['t887']]

bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue=='cerebellum' & meta$days>10000],'.bam')
nb = exon.birth.one[['t2397']]

nb = exon.birth.one[['t47']]
bams = paste0('processed/mapping/hisat2.s/',s,'/',meta[names(sort(born.exn.sajr$human$ir['hum.16953.s116',colnames(born.exn.sajr$human$ir) %in% rownames(meta)],decreasing = T)[1:30]),'fname'],'.bam')

if(all.anns.[[s]][nb[s,'useg_id'],'strand'] == 1){
	start = all.anns.[[s]][nb[s,'useg_id'],'start']
	stop = all.anns.[[s]][nb[s,'dseg_id'],'stop']
}else{
	start = all.anns.[[s]][nb[s,'dseg_id'],'start']
	stop = all.anns.[[s]][nb[s,'useg_id'],'stop']
}

cov=getReadCoverage(bams,all.anns.[[s]][nb[s,'useg_id'],'chr_id'],start - 100,stop+100,-all.anns.[[s]][nb[s,'useg_id'],'strand'])
maxcov = max(cov$cov,cov$juncs$score)
plotReadCov(cov,reverse = (all.anns.[[s]][nb[s,'useg_id'],'strand'] == -1),min.junc.cov = 5,bty='n',xlab=paste0('Chr ',all.anns.[[s]][nb[s,'useg_id'],'chr_id']),ylab='Reads',ylim=c(-0.2*maxcov,maxcov),plot.junc.only.within = T)


# look on t887 #####
library(msa)
blcl['t887']
exon.birth.one[['t887']]
h2te[h2te$seg_id=='hum.21203.s30',]
rm[4276702,]
paste(rev(c(A='T','T'='A',C='G',G='C')[strsplit('ACACAGTGGCTCATCCTTGTAATCCGAGCACTTTGGGAGGTGGAGGTGGGCAGATCACTTGAGTCCA','')[[1]]]),collapse='')

t887 = as.character(iei[paste0(rownames(species),'.iei.t887')])
names(t887) = rownames(species)
t887 = c(t887,chimp=tolower('CAGTGTACAAGAAGAAGACCCAAATGAAGAGCTTTCAAAAGATGAGTTCATTCTGAAGTTAAAGGCAGAAGTACAGCGTTTGCTGGGTAGCAACTCAATGAAGCGTCATCTGGTGTCTCAGTTACAAAATGACCTCAAAGACTGTCATAAGAAAATTGAAGATCTCCACCAAGTGAAGAAGGATGAAAAAAGCATTGAGGTTGAGGTCTGTAAAGCCTATAGAGTTAAAAAACCATTTTTAAAAGGCTCAATAAAGAAAACTTTTACACTTTTTTTTTTTTGAGACAGGGCATCACTCTGTCTTCCAGGCTGGAGGGCAGTGGTGGAAACATGGCTCCCTGTAGCCTCAACTCCCAGGCTAAAGCGATGCTCCCACCTCAGCCTCCTGAGTAGCTGGGACCACAGGCGCACGCCACTATGCCTGGCTAATTTTTAAATTTTTTGTAGAGACAGGTTCCCACAATGTTGCCCAGGCTGGTCTCAAACTCCTGGGCTCAAGCAATCCTCCCACCTCAATCTCTCAAGGTGTTGGATTACAGGTGTGAATCACCACATCCAGCCCTTTTACATGTTTTTAAATTAATAAACCTTTGGCTTTCAGGAGTTTGTTTGACATTTTTAGAAATAAAGTTTAAAGAAGGAAATACCATTTTTAAAGGGATACTGACAGTCCCTGACTTAGGATGGTTTGACTTACAATTTTTAGACTTTACAGTGGTGGAAAAGTTATGTGCATTCAGTAGAAACCATATTTCAAGTACCTGTAAAACAGTTCTGTTTTTTACTTTCAGTACAGTGTTCAATAAATTTCATGAGATATTCAACAGTTTATTATAAATGAAGCTTTGCATTAGATGACTTTGCATAACTGTTGGCTATTATAAGTGTTCTGAGCACATTTAAGGTAGGCCAGGCTAAGCTATGGTTTTAATAGGTTAGATGCATTAAATGCATTTTTGACTTAATGATATTTTCAACTTATAATGATTTCATTAGGTGATAAACCCATTGTAAGTTGAAAAGCATCTGTATGTTTTTTCTGGTTACTCTTGGCACCTAGTAAATGAGATACCAATTAGGTAAAATTCCTTTACCATAAATACAAGGAAGTATTCTCTGTAACAAGGAATCTGCATTTCTTTAAGATTGGCTTTTACCAATAGTGCTCTAAAATCTAGTTAAAATAAATCCTTCATAGAAATCCTGTAGTATTTTCTATGAGATTTTCCTGGCCCATGTTCTTAAAACAAAAGTCATTCTATATGAGGACATTAATTATATAGTGATGATGGTACAGATATTATATAGTGGTCGTGGTACTATGGGAAAGTACACAAAATCTGGTTATGTCAGTCATTGAAATCCACCTCTGTATTCTTCAGTTCGGAAAATGAAGCTGGAGCACTTTGTGACACTTTTGAATTTACAGAGCCAATTTCAGACAAATTGAAACTAAACACATGACTAAGTAAGTCAAAGATATTTCTGTAAGATCATGCTATCACAGATATCTTTAGGGCTGCTCTCCTTTATAGATCCTTTTTTTTTACCCTTCCATATTTCACAAATCATCAGCCATCTATACTCCCAATAGAAAATAATATATAAAACACCTAATAACTAAGAAAAGCAAAAGAAGAGGGAAAAGGAATGTGGGTAGTGTACTCCTTTTCCTCCTGGGTGGAGATAATTGTTTACTTCCTTCTTTAGTACTGTTTTATCTCCTATTCTTTTAAACTAGGAGTTGTGTTCATTAGAACAGTTCTGCTAAAGTATGTATAATATTCTGTGCTAATTATTGTGGAATTTTAATTGTGGTTAATTTTTATATTTCGGAAAAACAGTGTTTTTGCTATCAAAATCTTTTCCAGCCTGATTCCCTCTTCTGAGTAACTACATTTGATGACTAAAACATAGATTTTAATATTAACAACAAATGAGAGAGACTCCTCCCCGCCCTGCCCCCTGCACATTTAACATCTTTAAATCAAGCTGTTTTCATTGTTCCCCTACATTTCTTTAAAAGGATGGGAATGAAAGCTTGTTTCTCAATTTGTATTAAGACATTAGAGTTTATTAGTTTAAGCATTTTGGTGGGAGTGCTCTAGGATTTAGGATACTAGAAACCTGATTTTGAGTCTGAACCTTAGCTCAGATATCTTGGTGCTTTCTTTGTTCTGTATTAAACTTGCTAATGGCATTCGCTATACAGAAACACTTTATCATGTGTAGCAGTAGATAGTGGTCATGCATTGTGCTGTTTTTTAGCCTTCCCATTTTCCTCTGTAGTCTGGTTGATCAATAGATCAAAGGAATAACTATCTCATCTTTGAGTATTCTTTTTGTCCAGCTACTCTATTGCCATAGTGAATTCTCGGTAGTGTACATCTAAACTCCAATTTTTCTCAAAGTTACTTCAGCTGTGGACTCAAGTGATCTGCCCACCTCCACCTCCCAAAGTGCTCGGATTACAAGGATGAGCCACTGTGTTGTGGCAGAGCTCAGAGTCCAGGTACACAACAAGTGAAACGGAGTAGATATATTTTTAAAATCTTTGTAACCCCCATCAATGAATAAAAAGAAAAGTAGATAAGGAAAACTCTCAGAGGAAAATATGCAAAGTTAAATGAACAGATAAATCTGGGAAGGTGGAAAAACAAATGAACAATGAAGGTATATAAGAGGTTCCACTTTGGTGGTGACACAAGAAATACAAATGAAAGCAATTGTTAGATATTCTTTTCCCATCAAAGATAATTATTTAAAGATAATCACTGCCATTGGAGGCAAAGAGAGATGGGCAGACCCATGTACTACCTGACAGTATACATTTGGATGACATTTCTGGAAGGCAATTCGCCAAAGTCTATCATGAGTCTTAAAAATCGTATGTATTGAATAAAATATATGTCCATGAGTGAATGCTGACATAAATGAATACATACAGAAATAGAATAAATGCTTGATCTTAGTAGGATAACATATGATGAGAAATAGGATATTTGTATTGTCTCAAAGTATCCTATAAGATACTTATCATTTATTATAAAACAAAAAATAGTAAGTTTACAATGGAGATGACAGATCCCACCTTCAATTGATTATTACCAGTGGTGAGACATAGTGATATCATGTGTCTCCTGCAAAGAGCACATATCAGTTCTTGACATGCATGCCAAAGATGCATAACCTGATTTTAATAAACCTAAAAACTCCAAATTGAGGAACATCCTACAAAATAACGGGTTAGTATTCTTAGAAAACATCATGGGTTTGAAAGACTGAAGAGCTGTCCTAGATTGGAGGAGACTAATGAAACATGACAACTAAATGCCACTCAGAATCATGGATTTGATCTTGGACCAGAAAACAGACATTTGGGGCAAATTGGAGACATTTGTATAAGGTCTGTAGATTACTTGATAGTATTGCATCTGTGTTAATTTACTAGTTTTGATGATTGTGCTGTAGTTATATTGGGGGAAGATGAGTAAAAGCTATATGGAGATTATTTGTACTGTTTGCAATTTTTTTAAAGTATGGAATTATTTCAAAGTGAAACTTAGAATTTTTTAAAACTTGTAAAACTTTTGACTTTTTAATATTCAGACTAAAACAGATACCTCAGAAAAACCAAAGAATCAATTATGGCCTGAGTCTTCTACTTCTGATGTTGTCAGAGATGATATTCTGCTGCTTAAAAATGAAATTCAAGTTTTACAACAACAAAATCAG'))
al=msa(t887,method='Muscle',type='dna')
al.hq=msa(t887[1:2],method='Muscle',type='dna')
as.character(al.hq)
al.hqm=msa(t887[1:3],method='Muscle',type='dna')
as.character(al.hqm)
al.hpqm=msa(t887[c(1:3,8)],method='ClustalOmega',type='dna')
as.character(al.hpqm)
exon.birth.one[['t887']]
#in testis there are one more exon (missed in annotation): 49063417-49063504
# coverage: {humna*macaque*mouse}*{brain*testis}*{earliest*adult}
par(mfrow=c(4,3))
nb = exon.birth.one[['t887']]
covs = list()

s='human'
if(all.anns.[[s]][nb[s,'useg_id'],'strand'] == 1){
	start = all.anns.[[s]][nb[s,'useg_id'],'start']
	stop = all.anns.[[s]][nb[s,'dseg_id'],'stop']
}else{
	start = all.anns.[[s]][nb[s,'dseg_id'],'start']
	stop = all.anns.[[s]][nb[s,'useg_id'],'stop']
}
chr = all.anns.[[s]][nb[s,'useg_id'],'chr_id']
strand = all.anns.[[s]][nb[s,'useg_id'],'strand']
hcoor  = list(start=start,stop=stop,chr=chr,strand=strand,s=s)

hb1 = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue %in% c('brain','cerebellum') & meta$days<70],'.bam')
hb2 = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue %in% c('brain','cerebellum') & meta$days>15*365],'.bam')
ht1 = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue %in% c('testis') & meta$days<70],'.bam')
ht2 = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue %in% c('testis') & meta$days>15*365],'.bam')
covs$hb1 = getReadCoverage(hb1,chr,start - 100,stop+100,-strand)
covs$hb2 = getReadCoverage(hb2,chr,start - 100,stop+100,-strand)
covs$ht1 = getReadCoverage(ht1,chr,start - 100,stop+100,-strand)
covs$ht2 = getReadCoverage(ht2,chr,start - 100,stop+100,-strand)



s = 'macaque'
if(all.anns.[[s]][nb[s,'useg_id'],'strand'] == 1){
	start = all.anns.[[s]][nb[s,'useg_id'],'start']
	stop = all.anns.[[s]][nb[s,'dseg_id'],'stop']
}else{
	start = all.anns.[[s]][nb[s,'dseg_id'],'start']
	stop = all.anns.[[s]][nb[s,'useg_id'],'stop']
}
chr = all.anns.[[s]][nb[s,'useg_id'],'chr_id']
strand = all.anns.[[s]][nb[s,'useg_id'],'strand']
qcoor  = list(start=start,stop=stop,chr=chr,strand=strand,s=s)

qb1 = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue %in% c('brain','cerebellum') & meta$days<130],'.bam')
qb2 = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue %in% c('brain','cerebellum') & meta$days>10*365],'.bam')
qt1 = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue %in% c('testis') & meta$days<130],'.bam')
qt2 = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue %in% c('testis') & meta$days>10*365],'.bam')

covs$qb1 = getReadCoverage(qb1,chr,start - 100,stop+100,-strand)
covs$qb2 = getReadCoverage(qb2,chr,start - 100,stop+100,-strand)
covs$qt1 = getReadCoverage(qt1,chr,start - 100,stop+100,-strand)
covs$qt2 = getReadCoverage(qt2,chr,start - 100,stop+100,-strand)

s = 'mouse'
if(all.anns.[[s]][nb[s,'useg_id'],'strand'] == 1){
	start = all.anns.[[s]][nb[s,'useg_id'],'start']
	stop = all.anns.[[s]][nb[s,'dseg_id'],'stop']
}else{
	start = all.anns.[[s]][nb[s,'dseg_id'],'start']
	stop = all.anns.[[s]][nb[s,'useg_id'],'stop']
}
chr = all.anns.[[s]][nb[s,'useg_id'],'chr_id']
strand = all.anns.[[s]][nb[s,'useg_id'],'strand']
mcoor  = list(start=start,stop=stop,chr=chr,strand=strand,s=s)

mb1 = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue %in% c('brain','cerebellum') & meta$days<15],'.bam')
mb2 = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue %in% c('brain','cerebellum') & meta$days>40],'.bam')
mt1 = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue %in% c('testis') & meta$days<15],'.bam')
mt2 = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue %in% c('testis') & meta$days>40],'.bam')

covs$mb1 = getReadCoverage(mb1,chr,start - 100,stop+100,-strand)
covs$mb2 = getReadCoverage(mb2,chr,start - 100,stop+100,-strand)
covs$mt1 = getReadCoverage(mt1,chr,start - 100,stop+100,-strand)
covs$mt2 = getReadCoverage(mt2,chr,start - 100,stop+100,-strand)

rm(chr,s,start,stop,strand)

pdf('figures/CEP152.t887.human-specific.newborn.exon.pdf',w=12,h=12)
layout(matrix(1:15,ncol=3),heights = c(1,1,1,1,0.3))
par(tck=-0.02,mgp=c(1.3,0.2,0),mar=c(0.3,2.5,1,0),oma=c(0,0,0,1))
attach(hcoor)
plotReadCov(covs$hb1,reverse = (strand == -1),min.junc.cov = 2,bty='n',ylab='Reads',xaxt='n',plot.junc.only.within = T,main='Human brain < 10wpc',xlab='')
plotReadCov(covs$hb2,reverse = (strand == -1),min.junc.cov = 2,bty='n',ylab='Reads',xaxt='n',plot.junc.only.within = T,main='Human brain > 15wpc',xlab='')
plotReadCov(covs$ht1,reverse = (strand == -1),min.junc.cov = 2,bty='n',ylab='Reads',xaxt='n',plot.junc.only.within = T,main='Human testis < 10wpc',xlab='')
plotReadCov(covs$ht2,reverse = (strand == -1),min.junc.cov = 2,bty='n',ylab='Reads',xaxt='n',plot.junc.only.within = T,main='Human testis > 15ypc',xlab='')
xlim=range(covs$hb1$x)
if(strand==-1) xlim=rev(xlim)
par(mar=c(2.3,2.5,1.5,0))
plot(1,t='n',yaxt='n',ylab='',xlab=paste('Chr',chr),bty='n',ylim=c(-0.2,1.2),xlim=xlim)
abline(h=0.5)
seg=all.anns.[[s]][na.omit(unlist(nb[s,c('useg_id','seg_id','dseg_id')])),]
rect(seg$start,0,seg$stop,1,col=c('black','red','black'),border=NA)
detach(hcoor)

# macaque
attach(qcoor)
par(mar=c(0.3,2.5,1,0))
plotReadCov(covs$qb1,reverse = (strand == -1),min.junc.cov = 2,bty='n',ylab='Reads',xaxt='n',plot.junc.only.within = T,main='Macaque brain < 130dpc',xlab='')
plotReadCov(covs$qb2,reverse = (strand == -1),min.junc.cov = 2,bty='n',ylab='Reads',xaxt='n',plot.junc.only.within = T,main='Macaque brain > 10ypc',xlab='')
plotReadCov(covs$qt1,reverse = (strand == -1),min.junc.cov = 2,bty='n',ylab='Reads',xaxt='n',plot.junc.only.within = T,main='Macaque testis < 130dpc',xlab='')
plotReadCov(covs$qt2,reverse = (strand == -1),min.junc.cov = 2,bty='n',ylab='Reads',xaxt='n',plot.junc.only.within = T,main='Macaque testis > 10ypc ',xlab='')
xlim=range(covs$qb1$x)
if(strand==-1) xlim=rev(xlim)
par(mar=c(2.3,2.5,1.5,0))
plot(1,t='n',yaxt='n',ylab='',xlab=paste('Chr',chr),bty='n',ylim=c(-0.2,1.2),xlim=xlim)
abline(h=0.5)
seg=all.anns.[[s]][na.omit(unlist(nb[s,c('useg_id','seg_id','dseg_id')])),]
rect(seg$start,0,seg$stop,1,col='black')
detach(qcoor)


# mouse
attach(mcoor)
par(mar=c(0.3,2.5,1,0))
plotReadCov(covs$mb1,reverse = (strand == -1),min.junc.cov = 2,bty='n',ylab='Reads',xaxt='n',plot.junc.only.within = T,main='Mouse brain < 15dpc',xlab='')
plotReadCov(covs$mb2,reverse = (strand == -1),min.junc.cov = 2,bty='n',ylab='Reads',xaxt='n',plot.junc.only.within = T,main='Mouse brain > 40dpc',xlab='')
plotReadCov(covs$mt1,reverse = (strand == -1),min.junc.cov = 2,bty='n',ylab='Reads',xaxt='n',plot.junc.only.within = T,main='Mouse testis < 15dpc',xlab='')
plotReadCov(covs$mt2,reverse = (strand == -1),min.junc.cov = 2,bty='n',ylab='Reads',xaxt='n',plot.junc.only.within = T,main='Mouse testis > 40dpc',xlab='')
xlim=range(covs$mb1$x)
if(strand==-1) xlim=rev(xlim)
par(mar=c(2.3,2.5,1.5,0))
plot(1,t='n',yaxt='n',ylab='',xlab=paste('Chr',chr),bty='n',ylim=c(-0.2,1.2),xlim=xlim)
abline(h=0.5)
seg=all.anns.[[s]][na.omit(unlist(nb[s,c('useg_id','seg_id','dseg_id')])),]
rect(seg$start,0,seg$stop,1,col='black')
detach(mcoor)
dev.off()
# check coverage #####
id=9
id='t1008'
id ='t1020'
which(nchar(obs.sp) == 5)
blcl[id]
bl[[id]][order(sapply(strsplit(bl[[id]]$sseqid,'.',T),'[',1)),]
exon.birth.one[[id]]

s=sort(table(blast.sp))
o=sort(table(obs.sp))
s[nchar(names(s))==5]
o[nchar(names(o))==5]

# TE ####
#  repeats were d from ucsc
params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)

source('code/r.functions/paper.figures.4.F.R')
#h = readRDS('Rdata/human.as.u.all.Rdata')
anns = readRDS('Rdata/anns.Rdata')
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
border.stages = readRDS('Rdata/border.stages.Rdata')
alt.sp = readRDS('Rdata/paper.figures/alt.sp.Rdata')
age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)
age.dpsi$macaque = cbind(age.dpsi$macaque[,1:5],ovary=NaN,age.dpsi$macaque[,6,drop=FALSE])

ier.max = unlist(lapply(exon.birth.one,function(x)(max(x$int.cov/x$exon.cov,na.rm=T))))
len =     unlist(lapply(exon.birth.one,function(x)(min(x$length,na.rm=T))))
p = unlist(setNames(lapply(born.exn.sajr,function(x)apply(x$ir[,colnames(x$ir) %in% rownames(meta)],1,function(z)sum(z>0.2,na.rm=T))),NULL))
high.psi = sapply(exon.birth.one,function(x)sum(p[x$seg_id[!is.na(x$seg_id)]]>3))

obs.sp = sapply(exon.birth.one,function(x)paste(species$short[rownames(species) %in% rownames(x)[!is.na(x$seg_id)]],collapse = ''))
table(names(obs.sp) == names(blcl))
s = species[colnames(nb.stat)[1:7],'short']
s=!t(sapply(strsplit(obs.sp,''),function(x){s %in% x}))
fisher.test(table(s,nb.stat[,1:7]>0))
nb.missed.due.n = s & nb.stat[,1:7]>0

p = unlist(setNames(lapply(born.exn.sajr,function(x)apply(x$ir[,colnames(x$ir) %in% rownames(meta)],1,function(z)sum(z>0.2,na.rm=T))),NULL))
high.psi = sapply(exon.birth.one,function(x)sum(p[x$seg_id[!is.na(x$seg_id)]]>3))

pt = unlist(setNames(lapply(born.exn.sajr,function(x)apply(x$ir[,colnames(x$ir) %in% rownames(meta)[meta$tissue=='testis']],1,function(z)sum(z>0.2,na.rm=T))),NULL))
high.psi.t = sapply(exon.birth.one,function(x)sum(pt[x$seg_id[!is.na(x$seg_id)]]>1))
table(high.psi>0,testis=high.psi.t>0)

filters = list(all = rep(TRUE,length(exon.birth.one)),
							 orth                      = nb.stat$orthn == 7,
							 orth.adj                  = nb.stat$orthn == 7 & nb.stat$adj.exons==0,
							 orth.adj.bl               = nb.stat$orthn == 7 & nb.stat$adj.exons==0 & nb.stat$full.obs==1,
							 orth.adj.bl.int           = nb.stat$orthn == 7 & nb.stat$adj.exons==0 & nb.stat$full.obs==1 & ier.max<0.5,
							 orth.adj.bl.int.len       = nb.stat$orthn == 7 & nb.stat$adj.exons==0 & nb.stat$full.obs==1 & ier.max<0.5 & len < 500,
							 orth.adj.bl.int.len.psi   = nb.stat$orthn == 7 & nb.stat$adj.exons==0 & nb.stat$full.obs==1 & ier.max<0.5 & len < 500 & (high.psi > 0 | high.psi.t > 0),
							 orth.adj.bl.int.len.psi.n = nb.stat$orthn == 7 & nb.stat$adj.exons==0 & nb.stat$full.obs==1 & ier.max<0.5 & len < 500 & (high.psi > 0 | high.psi.t > 0) & apply(nb.missed.due.n,1,sum)==0)
sapply(filters,sum)


rm = read.table('raw/repeats.ucsc/hg19.ucsc.repeat.masker.gz',header = T,comment.char = '')
u = unique(rm$genoName)
names(u) = u

u = gsub('chr..?_','',u)
u = gsub('chr','',u)
u = gsub('_random','',u)
u = toupper(u)
u[u=='M'] = 'MT'
u[grep('GL',u)] = paste0(u[grep('GL',u)],'.1')
length(unique(u))

setdiff(u,h$seg$chr_id)
setdiff(h$seg$chr_id,u)

rm$genoName = u[rm$genoName]
table(rm$genoName)

h2rm = getAnnOverlap(h$seg,data.frame(chr_id=rm$genoName,start=rm$genoStart,stop=rm$genoEnd))
names(h2rm) = rownames(h$seg)
f = h$seg$sites=='ad' & (h$seg$type=='EXN' | rownames(h$seg) %in% rownames(psi.tsm$human))
tall = table(te=h2rm[f]!='-',h$seg$type[f])
fisher.test(tall)

age.segs = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,psi.thr = 0.2,border.stages,s))
names(age.segs) = rownames(species)
for(s in names(age.segs)) age.segs[[s]][is.na(per.tissue.age.qv[[s]][rownames(age.segs[[s]]),colnames(age.segs[[s]])])] = '-'
for(s in names(age.segs)) age.segs[[s]][age.segs[[s]] != '-' & per.tissue.age.qv[[s]][rownames(age.segs[[s]]),colnames(age.segs[[s]])]>0.05] = 'n'

z=lapply(unique(meta$tissue),function(t){
	f = anns$human$sites=='ad' & age.segs$human[,t] != '-'
	z=table(age.segs$human[f,t],h2rm[rownames(age.segs$human)[f]] != '-')[,c('TRUE','FALSE')]
	apply(z,1,function(x)c(x,my.binom.test(x)))
})

pdf('figures/newborn/TE.proportion.in.newborn.devAS.alternif.pdf',w=12,h=6)
par(mfrow=c(2,4),tck=-0.02,mgp=c(1.3,0.2,0),mar=c(3.5,2.5,1.5,0),oma=c(0,0,0,1))
for(i in names(filters)){
	f = filters[[i]]
	h2rm. = h2rm[born.seg.ids[f,1]]
	t = cbind(table(h2rm.!='-',sp.birth[f])[,c('h','hq','hqmrb','hqmrbo')],tall)[c('TRUE','FALSE'),]
	prob=apply(t,2,my.binom.test)
	colnames(prob)[5:6] = c('alt in h','cnst in h')
	plotArea(1:4,t(prob),new=T,col='red',xaxt='n',ylab='proportion of TE',xlim=c(1,6),main=paste0(i,' (',sum(f),')'))
	points(5:6,prob[1,5:6],col=c('orange','black'),pch=19,cex=2)
	segments(5:6,prob[2,5:6],5:6,prob[3,5:6],col=c('orange','black'),lwd=2)
	axis(1,1:6,colnames(prob))
}
par(mfrow=c(1,2),tck=-0.02,mgp=c(1.3,0.2,0),mar=c(3.5,2.5,1.5,0),oma=c(0,0,0,1))
names(z) = unique(meta$tissue)
cols=rep(params$tissue.col[names(z)],each=3)
b=barplot(sapply(z,function(x)x[3,]),col=cols,den=c(40,-1,40),angle=c(-45,0,45),beside = T,border = NA,ylim=range(0,sapply(z,function(x)x[4:5,])),ylab='proportion of TE',main='TE in devAS (qv < 0.05 & dPSI > 0.2)')
segments(b,sapply(z,function(x)x[4,]),b,sapply(z,function(x)x[5,]),col=cols)
legend('topleft',fill='gray',den=c(40,-1,40),angle=c(-45,0,45),legend=c('down','no-change','up'))

gain = c(species$short,'hq','mr','mrb','hqmrb','hqmrbo')
te.alt = apply(table(alt.sp,h2rm[names(alt.sp)]!='-')[gain,c('TRUE','FALSE')],1,my.binom.test)
b=barplot(te.alt[1,],ylim=range(0,te.alt),border=NA,las=3,ylab='proportion of TE',main='TE in alternificatin')
segments(b,te.alt[2,],b,te.alt[3,])

dev.off()


# look on TE types ####
# h2te = getAnnOverlap(h$seg,data.frame(chr_id=rm$genoName,start=rm$genoStart,stop=rm$genoEnd),T)$overlap
# h2te = as.data.frame(h2te)
# h2te = cbind(h2te,seg_id=rownames(h$seg)[h2te$from],h$seg[h2te$from,c('sites','type','length')])
# h2te$overlap = pmin(h$seg$stop[h2te$from],rm$genoEnd[h2te$to])-pmax(h$seg$start[h2te$from],rm$genoStart[h2te$to])+1
# h2te$te.class = rm$repClass[h2te$to]
# saveRDS(h2te,'Rdata/human.seg2TE.Rdata')
rownames(h2te) = NULL
range(h2te$overlap)
hist(h2te$overlap[h2te$sites=='ad'],1:max(h2te$overlap),xlim=c(0,60))
hist((h2te$overlap/h2te$length)[h2te$sites=='ad'],seq(0,1,length.out = 100))
table(h2te$overlap/h2te$length>0.1,h2te$sites)

table(h2te$te.class[h2te$sites=='ad'],h2te$overlap[h2te$sites=='ad']==1)
h2te[h2te$overlap==1 & h$seg$sites[h2te$from]=='ad',][1:2,]
h$seg[4525,]

rm[3299241,]

h2te[1:10,]
dim(h2te)
range(h2te[,2])
dim(h$seg)
table(table(h2te[h2te$overlap>10,1]))

# more than one overlaping TEs
t = table(h2te$seg_id[h2te$sites=='ad'])
table(t)
h2te[h2te$seg_id==names(t)[t>1][2],]
rm[c(3414532,3447700),]

te.sids= list(h     = born.seg.ids[sp.birth=='h',1],
						  hq    = born.seg.ids[sp.birth=='hq',1],
							#hqmrb = born.seg.ids[sp.birth=='hqmrb',1], #there are too few of them
						  alt=rownames(h$seg)[h$seg$sites=='ad' & h$seg$type=='ALT'],
						  exn=rownames(h$seg)[h$seg$sites=='ad' & h$seg$type=='EXN'])
sapply(te.sids,length)
te.sids = lapply(te.sids,intersect,y=h2te$seg_id)

f = function(s2c,sids,classes = unique(s2c$te.class)){
	sids = intersect(sids,s2c$seg_id)
	f = s2c$seg_id %in% sids
	r = table(factor(unlist(lapply(split(s2c$te.class[f],s2c$seg_id[f]),unique)),levels=classes))
	cbind(t(sapply(r,function(x)my.binom.test(x,length(sids)-x))),r)
}

ff = function(s2c,sids,classes = unique(s2c$te.class)){
	r = array(NA,dim=c(length(classes),4,length(sids)),dimnames=list(classes,c('p','ci.lower','ci.high','n'),names(sids)))
	for(n in names(sids))
		r[,,n] = f(s2c,sids[[n]],classes)
	r
}

cl = unique(unlist(seg2teclass))
z = ff(h2te,te.sids,cl)
z = z[order(apply(z[,'n',],1,sum),decreasing = T),,]

z10 = ff(h2te[h2te$overlap>9,],te.sids,cl)
z10 = z10[order(apply(z10[,'n',],1,sum),decreasing = T),,]


pdf('figures/newborn/TE.class.propo.pdf',w=8,h=6)
par(mfrow=c(1,2),tck=-0.02,mgp=c(1.3,0.2,0),mar=c(7.5,2.5,1.5,0),oma=c(0,0,0,1))
b=barplot(t(z[1:6,'p',]),beside = T,las=3,legend.text = T,ylab='exons_with_class/exons_with_TE',ylim=range(0,z[1:6,1:3,]),main='All TE overlaps')
segments(b,t(z[1:6,'ci.lower',]),b,t(z[1:6,'ci.high',]))

b=barplot(t(z10[1:6,'p',]),beside = T,las=3,legend.text = T,ylab='exons_with_class/exons_with_TE',ylim=range(0,z10[1:6,1:3,]),main='TE overlaps >= 10nt')
segments(b,t(z10[1:6,'ci.lower',]),b,t(z10[1:6,'ci.high',]))
dev.off()

# check whether lost exons are alternative #####
gain = c(species$short,'hq','mr','mrb','hqmrb','hqmrbo')
loss = sapply(gain,function(x)gsub(x,'','hqmrboc'))
altn = sapply(exon.birth.one,function(x)sum(x$type=='ALT',na.rm=T))


table(nchar(obs.sp),altn)
f = filters$orth.adj.bl.int.len.psi.n
table(nchar(obs.sp[f]),altn[f])

pdf('figures/newborn/alt.exons.in.lost-n-gained.exons.pdf',w=8,h=32)
par(mfrow=c(8,2),tck=-0.02,mgp=c(3.5,0.7,0),mar=c(2.5,4.5,1.5,0),oma=c(0,0,0,1))
for(n in names(filters)){
	f = filters[[n]]
	t = table(altn[f],obs.sp[f])[,rev(gain)]
	imageWithText(sweep(t,2,apply(t,2,sum),'/'),t,names.as.labs = T,xlab='# of ALT species',ylab='Observed species',main=paste0('Gain, ',n,' (',sum(f),')'))
	t = table(altn[f],obs.sp[f])[,rev(loss)]
	imageWithText(sweep(t,2,apply(t,2,sum),'/'),t,names.as.labs = T,xlab='# of ALT species',ylab='Observed species',main=paste0('Loss, ',n,' (',sum(f),')'))
}
dev.off()