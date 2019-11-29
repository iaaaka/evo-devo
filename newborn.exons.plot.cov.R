options(stringsAsFactors = FALSE)
setwd('~/projects/evo.devo')
source('code/r.functions/load.all.data.F.R')
source('../../r.code/util.R')
library(SAJR)
#library(seqinr)
library(reshape)
library(plyr)
library(doMC)
library(GenomicAlignments)

species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
#all.anns = readRDS('Rdata/all.anns.Rdata')
exon.birth.one = readRDS('Rdata/exon.birth.one.Rdata')
born.exn.sajr = readRDS('Rdata/born.exn.sajr.Rdata')
gene.descrs = readRDS('Rdata/ens.gene.descr.Rdata')
my.ge = readRDS('Rdata/my.ge.Rdata')

bl = read.table('processed/new.born.exons/same.exon.out.final.tab',sep='\t')
colnames(bl) = c('qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','qcovhsp','qcovs','qlen')
bl$map.frac = (bl$qend-bl$qstart+1)/bl$qlen

bl = bl[bl$map.frac> 0.5,]
bl$qid = sapply(strsplit(bl$qseqid,'.',T),'[',3)
bl = split(bl,bl$qid)

blast.sp = sapply(bl,function(x){paste(species$short[rownames(species) %in% sapply(strsplit(x$sseqid,'.',T),'[',1)],collapse = '')})
sort(table(blast.sp))

obs.sp = sapply(exon.birth.one,function(x)paste(species$short[rownames(species) %in% rownames(x)[!is.na(x$seg_id)]],collapse = ''))
blast.sp = blast.sp[names(obs.sp)]

# plot nb exons #####
seg2ens = readRDS('Rdata/seg2ens.Rdata')

sids = unique(unlist(lapply(exon.birth.one,function(x)x[,c('seg_id','useg_id','dseg_id')])))
#all.anns. = lapply(all.anns,function(x)x[rownames(x) %in% sids,])
#rm(all.anns)
#gc()
#saveRDS(all.anns.,'Rdata/anns.subset.for.newborn.cov.plot.Rdata',version = 2)
all.anns. = readRDS('Rdata/anns.subset.for.newborn.cov.plot.Rdata')
#all.bams = paste0('processed/mapping/hisat2.s/',meta$species,'/',meta$fname,'.bam')

#to restart
fls = list.files('figures/newborn.exon.cov/')
fls = sapply(strsplit(fls,'.',T),'[',2)
fls = setdiff(names(exon.birth.one),fls)
fls = match(fls,names(exon.birth.one))

registerDoMC(16)
l_ply(fls,function(i){
	print(i)
	#cat('\n',i)
	pdf(paste0('figures/newborn.exon.cov/',obs.sp[i],'-',blast.sp[i],'.',names(exon.birth.one)[i],'.pdf'),w=12,h=15)
	nb = exon.birth.one[[i]]
	par(tck=-0.01,mgp=c(1.1,0.2,0),mar=c(2.5,2.5,3.2,0),oma=c(0,0,0,1))
	layout(matrix(1:21,ncol=3,byrow = T),widths = c(5,1,1))
	cov = list()
	for(s in rownames(species)){
		cat(' ',s)
		bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s],'.bam')
		#bams = paste0('processed/mapping/hisat2.s/',s,'/',meta$fname[meta$species==s & meta$tissue=='testis' & meta$stage=='16wpc'],'.bam')
		#bams = all.bams[meta$species==s]
		
		if(all.anns.[[s]][nb[s,'useg_id'],'strand'] == 1){
			start = all.anns.[[s]][nb[s,'useg_id'],'start']
			stop = all.anns.[[s]][nb[s,'dseg_id'],'stop']
		}else{
			start = all.anns.[[s]][nb[s,'dseg_id'],'start']
			stop = all.anns.[[s]][nb[s,'useg_id'],'stop']
		}
		ensids = seg2ens[[s]][[nb[s,'useg_id']]]
		if(stop+100 <= 536870912 && start - 100 <= 536870912){
			cov[[s]]=getReadCoverage(bams,all.anns.[[s]][nb[s,'useg_id'],'chr_id'],start - 100,stop+100,-all.anns.[[s]][nb[s,'useg_id'],'strand'])
			maxcov = max(cov[[s]]$cov,cov[[s]]$juncs$score)
			plotReadCov(cov[[s]],reverse = (all.anns.[[s]][nb[s,'useg_id'],'strand'] == -1),min.junc.cov = 2,bty='n',xlab=paste0('Chr ',all.anns.[[s]][nb[s,'useg_id'],'chr_id']),ylab='Reads',
									main=paste0(s,'\n',paste(ensids,collapse=', '),'\n',paste(gene.descrs[[s]][ensids[ensids %in% rownames(gene.descrs[[s]])],2],collapse=', ')),ylim=c(-0.2*maxcov,maxcov))
			abline(h = -0.12*maxcov)
			seg = all.anns.[[s]][nb[s,'useg_id'],]
			rect(seg$start,-0.2*maxcov,seg$stop,-0.04*maxcov,col = 'black',border = NA)
			seg = all.anns.[[s]][nb[s,'dseg_id'],]
			rect(seg$start,-0.2*maxcov,seg$stop,-0.04*maxcov,col = 'black',border = NA)
		}else
			plot.new()
		
		if(!is.na(nb[s,'seg_id'])){
			#abline(v=c(nb[s,'start'],nb[s,'stop']),lty=3,col='blue')
			rect(nb[s,'start'],-0.2*maxcov,nb[s,'stop'],-0.04*maxcov,col = 'blue',border = NA)
		}	
		if(!is.na(nb[s,'seg_id']) && sum(!is.na(born.exn.sajr[[s]]$ir[nb[s,'seg_id'],colnames(born.exn.sajr[[s]]$ir) %in% rownames(meta)]))>0)
			plotTissueAgeProile(born.exn.sajr[[s]]$ir[nb[s,'seg_id'],],meta,age.axis = 'rank',ylab='PSI',main='AS')
		else
			plot.new()
		plotTissueAgeProile(my.ge[[s]]$rpkm[nb[s,'gene_id'],],meta,age.axis = 'rank',ylab='RPKM',main='GE')
	}
	dev.off()
	saveRDS(cov,paste0('Rdata/newborn.exon.cov/',obs.sp[i],'-',blast.sp[i],'.',names(exon.birth.one)[i],'.Rdata'))
	gc()
},.parallel = T)


# check time ####
# param = ScanBamParam(flag=scanBamFlag(isMinusStrand=FALSE),which=GRanges('5', IRanges(16226953, 16227875)))
# b = readGAlignments(bams,param = param)
# system.time({sapply(1:1000,function(i)ScanBamParam(flag=scanBamFlag(isMinusStrand=FALSE),which=GRanges('5', IRanges(16226953, 16227875))))})
# #0.01s
# system.time({sapply(1:100,function(i)readGAlignments(bams,param = param))})
#0.3s on laptop, 0.06 on cluster
