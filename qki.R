library(openxlsx)
library(SAJR)
source('code/r.functions/load.all.data.F.R')
source('~/skoltech/r.code/util.R')
source('code/r.functions/paper.figures.5.F.R')

species = readRDS('Rdata/species.Rdata')
ens.ge.cod = readRDS('Rdata/ens.ge.cod.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
all.anns = readRDS('Rdata/all.anns.Rdata')
seg2ens = readRDS('Rdata/seg2ens.Rdata')
gene.descrs = readRDS('Rdata/ens.gene.descr.Rdata')
age.segs = readRDS('Rdata/devAS.4patt.Rdata')
anns = readRDS('Rdata/anns.Rdata')
psi.tsm = readRDS('Rdata/psi.tsm.Rdata')


# QKI
qkis = read.csv('processed/mapping/QKI.RNAi/samples.csv')
qkis$id = paste0(substr(qkis$cell.line,1,1),rep(c('i','i','c','c'),times=2),rep(1:2,times=4))
qkis$full.names = paste0(qkis$cell.line,'-',rep(c('RNAi','RNAi','ctrl','ctrl'),times=2),'-',rep(1:2,times=4))
qkii.ge = readRDS('Rdata/qki/qkii.ge.Rdata')
qkii.as = readRDS('Rdata/qki/qkii.as.Rdata')
hqki.ints = readRDS('Rdata/qki/hqki.ints.Rdata')
mqki.ints = readRDS('Rdata/qki/mqki.ints.Rdata')
qkii.pv = readRDS('Rdata/qki/qki.pv.Rdata')
qki.isoforms = c('qki-5','qki-6','qki-7')

# try targets from papers ####
# qki.mus.targ = read.xlsx('input/QKI/Quaking and PTB control overlapping splicing regulatory networks during muscle cell differentiation suppl.xlsx',2)
# qki.mus.targ6 = read.xlsx('input/QKI/Quaking and PTB control overlapping splicing regulatory networks during muscle cell differentiation suppl.xlsx',6)
# 
# table(qki.mus.targ$eventName %in% qki.mus.targ6$eventName)
# qki.mus.targ$cooors = qki.mus.targ6$eventPosition[match(qki.mus.targ$eventName,qki.mus.targ6$eventName)]
# rm(qki.mus.targ6)
# qki.mus.targ[1:2,]
# t = strsplit(qki.mus.targ$cooors,'[:-]')
# qki.mus.targ$chr = gsub('chr','',sapply(t,'[',1))
# qki.mus.targ$start = as.numeric(sapply(t,'[',2))
# qki.mus.targ$stop = as.numeric(sapply(t,'[',3))
# qki.mus.targ$strand = ifelse(sapply(t,'[',4)=='+',1,-1)
# hist(qki.mus.targ$stop-qki.mus.targ$start)
# 
# sids = list()
# m = all.anns$mouse[all.anns$mouse$sites=='ad',]
# for(i in 1:nrow(qki.mus.targ)){
# 	cat('\r',i,'   ')
# 	sids[[i]] = rownames(m)[m$chr_id == qki.mus.targ$chr[i] & 
# 													m$strand == qki.mus.targ$strand[i] & 
# 													m$start >= qki.mus.targ$start[i] & 
# 													m$stop <= qki.mus.targ$stop[i]]
# }
# table(sapply(sids,length))

# look on ENCODE RNAi #####
# qkis$id = paste0(substr(qkis$cell.line,1,1),rep(c('i','i','c','c'),times=2),rep(1:2,times=4))
# qkii.ge = loadGData('processed/annotation/all.species/merged/human.sajr',paste0('processed/mapping/QKI.RNAi/sajr/',qkis$sra),qkis$id)
# qkii.ge$cpm = sweep(qkii.ge$cnts,2,apply(qkii.ge$cnts,2,sum),'/')*1e6
# qkii.ge[1:2,]
# gene.descrs$human[gene.descrs$human$gene.name=='QKI',] # ENSG00000112531
# qki.segs = names(seg2ens$human)[sapply(seg2ens$human,function(x) 'ENSG00000112531' %in% x)] # hum.57513
# saveRDS(qkii.ge,'Rdata/qki/qkii.ge.Rdata')

# qkii = loadSAData('processed/annotation/all.species/merged/human.sajr',paste0('processed/mapping/QKI.RNAi/sajr/',qkis$sra),qkis$id)
# table(rownames(qkii$seg) == rownames(all.anns$human))
# qkii$seg = all.anns$human
# qkii$ir[qkii$i+qkii$e < 10] = NA
# table(apply(!is.na(qkii$ir),1,sum),qkii$seg$type)
# table(apply(!is.na(qkii$ir[,1:4]),1,sum),qkii$seg$type)
# saveRDS(qkii,'Rdata/qki/qkii.as.Rdata')

hqki.ints.eclip=loadIntronCounts(paste0('processed/mapping/QKI.RNAi/sajr/',qkis$sra,'.intron'),qkis$id,'6:163984752-')
cbind(hqki.ints.eclip$intron,apply(hqki.ints.eclip$cnt,1,sum))

hqki.ints.eclip$intron$name = '-'
hqki.ints.eclip$intron['6:163984752-163987752:1','name'] = 'qki-5'
hqki.ints.eclip$intron['6:163984752-163986977:1','name'] = 'qki-6'
hqki.ints.eclip$intron['6:163984752-163985698:1','name'] = 'qki-7'
hqki.ints.eclip$freq = sweep(hqki.ints.eclip$cnt,2,apply(hqki.ints.eclip$cnt,2,sum),'/')


# plot MDS and QKI expression #3333
f = qkii.as$seg$sites=='ad' & qkii.as$seg$type=='ALT' & apply(!is.na(qkii.as$ir),1,sum) > 3
table(f)

pdf('figures/paper.figures/6/nature.review/qki.encode/qki.encode.overview.pdf',w=7,h=7)
par(mfrow=c(2,2),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(5,2.5,1.5,0),oma=c(0,0,0,1))
o = c(3,4,1,2,7,8,5,6)
b = barplot(qkii.ge$cpm['hum.57513',o],bty='n',main='QKI expression',ylab='CPM',names.arg = '',las=3,col=c('gray','gray','orange','orange'),border = NA)
text(b,-10,qkis$full.names[o],adj = c(0,0.5),srt=-45,xpd=NA)
mds=cmdscale(1-cor(qkii.ge$cnts,m='sp'),k=2)
mds = mds+runif(8,-0.03,0.03)
plot(mds*1.05,t='n',bty='n',xlab='Dim 1',ylab='Dim 2',main='GE MDS')
text(mds,rownames(mds),col=c('orange','orange','gray','gray'),font=2,cex=2,xpd=NA)

mds=cmdscale(1-cor(qkii.as$ir[f,],m='p',u='p'),k=2)
plot(mds*1.05,t='n',bty='n',xlab='Dim 1',ylab='Dim 2',main='AS MDS')
text(mds,rownames(mds),col=c('orange','orange','gray','gray'),font=2,cex=2,xpd=NA)
par(mar=c(5,2.5,1.5,5))
b=barplot(hqki.ints.eclip$freq[match(qki.isoforms,hqki.ints.eclip$intron$name),o],xaxt='n',border = NA,ylab='PSI',legend.text = qki.isoforms,args.legend = list(bty='n',max(b)*1.01,1,xpd=TRUE,xjust=0))
text(b,-0.01,qkis$full.names[o],adj = c(0,0.5),srt=-45,xpd=T)
dev.off()

# Make test #######
# hpv = fitSAGLM(qkii.as[f,1:4],formula = x ~ type,terms = qkis[1:4,],return.pv = TRUE)
# kpv = fitSAGLM(qkii.as[f,5:8],formula = x ~ type,terms = qkis[5:8,],return.pv = TRUE)
# #kpv. = fitSAGLM(qkii[f,c(5,7,6,8)],formula = x ~ type,terms = qkis[5:8,],return.pv = TRUE)

# hpv = cbind(hpv[,1:2] ,qv = p.adjust( hpv[,2],m='BH'),dpsi=apply(qkii.as$ir[f,1:2]-qkii.as$ir[f,3:4],1,mean))
# kpv = cbind(kpv[,1:2] ,qv = p.adjust( kpv[,2],m='BH'),dpsi=apply(qkii.as$ir[f,5:6]-qkii.as$ir[f,7:8],1,mean))
# #kpv.= cbind(kpv.[,1:2],qv = p.adjust(kpv.[,2],m='BH'),dpsi=apply(qkii$ir[f,c(5,7)]-qkii$ir[f,c(6,8)],1,mean,na.rm=T))
# qkii.pv = data.frame(hpv=hpv[,2],hqv=hpv[,3],hdpsi=hpv[,4],
# 										 kpv=kpv[,2],kqv=kpv[,3],kdpsi=kpv[,4])
# 
# z = apply(is.na(qkii.as$ir[f,1:4]),1,sum)>0
# qkii.pv$hpv[z] = qkii.pv$hdspi[z] = NA
# z = apply(is.na(qkii.as$ir[f,5:8]),1,sum)>0
# qkii.pv$kpv[z] = qkii.pv$kdspi[z] = NA
# qkii.pv$hqv = p.adjust(qkii.pv$hpv,m='BH')
# qkii.pv$kqv = p.adjust(qkii.pv$kpv,m='BH')
# 
# saveRDS(qkii.pv,'Rdata/qki/qki.pv.Rdata')

# look on results ####
png('figures/paper.figures/6/nature.review/qki.encode/qki.hepg2-vs-k562.png',width = 6,height = 6,units = 'in',res=300)
par(mfrow=c(1,1),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
plot2Sign(qkii.pv$hqv<0.05,qkii.pv$kqv<0.05,qkii.pv$hdpsi,qkii.pv$kdpsi,'HepG2','K562',bty='n')
dev.off()

# _actaac #####
fa = readRDS('Rdata/ad.alt.fa.Rdata')$human
fa = fa[names(fa) %in% rownames(qkii.pv)]
fa200 = sapply(fa,function(x)paste0(substr(x,1,200),' ',substr(x,nchar(x)-199,nchar(x))))
table(nchar(fa200))
actaac = gregexpr('actaac',fa200)
actaa  = gregexpr('actaa' ,fa200)
names(actaa) = names(actaac) = names(fa200)

# _eCLIP#####
heclip1 = getEClipSiteOccurence('processed/QKI.eCLIP/ENCFF454GLW.bed',ann = all.anns$human[rownames(qkii.pv),])
heclip2 = getEClipSiteOccurence('processed/QKI.eCLIP/ENCFF544QKV.bed',ann = all.anns$human[rownames(qkii.pv),])
keclip1 = getEClipSiteOccurence('processed/QKI.eCLIP/466_01.basedon_466_01.peaks.l2inputnormnew.bed.compressed.bed.narrowPeak.bed',ann = all.anns$human[rownames(qkii.pv),])
keclip2 = getEClipSiteOccurence('processed/QKI.eCLIP/466_02.basedon_466_02.peaks.l2inputnormnew.bed.compressed.bed.narrowPeak.bed',ann = all.anns$human[rownames(qkii.pv),])


pdf('figures/paper.figures/6/nature.review/qki.encode/qki.rnai.pdf',w=12,h=4)
par(mfrow=c(2,6),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(5,2.5,1.5,0),oma=c(0,0,1,6))

f = function(e,p,qv.thr = 0.05){
	e = cbind(e[,1:200],NA,e[,201:400])
	f =!is.na(p$qv) & !is.na(p$dpsi)
	list(n = apply(e[f & p$qv > qv.thr,],2,mean),
			 u = apply(e[f & p$qv < qv.thr & p$dpsi > 0,],2,mean),
			 d = apply(e[f & p$qv < qv.thr & p$dpsi < 0,],2,mean))
}

f0 = function(x,l){
	v = unlist(x)
	v = v[v != -1]
	r = rep(0,400)
	for(i in v)
		r[i:(i+l-1)] = r[i:(i+l-1)]+1
	r = c(r[1:200],NA,r[201:400])
	r/length(x)
}

f1 = function(o,p,l,qv.thr = 0.05){
	pv = pv[names(o),]
	f =!is.na(p$qv) & !is.na(p$dpsi)
	
	list(n = f0(o[f & p$qv > qv.thr],l),
			 u = f0(o[f & p$qv < qv.thr & p$dpsi > 0],l),
			 d = f0(o[f & p$qv < qv.thr & p$dpsi < 0],l))
}

for(qv.thr in c(0.05,0.2,0.5)[1]){
	for(c in c('HepG2','K562')){
		pv = qkii.pv[,substr(colnames(qkii.pv),1,1) == tolower(substr(c,1,1))]
		colnames(pv) = substr(colnames(pv),2,100)
		#ecl  = list(HepG2=hqki.tsm,K562=mqki.tsm)[[c]]
		e1 = f(get0(paste0(tolower(substr(c,1,1)),'eclip1')),pv,qv.thr = qv.thr)
		e2 = f(get0(paste0(tolower(substr(c,1,1)),'eclip2')),pv,qv.thr = qv.thr)
		
		plot.new()
		text(grconvertX(0.5,'nfc','user'),grconvertY(0.5,'nfc','user'),c,cex=2)
		barplot(c(up=sum(pv$qv<qv.thr & pv$dpsi>0,na.rm = TRUE),down=sum(pv$qv<qv.thr & pv$dpsi<0,na.rm = TRUE)),border=NA,col=c('red','blue'),ylab='# of exons')
		
		# motif
		x = c(1:200,NA,301:500)
		motif.occ = f1(actaa,pv,qv.thr = qv.thr,l = 5)
		plot(x,motif.occ$n,t='l',col='gray',lwd=3,xaxt='n',ylab='ACTAA freq',main='ACTAA',bty='n',ylim=range(0,unlist(motif.occ),na.rm=T),xlab='')
		y = grconvertY(c(-0.01,-0.05,-0.09),'npc','user')
		segments(1,y[2],450,y[2],xpd=T)
		rect(201,y[1],299,y[3],col='black',xpd=T)
		lines(x,motif.occ$d,col='blue',lwd=3)
		lines(x,motif.occ$u,col='red',lwd=3)
		lines(x,motif.occ$n,col='gray',lwd=3)
		
		# eCLIP
		plot(x,e1$n,t='l',col='gray',lwd=3,xaxt='n',ylab='QKI binding freq',main='eCLIP',bty='n',ylim=range(0,unlist(c(e1,e2)),na.rm=T),xlab='')
		y = grconvertY(c(-0.01,-0.05,-0.09),'npc','user')
		segments(1,y[2],450,y[2],xpd=T)
		rect(201,y[1],299,y[3],col='black',xpd=T)
		lines(x,e1$u,col='red',lwd=3)
		lines(x,e2$u,col='red',lwd=3)
		lines(x,e1$d,col='blue',lwd=3)
		lines(x,e2$d,col='blue',lwd=3)
		lines(x,e2$n,col='gray',lwd=3)
		lines(x,e1$n,col='gray',lwd=3)
		
		
		cmn = intersect(rownames(pv)[!is.na(pv$qv+pv$dpsi) & pv$qv<qv.thr],rownames(age.segs$human))
		b = table(age.segs$human[cmn,'brain'],eclip=sign(pv[cmn,'dpsi']))[c('d','u'),c('-1','1')]
		h = table(age.segs$human[cmn,'heart'],sign(pv[cmn,'dpsi']))[c('d','u'),c('-1','1')]
		colnames(b)=rownames(b)=colnames(h)=rownames(h)=c('down','up')
		
		#table(age.segs$human[f,'brain'],age.segs$human[f,'heart'])[c('u','d'),c('u','d')]
		#table(age.segs$human[cmn,'brain'],age.segs$human[cmn,'heart'])[c('u','d'),c('u','d')]
		
		
		bft = fisher.test(b)
		hft = fisher.test(h)
		par(mar=c(5,2.5,2.5,0))
		imageWithText(b,xlab='devAS',ylab='RNAi',main=paste0('brain\nodd=',round(bft$estimate,2),' pv=',format(bft$p.value,digits = 1)))
		imageWithText(h,xlab='devAS',ylab='RNAi',main=paste0('heart\nodd=',round(hft$estimate,2),' pv=',format(hft$p.value,digits = 1)))	
	}
	mtext(paste0('FDR < ',qv.thr),3,outer = T)
}
dev.off()



# _QKI.isoforms ######
o = read.csv('input/ens.orths.txt.gz')
qkio = o[o[,1]=='ENSG00000112531',]
colnames(qkio ) = rownames(species)
gene.descrs$rat[qkio$Rat.Ensembl.Gene.ID,]

# __human ######
# hens = loadSAData('processed/annotation/all.species/ensambl/Homo_sapiens.GRCh37.73.sajr.gz')
# hens = setSplSiteTypes(hens,'processed/annotation/all.species/ensambl/Homo_sapiens.GRCh37.73.sajr.gz')
# e = hens$seg[hens$seg$gene_id=='ENSG00000112531',]
# 
# hens=loadEnsGTF('processed/annotation/all.species/ensambl/Homo_sapiens.GRCh37.73.gtf.gz')
# hens[hens$gene_id=='ENSG00000112531',]
# plotTranscripts(hens[hens$gene_id=='ENSG00000112531' & hens$type=='protein_coding',],bty='n',xlim=c(163984000,163990000),cds.col='red')
# t = hens[hens$gene_id=='ENSG00000112531' & hens$type=='protein_coding',]
# t[t$feature=='exon' & t$start>163984000 & t$start<163990000,c(4,5,12)]

# m = meta[meta$species=='human',]
# hqki.ints=loadIntronCounts(paste0('processed/sajr/uniq/',m$fname,'.intron'),rownames(m),'6:163984752-')
# cbind(hqki.ints$intron,apply(hqki.ints$cnt,1,sum))
# 
# hqki.ints$intron$name = '-'
# hqki.ints$intron['6:163984752-163987752:1','name'] = 'qki-5'
# hqki.ints$intron['6:163984752-163986977:1','name'] = 'qki-6'
# hqki.ints$intron['6:163984752-163985698:1','name'] = 'qki-7'
# 
# hqki.ints$freq = sweep(hqki.ints$cnt,2,apply(hqki.ints$cnt,2,sum),'/')
# saveRDS(hqki.ints,'Rdata/qki/hqki.ints.Rdata')

m = meta[colnames(hqki.ints$cnt),]
hqki.tsm = calcMeanCols(hqki.ints$cnt,paste(m$species,m$tissue,m$stage),sum)
hqki.tsm = sweep(hqki.tsm,2,apply(hqki.tsm,2,sum),'/')

# mouse 
# mens=loadEnsGTF('processed/annotation/all.species/ensambl/Mus_musculus.GRCm38.84.gtf.gz')
# mens[1:2,]
# 
# mens$transcript_id[mens$gene_id=='ENSMUSG00000062078',]
# 
# plotTranscripts(mens[mens$gene_id=='ENSMUSG00000062078',],bty='n',cds.col='red',xlim=c(10220000,10205000))
# t = mens[mens$gene_id=='ENSMUSG00000062078',]
# t[t$feature=='exon' & t$start>10205000 & t$start<10220000,c(4,5,16)]

# __mouse #####
# m = meta[meta$species=='mouse',]
# mqki.ints=loadIntronCounts(paste0('processed/sajr/uniq/',m$fname,'.intron'),rownames(m),'17:\\d+-10216159:-1')
# cbind(mqki.ints$intron,apply(mqki.ints$cnt,1,sum))
# 
# mqki.ints$intron$name = '-'
# mqki.ints$intron['17:10213434-10216159:-1','name'] = 'qki-5'
# mqki.ints$intron['17:10214195-10216159:-1','name'] = 'qki-6'
# mqki.ints$intron['17:10215475-10216159:-1','name'] = 'qki-7'

# mqki.ints$freq = sweep(mqki.ints$cnt,2,apply(mqki.ints$cnt,2,sum),'/')
# saveRDS(mqki.ints,'Rdata/qki/mqki.ints.Rdata')

m = meta[colnames(mqki.ints$cnt),]
mqki.tsm = calcMeanCols(mqki.ints$cnt,paste(m$species,m$tissue,m$stage),sum)
mqki.tsm = sweep(mqki.tsm,2,apply(mqki.tsm,2,sum),'/')


pdf('figures/paper.figures/6/nature.review/qki.encode/qki.isoform.expression.pdf',w=18,h=6)
par(mfrow=c(2,7),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(5,2.5,1.5,0),oma=c(0,0,0,6))

for(s in c('human','mouse')){
	ints = list(human=hqki.ints,mouse=mqki.ints)[[s]]
	tsm  = list(human=hqki.tsm,mouse=mqki.tsm)[[s]]
	plot.new()
	plotPNG(paste0('figures/paper.figures/5/icons/',s,'.png'),0.5,0.5,0.9)
	plotTissueAgeProile(ens.ge.cod[[s]]$rpkm[qkio[1,s],],meta,age.axis = 'rank',main='GE',ylab='RPKM',bty='n')
	for(i in qki.isoforms)
		plotTissueAgeProile(ints$freq[ints$intron$name==i,],meta,age.axis = 'rank',main=i,bty='n',ylab='PSI',ylim=c(0,1))
	z = meta.tsm[meta.tsm$species==s & meta.tsm$tissue=='brain',]
	z = z[order(z$days),]
	b=barplot(tsm[match(qki.isoforms,ints$intron$name),rownames(z)],names.arg = z$stage,las=3,main='Brain')
	
	z = meta.tsm[meta.tsm$species==s & meta.tsm$tissue=='heart',]
	z = z[order(z$days),]
	barplot(tsm[match(qki.isoforms,ints$intron$name),rownames(z)],names.arg = z$stage,las=3,main='heart',legend.text = qki.isoforms,ylim=c(0,1),args.legend = list(max(b)*1.01,1,xjust=0,xpd=NA,bty='n'))
}
dev.off()


# QKI-targets ####
ff = function(d)paste0(ifelse(apply(d[,1:200],1,sum)>0,'u',''),ifelse(apply(d[,201:400],1,sum)>0,'d',''))
qki.targ = data.frame(heclip1 = ff(getEClipSiteOccurence('processed/QKI.eCLIP/ENCFF454GLW.bed',ann = anns$human[anns$human$sites=='ad',])),
											heclip2 = ff(getEClipSiteOccurence('processed/QKI.eCLIP/ENCFF544QKV.bed',ann = anns$human[anns$human$sites=='ad',])),
											keclip1 = ff(getEClipSiteOccurence('processed/QKI.eCLIP/466_01.basedon_466_01.peaks.l2inputnormnew.bed.compressed.bed.narrowPeak.bed',ann = anns$human[anns$human$sites=='ad',])),
											keclip2 = ff(getEClipSiteOccurence('processed/QKI.eCLIP/466_02.basedon_466_02.peaks.l2inputnormnew.bed.compressed.bed.narrowPeak.bed',ann = anns$human[anns$human$sites=='ad',])))

rownames(qki.targ) = rownames(anns$human)[anns$human$sites=='ad']
table(qki.targ$heclip1=='',qki.targ$heclip2=='')
table(qki.targ$keclip1=='',qki.targ$keclip2=='')
table(qki.targ$keclip1=='',qki.targ$heclip1=='')

fa = readRDS('Rdata/ad.alt.fa.Rdata')$human
fa = fa[names(fa) %in% rownames(qki.targ)]

qki.targ$actaac = qki.targ$actaa = NA
qki.targ[names(fa),'actaac'] = paste0(ifelse(grepl('actaac',substr(fa,1,200)),'u',''),
																			ifelse(grepl('actaac',substr(fa,nchar(fa)-199,nchar(fa))),'d',''))
qki.targ[names(fa),'actaa'] = paste0(ifelse(grepl('actaa',substr(fa,1,200)),'u',''),
																		 ifelse(grepl('actaa',substr(fa,nchar(fa)-199,nchar(fa))),'d',''))

qki.targ$hepg2 = qki.targ$k562 = '-'
table(pv=is.na(qkii.pv$kqv) ,is.na(qkii.pv$kdpsi))

f = !is.na(qkii.pv$hqv) & rownames(qkii.pv) %in% rownames(qki.targ)
qki.targ[rownames(qkii.pv)[f],'hepg2'] = ifelse(qkii.pv$hqv[f]<0.05,ifelse(qkii.pv$hdpsi[f]>0,'u','d'),'n')
f = !is.na(qkii.pv$kqv) & rownames(qkii.pv) %in% rownames(qki.targ)
qki.targ[rownames(qkii.pv)[f],'k562'] = ifelse(qkii.pv$kqv[f]<0.05,ifelse(qkii.pv$kdpsi[f]>0,'u','d'),'n')
table(qki.targ$k562,qki.targ$hepg2)

fisher.test(table(qki.targ$actaac,qki.targ$hepg2)[c('d','u'),c('d','u')])
f = qki.targ$hepg2 !='-'
fisher.test(table(qki.targ$actaac[f]!='',qki.targ$hepg2[f]!='n'))


table(qki.targ$actaac,qki.targ$heclip1,useNA='always')

fisher.test(table(qki.targ$actaac=='',qki.targ$heclip1==''))

dpat = age.segs$human[anns$human$sites=='ad',]
o = c('-','n','u','d','ud','du')
table(dpat[,'brain'],dpat[,'heart'])[o,o]

table(dpat[,'brain'],dpat[,'heart'],qki.targ$actaac)[o,o,]
table(dpat[,'brain'],qki.targ$actaac)[o,]
table(dpat[,'heart'],qki.targ$actaac)[o,]

# are QKI-targets more likely to be expressed in both brain and heart
table(qki.targ$hepg2,dpat[,'brain'])[c('d','u'),c('d','u')]
table(qki.targ$hepg2,dpat[,'heart'])[c('d','u'),c('d','u')]

table(dpat[,'brain']!='-',qki.targ$actaac!='')
f = dpat[,'brain']!='-' | dpat[,'heart']!='-'
fisher.test(table(dpat[f,'brain']!='-' & dpat[f,'heart']!='-',qki.targ$actaac[f]!='')) # eCLIP is assotiated with expression in both tissues, but not the motif

f = qki.targ$hepg2 != '-'
fisher.test(table(dpat[f,'heart']!='-',qki.targ$hepg2[f]!='n'))

f = (dpat[,'brain']!='-' | dpat[,'heart']!='-') & qki.targ$hepg2 != '-'
fisher.test(table(dpat[f,'brain']!='-' & dpat[f,'heart']!='-',qki.targ$hepg2[f]!='n')) # RNAi is assotiated with expression in both tissues, but not the motif


# check whether QKI-targets are more likely to be devAS in both tissues
f = dpat[,'brain']!='-' & dpat[,'heart']!='-' & (dpat[,'brain']!='n' | dpat[,'heart']!='n')
fisher.test(table(dpat[f,'brain']!='n' & dpat[f,'heart']!='n',qki.targ$actaa[f]!=''))

f = dpat[,'heart']!='-' & qki.targ$hepg2 != '-'
fisher.test(table(dpat[f,'heart']!='n',qki.targ$hepg2[f]!='n'))

f = dpat[,'brain']!='-' & dpat[,'heart']!='-' & (dpat[,'brain']!='n' | dpat[,'heart']!='n') & qki.targ$hepg2 != '-'
fisher.test(table(dpat[f,'brain']!='n' & dpat[f,'heart']!='n',qki.targ$hepg2[f]!='n'))

o = c('-','n','u','d')
o = c('-','n','u','d','ud','du')
table(brain=dpat[,'brain'],heart=dpat[,'heart'],qki.hepg2=qki.targ$hepg2)[o,o,c('u','d')]
write.xlsx(rbind(table(brain=dpat[,'brain'],heart=dpat[,'heart'],qki.hepg2=qki.targ$hepg2)[o,o,'u'],
			table(brain=dpat[,'brain'],heart=dpat[,'heart'],qki.hepg2=qki.targ$hepg2)[o,o,'d']),'figures/paper.figures/6/nature.review/qki.encode/hepg2.targets.devAS.xlsx')


# do QKI target have same direction?
f = dpat[,'brain']!='n' & dpat[,'heart']!='n' & dpat[,'brain']!='-' & dpat[,'heart']!='-'
fisher.test(table(dpat[f,'brain'] == dpat[f,'heart'],qki.targ$actaac[f]!=''))

f = dpat[,'brain']!='n' & dpat[,'heart']!='n' & dpat[,'brain']!='-' & dpat[,'heart']!='-' & qki.targ$hepg2 != '-'
fisher.test(table(dpat[f,'brain'] == dpat[f,'heart'],qki.targ$hepg2[f]!='n'))
table(dpat[f,'brain'] , dpat[f,'heart'],qki.targ$hepg2[f]!='n')

o = c('-','n','u','d','ud','du')
table(dpat[,'brain'],dpat[,'heart'],qki=apply(qki.targ!='',1,sum)>0)[o,o,]

table(dpat[,'brain'],qki.targ$actaac)[c('u','d'),c('u','d')]

# _mouse #####
mfa = readRDS('Rdata/ad.alt.fa.Rdata')$mouse
mdpat = age.segs$mouse[anns$mouse$sites=='ad' & rownames(age.segs$mouse) %in% names(mfa),]
dim(mdpat)
mfa = mfa[rownames(mdpat)]
mqki.targ = data.frame(actaac = paste0(ifelse(grepl('actaac',substr(mfa,1,200)),'u',''),
															ifelse(grepl('actaac',substr(mfa,nchar(mfa)-199,nchar(mfa))),'d','')))
mqki.targ$actaa = paste0(ifelse(grepl('actaa',substr(mfa,1,200)),'u',''),
																	ifelse(grepl('actaa',substr(mfa,nchar(mfa)-199,nchar(mfa))),'d',''))
rownames(mqki.targ) = names(mfa)
fisher.test(table(mdpat[,'heart'],mqki.targ$actaac)[c('u','d'),c('u','d')])


plotTissueAgeProile(apply(psi.tsm$mouse[rownames(mqki.targ)[mqki.targ$actaac=='u'],],2,mean,na.rm=T),meta.tsm,age.axis = 'rank')
plotTissueAgeProile(apply(psi.tsm$mouse[rownames(mqki.targ)[mqki.targ$actaac=='d'],],2,mean,na.rm=T),meta.tsm,age.axis = 'rank')


psi.t = calcMeanCols(psi.tsm$mouse,meta.tsm[colnames(psi.tsm$mouse),'tissue'])
# 
b = rep('-',nrow(psi.t))
b[psi.t[,'brain']-apply(psi.t[,-1:-2],1,max) >  0.1] = 'u'
b[psi.t[,'brain']-apply(psi.t[,-1:-2],1,min) < -0.1] = 'd'
names(b) = rownames(psi.t)

h = rep('-',nrow(psi.t))
h[psi.t[,'heart']-apply(psi.t[,-3],1,max) >  0.1] = 'u'
h[psi.t[,'heart']-apply(psi.t[,-3],1,min) < -0.1] = 'd'
names(h) = rownames(psi.t)

# or exlude heart for brain and brain for heart
colnames(psi.t)
b1 = rep('-',nrow(psi.t))
b1[psi.t[,'brain']-apply(psi.t[,-1:-3],1,max) >  0.1] = 'u'
b1[psi.t[,'brain']-apply(psi.t[,-1:-3],1,min) < -0.1] = 'd'
names(b1) = rownames(psi.t)

h1 = rep('-',nrow(psi.t))
h1[psi.t[,'heart']-apply(psi.t[,-1:-3],1,max) >  0.1] = 'u'
h1[psi.t[,'heart']-apply(psi.t[,-1:-3],1,min) < -0.1] = 'd'
names(h1) = rownames(psi.t)


fisher.test(table(mdpat[,'brain'],mqki.targ$actaac)[c('u','d'),c('u','d')])
fisher.test(table(b[rownames(mqki.targ)],mqki.targ$actaac)[c('u','d'),c('u','d')])
fisher.test(table(b1[rownames(mqki.targ)],mqki.targ$actaac)[c('u','d'),c('u','d')])


fisher.test(table(mdpat[,'heart'],mqki.targ$actaa)[c('u','d'),c('u','d')])
fisher.test(table(h[rownames(mqki.targ)],mqki.targ$actaa)[c('u','d'),c('u','d')])
fisher.test(table(h1[rownames(mqki.targ)],mqki.targ$actaa)[c('u','d'),c('u','d')])


table(h[rownames(mdpat)],mdpat[,'heart'])[c(1,3,2),o]
table(b[rownames(mdpat)],mdpat[,'brain'])[c(1,3,2),o]

table(h1[rownames(mdpat)],mdpat[,'heart'])[c(1,3,2),o]
table(b1[rownames(mdpat)],mdpat[,'brain'])[c(1,3,2),o]


f = anns$mouse$sites=='ad'
table(h[f],b[f])
table(h1[f],b1[f])

hb1 = h1
hb1[h1 != b1] = '-'


fisher.test(table(hb1[rownames(mqki.targ)],mqki.targ$actaac)[c('u','d'),c('u','d')])
plotMirroredMotFreq(list(mouse=mfa),list(mouse=cbind(brain=hb1,heart=hb1)),'actaac','heart',main='Heart, mean level')

table(mdpat[,'brain'],mdpat[,'heart'])[o,o]
hbp = mdpat[,'brain']
hbp[mdpat[,'brain']!=mdpat[,'heart']] = '-'
fisher.test(table(hbp[rownames(mqki.targ)],mqki.targ$actaa)[c('u','d'),c('u','d')])
plotMirroredMotFreq(list(mouse=mfa),list(mouse=cbind(brain=hbp,heart=hbp)),'actaa','heart',main='Heart, mean level')


par(mfrow=c(3,2),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(3,2.5,1.5,0),oma=c(0,0,0,1))
plotMirroredMotFreq(list(mouse=mfa),list(mouse=cbind(brain=b,heart=h)),'actaac','brain',main='Brain, mean level')
plotMirroredMotFreq(list(mouse=mfa),list(mouse=cbind(brain=b,heart=h)),'actaac','heart',main='Heart, mean level')

plotMirroredMotFreq(list(mouse=mfa),list(mouse=cbind(brain=b1,heart=h1)),'actaac','brain',main='Brain, mean level')
plotMirroredMotFreq(list(mouse=mfa),list(mouse=cbind(brain=b1,heart=h1)),'actaac','heart',main='Heart, mean level')

plotMirroredMotFreq(list(mouse=mfa),age.segs['mouse'],'actaac','brain',main='Brain, dev patterns')
plotMirroredMotFreq(list(mouse=mfa),age.segs['mouse'],'actaac','heart',main='Heart, dev patterns')


# qki5 exp by last exon #########
# hum.57513.s19, mou.19182.s9 (mouse exon is part of last exon annotated in ensembl, but I think it is fine
h = readRDS('Rdata/human.as.u.all.Rdata')
h$seg[h$seg$gene_id=='hum.57513',]
hrc = read.table('processed/mapping/hisat2.s/human/mapping.stat')
colnames(hrc) = c('lib.id','read.count','mapped.uniq','mapped.mult')
table(meta$species,meta$lib.id %in% hrc$lib.id)
hrc = hrc[hrc$lib.id %in% meta$lib.id,]
rownames(hrc) = rownames(meta)[match(hrc$lib.id,meta$lib.id)]
hrc = hrc[colnames(h$i),]
qki5.human.cpm = h$i['hum.57513.s19',]/(hrc$mapped.uniq)*1e6
rm(h);gc()

# m = readRDS('Rdata/mouse.as.u.all.Rdata')
# m$seg[m$seg$gene_id=='mou.19182',] 
# m$seg[m$seg$gene_id=='mou.19182' & m$seg$stop >=10202601 & m$seg$start <= 10209631, ] # 10,209,631-10,202,601	 (according to ensembl)
# plotTissueAgeProile(m$ir['mou.19182.s8',],meta,age.axis = 'rank',bty='n',ylab='CPM',ylim=c(0,1))
# mrc = read.table('processed/mapping/hisat2.s/mouse/mapping.stat')
# colnames(mrc) = c('lib.id','read.count','mapped.uniq','mapped.mult')
# table(meta$species,meta$lib.id %in% mrc$lib.id)
# mrc = mrc[mrc$lib.id %in% meta$lib.id,]
# rownames(mrc) = rownames(meta)[match(mrc$lib.id,meta$lib.id)]
# mrc = mrc[colnames(m$i),]
# qki5.mouse.cpm = m$i['mou.19182.s9',]/(mrc$mapped.uniq)*1e6
# rm(m);gc()
# saveRDS(qki5.mouse.cpm,'Rdata/qki5.mouse.cpm.Rdata')


pdf('figures/paper.figures/6/nature.review/reviewer.fig/R1.v2.qki-5.expression.by.last.exon.pdf',w=3,h=3)
par(mfrow=c(1,1),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
m = meta[meta$tissue %in% c('brain','heart') & meta$days <= 14515,]
#plotTissueAgeProile(qki5.human.cpm,m,age.axis = 'rank',bty='n',ylab='CPM',main='Human QKI-5',df = 4,xlab='Age',pch=19,cex=1)
plotTissueAgeProile(qki5.mouse.cpm,m,age.axis = 'rank',bty='n',ylab='CPM',main='Mouse QKI-5',df = 4,xlab='Age',pch=19,cex=1,ylim=range(0,qki5.mouse.cpm,na.rm=T))
dev.off()
