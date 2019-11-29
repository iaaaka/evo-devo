library(vioplot)
library(reshape)

plotProportionOfCannonicalSites = function(x,species,min.sam.no,more.or.equal.fun,...){
	get.freq = function(l,f){
		r = split(l,f)
		r = r[order(as.numeric(names(r)))]
		t = rev(cumsum(rev(sapply(r,length))))
		c = rev(cumsum(rev(sapply(r,sum))))
		list(y=c/t,x=as.numeric(names(t)))
	}
	canon.seq = x$introns$canonical
	self = x$sam.cnts[,species] >= min.sam.no
	sp.cnt = apply(x$sam.cnts >= min.sam.no,1,sum)
	sam.cnt = apply(x$sam.cnts,1,sum)
	canon.freq = list()
	sc = ncol(x$sam.cnts)
	for(i in 1:sc){
		ff = self & .Primitive(more.or.equal.fun)(sp.cnt,i)
		canon.freq[[i]] = get.freq(canon.seq[ff],sam.cnt[ff])
		ff = !self & .Primitive(more.or.equal.fun)(sp.cnt,i)
		if(sum(ff)>0)
			canon.freq[[sc+i]] = get.freq(canon.seq[ff],sam.cnt[ff])
	}
	ylim=range(sapply(canon.freq,function(z){range(z$y)}))
	xlim=range(sapply(canon.freq,function(z){range(z$x)}))
	plot(1,t='n',xlim=xlim,ylim=c(0.1,1),log='',xlab='# samples',ylab='% of canonical sites',...)
	col=getPal(c('#0000FF','#FF0000'),sc)
	col=c(col,paste(substr(col[-1],1,7),'60',sep=''))
	for(i in 1:length(canon.freq))
		lines(canon.freq[[i]],lwd=3,col=col[i],lty=(i>sc)+1)
	abline(h=0.95)
	legend('bottomright',col=col,lty=(1:length(col)>sc)+1,legend=c(paste(' self, sp',more.or.equal.fun,1:7,sep=''),paste('!self, sp',more.or.equal.fun,1:6,sep='')),ncol=2)
}

plotIntronSampleCntBySpecies = function(x,...){
	lens = split(log10(apply(x$sam.cnts,1,sum)+1),x$introns$species)
	nms = gsub('-','',names(lens),fixed = TRUE)
	
	lens = lens[order(-nchar(nms),-sapply(lens,sum))]
	ymax=max(unlist(lens))
	plot(1,t='n',xaxt='n',yaxt='n',ylim=c(0,ymax),xlim=c(0.5,length(lens)+0.5),...)
	axis(1,at=1:length(lens),names(lens),las=2)
	lab=10^(0:ymax)
	axis(2,at=log10(lab+1),labels = lab)
	col = nchar(gsub('-','',names(lens),fixed = TRUE))+1
	for(i in 1:length(lens))
		vioplot(lens[[i]],at = i,col=col[i],add=TRUE)
	#plot intron count
	cnt = log10(sapply(lens,length))
	lab = 10^(0:max(cnt))
	at = log10(lab)/max(cnt)*ymax
	cnt = cnt/max(cnt)*ymax
	points(1:length(cnt),cnt,pch=19,cex=2,type='b')
	axis(4,at=at,labels = lab)
	mtext('# of introns',side = 4,line = 2)
}

saveGoodSites = function(d,thr.canon,thr.other,file){
	cnt = apply(d$sam.cnts,1,max)
	f = (cnt >= thr.canon & d$intron$seqs %in% c('GTAG','CTAC')) | cnt >= thr.other
	cat('Saved',sum(f),'introns\n')
	write.table(d$intron[f,c('chr.id','start','stop','strand')],file,sep='\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
}

loadMergedLiftovered = function(f,species){
	species.r = setNames(names(species),species)
	sp.order = c('mouse','chicken','rabbit','rat','human','macaque','opossum') #it is order used in merge.lo.py
	species.r = species.r[sp.order]
	r = read.table(f,sep='\t')
	i = r[,1:6]
	colnames(i) =  c('chr.id','start','stop','strand','seq','doubled.sp')
	i$doubled.sp[!is.na(i$doubled.sp)] = sapply(strsplit(as.character(i$doubled.sp[!is.na(i$doubled.sp)]),''),function(x){paste(species.r[as.numeric(x)],collapse='')})
	iden = r[,7:13]
	s = r[,14:20]
	colnames(s) = colnames(iden) = sp.order
	s = s[,species]
	iden = iden[,species]
	i$species = apply(s,1,function(x){paste(ifelse(x>0,species.r,'-'),collapse='')})
	r = list(introns = i,sam.cnts=s,identity=iden)
	class(r) = c('sajr','list')
	r
}

loadLiftovered = function(original,species,path='~/iitp.disk/Solexa/ma.disk/mazin/evo.devo/mapping/junctions/',sp.names){
	stop("different introns can be liftovered into different positions!! Use merge.lo.py instead")
	original = original[,c('seq','sam.cnt','strand')]
	files = list.files(path,paste(species,'.out',sep=''))
	lo = vector('list',length(files))
	names(lo) = sapply(strsplit(files,'2',TRUE),'[',1)
	for(i in 1:length(files)){
		t = read.table(paste(path,files[i],sep='/'),skipNul = TRUE)[,c(1,4,5,7,9)] # skipNul = TRUE to evide 'embedded nul(s) found in input'... have no idea what it means
		lo[[i]] = setNames(t[,5],do.call(paste,c(t[,1:4],sep=':')))
	}
	
	all.ints = unique(c(rownames(original),unlist(lapply(lo,names))))
	#seqs = sam.cnts  = strands = matrix(nrow=length(all.ints),ncol=length(lo)+1)
	#rownames(seqs) = rownames(sam.cnts) =  rownames(strands) = all.ints
	#colnames(seqs) = colnames(sam.cnts) =  colnames(strands) =  c(species,sapply(strsplit(names(lo),'To',TRUE),'[',1))
	sam.cnts = matrix(nrow=length(all.ints),ncol=length(lo)+1)
	rownames(sam.cnts) = all.ints
	colnames(sam.cnts) = c(species,sapply(strsplit(names(lo),'To',TRUE),'[',1))
	
	#strands[rownames(original),1] = original[,'strand']
	#seqs[rownames(original),1] = original[,'seq']
	sam.cnts[rownames(original),1] = original[,'sam.cnt']
	lo = lapply(lo,function(x){strsplit(x,':',TRUE)})
	for(i in 1:length(lo)){
		print(i)
		#seqs[names(lo[[i]]),1+i] = sapply(lo[[i]],'[',6)
		sam.cnts[names(lo[[i]]),1+i] = as.numeric(sapply(lo[[i]],'[',7))
		#strands[names(lo[[i]]),1+i] = sapply(lo[[i]],'[',5)
	}
	introns = as.data.frame(do.call(rbind,strsplit(all.ints,':',TRUE)))
	colnames(introns) = c('chr.id','start','stop','strand')
	introns$start = as.integer(introns$start)
	introns$stop = as.integer(introns$stop)
	sam.cnts[is.na(sam.cnts)] = 0
	sp = sp.names[colnames(sam.cnts)]
	introns$species = apply(sam.cnts>0,1,function(x){paste(sort(sp[x]),collapse='')})
	r = list(introns = introns,sam.cnts=sam.cnts)#,seqs=seqs,strands=strands)
	class(r) = c('sajr','list')
	r
}

printIntronsAsGFF = function(i,min.sam.cnt,species,out){
	i = i[i$sam.cnt >= min.sam.cnt,]
	i = cbind(i$chr.id,'.','.',i$start,i$stop,'.',i$strand,'.',paste(species,rownames(i),i$seq,i$sam.cnt,sep=':'))
	write.table(i,file = out,col.names = FALSE,row.names = FALSE,quote = FALSE,sep='\t')
}


plotintronStatForSpecies = function(x,species){
	o = order(apply(sweep(x$seq.stat,1,apply(x$seq.stat,1,sum),'/')[,c('GTAG','CTAC')],1,sum))
	plotSeqStat(x$seq.stat[o,],xlab='Samples',main=paste(species,'. By samples.',sep=''),ylab='freq',ylim = c(0,1.5),yrange = c(0.5,1.5))
	plotSeqStat(x$seq.stat[o,],ord=c('ATAC','GCAG','CTGC','GTAT','others'),add=T,yrange = c(0,0.45))
	abline(h=0.475)
	
	
	t = table(x$sam.stat$sam.cnt,x$sam.stat$seq)
	plotSeqStat(t,xlab='# samples',main=paste(species,'. By number of samples where detected.',sep=''),ylab='freq',xlog=TRUE,yrange = c(0.5,1.5),ylim=c(0,1.5),plot.total = T)
	plotSeqStat(t,ord=c('ATAC','GCAG','CTGC','GTAT'),xlog=TRUE,add=T,yrange = c(0,0.45))
	mtext('# of introns',4,2.3,FALSE)
	abline(h=0.475)
}

loadJunctionStat = function(dir2sam,merged){
	r = read.table(merged,stringsAsFactors = FALSE)[,c(1,4,5,7,9)]
	colnames(r) = c('chr.id','start','stop','strand','seq')
	r$seq = sapply(strsplit(r$seq,':',fixed = TRUE),'[',6)
	rownames(r) = do.call(paste,c(r[,1:4],sep=':'))
	r$sam.cnt = 0
	sams = list.files(dir2sam,'.splicesites')
	unuq.seqs = unique(r$seq)
	seqs = matrix(nrow = length(sams),ncol=length(unuq.seqs))
	colnames(seqs) = unuq.seqs
	rownames(seqs) = sams
	for(i in 1:length(sams)){
		cat('\r',i,length(sams))
		t = read.table(paste(dir2sam,sams[i],sep='/'),stringsAsFactors = FALSE)
		t = do.call(paste,c(t[,1:4],sep=':'))
		seqs[i,] = table(factor(r[t,'seq'],levels = unuq.seqs))[unuq.seqs]
		r[t,'sam.cnt'] = r[t,'sam.cnt'] + 1
	}
	list(sam.stat=r,seq.stat = seqs)
}

plotSeqStat = function(t,...,add=F,yrange=c(0,1),ylim=NULL,xlog=FALSE,ord=c('ATAC','GCAG','GTAG','others','CTAC','CTGC','GTAT'),col=c('ATAC'='orange','GCAG'='yellow','GTAG'='red','others'='gray','CTAC'='blue','CTGC'='green','GTAT'='cyan'),plot.leg=!add,plot.total=FALSE){
	sam.cnts = suppressWarnings(as.numeric(rownames(t)))
	if(sum(is.na(sam.cnts)) > 0)
		sam.cnts = 1:nrow(t)
	scale = function(x,range,mn=min(x),mx=max(x)){
		(x - mn)/(mx-mn)*(range[2]-range[1]) + range[1]
	}
	xfun = xfun. = function(x){x}
	xat = c(1,seq(50,by = 50,to=max(51,nrow(t))))
	if(xlog){
		xfun = function(x,to=max(x)+1){
			log((to-1)*x/(to-x))
		}
		
		xfun. = function(x,to=max(x)+1){
			x = exp(y)
			to*x/(to-1+x)
		}
		xat = c(1:3,5,10,20,40)
		xat = c(xat,seq(100,by=100,to = max(101,max(sam.cnts)-max(xat))),max(sam.cnts)-xat+1)
	}

	seq.freq = t[,c('GTAG','CTAC','GCAG','CTGC','ATAC','GTAT')]
	seq.freq = cbind(seq.freq,others = apply(t[,!(colnames(t) %in% colnames(seq.freq))],1,sum))
	seq.freq = sweep(seq.freq,1,apply(seq.freq,1,sum),'/')

	seq.freq = seq.freq[,ord]
	col = rev(col[ord])
	seq.freq = t(apply(seq.freq,1,function(x){rev(cumsum(x))}))
	freq.range = range(0,seq.freq,na.rm=T)
	seq.freq = scale(seq.freq,yrange)
	if(is.null(ylim))
		ylim = c(0,max(apply(seq.freq,1,max)))
	seq.freq = cbind(seq.freq,yrange[1])
	
	x = xfun(c(sam.cnts,rev(sam.cnts)))
	if(!add){
		plot(1,t='n',xlim=c(0,max(x)*1.2),ylim=ylim,yaxt='n',xaxt='n',...)
		axis(1,at=xfun(xat),xat,las=3)
	}
	
	for(i in 2:ncol(seq.freq)){
		polygon(x,c(seq.freq[,i-1],rev(seq.freq[,i])),col=col[i-1],border = NA)
	}
	if(plot.leg)
		legend(xfun(max(sam.cnts))*1.01,yrange[2],fill=col,legend=names(col))
	if(plot.total){
		total = log(apply(t,1,sum))
		
		lab = c(1,3,10,30,100,300,1000,3000,10000,30000,100000,300000,1000000,3000000)
		lab = lab[log(lab)>min(total) & log(lab)<max(total)]
		at = log(lab)
		at = scale(at,yrange,min(total),max(total))
		
		total = scale(total,yrange)
		lines(xfun(sam.cnts),total,lwd=3)

		axis(4,at,lab,las=2,tck=0.02,mgp=c(1.1,0.1,0))
	}
	lab = seq(from=freq.range[1],to=freq.range[2],length.out = 5)
	n0 = as.character(lab)
	n0 = max(nchar(n0) - nchar(gsub('0\\.0*','',n0,perl=TRUE)))-1
	lab = unique(round(lab,n0))
	at = scale(lab,yrange)
	axis(2,at=at,labels = lab,las=2,tck=0.02,mgp=c(1.1,0.1,0))
}

