library(seqinr)
library(MyLogo)


getCisExplanedPSIonThr = function(p2s,thrs=0:19/20,sites=c('ad','da','dd','aa'),substractUnexplained=TRUE){
	p2s = p2s[p2s[,3] %in%  sites,]
	x=sapply(thrs,function(f){
		t = abs(p2s[,1]) > f
		(sum(t & (sign(p2s[,1]) == sign(p2s[,2])))- #number of explained
			sum(t & p2s[,2] != 0 & (sign(p2s[,1]) != sign(p2s[,2])))*substractUnexplained)/sum(t) # minus number of contr-explained. dss equal to zero is not not-explained and not contr-explained.
	})
	list(x=thrs,y=x)
}


getAllSpeciesSeq = function(sid,site.no,snp.ref.freq,rm.na.ref.freq){
	seg = ha1$seg[sid,]
	site.type = substr(seg$sites,site.no,site.no)
	
	if(site.type == 'a')
		inx = 3:27
	
	
	
	
	if(site.type == 'd'){
		r=getSNPInterSpeciesDiff(sid,ha1$seg,ra1$seg,site.no,don.ort.hs,don.ort.rm,sites,snp,sites2snp,24:32,
														 snp.ref.freq,rm.na.ref.freq,only.equal.in.species.ref=FALSE)[1,]
		h = toupper(getSegSiteSeq(ha1$seg[sid,],site.no,don.ort.hs))
		p = toupper(getSegSiteSeq(pa1$seg[sid,],site.no,don.ort.pt))
		t = toupper(getSegSiteSeq(ra1$seg[sid,],site.no,don.ort.rm))
	}else{
		r=getSNPInterSpeciesDiff(sid,ha1$seg,ra1$seg,site.no,acc.ort.hs,acc.ort.rm,sites,snp,sites2snp,3:27 ,
														 snp.ref.freq,rm.na.ref.freq,only.equal.in.species.ref=FALSE)[1,]
		h = toupper(getSegSiteSeq(ha1$seg[sid,],site.no,acc.ort.hs))
		p = toupper(getSegSiteSeq(pa1$seg[sid,],site.no,acc.ort.pt))
		t = toupper(getSegSiteSeq(ra1$seg[sid,],site.no,acc.ort.rm))
	}
	
	r = sapply(r,function(x){if(!is.na(x)){paste(substr(x,2,2),'->',substr(x,1,1),sep='')}else{x}})
	names(r) = NULL
	r=rbind(r,h,p,t)
	rownames(r) = c('hs.snp','hs','pt','rm')
	r
}


plot3spPSIhist = function(sid,breaks=0:10/10,main=sid,...){
	hst = rbind(hs=hist(ha1$ir[sid,],breaks,plot=F)$counts,
							pt=hist(pa1$ir[sid,],breaks,plot=F)$counts,
							rm=hist(ra1$ir[sid,],breaks,plot=F)$counts)
	barplot(hst,beside=T,col=c('red','blue','green'),xlab='PSI',ylab='# of sample',main=main)
	axis(1,0:5/5*40,0:5/5)
	legend('topright',fill=c('red','blue','green'),legend=c('human','chimpanzee','monkey'))
}

getSNPForSeg = function(sid,max.ref.freq=1,rm.na.ref.freq=FALSE){
	z = ha1$seg[sid,]
	ss = s2c(z$sites)
	z = c(paste(z$chr_id,z$strand,z$start,ifelse(z$strand==1,ss[1],ss[2]),sep='.'),paste(z$chr_id,z$strand,z$stop+1,ifelse(z$strand==1,ss[2],ss[1]),sep='.'))
	z = (z)
	r=snp[unique(sites2snp[sites2snp[,1] %in% z,3]),]
	r = r[is.na(r$ref.freq) | r$ref.freq<=max.ref.freq,]
	if(rm.na.ref.freq)
		r = r[!is.na(r$ref.freq),]
	r
}

plotChangePropExplainedByMutations = function(p2s,thrs,sites,main){
	r = t(sapply(thrs,function(x){c(thr=x,calcChangePropExplainedByMutations(p2s,x,sites))}))
	x = log(r[,1]+0.1)
	plot(x,r[,4],t='b',col='red',ylim=c(0,max(r[,3])),xaxt='n',xlab='splicing propensity change thr (bits)',ylab='# of changes',main=main)
	lines(x,r[,5],t='b',col='blue')
	lines(x,r[,4]-r[,5],t='b',col='green')
	lines(x,r[,3],t='b',col='black')
	axis(1,at=x,r[,1])
	mtext('proportion of significant changes',4,1.1)
	l = 0:20/20
	axis(4,l*r[1,2],l)
	abline(h=l*r[1,2],lty=2,col='#999999')
	legend('topright',lwd=2,col=c('black','red','green','blue'),legend=c('dPWM > thr','correct','explained','wrong'))
}


calcChangePropExplainedByMutations = function(x,thr=0,sites=c('ad','aa','dd','da')){
	x = x[x$sites %in% sites,]
	c(total=nrow(x),
		not0=sum(abs(x[,2]) > thr),
		expl=sum((sign(x[,1]) == sign(x[,2])) & (abs(x[,2]) > thr)),
		wrong=sum((sign(x[,1]) != sign(x[,2])) & (abs(x[,2]) > thr)))
}

getSNPInterSpeciesDiff = function(seg.inx,seg1,seg2,site.no,seq1,seq2,sts1,snps1,s2s1,snp.pos.inx,max.ref.freq=1,rm.na.freq=FALSE,only.equal.in.species.ref=TRUE){
	# findes positions wherer references are equal in both species, first species have SNPs
	snps1 = snps1[(!rm.na.freq & is.na(snps1$ref.freq)) | (!is.na(snps1$ref.freq) & snps1$ref.freq <= max.ref.freq),]
	seq1 = toupper(getSegSiteSeq(seg1[seg.inx,],site.no,seq1))
	seq2 = toupper(getSegSiteSeq(seg2[seg.inx,],site.no,seq2))
	snp.seq = getSNPforSeg(seg1[seg.inx,],sts1,snps1,s2s1,site.no)
	snp.seq$ref = snp.seq$ref[,snp.pos.inx]
	snp.seq$alt = snp.seq$alt[,snp.pos.inx]
	r = seq1
	r[] = NA
	if(sum(snp.seq$ref != seq1,na.rm=TRUE)>0)
		stop("ERROR: Some SNP ref aren't equal to expected ref!")
	print(table(ref.diff = seq1 != seq2,snp=!is.na(snp.seq$alt)))
	if(only.equal.in.species.ref)
		f = (seq1 == seq2) & !is.na(snp.seq$alt)
	else
		f = !is.na(snp.seq$alt)
	r[f] = paste(snp.seq$alt[f],snp.seq$ref[f],sep='')
	r
}

getSNPforSeg = function(segs,sts,snps,s2s,site.no){
	s = data.frame(chr_id=segs$chr_id,strand=segs$strand,pos=ifelse(!xor(segs$strand == 1,site.no==1),segs$start,segs$stop+1),site = substr(segs$sites,site.no,site.no))
	getSNPforSites(s,s2s,snps)
}


getSNPforSites = function(sites,s2s,snps){
	rownames(sites) = do.call(paste,c(sites,sep='.'))
	sites$inx = 1:nrow(sites)
	
	snps = snps[!is.na(snps$max.nonref.allele),]
	s2s = s2s[s2s[,1] %in% rownames(sites) & s2s[,3] %in% rownames(snps),]
	snp = snps[s2s[,3],]
	snp$dist = snp$V4 - sites[s2s[,1],'pos'] 
	snp$dist[sites[s2s[,1],'strand']==-1] = -snp$dist[sites[s2s[,1],'strand']==-1]-1
	snp$dist = snp$dist + 27
	
	ref = alt = matrix(NA,ncol=52,nrow=nrow(sites))
	if(nrow(s2s)>0)
		for(i in 1:nrow(s2s)){
			#cat('\r',i,snp$V9[i],snp$max.nonref.allele[i])
			site.inx = sites[s2s[i,1],'inx']
			if(sites[s2s[i,1],'strand'] == 1){
				ref[site.inx,snp$dist[i]] = snp$V9[i]
				alt[site.inx,snp$dist[i]] = snp$max.nonref.allele[i]
			}else{
				ref[site.inx,snp$dist[i]] = comp(snp$V9[i],forceToLower = FALSE)
				alt[site.inx,snp$dist[i]] = comp(snp$max.nonref.allele[i],forceToLower = FALSE)
			}
		}
	list(ref=ref,alt=alt)
}



readAndFilterSNPs = function(infile){
	snp = read.table(infile,sep='\t')
	#remove all not primary
	snp = snp[!grepl('_',snp$V2,TRUE),]
	# some "genomic single" have length greater than one. Remove them
	snp[snp$V4-snp$V3!=1,][1:10,]
	snp = snp[snp$V4-snp$V3==1,]
	# have no idea what 'H' means
	snp = snp[snp$V10 != 'G/H',]
	#SNP in some parts of X and Y chrs have the same ids
	rownames(snp) = paste(snp$V2,snp$V5,sep='.')
	
	#some alleles are longer than 1 nt, and I have no idea what '0' allele means
	#snp = snp[sapply(alleles,function(x){(length(x) == 0 | max(nchar(x))==1) & sum(x %in% c('N','0'))==0}),]
	
	observed = strsplit(snp[,10],'/',TRUE)
	#print(table(sapply(1:nrow(snp),function(i){snp$V9[i] %in% observed[[i]]}),snp$V7))
	observed = lapply(1:length(observed),function(i){
		r = observed[[i]]
		if(snp$V7[i] == '-')
			r = comp(r,forceToLower = FALSE)
		r
		})
	#print(table(sapply(1:nrow(snp),function(i){snp$V9[i] %in% observed[[i]]}),snp$V7))
	alleles = strsplit(snp[,23],',',TRUE)
	allele.freq = lapply(strsplit(snp[,25],',',TRUE),as.numeric)
	allele.freq = lapply(1:length(allele.freq),function(i){
		r = setNames(allele.freq[[i]],alleles[[i]])[observed[[i]]]
		names(r) = observed[[i]]
		r[is.na(r)] = 0
		if(sum(r)>0)
			r = r/sum(r)
		else
			r[] = NA
		r
		})
	#snp$max.freq = sapply(allele.freq,max)

	snp$ref.freq = sapply(1:length(allele.freq),function(i){allele.freq[[i]][snp$V9[i]]})
	
	snp$max.nonref.allele = sapply(1:nrow(snp),function(i){
		# take allele with max frequency if frequency of alt allele is defines
		a = allele.freq[[i]]
		a = a[names(a) != snp$V9[i]]
		a = sort(a,decreasing = TRUE,na.last=TRUE)
		r = NA
		if(!is.na(a[1]) || length(a) == 1)
				r = names(a)[1]
		r
		})
	snp
}

getSegSiteSeq = function(seg,site.no,seq){
	sites = paste(seg$chr_id,":",ifelse(seg$strand==1,'+','-'),":",ifelse(!xor(seg$strand == 1,site.no==1),seg$start-1,seg$stop),sep='')
	do.call(rbind,lapply(seq[sites,1],s2c))
}

findInterSpecDiff = function(seg.inx,seg1,seg2,site.no,seq1,seq2){
	seq1 = getSegSiteSeq(seg1[seg.inx,],site.no,seq1)
	seq2 = getSegSiteSeq(seg2[seg.inx,],site.no,seq2)
	r = seq1
	r[seq1 == seq2] = NA
	f = seq1 != seq2
	r[f] = paste(seq1[f],seq2[f],sep='')
	r
}


plotMutations = function(d1,d2,s,dpsi,logo1,logo2,main=''){
	par(tck=-0.02,mgp=c(1.4,0.2,0),mar=c(0.5,3,1.5,0),oma=c(0,0,1,1))
	pv.thrs = c(1,0.1,0.05,0.01,0.005,0.001,0)
	layout(matrix(c(1,2,3,3,4,4,5,6,7,7,8,8),ncol=2))
	plotNOofMutations(d1,NULL,s,'')
	plotLogo(logo1$seq,ic.scale = TRUE,order.by.freq = TRUE,xlabs=NA)
	par(mar=c(3.5,3,1.5,0))
	plotMutationMatrix(calcMutationEnrichment(d1,s & dpsi>0,!s),pv.thrs=pv.thrs,main='rm down',ylab='hs->rm')
	plotMutationMatrix(calcMutationEnrichment(d1,s & dpsi<0,!s),pv.thrs=pv.thrs,main='rm up',ylab='hs->rm')
	
	par(mar=c(0.5,3,1.5,0))
	plotNOofMutations(d2,NULL,s,'')
	plotLogo(logo2$seq,ic.scale = TRUE,order.by.freq = TRUE,xlabs=NA)
	par(mar=c(3.5,3,1.5,0))
	plotMutationMatrix(calcMutationEnrichment(d2,s & dpsi>0,!s),pv.thrs=pv.thrs,main='rm down',ylab='hs->rm')
	plotMutationMatrix(calcMutationEnrichment(d2,s & dpsi<0,!s),pv.thrs=pv.thrs,main='rm up',ylab='hs->rm')
	mtext(main,3,-0.5,outer=TRUE)
}

plotMutationMatrix = function(m,main='',pv.thrs = c(1,0.1,0.05,0.01,0.005,0.001,0),ylab=''){
	z   = t(m$pv)
	z[] = findInterval(-z,-pv.thrs)
	col. = c('white',rev(heat.colors(length(pv.thrs)-2)))
	col = col.[min(z):max(z)]
	image(x=1:nrow(z),y=1:ncol(z),z=z[,ncol(z):1],col=col,xaxt='n',yaxt='n',ylab=ylab,xlab='',main=main,xaxs='r')
	
	axis(2,nrow(m$odds):1,rownames(m$odds),las=2)
	text(rep(1:ncol(m$odds),each=nrow(m$odds)),rep(nrow(m$odds):1,times=ncol(m$odds)),round(m$odds,digits=0))
	#rect(0,0.5,ncol(m$odds)+0.5,nrow(m$odds)+0.5)
	par(xpd=TRUE)
	xx = (0:(length(pv.thrs)-1))/length(pv.thrs)*(nrow(z))+1.5
	yy = c(-1,0)
	rect(xx[-length(xx)],rep(yy[1],length(xx)-1),xx[-1],rep(yy[2],length(xx)-1),col=col.)
	text(0.5,-0.5,'pv>',c(0.5,0.5))
	text((xx[-length(xx)]+xx[-1])/2,rep(-0.5,length(xx)-1),pv.thrs[-length(pv.thrs)],c(0.5,0.5))
	par(xpd=FALSE)
}

calcMutationEnrichment = function(m,s,ns){
	m = toupper(m)
	mv = c('AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG')
	sgn = apply(m[ s,],2,function(x){table(factor(x[x %in% mv],levels = mv))})
	nsg = apply(m[ns,],2,function(x){table(factor(x[x %in% mv],levels = mv))})
	odd = pv = sgn
	sgn.t = sum( s)
	nsg.t = sum(ns)
	for(i in 1:nrow(sgn))
		for(j in 1:ncol(sgn)){
			t = matrix(c(sgn[i,j],sgn.t-sgn[i,j],nsg[i,j],nsg.t-nsg[i,j]),ncol=2)
			t = fisher.test(t,alternative = 'greater')
			odd[i,j] = t$estimate
			pv[i,j] = t$p.value
		}
	list(odds=odd,pv=pv,sgn=sgn,nsg=nsg)
}

compareMutEffects = function(d1,d2,sgn,dpsi,dir,pwm1,pwm2){
	d1 = dir[1]*calcMutationEffect(d1,pwm1)
	d2 = dir[2]*calcMutationEffect(d2,pwm2)
	d = cbind(d1,d2)
	f = function(d){
		table(sign(d))
	}
	rbind(up = f(d[sgn & dpsi > 0 ,]),
				no = f(d[!sgn,]),
				dn = f(d[sgn & dpsi < 0 ,]))
}

plotNOofMutations = function(d1,d2=NULL,s,main=''){
	if(is.null(d2)){
		sgn = c(apply(!is.na(d1[ s,]),2,sum))/sum( s)
		nsg = c(apply(!is.na(d1[!s,]),2,sum))/sum(!s)
	}else{
		sgn = c(apply(!is.na(d1[ s,]),2,sum),0,0,apply(!is.na(d2[ s,]),2,sum))/sum( s)
		nsg = c(apply(!is.na(d1[!s,]),2,sum),0,0,apply(!is.na(d2[!s,]),2,sum))/sum(!s)
	}
	barplot(sgn*100,col='red',main=main,ylab='% of segs')
	barplot(nsg*100,col='blue',density=10,add=T)
}

calcMutationEffect = function(mutations,pwm){
	r = matrix(0,ncol=ncol(mutations),nrow=nrow(mutations))
	for(i in 1:nrow(r)){
		for(j in 1:ncol(r)){
			if(!is.na(mutations[i,j])){
				m = s2c(toupper(mutations[i,j]))
				if(sum(m %in% rownames(pwm)) == length(m))
					r[i,j] = pwm[m[2],j] - pwm[m[1],j]
			}
		}
	}
	r
}

findInterSpecSiteSeqDiff = function(seg.inx,seg1,seg2,site.no,seq1,seq2){
	seg1 = seg1[seg.inx,]
	seg2 = seg2[seg.inx,]
	sites1 = paste(seg1$chr_id,":",ifelse(seg1$strand==1,'+','-'),":",ifelse(!xor(seg1$strand == 1,site.no==1),seg1$start-1,seg1$stop),sep='')
	sites2 = paste(seg2$chr_id,":",ifelse(seg2$strand==1,'+','-'),":",ifelse(!xor(seg2$strand == 1,site.no==1),seg2$start-1,seg2$stop),sep='')
	seq1 = do.call(rbind,lapply(seq1[sites1,1],s2c))
	seq2 = do.call(rbind,lapply(seq2[sites2,1],s2c))
	r = seq1
	r[seq1 == seq2] = NA
	f = seq1 != seq2
	r[f] = paste(seq1[f],seq2[f],sep='')
	r
}

loadSites = function(f){
	f = read.fasta(f,as.string = TRUE)
	f = data.frame(pos=names(f),seq = as.character(f))
	f = unique(f)
	rownames(f) = f$pos
	f[,'seq',drop=FALSE]
}

getSeqWeight = function(pwm,seq,from=1){
	if(nchar(seq)-from+1 < ncol(pwm) | from < 1)
		return(NA)
	seq = s2c(toupper(substr(seq,from,from+dim(pwm)[2]-1)))
	if(sum(seq %in% rownames(pwm)) != length(seq))
		return(NA);
	s = 0
	for(i in 1:length(seq))
		s = s + pwm[seq[i],i]
	s = s
	names(s) = NULL
	s
}


pwm = function(seq,prior=c(A=0.25, C=0.25, G=0.25, T=0.25)){
	prior = prior[order(names(prior))]
	seq = toupper(seq)
	cnt = as.matrix(as.data.frame(consensusMatrix(seq))[names(prior),])
	cnt[is.na(cnt)] = 0
	rownames(cnt) = names(prior)
	post = sweep(cnt,1,prior,'+')/(length(seq)+sum(prior))
	prior = prior/sum(prior)
	log2(sweep(post,1,prior,'/'))
}

calcSegSplPower = function(seg,a,d){
	rownames(a) = paste('a:',rownames(a),sep='')
	rownames(d) = paste('d:',rownames(d),sep='')
	s = rbind(a,d)
	
	sites = cbind(substr(seg$sites,1,1),substr(seg$sites,2,2))
	sites. = sites
	sites[,1] = paste(sites[,1],":",seg$chr_id,":",ifelse(seg$strand==1,'+','-'),":",ifelse(seg$strand== 1,seg$start-1,seg$stop),sep='')
	sites[,2] = paste(sites[,2],":",seg$chr_id,":",ifelse(seg$strand==1,'+','-'),":",ifelse(seg$strand==-1,seg$start-1,seg$stop),sep='')
	
	w5 = s[sites[,1],'weight']
	w3 = s[sites[,2],'weight']
	ifelse(sites.[,1]=='a',1,-1)*w5 + ifelse(sites.[,2]=='d',1,-1)*w3
}

plot2Dens = function(pwm,ir,main,min.pwm,show.reg=2,xlab='Site weight change',ylab='Splicing change'){
	tot = length(ir)
	to.use = !is.na(pwm) & abs(pwm) > min.pwm
	ir = ir[to.use]
	pwm = pwm[to.use]
	xlim=sd(pwm)*show.reg
	ylim=sd(ir)*show.reg
	d = kde2d(pwm,ir,n=100,lims=c(-xlim,xlim,-ylim,ylim))
	image(d,xlab=xlab,ylab=ylab,main=paste(main,' (',sum(sign(pwm)==sign(ir)),'/',length(pwm),'/',tot,')',sep=''),col=rev(heat.colors(100)))
	contour(d,add=T,xaxt='n',col='darkgray')
	points(pwm,ir,col='#FF0000',pch=16,cex=0.2)		
	lines(lowess(pwm,ir),lwd=2,col='blue')
	abline(v=-min.pwm,col='green',lty=2)
	abline(v= min.pwm,col='green',lty=2)
	abline(h=0,col='black',lty=2)
	abline(v=0,col='black',lty=2)
	text(-xlim,ylim,
			 labels=paste('rho=',round(cor(pwm,ir),3),sep=''),adj=c(0,1))
	box()
}