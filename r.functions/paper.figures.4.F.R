# common #####

my.binom.test = function(s,f=NULL){
	if(length(s)>1){
		f = s[2]
		s = s[1]
	}
	if(f+s == 0)
		return(c(NA,NA,NA))
	r=binom.test(s,s+f)
	c(r$estimate,r$conf.int)
}

getAgeASEns = function(psi,m,dPSI,bs,sp,sites='ad',cod=c('c','p')){
	m = getAgeASchanges(psi,m,dPSI,bs,sp)
	m.ens = apply(m,2,function(x){
		r=lapply(c('n','u','d'),function(d){
			unique(unlist(seg2ens[[sp]][rownames(m)[x==d & anns[[sp]]$sites %in% sites & anns[[sp]]$cod %in% cod]]))
		})
		names(r) = c('n','u','d')
		r$n = setdiff(r$n,c(r$u,r$d))
		r
	})
	m.ens
}

# figure 1 ######
getNoOfEvents = function(a,gene=FALSE){
	r = c()
	for(ss in c('ad','aa','dd','da')){
		for(c in c('c','p','n'))
			r[paste(ss,c)] = ifelse(gene,length(unique(a$gene_id[a$sites==ss & a$cod==c & a$type!='EXN'])),sum(a$sites==ss & a$cod==c & a$type!='EXN'))
	}
	r
}

# figure 2 #####
caclCor2Embryo = function(d,m,cor.m='p',use.mean.embryo=TRUE){
	cmn = intersect(rownames(m),colnames(d))
	d=d[,cmn]
	m=m[cmn,]
	days = unique(m[,c('days','stage')])
	days = days[order(days$days),]
	tissues = unique(m$tissue)
	r = array(NA,dim = c(length(tissues),nrow(days),3),dimnames = list(tissues,days$stage,c('rho','ci1','ci1')))
	if(use.mean.embryo)
		base = apply(d[,m$stage == days$stage[1]],1,mean,na.rm=T)
	
	for(t in tissues){
		if(!use.mean.embryo){
			mm = m[m$tissue==t,]
			base = d[,rownames(mm)[order(mm$days)[1]]]
		}
		for(s in days$stage){
			ff = m$tissue == t & m$stage == s
			if(sum(ff)!=0){
				tmp = cor.test(base,d[,ff],use='pair',method = cor.m)
				r[t,s,] = c(tmp$estimate,tmp$conf.int)
			}
		}
	}
	r
}

plotCor2Embryo = function(cor,lwd=1,area.transp=0.2,...){
	x = 1:dim(cor)[2]
	plot(1,t='n',xlim=range(x),ylim=range(cor,na.rm=T),xaxt='n',...)
	for(t in 1:dim(cor)[1]){
		plotArea(x,cor[t,,],col=params$tissue.col[dimnames(cor)[[1]][t]],lwd=lwd,area.transp=area.transp)
	}
	axis(1,x,dimnames(cor)[[2]],las=3)
}

plotDevASFreq = function(ens.ids,gids,gene.set.name='',legend.pos=NULL,filter=unique(unlist(ens.ids)),...){
	ens.ids = lapply(ens.ids,function(x){x$d = unique(x$u,x$d);x$u=NULL;list(n=intersect(filter,x$n),d=intersect(filter,x$d))})
	r = lapply(ens.ids,function(x){
		t = unlist(x)
		f = intersect(t,gids) %in% x$d
		r1 = my.binom.test(sum(f),sum(!f))
		f = setdiff(t,gids) %in% x$d
		r2 = my.binom.test(sum(f),sum(!f))
		list(ingenes=r1,outgenes=r2)
	})
	r = list(ingenes = sapply(r,'[[',1),outgenes = sapply(r,'[[',2))
	col=rep(params$tissue.col[colnames(r$ingenes)],each=2)
	b=barplot(rbind(r$ingenes[1,],r$outgenes[1,]),ylim=range(0,r$ingenes,r$outgenes),beside = T,col=col,den=c(-1,40),las=3,ylab='proportion of genes with devAS',border=NA,...)
	segments(b,rbind(r$ingenes[2,],r$outgenes[2,]),b,rbind(r$ingenes[3,],r$outgenes[3,]))#,col=col)
	if(!is.null(legend.pos)){
		legend(legend.pos,col='black',den=c(-1,40),legend=paste0(c('','not '),gene.set.name))
	}
	invisible(r)
}

plotASinGEPatterns = function(up,dw){
	c = rep(params$tissue.col[colnames(up)],each=2)
	b = barplot(up[c(1,4),],beside = T,col=c,den=c(-1,30),las=3,ylim=c(-max(dw[c(3,6),]),max(up[c(3,6),])),ylab='proportion of genes with devAS',main='DevAS in mouse GE clusters',yaxt='n')
	segments(b,up[c(2,5),],b,up[c(3,6),])
	b = barplot(-dw[c(1,4),],beside = T,col=c,den=c(-1,30),add=T,xaxt='n',yaxt='n')
	segments(b,-dw[c(2,5),],b,-dw[c(3,6),])
	abline(h=0)
	at=-1:3*0.05
	axis(2,at,abs(at))
	legend('topright',fill='black',den=c(30,-1),legend=c('Early genes','Late Genes'))
	text(b[2],par('usr')[4],'AS up',adj=c(-0.1,1.1))
	text(b[2],par('usr')[3],'AS down',adj=c(-0.1,-1.1))
}

# figure 4 ####
plotFTMotifSimMatrix = function(m,plot.nums=FALSE,pv.thr=0.01,pv.col=c(gray=1,yellow=1e-2,orange=1e-5,red=1e-20,none=0),or.col=c(blue=0,'#0000FFAA'=1/3,'#0000FF55'=1/1.5,gray=1,yellow=1.5,orange=3,red=Inf),main='',leg.cex=1,diag.col='white',diag.text=colnames(m)){
	plot(1,t='n',xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0,nrow(m)),ylim=c(0,nrow(m)),bty='o',xaxs='i',yaxs='i',main=main)
	for(i in 1:nrow(m)){
		rect(i-1,nrow(m)-i+1,i,nrow(m)-i,col=diag.col[i],border = NA)
		text(i-0.5,nrow(m)-i+0.5,diag.text[i],adj = c(0.5,0.5),font=2,cex=leg.cex)
		if(i == nrow(m)) break
		for(j in (i+1):nrow(m)){
			c = names(pv.col)[min(which(m[i,j] > pv.col))-1]
			rect(j-1,nrow(m)-i+1,j,nrow(m)-i,col=c,border = NA)
			if(plot.nums)
				text(j-0.5,nrow(m)-i+0.5,format(m[i,j],dig=2),adj = c(0.5,0.5))
			c = names(or.col)[min(which(m[j,i] <= or.col))]
			if(m[i,j] > pv.thr)
				c = 'gray'
			rect(i-1,nrow(m)-j+1,i,nrow(m)-j,col=c,border = NA)
			if(plot.nums){
				#text(i-0.5,nrow(m)-j+0.5,round(log2(m[j,i]),2),adj = c(0.5,0.5))
				text(i-0.5,nrow(m)-j+0.5,round(m[j,i],1),adj = c(0.5,0.5))
			}
		}
	}
}


plotMirroredMotFreq = function(f,sign,m,tissue,no.change.in.tiss=NULL,main='',leg=c('inclusion','exclusion'),plot.leg=TRUE,col=c('#FF000077','#0000FF77')){
	if(is.null(no.change.in.tiss)){
		fu = unlist(lapply(names(f),function(s){f[[s]][sign[[s]][names(f[[s]]),tissue]=='u']}))
		fd = unlist(lapply(names(f),function(s){f[[s]][sign[[s]][names(f[[s]]),tissue]=='d']}))
	}else{
		fu = unlist(lapply(names(f),function(s){f[[s]][sign[[s]][names(f[[s]]),tissue]=='u' & sign[[s]][names(f[[s]]),no.change.in.tiss] %in% c('n','-')]}))
		fd = unlist(lapply(names(f),function(s){f[[s]][sign[[s]][names(f[[s]]),tissue]=='d' & sign[[s]][names(f[[s]]),no.change.in.tiss] %in% c('n','-')]}))
	}
	
	p1u  = getPosForMot(fu,m,1,200,0:20*10,TRUE)
	p2u  = getPosForMot(fd,m,1,200,0:20*10,TRUE)
	p1d  = getPosForMot(fu,m,-1,-200,0:20*10,TRUE)
	p2d  = getPosForMot(fd,m,-1,-200,0:20*10,TRUE)
	
	p1eu  = getPosForMot(fu,m,201,250,0:5*10,TRUE)
	p2eu  = getPosForMot(fd,m,201,250,0:5*10,TRUE)
	p1ed  = getPosForMot(fu,m,-201,-250,0:5*10,TRUE)
	p2ed  = getPosForMot(fd,m,-201,-250,0:5*10,TRUE)
	
	space = c(rep(0,19),0.2,rep(0,4),0.2,rep(0,4),0.2,rep(0,20))
	ymax=max(p1u,p1eu,p1ed,p1d,p2u,p2eu,p2ed,p2d)
	b=barplot(c(p1u,p1eu,p1ed,p1d),main=main,space = space,border=NA,ylim=range(0,p1u,p1eu,p1ed,p1d,p2u,p2eu,p2ed,p2d),col=col[1])
	barplot(c(p2u,p2eu,p2ed,p2d),add=T,col=col[2],space = space,border=NA)
	if(plot.leg)
		legend('topright',fill=col,legend = leg,bty='n')
	segments(b[1],-0.05*ymax,b[50],-0.05*ymax)
	rect((b[19]+b[20])/2,-0.09*ymax,(b[29]+b[30])/2,-0.01*ymax,col='red')
	axis(1,c((b[1]+b[19])/2,(b[19]+b[20])/2,(b[29]+b[30])/2,(b[31]+b[50])/2),c("-100","exon","exon","+100"),las=1)
}

getPosForMot = function(fa,m,from,to,breaks,rate=FALSE){
	if(from>0)
		fa = substr(fa,from,to)
	else{
		l=nchar(fa)
		fa = substr(fa,l+to+1,l+from+1)
	}
	r = gregexpr(m,tolower(fa))
	r = unlist(r)
	r = r[r != -1]
	r = hist(r,breaks,plot=F)$counts
	if(rate) return(r/length(fa))
	r
}

plotQKICov = function(l){
	h = c(1,1,1,0.25,0.25,0.25,0.25,0.2,0.8)
	plot(1,t='n',xlim=qki.coor,ylim=c(0,sum(h)),xlab='chr6',ylab='',main='Human QKI splicing',yaxt='n')
	x=c(qki.coor[1],qki.coor[1]:qki.coor[2],qki.coor[2])
	plotCov(x,qki.brain,sum(h[-1]),h[1],qki.coor)
	plotCov(x,qki.hearte,sum(h[-(1:2)]),h[2],qki.coor)
	plotCov(x,qki.hearta,sum(h[-(1:3)]),h[3],qki.coor)
	
	text(qki.coor[1],rev(cumsum(rev(h)))[c(1:4,6)],
			 c('adult brain RNA-Seq','embryo heart RNA-Seq','adult heart RNA-Seq','ENCODE QKI-eCLIP (HEPG2)','ENCODE QKI-eCLIP (K562)'),
			 adj = c(0,1))
	
	param=ScanBamParam(which=GRanges("chr6", IRanges(qki.coor[1], qki.coor[2])))
	eclips = c('ENCFF179YDO.bam','ENCFF081ZJS.bam','ENCFF372LCI.bam','ENCFF494AKD.bam')
	for(i in 1:4){
		b = readGAlignments(paste0('processed/QKI.eCLIP/',eclips[i]),param = param)
		c = as.numeric(coverage(b)[['chr6']][qki.coor[1]:qki.coor[2]])
		lines(x,c(0,c,0)/max(c)*h[3+i]+sum(h[-(1:(3+i))]),lwd=2,col='magenta')
	}
	
	qki.sites = rbind(read.table('processed/QKI.eCLIP/466_01.basedon_466_01.peaks.l2inputnormnew.bed.compressed.bed.narrowPeak.bed'),
										read.table('processed/QKI.eCLIP/466_02.basedon_466_02.peaks.l2inputnormnew.bed.compressed.bed.narrowPeak.bed'),
										read.table('processed/QKI.eCLIP/ENCFF454GLW.bed'),
										read.table('processed/QKI.eCLIP/ENCFF544QKV.bed'))
	
	qki.sites = qki.sites[qki.sites$V1=='chr6' & qki.sites$V3>qki.coor[1] & qki.sites$V2<qki.coor[2],]
	rect(qki.sites$V2,h[9]+h[8]*0.2,qki.sites$V3,h[9]+h[8]*0.8,border = NA,col='magenta')
	text(qki.coor[1],h[9]+h[8]*0.5,'QKI-eCLIP peaks:',adj=c(0,0.5))
	# plot annotation
	abline(h=h[9]/2)
	for(i in 1:nrow(human.qki.gtf)){
		rect(human.qki.gtf$V4[i],h[9]*0.2,human.qki.gtf$V5[i],h[9]*0.8,border = ifelse(human.qki.gtf$V3[i]=='CDS',NA,'black'),col = ifelse(human.qki.gtf$V3[i]=='CDS','red','white'))
	}
	rect(163986978,0,163987752,sum(h),col=NA,lty=2)
	plotPanelLetter(l)
}

plotCov = function(x,cov,fromh,height,region,min.rel.junc.cov=0.05){
	y=c(0,as.numeric(cov$cov),0)
	maxcov = max(y)
	y=y/maxcov*height + fromh
	polygon(x,y,col = 'gray',border = NA)
	
	j = cov$juncs[cov$juncs$score >= min.rel.junc.cov*maxcov &
									((cov$juncs$start>=region[1] & cov$juncs$start<=region[2]) |
									 	(cov$juncs$end  >=region[1] & cov$juncs$end  <=region[2])),]
	if(nrow(j)>0)
		for(i in 1:nrow(j))
			plotArc(j$start[i],j$end[i],j$score[i]/maxcov*height,col='red',lwd=3,y.base=fromh)
}

# figure 3 #### 
caclSegOverlap = function(d){
	overlap = tounion = tomax = tomin = or = pv = matrix(NA,ncol=ncol(d),nrow=ncol(d),dimnames = list(colnames(d),colnames(d)))
	for(t1 in 1:ncol(d)){
		overlap[t1,t1] = sum(d[,t1],na.rm=TRUE)
		if(t1 < ncol(d))
			for(t2 in (t1+1):ncol(d)){
				f = !is.na(d[,t1]) & !is.na(d[,t2])
				overlap[t2,t1] = overlap[t1,t2] = sum(d[f,t1] & d[f,t2])
				tounion[t2,t1] = tounion[t1,t2] = overlap[t1,t2]/sum(d[f,t1] | d[f,t2])
				tomax[t2,t1]   = tomax[t1,t2]   = overlap[t1,t2]/max(sum(d[f,t1]), sum(d[f,t2]))
				tomin[t2,t1]   = tomin[t1,t2]   = overlap[t1,t2]/min(sum(d[f,t1]), sum(d[f,t2]))
				ft = fisher.test(table(d[f,t1], d[f,t2]))
				pv[t2,t1] = pv[t1,t2] = ft$p.value
				or[t2,t1] = or[t1,t2] = ft$estimate
			}
	}
	list(overlap = overlap, tounion = tounion, tomax = tomax, tomin = tomin, or = or,pv=pv)
}

plotASSegStat = function(x,y,pv,col,xax,lty,pv.at=x[-1:-2],pv.thr=c(' '=1,' '=0.05,'*'=0.01,'**'=0.001,'***'=0),...){
	ylim = range(y)
	ylim[2] = ylim[2]+(ylim[2]-ylim[1])*0.05
	plot(x,y[,2],col=col,xlab='',ylim=ylim,xaxt='n',...)
	segments(x,y[,1],x,y[,3],col=cols,lty=lty)
	abline(h=y[1,2],lty=3,col=col[1])
	abline(h=y[2,2],lty=3,col=col[2])
	axis(1,xax,names(xax),las=3)
	pv.thr = sort(pv.thr)
	pv.sym = sapply(pv,function(x)names(pv.thr)[findInterval(x,pv.thr,all.inside = T)])
	d = (ylim[2]-ylim[1])*0.05/3
	for(i in 1:7)
		if(pv[i] < max(pv.thr[pv.thr<1])){
			segments(pv.at[i*2-1]  ,ylim[2]-d,pv.at[i*2-1]  ,y[i*2+1,3]+d)
			segments(pv.at[i*2-1]  ,ylim[2]-d,pv.at[i*2],ylim[2]-d)
			segments(pv.at[i*2],ylim[2]-d,pv.at[i*2],y[i*2+2,3]+d)
			text(pv.at[i*2-1]/2+pv.at[i*2]/2,ylim[2]-d,pv.sym[i],adj=c(0.5,0))
		}
}

plotAgeSegOverlap = function(t,main='Overlap of ageAS exons across mouse tissues',xaxlab=NULL,yaxlab=NULL,...){
	d = t$tomin*100
	d[upper.tri(d)] = t$or[upper.tri(d)]
	d = round(d,1)
	diag(d) = sapply(strsplit(colnames(d),' '),function(x)paste0(substr(x[1],1,1),'\n',x[2]))
	pv.breaks = c(0,1e-20,1e-10,0.05)
	pv.cols = c('red','orange','yellow','gray')
	imageWithText(t$pv,t =d ,xaxt='n',yaxt='n',breaks=c(pv.breaks/7/13,2),col=pv.cols,bty='n',xlab='',ylab='',main=main,xaxlab=xaxlab,yaxlab=yaxlab,...)
	abline(v=ncol(d)/2+.5);abline(h=ncol(d)/2+.5)
	mtext('overlap/min (%)',1,adj = 1)
	mtext('odds ratio',3,adj = 0)
	mtext('inclusion',2,adj = 0.25)
	mtext('exclusion',2,adj = 0.75)
	legend(par('usr')[1],par('usr')[4],xpd=T,fill=pv.cols,legend = paste('<',c(pv.breaks[-1],1)),title = 'Bonf. pv',xjust=1)
}

plotTisUpDownCOns = function(u,d,panel.l,...){
	x = 1:ncol(u)
	plotArea(x,t(u[5:7,]),col='red',new=T,ylim=range(u[5:7],d[5:7,]),xaxt='n',xlab='',ylab='proportion of conserved',...)
	plotArea(x,t(d[5:7,]),col='blue')
	axis(1,x,colnames(u),las=3)
	plotPanelLetter(panel.l)
}


# figure 5 #####
calcBootstrapSpeciesDiv = function(d,sp,fun,age.al,N=100,qv=NULL,cor.meth='pearson',sign.both=FALSE){
	ts = unique(meta$tissue)
	if(N>0){
		perms = lapply(1:N,function(i)sample(nrow(d[[1]]),replace = TRUE))
		r = array(NA,dim=c(length(ts),nrow(age.al),N+1),dimnames = list(ts,age.al$mouse,c('obs',1:N)))
	}else
		r = array(NA,dim=c(length(ts),nrow(age.al),N+1),dimnames = list(ts,age.al$mouse,c('obs')))
	for(t in ts){
		if(!is.null(qv)){
			if(sign.both)
				sgn = qv[[sp[1]]][,t] <0.05 & qv[[sp[2]]][,t] < 0.05
			else
				sgn = qv[[sp[1]]][,t] <0.05 | qv[[sp[2]]][,t] < 0.05
			dd = lapply(d[sp],function(x)x[sgn,])
		}else
			dd = d[sp]
		for(s in 1:nrow(age.al)){
			if(sum(age.al[s,sp]=='') == 0){
				r[t,as.character(age.al$mouse[s]),1] = fun(getSpeciesCor(dd,age.al[s,sp],t,method = cor.meth))
			}
		}
	}
	if(N>0)
		for(i in 1:N){
			dd = lapply(d,function(x)x[perms[[i]],])
			for(t in ts){
				for(s in 1:nrow(age.al)){
					if(sum(age.al[s,sp]=='') == 0){
						cat('\r',i,t,s,'         ')
						r[t,as.character(age.al$mouse[s]),1+i] = fun(getSpeciesCor(dd,age.al[s,sp],t,method = cor.meth))
					}
				}
			}
		}
	r
}

getSpeciesCor = function(psi,stages,tissue,f=rep(TRUE,nrow(psi[[1]])),method='pearson'){
	m = matrix(NA,ncol=length(stages),nrow=sum(f))
	colnames(m) = names(stages)
	for(s in names(stages)){
		i = paste(s,tissue,stages[[s]])
		if(i %in% colnames(psi[[s]]))
			m[,s] = psi[[s]][f,i]
	}
	cor(m,u='p',method = method)
}


plotDivergenceOnAgeWithConf = function(d,col,ci.pv=0.05,lwd=3,straight.line=FALSE,plot.points=FALSE,...){
	x = 1:dim(d)[2]
	if(straight.line)
		r = lapply(1:dim(d)[1],function(t){
			v = d[t,,1]
			f = !is.na(v)
			xx = x[f]
			m = lm(v[f] ~ xx)
			r=predict(m,interval = 'conf',newdata=list(xx=x))
			r[!f,] = NA
			r
		})
	else
		r = lapply(1:dim(d)[1],function(t){
			cbind(d[t,,1],t(apply(d[t,,-1],1,quantile,probs=c(ci.pv/2,1-ci.pv/2),na.rm=TRUE)))
		})
	plot(1,t='n',xlim=c(1,dim(d)[2]),ylim=range(unlist(r),na.rm=T),xaxt='n',xlab='Mouse stage',...)
	
	for(t in 1:length(r)){
		f = !is.na(r[[t]][,1])
		plotArea(x[f],r[[t]][f,],new = FALSE,col=col[dimnames(d)[[1]][t]],lwd=lwd)
		if(plot.points)
			points(x[f],d[t,f,1],pch=19,col=col[dimnames(d)[[1]][t]])
	}
	axis(1,x,dimnames(d)[[2]],las=3)
}
plotSpeciesTree.3 = function(l,patt,mean,birt,lost){
	cnts = rbind(mean,patt,birt,lost)
	colnames(cnts) = names(birt)
	ff = apply(cnts=='',1,sum) == 0
	cnts = cnts[ff,]
	text.col=c('red','orange','blue','green')[ff]
	
	plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(1,11),ylim=c(0,7))
	title('Number of evolutionary events')
	segments(1,1,1,6)
	segments(2,2,2,6)
	segments(3,4,3,6)
	segments(4:7,5,4:7,6)
	segments(4.5,4,4.5,5)
	segments(3.75,3,3.75,4)
	segments(6.5,3,6.5,5)
	segments(5.125,2,5.125,3)
	segments(7.125/2,1,7.125/2,2)
	
	segments(1,1,7.125/2,1)
	segments(2,2,5.125,2)
	segments(3.75,3,6.5,3)
	segments(3,4,4.5,4)
	segments(4,5,5,5)
	segments(6,5,7,5)
	text(1:7,6,rev(rownames(species)),adj = c(0,0),srt=50)
	
	y = (0:4/4)[-5][1:nrow(cnts)]
	text(1.05,4-y,cnts[,'c'],col = text.col,adj = c(0,0))
	text(2.05,4.5-y,cnts[,'o'],col = text.col,adj = c(0,0))
	text(3.05,5.5-y,cnts[,'b'],col = text.col,adj = c(0,0))
	text(4.05,6-y,cnts[,'r'],col = text.col,adj = c(0,1))
	text(5.05,6-y,cnts[,'m'],col = text.col,adj = c(0,1))
	text(6.05,6-y,cnts[,'q'],col = text.col,adj = c(0,1))
	text(7.05,6-y,cnts[,'h'],col = text.col,adj = c(0,1))
	
	text(6.55,4.5-y,cnts[,'hq'],col = text.col,adj = c(0,1))
	text(4.55,5-y,cnts[,'mr'],col = text.col,adj = c(0,1))
	text(3.8,4-y,cnts[,'mrb'],col = text.col,adj = c(0,1))
	text(5.175,3-y,cnts[,'hqmrb'],col = text.col,adj = c(0,1))
	
	text(6,2-y,c('mean','pattern','birth','lost')[ff],col = text.col,adj = c(0,1))
	plotPanelLetter(l)
}

# figure 6 #####
getSgnFraqWithBootstrap = function(spec.sp,sps,sgn,N=100,s2s=setNames(rownames(species),species$short)){
	spec.sp.list = lapply(strsplit(spec.sp,''),function(x){intersect(s2s[x],colnames(sgn))})
	sgn.pr = sapply(1:nrow(sgn),function(i){
		mean(sgn[i,spec.sp.list[[i]]],na.rm=TRUE)
	})
	f = function(spec.sp,sgn.pr){
		sgn.pr = split(sgn.pr,spec.sp)
		sapply(sps,function(s)mean(unlist(sgn.pr[s]),na.rm=T))
	}
	r = f(spec.sp,sgn.pr)
	#bootstrap
	bs=sapply(1:N,function(i){
		t = sample(1:nrow(sgn),replace = TRUE)
		f(spec.sp[t],sgn.pr[t])
	})
	cnt = sapply(sps,function(s)sum(spec.sp %in% s))
	cbind(mean=r,cnt=cnt,t(apply(bs,1,quantile,probs=c(0.025,0.975),na.rm=T)))
}
