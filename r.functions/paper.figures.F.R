
calcDevVsNdevCor = function(s1,s2,maxPSI=0.9,per.tissue.sgn=TRUE,cor.meth='sp',psi.in.both=TRUE){
	a = orth.seg.ad[[s1]]$seg$type=='ALT' & orth.seg.ad[[s2]]$seg$type=='ALT'
	r = NULL
	for(t in unique(meta$tissue)){
		f1 = a & orth.seg.ad.tsm[[s1]][,paste(s1,t,border.stages[[s1]][t,2])] <= maxPSI
		f2 = a & orth.seg.ad.tsm[[s2]][,paste(s2,t,border.stages[[s2]][t,2])] <= maxPSI
		if(psi.in.both)
			f = f1 & f2
		else
			f = f1 | f2
		f[is.na(f)] = FALSE
		if(per.tissue.sgn){
			d1 = orth.per.tissue.age.qv[[s1]][,t] < 0.05 & abs(orth.age.dpsi[[s1]][,t]) > 0.2
			d2 = orth.per.tissue.age.qv[[s2]][,t] < 0.05 & abs(orth.age.dpsi[[s2]][,t]) > 0.2
		}else{
			d1 = apply(orth.per.tissue.age.qv[[s1]] < 0.05 & abs(orth.age.dpsi[[s1]]) > 0.2,1,sum,na.rm=T) > 0
			d2 = apply(orth.per.tissue.age.qv[[s2]] < 0.05 & abs(orth.age.dpsi[[s2]]) > 0.2,1,sum,na.rm=T) > 0
		}
		d1 = !is.na(d1) & d1
		d2 = !is.na(d2) & d2
		for(s in 1:nrow(age.al.i)){
			sid1 = paste(s1,t,age.al.i[s,s1])
			sid2 = paste(s2,t,age.al.i[s,s2])
			if(sid1 %in% colnames(orth.seg.ad.tsm[[s1]]) & sid2 %in% colnames(orth.seg.ad.tsm[[s2]])){
				dev = cor(orth.seg.ad.tsm[[s1]][f &  (d1 | d2),sid1],
									orth.seg.ad.tsm[[s2]][f &  (d1 | d2),sid2],u='p',m=cor.meth)
				ndv = cor(orth.seg.ad.tsm[[s1]][f & !(d1 | d2),sid1],
									orth.seg.ad.tsm[[s2]][f & !(d1 | d2),sid2],u='p',m=cor.meth)
				r = rbind(r,data.frame(tissue=t,mouse.stage=age.al.i$mouse[s],devAS.cor = dev,ndevAS.cor = ndv))
			}
		}
	}
	r$sid = paste('mouse',r$tissue,r$mouse.stage)
	r
}

getSpDivFromCorMatrix = function(c,a,m){
	ts = unique(m$tissue)
	r = list()
	for(t in ts){
		r[[t]] = list()
		for(i in 1:nrow(a)){
			t1 = rownames(m)[m$species==colnames(a)[1] & m$tissue==t & m$stage==a[i,1]]
			t2 = rownames(m)[m$species==colnames(a)[2] & m$tissue==t & m$stage==a[i,2]]
			if(length(t1)>0 & length(t2)>0){
				r[[t]][[a[i,1]]] = c[t1,t2]
			}else
				r[[t]][[a[i,1]]] = c()
		}
	}
	r
}

plotSpCor2Cor = function(c1,c2,a,m,same.lims=FALSE,...){
	f=function(x){
		do.call(rbind,lapply(names(x),function(z){
			do.call(rbind,lapply(names(x[[z]]),function(y){
				data.frame(as.numeric(x[[z]][[y]]),z,y)}))
		}))
	}
	d1  = f(getSpDivFromCorMatrix(c1,a,m))
	d2  = f(getSpDivFromCorMatrix(c2,a,m))
	cols = unique(m[,c('tissue','col')])
	cexs = unique(m[m$stage %in% d1[,3] & m$species == colnames(a)[1],c('stage','cex')])
	cols = setNames(cols[,2],cols[,1])[d1[,2]]
	cs = setNames(cexs[,2],cexs[,1])[d1[,3]]
	cs = (cs-min(cs))/(max(cs)-min(cs))*1+0.05
	if(same.lims){
		xlim=ylim=range(d1[,1],d2[,1])
	}else{
		xlim=range(d1[,1])
		ylim=range(d2[,1])
	}
	plot(d1[,1],d2[,1],col=cols,cex=cs*1.5,pch=19,xlim=xlim,ylim=ylim,...)
	invisible(cbind(d1,d2))
}




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

plotTisUpDownCOns = function(u,d,panel.l,...){
	x = 1:ncol(u)
	plotArea(x,t(u[5:7,]),col='red',new=T,ylim=range(u[5:7],d[5:7,]),xaxt='n',xlab='',ylab='proportion of conserved',...)
	plotArea(x,t(d[5:7,]),col='blue')
	axis(1,x,colnames(u),las=3)
	plotPanelLetter(panel.l)
}

plotASSegStat = function(x,y,pv,col,xax,lty,pv.at=x[-1],pv.thr=c(' '=1,' '=0.05,'*'=0.01,'**'=0.001,'***'=0),...){
	ylim = range(y)
	ylim[2] = ylim[2]+(ylim[2]-ylim[1])*0.05
	plot(x,y[,2],col=col,xlab='',ylim=ylim,xaxt='n',...)
	segments(x,y[,1],x,y[,3],col=cols,lty=lty)
	abline(h=y[1,2],lty=3,col=col[1])
	axis(1,xax,names(xax),las=3)
	pv.thr = sort(pv.thr)
	pv.sym = sapply(pv,function(x)names(pv.thr)[findInterval(x,pv.thr,all.inside = T)])
	d = (ylim[2]-ylim[1])*0.05/3
	for(i in 1:7)
		if(pv[i] < max(pv.thr[pv.thr<1])){
			segments(pv.at[i*2-1]  ,ylim[2]-d,pv.at[i*2-1]  ,y[i*2,3]+d)
			segments(pv.at[i*2-1]  ,ylim[2]-d,pv.at[i*2],ylim[2]-d)
			segments(pv.at[i*2],ylim[2]-d,pv.at[i*2],y[i*2+1,3]+d)
			text(pv.at[i*2-1]/2+pv.at[i*2]/2,ylim[2]-d,pv.sym[i],adj=c(0.5,0))
		}
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

calcBootstrapSpeciesDiv = function(d,sp,fun,age.al,N=100,qv=NULL){
	ts = unique(meta$tissue)
	if(N>0){
		perms = lapply(1:N,function(i)sample(nrow(d[[1]]),replace = TRUE))
		r = array(NA,dim=c(length(ts),nrow(age.al),N+1),dimnames = list(ts,age.al$mouse,c('obs',1:N)))
	}else
		r = array(NA,dim=c(length(ts),nrow(age.al),N+1),dimnames = list(ts,age.al$mouse,c('obs')))
	for(t in ts){
		if(!is.null(qv)){
			sgn = qv[[sp[1]]][,t] <0.05 | qv[[sp[2]]][,t] < 0.05
			dd = lapply(d[sp],function(x)x[sgn,])
		}else
			dd = d[sp]
		for(s in 1:nrow(age.al)){
			if(sum(age.al[s,sp]=='') == 0){
				r[t,as.character(age.al$mouse[s]),1] = fun(getSpeciesCor(dd,age.al[s,sp],t))
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
						r[t,as.character(age.al$mouse[s]),1+i] = fun(getSpeciesCor(dd,age.al[s,sp],t))
					}
				}
			}
		}
	r
}


getNoOfEvents = function(a,gene=FALSE){
	r = c()
	for(ss in c('ad','aa','dd','da')){
		for(c in c('c','p','n'))
			r[paste(ss,c)] = ifelse(gene,length(unique(a$gene_id[a$sites==ss & a$cod==c & a$type!='EXN'])),sum(a$sites==ss & a$cod==c & a$type!='EXN'))
	}
	r
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

plotOhnologsFreq = function(ens.ids,gids,filter=unique(unlist(ens.ids)),den=c(n=0,u=40,d=40),angle=c(45,45,-45),ylab='proportion of ohnologs',legend=FALSE,merge.up.down=FALSE,...){
	if(merge.up.down){
		ens.ids = lapply(ens.ids,function(x){x$d = unique(x$u,x$d);x$u=NULL;x})
		den=c(n=20,d=-1)
		angle = 45
	}
	r = lapply(ens.ids,function(x)lapply(x,function(z){
		t = intersect(z,filter) %in% gids
		r=my.binom.test(sum(t),sum(!t))
		names(r) = NULL
		r
	}))
	
	p = sapply(r,function(x)sapply(x,'[',1))[names(den),]
	ci1 = sapply(r,function(x)sapply(x,'[',2))[names(den),]
	ci2 = sapply(r,function(x)sapply(x,'[',3))[names(den),]
	
	
	col=rep(params$tissue.col[colnames(p)],each=nrow(p))
	tis = colnames(p)
	colnames(p) = NULL
	b=barplot(p,beside = T,ylab=ylab,col=col,border = col,den=den,angle = angle,ylim=range(0,ci2,na.rm=T),las=3,...)
	axis(1,b,rep(names(den),times=ncol(p)),tck=-0.01,mgp=c(0.9,0.0,0),cex.axis=0.5)
	axis(1,apply(b,2,mean),tis,tck=0,las=3,mgp=c(0,1.0,0))
	segments(b,ci1,b,ci2)
	if(legend){
		if(merge.up.down){
			legend('topleft',fill='black',den=den*2,angle=angle,legend=c('non-devAS','devAS'),title='AS direction')
		}else
			legend('topleft',fill='black',den=den*2,angle=angle,legend=c('no change','up','down'),title='AS direction')
	}
	invisible(r)
}

plotSpeciesDiv.fig3 = function(cors,species,tree.xlim,tree.ylim,max=1,trees.stages = list(tissue=c('brain','liver'),stages=c('11.5','9wpb')),...){
	trees = list(nj(max-cors[[trees.stages$tissue[1]]][[trees.stages$stages[1]]][species,species]),
							 nj(max-cors[[trees.stages$tissue[1]]][[trees.stages$stages[2]]][species,species]),
							 nj(max-cors[[trees.stages$tissue[2]]][[trees.stages$stages[1]]][species,species]),
							 nj(max-cors[[trees.stages$tissue[2]]][[trees.stages$stages[2]]][species,species]))
	
	plotDivergenceOnAge(cors,species=species,ylab='mean Pearson correlation',...)
	plotPanelLetter('A')
	p = par(mar=c(3,3,1.5,0))
	plot(trees[[1]],type = 'unrooted',x.lim=tree.xlim,y.lim=tree.ylim,main=paste(trees.stages$tissue[1],trees.stages$stages[1]))
	plotPanelLetter('B')
	plot(trees[[2]],type = 'unrooted',x.lim=tree.xlim,y.lim=tree.ylim,main=paste(trees.stages$tissue[1],trees.stages$stages[2]))
	plot(trees[[3]],type = 'unrooted',x.lim=tree.xlim,y.lim=tree.ylim,main=paste(trees.stages$tissue[2],trees.stages$stages[1]))
	plot(trees[[4]],type = 'unrooted',x.lim=tree.xlim,y.lim=tree.ylim,main=paste(trees.stages$tissue[2],trees.stages$stages[2]))
	par(p)
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
plotDivergenceOnAge = function(cor,fun=function(x)mean(x[upper.tri(x)]),species=colnames(cor[[1]][[1]]),ltype='l',lwd=1,...){
	r = sapply(cor,function(x)sapply(x,function(y)fun(y[species,species])))
	plot(1,t='n',ylim=range(r,na.rm=T),xaxt='n',xlab='Mouse age',xlim=c(1,nrow(r)),...)
	x = 1:nrow(r)
	for(t in colnames(r)){
		lines(x[!is.na(r[,t])],r[x[!is.na(r[,t])],t],col=params$tissue.col[t],t=ltype,lwd=lwd)
	}
	axis(1,x,rownames(r),las=3)
}


insertIntoText = function(t,i,each,max.chunks=NULL){
	s = seq(1,nchar(t),by = each)
	r = sapply(s,function(p)substr(t,p,p+each-1))
	if(!is.null(max.chunks) && length(r) > max.chunks){
		r = r[1:max.chunks]
		substr(r[length(r)],each-2,each) = '...'
	}
	paste(r,collapse=i)
}

wrapText = function(t,maxlen,max.chunks=Inf){
	t = strsplit(t,' ',fixed = TRUE)[[1]]
	if(length(t)==1)
		return(t)
	max.chunks = min(max.chunks,length(t))
	r = rep('',max.chunks)
	r[1] = paste0(t[1],' ')
	j = 1
	for(i in 2:length(t)){
		if(nchar(r[j]) + nchar(t[i])+1 <= maxlen)
			r[j] = paste0(r[j],t[i],' ')
		else{
			r[j] = substr(r[j],1,nchar(r[j])-1)
			j = j+1
			if(j <= length(r)){
				r[j] = paste0(t[i],' ')
			}else{
				r[j-1] = paste0(r[j-1],'...')
				break;
			}
		}
	}
	paste(r[r!=''],collapse='\n')
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

makeNJTreesByCorLastTopology = function(x){
	t1 = lapply(x,function(c){
		t = root(nj(as.dist(1-c)),outgroup='base',resolve.root=T)
		t$edge.length[t$edge.length<0]=0
		t
	})
	l = t1[[length(t1)]]
	
	t2 = lapply(names(t1),function(d){
		cmn = intersect(t1[[d]]$tip.label,l$tip.label)
		optim.phylo.wls(1-x[[d]][cmn,cmn],stree=drop.tip(l,setdiff(l$tip.label,cmn)),fixed=TRUE,collapse=FALSE)
	})
	
	names(t2) = names(x)
	t2
}

getSNPDistr = function(g,sids){
	g = g[g$alt_cnt>1,]
	t = lapply(sids,function(x)log10(g$freq[g$seg_id %in% x]))
	r = sapply(t,function(x){m=mean(x);s=sd(x)/sqrt(length(x));c(m,m-2*s,m+2*s)})
	r = rbind(r,pv=sapply(1:ncol(r),function(x)wilcox.test(t[[1]],t[[x]],a='g')$p.value),snp.cnt = sapply(t,length),seg.cnt = sapply(sids,length))
	r
}
plotSNPFreq = function(d,cols,ylim,pv.thr=0.05,...){
	x = 1:ncol(d)
	plot(x,d[1,],col=cols,pch=ifelse(d[4,]<pv.thr,19,1),ylab='log10(SNP freq)',xlab='',xaxt='n',ylim=ylim,...)
	segments(x,d[2,],x,d[3,],col=cols)
	axis(1,x,paste0(colnames(d),'\n',d[5,]),las=3)
}

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
getSpSpPerTissue = function(d,othr,rthr){
	r = matrix(NA,ncol=length(d),nrow=nrow(d[[1]]),dimnames = list(rownames(d[[1]]),names(d)))
	for(t in names(d)){
		r[!is.na(d[[t]]$out),t] = 'none'
		f = !is.na(d[[t]]$out) & d[[t]]$out >= othr & d[[t]]$out/d[[t]]$within > rthr
		r[f,t] = d[[t]]$species[f]
	}
	r
}

plotOddsRatioBars = function(d,filts,L,...){
	f = function(d,fs){
		apply(fs,2,function(f){
			t=fisher.test(d,f)
			c(t$estimate,t$conf.int,t$p.value,sum(f))
		})
	}
	t = f(d, filts)
	b=barplot(log2(t[1,]),ylim=log2(range(t[2:3,])),las=3,names.arg = paste0(colnames(t),'\n(',t[5,],')'),ylab='log2(odds ratio)',...)
	segments(b,log2(t[2,]),b,log2(t[3,]))
	plotPanelLetter(L)
}

caclCramersVPerCols = function(d,pvalue=FALSE){
	r = matrix(NA,ncol=ncol(d),nrow=ncol(d),dimnames = list(colnames(d),colnames(d)))
	for(i in 1:(ncol(d)-1))
		for(j in (i+1):ncol(d)){
			if(pvalue)
				r[i,j]=r[j,i] = chisq.test(d[,i],d[,j])$p.value
			else
				r[i,j]=r[j,i] = cramersV(d[,i],d[,j])
		}
	r
}

imageWithText = function(d,t=NULL,digits=2,text.col=NULL,xaxlab=rownames(d),yaxlab=colnames(d),las=1,...){
	if(is.null(t))
		t = round(d,digits = digits)
	x = 1:nrow(d)
	y = 1:ncol(d)
	image(x,y,d,xaxt='n',yaxt='n',...)
	if(is.null(text.col))
		text.col = 'black'
	text(rep(x,times=length(y)),rep(y,each=length(x)),t,col=text.col)
	if(!is.null(xaxlab))
		axis(1,x,xaxlab,las=las)
	if(!is.null(yaxlab))
		axis(2,y,yaxlab,las=las)
}

imageSpSpCramer = function(d,...){
	library(lsr)
	cv = caclCramersVPerCols(d,FALSE)
	pv = p.adjust(caclCramersVPerCols(d,TRUE),m='BH')
	imageWithText(cv,text.col=ifelse(pv<0.05,'black','gray'),yaxt='n',xaxt='n',xlab='',ylab='',...)
	axis(1,1:10,paste(rep(colnames(patt.sp.sp),times=2),rep(c('patt.','mean'),each=5)),las=2)
	axis(2,1:10,paste(rep(colnames(patt.sp.sp),times=2),rep(c('patt.','mean'),each=5)),las=2)
	invisible(list(cv=cv,pv=pv))
}

plotTreeLength = function(tree.lens,plot.leg=TRUE,...){
	plot(1,t='n',xlim=c(1,3),ylim=range(tree.lens),xlab='Age (days from conception)',ylab='Total tree length',xaxt='n',...)
	lines(1:3,tree.lens[,1],t='b',col=params$tissue.col[colnames(tree.lens)[1]],lwd=2)
	lines(1:3,tree.lens[,2],t='b',col=params$tissue.col[colnames(tree.lens)[2]],lwd=2)
	lines(1:3,tree.lens[,3],t='b',col=params$tissue.col[colnames(tree.lens)[3]],lwd=2)
	lines(2:3,tree.lens[2:3,4],t='b',col=params$tissue.col[colnames(tree.lens)[4]],lwd=2)
	axis(1,1:3,gsub('e','',rownames(tree.lens)))
	if(plot.leg)
		legend('bottomright',col=params$tissue.col[colnames(tree.lens.anc)],lwd=2,legend=colnames(tree.lens.anc),bty='n')
}


plotSpeciesTrees = function(data,main=''){
	spec.cor = getAllSpeciesCor(data)
	tree.len = matrix(NA,nrow=length(spec.cor),ncol=length(spec.cor[[1]]),dimnames=list(names(spec.cor),names(spec.cor[[1]])))
	par(mfcol=c(4,3),tck=-0.02,mgp=c(1.3,0.2,0),mar=c(0,0,0,0),oma=c(2,0,3,1))
	for(d in names(spec.cor)[])
		for(t in names(spec.cor[[d]])){
			tree = nj(as.dist(1-spec.cor[[d]][[t]]))
			if('c' %in% tree$tip.label)
				tree=root(tree,outgroup='c',resolve.root=T)
			else if('o' %in% tree$tip.label)
				tree=root(tree,outgroup='o',resolve.root=T)
			else if('h' %in% tree$tip.label)
				tree=root(tree,outgroup='h',resolve.root=T)
			cols = c(m='gray',r='black',h='red',b='brown',c='blue',o='green')
			tree.len[d,t] = sum(tree$edge.length)
			plot.phylo(tree,y.lim=c(-0.165,0.18),cex=1.4,tip.color=cols[tree$tip.label],direction='upwards',type='fan',use.edge.lengt=T)
			#box()
			if(d=='e10.5'){
				mtext(t,2,line = -2)
			}
			if(t=='brain')
				mtext(d,3)
			#axis(1)
		}
	mtext(main,outer = T,line = 1.5)
	tree.len
}

# getAllSpeciesCor = function(d){
# 	ts = c('brain','heart','liver','testis')
# 	r = list(e10.5=list(),e15.5=list(),e82=list())
# 	r$e10.5 = lapply(ts[1:3],function(t)getSpeciesCor(d,'hmrbo',t,10.5))
# 	r$e10.5$testis = getSpeciesCor(d,'rbmh','testis',10.5)
# 	r$e15.5 = lapply(ts,function(t)getSpeciesCor(d,'hmrboc',t,c(15.5,15.5,15.5,15.5,14.5,15.5)))
# 	r$e82 = lapply(ts,function(t)getSpeciesCor(d,'hmrboc',t,82))
# 	for(i in 1:length(r))
# 		names(r[[i]]) = ts
# 	r
# }
# 
# getSpeciesCor = function(d,sps,t,mdays){
# 	sps = setNames(rownames(species),species$short)[strsplit(sps,'')[[1]]]
# 	if(length(mdays)==1)
# 		mdays = rep(mdays,length(sps))
# 	r = NULL
# 	for(i in 1:length(sps)){
# 		r = cbind(r,d[,paste(sps[i],t,mdays[i])])
# 	}
# 	colnames(r) = species[sps,'short']
# 	cor(r,u='p')
# }

compareAnns = function(e,c,bty=NULL,min.int.cnt=0){
	if(is.null(bty)){
		f = TRUE
	}else{
		f = rep(FALSE,nrow(c))
		f[c$ens=='-'] = TRUE
		
		f[c$ens!='-'] = e[c$ens[c$ens!='-']] %in% bty
	}
	
	c = c[(c$my.ints>=min.int.cnt | c$my == '-') & (c$ens.ints>=min.int.cnt | c$ens == '-') & f,]
	my0 = sum(c$my=='-')
	ens0 = sum(c$ens=='-')
	c = c[c$class %in% c('=','c','j') | (c$class=='e' & c$my.ints==0 & c$ens.ints==0),]
	
	my.cnt = table(c$my)
	ens.cnt = table(c$ens)
	my.stat  = c('0'= my0,as.numeric(table(factor( my.cnt,levels = 1:max(my.cnt)))))
	ens.stat = c('0'=ens0,as.numeric(table(factor(ens.cnt,levels = 1:max(ens.cnt)))))
	one2one= sum(c$my %in% names(my.cnt)[my.cnt==1] & c$ens %in% names(ens.cnt)[ens.cnt==1])
	list(my.stat=my.stat,ens.stat=ens.stat,total=c(my.genes=sum(my.stat),ens.genes=sum(ens.stat),one2one=one2one))
}
plotAnnComp = function(x,...){
	plot(1,t='n',yaxt='n',xlab='# of genes',ylab='number of genes merged together',xlim=c(0,10),ylim=c(0,log(max(x$my.stat,x$ens.stat)+1)),...)
	at = c(0,5,20,100,500,2000,10000)
	axis(2,log(at+1),at)
	lines(0:(length(x$my.stat)-1),log(x$my.stat+1),col='red',pch=19,t='b',lwd=2)
	lines(0:(length(x$ens.stat)-1),log(x$ens.stat+1),col='blue',pch=19,t='b',lwd=2)
	l=legend('topright',title='one-to-one/present in other ann./total',
					 legend=c(paste0('my:  ',x$total[3],'/',x$total[1]-x$my.stat[1],'/',x$total[1]),
					 				 paste0('ens: ',x$total[3],'/',x$total[2]-x$ens.stat[1],'/',x$total[2])))
	legend(l$rect$left,l$rect$top-l$rect$h,pch=19,col=c('red','blue'),legend=c('Merged in my annotation','Splitted Ensambl genes'))
}

plotHexs = function(label){
	f = function(h,n,l){
		h = sort(h)
		r = ''
		for(i in 1:l){
			z=h[(i-1)*n+1:n]
			z=z[!is.na(z)]
			if(length(z)>0){
				if(nchar(r)>0)
					r = paste0(r,'\n',paste(z,collapse=' '))
				else
					r = paste(z,collapse=' ')
			}
		}
		r
	}
	
	minx=-0.6
	plot(1,bty='n',xaxt='n',yaxt='n',xlim=c(minx,5),ylim=c(1,8),t='n',xlab='',ylab='',xaxs='i')
	abline(h=4)
	abline(h=c(1:3,5:7),lty=2)
	rect(2.1,3.5,2.9,4.5,col='red')
	text(minx,6.5:4.5,c('early brain\n(pattern 3)','brain+heart\n(pattern 5)','heart\n(pattern 11)'),adj = c(0,0.5))
	text(minx,1.5:3.5,c('early brain\n(pattern 4)','brain+heart\n(pattern 6)','heart\n(pattern 12)'),adj = c(0,0.5))
	text(minx-.2,5.5,'inclusion',xpd=TRUE,srt=90,adj=c(0.5,0.5),cex=2)
	text(minx-.2,2.5,'exclusion',xpd=TRUE,srt=90,adj=c(0.5,0.5),cex=2)
	n=5
	l=3
	text(3,6.5,f(intersect(hex$upm$`4`,hex$dwm$`3`),n,l),adj=c(0,0.5),family="mono")
	text(2,1.5,f(intersect(hex$upm$`4`,hex$dwm$`3`),n,l),adj=c(1,0.5),family="mono")
	arrows(2,1.5,3,6.5,code=3,angle = 20,length = 0.12)
	
	text(3,2.5,f(intersect(hex$upm$`5`,hex$dwm$`6`),n,l),adj=c(0,0.5),family="mono")
	text(2,5.5,f(intersect(hex$upm$`5`,hex$dwm$`6`),n,l),adj=c(1,0.5),family="mono")
	arrows(2,5.5,3,2.5,code=3,angle = 20,length = 0.12)
	
	text(3,4.5,f(intersect(hex$upm$`12`,hex$dwm$`11`),n,l),adj=c(0,0.5),family="mono")
	text(2,3.5,f(intersect(hex$upm$`12`,hex$dwm$`11`),n,l),adj=c(1,0.5),family="mono")
	arrows(2,3.5,3,4.5,code=3,angle = 20,length = 0.12)
	
	plotPanelLetter(label)
}

getHexAnn = function(hexs,h2m){
	r=h2m[hexs,-2]	
	r=unique(r[r$V3!='',])
	r[order(r$V3),]
}



# old versions
# plotMirroredMotFreq = function(m,c1,c2,main='',leg=c('inclusion','exclusion'),l,plot.leg=TRUE,col=c('#FF000077','#0000FF77')){
# 	p1u  = getPosForMot(ad.fa,orth.ad.cl.reassign,c1,m,1,200,0:20*10,TRUE)
# 	p2u  = getPosForMot(ad.fa,orth.ad.cl.reassign,c2,m,1,200,0:20*10,TRUE)
# 	p1d  = getPosForMot(ad.fa,orth.ad.cl.reassign,c1,m,-1,-200,0:20*10,TRUE)
# 	p2d  = getPosForMot(ad.fa,orth.ad.cl.reassign,c2,m,-1,-200,0:20*10,TRUE)
# 	
# 	p1eu  = getPosForMot(ad.fa,orth.ad.cl.reassign,c1,m,201,250,0:5*10,TRUE)
# 	p2eu  = getPosForMot(ad.fa,orth.ad.cl.reassign,c2,m,201,250,0:5*10,TRUE)
# 	p1ed  = getPosForMot(ad.fa,orth.ad.cl.reassign,c1,m,-201,-250,0:5*10,TRUE)
# 	p2ed  = getPosForMot(ad.fa,orth.ad.cl.reassign,c2,m,-201,-250,0:5*10,TRUE)
# 	
# 	space = c(rep(0,19),0.2,rep(0,4),0.2,rep(0,4),0.2,rep(0,20))
# 	ymax=max(p1u,p1eu,p1ed,p1d,p2u,p2eu,p2ed,p2d)
# 	b=barplot(c(p1u,p1eu,p1ed,p1d),main=paste(main,' (',m,' ',c1,'-',c2,')',sep=''),space = space,border=NA,ylim=range(0,p1u,p1eu,p1ed,p1d,p2u,p2eu,p2ed,p2d),col=col[1])
# 	barplot(c(p2u,p2eu,p2ed,p2d),add=T,col=col[2],space = space,border=NA)
# 	if(plot.leg)
# 		legend('topright',fill=col,legend = leg,bty='n')
# 	segments(b[1],-0.05*ymax,b[50],-0.05*ymax)
# 	rect((b[19]+b[20])/2,-0.09*ymax,(b[29]+b[30])/2,-0.01*ymax,col='red')
# 	axis(1,c((b[1]+b[19])/2,(b[19]+b[20])/2,(b[29]+b[30])/2,(b[31]+b[50])/2),c("-100","exon","exon","+100"),las=1)
# 	plotPanelLetter(l)
# }
# 
# getPosForMot = function(fa,clust,cl,m,from,to,breaks,rate=FALSE){
# 	names(fa) = NULL
# 	fa = unlist(fa)
# 	if(from>0)
# 		fa = substr(fa,from,to)
# 	else{
# 		l=nchar(fa)
# 		fa = substr(fa,l+to+1,l+from+1)
# 	}
# 	s = intersect(names(fa),rownames(clust)[clust$orth.cl == cl])
# 	r = gregexpr(m,tolower(fa[s]))
# 	r = unlist(r)
# 	r = r[r != -1]
# 	r = hist(r,breaks,plot=F)$counts
# 	if(rate) return(r/length(s))
# 	c(length(s),r)
# }
# 




plotNoOfMots = function(ms,kn,main,plot.leg=FALSE,l){
	t = sapply(ms,function(x){table(factor(x %in% kn,levels = c(TRUE,FALSE)))})
	barplot(t[1,]+t[2,],xlab='cluster',ylab='# of hexamers',main=main,col=c(rep(c('#FFAAAA','#AAAAFF'),each=1,times=8),rep('#CCCCCC',20)))
	barplot(t['TRUE',],yaxt='n',xlab='cluster',ylab='# of hexamers',main=main,col=c(rep(c('#FF0000','#0000FF'),each=1,times=8),rep('#888888',20)),add=TRUE)
	if(plot.leg)
		legend('topright',fill=c('#FFAAAA','#AAAAFF','#CCCCCC','#888888'),ncol = 2,legend=c('inclusion','exclusion','unknown','known'),bty='n')
	plotPanelLetter(l)
}

getPolar = function(a,r,center=c(0,0)){
	list(x=center[1]+cos(a)*r,y=center[2]+sin(a)*r)
}

psichart = function(x,labels = names(x),radius=0.4,spiral.shift=0,init.angle=0,col=rainbow(length(x)),border=FALSE,main='',center=c(0.5,0.5),tck.length = radius/10,N=200,add=FALSE){
	if(!add){
		plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,1),ylim=c(0,1))
		title(main=main)
	}
	x = x/sum(x)*2*pi
	angle = init.angle
	for(i in 1:length(x)){
		a = c(seq(from=angle,to=angle+spiral.shift,length.out = N),
					seq(from=angle+spiral.shift,to=angle+spiral.shift+x[i],length.out = N),
					seq(from=angle+spiral.shift+x[i],to=angle+x[i],length.out = N))
		r = c(seq(from=0,to=radius,length.out = N),
					rep(radius,N),
					seq(from=radius,to=0,length.out = N))
		color = col[((i-1) %% length(col)) + 1]
		polygon(getPolar(a,r,center),border = border,col = color)
		if(!is.null(labels)){
			a = rep(angle+spiral.shift+x[i]/2,2)
			r = c(radius,radius + tck.length)
			tck = getPolar(a,r,center)
			lines(tck,col=color)
			t = getPolar(a[1],r[2]+tck.length/3,center)
			text(t$x,t$y,labels[i],adj=c(0.5-cos(a[1])/2,0.5-sin(a[1])/2))
		}
		angle = angle + x[i]
		
	}
}

plotPanelLetter = function(l){
	x=grconvertX(0,from='nfc',to='user')
	y=grconvertY(1,from='nfc',to='user')
	text(x=x,y=y,labels=l,adj=c(0,1.1),font=2,cex=1.2,xpd=TRUE)
}

plotMouseTrees = function(trees,main,l){
	max.len = max(unlist(lapply(trees,node.depth.edgelength)))*1.1
	max.root.len = max(sapply(trees,function(t){t$edge.length[t$edge[,2] == which(t$tip.label=='base')]}))
	
	plotTree(trees[[1]],setNames(params$tissue.col,substr(names(params$tissue.col),1,1)),T,
					 xmax=max.len,root.max=max.root.len,root.name='base',xaxs=T,lwd=2,
					 main=paste(names(trees)[1],'days'),direction='upwards',srt=-90,adj=0.5)
	legend('topleft',bty='n',legend = main)
	mtext('1 - rho',2,1.1,cex = 0.8)
	#title(ylab = '1 - rho')
	plotPanelLetter(l)
	for(i in c(3,4,9,12,14))
		plotTree(trees[[i]],setNames(params$tissue.col,substr(names(params$tissue.col),1,1)),T,
						 xmax=max.len,root.max=max.root.len,root.name='base',xaxs=F,lwd=2,
						 main=paste(names(trees)[i],'days'),direction='upwards',srt=-90,adj=0.5)
	
}



plotExonCount = function(l){
	aa = aanns[rownames(clust),]
	t = table(paste0(ifelse(aa$cod!='n','c','n'),ifelse(rownames(aa) %in% orth.seg.ad.ids,'o','n')),clust$species)[c('co','cn','no','nn'),rownames(species)]
	barplot(t,den=c(-1,20),col=rep(c('#FF9999','#9999FF'),each=2),las=3,ylab='# of exons',main='Tissue/age-variable exons',ylim=c(0,10000))
	legend('topright',fill=c('#FF9999','#9999FF','black','black'),den=c(-1,-1,-1,20),legend=c('coding','non-coding','orth.','non-orth.'),bty='n')
	plotPanelLetter(l)
}

plotOrthExonStat = function(l){
	s = apply(orth.seg.ad.ids,2,'%in%',rownames(aanns))
	s = s[apply(s,1,sum)>0,]
	tab = sapply(1:7,function(i){table(apply(s[s[,i],],1,sum))})[7:1,]
	rownames(tab) = paste(rownames(tab),'species')
	colnames(tab) = rownames(species)
	barplot(tab,las=3,ylab='# of orth. exons',legend.text = T,xlim=c(0,18),args.legend = list(title='alternative in:'))
	plotPanelLetter(l)
}

plotMDS = function(d,l,color.by.species=FALSE,.meta=NULL,cex=1,...){
	if(!is.null(.meta))
		meta = .meta
	meta = meta[rownames(d),]
	sp.age.range = sapply(split(meta$age.use,meta$species),range)
	norm.ages = (meta$age.use - sp.age.range[1,meta$species])/(sp.age.range[2,meta$species] - sp.age.range[1,meta$species])
	col = params$tissue.col[meta[rownames(d),'tissue']]
	if(color.by.species)
		col = params$species.col[meta[rownames(d),'species']]
	#col = rgb(t(col2rgb(col)),maxColorValue = 255,alpha=norm.ages*230 + 25)
	plot(d,col=col,pch=params$species.pch[meta[rownames(d),'species']],xlab='Dimension 1',ylab='Dimension 2',cex=cex*(0.3+norm.ages*0.7),...)
	plotPanelLetter(l)
}




plotClusters = function(l){
	c = clust[clust$species=='mouse',]
	m = meta.tsm[meta.tsm$species=='mouse',]
	xx = seq(min(m$age.use),max(m$age.use),length.out = 100)
	for(i in 1:28){
		p = psi.tsm$mouse[rownames(c)[c$orth.cl==i],rownames(m)]
		plot(1,t='n',xaxt='n',yaxt='n',xlim=range(m$age.use),ylim=c(0,1))
		for(t in names(params$tissue.col)){
			y = apply(p[,m$tissue==t],2,mean,na.rm=T)
			x = m$age.use[m$tissue==t]
			x = x[!is.na(y)]
			y = y[!is.na(y)]
			lines(predict(smooth.spline(x,y,df=5),xx),col=params$tissue.col[t],lwd=3)
		}
		text(min(xx),1,paste0('pattern ',i,' (',nrow(p),')'),adj = c(0,1))
		if(i %in% c(15,16,24,28)){
			axis(4,c(0.2,0.4,0.6,0.8))
			mtext('PSI',4,1.1,cex=0.7)
		}
		if(i > 20){
			d = c(12,14,25,40,60)
			axis(1,log(d),d)
			mtext('age (Days)',1,1.1,cex=0.7)
		}
		#if(i == 16) {plot.new();plot.new()}
		#if(i == 24) {plot.new()}
		if(i == 1) { 
			x=grconvertX(0,from='device',to='user')
			y=grconvertY(1,from='nfc',to='user')
			text(x=x,y=y,labels=l,adj=c(0,0),font=2,cex=1.2,xpd=NA)
		}
	}
}

plotPropInSameClust = function(l,prop.of.orth=FALSE){
	cls = list(tu=1,td=2,bhu=c(3,5,7,11),bhd=c(4,6,8,12),lu=9,ld=10,all=1:28)
	r = NULL
	for(sp in c('r','b','h','o','c')){
		o = orth.segs[[paste0(sp,'m')]]
		o = o[o[,1] %in% rownames(clust) & o[,2] %in% rownames(clust),]
		o = cbind(clust[o[,1],'orth.cl'],clust[o[,2],'orth.cl'])
		r. = setNames(rep(0,length(cls)),names(cls))
		for(n in names(cls)){
			if(prop.of.orth)
				r.[n] = sum(o[,2] %in% cls[[n]])/sum(clust$species=='mouse' & clust$orth.cl %in% cls[[n]])
			else
				r.[n] = sum(o[,2] %in% cls[[n]] & o[,1]==o[,2])/sum(o[,2] %in% cls[[n]])
			#r.[n] = sum(o[,2] %in% cls[[n]] & o[,1] %in% cls[[n]])/sum(o[,2] %in% cls[[n]])
		}
		r = rbind(r,r.)
	}
	r = r*100
	col=c('black','black','red','red','green','green','gray')
	lty=c(1,2,1,2,1,2,1)
	
	plot(1,t='n',xlim=c(1,nrow(r)),ylim=range(0,r),xaxt='n',ylab=ifelse(prop.of.orth,'% of orthologous','% of conserved'),main='Pattern conservaion',xlab='')
	axis(1,1:nrow(r),c('rat','rabbit','human','opossum','chicken'),las=2)
	for(i in 1:ncol(r))
		lines(1:nrow(r),r[,i],col=col[i],lty=lty[i],lwd=3)
	if(!prop.of.orth)
		legend('bottomleft',col=c('gray','red','green','black','black','black'),lty=c(1,1,1,1,1,2),legend=c('all','brain and heart','liver and kidney','testis','inclusion','exclusion'),bty='n')
	plotPanelLetter(l)
}



plotAltCntProp = function(l,col){
	z = mouse2rat.alt.cons$conserved/mouse2rat.alt.cons$total*100
	x = log(as.numeric(colnames(z)))
	plot(1,t='n',xlab='Age (days)',ylab='% of conserved in rat',main='Mouse to rat conservation',xlim=range(x),ylim=c(0,max(z)),xaxt='n')
	for(i in 1:nrow(z))
		lines(x,z[i,],col=col[rownames(z)[i]],lwd=3)
	days = c(12,20,40,80)
	axis(1,log(days),days)
	plotPanelLetter(l)
}


plotAltCnt = function(l){
	col = unique(meta[,c('tissue','col')])
	names = c(setNames(col$tissue,substr(col$tissue,1,1)),const='const.',bh='brain+heart','2ts'='2 tissues','3ts'='3 tissues','4ts'='4 tissues','5ts'='5 tissues','6ts'='all tissues')
	col = setNames(col$col,substr(col$tissue,1,1))
	col = c(col,const='gray',bh='orange','2ts'='#0000FFFF','3ts'='#0000FFAA','4ts'='#0000FF88','5ts'='#0000FF66','6ts'='#0000FF44')
	z = mouse.alt.cnt
	rownames(z) = names[rownames(z)]
	areaplot(z,x = log(as.numeric(colnames(z))),xaxt='n',col = col[rownames(mouse.alt.cnt)],xlab='Age (days)',ylab='# of segments',main='Dynamic of AS on age')
	days = c(12,20,40,80)
	axis(1,log(days),days)
	plotPanelLetter(l)
	x=grconvertX(1,from='npc',to='user')
	y=grconvertY(1,from='npc',to='user')
	text(x=x,y=y,labels='Alt. in:',adj=c(0,0),font=2,cex=1.2,xpd=TRUE)
	invisible(col)
}

plotMeanAltProfile = function(sps,main,sp.sp,l,ylim){
	sp.sp$spec[sp.sp$spec=='oc'] = 'hqmrb'
	f = good.orth & sp.sp$spec %in% sps
	f = sapply(1:7,function(s){
		f & grepl(species$short[s],sp.sp$spec)
	})
	ages = sort(unique(meta.tsm.orth$mouse.days))
	tissues = c('brain','cerebellum','heart','liver','kidney','testis','ovary')
	r = matrix(NA,nrow=length(tissues),ncol=length(ages))
	for(t in 1:nrow(r))
		for(a in 1:ncol(r)){
			tmp = c()
			for(s in 1:7){
				cols = meta.tsm.orth$tissue==tissues[t] & meta.tsm.orth$mouse.days==ages[a] & meta.tsm.orth$species==rownames(species)[s]
				tmp = c(tmp,psi.tsm.orth[f[,s],cols])
			}
			r[t,a] = mean(tmp,na.rm=T)
		}
	
	lens = sapply(orth.seg.ad,function(x){x$seg$length})[f]
	
	lens = table(lens %% 3)
	lens = lens/sum(lens)
	yr = ylim
	plot(1,t='n',xlab='Age (mouse days)',ylab='mean PSI',main=paste0(main,' (',sum(sum(f)),')'),log='x',xlim=range(ages,115),ylim=yr)
	for(t in 1:nrow(r))
		lines(ages,r[t,],col=params$tissue.col[tissues[t]],lwd=3)
	
	lenp = round(lens*100,0)
	lens=yr[1] + (yr[2]-yr[1])*cumsum(lens)*0.9
	rect(max(ages)+5,yr[1]  ,100,lens[1],col='#555555')
	rect(max(ages)+5,lens[1],100,lens[2],col='#888888')
	rect(max(ages)+5,lens[2],100,lens[3],col='#CCCCCC')
	
	text(105,c((lens[1]+yr[1])/2,(lens[2:3]+lens[1:2])/2),paste0(0:2,'\n',lenp,'%'),adj=c(0,0.5))
	text(max(ages)+5,max(lens),'Exon\nphase:',adj=c(0,-0.5))
	plotPanelLetter(l)
	invisible(r)
}


plotMeanBornProfile = function(sps,main,l){
	ids = born.seg.ids[sp.birth %in% sps,]
	ages = sort(unique(meta.tsm.orth$mouse.days))
	tissues = c('brain','cerebellum','heart','liver','kidney','testis','ovary')
	r = matrix(NA,nrow=length(tissues),ncol=length(ages))
	for(t in 1:nrow(r))
		for(a in 1:ncol(r)){
			tmp = c()
			for(s in 1:ncol(ids)){
				i = ids[,s]
				i = i[!is.na(i)]
				if(length(i)>0){
					cols = rownames(meta)[meta$species==rownames(species)[s] & meta$tissue==tissues[t] & meta$mouse.days==ages[a] & !is.na(meta$mouse.days)]
					tmp = c(tmp,born.exn.sajr[[s]]$ir[i,cols])
				}
			}
			r[t,a] = mean(tmp,na.rm=T)
		}
	lens = sapply(1:ncol(ids),function(s){
		i = ids[!is.na(ids[,s]),s]
		born.exn.sajr[[s]]$seg[i,'length']
	})
	lens = table(unlist(lens) %% 3)
	lens = lens/sum(lens)
	yr = c(0.08,0.5)
	plot(1,t='n',xlab='Age (mouse days)',ylab='mean PSI',main=paste0(main,' (',sum(!is.na(ids)),')'),log='x',xlim=range(ages,115),ylim=yr)
	for(t in 1:nrow(r))
		lines(ages,r[t,],col=params$tissue.col[tissues[t]],lwd=3)
	lenp = round(lens*100,0)
	lens=yr[1] + (yr[2]-yr[1])*cumsum(lens)*0.9
	rect(max(ages)+5,yr[1]  ,100,lens[1],col='#555555')
	rect(max(ages)+5,lens[1],100,lens[2],col='#888888')
	rect(max(ages)+5,lens[2],100,lens[3],col='#CCCCCC')
	
	text(105,c((lens[1]+yr[1])/2,(lens[2:3]+lens[1:2])/2),paste0(0:2,'\n',lenp,'%'),adj=c(0,0.5))
	text(max(ages)+5,max(lens),'Exon\nphase:',adj=c(0,-0.5))
	plotPanelLetter(l)
	invisible(r)
}

plotNonNormCor = function(l){
	ages = sort(unique(meta.tsm.orth$mouse.days))
	tissues = c('brain','cerebellum','heart','liver','kidney','testis','ovary')
	tab = getCorOnAge(sta.cor,'mouse','rabbit',tissues,ages)
	plot(1,t='n',xlim=range(ages),ylim=range(tab,na.rm=T),log='x',ylab='Mouse-rabbit correlation',xlab='Age (mouse days)',main='Non-normolized PSI')
	for(t in tissues){
		f = !is.na(tab[t,])
		lines(ages[f],tab[t,f],t='b',col=params$tissue.col[t],lwd=3)
	}
	plotPanelLetter(l)
}

plotNormCor = function(l){
	ages = sort(unique(meta.tsm.orth$mouse.days))
	tissues = c('brain','cerebellum','heart','liver','kidney','testis','ovary')
	tab = getCorOnAge(sta.cor.n,'mouse','rabbit',tissues,ages)
	plot(1,t='n',xlim=range(ages),ylim=range(tab,na.rm=T),log='x',ylab='Mouse-rabbit correlation',xlab='Age (mouse days)',main='Normolized PSI')
	for(t in tissues){
		f = !is.na(tab[t,])
		lines(ages[f],tab[t,f],t='b',col=params$tissue.col[t],lwd=3)
	}
	plotPanelLetter(l)
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

plotSpeciesTree = function(l){
	lins = c(species$short,'hq','mr','mrb','oc')
	patt = table(sp.sp.patt.disp.cor$spec[good.orth])[lins]
	mean = table(sp.sp.mean.disp.dst$spec[good.orth & (sp.sp.patt.disp.cor$spec %in% c('n','NA','s5'))])[lins]
	names(patt) = lins
	patt[is.na(patt)] = 0
	birt = c(species$short,'hq','mr','mrb','hqmrb')
	lost = sapply(birt,function(x){paste(species$short[!(species$short %in% strsplit(x,'')[[1]])],collapse = '')})
	birt = table(sp.birth)[birt]
	lost = table(sp.birth)[lost]
	
	cnts = rbind(mean,patt,birt,lost)
	colnames(cnts) = names(birt)
	text.col=c('red','orange','blue','green')
	
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
	
	y = (0:4/4)[-5]
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
	
	text(6,2-y,c('mean','pattern','birth','lost'),col = text.col,adj = c(0,1))
	b = length(sp.birth) - sum(birt+lost)
	m = sum(!(sp.sp.mean.disp.dst$spec %in% c(lins,'n','NA','s5')) & good.orth & (sp.sp.patt.disp.cor$spec %in% c(lins,'n','NA','s5')))
	a = sum(good.orth & !(sp.sp.patt.disp.cor$spec %in% c(lins,'n','NA','s5')))
	text(8,6-y,c('Complex:',paste(c('mean:','pattern:','born/lost:'),c(m,a,b))),col = c('black','red','orange','blue'),adj = c(0,1))
	plotPanelLetter(l)
}

plotClassOfOrthExons = function(l){
	a = table(nchar(orth.nalt[good.orth]))
	names(a) = paste0(names(a),' (',a,')')
	b = table(nchar(sp.birth))
	names(b) = paste0(names(b),' (',b,')')
	ac = 10:4/11
	bc = 9:4/10
	col=c(rgb(1,ac,ac),rgb(bc,bc,1))
	psichart(c(a,b),spiral.shift = 1,col=col,radius = 0.25,init.angle = -1,center = c(0.6,0.6),main='Orthologous exons')
	legend('bottomleft',bty='n',fill=col[c(7,13)],legend=c('exon is alternative','exon is observed'),title = 'Number of species where:')
	cnt = c(sum(good.all.orth)-sum(a),sum(a),sum(b))
	psichart(cnt,spiral.shift = 1,init.angle = 1.1,labels = paste(c('const.','alt.','born/lost'), '(',cnt,')'),add=T,col=c('#88FF88',col[c(7,13)]),center = c(0.2,0.85),radius = 0.13)
	plotPanelLetter(l)
}
