bootstrapDevAS = function(m,d,sps,replace=TRUE,N=100,ss=c('ad','aa','dd','da'),only.samples=FALSE){
	m = m[!is.na(m$mouse.stage) & m$species %in% sps, ]
	d = d[,colnames(d$ir) %in% rownames(m)]
	stat = apply(table(factor(m$species,levels = sps),m$mouse.stage,m$tissue),2:3,min)
	m = m[colnames(d$ir),]
	
	res = llply(1:N,function(i){
		#cat(i,'   ')
		# sample samples
		sids = c()
		for(t in colnames(stat)){
			for(s in rownames(stat)){
				sids = c(sids,sample(rownames(m)[m$mouse.stage==s & m$tissue==t],stat[s,t],replace = replace))
			}
		}
		m = m[sids,]
		if(only.samples)
			return(m)
		dt = d[,sids]
		colnames(dt$ir) = colnames(dt$i) = colnames(dt$e) = rownames(m) # to fix names of duplicated due to bootstrap samples
		# filter
		f = rep(FALSE,length(dt))
		na = is.na(dt$ir)
		for(tis in unique(m$tissue)){
			cinx = m$tissue==tis
			f = f | (apply(!na[,cinx],1,mean) > 0.6 & apply(dt$ir[,cinx],1,function(x){x=x[!is.na(x)];sum(x>0.1 & x < 0.9)>3}))
		}
		dt = dt[f,]
		na = is.na(dt$ir)
		# test devAS
		# I will use species ages.. but maybe mouse-scale is more appropriate
		tissues = unique(m$tissue)
		qv = sapply(tissues,function(tissue){testASAge(dt,m,tissue,min.cov.sams=0.6,.parallel=FALSE)})
		patt = matrix('-',ncol=length(tissues),nrow=length(dt),dimnames = list(rownames(dt$ir),tissues))
		for(t in colnames(qv)){
			# 'c' mean that exon passed coverage thershold
			cinx = m$tissue==t
			patt[apply(!na[,cinx],1,mean) > 0.6,t] = 'c'
			qv[apply(!na[,cinx],1,mean) <= 0.6,t] = NA # just to be sure. it have to be done in testASAge
			#
			patt[!is.na(qv[,t]),t] = 'n'
			f = !is.na(qv[,t]) & qv[,t] < 0.05
			if(sum(f)==0) next
			a = m$age.use[m$tissue==t]
			x = apply(dt$ir[f,m$tissue==t,drop=F],1,function(y)getDevASPattern1(a,y,4,1000))
			x = t(x)
			x = x[apply(is.na(x),1,sum)==0 & (x[,1]+x[,2])>0 & x[,'max']-x[,'min'] > 0.2,,drop=F]
			if(nrow(x) == 0 ) next
			u = x[,1]/(x[,1]+x[,2])
			tim = x[,3]/x[,1]-x[,4]/x[,2]
			thr = 0.3
			patt[rownames(x),t] = 'du'
			patt[rownames(x)[!is.na(tim) & tim < 0],t] = 'ud'
			patt[rownames(x)[u < thr],t] = 'd'
			patt[rownames(x)[u > (1 - thr)],t] = 'u'
		}
		gc()
		patt
	},.parallel = TRUE)
	res
}


loadIntronCounts = function(files,names,ints=NULL){
	data = list()
	for(i in 1:length(files)){
		cat('\r',i,' ',length(files),'       ')
		x = read.table(files[i],header = T)
		if(!is.null(ints))
			x = x[grepl(paste(ints,collapse = '|'),x[,1]),]
		data[[length(data)+1]] = x
	}
	if(!is.null(ints))
		ints = unique(unlist(lapply(data,'[',1)))
	r = matrix(0,ncol=length(names),nrow=length(ints),dimnames = list(ints,names))
	for(i in 1:length(names)){
		r[data[[i]][,1],i] = data[[i]][,2]
	}
	ints = strsplit(ints,'[:-]')
	r = list(intron=data.frame(chr_id=sapply(ints,'[',1),
														 start =as.numeric(sapply(ints,'[',2)),
														 stop  =as.numeric(sapply(ints,'[',3)),
														 strand=as.numeric(sapply(ints,'[',4))
	),
	cnt=r)
	r$intron$strand[is.na(r$intron$strand)] = -1
	rownames(r$intron) = rownames(r$cnt)
	class(r) = c('sajr','list')
	r
}

getExonNumberByGene = function(gtf){
	a = read.table(gtf,sep='\t')
	a = a[a$V3=='exon',]
	z = lapply(strsplit(a$V9,'; ',T),function(x){
		z = strsplit(x,' ',F)
		setNames(sapply(z,'[',2),sapply(z,'[',1))
	})
	y=t(sapply(z,'[',c('gene_id','transcript_id','exon_number')))
	y = as.data.frame(y)
	y$exon_number = as.numeric(y$exon_number)
	y = do.call(rbind,lapply(split(y,y$transcript_id),function(x)x[order(-x$exon_number)[1],]))
	split(y$exon_number,y$gene_id)
}

getDiamBySpline = function(x,y,df){
	f = !is.na(x) & !is.na(y)
	x = x[f]
	y = y[f]
	if(length(unique(x))<= df) return(NA)
	p = predict(smooth.spline(x,y,df=df))$y
	max(p)-min(p)
}

getAdjStageChange = function(as,ge,m,tissue,s2e,melt=TRUE){
	library(reshape)
	m = m[m$tissue==tissue,]
	dpsi.age = calcdPSIonAge(as,m,stages2use=NULL)
	dexp.age = calcdPSIonAge(ge,m,stages2use=NULL)
	c = nrow(dpsi.age)
	dpsi.age = dpsi.age[rownames(dpsi.age) %in% names(s2e)[sapply(s2e,length)==1],]
	gid = unlist(s2e[rownames(dpsi.age)])
	f = gid %in% rownames(dexp.age)
	cat('Ambig. exon filter. ',c,' -> ',nrow(dpsi.age),' -> ',sum(f),'\n')
	dpsi.age = dpsi.age[f,]
	gid = gid[f]
	dexp.age = dexp.age[gid,colnames(dpsi.age)]
	colnames(dpsi.age) = colnames(dexp.age) = sapply(strsplit(colnames(dpsi.age),' ',fixed = T),'[',3)
	if(melt){
		dpsi.age = melt(dpsi.age)
		dexp.age = melt(dexp.age)
		colnames(dpsi.age) = c('seg.id','stage','dpsi')
		return(cbind(seg.id=dpsi.age[,1],ens.id=dexp.age[,1],dpsi.age[,2:3],de=dexp.age[,3]))
	}
	return(list(as=dpsi.age,ge=dexp.age))
}

equalyseDistrs = function(x,y,bins=10,min.bin.size=10){
	#samples from y with distr as in x
	r = range(x,y,na.rm=T)
	bin.size = (r[2]-r[1])/bins
	xr = floor((x-r[1])/bin.size)
	yr = floor((y-r[1])/bin.size)
	xt = table(factor(xr,levels = 0:(bins-1)))
	yt = table(factor(yr,levels = 0:(bins-1)))
	f = min((yt/xt)[xt>=min.bin.size])
	unlist(lapply(0:(bins-1),function(b){if(yt[b+1]>0){sample(which(yr==b),min(yt[b+1],xt[b+1]*f))}else{NULL}}))
}

getAdjStageCor = function(s,t,dpsi.thr=0.2,cor.meth='p'){
	f=anns[[s]]$sites=='ad' & per.tissue.age.qv[[s]][,t] < 0.05 & abs(age.dpsi[[s]][,t]) > 0.2
	m = psi.tsm[[s]][!is.na(f) & f,]
	m = m[,meta.tsm[colnames(m),'tissue'] == t]# & meta.tsm[colnames(m),'stage'] %in% age.al.i[,s]]
	m = m[,order(meta.tsm[colnames(m),'days'])]
	if(cor.meth=='s')
		m = apply(m,2,rank,na.last = 'keep')
	r=t(sapply(1:(ncol(m)-1),function(i){r=cor.test(m[,i],m[,i+1],u='p',m='p');c(r$estimate,r$conf.int)}))
	rownames(r) = meta.tsm[colnames(m),'stage'][1:nrow(r)]
	r
}

getSplineDevdPSI = function(s,t,dpsi.thr=0.2,cor.meth='p',df=4,N=100,shuffle=F){
	f=anns[[s]]$sites=='ad' & per.tissue.age.qv[[s]][,t] < 0.05 & abs(age.dpsi[[s]][,t]) > dpsi.thr
	m = psi.tsm[[s]][!is.na(f) & f,]
	m = m[,meta.tsm[colnames(m),'tissue'] == t]# & meta.tsm[colnames(m),'stage'] %in% age.al.i[,s]]
	m = m[,order(meta.tsm[colnames(m),'days'])]
	x = meta.tsm[colnames(m),'age.rank']
	xp = seq(min(x),max(x),length.out = N)
	r = apply(m,1,function(y){
		f= !is.na(y)
		if(shuffle)
			x = sample(x)
		m = smooth.spline(x[f],y[f],df=df)
		predict(m,xp)$y
	})
	ru = rd = r[-1,] - r[1:(N-1),]
	ru[ru<0] = 0
	rd[rd>0] = 0
	mu = apply(abs(ru),1,mean)
	sdu = apply(abs(ru),1,sd)/sqrt(ncol(ru))
	
	md = apply(abs(rd),1,mean)
	sdd = apply(abs(rd),1,sd)/sqrt(ncol(ru))
	
	c = (N-1)/(max(xp)-min(xp))
	list(x=xp,up= cbind(mu,mu-2*sdu,mu+2*sdu)*c,dw= cbind(md,md-2*sdd,md+2*sdd)*c)
}

calcdPSIonAge = function(p,m,stages2use=NULL){
	if(!is.null(stages2use))
		m = m[m$stage %in% stages2use,]
	cmn = intersect(rownames(m),colnames(p))
	m = m[cmn,]
	p = p[,cmn]
	ds = unique(m[,c('days','stage')])
	ds = ds[order(ds$days),]
	do.call(cbind,lapply(unique(m$tissue),function(t){
		r = matrix(NA,ncol=nrow(ds)-1,nrow=nrow(p),dimnames = list(rownames(p),paste(m$species[1],t,ds$stage[-nrow(ds)])))
		for(i in 1:(nrow(ds)-1)){
			i1 = which(m$tissue==t & m$days == ds$days[i])
			i2 = which(m$tissue==t & m$days == ds$days[i + 1])
			if(length(i1) == 1 & length(i2) == 1)
				r[,i] = p[,i2] - p[,i1]
		}
		r
	}))
}



plotdPSIOnAge = function(sp,dpsi,stages2use=NULL){
	dpsi.age = calcdPSIonAge(psi.tsm[[sp]],meta.tsm,stages2use=stages2use)
	for(t in unique(meta$tissue)){
		f = grep(t,colnames(dpsi.age))
		up = apply(dpsi.age[anns[[sp]]$sites=='ad' & per.tissue.age.qv[[sp]][,t]<0.05,f]>  dpsi,2,sum,na.rm=T)
		dw = apply(dpsi.age[anns[[sp]]$sites=='ad' & per.tissue.age.qv[[sp]][,t]<0.05,f]< -dpsi,2,sum,na.rm=T)
		
		# up =  apply(dpsi.age[anns$mouse$sites=='ad' & per.tissue.age.qv$mouse[,t]<0.05,f],2,function(x){mean(x[!is.na(x) & x > 0])})
		# dw = -apply(dpsi.age[anns$mouse$sites=='ad' & per.tissue.age.qv$mouse[,t]<0.05,f],2,function(x){mean(x[!is.na(x) & x < 0])})
		
		
		ylim=range(0,up,dw,na.rm=T)
		#ylim=c(0,0.05)
		plotTissueAgeProile(up,meta.tsm,df = 0,age.axis = 'rank',tissues = t,lty=1,pch=19,cex=1,ylim=ylim,xlim=c(1,max(meta.tsm$age.rank[meta.tsm$species==sp])),main=paste(sp,t),ylab='# of exons')
		plotTissueAgeProile(dw,meta.tsm,df = 0,age.axis = 'rank',tissues = t,lty=2,pch=1,add=T,cex=1)
		
		m = meta.tsm[meta.tsm$species==sp,]
		b=m$age.rank[order(abs(species[sp,'gestation']-m$days))[1]]
		abline(v=b,lty=3)
	}
}

getBestCor = function(c,gr){
	r=sapply(1:ncol(c),function(i){
		max(c[i,gr == gr[i] & (1:ncol(c)) != i])
	})
	names(r) = rownames(c)
	r
}

ftHex = function(h,f){
	f = factor(f,levels = c(FALSE,TRUE))
	apply(h,2,function(x){r=fisher.test(factor(x,levels = c(FALSE,TRUE)),f,a='g');c(pv=r$p.value,or=r$estimate)})
}
ftHexAge = function(h,a,dir){
	pv = or = array(NA,dim=c(ncol(h[[1]]),ncol(a[[1]]),length(h)),dimnames=list(colnames(h[[1]]),colnames(a[[1]]),names(h)))
	for(t in colnames(a[[1]])){
		cat('\n',t,' ')
		for(s in names(h)){
			cat(s,' ')
			if(t %in% colnames(a[[s]])){
				cmn = intersect(rownames(h[[s]]),rownames(a[[s]])[a[[s]][,t]!='-'])
				r = ftHex(h[[s]][cmn,],a[[s]][cmn,t]==dir)
				pv[,t,s] = r[1,]
				or[,t,s] = r[2,]
			}
		}
	}
	list(pv=pv,or=or)
}

downSampledevAS = function(as,dirs,sites='ad',per.tissue=FALSE){
	r=lapply(rownames(species),function(s){
		x = as[[s]][anns[[s]]$sites==sites,]
		cnt = unlist(lapply(dirs,function(dir)apply(x==dir,2,sum)))
		mcnt = min(cnt[cnt>0])
		for(i in 1:ncol(x)){
			if(per.tissue)
				mcnt = min(sapply(dirs,function(dir)sum(x[,i]==dir)))
			for(dir in dirs){
				sids = which(x[,i]==dir)
				if(length(sids)>mcnt){
					torem=sample(sids,length(sids)-mcnt)
					x[torem,i] = 'n'
				}
			}
		}
		x
	})
	names(r) = rownames(species)
	r
}

plotExpAndPsiForNewExons = function(s,sps,born.psi,born.ids,exp,center=F,scale=F,max.na.prop = 1,ylab.fun.name=''){
	m=calcExpAndPsiForNewExons(s,sps,born.psi,born.ids,exp,center,scale,max.na.prop)
	for(i in 1:length(m$psis)){
		plotTissueAgeProile(m$psis[[i]],meta,main=paste0(sps[i],' (',m$stat[i,2],')'),ylab=paste0(ylab.fun.name,ifelse(ylab.fun.name=='','','('),'PSI',ifelse(ylab.fun.name=='','',')')))
		plotTissueAgeProile(m$exps[[i]],meta,main=paste0(sps[i],' (',m$stat[i,3],')'),ylab=paste0(ylab.fun.name,ifelse(ylab.fun.name=='','','('),'RPKM',ifelse(ylab.fun.name=='','',')')))
	}
	invisible(m)
}


calcExpAndPsiForNewExons = function(s,sps,born.psi,born.ids,exp,center=F,scale=F,max.na.prop = 1){
	psis = exps = list()
	seg.stat = NULL
	for(ss in sps){
		segs = born.ids[born.ids$species==ss,s]
		#PSI
		stat = c(all.segs = length(segs))
		x = born.psi[[s]]$ir[segs,]
		x = t(scale(t(x),center = center,scale =scale))
		x = x[apply(is.na(x),1,mean)<=max.na.prop,]
		stat = c(stat,filtered.segs = nrow(x))
		psis[[ss]] = apply(x,2,mean,na.rm=TRUE)
		#Exp
		x = exp[[s]]$rpkm[intersect(rownames(exp[[s]]$rpkm),born.psi[[s]]$seg[rownames(x),'gene_id']),]
		x = t(scale(t(x),center = center,scale =scale))
		stat = c(stat,genes = nrow(x))
		exps[[ss]] = apply(x,2,mean,na.rm=TRUE)
		seg.stat = rbind(seg.stat,stat)
	}
	rownames(seg.stat) = sps
	list(psis=psis,exps=exps,stat=seg.stat)
}

getASChangeCons = function(age.segs,dir,sps){
	res = lapply(colnames(age.segs[[1]]),function(t){
		r = sapply(sps,function(s){
			f = age.segs[[sps[1]]][,t] != '-' & age.segs[[s]][,t] != '-'
			c(both.exp = sum(f),age1=sum(age.segs[[sps[1]]][,t] == dir),age2=sum(age.segs[[s]][,t] == dir),age.both=sum(age.segs[[sps[1]]][,t] == dir & age.segs[[s]][,t] == dir))
		})
		apply(r,2,function(x){
			bt = binom.test(x[4],x[2])
			c(x,freq=x[4]/x[2],freq.conf = bt$conf.int)
		})
	})
	names(res) = colnames(age.segs[[1]])
	res
}

caclOrthPsiCor = function(psi,t,ages,f=TRUE,cor.method='pe'){
	rs = names(ages)[1]
	r = sapply(names(ages),function(s){
		r = cor.test(psi[[rs]][f,paste(rs,t,ages[rs])],psi[[s]][f,paste(s,t,ages[s])],method = cor.method)
		if(!is.null(r$conf.int))
			r = c(r$estimate,conf=r$conf.int)
		else
			r = c(r$estimate,conf1=r$estimate,conf2=r$estimate)
	})
	r
}

plotAgeChangeCons = function(d,inxs=5:7,...){
	x=1:ncol(d[[1]])
	plot(1,t='n',xlim=range(x),xaxt='n',...)
	grid(nx = 0,ny=NULL)
	for(t in names(d))
		plotArea(x,t(d[[t]][inxs,]),col=cols[t],lwd=3)
	axis(1,x,colnames(d[[1]]),las=2)
}


getAgeASchanges = function(psi,m,psi.thr,stages,sp,get.dPSI=FALSE){
	psi = psi[[sp]]
	stages = stages[[sp]]
	m = m[colnames(psi),]
	sapply(unique(m$tissue),function(t){
		dpsi = psi[,paste(sp,t,stages[t,2])] - psi[,paste(sp,t,stages[t,1])]
		if(get.dPSI)
			r = dpsi
		else{
			r = ifelse(dpsi>psi.thr,'u',ifelse(dpsi< -psi.thr,'d','n'))
			r[is.na(r)] = '-'
		}
		r
	})
}


firstToupper = function(t){
	paste0(toupper(substr(t,1,1)),substr(t,2,nchar(t)))
}

number2bin = function(v,n){
	o = order(v)
	j = 1
	for(i in 1:length(v)){
		if(j < i/length(v)*n)
			j = j + 1
		v[o[i]] = j
	}
	v
}

getGeneBiotype = function(gtf){
	d = read.table(gtf,sep='\t')[,9]
	d = gsub('"','',d)
	d = strsplit(d,'; ?',FALSE,perl = TRUE)
	d = sapply(d,function(x){x=do.call(rbind,lapply(strsplit(x,' ',TRUE),'[',1:2));rownames(x) = x[,1];x[c('gene_id','gene_biotype'),2]})
	d = unique(as.data.frame(t(d)))
	setNames(d$gene_biotype,d$gene_id)
}

calcGEtsm = function(d,meta){
	lapply(d,function(x){
		x = x[,colnames(x$rpkm) %in% rownames(meta)]
		m = meta[colnames(x$rpkm),]
		calcMeanCols(x$rpkm,paste(m$species,m$tissue,m$stage))
	})
}

addCDSPos2Ann = function(t){
	ids = rownames(t)
	t$cds.pos='n'
	t$seg_id = rownames(t)
	t = split(t,t$gene_id)
	t = lapply(t,function(x){
		if(sum(x$cod!='n')>0){
			cds.min = min(x$start[x$cod!='n'])
			cds.max = max(x$stop[x$cod!='n'])
			x$cds.pos = 'cds'
			x$cds.pos[(x$stop<cds.min & x$strand == 1) | (x$start>cds.max & x$strand ==-1)] = '5utr'
			x$cds.pos[(x$stop<cds.min & x$strand ==-1) | (x$start>cds.max & x$strand == 1)] = '3utr'
		}
		x
	})
	gc()
	t = as.data.frame(rbindlist(t))
	rownames(t) = t$seg_id
	t$seg_id = NULL
	t[ids,]
}


sampleSameDistr = function(sel,all,bin.count,min.bin.cnt=10){
	all = all[!is.na(all)]
	sel = sel[!is.na(sel)]
	all = all[setdiff(names(all),names(sel))]
	mn = min(all,sel)
	bin.size = (max(all,sel)-    mn)/bin.count
	
	all.bins = floor((all-mn)/bin.size)+1
	sel.bins = floor((sel-mn)/bin.size)+1
	
	all.tab = table(all.bins)
	sel.tab = table(sel.bins)
	all.tab = all.tab[names(sel.tab)]
	all.tab[is.na(all.tab)] = 0
	ratio = min(all.tab[all.tab>=min.bin.cnt]/sel.tab[all.tab>=min.bin.cnt])
	good = c()
	sel.inx = names(all) %in% names(sel)
	for(i in 1:bin.count){
		cnt = ratio*sum(sel.bins==i)
		all.obs = names(all)[all.bins==i]
		good = c(good,sample(all.obs,size=min(cnt,length(all.obs))))
	}
	good
}


plotPTBPTargets = function(targ,psi,m,center,scale){
	psi = psi[,colnames(psi) %in% rownames(m)]
	getBrain = function(psi,m,center,scale){
		m = m[m$tissue == 'brain' & rownames(m) %in% colnames(psi),]
		m = m[order(m$days),]
		days = sapply(split(m$days,m$stage),function(x)round(mean(x)))
		m$days = days[m$stage]
		psi = psi[,rownames(m),drop=FALSE]
		psi = calcMeanCols(psi,m$days)
		psi = psi[,order(as.numeric(colnames(psi))),drop=FALSE]
		if(center)
			psi = sweep(psi,1,apply(psi,1,mean,na.rm=TRUE),'-')
		if(scale)
			psi = sweep(psi,1,apply(psi,1,sd,na.rm=TRUE),'/')
		psi
	}
	for(i in 1:4){
		pd = getBrain(psi[targ[[(i-1)*2 + 1]],colnames(psi) %in% rownames(m),drop=F],m,center=center,scale=scale)
		pi = getBrain(psi[targ[[(i-1)*2 + 2]],colnames(psi) %in% rownames(m),drop=F],m,center=center,scale=scale)
		pv = sapply(1:ncol(pd),function(i)wilcox.test(pd[,i],pi[,i])$p.value)
		x = 1:ncol(pd)
		pda = cbind(apply(pd,2,mean,na.rm=T),apply(pd,2,function(x){x=x[!is.na(x)];sd(x)/sqrt(length(x))}))
		pia = cbind(apply(pi,2,mean,na.rm=T),apply(pi,2,function(x){x=x[!is.na(x)];sd(x)/sqrt(length(x))}))
		
		pda = cbind(pda[,1],pda[,1]-2*pda[,2],pda[,1]+2*pda[,2])
		pia = cbind(pia[,1],pia[,1]-2*pia[,2],pia[,1]+2*pia[,2])
		ylab = 'PSI'
		if(center)
			ylab = paste('centered',ylab)
		if(scale)
			ylab = paste('scaled',ylab)
		ylim=c(0,1)
		if(center || scale)
			ylim=range(pda,pia,na.rm=TRUE)
		plotArea(x,pda,'red',new=TRUE,ylim=ylim,xlab='Age',main=names(ptbp1.targ)[i*2],xaxt='n',ylab=ylab)
		plotArea(x,pia,'blue')
		
		p = par('usr')
		f = pv<0.1
		if(sum(f)>0){
			y = p[4] - (p[4]-p[3])*0.05
			points(x[f],rep(y,sum(f)),pch='*',cex=2)
		}
		f = pv<0.05
		if(sum(f)>0){
			y = p[4] - (p[4]-p[3])*0.1
			points(x[f],rep(y,sum(f)),pch='*',cex=2)
		}
		f = pv<0.005
		if(sum(f)>0){
			y = p[4] - (p[4]-p[3])*0.15
			points(x[f],rep(y,sum(f)),pch='*',cex=2)
		}
		legend('bottomright',col=c('red','blue'),lwd=2,legend=paste0(c('dependent','independent'),' (',c(length(targ[[(i-1)*2 + 1]]),length(targ[[(i-1)*2 + 2]])),')'))
		legend('bottomleft',legend=c('*: pv<0.1','**: pv<0.05','***: pv<0.05'))
		lab = as.numeric(colnames(pd)) - species[m[colnames(psi)[1],'species'],'gestation']
		alab = abs(lab)
		lab = paste0(ifelse(lab<0,'-',''),ifelse(alab<31,paste0(alab,'d'),ifelse(alab<366,paste0(round(alab/30),'m'),paste0(round(alab/365),'y'))))
		axis(1,x,lab)
	}
}


plotPhastGsnapOnID = function(sids,col){
	r = list()
	layout(rbind(c(1,1,1,1,1),2:6,7:11),heights=c(2,1,1))
	par(tck=-0.02,mgp=c(1.1,0.2,0),mar=c(5,2,1.5,0),oma=c(0,0,2,1))
	r$phast.prof = plotPhastConsProfilesByID(sids,phastcons,col=col)
	legend('topright',fill=col,legend=paste0(names(sids),' (',sapply(sids,length),')'))
	r$phast.distr = plotPhastConsDistrByID(sids,col=col)
	r$gnomad = plotGnomadFreqsByID(sids,col=col)
	invisible(r)
}

plotPhastConsDistrByID = function(sids,col,...){
	r = vector(mode = 'list',5)
	names(r) = c('u100.50','u50.0','exon','d0.50','d50.100')
	for(n in names(sids)){
		p = phastcons[sids[[n]]]
		r[['u100.50']][[n]] = sapply(p,function(x){mean(x[101:150])})
		r[['u50.0']][[n]] = sapply(p,function(x){mean(x[151:200])})
		r[['exon']][[n]] = sapply(p,function(x){mean(x[201:(length(x)-200)])})
		r[['d0.50']][[n]] = sapply(p,function(x){mean(x[(length(x)-199):(length(x)-150)])})
		r[['d50.100']][[n]] = sapply(p,function(x){mean(x[(length(x)-149):(length(x)-100)])})
	}
	for(n in names(r)){
		boxplot(r[[n]],col=col,nothc=T,ylab='phastcons',ylim=c(0,1),main=n,notch=TRUE,outline=F,las=3,...)
	}
	invisible(r)
}

plotPhastConsProfilesByID = function(sids,pc,cols,int.mar=100,exon.mar=25,ylim=NULL,...){
	profs = lapply(sids,function(x)getPhastconsProf(pc[x],int.mar,200+exon.mar))
	if(is.null(ylim))
		ylim = range(unlist(profs))
	x = 1:((int.mar+exon.mar+1)*2)
	plot(1,t='n',ylim=ylim,xlim=range(x),xaxt='n',ylab='phastcons',...)
	for(i in 1:length(profs))
		plotArea(x,profs[[i]],cols[i])
	axis(1,at=c(int.mar,int.mar+exon.mar*2),c('exon start','exon end'))
	invisible(profs)
}

plotGnomadFreqsByID = function(sids,col){
	r= vector(mode = 'list',5)
	names(r) = c('u100.50','u50.0','exon','d0.50','d50.100')
	for(n in names(sids)){
		g = gnomad[gnomad$seg_id %in% sids[[n]] & gnomad$alt_cnt>1,]
		r[['u100.50']][[n]] = log10(g$freq[g$dist2start>= -100 & g$dist2start< -50])
		r[['u50.0']][[n]] = log10(g$freq[g$dist2start>= -50 & g$dist2start< 0])
		r[['exon']][[n]] = log10(g$freq[g$dist2start>= 0 & g$dist2stop>= 0])
		r[['d0.50']][[n]] = log10(g$freq[g$dist2stop>= -50 & g$dist2stop< 0])
		r[['d50.100']][[n]] = log10(g$freq[g$dist2stop>= -100 & g$dist2stop< -50])
	}
	p = lapply(r,function(rr){
		sapply(rr,function(x){
			m = mean(x)
			s = sd(x)/sqrt(length(x)) 
			c(mean=m,lower = m-2*s,upper=m+2*s)
		})
	})
	ylim=range(unlist(p))
	x = 1:length(sids)
	for(n in names(p)){
		plot(x,p[[n]][1,],col=col,xaxt='n',ylab='log10(SNP freq)',main=n,ylim=ylim,xlab='')
		segments(x,p[[n]][2,],x,p[[n]][3,],col=col)
		axis(1,x,names(sids),las=3)
	}
	invisible(r)
}



plotBinomBar2 = function(s,t,...){
	p = s/t
	sd = p*(1-p)/sqrt(t)
	names(p) = paste0(names(p),'\n(',t,')')
	b = barplot(p,ylim=c(0,max(p+2*sd)),...)
	segments(b,pmin(1,p+2*sd),b,pmax(0,p-2*sd))
}


boxplotWithSgn = function(d,test.fun=wilcox.test,pv.thr=0.05,col.sgn.lines='magenta',from.top=TRUE,sgn.line.frac=0.1,col='white',add.group.sizes=F,...){
	n = names(d)
	if(add.group.sizes){
		n = paste0(n,'\n(',sapply(d,length),')')
	}
	boxplot(d,col=col,names=n,...)
	pv = matrix(1,ncol=3,nrow=length(d)*(length(d)-1)/2)
	n = 1
	for(i in 1:(length(d)-1))
		for(j in (i+1):length(d)){
			if(sum(!is.na(d[[i]])) > 0 & sum(!is.na(d[[j]])) > 0)
				pv[n,] = c(i,j,test.fun(d[[i]],d[[j]])$p.value)
			n = n + 1
		}
	#pv = pv[p.adjust(pv[,3],method = 'BH')<pv.thr,,drop=FALSE]
	pv = pv[pv[,3]<pv.thr,,drop=FALSE]
	if(nrow(pv)>0){
		yr = par('usr')[3:4]*0.95
		r = (yr[2]-yr[1])*sgn.line.frac
		yr = ifelse(rep(from.top,nrow(pv)),seq(yr[2],yr[2]-r,length.out = nrow(pv)),seq(yr[1],yr[1]+r,length.out = nrow(pv)))
		segments(pv[,1],yr,pv[,2],yr,col=col.sgn.lines)
	}
}


imageSpSpCramer = function(d,...){
	cv = caclCramersVPerCols(d,FALSE)
	pv = p.adjust(caclCramersVPerCols(d,TRUE),m='BH')
	imageWithText(cv,text.col=ifelse(pv<0.05,'black','gray'),yaxt='n',xaxt='n',xlab='',ylab='',...)
	axis(1,1:10,paste(rep(colnames(patt.sp.sp),times=2),rep(c('patt.','mean'),each=5)),las=2)
	axis(2,1:10,paste(rep(colnames(patt.sp.sp),times=2),rep(c('patt.','mean'),each=5)),las=2)
	invisible(list(cv=cv,pv=pv))
}

caclCramersVPerCols = function(d,pvalue=FALSE){
	library(lsr)
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

getSpSpPerTissue = function(d,othr,rthr){
	r = matrix(NA,ncol=length(d),nrow=nrow(d[[1]]),dimnames=list(rownames(d[[1]]),names(d)))
	for(t in names(d)){
		r[!is.na(d[[t]]$out),t] = 'none'
		f = !is.na(d[[t]]$out) & d[[t]]$out >= othr & d[[t]]$out/d[[t]]$within > rthr
		r[f,t] = d[[t]]$species[f]
	}
	r
}

plotPhastconsForSpSp = function(patt,mean,phast,f,tissue,s1,s2){
	no.change = apply(cbind(patt.sp.sp,mean.sp.sp)!='none',1,sum,na.rm=T)==0
	no.change[is.na(no.change)] = FALSE
	filts = list(
		exn.no.change = hmo.good &!hmo.anci & no.change & !is.na(patt.sp.sp[,tissue]),
		patt.mouse = f & !is.na(patt[,tissue]) & patt[,tissue] == 'mouse',
		patt.human = f & !is.na(patt[,tissue]) & patt[,tissue] == 'human',
		mean.mouse = f & !is.na(mean[,tissue]) & mean[,tissue] == 'mouse',
		mean.human = f & !is.na(mean[,tissue]) & mean[,tissue] == 'human',
		alt.no.change = hmo.good & hmo.anci & no.change & !is.na(patt.sp.sp[,tissue])
	)
	cols = c('black','red','orange','blue','cyan','green')
	x = 1:((s2-s1+1)*2)
	plot(1,t='n',xlim=range(x),ylim=c(0,1),main='',ylab='mean primate phastcons',xaxt='n',xlab='position (nt)')
	for(i in 1:length(filts)){
		plotArea(x,getPhastconsProf(phast[filts[[i]]],s1,s2),col=cols[i])
	}
	axis(1,c(1,200-s1+1,2*(s2-s1+1)-200+s1,length(x)),c(-(200-s1+1),'acc','don',(200-s1+1)))
	
	up = lapply(filts,function(f){sapply(phast[f],function(p){mean(p[s1:200])})})
	boxplotWithSgn(up,col=cols,notch = T,ylab='mean primate phastcons',main='upstream exon',las=3,col.sgn.lines='magenta')
	ex = lapply(filts,function(f){sapply(phast[f],function(p){mean(p[200:(length(p)-200)])})})
	boxplotWithSgn(ex,col=cols,notch = T,ylab='mean primate phastcons',main='Exon',las=3,from.top=FALSE,col.sgn.lines='magenta')
	up = lapply(filts,function(f){sapply(phast[f],function(p){mean(p[(length(p)-199):(length(p)-s1)])})})
	boxplotWithSgn(up,col=cols,notch = T,ylab='mean primate phastcons',main='downstream exon',las=3,col.sgn.lines='magenta')
	plot.new()
	legend('topleft',fill=cols,legend=paste0(names(filts),' (',sapply(filts,sum),')'))
}

plotPerTissueSpSpcounts = function(d,f,othrs,rthrs,...){
	p1 = sapply(d,function(x){table(factor(x$species[f & !is.na(x$out) & x$out>othrs[1] & x$out/x$within>rthrs[1]],levels=names(hmo.seg.ad)))})
	p2 = sapply(d,function(x){table(factor(x$species[f & !is.na(x$out) & x$out>othrs[2] & x$out/x$within>rthrs[2]],levels=names(hmo.seg.ad)))})
	barplot(p1,beside = T,ylab='# of segments',...)
	barplot(p2,beside = T,add=T,col='black',density = 20)
}

polorize3 = function(d){
	r = data.frame(species=rep(NA,dim(d)[1]))
	r$within = r$out = NA
	rownames(r) = dimnames(d)[[1]]
	s = dimnames(d)[[2]]
	for(i in 1:nrow(r)){
		diag(d[i,,]) = NA
		if(sum(is.na(d[i,,]))>dim(d)[2]) next
		close = which(d[i,,]==min(d[i,,],na.rm=T), arr.ind=T)[1,]
		r$species[i]=s[-close]
		r$within[i] = d[i,close[1],close[2]]
		r$out[i] = min(d[i,r$species[i],close])
	}
	r
}

make4DExpMatrix = function(d,m,age.field='mouse.days'){
	m = m[colnames(d),]
	ss = unique(m$species)
	ts = unique(m$tissue)
	as = sort(unique(m[,age.field]))
	r = array(NA,dim=c(nrow(d),length(ss),length(ts),length(as)),dimnames = list(rownames(d),ss,ts,as))
	for(s in ss)
		for(t in ts)
			for(a in as){
				inx=which(m$species==s & m$tissue==t & m[,age.field]==a)
				if(length(inx)!=0)
					r[,s,t,as.character(a)] = d[,inx]
			}
	r
}


caclInterSpeciesPerTissueDist = function(d,min.obs,mean.diff=FALSE){
	species = dimnames(d)[[2]]
	tissues = dimnames(d)[[3]]
	r = array(NA,dim=c(dim(d)[1],length(tissues),length(species),length(species)),dimnames=list(dimnames(d)[[1]],tissues,species,species))
	for(t in tissues){
		print(t)
		for(i in 1:nrow(d)){
			m = d[i,,t,]
			m = m[,apply(is.na(m),2,sum)==0,drop=FALSE]
			if(ncol(m)>=min.obs){
				if(!mean.diff)
					m = sweep(m,1,apply(m,1,mean))
				for(s1 in 1:(length(species)-1))
					for(s2 in (s1+1):length(species)){
						r[i,t,s1,s2] = r[i,t,s2,s1] = ifelse(mean.diff,abs(mean(m[s1,])-mean(m[s2,])),mean(abs(m[s1,]-m[s2,])))
					}
			}
		}
	}
	r
}

hgmdOverlapLocaction = function(h,s,from,to,upstream=TRUE){
	hstop = (h$pos_VCF_hg19+nchar(h$ref_VCF_hg19)-1)
	t2t = ifelse(s$strand == 1,h$pos_VCF_hg19-s$start,s$stop - hstop)
	p2t = ifelse(s$strand == 1,hstop         -s$start,s$stop - h$pos_VCF_hg19)
	
	p2p = ifelse(s$strand == 1,s$stop-h$pos_VCF_hg19,hstop         -s$start)
	t2p = ifelse(s$strand == 1,s$stop-hstop         ,h$pos_VCF_hg19-s$start)
	
	ifelse(rep(upstream,nrow(h)),to >= t2t & from <= p2t,to >= t2p & from <= p2p)
}

plotHGMDDen = function(h2s,mut,seg,cod,tags,int.mar=200,cl.names=c('-2'='const.','-1'='alt.'),col=setNames(c('green','orange','gray',rep(c('red','blue'),times=8),rep('gray',9)),-2:25),comp.to=NULL,cl.order=NULL){
	seg = seg[seg$cod %in% cod,]
	gid = unique(seg$gene_id)
	seg = seg[seg$gene_id %in% gid,]
	h2s = h2s[h2s$seg.id %in% rownames(seg) & mut[h2s$hgmd.id,'tag'] %in% tags,]
	cl.sizes = table(seg$cl)
	cl.nt.sizes = sapply(split(seg$length,seg$cl),sum)
	
	cl.mut.cnt = table(factor(seg[h2s$seg.id,'cl'],levels=names(col)),h2s$position)
	
	par(mfrow=c(3,1),tck=-0.01,mgp=c(1.7,0.4,0),mar=c(5,3,1.5,0),oma=c(0,0,2,1))
	plotBinomBar(cl.mut.cnt[,'u'],cl.sizes*int.mar,col=col[rownames(cl.mut.cnt)],comp.to=comp.to,cl.names=cl.names,main=paste0('[-',int.mar,',-1] upstream intron'),cl.sizes,cl.order=cl.order)
	plotBinomBar(cl.mut.cnt[,'d'],cl.sizes*int.mar,col=col[rownames(cl.mut.cnt)],comp.to=comp.to,cl.names=cl.names,main=paste0('[-',int.mar,',-1] downstream intron'),cl.sizes,cl.order=cl.order)
	plotBinomBar(cl.mut.cnt[,'i'],cl.nt.sizes,col=col[rownames(cl.mut.cnt)],comp.to=comp.to,cl.names=cl.names,main='Exon',cl.sizes,cl.order=cl.order)
	invisible(cl.mut.cnt)
}

getPhastconsProf = function(d,s1,s2){
	d = d[!sapply(d,is.null)]
	cbind(c(apply(sapply(d,function(x){x[s1:s2]}),1,mean),apply(sapply(d,function(x){x[(length(x)-s2):(length(x)-s1)]}),1,mean)),
				c(apply(sapply(d,function(x){x[s1:s2]}),1,sd  ),apply(sapply(d,function(x){x[(length(x)-s2):(length(x)-s1)]}),1,sd  ))/sqrt(length(d)))
}

plotArea = function(x,p,col,sd.mult=2,new=FALSE,ylim=NULL,xlim=range(x),area.transp=0.2,type='l',area.den=-1,...){
	#p should contain either mean and sd
	#or mean, lower and upper bounds
	o = order(x)
	x = x[o]
	p = p[o,]
	na = !is.na(p[,1])
	x = x[na]
	p = p[na,]
	if(ncol(p)==2)
		yp = c(p[,1]-p[,2]*sd.mult,rev(p[,1]+p[,2]*sd.mult))
	else
		yp = c(p[,2],rev(p[,3]))
	if(new){
		if(is.null(ylim))
			ylim = range(yp,na.rm=T)
		plot(1,t='n',xlim=xlim,ylim=ylim,...)
	}
	col.pol = col2rgb(col,alpha=TRUE)/255
	col.pol = rgb(col.pol[1],col.pol[2],col.pol[3],col.pol[4]*area.transp) 
	polygon(c(x,rev(x)),yp,col=col.pol,border=NA,den=area.den)
	lines(x,p[,1],col=col,type=type,...)
}

plotPhastConsProfs = function(p,s,s1,s2,cod=c('c','n','p')){
	getProf = function(d){getPhastconsProf(d,s1,s2)}
	
	gid = unique(s$gene_id[s$cod %in% cod & s$cl > -1])
	s = s[s$cod %in% cod & s$gene_id %in% gid,]
	p = p[rownames(s)]
	
	exn = getProf(p[s$cl==-2])
	alt = getProf(p[s$cl==-1])
	x = 1:((s2-s1+1)*2)
	for(i in 0:7){
		u = getProf(p[s$cl==1+i*2])
		d = getProf(p[s$cl==2+i*2])
		plot(1,t='n',xlim=range(x),ylim=c(0,1),main=paste0('Clusters ',1+i*2,'-',2+i*2),ylab='mean primate phastcons',xaxt='n',xlab='position (nt)')
		plotArea(x,u,col='red')
		plotArea(x,d,col='blue')
		plotArea(x,exn,col='green',lty=1)
		plotArea(x,alt,col='orange',lty=1)
		axis(1,c(1,200-s1+1,2*(s2-s1+1)-200+s1,length(x)),c(-(200-s1+1),'acc','don',(200-s1+1)))
		abline(v=c(200-s1+1,s2-s1+1.5,2*(s2-s1+1)-200+s1),lty=2)
		grid()
		legend('topleft',col=c('red','blue','green','yellow'),lwd=1,legend=c('up','down','const.','alt.'))
	}
	for(i in 17:25){
		u = getProf(p[s$cl==i])
		plot(1,xlim=range(x),t='n',ylim=c(0,1),main=paste0('Cluster ',i),ylab='mean primate phastcons',xaxt='n',xlab='position (nt)')
		plotArea(x,u,col='red')
		plotArea(x,exn,col='green',lty=1)
		plotArea(x,alt,col='orange',lty=1)
		axis(1,c(1,200-s1+1,2*(s2-s1+1)-200+s1,length(x)),c(-(200-s1+1),'acc','don',(200-s1+1)))
		abline(v=c(200-s1+1,s2-s1+1.5,2*(s2-s1+1)-200+s1),lty=2)
		grid()
	}
}

plotBinomBar = function(s,t,col,cl.names,cl.sizes,cl.order,pv.thr=0.05,comp.to,...){
	if(!is.null(cl.order)){
		col = col[cl.order[cl.order %in% names(s)]]
		s = s[cl.order[cl.order %in% names(s)]]
	}
	t = t[names(s)]
	t[is.na(t)] = 0
	p = as.numeric(s/t)
	p[s==0] = 0
	sd = sqrt(p*(1-p)/t)*2
	sd[p==0]=0
	b = 1:length(p)
	plot(b,p,col=col,pch=19,ylim=range(pmax(0,p-sd),p+sd),xlab='',xaxt='n',ylab='gnomAD SNP density (SNP per nt)',...)
	segments(b,pmax(0,p-sd),b,p+sd,col=col)
	y= par("usr")[3]
	l = names(s)
	l[l %in% names(cl.names)] = cl.names[l[l %in% names(cl.names)]]
	cl.sizes[is.na(cl.sizes)]=0
	l = paste0(l,' (',cl.sizes[names(s)],')')
	text(b,y,l,srt=45,xpd=T,adj = c(1,0.5))
	if(!is.null(comp.to)){
		y = range(p-sd,p+sd)
		y = seq(y[2],y[2]-(y[2]-y[1])*0.1,length.out=length(comp.to))
		for(i in 1:length(comp.to)){
			pv=sapply(1:length(t),function(j){
				if(min(c(t[comp.to[i]],t[j])) == 0)
					return(1)
				prop.test(c(s[comp.to[i]],s[j]),c(t[comp.to[i]],t[j]))$p.value
			})
			
			bi = b[p.adjust(pv,m='BH')<pv.thr]
			points(bi,rep(y[i],length(bi)),pch='*',col=col[comp.to[i]],cex=2)
		}
	}
}

plotSNPDen = function(g,s,title='',int.mar=100,cl.names=c('-2'='const.','-1'='alt.'),col=setNames(c('green','orange','gray',rep(c('red','blue'),times=8),rep('gray',12)),-2:28),cl.order=NULL,comp.to=NULL){
	cl.sizes = table(s$cl)
	cl.nt.sizes = sapply(split(s$length,s$cl),sum)
	up = table(s[g$seg_id[g$dist2start<0 & g$dist2start> -int.mar],'cl'])
	dw = table(s[g$seg_id[g$dist2stop<0 & g$dist2stop> -int.mar],'cl'])
	i = table(s[g$seg_id[g$dist2start>=0 & g$dist2stop>=0],'cl'])
	
	par(mfrow=c(3,1),tck=-0.01,mgp=c(1.7,0.4,0),mar=c(5,3,1.5,0),oma=c(0,0,2,1))
	plotBinomBar(up,cl.sizes*int.mar,col=col[names(up)],comp.to=comp.to,cl.names=cl.names,main=paste0('[-',int.mar,',-1] upstream intron'),cl.sizes,cl.order=cl.order)
	plotBinomBar(dw,cl.sizes*int.mar,col=col[names(up)],comp.to=comp.to,cl.names=cl.names,main=paste0('[1,',int.mar,'] downstream intron'),cl.sizes,cl.order=cl.order)
	plotBinomBar(i,cl.nt.sizes,col=col[names(up)],comp.to=comp.to,cl.names=cl.names,main='Exon',cl.sizes,cl.order=cl.order)
	mtext(title,outer = T)
}




plotSNPfreBoxplotByCl = function(f,col,cl.names,cl.order,comp.to,pv.thr=0.05,xlab='',all.cl=NULL,...){
	t = split(log10(f$freq),f$cl)
	if(!is.null(all.cl)){
		t$all=unlist(t[all.cl])
	}
	if(!is.null(cl.order))
		t = t[cl.order[cl.order %in% names(t)]]
	#boxplot(t,notch=T,col=col[names(t)],outline=F,ylab='log10(Allele freq.)',xaxt='n',...)
	#t = t[as.character(-2:25)]
	x = 1:length(t)
	y=sapply(t,mean)
	sd = sapply(t,sd)*1.96/sqrt(sapply(t,length))
	col=col[names(t)]
	plot(x,y,col=col,pch=19,ylab='log10(Allele freq.)',xaxt='n',ylim=range(y-sd,y+sd,na.rm=T),xlab=xlab,...)
	segments(x,y-sd,x,y+sd,col=col)
	l = names(t)
	l[l %in% names(cl.names)] = cl.names[l[l %in% names(cl.names)]]
	b = 1:length(t)
	yleg= par("usr")[3]
	text(b,yleg,l,srt=45,xpd=T,adj = c(1,0.5))
	
	if(!is.null(comp.to)){
		y = range(y-sd,y+sd,na.rm=T)
		y = seq(y[2],y[2]-(y[2]-y[1])*0.1,length.out=length(comp.to))
		for(i in 1:length(comp.to)){
			pv=sapply(1:length(t),function(j){
				wilcox.test(t[[comp.to[i]]],t[[j]])$p.value
			})
			pv=p.adjust(pv,m='BH')
			bi = x[pv<pv.thr]
			#print(pv)
			points(bi,rep(y[i],length(bi)),pch='*',col=col[comp.to[i]],cex=2)
		}
	}
	invisible(t)
}

plotSNPFreq = function(g,s,title='',int.mar=100,cl.names=c('-2'='const.','-1'='alt.'),col=setNames(c('green','orange','gray',rep(c('red','blue'),times=8),rep('gray',12)),-2:28),cl.order=NULL,comp.to=NULL){
	g = g[g$alt_cnt>1,]
	par(mfrow=c(3,1),tck=-0.01,mgp=c(1.7,0.4,0),mar=c(3,3,1.5,0),oma=c(0,0,2,1))
	plotSNPfreBoxplotByCl(g[g$dist2start<0 & g$dist2start> -int.mar,],comp.to=comp.to,cl.names=cl.names,col=col,main='[-100,-1] upstream intron',cl.order=cl.order)
	plotSNPfreBoxplotByCl(g[g$dist2stop<0 & g$dist2stop> -int.mar,],comp.to=comp.to,cl.names=cl.names,col=col,main='[1,100] downstream intron',cl.order=cl.order)
	plotSNPfreBoxplotByCl(g[g$dist2start>=0 & g$dist2stop>=0,],comp.to=comp.to,cl.names=cl.names,col=col,main='Exon',cl.order=cl.order)
	mtext(title,outer = T)
}

getReadCoverage = function(bams,chr,start,end,strand=NA,plot=FALSE,min.junc.cov=0,plot.junc.only.within=FALSE,ylim=NULL,reverse=FALSE,...){
	stop("deprecated, see util.R")
	if(start>end){
		t = start
		start = end
		end=t
	}
		
	require(GenomicAlignments)
	param = ScanBamParam(flag=scanBamFlag(isMinusStrand=strand==-1),which=GRanges(chr, IRanges(start, end)))
	r = NULL
	for(b in bams){
		bam = readGAlignments(b,param = param)
		cov=coverage(bam)[[chr]][start:end]
		juncs = as.data.frame(summarizeJunctions(bam))
		rownames(juncs)=paste(juncs$seqnames,juncs$start,juncs$end,sep='-')
		if(is.null(r))
			r=list(cov=cov,juncs=juncs)
		else{
			r$cov = r$cov + cov
			cmn = intersect(rownames(juncs),rownames(r$juncs))
			r$juncs[cmn,'score'] = r$juncs[cmn,'score'] + juncs[cmn,'score']
			r$juncs = rbind(r$juncs,juncs[setdiff(rownames(juncs),rownames(r$juncs)),])
		}
	}
	
	r$juncs = r$juncs[r$juncs$score>=min.junc.cov,]
	if(plot){
		x = start:end
		r$cov[c(1,length(r$cov))] = 0
		if(is.null(ylim))
			ylim = c(0,max(r$cov))
		xlim=range(x)
		if(reverse)
			xlim=rev(xlim)
		plot(x,r$cov,t='n',ylim=ylim,xlim=xlim,...)
		polygon(x,r$cov,col = 'gray',border=NA)
		for(i in 1:nrow(r$juncs))
			if(!plot.junc.only.within || (r$juncs$start[i] > start & r$juncs$end[i]< end))
				plotArc(r$juncs$start[i],r$juncs$end[i],r$juncs$score[i],col='red',lwd=3)
	}
	invisible(r)
}

plotArc = function(from,to,top,n=100,y.base=0,...){
	len = to - from
	x = seq(from=0,to=len,length.out = n)
	y = x*4*top/len - x^2*(4*top/len^2)
	lines(x+from,y+y.base,...)
}


plotAllTrees = function(ts,mouse.days,m){
	for(s in names(ts)){
		for(d in as.character(sort(mouse.days))){
			if(d %in% names(ts[[s]])){
				days = round(mean(m[m$species==s & m$mouse.days==d,'days'],na.rm=T),1)
				if(days > 365)
					days = paste(round(days/365,1),'years')
				else
					days = paste(days,'days')
				plotTree(ts[[s]][[d]],col,T,xmax=max.len,root.max=max.root.len,root.name='base',xaxs=T,lwd=2,main=paste(s,days),direction='upwards',srt=-90,adj=0.5)
			}else
				plot.new()
		}
	}
}

plotTrees = function(ts,cols,plot.steem,root.name='base',xaxs=FALSE,lwd=1,direction='rightwards',...){
	root.len = rep(0,length(ts))
	if(!is.na(root.name)){
		ts = lapply(ts,root,outgroup=root.name,resolve.root=T)
		root.len = sapply(ts,function(t){t$edge.length[t$edge[,2] == which(t$tip.label==root.name)]})
		ts = lapply(ts,drop.tip,tip=root.name)
	}
	if(!plot.steem)
		root.len[] = 0
	max.len = max(unlist(lapply(ts,node.depth.edgelength)))
	
	for(i in 1:length(ts)){
		plotTree(ts[[i]],cols,plot.steem,max.len*1.05,max(root.len),NULL,xaxs,lwd,root.len[i],main=paste('Day',names(ts)[i]),direction=direction,...)
	}
}

plotTree = function(t,col,plot.steem,xmax=NULL,root.max=NULL,root.name=NULL,xaxs=FALSE,lwd=2,root.len = 0,xmax.prop = 1.05,main='',direction='rightwards',...){
	if(!is.null(root.name)){
		t = root(t,outgroup=root.name,resolve.root=T)
		root.len = t$edge.length[t$edge[,2] == which(t$tip.label==root.name)]
		t = drop.tip(t,tip=root.name)
	}
	if(!plot.steem)
		root.len = 0
	if(is.null(root.max))
		root.max = root.len
	if(is.null(xmax))
		xmax = max(node.depth.edgelength(t))*xmax.prop
	edge.col = rep('gray',nrow(t$edge))
	f = t$edge[,2]<= length(t$tip.label)
	edge.col[f] = col[t$tip.label[t$edge[f,2]]]
	if(sum(c('b','c') %in% t$tip.label)==2){
		cb=t$edge[which(t$tip.label =='c')==t$edge[,2],1]
		edge.col[t$edge[,2]==cb] = col['b']
	}
	
	pars=list(x=t,x.lim=c(-root.len,xmax+root.max-root.len),root.edge=F,cex=1.5,main=main,type='p',tip.color=col[t$tip.label],edge.color=edge.col,edge.width=lwd,direction=direction,...)
	if(direction=='upwards')
		names(pars)[names(pars)=='x.lim']='y.lim'
	do.call(plot.phylo,pars)
	# n = chronos(t)
	# n$edge.length[n$edge.length<0] = 0
	# n=table(cutree(as.hclust(n),k=2))
	# n = n[1]+0.5
	n = length(t$tip.label)/2 + 1.5
	if(direction=='rightwards')
		lines(c(-root.len,0),c(n,n),lwd=lwd,col='gray')
	else
		lines(c(n,n),c(-root.len,0),lwd=lwd,col='gray')
	if(xaxs){
		at = seq(0,xmax+root.max,by = 0.05)
		axis(switch(direction,rightwards=1,upwards=2),at = at - root.len,labels=at)
	}
}


areaplot = function(d,x=1:ncol(d),col=gray(0:(nrow(d)-1)/nrow(d)),...){
	plot(1,t='n',xlim=range(x),ylim=range(0,apply(d,2,sum)),xaxs='i',yaxs='i',...)
	y = rep(0,ncol(d))
	x = c(x,rev(x))
	for(i in 1:nrow(d)){
		polygon(x,c(y,rev(y+d[i,])),border = NA,col=col[i])
		y = y + d[i,]
	}
	at =cumsum(d[,ncol(d)]) - d[,ncol(d)]/2
	#axis(4,at = at,labels = rownames(d))
	text(rep(max(x),length(at)),at,labels = rownames(d),adj=c(-0.1,0.5),xpd=TRUE)
}

getConservationOfAltExons = function(psi,m,s1,s2,tissues,psi.thr,summ=TRUE){
	m1 = m[m$tissue %in% tissues & m$species==s1,]
	m2 = m[m$tissue %in% tissues & m$species==s2,]
	st1 = paste(m1$tissue,m1$mouse.days)
	st2 = paste(m2$tissue,m2$mouse.days)
	st = intersect(st1,st2)
	s = sapply(strsplit(st,' ',T),'[',2)
	t = table(s)
	f = s %in% names(t)[t==length(tissues)]
	m1 = m1[st1 %in% st[f], ]
	m2 = m2[st2 %in% st[f], ]
	days = sort(unique(m1$mouse.days))
	r1 = r12 = list()
	for(d in days){
		p1 = psi[,rownames(m1)[m1$mouse.days==d]]
		t1 = substr(m1[colnames(p1),'tissue'],1,1)
		t1 = apply(p1>psi.thr & p1 < (1-psi.thr),1,function(x){paste(sort(t1[x]),collapse='')})
		
		p2 = psi[,rownames(m2)[m2$mouse.days==d]]
		t2 = substr(m2[colnames(p2),'tissue'],1,1)
		t2 = apply(p2>psi.thr & p2 < (1-psi.thr),1,function(x){paste(sort(t2[x]),collapse='')})
		r1[[as.character(d)]] = t1
		r12[[as.character(d)]] = t1[t1==t2]
	}
	list(total=getAltExonStat2table(r1,summ),conserved=getAltExonStat2table(r12,summ))
}



plotAltConsMouse.to = function(psi,m,psi.thr,main){
	r = list()
	r$rat     = getConservationOfAltExons(psi,m,'mouse','rat',c('brain','heart','liver','ovary','testis'),psi.thr)
	r$rabbit  = getConservationOfAltExons(psi,m,'mouse','rabbit',c('brain','heart','liver','ovary','testis'),psi.thr)
	r$human   = getConservationOfAltExons(psi,m,'mouse','human',c('brain','heart','liver','testis'),psi.thr)
	r$opossum = getConservationOfAltExons(psi,m,'mouse','opossum',c('brain','heart','liver','ovary','testis'),psi.thr)
	r$chicken = getConservationOfAltExons(psi,m,'mouse','chicken',c('brain','heart','liver','ovary','testis'),psi.thr)
	for(s in names(r)){
		x=as.numeric(colnames(r[[s]]$total))
		areaplot(r[[s]]$total                 ,x = x,col = col[rownames(r[[s]]$total)],log='x',xlab='Age (days)',ylab='# of segments',main='Mouse')
		areaplot(r[[s]]$conserved             ,x = x,col = col[rownames(r[[s]]$total)],log='x',xlab='Age (days)',ylab='# of segments',main=paste('Conserved in',s))
		y = r[[s]]$conserved/r[[s]]$total
		plot(1,t='n',ylim=c(0,max(y)),xlim=range(x),log='x',xlab='Age (days)',ylab='prob. of conserved',main=paste('Conserved in',s))
		for(i in 1:nrow(y)){
			lines(x,y[i,],col=col[rownames(y)[i]],lwd=3)
		}
	}
	mtext(paste(psi.thr,' < PSI < ',1-psi.thr,'; ',main,sep=''),3,outer=TRUE)
}

plotNumberOfAltExons = function(pv.thr,good.segs){
	z = getAltExonStat(psi.tsm$human[rownames(psi.tsm$human) %in% good.segs,],meta.tsm,pv.thr,tissues = c('brain','heart','liver','testis'))
	z = z[,apply(z==0,2,sum)==0]
	areaplot(z,x = as.numeric(colnames(z)),col = col[rownames(z)],log='x',xlab='Age (days)',ylab='# of segments',main='Human')
	
	z = getAltExonStat(psi.tsm$macaque[rownames(psi.tsm$macaque) %in% good.segs,],meta.tsm,pv.thr,tissues = c('brain','heart','kidney','liver','testis'))
	areaplot(z,x = as.numeric(colnames(z)),col = col[rownames(z)],log='x',xlab='Age (days)',ylab='# of segments',main='Macaque')
	
	z = getAltExonStat(psi.tsm$mouse[rownames(psi.tsm$mouse) %in% good.segs,],meta.tsm,pv.thr,tissues = c('brain','heart','kidney','liver','ovary','testis'))
	areaplot(z,x = as.numeric(colnames(z)),col = col[rownames(z)],log='x',xlab='Age (days)',ylab='# of segments',main='Mouse')
	
	z = getAltExonStat(psi.tsm$rat[rownames(psi.tsm$rat) %in% good.segs,],meta.tsm,pv.thr,tissues = c('brain','heart','liver','ovary','testis'))
	z = z[,-1]
	areaplot(z,x = as.numeric(colnames(z)),col = col[rownames(z)],log='x',xlab='Age (days)',ylab='# of segments',main='Rat')
	
	z = getAltExonStat(psi.tsm$rabbit[rownames(psi.tsm$rabbit) %in% good.segs,],meta.tsm,pv.thr,tissues = c('brain','heart','ovary','testis'))
	z = z[,apply(z==0,2,sum)==0]
	areaplot(z,x = as.numeric(colnames(z)),col = col[rownames(z)],log='x',xlab='Age (days)',ylab='# of segments',main='Rabbit')
	
	z = getAltExonStat(psi.tsm$opossum[rownames(psi.tsm$opossum) %in% good.segs,],meta.tsm,pv.thr,tissues = c('brain','heart','liver','ovary','testis'))
	z = z[,apply(z==0,2,sum)==0]
	areaplot(z,x = as.numeric(colnames(z)),col = col[rownames(z)],log='x',xlab='Age (days)',ylab='# of segments',main='Opossum')
	
	z = getAltExonStat(psi.tsm$chicken[rownames(psi.tsm$chicken) %in% good.segs,],meta.tsm,pv.thr,tissues = c('brain','heart','kidney','liver','ovary','testis'))
	areaplot(z,x = as.numeric(colnames(z)),col = col[rownames(z)],log='x',xlab='Age (days)',ylab='# of segments',main='Chicken')
}

plotLine = function(x,y,cor.method='pearson',line.col='red',leg.pos='topright',line.lwd=1,...){
	stop("use util.R")
	plot(x,y,...)
	abline(lm(y~x),col=line.col,lwd=line.lwd)
	c = cor.test(x,y,m=cor.method)
	ci = round(c$conf.int,2)
	leg=paste('rho=',round(c$estimate,2),' [',ci[1],',',ci[2],']; pv=',format(c$p.value,digits=2,scientific=TRUE),sep='')
	legend(leg.pos,lwd=line.lwd,col=line.col,legend=leg,bty='n',text.col=line.col)
}


classifyNonOrthEns = function(s2e,orth.segs,o,sp){
	#o is ens orth genes
	gids = unlist(s2e[[sp]][rownames(orth.segs[[sp]])[orth.segs$human$ens.orth.filter==1]])
	
	r = setNames(rep('',length(gids)),names(gids))
	d = rep(F,nrow(o))
	for(i in 1:ncol(o)){
		t = table(o[o[,i] != '',i])
		d = d | o[,i] %in% names(t)[t > 1]
	}
	r[gids %in% o[d,sp]] = 'duplicated'
	r[!(gids %in% o[,sp])] = 'human-lost'
	o. = o[!d & o[,sp] %in% gids,]
	rownames(o.) = o.[,sp]
	o. = apply(o.=='',1,sum)
	r[gids %in% names(o.)[o. == 0]] = 'contradict'
	r[gids %in% names(o.)[o. == 1]] = 'lost in one'
	r[gids %in% names(o.)[o. > 1]] = 'lost in more'
	r
}

plotOval = function(x0,y0,rx,ry,n=200,...){
	a = seq(0,2*pi,length.out=n)
	polygon(x0+rx*cos(a),y0+ry*sin(a),...)
}

plotClusterSSS = function(a,d,c,clusts = 1:25,main='',d.const=NULL,a.const=NULL){
	c = unlist(c)
	a = unlist(a)
	d = unlist(d)
	f = c %in% clusts
	a = a[f]
	d = d[f]
	c = c[f]
	ba = boxplot(a~c,plot=F,notch=T)
	bd = boxplot(d~c,plot=F,notch=T)
	const = NULL
	if(!is.null(d.const))
		const = boxplot(as.numeric(a.const),as.numeric(d.const),plot=F,notch=T)
	plot(1,t='n',xlim=range(ba$conf,const$conf[,1]),ylim=range(bd$conf,const$conf[,2]),main=main,xlab='acceptor strength',ylab='donor strength')
	x = ba$stats[3,]
	y = bd$stats[3,]
	for(i in clusts)
		plotOval(x[i],y[i],(ba$conf[2,i]-ba$conf[1,i])/2,(bd$conf[2,i]-bd$conf[1,i])/2,col='#00000010',border=NA)
	
	arrows(x[1+(0:7)*2],y[1+(0:7)*2],x[2+(0:7)*2],y[2+(0:7)*2])
	text(x,y,clusts,col=c(rep(c('red','blue'),times=8),rep('black',20)))
	if(!is.null(d.const)){
		plotOval(const$stats[3,1],const$stats[3,2],(const$conf[2,1]-const$conf[1,1])/2,(const$conf[2,2]-const$conf[1,2])/2,col='#00FF0020',border=NA)
		text(const$stats[3,1],const$stats[3,2],'const')	
	}
}

plotSSStrengthScatterPNG = function(sss,main,fname){
	png(fname,width=12,height=12,units='in',res=300)
	pairs(sss,pch='.',panel = function(x,y,...){points(x,y,...);abline(v=0,col='red');legend('topleft',bty='n',legend=paste('r=',round(cor(x,y,m='p',u='p'),3),sep=''),col='red',lwd=1,);abline(h=0,col='red');abline(a=0,b=1,col='red')},main=paste(main,' (',nrow(sss),')',sep=''))
	dev.off()	
}

plotSignPairs = function(sgn,means,f,amp.thrs = c(0.1,0.5,1),...){
	getSgnCnts = function(s,m,t1,t2){
		r = numeric(ncol(m)*ncol(m))
		names(r) = paste(rep(dimnames(sgn)[[1]],each=ncol(m)),rep(dimnames(sgn)[[1]],times=ncol(m)),sep='')
		for(i in 1:ncol(m))
			for(j in 1:ncol(m)){
				if(t1 > 0)
					cnt = sum(sgn[i,j,] & (means[,i] - means[,j]) > t1 & (means[,i] - means[,j]) < t2 ,na.rm=T)
				else
					cnt = sum(sgn[i,j,] & (means[,i] - means[,j]) < t1 & (means[,i] - means[,j]) > t2,na.rm=T)
				r[paste(dimnames(sgn)[[1]][i],dimnames(sgn)[[1]][j],sep='')] = cnt
			}
		r
	}
	
	sgn = sgn[,,f]
	means = means[f,]
	countsp = countsn = list()
	for(i in 1:(length(amp.thrs)-1))
		countsp[[paste('p',amp.thrs[i],sep='')]] = getSgnCnts(sgn,means, amp.thrs[i], amp.thrs[i+1])
	for(i in 1:(length(amp.thrs)-1))
		countsn[[paste('n',amp.thrs[i],sep='')]] = getSgnCnts(sgn,means,-amp.thrs[i],-amp.thrs[i+1])
	countsn = do.call(rbind,countsn)
	countsp = do.call(rbind,countsp)
	
	spaces = rep(c(rep(0,ncol(means)-1),2),times=ncol(means))
	spaces = c(0,spaces[1:(length(spaces)-1)])
	b=barplot( countsp,xaxt='n',ylim=c(-max(apply(countsn,2,sum)),max(apply(countsp,2,sum)))*1.1,space = spaces,ylab='# of exons',legend.text=amp.thrs[1:(length(amp.thrs)-1)],args.legend=list(x='topleft',title='dPSI'),...)
	barplot(-countsn,xaxt='n',ylim=c(-max(apply(countsn,2,sum)),max(apply(countsp,2,sum))),space = spaces,add=T)
	axis(1,b,substr(colnames(countsn),2,2),las=3)
	sp = sapply(split(b,rep(1:ncol(means),each=ncol(means))),mean)
	axis(1,sp,colnames(means),line=1,tick=F,lwd=0)
	invisible(list(countsn,countsp))
}


countNumberOfCassettesInAlt = function(a){
	r = numeric(nrow(a))
	r[0] = 0
	for(i in 1:nrow(a)){
		if(is.na(a$sites[i])) next
		ss = strsplit(a$sites[i],'',TRUE)[[1]]
		segs = strsplit(a$segs[i],';',TRUE)[[1]]
		segs = segs[!grepl('NA',segs)]
		if(length(segs)>0){
			segs = do.call(rbind,lapply(strsplit(segs,'-',TRUE),as.numeric))+1
			r[i] = sum(ss[segs[,1]] == 'a' & ss[segs[,2]] == 'd')
		}
	}
	names(r) = rownames(a)
	r
}

calcProportionOfLiftoveredExons = function(anns,o,cod.seg,s){
	require(data.table)
	ann = anns[[s]]
	o = o[,s]
	ann = ann[ann$sites=='ad' & ann$cod %in% cod.seg,]
	ann$has.orth = rownames(ann) %in% o
	data.frame(rbindlist(lapply(split(ann,ann$gene_id),function(x){
		data.frame(orth.freq=mean(x$has.orth),cod.gene=x$cod.gene[1],exn.cnt=nrow(x))
	})))
}

getDistanceToClosestSegs = function(d,min.psi){
	print('t\tgetDistanceToClosestSegs')
	d$seg$f = apply(d$ir,1,function(x){x=x[!is.na(x)];length(x)>1 & sort(x,decreasing = TRUE)[2]>=min.psi})
	d = d$seg
	gc()
	d$seg_id = rownames(d)
	order = rownames(d)
	d$up.dist = d$dw.dist = NA
	d = split(d,d$gene_id)
	k=1
	d = lapply(d,function(x){
		cat('\r',k)
		k <<- k + 1
		if(is.na(x$strand[1])) return(x)
		x = x[order(x$start),]
		#forward
		prev.pos = NA
		for(i in 1:nrow(x)){
			x$up.dist[i] = x$start[i] - prev.pos - 1
			if(x$f[i])
				prev.pos = x$stop[i]
		}
		#backward
		prev.pos = NA
		for(i in nrow(x):1){
			x$dw.dist[i] = prev.pos - x$stop[i] - 1
			if(x$f[i])
				prev.pos = x$start[i]
		}
		if(x$strand[1] == -1){
			x[,c('up.dist','dw.dist')] = x[,c('dw.dist','up.dist')]
		}
		x
	})
	d = as.data.frame(data.table::rbindlist(d))
	rownames(d) = d$seg_id
	d$seg_id = NULL
	d[order,]
}

getSegmentAltInfo = function(hsids){
	r = data.frame(seg.id=hsids,ens.id = NA,descr=NA,name=NA)
	for(i in 1:nrow(r)){
		eids = seg2ens$human[[hsids[i]]]
		if(length(eids)>0){
			r$descr[i] = paste(ens.descr[eids,'Description'],collapse='; ')
			r$name[i] = paste(ens.descr[eids,'Associated.Gene.Name'],collapse='; ')
		}
		r$ens.id[i] = paste(eids,collapse='; ')
	}
	r
}

plotSpSpMainTissueCount = function(e,cl,f,main,n=2){
	for(sp in c(species$short,'hq','mr','oc')){
		t = e[cl==sp & f,]
		t = table(factor(colnames(e)[unlist(apply(t,1,function(x){o = order(abs(x),decreasing = TRUE)[1:n];o[!is.na(x[o]) & abs(x[o])>0]}))],levels=colnames(e)))[colnames(e)]
		barplot(t,col=cols$tissue[colnames(e)],ylab=paste('Number of tissues in top ',n,sep=''),xlab='',main=sp,las=3)
	}
	mtext(main,3,0,TRUE)
}

plotSpSpEffectDistr = function(e,cl,f,main,ylim=c(-1,1)){
	for(sp in c(species$short,'hq','mr','oc')){
		t = e[cl==sp & f,]
		boxplot(t,col=cols$tissue[colnames(e)],ylab='species specific effect',xlab='',main=sp,ylim=ylim,las=3)
	}
	mtext(main,3,0,TRUE)
}

calcSpSpEffectPerTissue = function(psi,spec){
	m = do.call(rbind,strsplit(colnames(psi),' '))
	tissues = unique(m[,2])
	r = matrix(NA,nrow=nrow(psi),ncol=length(tissues),dimnames = list(rownames(psi),tissues))
	for(t in tissues){
		psi_ = psi[,m[,2]==t]
		colnames(psi_) = species[m[m[,2]==t,1],'short']
		for(s in 1:nrow(r)){
			if(spec[s] %in% c(species$short,'hq','mr','mrb','oc') & sum(is.na(psi_[s,])) == 0){
				g1 = strsplit(spec[s],'')[[1]]
				g2 = setdiff(species$short,g1)
				g1 = range(psi_[s,g1],na.rm=TRUE)
				g2 = range(psi_[s,g2],na.rm=TRUE)
				if(g1[1] <= g2[2] & g1[2] >= g2[1])
					r[s,t] = 0
				else if(g1[1] > g2[2])
					r[s,t] = g1[1]-g2[2]
				else 
					r[s,t] = g1[2]-g2[1]
			}
		}
	}
	r
}


countSpSpClasses =  function(sp,cl){
	ne = apply(cl!='e',1,sum)
	f = sp$spec != 'NA' & ne >= 1
	sp = sp[f,]
	ne = ne[f]
	class = sp$spec
	class[sp$spec.short %in% c('wrong','s12','s21')] = 'complex'
	class[sp$spec.short %in% c('n','s5')] = 'no change'
	colnames(cl) = species[colnames(cl),'short']
	group.order=c('no change',species$short,'hq','mr','mrb','oc','complex')
	return(table(ne==7,factor(class,levels = group.order))[,group.order])
}


plotOrthSegPSI = function(psi,meta,orth,id,...){
	ids = orth[which(apply(orth==id,1,sum)>0),]
	par(mfrow=c(3,3),tck=-0.02,mgp=c(1.3,0.4,0),mar=c(3,6,1.5,0),oma=c(0,0,2,1))
	for(i in 1:7){
		if(is.list(psi[[i]]))
			p = psi[[i]]$ir[ids[i],]
		else
			p = psi[[i]][ids[i],]
		plotTissueAgeProile(p,meta,main=ids[i],...)
	}
}

findOrthExonsByNeighbors = function(a,o,s='human',returnBad=FALSE){
	rownames(o) = o[,s]
	n = names(a)
	a = lapply(n,function(s){x=a[[s]][a[[s]]$sites =='ad',];x$has.orth = rownames(x) %in% o[,s];lapply(split(x,x$gene_id),function(y){y[order(y$exon.number),]})})
	names(a) = n
	bad = 0
	cases = 0
	res = new.env()
	bad.data = new.env()
	for(i in 1:length(a[[s]])){
		cat('\r',i,length(a[[s]]),'bad = ',bad,'; tested =',cases,'       ')
		g = a[[s]][[i]]
		loable = which(g$has.orth)
		if(length(loable) > 0 && length(loable) < (max(loable)-min(loable)+1)){
			for(e in 2:length(loable)){
				if(loable[e-1] + 1 < loable[e]){
					cases = cases + 1
					o1 = o[rownames(g)[loable[e-1]],]
					o2 = o[rownames(g)[loable[e  ]],]
					candidates = vector('list',ncol(o))
					names(candidates) = colnames(o)
					for(sp in colnames(o)){
						gid = unique(gsub('.s\\d+','',c(o1[sp],o2[sp])))
						if(length(gid) != 1)
							break
						gs = a[[sp]][[gid]]
						candidates[[sp]] = gs[gs$exon.number > gs[o1[sp],'exon.number'] & gs$exon.number < gs[o2[sp],'exon.number'],]
						
						if(sum(candidates[[sp]]$has.orth) > 0){ #skip if there is unexpected orth exon
							candidates[[sp]] = 	list(NULL)
							break
						}
					}
					if(sum(sapply(candidates,is.null))==0 && length(unique(sapply(candidates,nrow)))==1){
						res[[paste0('t',length(res))]] = candidates
					}else{
						bad = bad + 1
						if(sum(sapply(candidates,is.null))==0)
							bad.data[[paste0('t',length(bad.data))]] = candidates
					}
				}
			}
		}
	}
	cat('considered ',bad + length(res),' cases. ',length(res),' where ok\n',sep='')
	if(returnBad)
		return(as.list(bad.data))
	res = as.list(res)
	clnm = colnames(res[[1]][[length(res[[1]])]])
	res = sapply(res,function(x){
		lapply(1:nrow(x[[1]]),function(i){do.call(rbind,lapply(x,function(z){cbind(seg.id=rownames(z)[i],z[i,clnm],group.size=nrow(z))}))})
	})
	res = unlist(res,recursive = FALSE)
	do.call(rbind,lapply(res,function(x){x=data.frame(c(as.list(x$seg.id),range(x$length),length(unique(x$length%%3))==1,x$group.size[1]));colnames(x) = c(names(a),'length.min','length.max','same.frame','group.size');x}))
}


getAllPairsBySinglMate = function(pairs,o,all.pairs,anns){
	pairs$used = FALSE
	r = list()
	for(i in 1:nrow(pairs)){
		if(pairs$used[i]) next;
		t = getPairsBySinglMate(pairs[i,],o,all.pairs,anns)
		pairs$used[rownames(pairs) %in% paste(t[,1],t[,2])] = TRUE
		r[[length(r)+1]] = t
	}
	r
}
getPairsBySinglMate = function(pair,o,pairs,anns){
	inx = 1
	osids = which(apply(o==pair[1,1],1,sum)>0)
	if(length(osids)==0){
		osids = which(apply(o==pair[1,2],1,sum)>0)
		inx = 2
	}
	osids = o[osids,]
	r = do.call(rbind,lapply(1:length(pairs),function(x){pairs[[x]][pairs[[x]][,inx] == osids[x] & pairs[[x]]$exn.no.diff == pair$exn.no.diff,]}))
	rownames(r) = substr(r$seg.id1,1,3)
	rr = NULL
	for(i in 1:nrow(species)){
		s = substr(rownames(species)[i],1,3)
		if(s %in% rownames(r))
			rr = rbind(rr,r[s,])
		else{
			seg1 = anns[[i]][osids[i],]
			seg1$seg.id = rownames(seg1)
			seg2 = anns[[i]][anns[[i]]$gene_id == seg1$gene_id &  anns[[i]]$sites == 'ad' & anns[[i]]$exon.number == (seg1$exon.number + pair$exn.no.diff*ifelse(inx==2,-1,1)),]
			if(nrow(seg2)==0){
				seg2 = list(seg.id=NA,length=NA)
			}else
				seg2$seg.id = rownames(seg2)
			if(inx == 2){
				t = seg1
				seg1 = seg2
				seg2 = t
			}
			rr = rbind(rr,data.frame(
				seg.id1 = seg1$seg.id,
				seg.id2 = seg2$seg.id,
				cor = NA,exn.no.diff=pair$exn.no.diff, psi.sum.mean=NA, psi.sum.sd=NA, psi.sum.min=NA, psi.sum.max=NA, len1=seg1$length, len2=seg2$length)
			)
		}
	}
	rownames(rr) = rownames(species)
	rr
}

#' @param sgn is not null, only pairs where both segments are from sgn will be considered
getOrthMEXs = function(cors,orths,cor.thr,mean.psi.sum.dev,sgn=NULL){
	orth.mexs = list()
	no.orth = NULL
	cors = lapply(cors,function(x){x=cbind(x,used=FALSE);rownames(x)=paste(x$seg.id1,x$seg.id2);x})
	for(s in names(cors)){
		print(s)
		sp.mexs = cors[[s]][!is.na(cors[[s]]$cor) & cors[[s]]$cor <= cor.thr & abs(cors[[s]]$psi.sum.mean - 1) <= mean.psi.sum.dev & !cors[[s]]$used,]
		if(!is.null(sgn))
			sp.mexs = sp.mexs[sp.mexs$seg.id1 %in% sgn & sp.mexs$seg.id2 %in% sgn,]
		for(m in 1:nrow(sp.mexs)){
			has.orth = (sp.mexs$seg.id1[m] %in% orths[,s]) + (sp.mexs$seg.id2[m] %in% orths[,s])
			if(has.orth == 2){
				id1 = which(orths[,s] == sp.mexs$seg.id1[m])
				id2 = which(orths[,s] == sp.mexs$seg.id2[m])
				ids = orths[c(id1,id2),]
				same.genes = sum(!apply(ids,2,function(x){x=sapply(strsplit(x,'.',TRUE),'[',1);x[1]==x[2]}))==0
				if(same.genes){
					inxs = setNames(rep(NA,length(cors)),names(cors))
					for(s1 in names(cors)){
						idss = ids[,s1]
						#idss = idss[order(anns[[s1]][idss,'start'])]
						inxs[s1] = paste(idss,collapse=' ')
						if(is.na(cors[[s1]][inxs[s1],][,1]))
							inxs[s1] = NA
					}
					if(sum(is.na(inxs))==0){
						mex = NULL
						for(s1 in names(cors)){
							mex = rbind(mex,cors[[s1]][inxs[s1],])
							cors[[s1]][inxs[s1],'used'] = TRUE
						}
						orth.mexs[[length(orth.mexs)+1]] = mex
					}else
						no.orth = rbind(no.orth,cbind(sp.mexs[m,],has.orth=-1)) # orthologs are from same gene, but do not passed threshold in some of species
				}else
					no.orth = rbind(no.orth,cbind(sp.mexs[m,],has.orth=-2)) #has orthologs but they are in different genes
			}else{
				no.orth = rbind(no.orth,cbind(sp.mexs[m,],has.orth=has.orth))
			}
		}
	}
	list(orth=orth.mexs,no.orth=no.orth)
}

checkGeneClusterOverlapOnClDist = function(orth.cl,s,c1,c2,cl.dists){
	r = matrix(NA,ncol=2,nrow=length(cl.dists)-1,dimnames=list(paste(cl.dists[-length(cl.dists)],'-',cl.dists[-1],sep=''),c('estimate','p.value')))
	for(i in 1:(length(cl.dists)-1)){
		r[i,] = unlist(checkGeneClusterOverlap(orth.cl,s,c1,c2,'test',cl.dists[i+1],cl.dists[i])[c('estimate','p.value')])
	}
	r
}

plotGeneAllClustersOverlap = function(orth.cl,cl.max.dist){
	rr = NULL
	pv.col=c(gray=2,yellow=5e-2,orange=1e-4,red=1e-10,none=0)
	or.col=c(blue=0,'#0000FFAA'=1/3,'#0000FFAA'=1/1.5,'#0000FF66'=1,yellow=1.000000001,'#FF000066'=1.5,'#FF0000AA'=3,red=Inf)
	for(s in rownames(species)){
		r = checkGeneAllClustersOverlap(orth.cl,s,cl.max.dist = cl.max.dist)
		plotFTMotifSimMatrix(r,main=s,pv.col=pv.col,or.col=or.col,pv.thr = pv.col[2])
		abline(v=(1:7)*2)
		abline(h=28-(0:7)*2)
		r[is.infinite(r)] = 100
		r[r==0] = 0.01
		if(is.null(rr)) 
			rr = r
		else
			rr = rr * r
	}
	plotFTMotifSimMatrix(rr^(1/7),main='mean',pv.col=pv.col,or.col=or.col,pv.thr = pv.col[2])
	abline(v=(1:7)*2)
	abline(h=28-(0:7)*2)
}


plotSameGenePsiCor = function(c,sh){
	for(s in names(same.gene.psi.cor)){
		q = c[[s]]
		boxplot(log2(q$exn.no.diff) ~ I(floor(q$cor*20)/20),main=s,yaxt='n',ylab='Distance (segments)',xlab='Correlation')
		axis(2,1:10,2^(1:10))
		boxplot(q$psi.sum.mean ~ I(floor(q$cor*20)/20),main=s,ylab='Mean(PSI sum)',xlab='Correlation')
		boxplot(q$psi.sum.sd ~ I(floor(q$cor*20)/20),main=s,ylab='SD(PSI sum)',xlab='Correlation')
		
		hist(sh[[s]]$cor,-40:40/40,col='gray',border=NA,xlab='Correlation',main=s)
		hist(q$cor,-40:40/40,add=T,col='red',border=NA,den=40)
		legend('topleft',fil=c('gray','red'),den=c(-1,40),legend=c('shuffled','observed'))
	}
}

addExonNumber = function(segs){
	print('\taddExonNumber')
	require(data.table)
	ss = split(segs,segs$gene_id)
	k = 1
	r = lapply(ss,function(g){
		cat('\r',k,length(ss))
		k <<- k + 1
		g = g[order(g$start*g$strand),]
		g$exon.number = NA
		n = 1
		exn = TRUE
		for(i in 1:nrow(g)){
			if(exn){
				if(substr(g$sites[i],2,2) %in% c('d','.'))
					exn = FALSE
			}else{
				if(substr(g$sites[i],1,1) %in% c('a','.')){
					exn = TRUE
					n = n + 1
					if(substr(g$sites[i],2,2) %in% c('d','.'))
						exn = FALSE
				}
			}
			g$exon.number[i] = n
			if(substr(g$sites[i],1,1) %in% c('a','.') & g$position[i] == 'LAST'){ # to account for internal lasts
				n = n - 1
				exn = FALSE
			}
		}
		g$seg.id = rownames(g)
		g
	})
	r = data.frame(rbindlist(r))
	rownames(r) = r$seg.id
	r$seg.id = NULL
	r[rownames(segs),]
}

calcWithinGenePSICor = function(psi,seg,ids=NULL,shuffle=FALSE){
	if(!is.null(ids))
		seg = seg[ids,]
	psi = psi[rownames(seg),]
	
	if(shuffle)
		seg$gene_id = sample(seg$gene_id)
	psi = split.data.frame(psi,seg$gene_id)
	k=0
	r=lapply(psi,function(x){
		k <<- k + 1
		if(nrow(x)==1) return(NULL)
		cat('\r',k,length(psi),'       ')
		s = seg[rownames(x),]
		x = x[order(s$start*s$strand),]
		l = nrow(x)*(nrow(x)-1)/2
		r = data.frame(seg.id1 = rep(NA,l),seg.id2 = rep(NA,l),cor = rep(NA,l),exn.no.diff= rep(NA,l),psi.sum.mean= rep(NA,l),psi.sum.sd= rep(NA,l),psi.sum.min= rep(NA,l),psi.sum.max= rep(NA,l))
		k = 1
		for(i in 1:(nrow(x)-1))
			for(j in (i+1):nrow(x)){
				r[k,] = list(rownames(x)[i],rownames(x)[j],cor(x[i,],x[j,],u='pair'),abs(seg[rownames(x)[i],'exon.number']-seg[rownames(x)[j],'exon.number']),mean(x[i,] + x[j,],na.rm=TRUE),sd(x[i,] + x[j,],na.rm=TRUE),min(x[i,] + x[j,],na.rm=TRUE),max(x[i,] + x[j,],na.rm=TRUE))
				k = k + 1
			}
		r
	})
	require(data.table)
	p =  data.frame(rbindlist(r))
	p
}

checkGeneAllClustersOverlap = function(orth.cl,species,cl.max.dist=Inf,cl.min.dist=-Inf){
	clno = max(orth.cl$orth.cl)
	r = matrix(NA,ncol=clno,nrow=clno,dimnames = list(paste0('c',1:clno),paste0('c',1:clno)))
	for(c1 in 1:(clno-1))
		for(c2 in (c1+1):clno){
			t = checkGeneClusterOverlap(orth.cl,species,c1,c2,'test',cl.max.dist,cl.min.dist)
			r[c1,c2] = t$p.value
			r[c2,c1] = t$estimate
		}
	r
}

plotOrthClustCodProp = function(s){
	c = orth.ad.cl.reassign[orth.ad.cl.reassign$species==s,]
	cod = paste(anns[[s]][rownames(c),'cod'],ifelse(anns[[s]][rownames(c),'cod.gene'],'c','n'),sep='')
	t = table(cod,c$orth.cl)
	not.sgn = rownames(anns[[s]]) %in% rownames(c)
	not.sgn = paste(anns[[s]][not.sgn,'cod'],ifelse(anns[[s]][not.sgn,'cod.gene'],'c','n'),sep='')
	t = cbind('!s'=table(not.sgn)[rownames(t)],t)
	t = sweep(t,2,apply(t,2,sum),'/')
	t = t[c('cc','pc','nc','nn'),]
	barplot(t,col=c('red','red','red','blue'),den=c(-1,60,30,30),xlab='Cluster',ylab='% of exons',main=s)
}


r.comp = function(s){
	sapply(s,r.comp.)
}

r.comp. = function(s){
	if(s == '')
		return('')
	t = rev(comp(s2c(s)))
	t[is.na(t)] = 'n'
	c2s(t)
}



plotSpSpecificCounts = function(obs,exp,tissue.cols,tinxs=1:dim(obs)[4],main='',thrs=c(1,0.1)){
	r2 = r2.sh = r1 = r1.sh = NULL
	thrs = c(1,0.1)#c(1,0.1)
	for(t in unique(meta$tissue)){
		r1    = rbind(r1   ,table(annotateSpSpec(obs,t,1,thrs = thrs)))
		r1.sh = rbind(r1.sh,table(annotateSpSpec(exp,t,1,thrs = thrs)))
		r2    = rbind(r2   ,table(annotateSpSpec(obs,t,2,thrs = thrs)))
		r2.sh = rbind(r2.sh,table(annotateSpSpec(exp,t,2,thrs = thrs)))
	}
	barplot(r1[tinxs,-1],col=tissue.cols[unique(meta$tissue)[tinxs]],beside = TRUE,xlab='Species',ylab='# of segments')
	barplot(r1.sh[tinxs,-1],col='gray',beside = TRUE,den=30,add=TRUE,border=NA)
	
	inxs = c(2,13,22,3:12,14:21)
	barplot(r2[tinxs,inxs],col=tissue.cols[unique(meta$tissue)[tinxs]],beside = TRUE,xlab='Species pairs',ylab='# of segments')
	barplot(r2.sh[tinxs,inxs],col='gray',beside = TRUE,den=30,add=TRUE,border=NA)
	mtext(main,3,outer = TRUE)
}


plotOneSpecies2DvarDens = function(obs,exp,main='',thr=c(1,0.1)){
	to1 = te1 = NULL
	for(t in unique(meta$tissue)){
		for(s1 in 1:nrow(species)){
			o. =  getSpSpecSegsFromDivMatrix(obs,rownames(species)[s1],rownames(species)[-s1],t)
			e. =  getSpSpecSegsFromDivMatrix(exp,rownames(species)[s1],rownames(species)[-s1],t)
			o =  get2dFreq(o.,n=100)
			e =  get2dFreq(e.,n=100)
			o.cnt = sum(o.[,2] > thr[1]*o.[,1] & o.[,2] > thr[2],na.rm=T)
			e.cnt = sum(e.[,2] > thr[1]*e.[,1] & e.[,2] > thr[2],na.rm=T)
			if(is.null(to1)){
				to1 = o
				te1 = e
			}else{
				to1 = to1 + o
				te1 = te1 + e
			}
			
			plotPattDivQV(o,e,T,xlab='max internal distance',ylab='min external distance',main=paste(t,' ',species$short[s1],' (',o.cnt,' ',round(e.cnt/o.cnt*100,1),'%)',sep=''))
			#points(e.[,1]-0.005,e.[,2]+0.005,pch=16,col='blue',cex=1)
			#points(o.[,1]-0.005,o.[,2]+0.005,pch=16,col='magenta',cex=0.5)
		}
	}
	plotPattDivQV(to1,te1,T,xlab='max internal distance',ylab='min external distance',main='All together')
	mtext(main,3,outer = TRUE)
}

plotTwoSpecies2DvarDens = function(obs,exp,main='',thr=c(1,0.1)){
	sps = combn(rownames(species),2,simpl=F)
	names(sps) = sapply(sps,function(x){paste(species[x,'short'],collapse='')})
	pairs = c('hq','mr','oc')
	for(t in unique(meta$tissue)){
		o = lapply(sps,function(s){getSpSpecSegsFromDivMatrix(obs,s,setdiff(rownames(species),s),t)})
		e = lapply(sps,function(s){getSpSpecSegsFromDivMatrix(exp,s,setdiff(rownames(species),s),t)})
		o = c(o[pairs],list(exp=do.call(rbind,o[setdiff(names(o),pairs)])))
		for(p in pairs){
			o.cnt = sum(o[[p]][,2] > thr[1]*o[[p]][,1] & o[[p]][,2] > thr[2],na.rm=T)
			e.cnt = sum(e[[p]][,2] > thr[1]*e[[p]][,1] & e[[p]][,2] > thr[2],na.rm=T)
			ot =  get2dFreq(o[[p]],n=100)
			et =  get2dFreq(e[[p]],n=100)
			plotPattDivQV(ot,et,T,xlab='max internal distance',ylab='min external distance',main=paste('Comp to sh; ',t,' ',p,' (',o.cnt,' ',round(e.cnt/o.cnt*100,1),'%)',sep=''))
		}
		for(p in pairs){
			o.cnt = sum(o[[p]][,2] > thr[1]*o[[p]][,1] & o[[p]][,2] > thr[2],na.rm=T)
			e.cnt = sum(o$exp[,2] > thr[1]*o$exp[,1] & o$exp[,2] > thr[2],na.rm=T)/(length(sps)-length(pairs))
			ot =  get2dFreq(o[[p]],n=100)
			et =  get2dFreq(o$exp,n=100)
			plotPattDivQV(ot,et,T,xlab='max internal distance',ylab='min external distance',main=paste('Comp to pairs; ',t,' ',p,' (',o.cnt,' ',round(e.cnt/o.cnt*100,1),'%)',sep=''))
		}
	}
	mtext(main,3,outer = TRUE)
}

annotateSpSpec = function(d,t,sp.no,thrs=c(1,0.1),na.rm=FALSE){
	sps = combn(rownames(species),sp.no,simpl=F)
	r = rep('',dim(d)[1])
	names(r) = dimnames(d)[[1]]
	for(sp in sps){
		tmp = getSpSpecSegsFromDivMatrix(d,sp,setdiff(rownames(species),sp),t,na.rm = na.rm)
		f = tmp[,2] >= tmp[,1]*thrs[1] & tmp[,2] >= thrs[2]
		f[is.na(f)] = FALSE
		if(sum(r[f]!='')>0) browser()
		r[f] = paste(species[sp,'short'],collapse = '')
	}
	
	factor(r,levels=c('',sapply(sps,function(x){paste(species[x,'short'],collapse = '')})))
}

getCorOnAge = function(c,s1,s2,ts,as){
	r = matrix(NA, ncol=length(as),nrow=length(ts),dimnames = list(ts,as))
	for(t in ts)
		for(a in as){
			i1 = paste(s1,t,a)
			i2 = paste(s2,t,a)
			if(sum(c(i1,i2) %in% colnames(c))==2)
				r[t,as.character(a)] = c[i1,i2]
		}
	colnames(r) = round(as.numeric(colnames(r)),5)
	r
}

getCorOnSpecies = function(c,s1,a,ss,ts){
	r = matrix(NA, ncol=length(ss),nrow=length(ts),dimnames = list(ts,ss))
	for(t in ts)
		for(s2 in ss){
			i1 = paste(s1,t,a)
			i2 = paste(s2,t,a)
			if(sum(c(i1,i2) %in% colnames(c))==2)
				r[t,s2] = c[i1,i2]
		}
	r
}

getAltExonStat = function(psi,m,psi.thr,tissues = unique(m$tissue),summ=T,na.as.cnst=TRUE){
	m = m[rownames(m) %in% colnames(psi) & m$tissue %in% tissues ,]
	m = m[order(m$days,m$tissue),]
	psi = psi[,rownames(m)]
	stages = unique(m[,c('stage','days')])
	stages = stages[order(stages$days),]
	res = vector('list',length = nrow(stages))
	for(i in 1:nrow(stages)){
		p = psi[,m$stage==stages$stage[i],drop=F]
		if(!na.as.cnst)
			p = p[apply(is.na(p),1,sum)==0,]
		p = !is.na(p) & p > psi.thr & p < (1-psi.thr)
		ts = substr(m[colnames(p),'tissue'],1,1)
		res[[i]] = apply(p,1,function(x)paste(ts[x],collapse=''))
	}
	names(res) = stages$days
	getAltExonStat2table(res,summ)
}

getAltExonStat2table = function(res,summ){
	all = unique(unlist(res))
	r = sapply(res,function(x)table(factor(x,levels = all)))
	r = r[order(nchar(rownames(r)),rownames(r)),]
	if(summ){
		bh = r['bh',]
		r = r[rownames(r) != 'bh',]
		r. = NULL
		for(i in 2:max(nchar(rownames(r)))){
			r. = rbind(r.,apply(r[nchar(rownames(r)) == i,,drop=F],2,sum))
		}
		rownames(r.) = paste(2:max(nchar(rownames(r))),'ts',sep='')
		f = nchar(rownames(r)) ==0 | rownames(r) %in% c('b','h')
		r = rbind(r[f,],bh,r[!f & nchar(rownames(r)) ==1,],r.)
	}
	
	rownames(r)[rownames(r)==''] = 'const'
	colnames(r) = names(res)
	r
}


plotStageTissueMatrix = function(m,species){
	t = table(m$tissue[m$species==species],round(m$days[m$species==species],1))
	t = t[order(rownames(t),decreasing = T),]
	image(x=1:ncol(t),y=1:nrow(t),t(t),col=c('white','black'),xaxt='n',yaxt='n',xlab='Stages (days)',ylab='',main=species,breaks=c(0,0.5,1e10))
	axis(1,1:ncol(t),colnames(t))
	axis(2,1:nrow(t),rownames(t),las=2)
	invisible(t)
}


plotStageMatrix = function(m,tissue){
	m$a = factor(round(exp(m$age.use),1))
	t = table(m$species[m$tissue==tissue],m$a[m$tissue==tissue])
	image(x=1:ncol(t),y=1:nrow(t),t(t),col=c('white','black'),xaxt='n',yaxt='n',xlab='Stages',ylab='',main=tissue)
	axis(1,1:ncol(t),colnames(t))
	axis(2,1:nrow(t),rownames(t),las=2)
	invisible(t)
}

plotPattDivQV = function(o,e,cum=TRUE,col=c('gray','red','orange','yellow','white'),breaks = c(-1,-1e-20,0.00001,0.2,0.9,1e50),xlim=NULL,ylim=NULL,...){
	if(cum){#qv is calculated for each in-veriability bin
		o = t(apply(o,1,function(x){x=rev(cumsum(rev(x)));x/max(1,x)}))
		e = t(apply(e,1,function(x){x=rev(cumsum(rev(x)));x/max(1,x)}))
	}else{
		o = o/sum(o)
		e = e/sum(e)
	}
	q = e/o
	q[o==0] = 1
	q[o == 0 & e == 0] = -1
	if(is.null(xlim))
		xlim=c(-1,1+max(1,which(apply(o+e>0,1,sum)>0)))/nrow(o)
	if(is.null(ylim))
		ylim=c(-1,1+max(1,which(apply(o+e>0,2,sum)>0)))/nrow(o)
	image(q,col=col,breaks = breaks,xlim=xlim,ylim=ylim,xaxt='n',yaxt='n',...)
	axis(1,0:20/20,0:20/10)
	axis(2,0:20/20,0:20/10)
	abline(a=0,b=1)
	abline(h=0.05)
	invisible(q)
}


get2dFreq = function(x,y=NULL,n=20,rangex=c(0,2),rangey=c(0,2)){
	if(is.matrix(x)){
		y = x[,2]
		x = x[,1]
	}
	x = floor((x - rangex[1])/(rangex[2] - rangex[1])*n + 1)
	y = floor((y - rangey[1])/(rangey[2] - rangey[1])*n + 1)
	r = matrix(0,ncol=n,nrow=n)
	for(i in 1:length(x))
		r[x[i],y[i]] = r[x[i],y[i]] + 1
	r
}

getSpSpecSegsFromDivMatrix = function(m,ss1,ss2,tissue,na.rm=FALSE){
	r = cbind(rep(-Inf,dim(m)[1]),rep(Inf,dim(m)[1]))
	not.na.cnt = rep(0,dim(m)[1])
	for(s1 in ss1)
		for(s2 in ss2){
			r[,2] = pmin(r[,2] , m[,s1,s2,tissue],na.rm = na.rm)
			not.na.cnt = not.na.cnt + !is.na(m[,s1,s2,tissue])
		}
	if(na.rm) #to make classification unambiguous
		r[not.na.cnt <= min(length(ss1),length(ss2)),2] = NA
	for(s1 in ss1)
		for(s2 in setdiff(ss1,s1))
			r[,1] = pmax(r[,1] , m[,s1,s2,tissue],na.rm = na.rm)
	
	for(s1 in ss2)
		for(s2 in setdiff(ss2,s1))
			r[,1] = pmax(r[,1] , m[,s1,s2,tissue],na.rm = na.rm)
	r[is.infinite(r)] = NA
	r
}


#' Title
#'
#' @param psi table for all species
#' @param m sample information
#' @param substr.species should species means be substracted. Means are calculated based only on tissue/stages present in all species
#' @param spline.df number of df to be used with spline. set to NULL to use raw data
#' @param spline.n number of points to be used for approximation
#' @param min.obs results will be NA, if number of not NAs in given tissue-species is lower than min.obs
#' @param shuffle should raw data be shuffled (all not NA values will be shuffled per stage (within given tissue))
#' @param min.cmn.obs minimal allowed number of tissue/stages present in all species. Only relevant if substr.species is true.
#'
#' @return 4-way array: segment*species1*species2*tissue
#'
#' @examples
caclPatternDiv = function(psi,m,substr.species,spline.df=NULL,spline.n=500,min.obs=0,shuffle=FALSE,min.cmn.obs=5){
	sp = rownames(species)
	ts = unique(m$tissue)
	ages.d = sort(unique(m$mouse.days))
	ages = as.character(ages.d)
	psi = psi[,rownames(m),drop=FALSE]
	r = array(NA,c(nrow(psi),length(sp),length(sp),length(ts)),list(rownames(psi),sp,sp,ts))
	all.tissue.stages = paste(m$tissue,m$age.use)
	for(seg in 1:nrow(psi)){
		cat('\r',seg)
		p = psi[seg,]
		if(substr.species){ #calculate per-species mean but only for tissue-stages that are present in all species
			ats = all.tissue.stages
			for(s in sp)
				ats = intersect(ats,all.tissue.stages[m$species == s & !is.na(p)])
			if(length(ats) >= min.cmn.obs){
				for(s in sp)
					p[m$species == s] = p[m$species == s] - mean(p[m$species == s & all.tissue.stages %in% ats],na.rm=TRUE)
			}else
				next;
		}
		for(t in ts){
			#prepare PSIs to be arranged in the same age scale
			sp.psi = matrix(NA,nrow=length(sp),ncol=length(ages))
			for(i in 1:length(sp)){
				sp.psi[i,] = p[paste(sp[i],t,ages)]
				if(sum(!is.na(sp.psi[[i]])) < min.obs)
					sp.psi[i,] = NA
			}
			#shuffle
			if(shuffle){
				for(c in 1:ncol(sp.psi))
					sp.psi[!is.na(sp.psi[,c]),c] = sample(sp.psi[!is.na(sp.psi[,c]),c])
			}
			#approximate PSI
			if(!is.null(spline.df)){
				sp.psi. = matrix(NA,nrow=length(sp),ncol=spline.n)
				for(i in 1:length(sp)){
					if(!is.null(sp.psi[[i]])){
						a =      ages.d[!is.na(sp.psi[i,])]
						y =  sp.psi[i,][!is.na(sp.psi[i,])]
						if(length(a) >= max(spline.df,4)){
							pred = predict(smooth.spline(a,y,df=spline.df),seq(ages.d[1],ages.d[length(ages.d)],length.out = spline.n))
							pred$y[pred$x < a[1] | pred$x > a[length(a)]] = NA
							pred$y = pmin(1,pmax(-1,pred$y))
							pred$y = pmin(1,pmax(-1,pred$y))
							sp.psi.[i,] = pred$y
						}
					}
				}
				sp.psi = sp.psi.
			}
			for(s1 in 1:(length(sp)-1)){
				for(s2 in (s1+1):length(sp)){
					r[seg,s1,s2,t] = r[seg,s2,s1,t] = mean(abs(sp.psi[s1,] - sp.psi[s2,]),na.rm=TRUE)
				}
			}
		}
	}
	r
}


heatmap3.my = function(d,m,cols,f=rep(TRUE,nrow(d)),dist=FALSE,...){
	m = m[f,]
	d = d[rownames(m),rownames(m)]
	colsleg = NULL
	for(n in names(cols)){
		cols[[n]] = cols[[n]][names(cols[[n]]) %in% as.character(m[[n]])]
		colsleg = cbind(colsleg,cols[[n]][as.character(m[[n]])])
	}
	colnames(colsleg) = names(cols)
	col = getPal()
	if(dist) col = rev(col)
	heatmap3(d,col=col,showRowDendro=F,symm=TRUE,distfun = function(x){if(!dist){x=1-x};as.dist(x)},scale = 'none',labRow = NA,labCol = NA,RowSideColors = colsleg,...)
	legend(grconvertX(0,'ndc'),grconvertY(0.7,'ndc'),xpd=TRUE,fill=c(cols$species,cols$tissue),legend=c(names(cols$species),names(cols$tissue)),ncol = 1,title='species/tissues')
}

calcPatternDistance = function(ir,m,substr.mean=FALSE){
	st = unique(m[,c('species','tissue')])
	st = st[order(st$species,st$tissue),]
	r = matrix(NA,ncol=nrow(st),nrow=nrow(st))
	for(i in 1:nrow(st)){
		r[i,i] = 0
		if(i < nrow(st))
			for(j in (i+1):nrow(st)){
				r[i,j] = r[j,i] = calcDistanceBetweenSpeciesTissue(ir,m,st$species[i],st$species[j],st$tissue[i],st$tissue[j],substr.mean)
			}
	}
	colnames(r)=rownames(r) = paste(st$species,st$tissue)
	r
}

calcDistanceBetweenSpeciesTissue = function(ir,m,s1,s2,t1,t2,substr.mean){
	common.stages = intersect(m$age.use[m$species==s1 & m$tissue==t1],m$age.use[m$species==s2 & m$tissue==t2])
	ids1 = rownames(m)[m$species==s1 & m$tissue==t1 & m$age.use %in% common.stages]
	ids2 = rownames(m)[m$species==s2 & m$tissue==t2 & m$age.use %in% common.stages]
	ids1 = ids1[order(m[ids1,'age.use'])]
	ids2 = ids2[order(m[ids2,'age.use'])]
	ir1 = ir[,ids1,drop=F]
	ir2 = ir[,ids2,drop=F]
	if(substr.mean){
		ir1 = sweep(ir1,1,apply(ir1,1,mean,na.rm=TRUE))
		ir2 = sweep(ir2,1,apply(ir2,1,mean,na.rm=TRUE))
	}
	#sqrt(mean((ir1-ir2)^2,na.rm=T))
	#cor(as.numeric(ir1),as.numeric(ir2),u='p')
	t = apply(abs(ir1-ir2),1,max,na.rm=T)
	t = t[!is.infinite(t)]
	mean(t,na.rm=T)
}

testConsGOEnrich = function(go,ont){
	cmn = go[[1]][[1]]
	cmn = cmn$category[cmn$ontology == ont]
	for(s in 2:length(go)){
		cmn = intersect(cmn,go[[s]][[1]]$category[go[[s]][[1]]$ontology == ont])
	}
	ranks = list()
	for(s in names(go)){
		ranks[[s]] = matrix(NA,ncol=length(go[[1]]),nrow=length(cmn))
		rownames(ranks[[s]]) = cmn
		for(c in 1:length(go[[1]])){
			rownames(go[[s]][[c]]) = go[[s]][[c]][,'category']
			ranks[[s]][,c] = rank(go[[s]][[c]][cmn,'pv.over'],tie='random')/length(cmn)
		}
	}
	applyIrwin.hall(ranks)
}

getGO = function(sel,bkg,min.go.size,method = c('Hypergeometric',"Wallenius",'Sampling')[1],genome='hg19',id='ensGene',gene2cat=NULL) {
	if(is.character(bkg))
		bkg = setNames(rep(0,length(bkg)),bkg)
	sel.genes = names(bkg) %in% sel
	names(sel.genes) = names(bkg)
	if(method != "Hypergeometric")
		pwf=nullp(DEgenes=sel.genes,bias.data=bkg,plot.fit=T)
	else
		pwf = data.frame(DEgenes=sel.genes,row.names=names(bkg),bias.data=1,pwf=1)
	suppressMessages(go <- goseq(pwf,genome=genome,id=id,method=method,gene2cat=gene2cat))
	colnames(go) = c('category','pv.over','pv.under','sel.size','size','term','ontology')
	go = go[!is.na(go$ontology) & go$size >= min.go.size,] # to remove deprecated terms
	terms = unique(go$ontology)
	for(t in terms){
		f = go$ontology == t
		go$qv.over[f]  = p.adjust(go$pv.over [f],method='BH')
		go$qv.under[f] = p.adjust(go$pv.under[f],method='BH')
	}
	go = go[,c('category','pv.over','pv.under','qv.over','qv.under','ontology','term','size','sel.size')]
	go = go[order(go$pv.over,go$pv.under),]
	rownames(go) = go$category
	go
}


loadSpeciesGoAnn = function(f,go.terms){
	go.ann = read.csv(f)
	go.ann[,2] = go.terms[go.ann[,2],'current.go.id']
	go.ann = split(go.ann[,2],go.ann[,1])
	go.ann = lapply(go.ann,function(x){unique(x[x!=''])})
	go.ann
}


makeSeg2Ens = function(s){
	seg = readRDS(paste('Rdata/',s,'.as.u.all.Rdata',sep=''))$seg
	comp = read.table(list.files('processed/annotation/all.species/merged/',paste(s,'.*comp',sep=''),full.names=TRUE))
	comp = comp[comp[,3] %in% c('=','c','j','e'),]
	ens = read.table(list.files('processed/annotation/all.species/ensambl/',paste(gsub(' ','_',species[s,'name'],fixed = TRUE),'.*gtf.gz',sep=''),full.names=TRUE),sep='\t')
	#prepare comp
	comp = lapply(split(comp,comp$V1),function(x){
		for(c in list('=',c('c','j'),'e')){
			if(sum(x$V3 %in% c)>0)	return(x[x$V3 %in% c,])
		}
		NULL})
	#prepare ens
	ens = ens[ens$V3=='exon',]
	gids = gsub(';','',gsub('gene_id ','',my.regmatches(ens$V9,'gene_id [^;]+;')))
	ens.gene.coors = t(sapply(split(ens[,4:5],gids),function(x){c(min(x),max(x))}))
	#annotate segments
	s2e = comp[seg$gene_id]
	s2e = lapply(1:length(s2e),function(i){
		if(is.null(s2e[[i]])) return(character(0))
		if(nrow(s2e[[i]])==1) return(s2e[[i]]$V2)
		s = seg[i,]
		dist = apply(ens.gene.coors[s2e[[i]]$V2,],1,function(g){min(g[2],s$stop)-max(g[1],s$start)})
		s2e[[i]]$V2[dist == max(dist)]
	})
	names(s2e) = rownames(seg)
	s2e
}

getSpSpMAtrixFromList = function(x,i){
	r = array(NA,dim = c(nrow(species),nrow(species),nrow(x[[1]])),dimnames = list(species$short,species$short,rownames(x[[1]])))
	for(s1 in 1:(nrow(species)-1))
		for(s2 in (s1+1):nrow(species))
			r[species$short[s1],species$short[s2],] = r[species$short[s2],species$short[s1],] = x[[paste(species$short[s1],species$short[s2],sep='')]][,i]
		r
}

getAllSpecificSpeciesGroups = function(m,clust = NULL){
	if(!is.null(clust)){
		colnames(clust) = species[colnames(clust),'short']
		clust = clust[,dimnames(m)[[1]]]
		clust.m = m
		for(s in 1:nrow(clust)){
			for(p in 1:ncol(clust))
				clust.m[p,,s] = clust[s,] != clust[s,p]
		}
		m = m & clust.m
	}
	r = vector('list',dim(m)[3])
	names(r) = dimnames(m)[[3]]
	spgroups = list()
	used.groups = c()
	for(i in 1:floor((dim(m)[1]/2)))
		spgroups = c(spgroups,combn(dimnames(m)[[1]],i,simplify = FALSE))
	for(g in 1:length(spgroups)){
		g1 = spgroups[[g]]
		cat(g,paste(g1,collapse=''),length(spgroups),'\n')
		used.groups[length(used.groups)+1]=paste(g1,collapse='')
		g2 = setdiff(dimnames(m)[[1]],g1)
		if(paste(g2,collapse='') %in% used.groups) next
		for(s in 1:dim(m)[3]){
			if(sum(m[g1,g2,s],na.rm=TRUE)==length(g1)*length(g2)){
				if(length(r[[s]])==0)
					r[[s]] = list()
				if(sum(g1 %in% unlist(r[[s]]))<length(g1))
					r[[s]][[length(r[[s]])+1]] = g1
			}
		}
	}
	sapply(r,function(x){paste(sapply(x,paste,collapse=''),collapse=',')})
}

classifySpSPMatrix = function(m,max.outmax,clust = NULL,cor.sgn=NULL){
	pv.thr = max(m[p.adjust(m,m='BH')<0.05],na.rm=T) #used only to mark unclassified segments
	if(!is.null(clust)){
		colnames(clust) = species[colnames(clust),'short']
		clust = clust[,dimnames(m)[[1]]]
		clust.m = m
		for(s in 1:nrow(clust)){
			for(p in 1:ncol(clust))
				clust.m[p,,s] = clust[s,] != clust[s,p]
		}
		if(is.null(cor.sgn))
			sgn.cnt = apply(m<pv.thr & clust.m,3,sum,na.rm=TRUE)/2
		else
			sgn.cnt = apply(m<pv.thr & clust.m & cor.sgn,3,sum,na.rm=TRUE)/2
	}else{
		if(is.null(cor.sgn))
			sgn.cnt = apply(m<pv.thr,3,sum,na.rm=TRUE)/2
		else
			sgn.cnt = apply(m<pv.thr & cor.sgn,3,sum,na.rm=TRUE)/2
	}
	r = data.frame(spec=character(dim(m)[3]),inmin=NA,outmax=NA)
	rownames(r) = dimnames(m)[[3]]
	r$spec = 'n'
	n=dim(m)[1]
	r$spec[sgn.cnt > 0 & sgn.cnt < n-1] = 's5' #too few sign
	nn = floor(n/2)
	r$spec[sgn.cnt >= n-1 & sgn.cnt <= nn*(n-nn)] = 's12'
	r$spec[sgn.cnt > nn*(n-nn)] = 's21' #too many sign
	r$spec[apply(is.na(m),3,sum) > dim(m)[1]] = 'NA'
	if(!is.null(cor.sgn)){
		r$spec[apply(is.na(cor.sgn),3,sum) > dim(m)[1]] = 'NA'
	}
	spgroups = list()
	used.groups = c()
	for(i in 1:nn)
		spgroups = c(spgroups,combn(dimnames(m)[[1]],i,simplify = FALSE))
	for(g in 1:length(spgroups)){
		g1 = spgroups[[g]]
		cat(g,paste(g1,collapse=''),length(spgroups),'\n')
		used.groups[length(used.groups)+1]=paste(g1,collapse='')
		g2 = setdiff(dimnames(m)[[1]],g1)
		if(paste(g2,collapse='') %in% used.groups) next
		for(s in 1:dim(m)[3]){
			if(sum(is.na(m[,,s]))>dim(m)[1]) next
			inmin = sort(c(m[g1,g1,s],m[g2,g2,s]))[1]
			outmax = sort(m[g1,g2,s],decreasing = TRUE)[1]
			if(inmin > outmax & max.outmax > outmax & length(intersect(clust[s,g1],clust[s,g2]))==0 & (is.null(cor.sgn) || sum(cor.sgn[g1,g2,s],na.rm=TRUE)==length(g1)*length(g2))){
				r[s,2:3] = c(inmin,outmax)
				if(!(r[s,1] %in% c('n','s5','s12','s21','NA'))) stop(paste('abmiguous classification for exon ',s,', group "',paste(g1,collapse=''),'": "',r[s,1],'"',sep=''))
				r[s,1] = paste(g1,collapse='')
			}
		}
	}
	r$spec.short = r$spec
	r$spec.short[r$spec %in% species$short] = 'one'
	r$spec.short[!(r$spec.short %in% c('one','hq','mr','oc','mrb','s5','s12','s21','NA','n'))] = 'wrong'
	r
}

testTwoSpecies = function(os,s1,s2,m,disp.param=NULL,return.pv = FALSE){
	m = m[m$species %in% c(s1,s2),]
	os = os[,rownames(m)]
	pv = fitSAGLM(os,terms(x ~ s+t+a1+a2+a3 + t:a1+t:a2+t:a3 + s:a1+s:a2+s:a3 + s:t + s:t:a1+s:t:a2+s:t:a3,keep.order=TRUE),list(s=m$species,a1=m$age.use,a2=m$age.use^2,a3=m$age.use^3,t=m$tissue),pseudocount = 0.05,.parallel = TRUE,return.pv = TRUE,disp.param=disp.param)
	if(return.pv)
		qv = pv[,-1]
	else
		qv = apply(pv[,-1],2,p.adjust,method='BH')
	cbind(species=qv[,1],species.age=apply(qv[,9:15],1,min))
}

compASinTwoSpecies = function(os,s1,s2,m){
	m = m[m$species %in% c(s1,s2),]
	os = os[,rownames(m)]
	pv = fitSAGLM(os,terms(x ~ s + a + a2 + t + a:t + a2:t + s:a + s:a2 + s:t + s:t:a + s:t:a2,keep.order=TRUE),list(s=m$species,t=m$tissue,a=m$age.use,a2=m$age.use^2),pseudocount = 0.05,.parallel = TRUE,return.pv = TRUE)
	apply(pv[,-1],2,p.adjust,method='BH')
}

testSpeciesSpecASintissue = function(os,s,t,m){
	m = m[m$tissue==t,]
	os = os[,rownames(m)]
	s = factor(m$species == s)
	a = m$age.use
	pv = fitSAGLM(os,terms(x ~ s + a + I(a^2) + s:a + s:I(a^2),keep.order=TRUE),list(s=s,a=a),pseudocount = 0.05,.parallel = TRUE,return.pv = TRUE)
	apply(pv[,-1],2,p.adjust,method='BH')
}


my.regmatches = function(s,p,from=1,to=nchar(s),return.coor=FALSE,...){
	if(from<0){
		from. = from
		from = nchar(s) + to + 1
		to = nchar(s) + from. + 1
	}
	s = substr(s,from,to)
	rg = gregexpr(p,s,...)
	if(return.coor)
		return(rg)
	regmatches(s,rg)
}

plotQKIgene = function(from,to,h.ann,main=''){
	library(GenomicAlignments)
	# K562
	qki.clip1 = read.table('input/QKI.eCLIP/466_01.basedon_466_01.peaks.l2inputnormnew.bed.compressed.bed.narrowPeak.bed')
	qki.clip1 = qki.clip1[qki.clip1$V1=='chr6' & qki.clip1$V3>from & qki.clip1$V2<to,]
	qki.clip2 = read.table('input/QKI.eCLIP/466_02.basedon_466_02.peaks.l2inputnormnew.bed.compressed.bed.narrowPeak.bed')
	qki.clip2 = qki.clip2[qki.clip2$V1=='chr6' & qki.clip2$V3>from & qki.clip2$V2<to,]
	
	# HEPG2
	qki.clip3 = read.table('input/QKI.eCLIP/ENCFF454GLW.bed')
	qki.clip3 = qki.clip3[qki.clip3$V1=='chr6' & qki.clip3$V3>from & qki.clip3$V2<to,]
	qki.clip4 = read.table('input/QKI.eCLIP/ENCFF544QKV.bed')
	qki.clip4 = qki.clip4[qki.clip4$V1=='chr6' & qki.clip4$V3>from & qki.clip4$V2<to,]
	
	param = ScanBamParam(which=GRanges("chr6", IRanges(from, to)))
	# K562
	bam1 = readGAlignments('input/QKI.eCLIP/ENCFF372LCI.bam',param = param)
	c1 = coverage(bam1)
	bam2 = readGAlignments('input/QKI.eCLIP/ENCFF494AKD.bam',param = param)
	c2 = coverage(bam2)
	mock = readGAlignments('input/QKI.eCLIP/ENCFF849FKQ.mock.bam',param = param)
	mo = coverage(mock)
	
	# HEPG2
	bam3 = readGAlignments('input/QKI.eCLIP/ENCFF179YDO.bam',param = param)
	c3 = coverage(bam3)
	bam4 = readGAlignments('input/QKI.eCLIP/ENCFF081ZJS.bam',param = param)
	c4 = coverage(bam4)
	mock = readGAlignments('input/QKI.eCLIP/ENCFF015GLL.mock.bam',param = param)
	mk = coverage(mock)
	
	x=from:to
	
	qki.fa = read.fasta('input/qki.human.6:163831977-164013409.fa',as.string = TRUE)[[1]]
	
	qki.seg = h.ann[h.ann$gene_id=='hum_148187',]
	start = my.ge$human$gene['hum_148187','start']
	
	plot(c(start,my.ge$human$gene['hum_148187','stop']),c(0.4,0.4),ylim=c(0,3),t='l',xlim=c(from,to),main=main,xlab='genome position',ylab='')
	rect(qki.seg$start,sapply(qki.seg$type,switch,INT=0.2,EXN=0,ALT=0),qki.seg$stop,sapply(qki.seg$type,switch,INT=0.6,EXN=0.8,ALT=0.8),col=sapply(qki.seg$type,switch,INT='blue',EXN='green',ALT='red'))
	mot.pos = gregexpr('actaa',qki.fa)[[1]]+start
	rect(mot.pos,0.95,mot.pos,0.85,col='orange',border = 'orange')
	
	lines(x,mo$chr6[x]/max(mo$chr6[x])*0.7+1,col='gray')	
	lines(x,c1$chr6[x]/max(c1$chr6[x])*0.7+1,col='brown')
	lines(x,c2$chr6[x]/max(c2$chr6[x])*0.7+1,col='magenta')
	
	rect(qki.clip1$V2,1.75,qki.clip1$V3,1.85,col='brown',border = 'brown')
	rect(qki.clip2$V2,1.88,qki.clip2$V3,1.98,col='magenta',border = 'magenta')
	
	
	lines(x,mk$chr6[x]/max(mk$chr6[x])*0.7+1,col='gray')	
	lines(x,c3$chr6[x]/max(c3$chr6[x])*0.7+2,col='brown')
	lines(x,c4$chr6[x]/max(c4$chr6[x])*0.7+2,col='magenta')
	
	rect(qki.clip3$V2,2.75,qki.clip3$V3,2.85,col='brown',border = 'brown')
	rect(qki.clip4$V2,2.88,qki.clip4$V3,2.98,col='magenta',border = 'magenta')
	
}

plotEnrichHexClustOverlap = function(hexs,all.hex,...){
	f = function(m){
		r = matrix(0,ncol=nrow(m),nrow=nrow(m),dimnames = list(rownames(m),rownames(m)))
		for(i in 1:(nrow(m)-1))
			for(j in (i+1):nrow(m)){
				r[i,j]=r[j,i] = 1-sum(m[i,] & m[j,])/sqrt((sum(m[i,])*sum(m[j,])))
			}
		as.dist(r)
	}
	t = sapply(hexs,function(x){all.hex %in% x})
	t = t[apply(t,1,sum)>0,]
	t = t[,apply(t,2,sum)>0]
	heatmap(t+0,scale='none',distfun = f,...)
}

plotHexClustEnrichPv = function(m,main=''){
	image(1:nrow(m),1:ncol(m),-log(m[,ncol(m):1]),col=getPal(),xaxt='n',yaxt='n',xlab='motif',ylab='cluster',main=main)
	axis(1,1:nrow(m),rownames(m),las=2)
	axis(2,ncol(m):1,1:ncol(m),las=2)
	t = round(-log10(m))
	text(rep(1:nrow(m),times=ncol(m)),rep(ncol(m):1,each=nrow(m)),t)
}


#' @param o list with pair orths
#' @param id seg id
#' @param s1 species of id (single letter)
#' @param s2 species for that all pairs are avaliable
#' @param ss all other pseices
getOrthsByPairs = function(o,id,s1,s2,ss){
	res = rep(NA,length(ss)+2)
	names(res) = c(s1,s2,ss)
	res[1] = id
	x = o[[paste(s1,s2,sep='')]][o[[paste(s1,s2,sep='')]][,1] == id,2]
	if(length(x) != 0)
		res[2] = x
	if(!is.na(res[2]))
		for(s in ss){
			x = o[[paste(s,s2,sep='')]][o[[paste(s,s2,sep='')]][,2] == res[2],1]
			if(length(x) != 0)
				res[s] = x
		}
	res
}

plotAllSpesiecCLusterNtFreqTrend = function(t,n,r,substr.mean=FALSE){
	t = lapply(t,function(x){x[[r]][[n]]})
	if(substr.mean)
		t = lapply(t,scale,scale=FALSE)
	plot(1,t='n',xlim=c(0,28),ylim=range(unlist(t)),xlab='cluster',ylab=paste(n,'frequency'),main=paste(n,'; sequence region',r))
	for(s in rownames(species))
		lines(1:28,t[[s]],t='b',pch=19,col = c(rep(c('red','blue'),times=8),rep('gray',12)))
	if(substr.mean)
		abline(h=0,lty=2)
}

plotNtFreqByClust = function(s,n,min.len,from,to,seg.cod.in = c('c','n','p')){
	main = paste(s,' ',n,' [',from,',',to,']',sep='')
	n=toupper(n)
	seq = ad.fa[[s]]
	seq = seq[nchar(seq)>=min.len]
	seq = seq[anns[[s]][names(seq),'cod'] %in% seg.cod.in]
	if(from<0){
		from = nchar(seq) + from
		to = nchar(seq) + to
	}
	seq = toupper(substr(seq,from,to))
	seq = sapply(seq,function(s){mean(strsplit(s,'')[[1]]==n)})
	r = split(seq,orth.ad.cl.reassign[names(seq),'orth.cl'])
	r = r[as.character(1:28)]
	boxplot(r,outline=FALSE,col=c(rep(c('red','blue'),times=8),rep('gray',12)),main=main,xlab='cluster',ylab=paste(n,'frequency'))
	invisible(r)
}

plotSpeciesAgeProile = function(d,m,tissue,ylim=range(d,na.rm=TRUE),df=4,remove.mean=FALSE,...){
	m = m[m$tissue == tissue,]
	d = d[rownames(m)]
	if(remove.mean)
		for(s in unique(m$species))
			d[m$species==s] = d[m$species==s] - mean(d[m$species==s],na.rm=TRUE)
	plot(m$age.use,d,col=m$species.col,pch=m$pch,cex=m$cex,xaxt='n',xlab='Age (days)',ylim=ylim,...)
	
	for(s in unique(m$species)){
		x = m$age.use[m$species == s]
		y = d[m$species == s]
		x = x[!is.na(y)]
		y = y[!is.na(y)]
		df. = min(df,length(unique(y))-2)
		if(length(unique(y))<4 | length(unique(x))<4){
			points(x,y,t='b',col=m$species.col[m$species == s][1])
		}else{
			l=predict(smooth.spline(x,y,df=df.),seq(from=min(x),to=max(x),length.out=50))
			lines(l,lwd=3,col=m$species.col[m$species == s][1])
		}
	}
	plotAgeAxis(log)
}


testOrthTissueAS = function(d,m,species,tissue,glm=TRUE){
	d = d[,m[colnames(d$i),'tissue'] == tissue]
	m = m[colnames(d$i),]
	species.f = factor(m$species %in% species)
	if(glm)
		p = fitSAGLM(d,terms(x ~ ss + a + I(a^2) + I(a^3) + s:a + s:I(a^2) + s:I(a^3),keep.order = TRUE),list(a=m$age.use,ss=m$species,s=species.f),return.pv=TRUE,pseudocount=0.05,.parallel = FALSE)
	else{
		p = matrix(NA,nrow=length(d),ncol=8)
		data = data.frame(a=m$age.use,ss=m$species,s=species.f)
		for(i in 1:length(d)){
			data$x = d$ir[i,]
			try({
				p[i,-1] = anova(lm(terms(x ~ ss + a + I(a^2) + I(a^3) + s:a + s:I(a^2) + s:I(a^3),keep.order = TRUE),data=data))[1:7,5]
			},TRUE)
		}
	}
	for(i in 2:ncol(p))
		p[,i] = p.adjust(p[,i])
	r = cbind(species=p[,2],age=apply(p[,3:5],1,min,na.rm=T),species.age=apply(p[,6:8],1,min,na.rm=T))
	r[is.na(r) | is.infinite(r)] = 1
	rownames(r) = rownames(d$seg)
	r
}

plotNmerCount = function(t,col,add=F){
	x = log(apply(t,2,function(x){sum(p.adjust(x,m='BH')<0.05)})+1)
	if(!add)
		plot(x,yaxt='n',xaxt='n',pch=19,col=col,xlab='cluster',ylab='# of 6mers')
	else
		points(x,yaxt='n',pch=19,col=col)
	arrows(0:7*2+1,x[0:7*2+1],0:7*2+2,x[0:7*2+2],col=col)
}
calcAllPairsFT = function(m){
	r = matrix(NA,ncol=ncol(m),nrow=ncol(m),dimnames=list(colnames(m),colnames(m)))
	for(i in 1:(nrow(r)-1))
		for(j in (i+1):nrow(r)){
			t = fisher.test(factor(m[,i],levels = c(FALSE,TRUE)),factor(m[,j],levels = c(FALSE,TRUE)))
			r[i,j] = t$p.value
			r[j,i] = t$estimate
		}
	r
}


plotFTMotifSimMatrix = function(m,plot.nums=FALSE,pv.thr=0.01,pv.col=c(gray=1,yellow=1e-2,orange=1e-5,red=1e-20,none=0),or.col=c(blue=0,'#0000FFAA'=1/3,'#0000FF55'=1/1.5,gray=1,yellow=1.5,orange=3,red=Inf),main='',leg.cex=1,diag.col='white',diag.text=colnames(m)){
	stop('function is in paper.figures.5.F.R file')
}


plotMotFreqInTwoClusts = function(m,c1,c2,main='',leg=c(c1,c2),col=c(c('#FF000077','#0000FF77'))){
	p1u  = getPosForMot(ad.fa,orth.ad.cl.reassign,c1,m,1,200,0:20*10,TRUE)
	p2u  = getPosForMot(ad.fa,orth.ad.cl.reassign,c2,m,1,200,0:20*10,TRUE)
	p1d  = getPosForMot(ad.fa,orth.ad.cl.reassign,c1,m,-1,-200,0:20*10,TRUE)
	p2d  = getPosForMot(ad.fa,orth.ad.cl.reassign,c2,m,-1,-200,0:20*10,TRUE)
	
	p1eu  = getPosForMot(ad.fa,orth.ad.cl.reassign,c1,m,201,250,0:5*10,TRUE)
	p2eu  = getPosForMot(ad.fa,orth.ad.cl.reassign,c2,m,201,250,0:5*10,TRUE)
	p1ed  = getPosForMot(ad.fa,orth.ad.cl.reassign,c1,m,-201,-250,0:5*10,TRUE)
	p2ed  = getPosForMot(ad.fa,orth.ad.cl.reassign,c2,m,-201,-250,0:5*10,TRUE)
	
	space = c(rep(0,19),0.2,rep(0,4),0.2,rep(0,4),0.2,rep(0,20))
	b=barplot(c(p1u,p1eu,p1ed,p1d),main=paste(main,' (',m,' ',c1,'-',c2,')',sep=''),space = space,border=NA,ylim=range(0,p1u,p1eu,p1ed,p1d,p2u,p2eu,p2ed,p2d),col=col[1])
	barplot(c(p2u,p2eu,p2ed,p2d),add=T,col=col[2],space = space,border=NA)
	legend('topleft',fill=col,legend = leg)
	axis(1,c((b[1]+b[19])/2,(b[19]+b[20])/2,(b[29]+b[30])/2,(b[31]+b[50])/2),c("5'-intron","exon","exon","3'-intron"),las=1)
}



applyIrwin.hall2or = function(x){
	rnks = lapply(x$ft,'[[','or')
	rnks = lapply(rnks,function(x){apply(1/x,2,rank,tie='random')})
	rnks = lapply(rnks,function(x){sweep(x,2,apply(x,2,max),'/')})
	applyIrwin.hall(rnks)
}

testConsNmerEnrich = function(fa,orth.cl,N=6,removeN=TRUE){
	print('Count Nmers:')
	up6mers = fa
	if(removeN)
		up6mers = lapply(up6mers,function(s){s[!grepl('n',s,TRUE)]})
	up6mers = lapply(names(up6mers),function(s){cat('\t',s,'\n');countNmers(up6mers[[s]],N)})
	names(up6mers) = names(fa)
	
	print('Make Fisher Test:')
	up6mers.ft = lapply(names(up6mers),function(s){cat('\t',s,'\n');testMotifInClustEnrich(up6mers[[s]],orth.cl[rownames(up6mers[[s]]),'orth.cl'],TRUE)})
	names(up6mers.ft) = names(up6mers)
	
	rnks = lapply(up6mers.ft,'[[','pv')
	rnks = lapply(rnks,function(x){apply(x,2,rank,tie='random')})
	rnks = lapply(rnks,function(x){sweep(x,2,apply(x,2,max),'/')})
	r = applyIrwin.hall(rnks)
	list(cons.pv=r,ft=up6mers.ft)#,nmers=up6mers)
}

applyIrwin.hall = function(rnks){
	r  = rnks[[1]]
	r[] = NA
	for(m in 1:nrow(r)){
		cat('\r',m,'         ')
		for(c in 1:ncol(r)){
			sum.rank = 0
			for(s in 1:length(rnks)) sum.rank = sum.rank + rnks[[s]][m,c]
			r[m,c] = dirwin.hall(sum.rank,length(rnks))
		}
	}
	r
}

dirwin.hall = function(q,n){
	#it is probability, not density
	r = 0
	for(i in 0:q){
		r = r + ((-1)^i)*choose(n,i)*((q-i)^n)
	}
	r/factorial(n)
}

testMotifInClustEnrich = function(occ,cl,f){
	occ = occ[f & cl > 0,]
	cl = cl[f & cl>0]
	r = matrix(nrow=ncol(occ),ncol=max(cl))
	rownames(r) = colnames(occ)
	colnames(r) = 1:max(cl)
	tot = apply(occ,2,sum)
	or = fr = r
	tot.mot = apply(occ,2,sum)
	for(c in 1:max(cl)){
		#cat('\r',c)
		cl.size = sum(cl==c)
		for(m in 1:ncol(occ)){
			motINcl = sum(cl==c & occ[,m])
			t = cbind(c(nrow(occ)-tot.mot[m]-cl.size+motINcl,tot.mot[m]-motINcl),c(cl.size-motINcl,motINcl))
			f=fisher.test(t,a='g')
			r[m,c] = f$p.value
			or[m,c] = f$estimate
			fr[m,c] = motINcl/cl.size
		}
	}
	list(pv=r,or=or,freq=fr)
}

countNmers = function(seq,n){
	seq = tolower(seq)
	nmers = sapply(1:(4^n),function(i){c('a','t','g','c')[((i-1)  %% (4^(1:n)) %/% (4^(0:(n-1)))) + 1]})
	nmers = apply(nmers,2,paste,collapse='')
	r = matrix(0,ncol=length(nmers),nrow=length(seq))
	rownames(r) = names(seq)
	colnames(r) = nmers
	for(i in 1:length(seq)){
		if(i %% 1000 == 1) cat('\r',i)
		s = seq[i]
		t = table(sapply(1:(nchar(s)-n+1),function(j)substr(s,j,j+n-1)))
		#	t = t[names(t) %in% nmers]
		r[i,names(t)] = t
	}
	gc()
	r
}

getEns2My = function(f){
	t = read.table(f,sep='\t')
	t = t[t$V3 %in% c('=','c','j','e'),]
	t$V4 = sapply(strsplit(t$V4,',',T),function(x){as.numeric(x[3])}) 
	t$V3 = c('='=1,'c'=2,'j'=3,'e'=4)[t$V3]
	t = split(t,t$V2)
	t = sapply(t,function(x){
		x[order(x$V3,-x$V4)[1],1]
	})
	t
}

checkGeneClusterOverlap = function(orth.cl,sp,c1,c2,value=c('table','test','genes'),cl.max.dist=Inf,cl.min.dist=-Inf){
	orth.cl = orth.cl[orth.cl$species == sp & orth.cl$dist <= cl.max.dist & orth.cl$dist >= cl.min.dist,]
	t = unique(anns[[sp]][rownames(orth.cl),'gene_id'])
	g1 = unique(anns[[sp]][rownames(orth.cl)[orth.cl$orth.cl == c1],'gene_id'])
	g2 = unique(anns[[sp]][rownames(orth.cl)[orth.cl$orth.cl == c2],'gene_id'])
	value= value[1]
	t = table(factor(t %in% g1,levels = c(FALSE,TRUE)),factor(t %in% g2,levels = c(FALSE,TRUE)))
	if(value == 'table')
		return(t)
	if(value == 'test')
		return(fisher.test(t))
	if(value == 'genes')
		return(intersect(g1,g2))
}

plotClustOrthEnrich = function(sps,sp.to,cls,FUN,seg=TRUE,orth.cl,main=''){
	pv = ov = fq = matrix(NA,ncol=length(sps),nrow=length(cls),dimnames = list(cls,sps))
	for(cl in cls){
		for(s1 in sps){
			m = FUN(cl,s1,sp.to)
			pv[cl,s1] = fisher.test(m,alternative = 'greater')$p.value
			orth.size = min(m[1,2],m[2,1]) + m[2,2]
			ov[cl,s1] = m[2,2]/orth.size
			if(seg){
				#part of cluster segments that have orthologs in other species
				fq[cl,s1] = orth.size/min(sum(orth.cl$species==s1 & orth.cl$orth.cl == cl),sum(orth.cl$species==sp.to & orth.cl$orth.cl == cl))
			}else{
				s2 = sp.to
				e1 = unique(unlist(anns[[s1]][rownames(orth.cl)[orth.cl$species==s1 & orth.cl$orth.cl == cl],'ens.id']))
				e2 = unique(unlist(anns[[s2]][rownames(orth.cl)[orth.cl$species==s2 & orth.cl$orth.cl == cl],'ens.id']))
				fq[cl,s1] = orth.size/min(sum(!is.na(e1)),sum(!is.na(e2)))
			}
		}
	}
	image(1:6,1:28,t(pv[rev(cls),]),xaxt='n',yaxt='n',xlab='species',ylab='clusters',col=rev(c('white','yellow','orange','red')),breaks=c(0,1e-60,1e-15,5e-2,28*6)/28/6,main=main)
	text = paste(round(fq*100),'/',round(ov*100),sep='')
	text(rep(1:6,each=length(cls)),rep(28:1,times=length(sps)),text,c(0.5,0.5))
	axis(1,1:6,sps,las=2)
	axis(2,1:28,28:1,las=2)
	abline(h=26.5 - (0:7)*2)
}

compOrthSegsByFT = function(orth.cl,orth.segs,pair=TRUE,cl,s1,s2,value=c('matrix','test','rate')){
	if(pair){
		o = orth.segs[[paste(species[s1,'short'],species[s2,'short'],sep='')]]
	}else
		o = orth.segs[['hqmrboc']][,c(s1,s2)]
	o = o[o[,1] %in% rownames(orth.cl) & o[,2] %in% rownames(orth.cl),] # I take segments that are significant in both species
	value = value[1]
	incl1  = orth.cl[o[,1],'orth.cl'] == cl
	incl2  = orth.cl[o[,2],'orth.cl'] == cl 
	if(value=='rate'){
		return(sum(incl1 & incl2)/min(sum(incl1),sum(incl2)))
	}
	r = table(factor(incl1,levels = c(FALSE,TRUE)),factor(incl2,levels = c(FALSE,TRUE)))
	if(value=='matrix')
		return(r)
	r = fisher.test(r)
	r
}

compOrthGenesByFT = function(cl.gids,cl,s1,s2,value=c('matrix','test','rate')){
	t = intersect(cl.gids[[s1]]$total,cl.gids[[s2]]$total)
	value = value[1]
	if(value=='rate'){
		a = intersect(t,cl.gids[[s1]][[cl]])
		b = intersect(t,cl.gids[[s2]][[cl]])
		return(length(intersect(a,b))/min(length(a),length(b)))#length(union(a,b)))
	}
	r = table(t %in% cl.gids[[s1]][[cl]],t %in% cl.gids[[s2]][[cl]])
	if(value=='matrix')
		return(r)
	r = fisher.test(r)
	r
}


plotOrthClusters = function(pams,sites,sps,orth.cl,m){
	ncl = max(pams[[1]][[sites]]$clustering[,1])
	orth.cl = split(orth.cl,orth.cl$orth.cl)
	
	orth.cl = orth.cl[order(sapply(orth.cl,function(x){length(unique(x$species))}),sapply(orth.cl,function(x){-min(x$sp.cl)}),decreasing=TRUE)]
	cat("nrow = ",sum(sapply(orth.cl,function(x){max(table(x$species))})),'\n')
	for(ocn in names(orth.cl)){
		oc = orth.cl[[ocn]]
		nrow = max(table(oc$species))
		for(cc in 1:nrow){
			for(s in sps){
				sp.cl.inx = sort(oc[oc$species==s,'sp.cl'])
				if(length(sp.cl.inx) < cc)
					plot.new()
				else{
					cnt = sum(pams[[s]][[sites]]$clustering[,1] == sp.cl.inx[cc])
					plotTissueAgeProile(pams[[s]][[sites]]$means[sp.cl.inx[cc],],m,main=paste(s,', cl',sp.cl.inx[cc],' (',cnt,')',sep=''),ylab='')
				}
				if(s == sps[1])
					mtext(paste('ocl',ocn,sep=''),2,1,outer = FALSE)
			}
		}
		abline(h=grconvertY(0,from='nfc','user'),xpd=NA)
	}
}

addEnsIDToAnn = function(comp,ens.gtf,seg){
	c = read.table(comp,comment.char = '#')
	c = c[c[,3] %in% c('=','c','j'),]
	m2e = split(c[,2],c[,1])
	e = read.table(ens.gtf,sep='\t',comment.char = '#')
	e = e[e[,3]=='exon',]
	gid = sapply(regmatches(e[,9],regexec('gene_id "?([^";]+)"?;',e[,9])),'[',2)
	e = lapply(split(e,gid),function(x){
		data.frame(x[1,1],x[1,7],min(x[,4:5]),max(x[,4:5]))
	})
	e = do.call(rbind,e)
	
	eidss = m2e[seg$gene_id]
	seg$ens.id = vector('list',nrow(seg))
	for(i in 1:nrow(seg)){
		eids = eidss[[i]]
		if(is.null(eids)) next
		r = list()
		for(j in 1:length(eids)){
			eid = eids[j]
			if(seg[i,'start'] < e[eid,3] || seg[i,'stop'] > e[eid,4])
				eids[j] = NA
		}
		seg$ens.id[[i]] = as.list(eids[!is.na(eids)])
	}
	seg
}



#' @param adj adj matrix. should be positive.
lookForBestHits = function(adj,groups,lookForInhits = FALSE,hits=1){
	diag(adj) = 0
	r = adj
	r[] = 0
	
	for(g1. in 1:length(groups)){
		g1 = groups[[g1.]]
		if(lookForInhits)
			for(e in g1){
				t = order(adj[e,g1],decreasing = TRUE)[1]
				if(sum(adj[e,g1[t]] < adj[e,-g1])==0)
					r[e,g1[t]] = adj[e,g1[t]]
			}
		for(g2. in setdiff(1:length(groups),g1.)){
			g2 = groups[[g2.]]
			for(e in g1){
				t = order(adj[e,g2],decreasing = TRUE)[1:hits]
				r[e,g2[t]] = adj[e,g2[t]]
			}
		}
	}
	r
}

symmGraph = function(m){
	for(i in 1:nrow(m))
		for(j in i:ncol(m))
			if(m[i,j] != m[j,i])
				m[i,j] = m[j,i] = 0
			m
}


getMyLayout = function(n=30){
	function(g){
		a = as_adjacency_matrix(g,sparse = FALSE)
		link.comps = clusters(g,mode = "weak")$membership
		
		r = matrix(-1,ncol=2,nrow=ncol(a))
		r[,2] = -rep(1:(ncol(a)/n),each=n)
		group.count = nrow(a)/n
		
		for(gr in 1:group.count){
			gr.inx = ((gr-1)*n+1):(gr*n)
			x = max(r[min(gr.inx):nrow(r),1],na.rm=TRUE) + 1
			for(v in gr.inx){
				if(is.na(link.comps[v])) next #this vertex is already placed
				r[!is.na(link.comps) & link.comps == link.comps[v],1] = x
				link.comps[!is.na(link.comps) & link.comps == link.comps[v]] = NA
				x = x + 1
			}
		}
		# move overlapped nodes
		for(i in 1:nrow(r)){
			p = which(r[,1] == r[i,1] & r[,2] == r[i,2])
			if(length(p)>1){
				m = mean(r[p,2])
				r[p,2] = seq(m-0.3,m+0.3,length.out = length(p))
				#m = mean(r[p,1])
				#r[p,1] = seq(m-0.3,m+0.3,length.out = length(p))
			}
		}
		r*2
	}
}

approxPSI = function(psi,m,ages,age2pred,tissues,df=4){
	m$ages = ages
	m = m[colnames(psi),]
	a = m$ages
	N = length(age2pred)
	r = matrix(NA,nrow=nrow(psi),ncol=N*length(tissues))
	for(c in 1:nrow(psi))
		for(t in 1:length(tissues)){
			a.tmp = a[m$tissue==tissues[t]]
			if(length(a.tmp) < df) next
			ss = smooth.spline(a.tmp,psi[c,m$tissue==tissues[t]],df=df)
			p = predict(ss,age2pred)
			p$y[p$x < min(a.tmp) | p$x > max(a.tmp)] = NA
			r[c,(t-1)*N + (1:N)] = p$y
		}
	r
}



tranfromAges = function(m,suff,mult,addGestation){
	l = nchar(m$age)
	f = substr(m$age,l-nchar(suff)+1,l) == suff
	cat(suff,': ',sum(f)," values set, ",sum(f & !is.na(m$days)),' where already defined',"\n",sep='')
	m$days[f] = addGestation*species[m$species[f],'gestation'] + mult*as.numeric(gsub(suff,'',m$age[f]))
	m
}

getAnnOverlap = function(a1,a2,return.matches=FALSE){
	t=Sys.time()
	seqinf = Seqinfo(unique(c(a1$chr_id,a2$chr_id)))
	if(is.numeric(a2$strand))
		a2$strand = ifelse(is.na(a2$strand),'*',ifelse(a2$strand==1,"+","-"))
	a2 = GRanges(a2$chr_id ,IRanges(a2$start,a2$stop),a2$strand,seqinfo = seqinf)
	if(!return.matches)
		a2 = reduce(a2)
	
	a1 = GRanges(a1$chr_id,IRanges(a1$start,a1$stop),ifelse(is.na(a1$strand),'*',ifelse(a1$strand==1,"+","-")),seqinfo = seqinf)
	overlap = findOverlaps(a1,a2,type='any',ignore.strand = FALSE,select='all')
	cat('map1: ',as.numeric(Sys.time()-t, units = "secs"),"\n")
	include = findOverlaps(a1,a2,type="within",ignore.strand = FALSE,select='all')
	cat('map2: ',as.numeric(Sys.time()-t, units = "secs"),"\n")
	if(return.matches)
		return(list(overlap=cbind(from=overlap@from,to=overlap@to),include=cbind(from=include@from,to=include@to)))
	overlap = overlap@from
	include = include@from
	r = rep('-',length(a1))
	r[overlap] = 'o'
	r[include] = 'i'
	r
}


addIsCogingByEnsGTF = function(ens.gtf,as){
	t=Sys.time()
	
	cds = read.table(ens.gtf,sep = '\t',header=FALSE,comment.char = '#')[,c(1,3:5,7)]
	cds = cds[cds[,2] == 'CDS',-2]
	colnames(cds) = c('chr_id','start','stop','strand')
	
	seqinf = Seqinfo(unique(c(cds$chr_id,as$seg$chr_id)))
	cds = GRanges(cds$chr_id ,IRanges(cds$start,cds$stop),cds$strand,seqinfo = seqinf)
	cds = reduce(cds)
	
	seg.r = GRanges(as$seg$chr_id,IRanges(as$seg$start,as$seg$stop),ifelse(is.na(as$seg$strand),'*',ifelse(as$seg$strand==1,"+","-")),seqinfo = seqinf)
	overlap = findOverlaps(seg.r,cds,type='any',ignore.strand = FALSE,select='all')@from
	cat('map1: ',as.numeric(Sys.time()-t, units = "secs"),"\n")
	include = findOverlaps(seg.r,cds,type="within",ignore.strand = FALSE,select='all')@from
	cat('map2: ',as.numeric(Sys.time()-t, units = "secs"),"\n")
	as$seg$cod = 'n'
	as$seg$cod[overlap] = 'p'
	as$seg$cod[include] = 'c'
	cod.genes = as$seg$gene_id[as$seg$cod != 'n']
	as$seg$cod.gene = as$seg$gene_id %in% cod.genes
	as
}

plotAllMDSsforSubset = function(dist,meta,tissues=NULL,main=''){
	meta = meta[meta$tissue %in% tissues & rownames(meta) %in% rownames(dist),]
	dist = dist[rownames(meta),rownames(meta)]
	mds=isoMDS(dist,k=3)$points
	
	plot(mds[,1],mds[,2],pch=meta$pch,col=meta$col,cex=meta$cex,xlab='Dimension 1',ylab='Dimension 2',main=main)
	lg = do.call(rbind,lapply(split(meta,meta$tissue.orig),function(x){x[1,c('tissue.orig','col','pch')]}))
	lg = lg[order(lg$col),]
	colnames(lg)[1] = 'tissue'
	legend('topleft',pch=lg[,3],col=lg[,2],legend = lg[,1])
	plot(mds[,2],mds[,3],pch=meta$pch,col=meta$col,cex=meta$cex,xlab='Dimension 2',ylab='Dimension 3',main=main)
	
	tissues = unique(meta[,c('tissue','col')])
	for(d in 1:ncol(mds)){
		plot(meta$age.use,mds[,d],pch=meta$pch,col=meta$col,xlab='Age (days)',ylab=paste('Dimension',d),xaxt='n')
		plotAgeAxis()
		for(t in 1:nrow(tissues)){
			ft = meta$tissue == tissues[t,1]
			age = meta$age.use[ft]
			ages2pred = seq(from=min(age),to=max(age),length.out = 100)
			if(length(unique(age))>4){
				ss=smooth.spline(age,mds[ft,d],df=3)
				p = predict(ss,ages2pred)
				lines(p,col=tissues[t,2],lwd=2)
			}
		}
	}
	invisible(mds)
}

plotDivOnAge = function(d,m,t,FUN=function(x,y){mean(abs(x-y),na.rm=T)},ylab='Distance',main='',ylim=NULL){
	m = m[colnames(d),]
	
	tissues = setdiff(unique(m$tissue),t)
	ages = table(m$age.use)
	ages = sort(as.numeric(names(ages)[ages > 1]))
	
	r = matrix(ncol=length(ages),nrow=length(tissues))
	rownames(r) = tissues
	colnames(r) = ages
	
	this.d = d[,m$tissue==t]
	this.d = this.d[,order(m[colnames(this.d),'age.use'])]
	this.a = m[colnames(this.d),'age.use']
	r[] = NA
	for(tissue in tissues){
		cmn.age = sort(intersect(this.a,m$age.use[m$tissue==tissue]))
		curr.d  = d[,m$tissue==tissue & m$age.use %in% cmn.age]
		curr.d  = curr.d[,order(m[colnames(curr.d),'age.use'])]
		for(i in 1:length(cmn.age)){
			r[tissue,as.character(cmn.age[i])] = FUN(this.d[,this.a %in% cmn.age][,i],curr.d[,i])
		}
	}
	if(is.null(ylim)) ylim = range(r,na.rm=T)
	plot(1,t='n',xlim=range(ages),ylim=ylim,xlab='Age',xaxt='n',ylab=ylab,main=main)
	plotAgeAxis()
	for(tissue in tissues){
		mm = m[m$tissue==tissue,]
		mm = mm[order(mm$days),]
		f = !is.na(r[tissue,])
		lines(ages[f],r[tissue,f],t='b',lwd=2,pch=mm$pch,cex=mm$cex,col=mm$col)
	}
	invisible(r)
}

calcMeanCols = function(d,f,FUN=base::mean){
	u = unique(as.character(f))
	r = matrix(ncol=length(u),nrow=nrow(d))
	colnames(r) = u
	rownames(r) = rownames(d)
	for(i in u){
		print(i)
		r[,i] = apply(d[,f==i,drop=F],1,FUN,na.rm=TRUE)
	}
	r
}


plotLengthDistr = function(l,max=max(l),log=FALSE,...){
	l = l[l<=max]
	h=hist(l,0:max(l),plot=FALSE)$counts
	if(log){
		barplot(log(h+1),col=c('gray','darkgray','olivedrab'),border = NA,space = 0,xlab='Segment length',ylab='# of segments',yaxt='n',...)
		at = c(0,5,10,20,50,100,500,1000)
		axis(2,at=log(at+1),labels = at)
	}else
		barplot(h,col=c('gray','darkgray','olivedrab'),border = NA,space = 0,xlab='Segment length',ylab='# of segments',...)
	legend('topright',fill=c('gray','darkgray','olivedrab'),legend=c(1,2,0),title = 'Phase')
	axis(1,at = seq(from=0,to=max,by = 25))
	invisible(h)
}


plotPhaseDistr = function(l,max=max(l),...){
	l = l[l<=max]
	p = l %% 3
	l = floor(l / 3)
	t=table(p,l)
	t = sweep(t,2,apply(t,2,sum),'/')
	barplot(t,col=c('olivedrab','darkgray','gray'),xaxt='n',border = NA,space = 0,xlab='Segment length (codons)',ylab='proportion of segments',...)
	abline(h=c(1/3,2/3),lty=2)
	axis(1,at = seq(from=0,to=max(l),by=20))
	invisible(t)
}

calcInterSpeciesCorr = function(psi,method='sp',min.obs=10,noise=0.01,seed=1234){
	set.seed(seed)
	r = array(NA,dim=c(length(psi),length(psi),nrow(psi[[1]])),dimnames=list(names(psi),names(psi),rownames(psi[[1]])))
	for(i in 1:dim(r)[3]){
		cat('\r',i,dim(r)[3])
		for(s in 1:length(psi))
			psi[[s]][i,] = psi[[s]][i,] + rnorm(ncol(psi[[s]]),sd=noise)
		for(s1 in 1:(length(psi)-1))
			for(s2 in (s1+1):length(psi)){
				f = !is.na(psi[[s1]][i,]) & !is.na(psi[[s2]][i,])
				if(sum(f)>= min.obs)
					r[s1,s2,i] = r[s2,s1,i] = cor(psi[[s1]][i,f],psi[[s2]][i,f],method=method)
			}
	}
	r
}

plotTissueAgeProile = function(d,m,df=5,error.bar=NULL,ylim=NULL,xlim=NULL,by.species=FALSE,age.axis='log',tissues=unique(m$tissue),col=NULL,add=FALSE,pch=NULL,lty=1,cex=NULL,xlab='Age (days)',plot.xaxt=TRUE,lwd=3,age=NULL,...){
	if(!is.null(pch))
		m$pch = pch
	if(!is.null(col))
		m$col = col
	if(!is.null(cex))
		m$cex = cex
	if(age.axis == 'rank')
		m$age.use = m$age.rank
	m = m[m$tissue %in% tissues,]
	if(sum(!is.na(d)) == 0){
		plot(1,xaxt='n',t='n',xlab=xlab,...)
		return()
	}
	if(by.species){
		m$tissue=m$species
		m$col=m$species.col
	}
	if(!is.null(age))
		m$age.use = age
	cmn = intersect(rownames(m)[!is.na(m$age.use)],names(d))
	m = m[cmn,]
	d = d[cmn]
	if(is.null(ylim))
		ylim = range(d,error.bar,na.rm=T)
	if(is.null(xlim))
		xlim = range(m$age.use)
	if(add)
		points(m$age.use,d,col=m$col,pch=m$pch,cex=m$cex)
	else
		plot(m$age.use,d,col=m$col,pch=m$pch,cex=m$cex,xaxt='n',xlab=xlab,ylim=ylim,xlim=xlim,...)
	if(!is.null(error.bar)){
		segments(m$age.use,error.bar[,1],m$age.use,error.bar[,2],col=m$col)
	}
	for(t in unique(m$tissue)){
		x = m$age.use[m$tissue == t]
		y = d[m$tissue == t]
		nna = !is.na(y) & !is.na(x)
		x = x[nna]
		y = y[nna]
		o = order(x)
		x = x[o]
		y = y[o]
		df. = min(df,length(unique(x))-2)
		if(df.<=1){
			points(x,y,t='b',col=m$col[m$tissue == t][1],lwd=3,lty=lty,cex=0)
		}else{
			l=predict(smooth.spline(x,y,df=df.,tol=1e-10),seq(from=min(x),to=max(x),length.out=50))
			lines(l,lwd=lwd,col=m$col[m$tissue == t][1],lty=lty)
		}
	}
		if(plot.xaxt){
			if(age.axis == 'log')
				plotAgeAxis(log)
			else if(age.axis == 'rank'){
				ds = unique(m[,c('stage','age.use')])
				ds = ds[order(ds$age.use),]
				axis(1,ds$age.use,ds$stage)
				# d = unique(m[,c('days','age.use')])
				# d = unique(do.call(rbind,lapply(split(d,d$age.use),function(x){x$days = round(mean(x$days));x})))
				# d = d[order(d$age.use),]
				# d$days = d$days-species[m$species[1],'gestation']
				# ad = abs(d$days)
				# d$days = ifelse(ad<30,paste0(d$days,'d'),ifelse(ad<365,paste0(round(d$days/30),'m'),paste0(round(d$days/365),'y')))
				# axis(1,d$age.use,d$days)
			}
	}
	invisible(list(x=m$age.use,y=d))
}

plotAgeAxis = function(fun=function(x){x^0.25}){
	at = c(11,20,45,80,150,365,365*2,365*5,365*10,365*30,365*60)
	labs = paste(ifelse(at<365,at,at/365),ifelse(at<365,'d','y'),sep='')
	axis(1,at=fun(at),labels = labs)
}

testASTissueAge = function(d,m){
	m = m[colnames(d$ir),]
	p = fitSAGLM(d,terms(x ~ t + a + I(a^2) + I(a^3) + t:a + t:I(a^2) + t:I(a^3),keep.order = TRUE),list(a=m$age.use,t=m$tissue),return.pv=TRUE,pseudocount=0.05,.parallel = TRUE)
	for(i in 2:ncol(p))
		p[,i] = p.adjust(p[,i])
	cbind(tissue=p[,2],age=apply(p[,3:5],1,min),tissue.age=apply(p[,6:8],1,min))
}

testASAge = function(d,m,t,min.cov.sams=0,return.pv=FALSE,.parallel=TRUE){
	m = m[m$tissue==t & rownames(m) %in% colnames(d$ir),]
	d = d[,rownames(m)]
	p = fitSAGLM(d,formula(x ~ a + I(a^2) + I(a^3)),list(a=m$age.use),0.05,.parallel = .parallel,return.pv=TRUE)
	p[apply(!is.na(d$ir),1,mean) < min.cov.sams,] = NA
	if(return.pv)
		return(p)
	for(i in 2:ncol(p))
		p[,i] = p.adjust(p[,i],method = 'BH')
	apply(p[,-1],1,min)
}

plotAllDistHM = function(scls,add.scales=NULL,dist.pow=1){
	for(s in names(scls)){
		m=scls[[s]]$distMatrix^dist.pow
		image(1:nrow(m),1:ncol(m),m,xlab=s,ylab='mouse',main=s,col=terrain.colors(500))
		a=scls[[s]]$align
		lines(a$q,a$r,lwd=5,col='blue')
		if(!is.null(add.scales)){
			a=add.scales[[s]]$align
			lines(a$q,a$r,lwd=3,col='brown',lty=2)
		}
	}
}

loadInfoForOrths = function(o){
	info = vector('list',ncol(o))
	names(info) = colnames(o)
	for(s in colnames(o)){
		print(s)
		info[[s]] = readRDS(paste('Rdata/',s,'.as.u.all.Rdata',sep=''))[o[!is.na(o[,s]),s],]
		gc(verbose = FALSE)
	}
	info
}

loadAltOrthSegs = function(fs,only.filtered=TRUE){
	orth = NULL
	for(f in fs)
		orth = rbind(orth,read.table(f,comment.char = '#',na.strings = 'None',sep='\t'))
	if(is.numeric(orth[,1]))
		orth = orth[,-1]
	orth = unique(orth)
	orth = apply(orth,2,function(x){x[grepl('_',substr(x,5,100),fixed=TRUE)] = NA;x})
	orth = orth[apply(is.na(orth),1,sum)==0,]
	colnames(orth) = setNames(rownames(species),substr(rownames(species),1,3))[substr(unlist(orth[1,]),1,3)]
	specs = colnames(orth)#[-ncol(orth)]
	if(only.filtered){
		for(s in specs){
			orth = orth[is.na(orth[,s]) | orth[,s] %in% rownames(anns[[s]]),]
		}
	}
	orth
}

myPAM = function(d,k,meta,dist.fun=c('cor','pearson'),core.set.max.nas=0,max.dist.to.medoids=0.35,min.no.tissue.sam=0,seed=NULL){
	if(!is.null(seed)) set.seed(seed)
	#set PSI to NA in all segments/tissues where number of defined PSIs are less than min.no.tissue.sam
	meta = meta[colnames(d),]
	for(t in unique(meta$tissue)){
		tf = meta$tissue == t
		d[apply(!is.na(d[,tf]),1,sum)<min.no.tissue.sam,tf] = NA
	}
	
	f = apply(is.na(d),1,sum) <= core.set.max.nas #select core set to evide NAs in distance matrix
	dist = calcPairDistance(d[f,],fun.name = dist.fun[1],method = dist.fun[2])
	core.size1 = nrow(dist)
	dist = cleanNADistMatrix(dist)
	core.size2 = nrow(dist)
	cat('core consists of ',core.size2,' segments out of ',length(f),'; ',(core.size1-core.size2),' elements were removed from core due to NAs in dist matrix',"\n",sep='')
	
	pam.cl = pam(dist,k = k,diss = TRUE)
	cl = matrix(NA,ncol=2,nrow=nrow(d))
	
	colnames(cl) = c('cluster','distance')
	rownames(cl) = rownames(d)
	dist2meds = calcPairDistance(d,d[pam.cl$medoids,],fun.name = dist.fun[1],method = dist.fun[2])
	for(i in 1:nrow(cl)){
		o = order(dist2meds[i,])[1]
		cl[i,] = c(o,dist2meds[i,o])
	}
	cl[cl[,2]>max.dist.to.medoids,1] = 0
	
	cl[cl[,1] != 0,1] = reorderClustersBySize(cl[cl[,1] != 0,1])
	means = matrix(NA,nrow=k,ncol=ncol(d))
	colnames(means) = colnames(d)
	for(i in 1:k)
		means[i,] = apply(d[cl[,1]==i,,drop=FALSE],2,mean,na.rm=TRUE)
	
	list(original.pam = pam.cl,clustering=cl,medoids=pam.cl$medoids[order(cl[pam.cl$medoids,1])],means=means)
}

cleanNADistMatrix = function(d){
	nas = apply(is.na(d),1,sum)
	ord = (1:ncol(d))[order(nas,decreasing = TRUE)]
	if(nas[ord[1]] == 0) return(d)
	for(i in 1:length(ord)){
		if(sum(is.na(d[-ord[1:i],-ord[1:i]]))==0)
			return(d[-ord[1:i],-ord[1:i]])
	}
}


calcPairDistance = function(a,b=NULL,min.obs=1,fun.name='cor',method='pearson'){
	if(is.null(b)){
		if(fun.name == 'cor')
			r = 0.5 - cor(t(a),u='p',method=method)/2
		else
			r = as.matrix(dist(a,method=method))
	}else{
		r = matrix(NA,nrow=nrow(a),ncol=nrow(b))
		rownames(r) = rownames(a)
		colnames(r) = rownames(b)
		for(j in 1:nrow(b)){
			for(i in 1:nrow(a)){
				if(sum(!is.na(a[i,]) & !is.na(b[j,]))>= min.obs)
					if(fun.name == 'cor')
						r[i,j] = 0.5 - cor(a[i,],b[j,],u='p',method=method)/2
					else
						r[i,j] = dist(rbind(a[i,],b[j,]),method=method)[1]
			}
		}
	}
	# if(fun.name != 'cor'){
	# 	if(method == 'euclidian')
	# 		r = r / sqrt(ncol(r))
	# 	else if(method == 'manhattan')
	# 		r = r / ncol(r)
	# 	else if(method != 'maximum')
	# 		stop('Unexpected method!')
	# }
	r
}

plotClustering = function(pam.cl,psi,meta,norm=FALSE,ylab='mean(PSI)',mtitle='',meanORmedoids = 'both',vert.line=NULL,age.axis='log',plot.mtext=TRUE,plots.with.yaxis=NULL,...){
	clnames = c()
	r = NULL
	if(!is.list(pam.cl))
		pam.cl = list(clustering=pam.cl)
	if(!is.matrix(pam.cl$clustering))
		pam.cl$clustering = cbind(pam.cl$clustering,pam.cl$clustering)
	if(meanORmedoids %in% c('means','both'))
		for(i in 1:max(pam.cl$clustering[,1])){
			f = pam.cl$clustering[,1]==i
			t = psi[f,,drop=FALSE]
			if(norm)
				t = normRows(t)
			x = apply(t,2,mean,na.rm=T)
			r = rbind(r,x)
			sd = apply(t,2,sd,na.rm=T)/sqrt(nrow(t))*2
			clnames[i] = paste('c',i,' (',sum(f),')',sep='')
			plotTissueAgeProile(x,meta,error.bar = cbind(x-sd,x+sd),main=clnames[i],ylab=ylab,age.axis=age.axis,yaxt='n',...)
			if(is.null(plots.with.yaxis) || i %in% plots.with.yaxis)
				axis(2)
			if(!is.null(vert.line)) abline(v=vert.line,lty=2)
			if(i == 1 & plot.mtext) mtext(paste(mtitle,'. Cluster average.',sep=''),3,outer = TRUE)
		}
	if(meanORmedoids == 'both') par(mfrow=par('mfrow'))
	
	if(meanORmedoids %in% c('medoids','both'))
		for(i in 1:max(pam.cl$clustering[,1])){
			f = pam.cl$clustering[,1]==i
			clnames[i] = paste('c',i,' (',sum(f),')',sep='')
			r = rbind(r,psi[pam.cl$medoids[i],])
			plotTissueAgeProile(psi[pam.cl$medoids[i],],meta,ylab='PSI',main=clnames[i],plot.age.axis=plot.age.axis)
			if(!is.null(vert.line)) abline(v=vert.line,lty=2)
			if(i == 1 & plot.mtext) mtext(paste(mtitle,'. Cluster medoids'),3,outer = TRUE)
		}
	r
}

reorderClustersBySize = function(clusts){
	t = sort(table(clusts),decreasing=T)
	clusts. = clusts
	for(c in 1:length(t))
		clusts[clusts.==names(t)[c]] = c
	clusts
}


caclCurveDivergence = function(x1,y1,x2=NULL,y2=NULL,max.df=4,n=100){
	# approximates curve by spline with df equal to half of number of data points but not more than max.df
	na1 = is.na(x1) | is.na(y1)
	x1 = x1[!na1]
	y1 = y1[!na1]
	if(!is.null(x2)){
		na2 = is.na(x2) | is.na(y2)
		x2 = x2[!na2]
		y2 = y2[!na2]
	}
	uniq.cnt = length(unique(x1))
	if(uniq.cnt < 4)
		return(NA)
	x = seq(max(x1,x2),min(x1,x2),length.out = n)
	yp1 = predict(smooth.spline(x1,y1,df=min(max.df,uniq.cnt/2)),x)$y
	if(!is.null(x2))
		yp2 = predict(smooth.spline(x2,y2,df=min(max.df,length(x1)/2)),x)$y
	else
		yp2 = rep(0,n)
	# plot(x1,y1,ylim=range(y1,yp1))
	# lines(x,yp1)
	yp1 = yp1 - mean(yp1)
	yp2 = yp2 - mean(yp2)
	mean(abs(yp1-yp2))
}

plotPSIchange = function(dp,s,f,thr=0.1,breaks,...){
	dp = dp[f]
	s = s[f]
	#breaks = seq(from = 0,to=1,length.out=breaks)
	if(sum(!is.na(dp))>0){
		hist(dp[!s],breaks = breaks,col='#66666666',freq = FALSE,border = '#66666666',...)
		hist(dp[ s],breaks = breaks,col='#FF000066',add=TRUE,freq = FALSE,border='#FF000066')
	}else
		plot(1,t='n',...)
	abline(v=thr,lty=2)
	legend('topright',col=c('gray','red','white'),lwd=2,legend=c(
		paste("!sign: ",sum(!s & dp >= thr,na.rm=T),'/',sum(!s,na.rm=T)),
		paste(" sign: ",sum( s & dp >= thr,na.rm=T),'/',sum( s,na.rm=T)),
		paste(" sign but NA: ",sum( s & is.na(dp),na.rm=T),'/',sum( s,na.rm=T))))
}

plotASStat = function(sgn=NULL,seg=NULL,matrix=NULL,frac=FALSE,col=c(ad='red',da='blue',aa='green',dd='magenta'),den=c(c=-1,p=30,n=15),tissue.order=c('testis','brain','cerebellum','heart','liver','kidney','ovary'),plot=TRUE,...){
	if(is.null(matrix)){
		f = paste(seg$sites,seg$cod)
		g = seg$sites %in% names(col)
		sgn = sgn[g,,drop=FALSE]
		o = paste(rep(names(col),each=length(den)),rep(names(den),times=length(col)))
		f = factor(f[g],levels = o)
		matrix = apply(sgn,2,function(s){table(f[s])})[o,]
		if(!is.null(tissue.order))
			matrix = matrix[,tissue.order]
		else
			matrix = matrix[,order(-apply(matrix,2,sum))]
	}
	
	if(frac)
		matrix = sweep(matrix,2,apply(matrix,2,sum),'/')
	col=rep(col,each=length(den))
	density = rep(den,times=length(col))
	if(plot){
		barplot(matrix,col=col,density=density,las=2,...)
		legend('topright',fill=col,density=density,legend=rownames(matrix),ncol = 2)
	}
	invisible(matrix)
}

plotIntronExonStructure = function(a,new=TRUE,ylim=c(0,1),xlim=NULL,...){
	gstart = min(a$start)
	gstop = max(a$stop)
	if(all(a$strand==-1)){
		gstop = min(a$start)
		gstart = max(a$stop)
	}
	if(is.null(xlim))
		xlim=c(gstart,gstop)
	if(new)
		plot(1,t='n',xlim=xlim,ylim=ylim,xlab=paste0('Chr ',a$chr_id[1]),yaxt='n',ylab='',...)
	
	segments(gstart,mean(ylim),gstop,mean(ylim))
	
	f = a$type=='EXN'
	if(sum(f)>0)
		rect(a$start[f],ylim[1],a$stop[f],ylim[2],border = NA,col='green')
	f = a$type!='EXN' & a$sites %in% c('ad','sd','ae')
	if(sum(f)>0)
		rect(a$start[f],ylim[1],a$stop[f],ylim[2],border = NA,col='red')
	d = ylim[2]-ylim[1]
	f = a$type!='EXN' & !f & a$sites != 'da'
	if(sum(f)>0)
		rect(a$start[f],ylim[1]+d/4,a$stop[f],ylim[2]-d/4,border = NA,col='orange')
	f = a$sites == 'da'
	if(sum(f)>0)
		rect(a$start[f],ylim[1]+d/8,a$stop[f],ylim[2]-d/8,border = NA,col='blue')
}

