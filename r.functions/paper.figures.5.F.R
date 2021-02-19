getAltSp = function(p,thr,n){
	sps =species[names(p),'short']
	apply(sapply(p,function(x){apply(x<thr,1,sum,na.rm=T)>n}),1,function(x)paste(sps[x],collapse=''))
}


# 1 ######
getSegTestDevAsStat = function(patt,ann,orth=NULL){
	astypes = c(CE='ad',AA='aa',AD='dd',RI='da')
	tested = devas = list()
	for(ss in astypes){
		tested[[ss]] = devas[[ss]] = matrix(0,ncol=length(patt),nrow = ncol(patt[[1]]),dimnames=list(colnames(patt[[1]]),names(patt)))
		for(s in names(patt)){
			f = ann[[s]]$sites==ss
			if(!is.null(orth))
				f = f & rownames(patt[[s]]) %in% orth[,s]
			tested[[ss]][,s] = apply(patt[[s]][f,]!='-',2,sum)[rownames(tested[[ss]])]
			devas[[ss]][,s]  = apply(patt[[s]][f,]!='-' & patt[[s]][f,]!='n',2,sum)[rownames(tested[[ss]])]
		}
	}
	devasf = setNames(lapply(names(tested),function(ss)devas[[ss]]/tested[[ss]]),names(tested))
	list(tested=tested,devasn=devas,devasf=devasf)
}

plot1A = function(l){
	# it is not ready
	tissues = sort(unique(meta$tissue))
	plot(phyl.tree,x.lim=c(0,312*4),show.tip.label=F)
	at = c(300,180,90,25)
	axis(1,312-at,at,las=2)
	grconvertY(0:1,'npc','user')
	
	x = grconvertX(312,'user','npc')
	par(xpd=NA)
	nbx  = grconvertX(c(0.40,0.67),'npc','user')
	altx = grconvertX(c(0.73,1),'npc','user')
	for(i in 1:7){
		plotPNG(paste0('figures/paper.figures/5/icons/',phyl.tree$tip.label[i],'.png'),x+0.06,grconvertY(8-i,'user','npc'),0.1)
		plotPNG(paste0('figures/paper.figures/5/icons/',tissues[i],'.png'),2*x+0.06,grconvertY(8-i,'user','npc'),0.1)
	}
	plotPanelLetter(l,lab.cex)
}


plotAsEventCount = function(d,pchs=c(ad=19),by.tissue=FALSE,sps=NULL,tiss=NULL,ylim=NULL,xlab.cex=0.8,cil=NULL,cih=NULL,short.species.lab=FALSE,...){
	if(is.matrix(d)){
		d = list(ad=d)
		if(!is.null(cil)) cil = list(ad=cil)
		if(!is.null(cil)) cih = list(ad=cih)
	}
	if(!is.null(sps))
		d = lapply(d,function(x)x[,sps])
	if(!is.null(tiss))
		d = lapply(d,function(x)x[tiss,])
	
	
	dim = dim(d[[1]])
	if(by.tissue)
		dim = rev(dim)
	x = 1:length(d[[1]]) + rep((1:dim[2])*3,each=dim[1]) - 1
	if(is.null(ylim))
		ylim = range(0,unlist(list(d,cil,cih)),na.rm=T)*1.05
	plot(1,t='n',xaxt='n',xlab='',xlim=range(x),ylim=ylim,yaxs = 'i',...)
	for(ss in names(pchs)){
		d[[ss]][d[[ss]]==0] = NA
		if(by.tissue){
			points(x,t(d[[ss]]),pch=pchs[ss],col=rep(params$tissue.col[rownames(d[[ss]])],each=ncol(d[[ss]])))
			if(!is.null(cil))
				arrows(x,t(cil[[ss]]),x,t(cih[[ss]]),col=rep(params$tissue.col[rownames(d[[ss]])],each=ncol(d[[ss]])),angle=90,code=3,length=0.03)
		}else{
			points(x,  d[[ss]] ,pch=pchs[ss],col=rep(params$tissue.col[rownames(d[[ss]])],times=ncol(d[[ss]])))
			if(!is.null(cil))
				arrows(x,cil[[ss]],x,cih[[ss]],col=rep(params$tissue.col[rownames(d[[ss]])],times=ncol(d[[ss]])),angle=90,code=3,length=0.03)
		}
	}
	if(by.tissue){
		
		segments(x,grconvertY(0.0,'npc','user'),x,grconvertY(-0.02,'npc','user'),xpd=T)
		xx = x[1:ncol(d[[1]])]
		xx = 2.2*(xx - mean(xx))/7*nrow(d[[1]]) + mean(xx)*min(1,0.7/nrow(d[[1]])*7)
		
		# adjust q and r
		if(ncol(d[[1]]) == 7){
			xx[2] = xx[2] - (max(xx)-min(xx))/50
			xx[4] = xx[4] + (max(xx)-min(xx))/50
		}
		
		segments(x[1:ncol(d[[1]])],grconvertY(-0.02,'npc','user'),xx,grconvertY(-0.03,'npc','user'),xpd=T)
		text(xx,grconvertY(-0.03,'npc','user'),species[colnames(d[[1]]),'short'],adj = c(0.5,1),xpd=T,cex=1)
		# text(x,0,species[colnames(d[[1]]),'short'],adj = c(0.5,1),xpd=T,cex=0.5)
		# xx = sapply(split(x,rep(1:7,each=7)),mean)
		# text(xx,0,rownames(d[[1]]),adj = c(0.5,2.3),xpd=T,cex=0.7)
	}else{
		# text(x,0,substr(rownames(d[[1]]),1,1),adj = c(0.5,1),xpd=T,cex=0.5)
		xx = sapply(split(x,rep(1:ncol(d[[1]]),each=nrow(d[[1]]))),mean)
		slab = colnames(d[[1]])
		if(short.species.lab)
			slab = species[slab,'short']
		text(xx,0,slab,adj = c(0.5,2.3),xpd=T,cex=xlab.cex)
	}
}

plotASTypes = function(pchs,cnst.col='gray',alt.col='orange',cex=0.8){
	plot(1,t='n',bty='n',xaxt='n',yaxt='n',xlab='',ylab='',main='',xlim=c(0,100),ylim=c(0,100),xaxs = 'i',yaxs = 'i')
	h = 1.2 #half of exon height
	top = 90-h # top line
	x1 = 30 # left coor of cschems
	x2 = 95 # right coor of cschems
	w = 10 # with of exons
	f = function(ss,ssn,dd,t=''){
		points(x1-13,top-dd-h*0.75,pch=pchs[ss],cex=1)
		#text(14,top-dd-h*0.75,adj = c(0,0.5),ssn,cex=1)
		
		segments(x1,top-dd,x2,top-dd)
		segments(x1,top-dd-h*8/3,x2,top-dd-h*8/3)
		rect(c(x1,x2-w),top-dd+h,c(x1+w,x2),top-dd-h,col=cnst.col,border=cnst.col)
		rect(c(x1,x2-w),top-dd-h*5/3,c(x1+w,x2),top-dd-h*11/3,col=cnst.col,border=cnst.col)
		if(t=='')
			return(dd+h*6)
		text(10,top-dd+1.8*h,adj = c(0,0),t,cex=cex,xpd=NA)
		return(dd+h*12)
	}
	
	# CE
	dd0=f('ad','CE',h*0,'Cassette\nExon')
	rect((x1+x2-w/2)/2,top+h,(x1+x2+w/2)/2,top-h,col=alt.col,border=alt.col)
	# AA
	dd1=f('dd','AD',dd0,'Alternative\nDonor')
	rect(x1+w,top+h-dd0,x1+w*1.5,top-h-dd0,col=alt.col,border=alt.col)
	dd0 = dd1
	# DD
	dd1 = f('aa','AA',dd0,'Alternative\nAcceptor')
	rect(x2-w*1.5,top+h-dd0,x2-w,top-h-dd0,col=alt.col,border=alt.col)
	dd0 = dd1
	# DA
	dd1=f('da','RI',dd0,'Retained\nIntron')
	rect(x1+w,top+h-dd0,x2-w,top-h-dd0,col=alt.col,border=alt.col)
	top - dd0-h*3
}

# 2 #####
caclCor2EmbryoAllForTissue = function(v,m,t,cor.m='p',boot=0,ci=0.95){
	cat(t,': ')
	m = m[colnames(v),]
	v = v[,m$tissue==t]
	m = m[colnames(v),]
	stages = sort(sapply(split(m$days,m$stage),mean))
	s1 = which(m$stage == names(stages)[1])
	r=lapply(names(stages),function(s){
		cat(s,' ')
		s2 = which(m$stage == s)
		r = matrix(NA,ncol=3,nrow=length(s1)*length(s2))
		j = 1
		for(i1 in s1){
			for(i2 in s2){
				p1 = v[,i1]
				p2 = v[,i2]
				f = !is.na(p1) & !is.na(p2)
				p1 = p1[f]
				p2 = p2[f]
				r[j,1] = cor(p1,p2,meth=cor.m)
				if(boot == 0)
					r[j,2:3] = r[j,1]
				else{
					bs = laply(1:boot,function(i){
						s = sample(length(p1),replace = TRUE)
						cor(p1[s],p2[s],meth=cor.m)
						},.parallel=TRUE)
					r[j,2:3] = quantile(bs,prob=c((1-ci)/2,0.5+ci/2))
				}
				j = j + 1
			}
		}
		r
	})
	names(r) = names(stages)
	cat('\n')
	r
}

caclCor2EmbryoAllStat = function(v,m,cor.m='p',stat.fun=function(x)c(median(x[,1]),min(x[,2]),max(x[,3])),boot=0,ci=0.95){
	m = m[colnames(v),]
	r = lapply(unique(m$tissue),function(t)caclCor2EmbryoAllForTissue(v,m,t,cor.m,boot=boot,ci=ci))
	names(r) = unique(m$tissue)
	stages = sort(sapply(split(m$days,m$stage),mean))
	rt = array(NA,dim=c(7,length(stages),3),dimnames = list(unique(meta$tissue),names(stages),c('rho','ci1','ci2')))
	for(t in unique(m$tissue)){
		for(s in names(stages))
			if(!is.null(r[[t]][[s]]))
				rt[t,s,] = stat.fun(r[[t]][[s]])
	}
	rt
}


caclCor2Embryo = function(d,m,cor.m='p',use.mean.embryo=TRUE,use.bootstrap=0){
	cmn = intersect(rownames(m),colnames(d))
	d=d[,cmn]
	m=m[cmn,]
	days = unique(m[,c('days','stage')])
	days = days[order(days$days),]
	tissues = unique(m$tissue)
	r = array(NA,dim = c(length(tissues),nrow(days),3),dimnames = list(tissues,days$stage,c('rho','ci1','ci2')))
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
				if(use.bootstrap==0){
					tmp = cor.test(base,d[,ff],use='pair',method = cor.m)
					r[t,s,] = c(tmp$estimate,tmp$conf.int)
				}else{
					y = d[,ff]
					boots = sapply(1:use.bootstrap,function(i){
						s = sample(length(y),replace = TRUE)
						cor(base[s],y[s],use='pair',method = cor.m)
						})
					r[t,s,] = c(cor(base,y,use='pair',method = cor.m),quantile(boots,prob=c(0.025,0.975)))
				}
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

getsPSIbyEnsID = function(dpsi,s2e,use.random=FALSE){
	#getsPSIbyEnsID_ is an old version
	dpsi = dpsi[!is.na(dpsi)]
	r=do.call(rbind,lapply(names(dpsi),function(sid){
		eids = s2e[[sid]]
		if(is.null(eids))
			return(NULL)
		data.frame(eids,rep(dpsi[sid],length(eids)))
	}))
	sapply(split(r[,2],r[,1]),function(x){
		if(use.random){
			t = sample(length(x),1)
			c(up=x[t],down=x[t])
		}else
			c(up=max(x),down=min(x)) # now down doesn't make any sence
	})
}

getsPSIbyEnsID_ = function(psi,stages,tissue,s2e,sp,use.random=FALSE){
	stop("It is outdated version that uses ddirected dPSI")
	psi = psi[[sp]]
	stages = stages[[sp]][tissue,]
	s2e = s2e[[sp]]
	dpsi = psi[,paste(sp,tissue,stages[2])] - psi[,paste(sp,tissue,stages[1])]
	dpsi = dpsi[!is.na(dpsi)]
	r=do.call(rbind,lapply(names(dpsi),function(sid){
		eids = s2e[[sid]]
		if(is.null(eids))
			return(NULL)
		data.frame(eids,rep(dpsi[sid],length(eids)))
	}))
	sapply(split(r[,2],r[,1]),function(x){
		if(use.random){
			t = sample(length(x),1)
			c(up=x[t],down=x[t])
		}else
			c(up=max(x),down=min(x))
	})
}

plotAsInTauDistr = function(ge.tau,dpsi,dpsi.thr=0.2,cols=c('darkgray','lightgray','orange'),norm.by.col=F,...){
	ge.tau$as = ifelse(rownames(ge.tau) %in% rownames(dpsi),'y','!AS')
	ge.tau$as[ge.tau$as=='y'] = ifelse(apply(abs(dpsi[rownames(ge.tau)[ge.tau$as=='y'],]),1,max)>dpsi.thr,'devAS','!devAS')
	ge.tau$TissueTau[!is.na(ge.tau$TissueTau) & ge.tau$TissueTau==1] = max(ge.tau$TissueTau[!is.na(ge.tau$TissueTau) & ge.tau$TissueTau<1])
	r=table(ge.tau$as,floor(ge.tau$TissueTau*10)/10)[c('!AS','!devAS','devAS'),]
	if(norm.by.col)
		r = sweep(r,2,apply(r,2,sum),'/')
	barplot(r,names.arg = seq(0.05,0.95,by=0.1)[1:ncol(r)],col=cols,ylab='# genes',xlab='TissueTau',...)
	invisible(ge.tau)
}

# 3 #####
getDevASCons = function(s1,s2,t,q,d,f,thr=0.2){
	sg1 = q[[s1]][,t]<0.05 & abs(d[[s1]][,t]) > thr
	sg1[is.na(sg1)] = FALSE
	
	sg2 = q[[s2]][,t]<0.05 & abs(d[[s2]][,t]) > thr
	sg2[is.na(sg2)] = FALSE
	
	ff = f & sg1 & sg2
	same.dir=sign(d[[s1]][ff,t])==sign(sign(d[[s2]][ff,t]))
	bn = my.binom.test(sum(ff),sum(sg1 & f)-sum(ff))
	data.frame(s1=s1,s2=s2,t=t,n1=sum(sg1 & f),n12=sum(ff),same.dir.cnt=sum(same.dir),also.sgn=mean(sg2[sg1 & f]),same.dir=mean(same.dir),ci1=bn[2],ci2=bn[3])
}

calcDevASCor = function(psi,age.al,s1,s2,t,fseg,cor.meth='pearson',return.stat=TRUE){
	c1 = paste(s1,t,age.al[,s1])
	c2 = paste(s2,t,age.al[,s2])
	f = c1 %in% colnames(psi[[s1]]) & c2 %in% colnames(psi[[s2]])
	c1 = c1[f]
	c2 = c2[f]
	p1 = psi[[s1]][fseg,c1]
	p2 = psi[[s2]][fseg,c2]
	r = sapply(1:nrow(p1),function(i)cor(p1[i,],p2[i,],m=cor.meth,u='p'))
	if(return.stat){
		sd = sd(r,na.rm=T)/sqrt(length(r))
		m=mean(r,na.rm=T)
		return(c(mean=m,ci1=m-2*sd,ci2=m+2*sd))
	}
}

plotExampleDLG3 = function(fig=c(0,1,0,1)){
	fig[3] = fig[3] + grconvertY(1,'lines','ndc')
	sid1 = 12502 ##Dlg3 - looks nice + Intellectual disability (HGMD)
	sid2 = 8396
	sid3 = 12506
	h = fig[4]-fig[3]
	w = fig[2]-fig[1]
	rows=6
	ad = dlg3.mdata$ann[!(dlg3.mdata$ann$sites %in% c('aa','dd','da')),]
	alt.seg = rownames(ad) %in% rownames(orth.age.dpsi$mouse)[c(sid1,sid2,sid3)]
	alt.seg.coor = c(min(dlg3.mdata$ann$start[alt.seg]),max(dlg3.mdata$ann$stop[alt.seg]))
	par(fig=c(fig[1],fig[2],fig[4]-1/rows*h,fig[4]),new = TRUE)
	plotReadCov(dlg3.mdata$cov,reverse = (dlg3.mdata$ann$strand[1] == -1),min.junc.cov = 5,bty='n',xlab=paste0('Chr ',dlg3.mdata$ann$chr_id[1]),ylab='',junc.col='black',junc.lwd=1,axes=F)
	y=grconvertY(-0.1,'nfc','user')
	segments(min(dlg3.mdata$ann$stop),y/2,max(dlg3.mdata$ann$start),y/2,xpd=NA)
	rect(ad$start,y*0.8,ad$stop,y*0.2,xpd=NA,border=NA,col=ifelse(alt.seg,'red','black'))
	alt.seg.y = grconvertY(y,'user','ndc')
	#plotPanelLetter(l,lab.cex)
	
	zoom.coor.ndc = grconvertX(dlg3.mdata$zoom.coor,'user','ndc')
	
	f = ad$start>=dlg3.mdata$zoom.coor[1] & ad$stop<=dlg3.mdata$zoom.coor[2]
	ad = ad[f,]
	alt.seg = alt.seg[f]
	par(fig=c(fig[1],fig[2],fig[4]-2/rows*h,fig[4]-1/rows*h),new = TRUE)
	plotReadCov(dlg3.mdata$cov,reverse = (dlg3.mdata$ann$strand[1] == -1),min.junc.cov = 5,bty='n',xlab=paste0('Chr ',dlg3.mdata$ann$chr_id[1]),ylab='',junc.col='black',junc.lwd=1,axes=F,xlim=dlg3.mdata$zoom.coor,plot.junc.only.within = T)
	segments(dlg3.mdata$zoom.coor[1],y/2,dlg3.mdata$zoom.coor[2],y/2,xpd=NA)
	rect(ad$start,y*0.8,ad$stop,y*0.2,xpd=NA,border=NA,col=ifelse(alt.seg,'red','black'))

	alt.seg.coor0 = grconvertX(zoom.coor.ndc,'ndc','user')
	
	alt.seg.y = grconvertY(alt.seg.y,'ndc','user')
	y = 0.9*max(dlg3.mdata$cov$cov[dlg3.mdata$cov$x>=dlg3.mdata$zoom.coor[1] & dlg3.mdata$cov$x<=dlg3.mdata$zoom.coor[2]])
	segments(alt.seg.coor0,c(alt.seg.y,alt.seg.y),dlg3.mdata$zoom.coor,rep(y,2),xpd=NA,lty=2)
	
	for(j in 1:6){
		s = rownames(species)[-2][j]
		par(fig=c(fig[1]+w/6*(j-1),fig[1]+w/6*j,fig[4]-3/rows*h,fig[4]-2/rows*h),mar=c(0,0,0,0),new = TRUE,xpd=NA)
		plot.new()
		plotPNG(paste0("figures/paper.figures/5/icons/",s,".png"),0.6,0.25,0.6)
		par(fig=c(fig[1]+w/6*(j-1),fig[1]+w/6*j,fig[4]-4/rows*h,fig[4]-3/rows*h),mar=c(1,1,0.5,0),new = TRUE)
		plotTissueAgeProile(orth.seg.ad.tsm[[s]][sid1,],meta.tsm,age.axis = 'rank',yaxt=ifelse(j==1,'s','n'),bty='n',ylim=c(0,1),xlab='',ylab=ifelse(j==1,'PSI',''),plot.xaxt=F,pch=19,cex=0.3,lwd=1,age=meta.tsm$corr.age.rank,xlim=c(0,14))
		par(fig=c(fig[1]+w/6*(j-1),fig[1]+w/6*j,fig[4]-5/rows*h,fig[4]-4/rows*h),new = TRUE)
		plotTissueAgeProile(orth.seg.ad.tsm[[s]][sid2,],meta.tsm,age.axis = 'rank',yaxt=ifelse(j==1,'s','n'),bty='n',ylim=c(0,1),xlab='',ylab=ifelse(j==1,'PSI',''),plot.xaxt=F,pch=19,cex=0.3,lwd=1,age=meta.tsm$corr.age.rank,xlim=c(0,14))
		if(j==1) title(ylab='PSI',xpd=NA)
		par(fig=c(fig[1]+w/6*(j-1),fig[1]+w/6*j,fig[4]-6/rows*h,fig[4]-5/rows*h),new = TRUE) 
		plotTissueAgeProile(orth.seg.ad.tsm[[s]][sid3,],meta.tsm,age.axis = 'rank',yaxt=ifelse(j==1,'s','n'),bty='n',ylim=c(0,1),xlab='',ylab=ifelse(j==1,'PSI',''),plot.xaxt=F,pch=19,cex=0.3,lwd=1,age=meta.tsm$corr.age.rank,xlim=c(0,14))
		axis(1,labels = NA)
	}
	x = grconvertX(mean(fig[1:2]),'ndc','user')
	y = grconvertY(fig[3],'ndc','user')
	text(x,y,'Development',xpd=NA,adj=c(0.5,0))
	y = grconvertY(fig[4],'ndc','user')
	text(x,y,substitute(paste(italic('DLG3'), " - membrane-associated guanylate kinase")),adj=c(0.5,1.4),xpd=NA)
}

# 4 ######
calcCor2TisOnDev = function(p,tissue,m,cor.meth='p'){
	m = m[colnames(p),]
	r = NULL;
	stages = unique(m$stage)
	tissues = setdiff(unique(m$tissue),tissue)
	for(t in tissues)
		for(s in stages){
			s1 = which(m$tissue == tissue & m$stage==s)
			s2 = which(m$tissue == t & m$stage==s)
			if(length(s1)==1 & length(s2)==1){
				r = rbind(r,data.frame(id=paste(m$species[1],t,s),cor=cor(p[,s1],p[,s2],m=cor.meth,u='p')))
			}
		}
	setNames(r$cor,r$id)
}

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



plotLine1 = function(x,y,col,area.alpha=0.2,lwd=3,leg.pos=c(0.01,0.01),leg.adj=c(0,0),leg.cex=0.7,add=T,...){
	o = order(x)
	x = x[o]
	y = y[o]
	if(!add)
		plot(x,y,...)
	f = !is.na(y)
	p = predict(lm(y[f] ~ x[f]),interval='conf')
	col.area =  col2rgb(col)[,1]
	col.area = rgb(col.area[1],col.area[2],col.area[3],area.alpha*255,maxColorValue = 255)
	polygon(c(x[f],rev(x[f])),c(p[,2],rev(p[,3])),border=NA,col=col.area,xpd=FALSE)
	lines(x[f],p[,1],col=col,lwd=lwd)
	c = cor.test(x[f],y[f])
	#text=paste0('r^2 = ',round(c$estimate^2,3),'\n','P = ',format(c$p.value,scientific = T,digits = 1))
	r2=round(c$estimate^2,3)
	p=gsub("^-0+","-",strsplit(format(c$p.value,scientific = T,digits = 1),'e')[[1]],perl = T)
	t1 = bquote("r"^"2" ~ "=" ~ .(r2))
	t2 = bquote(italic("P")~ "=" ~ .(p[1])  %*% "10"^~.(p[2]))
	text(grconvertX(leg.pos[1],'npc','user'),grconvertY(leg.pos[2],'npc','user'),t1,adj=leg.adj,cex=leg.cex)
	text(grconvertX(leg.pos[1],'npc','user'),grconvertY(leg.pos[2],'npc','user')+(grconvertY(1,'lines','user')-grconvertY(0,'lines','user'))*leg.cex,t2,adj=leg.adj,cex=leg.cex)
}

plot4C.DivergenceOnAge = function(s1,s2,as.mr.cor,l,yrange=c(0.25,0.5),plot.xlab=T,ylab.pos=0.5,panel.lab.yadj=-0.5,plot.tissue.lab=TRUE){
	#just to organize code
	x = 1:ncol(as.mr.cor)
	left = grconvertX(1,'lines','ndc')
	plot.pos = seq(left,1,length.out = nrow(as.mr.cor)+1)
	for(i in 1:nrow(as.mr.cor)){
		t = rownames(as.mr.cor)[i]
		par(fig=c(plot.pos[i],plot.pos[i+1],yrange[1],yrange[2]),new=T)
		y = as.mr.cor[t,]
		plot(x,y,pch=19,col=params$tissue.col[t],axes=F,main='',xlab='',ylab='',ylim=range(as.mr.cor,na.rm=T))
		if(plot.tissue.lab)
			text(grconvertX(0.5,'npc','user'),grconvertY(1,'nfc','user'),t,adj=c(0.5,1),xpd=NA)
		axis(1,x,colnames(as.mr.cor))
		if(i==1) axis(2)
		plotLine1(x,y,col=params$tissue.col[t],leg.pos=c(0.05,0.05))
		if(i == 1){
			plotPNG(paste0("figures/paper.figures/5/icons/",s1,".png"),0.2,0.3,0.3)
			plotPNG(paste0("figures/paper.figures/5/icons/",s2,".png"),0.7,0.32,0.3)
			arrows(grconvertX(0.37,'npc','user'),grconvertY(0.3,'npc','user'),grconvertX(0.5,'npc','user'),grconvertY(0.3,'npc','user'),lwd=1,length = 0.05)
			
		}
	}
	text(x=grconvertX(0,'ndc','user'),y=grconvertY(yrange[2],'ndc','user'),labels=l,adj=c(0,panel.lab.yadj),font=2,cex=lab.cex,xpd=NA)
	text(grconvertX(left/10,'ndc','user'),grconvertY(ylab.pos,'npc','user'),"Correlation",srt=90,xpd=NA,adj=c(0.5,1))
	if(plot.xlab)
		text(grconvertX(mean(plot.pos),'ndc','user'),grconvertY(yrange[1],'ndc','user'),'Developmental stage',xpd=NA,adj=c(0.5,-0.1))
	
}


plot4D.PeakChange = function(sp,peak.changes,l,yrange = c(0,0.25),plot.tissue.lab=FALSE,plot.xlab=TRUE,plot.inset=TRUE,letter.y.pos=1.7/6.2,sp.icon.y=0.7){
	left = grconvertX(1,'lines','ndc')
	
	plot.pos = seq(left,1,length.out = length(peak.changes)+1)
	f = function(x)log(x+1)
	atx = c(0,10,100,300,1000,2000,5000)
	aty = c(50,100,150,200,300,400)
	for(t in 1:length(peak.changes)){
		par(fig=c(plot.pos[t],plot.pos[t+1],yrange[1],yrange[2]),new=T)
		x = f(peak.changes[[t]]$ge)
		y = f(peak.changes[[t]]$as)
		plot(x,y,pch=19,col=params$tissue.col[t],axes=F,xlab='',ylab='',xpd=NA,ylim=range(f(30),y))
		at. = atx[atx <= max(peak.changes[[t]]$ge)]
		axis(1,f(at.),at.)
		at. = aty[aty <= max(peak.changes[[t]]$as) & aty >= min(peak.changes[[t]]$as)]
		axis(2,f(at.),at.)
		plotLine1(x,y,col=params$tissue.col[t],leg.pos=c(0.95,0.05),leg.adj = c(1,0))
		if(plot.tissue.lab)
			text(grconvertX(0.5,'npc','user'),grconvertY(1,'nfc','user'),names(peak.changes)[t],adj=c(0.5,1),xpd=NA)
		if(t==1){
			par(xpd=T)
			plotPNG(paste0("figures/paper.figures/5/icons/",sp,".png"),0.22,sp.icon.y,0.4)
			par(xpd=F)
			if(plot.inset){
				ff = function(x){x/max(x)}
				ppar=par(fig=c(grconvertX(c(-0.05,0.45),'npc','ndc'),grconvertY(c(0.85,1.15),'nfc','ndc')),new=TRUE,cex=0.4,mar=par('mar')*0.2,mgp=par('mgp')*0,xpd=NA)
				plot(ff(peak.changes$heart$ge),t='l',bty='n',xlab='Development',ylab='# events',axes=F)
				lines(ff(peak.changes$heart$as),col='red')
				text(grconvertX(c(1,1),'npc','user'),grconvertY(c(0.8,0.65),'nfc','user'),c("Exons","Genes"),col=c('red','black'),adj=c(0,1))
				rect(grconvertX(-0.2,'nfc','user'),grconvertY(-0.1,'nfc','user'),grconvertX(1.7,'nfc','user'),grconvertY(1,'nfc','user'),col=NA)
				par(cex=ppar$cex,mar=ppar$mar,mgp=ppar$mgp,xpd=T)
			}
		}
		#if(t==2)
			#plotPNG(paste0("figures/paper.figures/5/icons/",sp,".png"),0.3,0.9,0.4)
	}
	text(x=grconvertX(0,'ndc','user'),y=grconvertY(letter.y.pos,'ndc','user'),labels=l,adj=c(0,1),font=2,cex=lab.cex,xpd=NA)
	text(grconvertX(left/10,'ndc','user'),grconvertY(0.5,'npc','user'),'# exons',srt=90,xpd=NA,adj=c(0.5,1))
	if(plot.xlab)
		text(grconvertX(mean(plot.pos),'ndc','user'),grconvertY(yrange[1],'ndc','user'),'# genes',xpd=NA,adj=c(0.5,-0.4))
}

plot4B.legend = function(){
	y = grconvertY(1,'npc','user')
	top = grconvertY(0.85,'npc','user')
	x = grconvertX(2/5,'npc','user')
	h = y/5

	
	w=grconvertX(grconvertY(c(0,h),'user','inches'),'inches','user')
	w=w[2]-w[1]
	
	yy = (1+cos(seq(0,pi,length.out = 30)))/2*0.84+0.08
	xx = seq(x+0.08*w,x+.92*w,length.out = 30)
	
	#rect(x-w,y-h,x+w,y+h,border = NA,col='#000000')
	segments(x,top-h,x+w,top-h)
	segments(x,top-h,x,top)
	lines(xx,top - h + yy*h,lwd=3,col='#00000030')
	text(mean(xx),top+h/4,'Early\ngenes',adj = c(0.5,0.5))
	rect(x-w/8,top+h/4-h/8,x+w/8,top+h/4+h/8,border = NA,col='#00000030')
	text(x,top-h/2,'Gene\nexpression',srt=90,adj=c(0.5,-0.4))
	
	x = x + w*1.5
	xx = xx + w*1.5
	yy = 1-yy
	segments(x,top-h,x+w,top-h)
	segments(x,top-h,x,top)
	lines(xx,top - h + yy*h,lwd=3,col='#000000')
	text(mean(xx),top+h/4,'Late\ngenes',adj = c(0.5,0.5))
	rect(x-w/8,top+h/4-h/8,x+w/8,top+h/4+h/8,border = NA,col='#000000')
	text(x-.25*w,top-h,'Development',adj=c(0.5,1.5))
	
	#plotPNG("figures/paper.figures/5/icons/mouse.png",0.45,0.9,0.15)
}

# _5 ######
getHexStat = function(hex.ups,hex.dws,pv.thr){
	hex.stat = rbind(apply(hex.ups$up$ih.qv<pv.thr,2,sum),
									 apply(hex.dws$up$ih.qv<pv.thr,2,sum),
									 apply(hex.ups$dw$ih.qv<pv.thr,2,sum),
									 apply(hex.dws$dw$ih.qv<pv.thr,2,sum))
	
	k = rownames(hex.ups$up$ih.qv) %in% rownames(hex2mot)[hex2mot$V2!='']
	hex.stat.known = rbind(apply(hex.ups$up$ih.qv[k,]<pv.thr,2,sum),
												 apply(hex.dws$up$ih.qv[k,]<pv.thr,2,sum),
												 apply(hex.ups$dw$ih.qv[k,]<pv.thr,2,sum),
												 apply(hex.dws$dw$ih.qv[k,]<pv.thr,2,sum))
	list(hex.stat=hex.stat,hex.stat.known=hex.stat.known)
}

plotHexStat = function(hex.stat,hex.stat.known,panles,plot.leg=TRUE,leg.h=450){
	den=50
	barplotWithText(hex.stat[c(1,3),],col=paste0(rep(params$tissue.col,each=2),'55'),las=3,main='upstream',beside = T,border=NA,den=c(-1,den),srt = 90,adj=c(-0.1,0.5),ylim=c(0,max(hex.stat[c(1,3),])*1.2),xaxt='n',ylab='# hexamers')
	b=barplot(hex.stat.known[c(1,3),],col=rep(params$tissue.col,each=2),add=T,beside = T,xaxt='n',yaxt='n',den=c(-1,den),border=NA)
	
	text(apply(b,2,mean),rep(0,7),colnames(hex.stat),xpd=NA,srt=-45,adj=c(0,1))
	plotExonScheme(TRUE)
	plotPanelLetter(panles[1],lab.cex)
	
	barplotWithText(hex.stat[c(2,4),],col=paste0(rep(params$tissue.col,each=2),'55'),las=3,main='downstream',beside = T,border=NA,den=c(-1,den),srt = 90,adj=c(-0.1,0.5),ylim=c(0,max(hex.stat[c(2,4),])*1.2),xaxt='n',ylab='# hexamers')
	b=barplot(hex.stat.known[c(2,4),],col=rep(params$tissue.col,each=2),add=T,beside = T,xaxt='n',yaxt='n',den=c(-1,den),border=NA)
	text(apply(b,2,mean),rep(0,7),colnames(hex.stat),xpd=NA,srt=-45,adj=c(0,1))
	if(plot.leg){
		l=legend(10,leg.h,fill=c('black','black','black','#00000055'),den=c(-1,den,-1,-1),legend=c('up','down','known','unknown'),bty = 'n',border = NA)
		rect(l$rect$left+l$rect$w*0.05,l$rect$top-l$rect$h*0.06,l$rect$left+l$rect$w*0.98,l$rect$top-l$rect$h*0.94)
	}
	
	plotExonScheme(FALSE)
	plotPanelLetter(panles[2],lab.cex)
}

getDevASpattern = function(psi,sgn,dpsi,m,f,at.least.two.sgn=FALSE){
	psi = psi[f,colnames(psi) %in% rownames(m)]
	sgn = sgn[f,colnames(sgn) %in% m$tissue]
	dpsi=dpsi[f,colnames(dpsi) %in% m$tissue]
	
	sgn[is.na(sgn)] = FALSE
	m = m[colnames(psi),]
	if(at.least.two.sgn){
		f = apply(sgn,1,sum)>1
		psi = psi[f,]
		sgn = sgn[f,]
		dpsi=dpsi[f,]
	}
	r = data.frame(dir=character(nrow(psi)),mean.psi=NA,ntissue=0)
	rownames(r) = rownames(psi)
	for(i in 1:nrow(psi)){
		r$ntissue[i] = sum(sgn[i,])
		dir = unique(sign(dpsi[i,sgn[i,]]))
		if(length(dir)==1)
			r$dir[i] = ifelse(dir>0,'u','d')
		else if(length(dir)==0)
			r$dir[i] = 'n'
		else
			r$dir[i] = 'b'
		r$mean.psi[i] = mean(psi[i,!(m$tissue %in% colnames(sgn)[sgn[i,]])],na.rm=T)
	}
	r
}

getRandomCurve = function(x,y,sd,ylim=NULL,df=5){
	y = rnorm(length(y),y,sd)
	xx = seq(min(x),max(x),length.out=100)
	r=predict(smooth.spline(x,y,df=df),xx)
	if(!is.null(ylim))
		r$y = pmax(pmin(r$y,ylim[2]),ylim[1])
	r
}

plotDevAS4Pattern = function(p,main=p,s='mouse',t='brain',sites='ad',...){
	m = meta.tsm
	col = params$tissue.col
	col[names(col)!=t] = paste0(col[names(col)!=t],'40')
	m$col = col[m$tissue]
	f = anns[[s]]$sites==sites & age.segs[[s]][,t] == p
	y = apply(psi.tsm[[s]][f,],2,mean,na.rm=T)
	y = (y-min(y))/(max(y)-min(y))
	plotTissueAgeProile(y,m,age.axis = 'rank',bty='n',xlab='',ylab='',yaxt='n',plot.xaxt = F,main=paste0(main,' (',sum(f),')'),...)
	axis(1,labels = NA)
	axis(2,labels = NA)
}

plot5A = function(l){
	par(mar=c(0,0,0,0))
	t = round(sweep(all02,2,all02['total',],'/')*100,0)
	t = setNames(paste0(apply(t,1,min),'-',apply(t,1,max)),rownames(t))
	text = paste(t[c('u0','d1','u1','d0')],'%')
	plot(1,t='n',axes=F,xlim=c(0,100),ylim=c(0,110),xlab='',ylab='')
	x = c(15,60)
	y = c(15,65)
	w = 35
	segments(c(x[1],x[1]),y,c(x[1],x[1]),y+w)
	segments(c(x[2],x[2]),y,c(x[2],x[2]),y+w)
	
	segments(x,c(y[1],y[1]),x+w,c(y[1],y[1]))
	segments(x,c(y[2],y[2]),x+w,c(y[2],y[2]))
	text(50,0,'Development',adj=c(0.5,0),cex=1.5)
	text(5,50,'PSI',adj=c(0,0.5),srt=90,cex=1.5)
	text(rep(c(x+w/2),times=2),c(y[2],y[2],y[1],y[1])-1,text,adj=c(0.5,1),cex=1.5,font=2)
	text(rep(c(x+w/2),times=2),c(y[2],y[2],y[1],y[1])+w+5,c('Inclusion','Exclusion','',''),adj=c(0.5,0),cex=1.2,font=2)

		set.seed(1551)
	
	xx = seq(x[1]+2,x[1]+w-2,length.out=20)
	yflat = rep(y[2]+3,20)
	yup = (xx-min(xx))/(max(xx)-min(xx))
	ydw = 2+exp(-6*yup)*(w-4)
	yup = 2+(1-exp(-6*yup))*(w-4)
	sd=3
	
	for(c in params$tissue.col[-2:-2])	lines(getRandomCurve(xx,yflat,sd,c(y[2],y[2]+w)),lwd=2,col=c)
	lines(getRandomCurve(xx,yup+y[2],2),lwd=3,col=params$tissue.col['brain'])
	
	for(c in params$tissue.col[-2:-2])	lines(getRandomCurve(xx-x[1]+x[2],yflat+w-6,sd,c(y[2],y[2]+w)),lwd=2,col=c)
	lines(getRandomCurve(xx-x[1]+x[2],ydw+y[2],2),lwd=3,col=params$tissue.col['brain'])
	
	for(c in params$tissue.col[-2:-2])	lines(getRandomCurve(xx,yflat-y[2]+y[1]+w-6,sd,c(y[1],y[1]+w)),lwd=2,col=c)
	lines(getRandomCurve(xx,yup+y[1],2),lwd=3,col=params$tissue.col['brain'])
	
	for(c in params$tissue.col[-2:-2])	lines(getRandomCurve(xx-x[1]+x[2],yflat-y[2]+y[1],sd,c(y[1],y[1]+w)),lwd=2,col=c)
	lines(getRandomCurve(xx-x[1]+x[2],ydw+y[1],2),lwd=3,col=params$tissue.col['brain'])
	
	plotPanelLetter(l,lab.cex)
}

plotASSegStat = function(x,y,pv,col,xax,lty,skip.for.pv=1:2,pv.thr=c(' '=1,' '=0.05,'*'=0.01,'**'=0.001,'***'=0),first2horiz=TRUE,...){
	ylim = range(y,na.rm=TRUE)
	ylim[2] = ylim[2]+(ylim[2]-ylim[1])*0.05
	plot(x,y[,2],col=col,xlab='',ylim=ylim,xaxt='n',...)
	arrows(x,y[,1],x,y[,3],col=cols,lty=1,angle=90,code=3,length=0.03,xpd=T)
	if(first2horiz){
		abline(h=y[1,2],lty=3,col=col[1])
		abline(h=y[2,2],lty=3,col=col[2])
	}
	axis(1,xax,FALSE)
	text(xax,grconvertY(-0.05,'npc','user'),names(xax),adj=c(0,0),srt=-45,xpd=NA)
	
	
	pv.thr = sort(pv.thr)
	pv.sym = sapply(pv,function(x)names(pv.thr)[findInterval(x,pv.thr,all.inside = T)])
	d = (ylim[2]-ylim[1])*0.05/3
	if(!is.null(skip.for.pv)){
		pv.at = x[-skip.for.pv]
		y = y[-skip.for.pv,]
	}else
		pv.at = x
	for(i in 1:7)
		if(pv[i] < max(pv.thr[pv.thr<1])){
			segments(pv.at[i*2-1]  ,ylim[2]-d,pv.at[i*2-1]  ,y[i*2-1,3]+d)
			segments(pv.at[i*2-1]  ,ylim[2]-d,pv.at[i*2],ylim[2]-d)
			segments(pv.at[i*2],ylim[2]-d,pv.at[i*2],y[i*2,3]+d)
			text(pv.at[i*2-1]/2+pv.at[i*2]/2,ylim[2]-d,pv.sym[i],adj=c(0.5,0))
		}
}

plotMirroredMotFreq = function(f,sign,m,tissue,no.change.in.tiss=NULL,main='',leg=c('up','down'),plot.leg=TRUE,col=c('#FF000077','#0000FF77')){
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
	b=barplot(c(p1u,p1eu,p1ed,p1d),main=main,space = space,border=NA,ylim=range(0,p1u,p1eu,p1ed,p1d,p2u,p2eu,p2ed,p2d),col=col[1],ylab='fraction of exons with the motif')
	barplot(c(p2u,p2eu,p2ed,p2d),add=T,col=col[2],space = space,border=NA,ylab='fraction of exons with the motif')
	if(plot.leg)
		legend('topright',fill=col,legend = leg,bty='n',border=NA)
	segments(b[1],-0.05*ymax,b[50],-0.05*ymax)
	rect((b[19]+b[20])/2,-0.09*ymax,(b[29]+b[30])/2,-0.01*ymax,col='red')
	
	exon = c((b[1]+b[19])/2,(b[19]+b[20])/2,(b[29]+b[30])/2,(b[31]+b[50])/2)
	y = grconvertY(c(-0.07,-0.04,-0.01),'npc','user')
	segments(min(b),y[2],max(b),y[2],xpd=T)
	rect(exon[2],y[1],exon[3],y[3],col='gray',border=NA,xpd=T)
	ysh = grconvertY(c(0,-0.09),'npc','line')
	axis(1,exon,c("-100","acc.","don.","+100"),las=1,mgp=c(1,0.1,0)+ysh[1]-ysh[2])
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



plotExonScheme = function(up=TRUE){
	y = grconvertY(0.9,'npc','user')
	x = grconvertX(c(0.55,0.6,0.7,0.75,0.85,0.9),'npc','user')
	h = grconvertY(grconvertX(x[1:2],'user','inch'),'inch','user')
	h = h[2]-h[1]
	w = x[2]-x[1]
	segments(x[1],y,x[6],y)
	rect(x[1],y-h/2,x[2],y+h/2,col='gray',border='NA')
	rect(x[3],y-h/2,x[4],y+h/2,col='orange',border='NA')
	rect(x[5],y-h/2,x[6],y+h/2,col='gray',border='NA')
	if(up){
		arrows(x[3]-w/2,y+1.1*h,x[3]-w/2,y+0.2*h,col='red',angle=45,length = 0.05,lwd=3)
		#text(x[3]-0.1*w,y-0.1*h,'6mer',cex=0.7,adj=c(1,1))
	}
	else{
		arrows(x[4]+w/2,y+1.1*h,x[4]+w/2,y+0.2*h,col='red',angle=45,length = 0.05,lwd=3)
		#text(x[4]+0.1*w,y-0.1*h,'6mer',cex=0.7,adj=c(0,1))
	}
}

plot5G3.old = function(){
	w = 10
	plot(1,t='n',axes=F,xlim=c(0,19)*w,ylim=c(0,12)*w,xlab='',ylab='')
	x = c(0,1,5,6,10,11,16)*w
	y = c(10,6)*w
	segments(rep(x[1],2),y,rep(x[6],2),y)
	rect(rep(x[1],2),y+w/2,rep(x[2],2),y-w/2,col='gray',border=NA)
	rect(rep(x[5],2),y+w/2,rep(x[6],2),y-w/2,col='gray',border=NA)
	rect(rep(x[3],2),y+w/2,rep(x[4],2),y-w/2,col='orange',border=NA)
	
	rect(rep(x[7],2),y+w/2,rep(x[7]+w,2),y-w/2,col='gray',border=NA)
	rect(rep(x[7]+w,2),y[1]+w/2,rep(x[7]+2*w,2),y[1]-w/2,col='orange',border=NA)
	rect(rep(x[7]+w,2),y[2]+w/2,rep(x[7]+2*w,2),y[2]-w/2,col='gray',border=NA)
	rect(rep(x[7]+2*w,2),y[1]+w/2,rep(x[7]+3*w,2),y[1]-w/2,col='gray',border=NA)
	
	sh = 0.5
	arrows(x[6]+w*sh,y[1],x[7]-w*sh,y[1],lwd=2.5,col=params$tissue.col['brain'],angle = 15,length=0.1)
	arrows(x[6]+w*sh,y[2],x[7]-w*sh,y[2],lwd=2.5,col=params$tissue.col['brain'],angle = 15,length=0.1)
	
	arrows(x[6]+w*sh,y[1]-w*sh,x[7]-w*sh,y[2]+w*sh,lwd=2.5,col=params$tissue.col['heart'],angle = 15,length=0.1)
	arrows(x[6]+w*sh,y[2]+w*sh,x[7]-w*sh,y[1]-w*sh,lwd=2.5,col=params$tissue.col['heart'],angle = 15,length=0.1)
	
	plotOval(x[3]-1.4*w,y[1]+w/2,w*1.4,w/2,col='#7FC97F',border=NA)
	text(x[3]-1.4*w,y[1]+w/2,'QKI')
	text(x[3]-w/10,y[1]-w/10,'ACTAAC',adj=c(1,1),cex=0.7)
	
	plotOval(x[4]+1.4*w,y[2]+w/2,w*1.4,w/2,col='#7FC97F',border=NA)
	text(x[4]+1.4*w,y[2]+w/2,'QKI')
	text(x[4]+w/10,y[2]-w/10,'ACTAAC',adj=c(0,1),cex=0.7)
	
	text(x[7],y[1]+0.6*w,'Inclusion',adj=c(0,0),xpd=NA)
	text(x[7],y[2]+0.6*w,'Exclusion',adj=c(0,0),xpd=NA)
	plotPNG("figures/paper.figures/5/icons/brain.png",0.35,0.3,0.25)
	plotPNG("figures/paper.figures/5/icons/heart.png",0.65,0.3,0.2)
}

plot5G3 = function(){
	w = 10
	plot(1,t='n',axes=F,xlim=c(0,25.3)*w,ylim=c(0,12)*w,xlab='',ylab='')
	x = c(0,1,5,6,10,11,15)*w
	y = c(10,6)*w
	segments(rep(x[1],2),y,rep(x[6],2),y)
	rect(rep(x[1],2),y+w/2,rep(x[2],2),y-w/2,col='gray',border=NA)
	rect(rep(x[5],2),y+w/2,rep(x[6],2),y-w/2,col='gray',border=NA)
	rect(rep(x[3],2),y+w/2,rep(x[4],2),y-w/2,col='orange',border=NA)
	
	
	sh = 0.5
	arrows(x[6]+w*sh,y[1],x[7]-w*sh,y[1],lwd=2.5,angle = 15,length=0.1)
	arrows(x[6]+w*sh,y[2],x[7]-w*sh,y[2],lwd=2.5,angle = 15,length=0.1)
	

	plotOval(x[3]-1.4*w,y[1]+w/2,w*1.4,w/2,col='#7FC97F',border=NA)
	text(x[3]-1.4*w,y[1]+w/2,'QKI')
	text(x[3]-w/10,y[1]-w/10,'ACTAAC',adj=c(1,1),cex=0.6)
	
	plotOval(x[4]+1.4*w,y[2]+w/2,w*1.4,w/2,col='#7FC97F',border=NA)
	text(x[4]+1.4*w,y[2]+w/2,'QKI')
	text(x[4]+w/10,y[2]-w/10,'ACTAAC',adj=c(0,1),cex=0.6)
	
	# text(x[7],y[1]+0.6*w,'Inclusion',adj=c(0,0),xpd=NA)
	# text(x[7],y[2]+0.6*w,'Exclusion',adj=c(0,0),xpd=NA)
	plotPNG("figures/paper.figures/5/icons/brain.png",0.3,0.25,0.25)
	plotPNG("figures/paper.figures/5/icons/heart.png",0.7,0.25,0.2)
	
	xqki = grconvertX(0.52,'npc','user')
	yqki = grconvertX(c(0.09,0.17),'npc','user')
	arrows(xqki-w/2,yqki[2],xqki-w/2,yqki[1],col=params$tissue.col['brain'],lwd=2.5,angle = 15,length=0.1)
	arrows(xqki+w/2,yqki[1],xqki+w/2,yqki[2],col=params$tissue.col['heart'],lwd=2.5,angle = 15,length=0.1)
	text(xqki,yqki[1],'Likely\nQKI activity',adj=c(0.5,1.2))
	
	#make plots for QKI regulated exons
	ww = grconvertX(c(0,4*w),'user','inches')
	h = grconvertY(ww,'inches','user')
	h = h[2]-h[1]
	
	hline = grconvertY(c(0,1),'lines','user')
	hline = hline[2]-hline[1]
	#browser()
	xs = c(x[7]+c(0,4,5,9)*w) + 1.3*w #grconvertX(c(x[8]+c(0,3,5,8)*w),'user','ndc')
	ys = c(y[1]+c(-h/2,h/2),y[2]+c(-h/2,h/2)) #grconvertY(c(y[1]+c(-h/2,h/2+hline),y[2]+c(-h/2,h/2+hline)),'user','ndc')
	linex = seq(-pi/2,pi/2,length.out = 100)
	

	f = function(xc,yc,x,y,main,fr=0.9,...){
		x = scaleTo(x,xc[1],xc[2],fr)
		y = scaleTo(y,yc[1],yc[2],fr)
		segments(xc[1],yc[1],xc[2],yc[1])
		segments(xc[1],yc[1],xc[1],yc[2])
		lines(x,y,...)
		text(mean(xc),yc[2],main,adj=c(0.5,-0.1))
	}
	
	f(xs[1:2],ys[1:2],linex,sin(linex),main='Brain',lwd=3,col=params$tissue.col['brain'])
	f(xs[3:4],ys[1:2],linex,-sin(linex),main='Heart',lwd=3,col=params$tissue.col['heart'])
	f(xs[1:2],ys[3:4],linex,-sin(linex),main='Brain',lwd=3,col=params$tissue.col['brain'])
	f(xs[3:4],ys[3:4],linex,sin(linex),main='Heart',lwd=3,col=params$tissue.col['heart'])
	text(mean(xs[2:3]),ys[1],'Development',adj=c(0.5,1.2))
	text(mean(xs[2:3]),ys[3],'Development',adj=c(0.5,1.2))
	text(xs[1],mean(ys[1:2]),'PSI',adj=c(0.5,-0.2),srt=90)
	text(xs[1],mean(ys[3:4]),'PSI',adj=c(0.5,-0.2),srt=90)
}

# _6 ####
plotPhast = function(phast,ss,col=c('orange','olivedrab3','black'),ylab='human 46 primates PhastCons',...){
	ph = lapply(ss,function(s)t(sapply(phast[intersect(s,names(phast))],function(x)c(x[1:200],rep(NA,30),x[length(x)-(199:0)]))))
	ph.p = lapply(ph,function(p)cbind(apply(p,2,mean,na.rm=T),apply(p,2,sd,na.rm=T)/sqrt(nrow(p))))
	x = 1:nrow(ph.p$mi)
	
	plotArea(x,ph.p[[1]],col=col[1],xaxt='n',new = T,ylab=ylab,...)
	for(i in 2:length(ss))
		plotArea(x,ph.p[[i]],col=col[i],xaxt='n',new = F)
	rect(200,0,230,1,col = 'white',border = NA)
	
	
	y = grconvertY(c(0.07,0.04,0.01),'npc','user')
	segments(0,y[2],430,y[2],xpd=T)
	rect(200,y[1],230,y[3],col='gray',border=NA,xpd=T)
	
	axis(1,c(100,200,230,330),c('-100','acc.','don.','+100'),las=2)
	legend('topright',fill=col,legend = paste0(names(ss),' (',sapply(ss,length),')'),bty='n')
}

# _7 #####
getNewExonDevPref = function(s,t,thr=0){
	sps = c(species$short,'hq','mr','mrb','hqmrb')
	birth = meta.tsm$days[meta.tsm$species==s & meta.tsm$tissue==t & meta.tsm$stage== age.al.i[10,s]]
	bb = meta.tsm$stage[meta.tsm$species==s & meta.tsm$tissue==t & meta.tsm$days <= birth]
	ab = meta.tsm$stage[meta.tsm$species==s & meta.tsm$tissue==t & meta.tsm$days > birth]
	sids = born.seg.ids[nbf & sp.birth %in% sps[grep(species[s,'short'],sps)],s]
	psi = born.exn.tsm[[s]][sids,]
	
	bbp = apply(psi[,colnames(psi) %in% paste(s,t,bb)],1,mean,na.rm=TRUE)
	abp = apply(psi[,colnames(psi) %in% paste(s,t,ab)],1,mean,na.rm=TRUE)
	cnt = c(before=sum(bbp > abp + thr,na.rm=T),after=sum(sum(bbp + thr < abp,na.rm=T)))
	c(my.binom.test(cnt),cnt)
}



plotNewExonDevPref = function(t,sps=c('human','mouse','rat','rabbit'),plot.species.icons=FALSE,...){
	v1 = sapply(sps,getNewExonDevPref,t=t,thr=0.1)
	v2 = 1-v1
	b=barplot(rbind(v1[1,],v2[1,])*100,col=c(paste0(params$tissue.col[t],'55'),params$tissue.col[t]),border=NA,beside = T,ylim=range(0,v1[1:3,],v2[1:3,])*100,...)
	arrows(b,rbind(v1[2,],v2[2,])*100,b,rbind(v1[3,],v2[3,])*100,angle=90,code=3,length=0.03)
	if(plot.species.icons){
		b = grconvertX(apply(b,2,mean),'user','npc')
		for(s in 1:ncol(v1))
			plotPNG(paste0('figures/paper.figures/5/icons/',colnames(v1)[s],'.png'),b[s],-0.25,0.15)
	}
	p = binom.test(c(sum(v1[4,]),sum(v1[5,])))$p.value
	
	p=gsub("^-0+","-",strsplit(format(p,scientific = T,digits = 1),'e')[[1]],perl = T)
	p = bquote(italic("P")~ "=" ~ .(p[1])  %*% "10"^~.(p[2]))
	text(grconvertX(1,'npc','user'),grconvertY(1.05,'npc','user'),p,adj=c(1,0),xpd=TRUE,cex=0.7)
	invisible(list(v1=v1,v2=v2))
}

plot7A = function(l,nb.sp,alt.sp,nb.text,alt.text){
	plot(phyl.tree,x.lim=c(0,312*5),show.tip.label=F)
	text(grconvertX(0.05,'nfc','user'),grconvertY(0,'nfc','user'),'Million years ago',cex=1,adj=c(0,-4.8),xpd=NA)
	at = c(300,180,90,25)
	axis(1,312-at,at,las=2)
	grconvertY(0:1,'npc','user')
	
	x = grconvertX(312,'user','npc')
	par(xpd=NA)
	nbx  = grconvertX(c(0.40,0.67),'npc','user')
	altx = grconvertX(c(0.73,1),'npc','user')
	for(i in 1:7){
		plotPNG(paste0('figures/paper.figures/5/icons/',phyl.tree$tip.label[i],'.png'),x+0.08,grconvertY(8-i,'user','npc'),0.15)
		plotCE(nbx[1] ,8-i-0.3,nbx[2] ,8-i+0.3,lty=c(1,2),has.exon = grepl(species[phyl.tree$tip.label[i],'short'],nb.sp))
		plotCE(altx[1],8-i-0.3,altx[2],8-i+0.3,lty=c(ifelse(grepl(species[phyl.tree$tip.label[i],'short'],alt.sp),2,0),1),has.exon = T)
	}
	
	#text(312-90,grconvertY(0,'nfc','user'),'Million years',adj=c(0.5,0))
	# text(mean(nbx),7.5,nb.text,adj=c(0.5,0),cex=1)
	# text(mean(altx),7.5,alt.text,adj=c(0.5,0),cex=1)
	
	lw=grconvertY(0:1,'line','user')
	
	text(mean(nbx) ,seq(grconvertY(0.95,'nfc','user'),by=(lw[1]-lw[2]),length.out=length(nb.text)) ,nb.text,xpd=NA,adj=c(0.5,0.5),cex=1)
	text(mean(altx),seq(grconvertY(0.95,'nfc','user'),by=(lw[1]-lw[2]),length.out=length(alt.text)),alt.text,xpd=NA,adj=c(0.5,0.5),cex=1)
	plotPanelLetter(l,lab.cex)
}


getMaxStage = function(psi,ts,s,age.al,FUN=max){
	# select only stages were all tissues are present accourding to colnames(psi)
	t2t = unique(meta$tissue)
	t2t = setNames(t2t,substr(t2t,1,1))
	f = rep(T,nrow(age.al))
	ts = t2t[strsplit(ts,'')[[1]]]
	for(t in ts)
		f = f & paste(s,t,age.al[,s]) %in% colnames(psi)
	psi = psi[,paste(s,rep(ts,each=sum(f)), rep(age.al[f,s],times=length(ts)))]
	
	max.ts=apply(psi,1,function(x){
		m = FUN(x,na.rm=T)
		m=which(x==m)
		if(length(m)!=1)
			return(NA)
		if(length(m)>1)
			m = sample(m,1)
		colnames(psi)[m]
	})
	list(max.stages=data.frame(max.ts=max.ts,na.cnt=apply(is.na(psi),1,sum)),species=s,tissues=ts,stages=age.al[f,s],all=colnames(psi))
}

getMaxStageStat = function(m,bids,na.thr=Inf){
	max.stages = m$max.stages
	max.stages = max.stages[max.stages$na.cnt <= na.thr,]
	bids = bids[!is.na(bids[,m$species]),]
	rownames(bids) = bids[,m$species]
	bids = bids[,rownames(species)]
	max.stages = max.stages[bids[,m$species],]
	sb = apply(!is.na(bids),1,function(x)paste(species$short[x],collapse=''))
	st = strsplit(max.stages$max.ts,' ')
	r=list(tissue=table(sb,factor(sapply(st,'[',2),levels=m$tissues)),stage=table(sb,factor(sapply(st,'[',3),levels=m$stages)),all = table(sb,factor(max.stages$max.ts,levels=m$all)))
}

getMinStageStat = function(m,alt.sp,na.thr=Inf){
	max.stages = m$max.stages
	max.stages = max.stages[max.stages$na.cnt <= na.thr & rownames(max.stages) %in% names(alt.sp),]
	alt.sp = alt.sp[rownames(max.stages)]
	st = strsplit(max.stages$max.ts,' ')
	r=list(tissue=table(alt.sp,factor(sapply(st,'[',2),levels=m$tissues)),
				 stage= table(alt.sp,factor(sapply(st,'[',3),levels=m$stages)),
				 all =  table(alt.sp,factor(max.stages$max.ts,levels=m$all)))
}


plot7B = function(d,t,sps=c('human','mouse','rat','rabbit'),plot.species.icons=FALSE,...){
	v1 = sapply(sps,function(s){
		cnt = c(tissue=d[[s]]$tissue[species[s,'short'],t],other=sum(d[[s]]$tissue[species[s,'short'],]) - d[[s]]$tissue[species[s,'short'],t])
		c(100*my.binom.test(cnt[1],cnt[2]),cnt)
	})
	v2 = sapply(sps,function(s){
		cnt=c(tissue=d[[s]]$tissue['hqmrb',t],other=sum(d[[s]]$tissue['hqmrb',]) - d[[s]]$tissue['hqmrb',t])
		c(100*my.binom.test(cnt[1],cnt[2]),cnt)
	})
	b=barplot(rbind(v1[1,],v2[1,]),col=c(paste0(params$tissue.col[t],'55'),params$tissue.col[t]),border=NA,beside = T,ylim=range(0,v1[1:3,],v2[1:3,]),...)
	arrows(b,rbind(v1[2,],v2[2,]),b,rbind(v1[3,],v2[3,]),angle=90,code=3,length=0.03)
	if(plot.species.icons){
		b = grconvertX(apply(b,2,mean),'user','npc')
		for(s in 1:ncol(v1))
			plotPNG(paste0('figures/paper.figures/5/icons/',colnames(v1)[s],'.png'),b[s],-0.21,0.14)
	}
	p = prop.test(c(sum(v1[4,]),sum(v2[4,])),c(sum(v1[4:5,]),sum(v2[4:5,])))$p.value
	
	p=gsub("^-0+","-",strsplit(format(p,scientific = T,digits = 1),'e')[[1]],perl = T)
	p = bquote(italic("P")~ "=" ~ .(p[1])  %*% "10"^~.(p[2]))
	text(grconvertX(1,'npc','user'),grconvertY(1.05,'npc','user'),p,adj=c(1,0),xpd=TRUE,cex=0.7)
	invisible(list(v1=v1,v2=v2))
}


plotCE = function(x1,y1,x2,y2,col=c('gray','orange'),lty,has.exon){
	#rect(x1,y1,x2,y2)
	h = (y2-y1)/4
	w = (x2-x1)/10
	#segments(x1+w/1,y1+h,x2-w/1,y1+h)
	rect(x1,y1+0.2*h,x1+2*w,y1+1.8*h,col=col[1],border=NA)
	rect(x2,y1+0.2*h,x2-2*w,y1+1.8*h,col=col[1],border=NA)
	if(lty[1] != 0)
		lines(c(x1+2*w,(x1+x2)/2,x2-2*w),c(y1+h,y1+3.8*h,y1+h),lty=lty[1])
	
	if(has.exon){
		rect(x1+4.5*w,y1+0.2*h,x1+5.5*w,y1+1.8*h,col=col[2],border=NA)
		lines(c(x1+2*w,x1+6.5/2*w,x1+4.5*w),c(y1+h,y1+1.8*h,y1+h),lty=lty[2])
		lines(c(x1+2*w,x1+6.5/2*w,x1+4.5*w)+3.5*w,c(y1+h,y1+1.8*h,y1+h),lty=lty[2])
	}
}

plotCE1 = function(x1,y1,x2,y2,col=c('gray','orange'),has.exon=0){
	#rect(x1,y1,x2,y2)
	h = (y2-y1)/4
	w = (x2-x1)/10
	segments(x1+w/1,y1+h,x2-w/1,y1+h)
	segments(x1+w/1,y1+3*h,x2-w/1,y1+3*h)
	rect(x1,y1+0.2*h,x1+2*w,y1+1.8*h,col=col[1],border=NA)
	rect(x1,y1+2.2*h,x1+2*w,y1+3.8*h,col=col[1],border=NA)
	
	rect(x2,y1+0.2*h,x2-2*w,y1+1.8*h,col=col[1],border=NA)
	rect(x2,y1+2.2*h,x2-2*w,y1+3.8*h,col=col[1],border=NA)
	if(has.exon > 0)
		rect(x1+4.5*w,y1+2.2*h,x1+5.5*w,y1+3.8*h,col=col[2],border=NA)
	if(has.exon == 2)
		rect(x1+4.5*w,y1+0.2*h,x1+5.5*w,y1+1.8*h,col=col[2],border=NA)
}

# suppl ####
# _2 ####
plotSegDef = function(l,cols,sp,wd){
	sites = data.frame(pos=c(10,15,18,26,31,55,58,63,85,90),
										 sites=c('a','d','d','a','d','a','a','d','a','d'))
	top=90
	ys = seq(top-wd/2,by=-(wd+sp),length.out = 5)
	ys[length(ys)] = ys[length(ys)] - sp
	
	plot(1,t='n',bty='n',xaxt='n',yaxt='n',xlab='',ylab='',main='',xlim=c(0,100),ylim=c(0,100),xaxs = 'i',yaxs = 'i')
	segments(sites$pos[1],ys,sites$pos[nrow(sites)],ys)
	#const
	rect(sites$pos[1],ys-wd/2,sites$pos[2],ys+wd/2,col=cols$const,border = cols$const)
	rect(sites$pos[7],ys-wd/2,sites$pos[8],ys+wd/2,col=cols$const,border = cols$const)
	rect(sites$pos[9],ys-wd/2,sites$pos[10],ys+wd/2,col=cols$const,border = cols$const)
	
	rect(sites$pos[4],ys[c(1,5)]-wd/2,sites$pos[5],ys[c(1,5)]+wd/2,col=cols$alt,border = cols$alt)
	rect(sites$pos[2],ys[c(2,5)]-wd/2,sites$pos[3],ys[c(2,5)]+wd/2,col=cols$alt,border = cols$alt)
	rect(sites$pos[6],ys[c(3,5)]-wd/2,sites$pos[7],ys[c(3,5)]+wd/2,col=cols$alt,border = cols$alt)
	rect(sites$pos[8],ys[c(4,5)]-wd/2,sites$pos[9],ys[c(4,5)]+wd/2,col=cols$alt,border = cols$alt)
	segments(sites$pos,ys[5]-wd/2,sites$pos,ys[1]+wd/2+sp,lty=3,col='gray')
	text(sites$pos,ys[1]+wd/2+sp,sites$sites,adj = c(0.5,0),cex=0.8)
	text(sites$pos[1],mean(ys[2:3]),"transcripts",srt=90,adj=c(0.5,-0.3))
	#seg.titles = c('DD','CE','Cnst. intron','AA','RI','Cnst. exon')
	seg.titles.x = c(mean(sites$pos[2:3]),mean(sites$pos[4:5]),mean(sites$pos[5:6]),mean(sites$pos[6:7]),mean(sites$pos[8:9]),mean(sites$pos[9:10]))
	#text(seg.titles.x,ys[length(ys)]-wd/2-.2*sp,seg.titles,adj=c(0.5,1))
	text(seg.titles.x,ys[length(ys)]-wd/2-.2*sp,c('Alternative donor','Cassette exon','Constitutive intron','Alternative acceptor','Retained intron','Constitutive exons'),adj=c(1,0.5),srt=90)
	#segments(seg.titles.x,ys[length(ys)]-wd/2-0.1*sp,seg.titles.x,ys[length(ys)]-wd/2-1.9*sp)
	plotPanelLetter(l)
}

plotEJread = function(x1,x2,l1,l2,y,h,...){
	plotArc(x1,x2,h,y.base = y,...)
	segments(x1,y,x1-l1,y,...)
	segments(x2,y,x2+l2,y,...)
}
plotSAJR.AS.Q = function(l,wd){
	plot(1,t='n',bty='n',xaxt='n',yaxt='n',xlab='',ylab='',main='',xlim=c(0,100),ylim=c(0,100),xaxs = 'i',yaxs = 'i')
	x = c(10,20,45,55,80,90)
	y = 71
	c = c(cols$const,cols$alt,cols$const)
	segments(x[1],y,x[6],y)
	rect(x[c(1,3,5)],y-wd/2,x[c(2,4,6)],y+wd/2,col=c,border = c)
	
	plotEJread(x[2],x[3],2,3,y+wd*0.7,wd*2,lwd=3)
	plotEJread(x[4],x[5],1,4,y+wd*0.7,wd*2,lwd=3)
	segments(x[3]+1,y+wd*0.9,x[3]+6,y+wd*0.9,lwd=3)
	plotEJread(x[2],x[5],3,2,y-wd*0.7,-wd*2,lwd=3)
	text(50,y+wd*3.1,'inclusion read (i)',adj=c(0.5,0))
	text(50,y-wd*2.9,'exclusion read (i)',adj=c(0.5,1))
	
	text(50,35,expression(psi == frac(i/(ls + lr - 1),i/(ls + lr - 1) + e/(lr - 1))),cex=1)
	text(50,10,c('ls is the length of the segment\nlr is the length of the reads'))
	plotPanelLetter(l)
}

# _10 ####
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

# _18 ######
plotS18A = function(){
	plot(meta.example$age.use,psi.example,pch=19,cex=2,col=params$tissue.col[t],bty='n',ylab='PSI',xlab='Age, log(dpc)',main='dPSI definition',xlim=c(2.3,5),ylim=c(0,1))
	lines(age1000,psi.pred1000,lwd=4,col=params$tissue.col[t])
	wmx = which.max(psi.pred$y)
	wmn = which.min(psi.pred$y)
	xto = 4.7
	segments(psi.pred$x[wmx],psi.pred$y[wmx],xto,psi.pred$y[wmx],lty=2)
	segments(psi.pred$x[wmn],psi.pred$y[wmn],xto,psi.pred$y[wmn],lty=2)
	arrows(xto,psi.pred$y[wmx]-0.01,xto,psi.pred$y[wmn]+0.01,code=3,length=0.1)
	text(xto,psi.pred$y[wmx]/2+psi.pred$y[wmn]/2,'dPSI',srt=90,adj=c(0.5,1.3),cex=1.4,font=2)
}

plotS18B = function(){
	d = age.dpsi[[s]][,'brain']
	br = 0:100/100
	hist(d[anns[[s]]$sites=='ad' & !is.na(per.tissue.age.qv[[s]][,'brain']) & per.tissue.age.qv[[s]][,'brain'] > 0.05],br,border =NA,col='#00000030',xlab='dPSI',main='')
	hist(d[anns[[s]]$sites=='ad' & !is.na(per.tissue.age.qv[[s]][,'brain']) & per.tissue.age.qv[[s]][,'brain'] < 0.05],br,border =NA,col=paste0(params$tissue.col[t],'90'),ylab='',xlab='',add=T)
	abline(v=0.2,col='red',lwd=3)
	legend('topright',fill=c('#00000030',paste0(params$tissue.col[t],'90')),border='NA',bty='n',legend=c('pv.adj > 0.05','pv.adj < 0.05'))
	text(0.6,600,"devAS",font=2,cex=2,col=params$tissue.col[t])
}

plotS18C = function(){
	plot(meta.example$age.use,psi.example,pch=19,cex=2,col=params$tissue.col[t],bty='n',ylab='PSI',xlab='Age, log(dpc)',main='Four patterns',xlim=c(2.3,5),ylim=c(0,1))
	up = psi.pred1000[-1] > psi.pred1000[-length(psi.pred1000)]
	up = c(up[1],up)
	lines(age1000[up],psi.pred1000[up],lwd=4,col='red')
	lines(age1000[!up],psi.pred1000[!up],lwd=4,col='blue')
	
	wmx = which.max(psi.pred$y)
	wmn = which.min(psi.pred$y)
	xto = 4.8
	
	segments(psi.pred$x[wmx],psi.pred$y[wmx],xto,psi.pred$y[wmx],lty=2,col='red')
	segments(psi.pred$x[wmn],psi.pred$y[wmn],xto,psi.pred$y[wmn],lty=2,col='red')
	arrows(xto,psi.pred$y[wmx]-0.01,xto,psi.pred$y[wmn]+0.01,code=3,length=0.1,col='red')
	text(xto+0.05,psi.pred$y[wmx]/2+psi.pred$y[wmn]/2,'up',adj=c(0,0.5),cex=1.4,col='red',font=2,xpd=T)
	
	
	xto = 4.7
	end = length(psi.pred$x)
	segments(psi.pred$x[wmx],psi.pred$y[wmx],xto,psi.pred$y[wmx],lty=3,col='blue')
	segments(psi.pred$x[end],psi.pred$y[end],xto,psi.pred$y[end],lty=3,col='blue')
	arrows(xto,psi.pred$y[wmx]-0.01,xto,psi.pred$y[end]+0.01,code=3,length=0.1,col='blue')
	text(xto-0.05,psi.pred$y[wmx]/2+psi.pred$y[wmn]/2,'down',adj=c(1,0.5),cex=1.4,col='blue',font=2)
	
	timing = devAS.change.stat[[s]][[t]][sid,]
	timing = timing[3:4]/timing[1:2]
	
	segments(timing[1],predict(spline.fit,timing[1])$y,timing[1],0,col='red')
	text(timing[1]+0.01,-0.01,'up_timing',adj=c(0,0),font=2,col='red')
	
	segments(timing[2],predict(spline.fit,timing[2])$y,timing[2],0,col='blue')
	text(timing[2]+0.01,-0.01,'down_timing',adj=c(0,0),font=2,col='blue')
}

plotS18D = function(){
	dstat = devAS.change.stat[[s]][[t]]
	dstat = dstat[anns[[s]][rownames(dstat),'sites']=='ad',]
	upfr = dstat[,1]/(dstat[,1]+dstat[,2])
	hist(upfr,0:20/20,col='gray',border=NA,xlab='up/(up+down)',main='')
	abline(v=0.3,col='blue',lwd=3)
	abline(v=0.7,col='red',lwd=3)
	y=800
	text(0.15,y,'down',font=2,col='blue')
	text(0.85,y,'up',font=2,col='red')
	text(0.5,y,'up-down\nor\ndown-up',font=2,col='black')
}


plotS18E = function(){
	dstat = devAS.change.stat[[s]][[t]]
	dstat = dstat[anns[[s]][rownames(dstat),'sites']=='ad',]
	upfr = dstat[,1]/(dstat[,1]+dstat[,2])
	
	timing = dstat[,3]/dstat[,1] - dstat[,4]/dstat[,2]
	timing = timing[upfr > 0.3 & upfr < 0.7]
	
	hist(timing,20,col='gray',border=NA,xlab='up_timing - down_timing',main='')
	abline(v=0,lwd=3)
	
	y=60
	text(median(timing[timing<0]),y,'up-down',font=2)
	text(median(timing[timing>0]),y,'down-up',font=2)
}

# Nature review ####
plotAsEventCount.boot = function(data,pchs,...){
	mean   = setNames(lapply(names(pchs),function(ss)sapply(data,function(x)sapply(x,function(y)mean(y[,ss])))),names(pchs))
	q0.025 = setNames(lapply(names(pchs),function(ss)sapply(data,function(x)sapply(x,function(y)quantile(y[,ss],p=0.025)))),names(pchs))
	q0.975 = setNames(lapply(names(pchs),function(ss)sapply(data,function(x)sapply(x,function(y)quantile(y[,ss],p=0.975)))),names(pchs))
	plotAsEventCount(mean,cil=q0.025,cih=q0.975,pchs=pchs,ylim=range(0,q0.975),...)	
}

plotSampleTable1 = function(meta_,sps.order,sspace=7,max.sams = 6,cex=0.7,...){
	for(s in rownames(species)){
		s2s = unique(meta[meta$species==s,c('paper.stages','stage')])
		s2s = setNames(s2s$paper.stages,s2s$stage)
		age.al.i[age.al.i[,s]!='',s] = s2s[age.al.i[age.al.i[,s]!='',s]]
	}
	
	meta_$stage = meta_$paper.stages
	sstages = lapply(rownames(species),function(s){
		t = meta_[meta_$species==s,]
		sort(rank(sapply(split(t$days,t$stage),mean)))
	})
	names(sstages) = rownames(species)
	age.al.i = age.al.i[age.al.i$mouse %in% names(sstages$mouse),]
	r = list()

	for(s in rownames(species)){
		r[[s]] = list()
		for(ms in 1:nrow(age.al.i)){
			c=age.al.i[ms,s]
			if(c==''){
				r[[s]][[ms]] = NULL
			}else if(ms == 1 || is.na(sstages[[s]][age.al.i[ms-1,s]])){
				r[[s]][[ms]] = sstages[[s]][sstages[[s]]<=sstages[[s]][c]]
			}else{
				r[[s]][[ms]] = sstages[[s]][sstages[[s]]<=sstages[[s]][c] & sstages[[s]]>sstages[[s]][age.al.i[ms-1,s]]]
			}
		}
		r[[s]][[ms+1]] = sstages[[s]][sstages[[s]]>sstages[[s]][c]]
	}
	xmax = length(sstages$mouse)*3+3
	plot(1,xlim=c(-0.2*xmax,xmax),ylim=c(length(sps.order)*length(unique(meta$tissue))+(length(sps.order)-1)*sspace+1,0),xlab='',ylab='',bty='n',xaxt='n',yaxt='n',t='n',yaxs='i',xaxs='i',...)
	l = 0
	abline(v=seq(0,by=3,length.out = nrow(age.al.i)+1),col='gray',lty=3)
	for(s in sps.order){
		stage.coors = c()
		plotPNG(paste0("figures/paper.figures/5/icons/",s,".png"),grconvertX(-0.1*xmax,'user','npc'),grconvertY(l+5,'user','npc'),0.09)
		for(t in sort(unique(meta$tissue))){
			text(-0.7,l+0.5,substr(t,1,1),adj=c(0.5,0.5),cex=0.7)
			for(ms in 0:nrow(age.al.i)){
				stages = r[[s]][[ms+1]]
				if(length(stages)==0) next
				for(stg in 1:length(stages)){
					sams = sum(meta_$species==s & meta_$tissue==t & meta_$stage == names(stages)[stg],na.rm = T)
					cellen = 3*max(1,sum(age.al.i[,s] %in% names(stages)))/length(stages)
					if(sams>0){
						col = makeTransparent(params$tissue.col[t],alpha = sams/max.sams)
						if(cellen>3)
							rect(ms*3+(stg-1)*cellen+1,l,ms*3+stg*cellen-1,l+1,col='white',border = NA)
						rect(ms*3+(stg-1)*cellen,l,ms*3+stg*cellen,l+1,col=col,border = NA)
						text(ms*3+(stg-0.5)*cellen,l+0.5,sams,adj = c(0.5,0.5),cex=cex)
						stage.coors[names(stages)[stg]] = ms*3+(stg-0.5)*cellen
					}
				}
			}
			l = l +1
		}
		text(stage.coors,l+0.2,names(stage.coors),adj = c(0,0.5),srt=-65,cex=ifelse(s==sps.order[length(sps.order)],1,cex),xpd=T)
		l = l + sspace
		abline(h=l-0.5,col='gray',lty=3)
	}
	#axis(1,(1:nrow(age.al.i))*3-1.5,labels = age.al.i$mouse,las=3)
}

plotSampleTable = function(exact.stage=TRUE,main='',cex=0.7){
	sams = list()
	for(s in rownames(species)){
		sams[[s]] = matrix(0,ncol=nrow(age.al.i),nrow=7,dimnames = list(sort(unique(meta$tissue)),age.al.i$mouse))
		for(t in rownames(sams[[s]]))
			for(a in 1:nrow(age.al.i)){
				if(age.al.i[a,s] == '') next
				# just matched samples, some could be counted twice, some could be not counted at all
				if(exact.stage)
					sams[[s]][t,a] = sum(meta$species == s & meta$tissue == t & meta$stage==age.al.i[a,s])
				else{
					days = max(meta$days[meta$species == s & meta$stage==age.al.i[a,s]])
					f = meta$tissue == t & meta$species == s & meta$days <= days
					if(a>1 && age.al.i[a-1,s] != ''){
						pdays = max(meta$days[meta$species == s & meta$stage==age.al.i[a-1,s]])
						f = f & meta$days > pdays
					}
					sams[[s]][t,a] = sum(f)
				}
			}
	}
	
	
	nr = sapply(sams,nrow)
	max = max(sapply(sams,max))
	plot(1,xlim=c(-3,nrow(age.al.i)),ylim=c(sum(nr)+length(nr)-1,0),main=main,xlab='Mouse stage',ylab='',bty='n',xaxt='n',yaxt='n',t='n')
	l = 0
	for(s in 1:length(sams)){
		plotPNG(paste0("figures/paper.figures/5/icons/",names(sams)[s],".png"),grconvertX(-1.5,'user','npc'),grconvertY(l+3.5,'user','npc'),0.1)
		for(t in 1:nrow(sams[[s]])){
			for(a in 1:ncol(sams[[s]])){
				if(sams[[s]][t,a]>0){
					col = makeTransparent(params$tissue.col[rownames(sams[[s]])[t]],alpha = sams[[s]][t,a]/max)
					rect(a-1,l,a,l+1,col=col,border = NA)
					text(a-0.5,l+0.5,sams[[s]][t,a],adj = c(0.5,0.5),cex=cex)
				}
			}
			l = l + 1
		}
		if(s != length(sams)) abline(h=l+0.5,col='gray',lty=3)
		l = l + 1
	}
	axis(1,1:nrow(age.al.i)-0.5,labels = age.al.i$mouse,las=3)
	invisible(sams)
}


getEClipSiteOccurence = function(f,ann){
	keclip1 = read.table(f,sep='\t')
	colnames(keclip1)[c(1,2,3,6)]=c('chr_id','start','stop','strand')
	keclip1$chr_id = gsub('chr','',keclip1$chr_id)
	table(keclip1$chr_id,keclip1$chr_id %in% ann$chr_id)
	keclip1 = keclip1[keclip1$chr_id %in% ann$chr_id,]
	keclip1[1:2,]
	
	#seqinf = Seqinfo(unique(ann$human$chr_id))
	
	gre = GRanges(keclip1$chr_id ,IRanges(keclip1$start,keclip1$stop),keclip1$strand)#,seqinfo = seqinf)
	gra = GRanges(ann$chr_id ,IRanges(ann$start,ann$stop),ifelse(ann$strand==1,'+','-'))#,seqinfo = seqinf)
	o = findOverlaps(gra,gre,type = 'any',maxgap = 200,ignore.strand=FALSE)
	o = cbind(o@from,o@to)

	sites = matrix(0,ncol=400,nrow=nrow(ann))
	for(i in 1:nrow(o)){
		a = ann[o[i,1],c('start','stop')]
		e = keclip1[o[i,2],c('start','stop')]
		#upstream
		if(e$start<a$start){
			inx = c(max(1,200 + e$start-a$start),min(200,200 + e$stop-a$start))
			if(inx[2]>=inx[1]){
				inx = inx[1]:inx[2]
				sites[o[i,1],inx] = sites[o[i,1],inx] + 1
			}
		}
		# downstream
		if(e$stop>a$stop){
			inx = 200 + c(max(1,e$start -a$stop),min(200,e$stop -a$stop))
			if(inx[2]>=inx[1]){
				inx = inx[1]:inx[2]
				sites[o[i,1],inx] = sites[o[i,1],inx] + 1
			}
		}
	}
	sites[ann$strand==-1,] = sites[ann$strand==-1,ncol(sites):1]
	sites
}


# additioanl ######
devASSameDir = function(d){
	r = matrix(NA,ncol=ncol(d),nrow=ncol(d),dimnames=list(colnames(d),colnames(d)))
	for(i in 1:(ncol(r)-1))
		for(j in (i+1):ncol(r)){
			f = d[,i] %in% c('u','d') & d[,j] %in% c('u','d')
			r[i,j] = mean(d[f,i] == d[f,j])
			r[j,i] = sum(f)
		}
	r
}

getDevASPattern = function(x,y,df,N=1000){
	stop('use getDevASPattern1: it just adds min/max and dPSI estimation')
	f = !is.na(x) & !is.na(y)
	x = x[f]
	y = y[f]
	if(length(unique(x))<= df) return(c(up=NA,down=NA,up.time=NA,down.time=NA))
	xx = seq(min(x),max(x),length.out = N)
	p = predict(smooth.spline(x,y,df=df),x=xx)$y
	dp = p[-1]-p[-N]
	dpx = dp*xx[-1]
	x = x[-1]
	c(up   = sum(dp[dp>0]),
		down = sum(-dp[dp<0]),
		up.time = sum(dpx[dp>0]),
		down.time = sum(-dpx[dp<0]))
}

getDevASPattern1 = function(x,y,df,N=1000){
	f = !is.na(x) & !is.na(y)
	x = x[f]
	y = y[f]
	if(length(unique(x))<= df | IQR(x) == 0) return(c(up=NA,down=NA,up.time=NA,down.time=NA,min=NA,max=NA))
	xx = seq(min(x),max(x),length.out = N)
	p = predict(smooth.spline(x,y,df=df),x=xx)$y
	dp = p[-1]-p[-N]
	dpx = dp*xx[-1]
	x = x[-1]
	c(up   = sum(dp[dp>0]),
		down = sum(-dp[dp<0]),
		up.time = sum(dpx[dp>0]),
		down.time = sum(-dpx[dp<0]),
		min=min(p),
		max=max(p))
}


rarefy2uniform = function(x,n,min.cnt=1,seed=NULL){
	if(!is.null(seed))
		set.seed(seed)
	bins = seq(from=min(x),to=max(x),length.out = n+1)
	x2bin = findInterval(x,bins,rightmost.closed = TRUE)
	cnt = table(x2bin)
	cnt = max(min(cnt[cnt>0]),min.cnt)
	f = rep(FALSE,length(x))
	for(i in 1:n){
		inx = which(x2bin == i)
		if(length(inx) > 1)
			inx = sample(inx,cnt)
		f[inx] = TRUE
	}
	f
}


getSpeciesCorForTissue = function(s1,s2,t,age.al){
	age.al = age.al[,c(s1,s2)]
	age.al = age.al[apply(age.al=='',1,sum)==0,]
	age.al = age.al[paste(s1,t,age.al[,s1]) %in% colnames(orth.seg.ad.tsm[[s1]]) & 
										paste(s2,t,age.al[,s2]) %in% colnames(orth.seg.ad.tsm[[s2]]),]
	if(nrow(age.al) < 3) return(NULL)
	psi1 = orth.seg.ad.tsm[[s1]][,paste(s1,t,age.al[,s1])]
	psi2 = orth.seg.ad.tsm[[s2]][,paste(s2,t,age.al[,s2])]
	cor = sapply(1:nrow(psi1),function(i){
		cor(psi1[i,],psi2[i,],u='p')
	})
	diff = apply(abs(psi1-psi2),1,mean,na.rm=T)
	r = data.frame(species1 = s1,species2=s2,tissue=t,sid1=rownames(psi1),sid2=rownames(psi2),cor=cor,diff=diff,nnna = apply(!is.na(psi1) & !is.na(psi2),1,sum))
	cmn1 = intersect(rownames(age.segs[[s1]]),r$sid1)
	cmn2 = intersect(rownames(age.segs[[s2]]),r$sid2)
	r$p2 = r$p1 = '-'
	r$dpsi2  = r$dpsi1 = NA
	r[match(cmn1,r$sid1),'p1'] = age.segs[[s1]][cmn1,t]
	r[match(cmn2,r$sid2),'p2'] = age.segs[[s2]][cmn2,t]
	r[match(cmn1,r$sid1),'dpsi1'] = age.dpsi[[s1]][cmn1,t]
	r[match(cmn2,r$sid2),'dpsi2'] = age.dpsi[[s2]][cmn2,t]
	r
}
