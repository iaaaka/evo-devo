getTestedSids = function(s,t,MEL=500,MAS=5,sites='ad',max.nn.frac=0.5){
	f=anns[[s]]$length<=MEL & anns[[s]]$alt.size<=MAS & anns[[s]]$sites %in% sites
	psi = psi.tsm[[s]][,meta.tsm[colnames(psi.tsm[[s]]),'tissue']==t]
	rownames(anns[[s]])[f & apply(is.na(psi),1,mean)<=max.nn.frac]
}

getASPropPerTissue = function(d,f,as,thr,m){
	d = d[f,] >= thr
	m = m[colnames(d),]
	as = as[f,]
	do.call(cbind,lapply(unique(m$tissue),function(t){
		a = as[,paste0(t,'.age.as')]
		s = as[,paste0(t,'.tested')]
		dd = d[,m$tissue==t]
		apply(dd,2,function(x){expressed=c(sum(x),age=mean(a[x]),tested=mean(s[x]))})
	}))
}

getASProp = function(d,f,as,thr,m){
	d = d[f,] >= thr
	m = m[colnames(d),]
	as = as[f]
	apply(d,2,function(x)mean(as[x]))
}


plotGeneASFreqOnExprThrs = function(exp,filter,gene.info,thrs,meta){
	for(t in thrs){
		x=getASPropPerTissue(exp,filter,gene.info,t,meta)
		plotTissueAgeProile(x[1,],meta,ylab='# of expressed genes',age.axis = 'rank',main=paste0('# genes, RPKM > ',t))
		plotTissueAgeProile(x[3,]*100,meta,ylab='% of tested',age.axis = 'rank',main='% of tested (this tissue)')
		plotTissueAgeProile(x[2,]*100,meta,ylab='% of age-AS',age.axis = 'rank',main='% of age-AS (this tissue)')
		plotTissueAgeProile(x[2,]/x[3,]*100,meta,ylab='% of age-AS (in tested)',age.axis = 'rank',main='% of age-AS (in tested)')
		
		x=getASProp(exp,filter,apply(gene.info[,paste0(unique(meta$tissue),'.tested')],1,sum)>0,t,meta)
		plotTissueAgeProile(x*100,meta,ylab='% of tested',age.axis = 'rank',main='% of tested (any tissue)')
		x=getASProp(exp,filter,apply(gene.info[,paste0(unique(meta$tissue),'.age.as')],1,sum)>0,t,meta)
		plotTissueAgeProile(x*100,meta,ylab='% of age-AS',age.axis = 'rank',main='% of age-AS (any tissue)')
	}
}

saveGO2xlsx = function(go,names,fname,qv.thr=0.1){
	file.remove(fname)
	for(t in names(go)) 
		for(n in names){
			tmp = go[[t]][[n]]
			tmp = tmp[tmp$qv.over<qv.thr,]
			tmp = tmp[order(tmp$ontology,tmp$pv.over),]
			write.xlsx2(tmp, file=fname, sheetName=paste(t,n), row.names=FALSE,append=TRUE)
		}
}

testGoByTIssues = function(tested.segs,meta,genome,max.na.freq=0.5){
	m = meta
	h = tested.segs
	human.go = list()
	for(t in unique(meta$tissue)){
		sel = intersect(rownames(h$seg),rownames(sp.tis.clust[[s]][[t]][[1]]))
		selu = sel[cor2agerank[[s]][sel,t]>0]
		seld = sel[cor2agerank[[s]][sel,t]<0]
		bkg = rownames(h$seg)[apply(is.na(psi.tsm[[s]][rownames(h$seg),meta.tsm[colnames(psi.tsm[[s]]),'tissue']==t]),1,mean)<=max.na.freq]
		bkg = log(apply(h$i[bkg,m$tissue==t] + h$e[bkg,m$tissue==t],1,sum))
		selu = unique(unlist(seg2ens[[s]][selu]))
		seld = unique(unlist(seg2ens[[s]][seld]))
		bkg = sapply(revList(seg2ens[[s]][names(bkg)]),function(x){max(bkg[x])})
		cat(t,length(selu),length(seld),length(bkg),'\n')
		human.go[[t]] = list(all=getGO(union(seld,selu),bkg,4,'Wallenius',genome=genome,id = 'ensGene',gene2cat = NULL),
													up=getGO(selu,bkg,4,'Wallenius',genome=genome,id = 'ensGene',gene2cat = NULL),
												  dw=getGO(seld,bkg,4,'Wallenius',genome=genome,id = 'ensGene',gene2cat = NULL))
		human.go[[t]]$up2dw = getGO(selu,union(seld,selu),4,'Hypergeometric',genome=genome,id = 'ensGene',gene2cat = NULL)
		human.go[[t]]$dw2up = getGO(seld,union(seld,selu),4,'Hypergeometric',genome=genome,id = 'ensGene',gene2cat = NULL)
	}
	human.go
}


exontSel = function(d,qv=0.2,n=NULL,or=1,or.fun=`>`){
	if(is.character(or.fun)){
		u = exontSel(d,qv,n,or,`>`)
		d = exontSel(d,qv,n,or,`<`)
		t = d[1,]
		t[1,] = rownames(t) = 0
		rbind(u,t,d)
	}
	else if(is.null(n))
		d[d$q.value<=qv & or.fun(d$odds.ratio,or),]
	else
		d[or.fun(d$odds.ratio,or),][1:n,]
}

calcExontEnrichmentWithSelfControl = function(test,contr,onth,remove.unannotated,exont,cod='c',Ncontr=2,min.size=4,add.parents=TRUE){
	exont.gr = GRanges(exont$chr.id,IRanges(exont$start,exont$end),ifelse(exont$strand==1,'+','-'))
	t = makeControlExons(test,contr,cod,Ncontr)
	a = rbind(t$test[,1:5],t$contr[,1:5])
	seg.gr = GRanges(a$chr_id,IRanges(a$start,a$stop),ifelse(is.na(a$strand),'*',ifelse(a$strand== 1,'+','-')))
	s2e = findOverlaps(exont.gr,seg.gr,maxgap=0,type='any',select='all',ignore.strand=FALSE)
	s2e = sapply(split(exont$exont.id[s2e@from],rownames(a)[s2e@to]),unique)
	exontEnrich(onth,rownames(t$test),rownames(t$contr),s2e,min.size=min.size,add.parents=add.parents,remove.unannotated=remove.unannotated)
}

makeControlExons = function(test,contr,cod,n=2){
	test = test[test$cod %in% cod,]
	test$seg.id = rownames(test)
	contr$seg.id = rownames(contr)
	contr = contr[contr$type=='EXN' & contr$gene_id %in% test$gene_id & contr$cod %in% cod,]
	contr$selected = FALSE
	test = split(test,test$gene_id)
	contr = split(contr,contr$gene_id)
	nalt = sapply(test,nrow)
	nctr = sapply(contr,nrow)
	cmn = intersect(names(nalt),names(nctr))
	cmn = cmn[nalt[cmn]*n <= nctr[cmn]]
	cat(sum(sapply(test,nrow))-sum(sapply(test[cmn],nrow)),'from',sum(sapply(test,nrow)),'exons do not have controls and were removed\n')
	test  = test[cmn]
	contr = contr[cmn]
	
	for(i in 1:length(test)){
		t = test[[i]]
		for(k in 1:n){
			for(j in 1:nrow(t)){
				f = which(!contr[[i]]$selected & contr[[i]]$length>=t$length[j])
				if(length(f)>0){#cut part of longest exon
					if(length(f)>1)
						id = sample(f,1)
					else
						id =f
					#if(i == 126) cat(1,id,paste(f,collapse=','),'\n')
					contr[[i]]$selected[id] = TRUE
					diff = contr[[i]]$length[id] - t$length[j]
					start.pos = sample(diff,1)
					contr[[i]]$start[id] = contr[[i]]$start[id] + start.pos
					contr[[i]]$stop[id] = contr[[i]]$stop[id] - diff + start.pos
					contr[[i]]$length[id] = contr[[i]]$stop[id] - contr[[i]]$start[id] + 1
				}else{#take closest by size
					f = which(!contr[[i]]$selected)
					id = f[order(t$length[j] - contr[[i]]$length[f])[1]]
					#	if(i == 126) cat(1,id,paste(f,collapse=','),'\n')
					contr[[i]]$selected[id] = TRUE
				}
			}
		}
	}
	contr=do.call(rbind,lapply(contr,function(x)x[x$selected,]))
	contr$selected = NULL
	test = do.call(rbind,test)
	rownames(test) = test$seg.id
	rownames(contr) = contr$seg.id
	test$seg.id = contr$seg.id = NULL
	list(test=test,contr=contr)
}

exontEnrich = function(o,s,b,s2e,min.size=0,add.parents=TRUE,remove.unannotated=TRUE){
	if(remove.unannotated){
		s = s[s %in% names(s2e)]
		b = b[b %in% names(s2e)]
	}
	b = setdiff(b,s)
	s2e = s2e[intersect(names(s2e),union(s,b))]
	if(add.parents)
		s2e = lapply(s2e,function(x)get_ancestors(o,x))
	r = data.frame(exont.id = unique(unlist(s2e)),sel.cnt=0,total.cnt=0,p.value=NA,q.value=NA,odds.ratio=NA)
	rownames(r) = r$exont.id
	sel = table(unlist(s2e[intersect(names(s2e),s)]))
	bkg = table(unlist(s2e[intersect(names(s2e),b)]))
	r[names(sel),'sel.cnt'] = sel
	r[names(bkg),'total.cnt'] = bkg
	r = r[r$total.cnt >=min.size,]
	tc = length(b)
	sc = length(s)
	for(i in 1:nrow(r)){
		cat('\r',i,nrow(r),'      ')
		t = fisher.test(cbind(c(tc-r$total.cnt[i],sc-r$sel.cnt[i]),c(r$total.cnt[i],r$sel.cnt[i])))
		r$p.value[i] = t$p.value
		r$odds.ratio[i] = t$estimate
	}
	r = r[r$exont.id !='EXONT:000001',] # there are terms that are not children of root...
	r$total.cnt = r$total.cnt + r$sel.cnt
	r$q.value = p.adjust(r$p.value,m='BH')
	r = r[order(r$p.value),]
	r$def = o$name[r$exont.id]
	r
}


plotGeneExprPropertiesForAS = function(ge.info,tissue,field,amp.thr,species,as.data,dir=c(pmax,`>`,`<`)[[1]],cod='c',n=10,max.nn.frac=0.5,pv.thr=0.05,...){
	all.ids = rownames(anns[[species]])[anns[[species]]$sites=='ad' & anns[[species]]$length <= MEL & anns[[species]]$alt.size<= MAS & anns[[species]]$cod %in% cod & apply(is.na(psi.tsm[[species]][,meta.tsm[colnames(psi.tsm[[species]]),'tissue']==tissue]),1,mean)<=max.nn.frac]
	sgn.ids = all.ids[per.tissue.age.qv[[species]][all.ids,tissue] < 0.05 & per.tissue.age.ampl[[species]][all.ids,tissue] >= amp.thr & as.logical(dir(cor2agerank[[species]][all.ids,tissue],0))]
	sgn.ens = unique(unlist(seg2ens[[species]][sgn.ids]))
	all.ens = unique(unlist(seg2ens[[species]][all.ids]))
	f = meta[colnames(as.data$i),'tissue'] == tissue
	cov = apply(as.data$i[all.ids,f]+as.data$e[all.ids,f],1,sum)
	ens.cov = sapply(revList(seg2ens[[species]][all.ids]),function(sids){max(cov[sids])})
	all.ens = intersect(all.ens,rownames(ge.info))
	sgn.ens = intersect(sgn.ens,rownames(ge.info))
	nsg.ens = setdiff(all.ens,sgn.ens)
	ens.cov.inx = number2bin(ens.cov,n)
	f = function(x){
		x = x[!is.na(x)]
		c(mean(x),sd(x)/sqrt(length(x)))
	}
	sgn.list = split(ge.info[sgn.ens,field],factor(ens.cov.inx[sgn.ens],levels = 1:n))
	nsg.list = split(ge.info[nsg.ens,field],factor(ens.cov.inx[nsg.ens],levels = 1:n))
	sgn.stat = sapply(sgn.list,f)
	nsg.stat = sapply(nsg.list,f)
	pv = sapply(1:n,function(i){
		if(sum(!is.na(sgn.list[[i]])) > 0 && sum(!is.na(nsg.list[[i]])) > 0)
			wilcox.test(sgn.list[[i]],nsg.list[[i]])$p.value
		else
			NA
	})
	qv = p.adjust(pv,m='BH') < pv.thr
	ylim = range(sgn.stat[1,],nsg.stat[1,],na.rm=TRUE)
	d = ylim[2] - ylim[1]
	ylim = c(ylim[1] - d*0.1,ylim[2] + d*0.1)
	plot(1:n,sgn.stat[1,],col='red',t='b',ylim=ylim,xlab='coverage bin',pch=ifelse(qv,19,1),...)
	lines(1:n,nsg.stat[1,],col='blue',t='b',pch=ifelse(qv,19,1))
	arrows(1:n,sgn.stat[1,]-sgn.stat[2,]/2,1:n,sgn.stat[1,]+sgn.stat[2,]/2,col='red',angle=90,code=3,length=0.1)
	arrows(1:n,nsg.stat[1,]-nsg.stat[2,]/2,1:n,nsg.stat[1,]+nsg.stat[2,]/2,col='blue',angle=90,code=3,length=0.1)
	text(1:n,rep(ylim[2],n),paste0(sapply(sgn.list,length),'\n',table(ens.cov.inx)),adj=c(0.5,1))
	mod = lm(ge.info[all.ens,field] ~ factor(ens.cov.inx[all.ens]) + factor(I(all.ens %in% sgn.ens)))
	mod.pv = anova(mod)[2,5]
	
	pv = pv[!is.na(pv)]
	ih.pv = dirwin.hall(sum(pv),length(pv))
	res = c(eff=mod$coefficients[length(mod$coefficients)],anova.pv=mod.pv,ih.pv=ih.pv)
	legend('bottomleft',legend=paste0(c('effect=','pv=','ih.pv='),format(res,digits=3,scientific=TRUE)))
	invisible(res)
}

revList = function(l){
	s = sapply(l,length)
	n = lapply(1:length(l),function(i){rep(names(l)[i],s[i])})
	split(unlist(n),unlist(l))
}


plotSNPfreqPerTissueWithCont = function(c,x,ylim=NULL,...){
	pv = p.adjust(sapply(1:length(x),function(i){wilcox.test(x[[i]],c[[i]],a='l')$p.value}),m='BH')<0.05
	mx = sapply(x,function(x){mean(x)})
	mc = sapply(c,function(x){mean(x)})
	
	sx = sapply(x,function(x){sd(x)/sqrt(length(x))})*2
	sc = sapply(c,function(x){sd(x)/sqrt(length(x))})*2
	if(is.null(ylim))
		ylim=range(mx-sx,mx+sx,mc-sc,mc+sc)
	plot(1,t='n',ylab='SNP log10(freq)',xlab='',xlim=c(1,length(x)*2),ylim=ylim,xaxt='n',...)
	for(i in 1:length(mx)){
		x = (i-1)*2 + 1
		points(x,mc[i],pch=19,col='black',cex=2)
		segments(x,mc[i] - sc[i],x,mc[i] + sc[i])
		col = params$tissue.col[names(mx)[i]]
		x = x + 1
		if(!pv[i]){
			col=col2rgb(col)
			col=rgb(col[1],col[2],col[3],alpha = 50,maxColorValue = 255)
		}
		points(x,mx[i],pch=ifelse(pv[i],19,1),col=col,cex=2)
		segments(x,mx[i] - sx[i],x,mx[i] + sx[i],col=col)
	}
	axis(1,1:length(mx)*2-0.5,names(mx),las=3)
}

plotSgnPropOnCov = function(s,n=20,ampl.thr=0.05){
	d = readRDS(paste0('Rdata/',s,'.as.u.filtered.Rdata'))
	d = d[d$seg$length <= MEL & anns[[s]]$alt.size <= MAS & d$seg$sites == 'ad',]
	m = meta[colnames(d$i),]
	d = calcMeanCols(d$i + d$e,m$tissue)
	a = per.tissue.age.ampl[[s]][rownames(d),]
	q = per.tissue.age.qv[[s]][rownames(d),] <= 0.05
	q[is.na(q) | is.na(a)] = FALSE
	tissues = unique(meta$tissue)
	f1 = function(x){x[,4]/x[,5]*100}

	
	tissue1 =  list()
	
	for(t in colnames(d)){
		tcr = factor(number2bin(d[,t],n))
		amp = sapply(split(a[q[,t],t],tcr[q[,t]]),quantile,probs=c(0.1,0.5,0.9))
		if(sum(q[,t])>0)
			plot(1:n,amp[2,],ylim=range(amp,na.rm=TRUE),xlab='coverage bin',ylab='PSI amplitude',main=paste0(t,' (',sum(q[,t]),')'),t='b')
		else
			plot(1,t='n',xlab='coverage bin',ylab='PSI amplitude',main=paste0(t,' (',sum(q[,t]),')'))
		segments(1:n,amp[1,],1:n,amp[3,])
		abline(h=ampl.thr,col='red',lty=2)
		legend('topright',lwd=1,legend='80% conf. interval')
		tissue1[[t]] = tissue2[[t]] = list()
		
		
		tissue1[[t]]$s2c1 = f1(getFreqOfSignOnCov(d[,t],q[,t],n))
		tissue1[[t]]$s2c2 = f1(getFreqOfSignOnCov(d[,t],q[,t] & a[,t] >= ampl.thr,n))

	}
	ylim = c(0,40)#c(0,max(unlist(tissue1)))
	for(t in colnames(d)){
		if(sum(q[,t])>0){
			plot(1:n,tissue1[[t]]$s2c1,col='red',t='l',xlab='coverage bin',ylab='% of significant',main=paste0(t,' (',nrow(d),')'),ylim=ylim)
			lines(1:n,tissue1[[t]]$s2c2,col='blue')
			legend('topleft',col=c('red','blue'),lwd=1,legend=paste0(c('only q-value','with ampl. thr'),' (',c(sum(q[,t]),sum(q[,t] & a[,t] >= ampl.thr)),')'))
		}else
			plot(1,t='n',xlab='coverage bin',ylab='% of significant',main=paste0(t,' (',nrow(d),')'))
	}
	plotT2T = function(d,a,q,n,ampl.thr,dir,same.dir){
		tissue2 = list()
		f2 = function(x){x[,5]/pmin(x[,3],x[,4])*100}
		for(t in colnames(d)){
			for(t2 in colnames(d)){
				if(t != t2){
					if(is.null(dir))
						dir. = NULL
					else
						dir. = dir[,c(t,t2)]
					tissue2[[t]][[t2]] = f2(getFreqOfSignOnCov21(d[,t],d[,t2],q[,t] & a[,t] >= ampl.thr,q[,t2] & a[,t2] >= ampl.thr,n,dir.,same.dir))
				}
			}
		}
		
		ylim = c(0,ifelse(same.dir,90,25))#c(0,max(unlist(tissue2),na.rm=T))
		for(t in colnames(d)){
			plot(1:n,t='n',xlab='coverage bin',ylab=ifelse(is.null(dir),'intersect/min',ifelse(same.dir,'(intersect same dir)/min','(intersect opposite dir)/min')),main=t,ylim=ylim)
			grid()
			if(sum(q[,t])>0){
				for(t2 in setdiff(colnames(d),t)){
					lines(1:n,tissue2[[t]][[t2]],col=params$tissue.col[t2])
				}
			}
		}
	}
	plotT2T(d,a,q,n,ampl.thr,sign(cor2agerank[[s]][rownames(d),]),TRUE)
	plotT2T(d,a,q,n,ampl.thr,sign(cor2agerank[[s]][rownames(d),]),FALSE)
	mtext(s,outer=T)
}

getFreqOfSignOnCov2 = function(c1,c2,s1,s2,n=20){
	s1 = names(c1) %in% s1 
	s2 = names(c1) %in% s2
	r1 = number2bin(c1,n)
	r2 = number2bin(c2,n)
	sgn = total = matrix(NA,ncol=n,nrow=n)
	cov1 = cov2 = matrix(NA,ncol=5,nrow=n,dimnames = list(NULL,c('min','max','mean.log10','sgn.count','total')))
	for(i in 1:n){
		cov1[i,] = c(min(c1[r1==i]),max(c1[r1==i]),mean(log10(c1[r1==i])),sum(s1[r1==i]),sum(r1==i))
		cov2[i,] = c(min(c2[r2==i]),max(c2[r2==i]),mean(log10(c2[r2==i])),sum(s2[r2==i]),sum(r2==i))
		for(j in 1:n){
			f = r1==i & r2==j
			total[i,j] = sum(f)
			sgn[i,j] = sum(s1[f] & s2[f])
		}
	}
	list(cov1=cov1,cov2=cov2,sgn=sgn,total=total)
}

getFreqOfSignOnCov21 = function(c1,c2,s1,s2,n=20,dir=NULL,same.dir=TRUE){
	if(is.character(s1))
		s1 = names(c1) %in% s1 
	if(is.character(s2))
		s2 = names(c1) %in% s2
	covr = number2bin(pmin(c1,c2),n)
	r = matrix(NA,ncol=6,nrow=n,dimnames = list(NULL,c('log10.mean1','log10.mean2','sgn.count1','sgn.count2','sgn.count12','total')))
	for(i in 1:n){
		f = covr==i
		if(is.null(dir))
			both = sum(s1[f] & s2[f])
		else
			both = sum(s1[f] & s2[f] & xor(dir[f,1] == dir[f,2],!same.dir))
		r[i,] = c(mean(log10(c1[f])),mean(log10(c2[f])),sum(s1[f]),sum(s2[f]),both,sum(f))
	}
	r
}

getFreqOfSignOnCov = function(cov,sgn,n=20){
	if(is.character(sgn))
		sgn = names(cov) %in% sgn
	covr = number2bin(cov,n)
	t(sapply(split(data.frame(cov=cov,sgn=sgn),covr),function(x){
		x = x[x$cov > 0,]
		c(min.cov=min(x$cov),max.cov=max(x$cov),mean.log10.cov=mean(log10(x$cov)),sgn.count=sum(x$sgn),total=nrow(x))
	}))
}



plotAllAgeRelExonProfiles = function(psi,m,cls,base,s2e,descr,ncol=6,nrow=6){
	s2a=unique(m[,c('stage','age.use')])
	s2a = setNames(rank(s2a$age.use),s2a$stage)
	m$age.use = s2a[m$stage]
	birth=13.5
	
	vl = m
	vl=vl$age.use[order(abs(species[m$species[1],'gestation']-vl$days))[1]] + 0.5
	
	for(t in names(cls)){
		print(t)
		pdf(paste0(base,t,'.pdf'),h=nrow*2*1.2,w=ncol*2)
		sids = rownames(cls[[t]])[order(cls[[t]][,1],-cls[[t]][,2])]
		c = 0
		for(i in 1:length(sids)){
			cat('\r',i)
			sid = sids[i]
			if(cls[[t]][sid,1] != c){
				c = cls[[t]][sid,1]
				par(mfrow=c(nrow,ncol),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(0,0,5,1),oma=c(2,3,1,1),cex.main=0.8)
			}
			d = paste0(descr[s2e[[sid]],'Description'],collapse=', ')
			ll = 40
			d=paste0(sapply(1:ceiling(nchar(d)/ll),function(i){substr(d,(i-1)*ll+1,i*ll)}),collapse='\n')
			plotTissueAgeProile(psi[sid,],m,age.axis=ifelse((i-1)%%(ncol*nrow)>=(ncol*(nrow-1)),'rank','none'),ylim=c(0,1),yaxt=ifelse(i%%ncol==1,'s','n'),
													main=paste0(sid,', c',c,'\n',
																			paste0(s2e[[sid]],collapse=', '),'\n',
																			d,'\n'))
			abline(v=vl,lty=2)
		}
		dev.off()
	}
}

getRowsHGMDCounts = function(hgm.in.tissue,dir,col,total,property,values=NULL){
	sapply(hgm.in.tissue,function(x){
		x = x[[dir]]
		tot = x[total,col]
		x = x[grepl(paste0(property,':'),rownames(x)),col]
		names(x) = gsub(paste0(property,':'),'',names(x))
		if(!is.null(values))
			x = x[values]
		c(total=tot,x)
	})
}

plotTissue2Concept = function(hgm.in.tissue){
	l = list()
	l$up.alt    = getRowsHGMDCounts(hgm.in.tissue,'up',1,'has.concept','concept',concepts)
	l$up.counst = getRowsHGMDCounts(hgm.in.tissue,'up',2,'has.concept','concept',concepts)
	l$dw.alt    = getRowsHGMDCounts(hgm.in.tissue,'dw',1,'has.concept','concept',concepts)
	l$dw.counst = getRowsHGMDCounts(hgm.in.tissue,'dw',2,'has.concept','concept',concepts)
	for(i in 1:length(l)){
		d = l[[i]]
		colnames(d) = paste0(colnames(d),'\n',d[1,])
		d = d[,ncol(d):1]
		c = d[-1,]
		f = sweep(d[-1,],2,d[1,],'/')
		imageWithText(apply(f,2,function(x){x/max(x)}),paste0(round(f*100),'%; ',c),names.as.labs=T,xlab='',ylab='',col=getPal(c('white','yellow','red'),100),main=names(l)[i])
	}
}

plotHGMDFreqForTissue = function(d,row,tot='with.mut',ylim=NULL,ylab='% exons with mutations',col=c('red','blue','gray'),...){
	tots = t(sapply(d,function(x)x[tot,]))[,2:1]
	vals = t(sapply(d,function(x)x[row,]))[,2:1]
	muts = t(sapply(d,function(x)x[row,]))[,4:3]
	f = vals/tots
	sd = sapply(1:length(tots),function(i){prop.test(vals[i],tots[i],correct=FALSE)$conf.int})
	if(is.null(ylim))
		ylim = c(0,max(sd)*120)
	b = barplot(f*100,beside = T,names.arg = c('const. same gene','alternative'),ylab=ylab,ylim=ylim,col=col,...)
	segments(b,sd[1,]*100,b,sd[2,]*100)
	text(b,sd[2,]*100,paste0('e:',vals,'\nm:',muts),adj = c(0.5,-0.1))
	text(b,0,tots,adj=c(0.5,1.1),xpd=T)
}

getNumberOfExonsWithHGM = function(sids,const,h2s,cod='c',locs=c(exon='e',ss='a|d|A|D',intron='u|p|w'),ages=c(congenital='c',postnatal='p',adult='a',all='l')){
	f = function(seg1,seg2,sids){
		c(sum(seg1 %in% sids),sum(seg2 %in% sids),sum(sids %in% seg1),sum(sids %in% seg2))
	}
	sids = sids[anns$human[sids,'cod'] %in% cod]
	const = rownames(const)[const$cod %in% cod & const$gene_id %in% anns$human[sids,'gene_id']]
	r = matrix(NA,ncol=4,nrow=0,dimnames = list(c(),c('exn.alt','exn.const','mut.alt','mut.const')))
	r = rbind(r,total=c(length(sids),length(const),NA,NA))
	
	r = rbind(r,with.mut=f(sids,const,h2s$seg.id))
	
	for(l in 1:length(locs)){
		r = rbind(r,f(sids,const,h2s$seg.id[grepl(locs[l],h2s$loc)]))
		rownames(r)[nrow(r)] = paste0('loc:',names(locs)[l])
	}
	
	r = rbind(r,has.age=f(sids,const,h2s$seg.id[h2s$hgmd.id %in% names(h2a)]))
	
	for(a in 1:length(ages)){
		hids = names(h2a)[sapply(h2a,function(x){ages[a] %in% x})]
		r = rbind(r,f(sids,const,h2s$seg.id[h2s$hgmd.id %in% hids]))
		rownames(r)[nrow(r)] = paste0('age:',names(ages)[a])
	}
	
	r = rbind(r,has.concept=f(sids,const,h2s$seg.id[h2s$hgmd.id %in% names(h2c)]))
	concepts = names(sort(table(unlist(h2c)),decreasing = T))
	for(c in 1:length(concepts)){
		hids = names(h2c)[sapply(h2c,function(x){concepts[c] %in% x})]
		r = rbind(r,f(sids,const,h2s$seg.id[h2s$hgmd.id %in% hids]))
		rownames(r)[nrow(r)] = paste0('concept:',concepts[c])
	}
	r
}

calcCorToAgeRank = function(psi,m){
	m = m[colnames(psi),]
	m = m[order(m$tissue,m$age.use),]
	psi = psi[,rownames(m)]
	tissues = unique(m$tissue)
	r = matrix(NA,ncol=length(tissues),nrow=nrow(psi),dimnames = list(rownames(psi),tissues))
	for(t in tissues){
		p = psi[,m$tissue==t]
		a = 1:ncol(p)
		r[,t] =  apply(p,1,function(x)cor(x,a,u='p'))
	}
	r
}

plotClusterPhastcons = function(p,cl,cnst,ann,s1,s2,ylim=c(0,1),cod=c('c','n','p')){
	cl = lapply(1:max(cl[,1]),function(i){sid=rownames(cl)[cl[,1]==i];sid[ann[sid,'cod'] %in% cod]})
	cnst = cnst[cnst$cod %in% cod,]
	x = 1:((s2-s1+1)*2)
	for(i in 1:length(cl)){
		cnstp = rownames(cnst)[cnst$gene_id %in% ann[cl[[i]],'gene_id']]
		cnstp = getPhastconsProf(phastcons[cnstp],s1,s2)
		p=getPhastconsProf(phastcons[cl[[i]]],s1,s2)
		plotArea(x,cnstp,col='black',new = TRUE,ylim=ylim,xaxt='n',ylab='mean phastcons')
		grid(1,10,col = '#999999')
		plotArea(x,p,col='red',new = FALSE,ylim=)
		axis(1,c(1,200-s1+1,2*(s2-s1+1)-200+s1,length(x)),c(-(200-s1+1),'acc','don',(200-s1+1)))
	}
}

getOrthSegClids = function(oids,cl){
	if(is.matrix(cl))
		cl = setNames(cl[,1],rownames(cl))
	r = oids
	r[,] = NA
	for(i in 1:ncol(r)){
		rownames(r) = oids[,i]
		cmn = intersect(oids[,i],names(cl))
		r[cmn,i] = cl[cmn]
	}
	rownames(r) = oids[,1]
	r = r[apply(!is.na(r),1,sum)>0,]
	r
}

calcOrthClusteringAgrement = function(oids,cl){
	library(lsr)
	ocl = getOrthSegClids(oids,cl)
	r = list()
	clusters = unique(cl)
	for(s1 in 1:(ncol(oids)-1)){
		t = c()
		for(s2 in (s1+1):ncol(oids)){
			t1 = factor(ocl[apply(is.na(ocl[,c(s1,s2)]),1,sum)==0,s1],levels=clusters)
			t2 = factor(ocl[apply(is.na(ocl[,c(s1,s2)]),1,sum)==0,s2],levels=clusters)
			t['cmn.sign']       = length(t1)
			t['agree']          = mean(t1==t2)
			t['agree.expected'] = sum(table(t1)/length(t1)*table(t2)/length(t2))
			t['cramerV']        = cramersV(t1,t2)
			t['chisq.pv']       = chisq.test(t1,t2)$p.value
			r[[paste0(colnames(oids)[s1],'-',colnames(oids)[s2])]] = t
		}
	}
	do.call(rbind,r)
}


getNAStatInTissueAgeRelated = function(t,s,amp.thr=0.05,MEL=500,MAS=5,sites='ad'){
	f=per.tissue.age.qv[[s]][,t]<=0.05 & per.tissue.age.ampl[[s]][,t]>=amp.thr & anns[[s]]$length<=MEL & anns[[s]]$alt.size<=MAS & anns[[s]]$sites %in% sites
	f[is.na(f)] = FALSE
	p = psi.tsm[[s]][f,meta.tsm[colnames(psi.tsm[[s]]),'tissue']==t]
	c(sample.no=ncol(p),table(apply(is.na(p),1,sum)))
}

getSignSegOverlap = function(s,only.sgn=TRUE,amp.thr=0.05,MEL=500,MAS=5,sites='ad',max.nn.frac=1,with.dir=FALSE){
	f=per.tissue.age.qv[[s]]<=0.05 & per.tissue.age.ampl[[s]]>=amp.thr
	f[is.na(f)] = FALSE
	f = f[anns[[s]]$length<=MEL & anns[[s]]$alt.size<=MAS & anns[[s]]$sites %in% sites,]
	psi = psi.tsm[[s]][rownames(f),]
	for(t in colnames(f)){
		p = psi[,meta.tsm[colnames(psi),'tissue']==t]
		f[,t] = f[,t] & apply(is.na(p),1,mean)<=max.nn.frac
	}
	if(only.sgn)
		f = f[apply(f,1,sum)>0,]
	if(with.dir){
		dir = cor2agerank[[s]][rownames(f),colnames(f)]
		ts = colnames(f)
		f = cbind(f & dir > 0,f & dir < 0)
		colnames(f) = c(paste0(ts,' up'),paste0(ts,' down'))
	}
	pv = ol = un = matrix(NA,ncol=ncol(f),nrow=ncol(f),dimnames = list(colnames(f),colnames(f)))
	for(i in 1:(ncol(f)-1)){
		for(j in (i+1):ncol(f)){
			ol[i,j]=ol[j,i] = sum(f[,i] & f[,j])
			un[i,j]=un[j,i] = sum(f[,i] | f[,j])
			pv[i,j]=pv[j,i] = fisher.test(factor(f[,i],levels=c(F,T)) , factor(f[,j],levels=c(F,T)),a='g')$p.value
		}
	}
	list(pv=pv,ol=ol,un=un,totals = apply(f,2,sum),total=nrow(f))
}

getSelfClusterCor = function(cl){
	r = rep(NA,length(cl$clustering))
	for(c in 1:max(cl$clustering)){
		f = cl$clustering == c
		r[f] = 1-cl$data[cl$id.med[c],f]
	}
	r
}

reassignToClusters = function(psi,cl){
	meds = psi[colnames(cl$data)[cl$id.med],]
	r = matrix(NA,nrow=nrow(psi),ncol=2)
	rownames(r) = rownames(psi)
	for(i in 1:nrow(psi)){
		cor = apply(meds,1,function(m)cor(m,psi[i,],use='pair'))
		r[i,1] = order(cor,decreasing = T)[1]
		r[i,2] = cor[r[i,1]]
		if(is.na(r[i,2]))
			r[i,1] = NA
	}
	r
}

clusterPSIs.by.cor.pam = function(p,n){
	pf = p[apply(is.na(p),1,sum)==0,]
	c = cor(t(pf))
	cl = pam(1-c,n)
	r=reassignToClusters(p,cl)
	r[,1] = reorderClustersBySize(r[,1])
}

clusterTissueAgeRelExons = function(s,t,n,amp.thr=0.05,MEL=500,MAS=5,max.nn.frac=0.5,max.nn.frac.core=0,sites='ad'){
	f=per.tissue.age.qv[[s]][,t]<=0.05 & per.tissue.age.ampl[[s]][,t]>=amp.thr & anns[[s]]$length<=MEL & anns[[s]]$alt.size<=MAS & anns[[s]]$sites %in% sites
	f[is.na(f)] = FALSE
	psi = psi.tsm[[s]][f,meta.tsm[colnames(psi.tsm[[s]]),'tissue']==t]
	psi = psi[apply(is.na(psi),1,mean)<=max.nn.frac,]
	psi.core = psi[apply(is.na(psi),1,mean)<=max.nn.frac.core,]
	if(length(psi.core) == 0)
		return(NULL)
	c = cor(t(psi.core),use = 'pair')
	r=lapply(n,function(i){
		cat(i)
		r = reassignToClusters(psi,pam(1-c,i))
		r[,1] = reorderClustersBySize(r[,1])
		r
	})
	names(r) = n
	r
}

getAlignedSgnAgePSIs = function(psis,m,t,amp.thr=0.05,MEL=500,MAS=5,sites='ad'){
	m = m[m$tissue==t,]
	r = NULL
	ages = sort(unique(m$mouse.days))
	ff = rep(TRUE,length(ages))
	for(s in names(psis)){
		print(s)
		sm = m[m$species==s,]
		f=per.tissue.age.qv[[s]][,t]<=0.05 & per.tissue.age.ampl[[s]][,t]>=amp.thr & anns[[s]]$length<=MEL & anns[[s]]$alt.size<=MAS & anns[[s]]$sites %in% sites
		f[is.na(f)] = FALSE
		p = psis[[s]]$ir[f,rownames(sm)]
		tmp = matrix(NA,ncol=length(ages),nrow=nrow(p))
		rownames(tmp) = rownames(p)
		for(i in 1:length(ages)){
			if(sum(sm$mouse.days==ages[i])>0){
				tmp[,i] = apply(p[,sm$mouse.days==ages[i],drop=FALSE],1,mean,na.rm=TRUE)
			}else
				ff[i] = FALSE
		}
		r = rbind(r,tmp)
	}
	colnames(r) = ages
	r[,ff]
}