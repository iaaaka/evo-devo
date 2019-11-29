getsPSIbyEnsID = function(psi,stages,tissue,s2e,sp,use.random=FALSE){
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

makeASinGEpatterns = function(sp,use.random.seg=FALSE){
	as.in.ge = list()
	mc = read.csv(paste0('processed/GE.from.marg/',firstToupper(sp),'Clusters.csv'),row.names = 1)
	colnames(mc) = tolower(colnames(mc))
	for(tis in unique(meta$tissue)){
		cat(tis)
		psi=psi.tsm[sp]
		psi[[sp]] = psi[[sp]][anns[[sp]]$sites=='ad' & anns[[sp]]$cod!='n',]
		t = getsPSIbyEnsID(psi,border.stages,tis,seg2ens,sp,use.random=use.random.seg)
		t = t(t)
		t = t[intersect(rownames(t),rownames(mc)[!is.na(mc[,paste0(tis,'pattern')])]),]
		as.in.ge[[tis]] = cbind(data.frame(t),ge.pattern=mc[rownames(t),paste0(tis,'pattern')])
	}
	as.in.ge
}


plotBinomWithConf = function(tab,col='black',add.totals.to.lab=FALSE,...){
	#first column success, second = failure 
	f = t(apply(tab,1,function(x){r=binom.test(x[1],x[1]+x[2]);c(r$estimate,r$conf.int)}))
	x=1:nrow(f)
	plot(x,f[,1],pch=19,ylim=range(f),xaxt='n',xlab='',col=col,cex=2,...)
	segments(x,f[,2],x,f[,3],col=col,lwd=2)
	labs=rownames(f)
	if(add.totals.to.lab)
		labs = paste0(labs,'\n(',tab[,1]+tab[,2],')')
	axis(1,x,labs,las=3)
}

getdPSI2lfcCorOnAge = function(psi,exp,t1,t2,m,sp,cor.method,age.dpsi.thr=0){
	psi = psi[[sp]]
	exp = log2(exp[[sp]]+0.1)
	stages = intersect(m[m$species==sp & m$tissue == t1,'stage'],m[m$species==sp & m$tissue == t2,'stage'])
	r = unique(m[m$species==sp & m$stage %in% stages,c('stage','days','age.rank')])
	r = r[order(r$days),]
	r$rho = r$p.value = r$conf1 = r$conf2 = NA
	f = abs(psi[,paste(sp,t1,r$stage[1])]-psi[,paste(sp,t1,r$stage[nrow(r)])])>=age.dpsi.thr | 
		abs(psi[,paste(sp,t2,r$stage[1])]-psi[,paste(sp,t2,r$stage[nrow(r)])])>=age.dpsi.thr
	f[is.na(f)] = FALSE
	for(i in 1:nrow(r)){
		dpsi = abs(psi[f,paste(sp,t1,r$stage[i])]-psi[f,paste(sp,t2,r$stage[i])])
		lfc  = abs(exp[f,paste(sp,t1,r$stage[i])]-exp[f,paste(sp,t2,r$stage[i])])
		cor = cor.test(dpsi,lfc,method=cor.method,use='pair')
		r$rho[i] = cor$estimate
		r$p.value[i] = cor$p.value
		if(is.null(cor$conf.int))
			cor$conf.int = c(NA,NA)
		r$conf1[i] = cor$conf.int[1]
		r$conf2[i] = cor$conf.int[2]
	}
	r
}

getAS.GE.change = function(psi,rpkm,tissue,m,by.gene=FALSE,stages=NULL,rpkm.pseudo=0.1){
	m = m[colnames(psi),]
	m = m[m$tissue==tissue,]
	if(is.null(stages))
		stages = rownames(m)[order(m$days)[c(1,nrow(m))]]
	else
		stages = c(rownames(m)[m$stage==stages[1]],rownames(m)[m$stage==stages[2]])
	r = data.frame(psi1=psi[,stages[1]],psi2=psi[,stages[2]],exp1=log2(rpkm[,stages[1]]+rpkm.pseudo),exp2=log2(rpkm[,stages[2]]+rpkm.pseudo),gene.id = rownames(rpkm))
	r$dpsi = r$psi2-r$psi1
	r$log2fc=r$exp2-r$exp1
	if(by.gene)
		r = do.call(rbind,lapply(split(r,r$gene.id),function(x){x[order(-abs(x$dpsi))[1],]}))
	r
}

plotProportionOFSegOnExp = function(dpsi,exp,bins=20,dpsi.bins=c(0.2,0.3,0.5,1),col=c('yellow','orange','red'),cor.by.dPSI=TRUE,...){
	f = !is.na(dpsi) & !is.na(exp)
	dpsi = dpsi[f]
	exp = exp[f]
	exp.inx = number2bin(exp,bins)
	exp.names = round(sapply(split(exp,exp.inx),mean),1)
	frequ = freqd = sdu = sdd = list()
	min = Inf
	max = -Inf
	for(i in 1:(length(dpsi.bins)-1)){
		frequ[[i]] = sapply(split(dpsi > dpsi.bins[i] & dpsi <= dpsi.bins[i+1],exp.inx),mean,na.rm=T)
		sdu[[i]] = sapply(split(dpsi > dpsi.bins[i] & dpsi <= dpsi.bins[i+1],exp.inx),function(x){x=x[!is.na(x)];sd(x)/sqrt(length(x))})
		max = max(max,frequ[[i]]+sdu[[i]])
		
		freqd[[i]] = sapply(split(dpsi < -dpsi.bins[i] & dpsi >= -dpsi.bins[i+1],exp.inx),mean,na.rm=T)
		sdd[[i]] = sapply(split(dpsi < -dpsi.bins[i] & dpsi >= -dpsi.bins[i+1],exp.inx),function(x){x=x[!is.na(x)];sd(x)/sqrt(length(x))})
		min = min(min,-freqd[[i]]-sdd[[i]])
	}
	x = 1:bins
	plot(1,t='n',xlim=c(0,bins),ylim=c(min,max*1.1),ylab='Proportion of changed exons',xaxt='n',yaxt='n',...)
	axis(1,x,exp.names)
	axis(2,-20:20/20,abs(-20:20/20))
	for(i in 1:(length(dpsi.bins)-1)){
		lines(x,frequ[[i]],col=col[i],lwd=3)
		segments(x,frequ[[i]]-sdu[[i]],x,frequ[[i]]+sdu[[i]],col=col[i])
		
		lines(x,-freqd[[i]],col=col[i],lwd=3)
		segments(x,-freqd[[i]]-sdd[[i]],x,-freqd[[i]]+sdd[[i]],col=col[i])
	}
	abline(h=0,lty=2)
	if(cor.by.dPSI){
		rho = cor.test(abs(dpsi),exp,m='p')
		cor = paste0('Pearson(|dPSI|,x)=',round(rho$estimate,4),'; pv=',format(rho$p.value,digits=1,scientific=TRUE))
	}else{
		rho = cor.test(x,sapply(split(abs(dpsi)>0.5,exp.inx),mean),m='sp')
		cor = paste0('Spearman(freq(dPSI>0.5),x)=',round(rho$estimate,4),'; pv=',format(rho$p.value,digits=1,scientific=TRUE))
	}
	
	legend('topleft',lwd=3,col=c(col,NA),legend=c(paste0('dPSI in [',dpsi.bins[-length(dpsi.bins)],',',dpsi.bins[-1],']'),cor),bty = 'n')
	text(0, min/2,labels = 'PSI decrease',srt=90,adj = c(0.5,0))
	text(0, max/2,labels = 'PSI increase',srt=90,adj = c(0.5,0))
	invisible(c(rho$estimate,rho$p.value))
}

plotAllProportionOFSegOnExp = function(psi,rpkm,by.gene,bins,stages,cor.by.dPSI){
	r = array(NA,dim = c(7,3,2),dimnames=list(unique(meta$tissue),1:3,1:2))
	for(t in  unique(meta$tissue)){
		x = getAS.GE.change(psi,rpkm,t,meta.tsm,by.gene = by.gene,stages = stages[t,])
		main = paste0(t,' (',stages[t,1],'-',stages[t,2],')')
		r[t,1,]=plotProportionOFSegOnExp(x$dpsi,x$exp1,xlab='pre-natal log2(RPKM)',main=main,bins=bins,cor.by.dPSI = cor.by.dPSI)
		r[t,2,]=plotProportionOFSegOnExp(x$dpsi,x$exp2,xlab='adult log2(RPKM)',main=main,bins=bins,cor.by.dPSI = cor.by.dPSI)
		r[t,3,]=plotProportionOFSegOnExp(x$dpsi,x$log2fc,xlab='log2FC',main=main,bins=bins,cor.by.dPSI = cor.by.dPSI)
	}
	plotCorStat(r)
	invisible(r)
}

plotCorStat = function(r,species=NULL,pv.thr=0.001,...){
	col = unique(meta[,c('tissue','col')])
	rownames(col) = col$tissue
	col = col[dimnames(r)[[1]],'col']
	n = c('pre-natal log2(RPKM)','adult log2(RPKM)','log2FC')
	for(i in 1:3){
		if(is.null(species))
			main = n[i]
		else
			main = paste0(species,' (',n[i],')')
		barplot(r[,i,1],col=paste0(col,ifelse(r[,i,2]<pv.thr,'FF','15')),ylab='Pearson cor',xlab='tissue',main=main,...)
	}
}