getASComplexity = function(psi,sgn=NULL,type=c('alt','dpsi'),alt.thr=0.1,plot=TRUE,meta=meta.tsm,main='',...){
	type = type[1]
	altF = function(p,t,thr){
		apply(p,2,function(x){
			if(t=='alt')
				mean(x>=thr & x<=1-thr,na.rm=T)
			else
				mean(0.5-abs(x-0.5),na.rm=T)
		})
	}
	if(is.null(sgn)){
		main = paste0(main,' (',nrow(psi),')')
		r=altF(psi,type,alt.thr)
	}else if(!is.matrix(sgn)){
		main = paste0(main,' (',nrow(psi[sgn,]),')')
		r=altF(psi[sgn,],type,alt.thr)
	}else{
		r=unlist(lapply(colnames(sgn),function(t){altF(psi[sgn[,t],meta[colnames(psi),'tissue']==t],type,alt.thr)}))
	}
	if(plot)
		plotTissueAgeProile(r,meta.tsm,main=main,...)
	invisible(r)
}

plotSpeciesADAltProp = function(s,sgn=FALSE,main,col,psi.thr=0.1){
	m = getAgeASchanges(psi.tsm,meta.tsm,0.5,border.stages,s)
	mm = apply(m,2,function(x){x %in% c('u','d') & anns[[s]]$sites=='ad'})
	f = TRUE
	if(sgn)
		f = apply(mm,1,sum)>0
	
	getASComplexity(psi.tsm[[s]][anns[[s]]$sites=='ad' & f,],type='alt',age.axis='rank',main=main,ylab='proportion PSI in [0.1,0.9]',alt.thr = psi.thr)
	plotTissueAgeProile(apply(!is.na(psi.tsm[[s]][anns[[s]]$sites=='ad' & f,]),2,sum),meta.tsm,ylab='# of expressed exons',main='Expressed exons')
	tissues = c('brain','heart','liver','ovary','testis')
	if(s == 'human')
		tissues = c('brain','heart','liver','testis')
	z=getAltExonStat(psi.tsm[[s]][anns[[s]]$sites=='ad' & f,],meta.tsm,psi.thr,tissues = tissues,na.as.cnst=F)
	z = z[,apply(z==0,2,sum)==0]
	mar=par(mar=c(4,2,1.5,5))
	areaplot(z,col=col[rownames(z)],ylab='# of exons',main=main,xaxt='n',xlab='Age')
	axis(1,1:ncol(z),colnames(z))
	z = sweep(z,2,apply(z,2,sum),'/')
	areaplot(z,col=col[rownames(z)],ylab='proportion of exons',main=main,xaxt='n',xlab='Age')
	axis(1,1:ncol(z),colnames(z))
	par(mar=mar)
}
