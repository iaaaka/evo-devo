# 1 ######
plotAsEventCount = function(d,pchs,by.tissue=FALSE,...){
	x = 1:length(d[[1]]) + rep((0:6)*3,each=7)
	plot(1,t='n',xaxt='n',xlab='',xlim=range(x),ylim=range(0,unlist(d),na.rm=T)*1.05,yaxs = 'i',...)
	for(ss in names(d)){
		d[[ss]][d[[ss]]==0] = NA
		if(by.tissue)
			points(x,t(d[[ss]]),pch=pchs[ss],col=rep(params$tissue.col[rownames(d[[ss]])],each=ncol(d[[ss]])))
		else
			points(x,d[[ss]],pch=pchs[ss],col=rep(params$tissue.col[rownames(d[[ss]])],times=ncol(d[[ss]])))
	}
	if(by.tissue){
		text(x,0,species[colnames(d[[1]]),'short'],adj = c(0.5,1),xpd=T,cex=0.6)
		xx = sapply(split(x,rep(1:7,each=7)),mean)
		text(xx,0,rownames(d[[1]]),adj = c(0.5,2.3),xpd=T,cex=0.8)
	}else{
		text(x,0,substr(rownames(d[[1]]),1,1),adj = c(0.5,1),xpd=T,cex=0.6)
		xx = sapply(split(x,rep(1:7,each=7)),mean)
		text(xx,0,colnames(d[[1]]),adj = c(0.5,2.3),xpd=T,cex=0.8)
	}
}

plotASTypes = function(pchs,cnst.col='gray',alt.col='orange'){
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
		text(50,top-dd+h*3,adj = c(0.5,0.5),t,cex=0.4,xpd=NA)
		return(dd+h*10)
	}
	
	# CE
	dd0=f('ad','CE',h*0,'Cassette Exon')
	rect((x1+x2-w/2)/2,top+h,(x1+x2+w/2)/2,top-h,col=alt.col,border=alt.col)
	# AA
	dd1=f('dd','AD',dd0,'Alternative Donor')
	rect(x1+w,top+h-dd0,x1+w*1.5,top-h-dd0,col=alt.col,border=alt.col)
	dd0 = dd1
	# DD
	dd1 = f('aa','AA',dd0,'Alternative Acceptor')
	rect(x2-w*1.5,top+h-dd0,x2-w,top-h-dd0,col=alt.col,border=alt.col)
	dd0 = dd1
	# DA
	dd1=f('da','RI',dd0,'Retained Intron')
	rect(x1+w,top+h-dd0,x2-w,top-h-dd0,col=alt.col,border=alt.col)
	top - dd0-h*3
}

# 2 #####
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

plotExampleDLG3 = function(fig=c(0,1,0,1),l){
	fig[3] = fig[3] + grconvertY(1,'lines','ndc')
	sid1 = 12502 ##Dlg3 - looks nice + Intellectual disability (HGMD)
	sid2 = 8396
	sid3 = 12506
	h = fig[4]-fig[3]
	w = fig[2]-fig[1]
	rows=6
	ad = dlg3.mdata$ann[dlg3.mdata$ann$sites=='ad',]
	alt.seg = rownames(ad) %in% rownames(orth.age.dpsi$mouse)[c(sid1,sid2,sid3)]
	alt.seg.coor = c(min(dlg3.mdata$ann$start[alt.seg]),max(dlg3.mdata$ann$stop[alt.seg]))
	par(fig=c(fig[1],fig[2],fig[4]-1/rows*h,fig[4]),new = TRUE)
	plotReadCov(dlg3.mdata$cov,reverse = (dlg3.mdata$ann$strand[1] == -1),min.junc.cov = 5,bty='n',xlab=paste0('Chr ',dlg3.mdata$ann$chr_id[1]),ylab='',junc.col='black',junc.lwd=1,axes=F)
	y=grconvertY(-0.1,'nfc','user')
	segments(min(dlg3.mdata$ann$start),y/2,max(dlg3.mdata$ann$stop),y/2,xpd=NA)
	rect(ad$start,y*0.8,ad$stop,y*0.2,xpd=NA,border=NA,col=ifelse(alt.seg,'red','black'))
	alt.seg.y = grconvertY(y,'user','ndc')
	plotPanelLetter(l,lab.cex)
	
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
	y = max(dlg3.mdata$cov$cov[dlg3.mdata$cov$x>=dlg3.mdata$zoom.coor[1] & dlg3.mdata$cov$x<=dlg3.mdata$zoom.coor[2]])
	segments(alt.seg.coor0,c(alt.seg.y,alt.seg.y),dlg3.mdata$zoom.coor,rep(y,2),xpd=NA,lty=2)
	
	for(j in 1:6){
		s = rownames(species)[-2][j]
		par(fig=c(fig[1]+w/6*(j-1),fig[1]+w/6*j,fig[4]-3/rows*h,fig[4]-2/rows*h),mar=c(0,0,0,0),new = TRUE)
		plot.new()
		plotPNG(paste0("figures/paper.figures/5/icons/",s,".png"),0.5,0.5,0.6)
		par(fig=c(fig[1]+w/6*(j-1),fig[1]+w/6*j,fig[4]-4/rows*h,fig[4]-3/rows*h),mar=c(1,1,0.5,0),new = TRUE)
		plotTissueAgeProile(orth.seg.ad.tsm[[s]][sid1,],meta.tsm,age.axis = 'rank',yaxt=ifelse(j==1,'s','n'),bty='n',ylim=c(0,1),xlab='',ylab=ifelse(j==1,'PSI',''),plot.xaxt=F,pch=19,cex=0.3,lwd=1)
		par(fig=c(fig[1]+w/6*(j-1),fig[1]+w/6*j,fig[4]-5/rows*h,fig[4]-4/rows*h),new = TRUE)
		plotTissueAgeProile(orth.seg.ad.tsm[[s]][sid2,],meta.tsm,age.axis = 'rank',yaxt=ifelse(j==1,'s','n'),bty='n',ylim=c(0,1),xlab='',ylab=ifelse(j==1,'PSI',''),plot.xaxt=F,pch=19,cex=0.3,lwd=1)
		if(j==1) title(ylab='PSI',xpd=NA)
		par(fig=c(fig[1]+w/6*(j-1),fig[1]+w/6*j,fig[4]-6/rows*h,fig[4]-5/rows*h),new = TRUE) 
		plotTissueAgeProile(orth.seg.ad.tsm[[s]][sid3,],meta.tsm,age.axis = 'rank',yaxt=ifelse(j==1,'s','n'),bty='n',ylim=c(0,1),xlab='',ylab=ifelse(j==1,'PSI',''),plot.xaxt=F,pch=19,cex=0.3,lwd=1)
	}
	x = grconvertX(mean(fig[1:2]),'ndc','user')
	y = grconvertY(fig[3],'ndc','user')
	text(x,y,'Development',xpd=NA,adj=c(0.5,0))
	y = grconvertY(fig[4],'ndc','user')
	text(x,y,"DLG3 - membrane-associated guanylate kinase",adj=c(0.5,1.4),xpd=NA)
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



plotLine1 = function(x,y,col,area.alpha=0.2,lwd=3,leg.pos=c(0.01,0.01),leg.adj=c(0,0),leg.cex=0.5,add=T,...){
	o = order(x)
	x = x[o]
	y = y[o]
	if(!add)
		plot(x,y,...)
	f = !is.na(y)
	p = predict(lm(y[f] ~ x[f]),interval='conf')
	col.area =  col2rgb(col)[,1]
	col.area = rgb(col.area[1],col.area[2],col.area[3],area.alpha*255,maxColorValue = 255)
	polygon(c(x[f],rev(x[f])),c(p[,2],rev(p[,3])),border=NA,col=col.area)
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


plot4C.PeakChange = function(sp,peak.changes,l,yrange = c(0,0.25),plot.tissue.lab=FALSE,plot.xlab=TRUE){
	left = grconvertX(1,'lines','ndc')
	
	plot.pos = seq(left,1,length.out = length(peak.changes)+1)
	f = function(x)log(x+1)
	atx = c(0,10,100,300,1000,2000,5000)
	aty = c(50,100,150,200,300,400)
	for(t in 1:length(peak.changes)){
		par(fig=c(plot.pos[t],plot.pos[t+1],yrange[1],yrange[2]),new=T)
		x = f(peak.changes[[t]]$ge)
		y = f(peak.changes[[t]]$as)
		plot(x,y,pch=19,col=params$tissue.col[t],axes=F,xlab='',ylab='',xpd=NA)
		at. = atx[atx <= max(peak.changes[[t]]$ge)]
		axis(1,f(at.),at.)
		at. = aty[aty <= max(peak.changes[[t]]$as) & aty >= min(peak.changes[[t]]$as)]
		axis(2,f(at.),at.)
		plotLine1(x,y,col=params$tissue.col[t],leg.pos=c(0.95,0.05),leg.adj = c(1,0))
		if(plot.tissue.lab)
			text(grconvertX(0.5,'npc','user'),grconvertY(1,'nfc','user'),names(peak.changes)[t],adj=c(0.5,1),xpd=NA)
		if(t==1){
			ff = function(x){x/max(x)}
			ppar=par(fig=c(grconvertX(c(-0.09,0.45),'npc','ndc'),grconvertY(c(0.75,1.05),'nfc','ndc')),new=TRUE,cex=0.4,mar=par('mar')*0.2,mgp=par('mgp')*0,xpd=NA)
			plot(ff(peak.changes$heart$ge),t='l',bty='n',xlab='Development',ylab='# events',axes=F)
			lines(ff(peak.changes$heart$as),col='red')
			text(grconvertX(c(1,1),'npc','user'),grconvertY(c(0.8,0.65),'nfc','user'),c("Exons","Genes"),col=c('red','black'),adj=c(0,1))
			par(cex=ppar$cex,mar=ppar$mar,mgp=ppar$mgp,xpd=T)
		}
		if(t==2)
			plotPNG(paste0("figures/paper.figures/5/icons/",sp,".png"),0.3,0.9,0.4)
	}
	text(x=grconvertX(0,'ndc','user'),y=grconvertY(yrange[2],'ndc','user'),labels=l,adj=c(0,1),font=2,cex=lab.cex,xpd=NA)
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
plotASSegStat = function(x,y,pv,col,xax,lty,pv.at=x[-1:-2],pv.thr=c(' '=1,' '=0.05,'*'=0.01,'**'=0.001,'***'=0),...){
	ylim = range(y)
	ylim[2] = ylim[2]+(ylim[2]-ylim[1])*0.05
	plot(x,y[,2],col=col,xlab='',ylim=ylim,xaxt='n',...)
	segments(x,y[,1],x,y[,3],col=cols,lty=lty)
	abline(h=y[1,2],lty=3,col=col[1])
	abline(h=y[2,2],lty=3,col=col[2])
	axis(1,xax,FALSE)
	text(xax,grconvertY(-0.05,'npc','user'),names(xax),adj=c(0,0),srt=-45,xpd=NA)
	
	
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
