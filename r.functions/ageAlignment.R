library(dtw)
if(pr_DB$entry_exists("correlation.pair"))
	pr_DB$delete_entry("correlation.pair")
pr_DB$set_entry(FUN = function(x,y){((1-cor(x,y,use='pair'))/2)^1}, names = "correlation.pair")

calcAllScales = function(all.tissues,all.orth,local,step.pattern,removeQambig,norm.psi,sgn){
	scales = list()
	if(all.orth){
		orth = loadAltOrthSegs(f = 'output/hqmrboc.orth.segs')
		orth = orth[anns$mouse[orth[,'mouse'],'sites'] == 'ad',]
	}
	for(s in c('human','macaque','rat','rabbit','opossum','chicken')){
		if(!all.orth){
			orth = loadAltOrthSegs(f = paste('output/',species[s,'short'],'m.orth.segs',sep=''))
			orth = orth[anns$mouse[orth[,'mouse'],'sites'] == 'ad',]
		}
		if(all.tissues)
			scales[[s]]= alignAges(orth,NULL,s,'mouse',psi.tsm,meta.tsm,oneBYone = F,norm.psi = norm.psi,step.pattern = step.pattern,open.end=local,open.begin=local,keep.internals=TRUE,removeQambig=removeQambig,sgn)
		else
			scales[[s]]= lapply(unique(meta.tsm$tissue),function(ts){alignAges(orth,ts,  s,'mouse',psi.tsm,meta.tsm,oneBYone = F,norm.psi = norm.psi,step.pattern = step.pattern,open.end=local,open.begin=local,keep.internals=TRUE,removeQambig=removeQambig)})
	}
	scales
}

alignAges = function(o,t=NULL,s1,s2,psi,m,oneBYone=FALSE,norm.psi=TRUE,sgn=NULL,step.pattern=symmetric1,dist.method="correlation.pair",removeQambig = FALSE,...){
	o = o[,c(s1,s2)]
	if(!is.null(sgn))
		o = o[sgn[[s1]][o[,1]] & sgn[[s2]][o[,2]],]
	if(!is.null(t)){
		m1 = m[m$species == s1 & m$tissue == t,]
		m2 = m[m$species == s2 & m$tissue == t,]
		m1 = m1[order(m1$days),]
		m2 = m2[order(m2$days),]
	}else{
		m1 = m[m$species == s1,]
		m2 = m[m$species == s2,]
	}
	
	psi1 = psi[[s1]][o[,1],rownames(m1)]
	psi2 = psi[[s2]][o[,2],rownames(m2)]
	
	if(norm.psi){
		psi1 = normRows(psi1)
		psi2 = normRows(psi2)
	}
	
	f = function(p,m,tissues){
		stages = unique(m[,c('days','stage')])
		stages = stages$stage[order(stages$days)]
		r = matrix(NA,ncol=length(stages),nrow=nrow(p)*length(tissues))
		colnames(r) = paste(m$species[1],stages)
		rownames(r) = paste(rep(rownames(p),times=length(tissues)),rep(tissues,each=nrow(p)))
		for(t in 1:length(tissues)){
			tmp = p[,rownames(m)[m$tissue==tissues[t]]]
			colnames(tmp) = m$stage[m$tissue==tissues[t]]
			tissue.stages =  m$stage[m$tissue==tissues[t]]
			tissue.stages = stages[stages %in% tissue.stages]
			r[((t-1)*nrow(p)+1):(t*nrow(p)),stages %in% tissue.stages] = tmp[,tissue.stages]
		}
		r
	}
	
	if(is.null(t)){ #consider tissues as sepatare segments
		tissues = unique(m$tissue)
		psi1 = f(psi1,m1,tissues)
		psi2 = f(psi2,m2,tissues)
	}
	
	if(oneBYone){
		f = which(apply(is.na(psi1),1,sum) == 0 & apply(is.na(psi2),1,sum) == 0)
		r = vector('list',length(f))
		names(r) = rownames(psi1)[f]
		for(i in f){
			r[[rownames(psi1)[i]]] = my.dtw(psi1[i,],psi2[i,],step.pattern = step.pattern,dist.method='Euclidean',removeQambig=removeQambig,...)$align
		}
	}else{
		r = my.dtw(t(psi1),t(psi2),step.pattern = step.pattern,dist.method = dist.method,removeQambig=removeQambig,...)
	}
	return(r)
}

my.dtw = function(v1,v2,step.pattern,dist.method,removeQambig,...){
	dist = dist(v1,v2,method=dist.method)
	d = dtw(dist,step.pattern = step.pattern,...)
	r = data.frame(q=d$index1,r=d$index2,qname=rownames(v1)[d$index1],rname=rownames(v2)[d$index2])
	r$dist = NA
	for(i in 1:nrow(r))
		r$dist[i] = dist[r[i,1],r[i,2]]
	if(removeQambig)
		r = filterAmbigQPaths(r)
	d$align=r
	d$distMatrix = dist
	d$v1 = v1
	d$v2 = v2
	d
}

filterAmbigQPaths = function(i){
	curr.q = i$q[1]
	min.v = i$dist[1]
	min.p = 1
	r = NULL
	for(j in 1:nrow(i)){
		if(curr.q == i$q[j]){
			if(min.v > i$dist[j]){
				min.v = i$dist[j]
				min.p = j
			}
		}
		if(curr.q != i$q[j]){
			r = rbind(r,i[min.p,])
			min.v = i$dist[j]
			min.p = j
			curr.q = i$q[j]
		}
		if(j == nrow(i))
			r = rbind(r,i[min.p,])
	}
	r
}

NW = function(m,FUN=function(a,b,l){a+b}){
	stop("deprecated. use dtw package")
	
	pr = pc = l = m
	pr[] = pc[] = NA
	l[] = 0
	pr[1,1] = pc[1,1] = 0
	for(r in 1:nrow(m)){
		for(c in 1:ncol(m)){
			if(r == 1){
				if(c == 1) next
				l[r,c] = l[r,c-1] + 1
				m[r,c] = FUN(m[r,c] , m[r,c-1],l[r,c])
				pr[r,c] = r
				pc[r,c] = c-1
			}else{
				if(c == 1){
					l[r,c] = l[r-1,c] + 1
					m[r,c] = FUN(m[r,c] , m[r-1,c],l[r,c])
					pr[r,c] = r-1
					pc[r,c] = c
				}else{
					v = order(c(m[r,c-1],m[r-1,c-1],m[r-1,c]),c(l[r,c-1],l[r-1,c-1],l[r-1,c]),decreasing = FALSE)[1]
					if(v == 1){
						l[r,c] = l[r,c-1] + 1
						m[r,c] = FUN(m[r,c] , m[r,c-1],l[r,c])
						pr[r,c] = r
						pc[r,c] = c-1
					}else if(v == 2){
						l[r,c] = l[r-1,c-1] + 1
						m[r,c] = FUN(m[r,c] , m[r-1,c-1],l[r,c])
						pr[r,c] = r-1
						pc[r,c] = c-1
					}else{
						l[r,c] = l[r-1,c] + 1
						m[r,c] = FUN(m[r,c] , m[r-1,c],l[r,c])
						pr[r,c] = r-1
						pc[r,c] = c
					}
				}
			}
		}
	}
	list(m=m,pr=pr,pc=pc,l=l)
}

NWp = function(nw){
	stop("deprecated. use dtw package")
	i = nrow(nw$m)
	j = ncol(nw$m)
	r = i
	c = j
	repeat{
		ni = nw$pr[i,j]
		nj = nw$pc[i,j]
		r = c(r,ni)
		c = c(c,nj)
		i = ni
		j = nj
		if(i == 1 & j== 1) break
	}
	data.frame(r=rev(r),c=rev(c))
}
