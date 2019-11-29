plotScaleHist = function(d,main,from=0.7,to=2.5){
	total = dim(d)[1]
	d = d[d[,1]!=from & d[,1]!=to,1]
	b = (0:20)/20*(to-from)+from
	h = hist(d,breaks=b,xlab='Scale',main=paste(main,' (',length(d),'/',total,')',sep=''),ylab='Segment count')
	dns = density(d,from=from,to=to)
	dns$y = dns$y/max(dns$y)*max(h$counts)
	lines(dns,col='blue')
	max = dns$x[order(dns$y,decreasing=T)[1]]
	abline(v=max,col='red')
	legend('topright',bty='n',col='red',lwd=1,legend=paste('max=',round(max,4),sep=''))
}

#' Looks for the best age scaling coefficient (x1_ = x1*coef) that minimizes distance between curves
#' 
#' Distance calculated in shared (after scaling) age range after mean substraction
#' 
#' @param ir1 matrix of expressions in species 1
#' @param ir2 matrix of expressions (of the same genes) in species 2
#' @param a1 age in species 1
#' @param a2 age in species 2
#' @param from minimal possible scaling coefficient
#' @param to maximal possible scaling coefficient
#' @param n number of points used in approximation
#' @param df number of spline degrees of freedom used in approximation
findScales = function(ir1,ir2,a1,a2,from,to,n=100,df=4){
	ir1 = as.matrix(ir1)
	ir2 = as.matrix(ir2)
	r = matrix(nrow=dim(ir1)[1],ncol=2)
	colnames(r) = c('min.scale','min.dist')
	rownames(r) = rownames(ir1)
	for(i in 1:dim(ir1)[1]){
		r[i,] = optScale(a1,ir1[i,],a2,ir2[i,],from,to,n,df)
	}
	r
}

optScale = function(x1,y1,x2,y2,from,to,n,df){
	x1 = x1[!is.na(y1)]
	x2 = x2[!is.na(y2)]
	y1 = y1[!is.na(y1)]
	y2 = y2[!is.na(y2)]
	
	x1. = sort(x1)[c(2,length(x1)-1)]
	x2. = sort(x2)[c(2,length(x2)-1)]
	
	from = max(from,x2.[1]/x1.[2])
	to = min(to,x2.[2]/x1.[1])
	
	f1 = smooth.spline(x1,y1,df=df)
	f2 = smooth.spline(x2,y2,df=df)
	min = optimise(function(s){fitDist(f1,f2,range(x1),range(x2),s,n)},c(from,to),maximum=F)
	r = c(min$minimum,min$objective)
	names(r) = c('min.scale','min.dist')
	v = fitDist(f1,f2,range(x1),range(x2),from,n)
	s = from
	to.v = fitDist(f1,f2,range(x1),range(x2),to,n)
	if(to.v<v){
		v = to.v
		s = to
	}
	if(v<r[2]){
		r[1] = s 
		r[2] = v
	}
	r
}

getEqD = function(mn,mx,n){
	n = n-1
	(0:n)/n*(mx-mn)+mn
}

fitDist = function(f1,f2,r1,r2,scale,n){
	r1 = r1*scale
	r = c(max(r1[1],r2[1]),min(r1[2],r2[2]))
	if(r[1]>=r[2]) {
		cat('bad scale:',scale,'\n')
		print(r1)
		print(r2)
		return(NA)
	}
	a = getEqD(r[1],r[2],n) 
	p1 = predict(f1,a/scale)$y
	p2 = predict(f2,a)$y
	p1 = p1-mean(p1)
	p2 = p2-mean(p2)
	r = mean(abs((p1-p2)))
	r
}
