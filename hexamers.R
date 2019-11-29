library(seqinr)
library(MyLogo)
library(SAJR)
library(plyr)
library(doMC)
source('code/r.functions/ssPWM.F.R')
source('code/r.functions/load.all.data.F.R')
don.pwm = readRDS('Rdata/old/don.pwm.Rdata')
acc.pwm = readRDS('Rdata/old/acc.pwm.Rdata')
options(stringsAsFactors = FALSE)
#all.anns = readRDS('Rdata/all.anns.Rdata')
anns = readRDS('Rdata/anns.Rdata')
species = readRDS('Rdata/species.Rdata')
psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
border.stages = readRDS('Rdata/border.stages.Rdata')
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')
orth.per.tissue.age.qv = readRDS('Rdata/orth.per.tissue.age.qv.Rdata')
#here I didn't use stat tests 
# hex.dws.age03 = readRDS('Rdata/hex.dws.age03.Rdata') 
# hex.ups.age03 = readRDS('Rdata/hex.ups.age03.Rdata')

params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)


# 
# mar = 200
# for(s in rownames(species)){
# 	ad = all.anns[[s]][all.anns[[s]]$sites=='ad',]
# 	write.table(cbind(ad$chr_id,'.','.',as.integer(ad$start-mar),as.integer(ad$stop+mar),'.',ifelse(ad$strand==1,'+','-'),'.',rownames(ad)),sep='\t',quote = F,row.names = FALSE,col.names = FALSE,file=paste('processed/extract.seq/new.ad.200/',s,'200.gtf',sep=''))
# }
# run qsub -t 1-7 extract.ad.seq.sh

#add 2nt
# for(s in rownames(species)){
# 	print(s)
# 	f = read.table(paste0('processed/extract.seq/new.ad.200/',s,'200.fa'),sep='\t')
# 	all.anns[[s]]$dint = all.anns[[s]]$upweight = all.anns[[s]]$dwweight = NA
# 	all.anns[[s]][f[,1],'dint'] = paste(substr(f[,2],199,200),substr(f[,2],nchar(f[,2])-199,nchar(f[,2])-198))
# 	
# 	all.anns[[s]][f[,1],'upweight'] = sapply(f[,2],function(s){getSeqWeight(acc.pwm,s,177)})
# 	all.anns[[s]][f[,1],'dwweight'] = sapply(f[,2],function(s){getSeqWeight(don.pwm,s,nchar(s)-202)})
# 	
# 	anns[[s]] = all.anns[[s]][rownames(anns[[s]]),]
# }
# 
# saveRDS(all.anns,'Rdata/all.anns.Rdata')
# saveRDS(anns,'Rdata/anns.Rdata')

#use only without 'n' and AG-GT
# fa = lapply(rownames(species),function(s)read.table(paste0('processed/extract.seq/new.ad.200/',s,'200.fa'),sep='\t'))
# names(fa) = rownames(species)
# fa = lapply(fa,function(x){setNames(tolower(x[,2]),x[,1])})

# orth.fa = list()
# for(s in names(fa)){
# 	orth.fa[[s]] = fa[[s]][rownames(orth.seg.ad[[s]]$seg)]
# 	orth.seg.ad[[s]]$seg$dint = all.anns[[s]][rownames(orth.seg.ad[[s]]$seg),'dint']
# }

# for(s in rownames(species)) fa[[s]] = fa[[s]][rownames(anns[[s]])[!is.na(anns[[s]]$dint) & anns[[s]]$dint=='AG GT']]
# sapply(fa,function(x)table(grepl('n',x)))
#saveRDS(fa,'Rdata/ad.alt.fa.Rdata')
# gc()
#saveRDS(orth.seg.ad,'Rdata/orth.seg.ad.Rdata')
#saveRDS(orth.fa,'Rdata/ad.orth.fa.Rdata')
fa = readRDS('Rdata/ad.alt.fa.Rdata')

age.segs = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,0.2,border.stages,s))
names(age.segs) = rownames(species)
for(s in names(age.segs)) age.segs[[s]][age.segs[[s]] != '-' & (is.na(per.tissue.age.qv[[s]]) | per.tissue.age.qv[[s]]>0.05)[rownames(age.segs[[s]]),colnames(age.segs[[s]])]] = 'n'

up.hex = lapply(fa,function(x){print(':');r=substr(x,1,200);countNmers(r[!grepl('n',r)],6)>0})
up = ftHexAge(up.hex,age.segs,'u')
dw = ftHexAge(up.hex,age.segs,'d')

up$ih.pv = apply(up$pv[,,-2],1:2,function(x)dirwin.hall(sum(x),length(x)))
up$ih.qv = apply(up$ih.pv,2,p.adjust,m='BH')

dw$ih.pv = apply(dw$pv[,,-2],1:2,function(x)dirwin.hall(sum(x),length(x)))
dw$ih.qv = apply(dw$ih.pv,2,p.adjust,m='BH')

hex.ups.age02sgn = list(up=up,dw=dw)
saveRDS(hex.ups.age02sgn,'Rdata/hex.ups.age02sgn.Rdata')


dw.hex = lapply(fa,function(x){print(':');l=nchar(x);r=substr(x,l-199,l);countNmers(r[!grepl('n',r)],6)>0})

up = ftHexAge(dw.hex,age.segs,'u')
dw = ftHexAge(dw.hex,age.segs,'d')

up$ih.pv = apply(up$pv[,,-2],1:2,function(x)dirwin.hall(sum(x),length(x)))
up$ih.qv = apply(up$ih.pv,2,p.adjust,m='BH')

dw$ih.pv = apply(dw$pv[,,-2],1:2,function(x)dirwin.hall(sum(x),length(x)))
dw$ih.qv = apply(dw$ih.pv,2,p.adjust,m='BH')

hex.dws.age02sgn = list(up=up,dw=dw)
saveRDS(hex.dws.age02sgn,'Rdata/hex.dws.age02sgn.Rdata')




apply(up$ih.qv<0.05,2,sum)
up6mers$cons.pv[1:2,]

# comp to hexamers in clusters
print(load('Rdata/old/ad.6mer.cl.enrich.Rdata'))
z=cor(cbind(hex.ups.age03$up$ih.pv,hex.ups.age03$dw$ih.pv),up6mers$cons.pv,m='sp')
z=cor(cbind(hex.dws.age03$up$ih.pv,hex.dws.age03$dw$ih.pv),dw6mers$cons.pv,m='sp')
rownames(z) = paste(rownames(z),rep(c('u','d'),each=7))
apply(z,2,function(x){r=setNames(x,rownames(z));r=sort(r,decreasing = T);r=r[r>0.2];if(length(r)>3)r=r[1:3];r})




z = cor(cbind(hex.ups.age03$up$ih.pv,hex.ups.age03$dw$ih.pv),cbind(hex.dws.age03$up$ih.pv,hex.dws.age03$dw$ih.pv),m='sp')
rownames(z) = colnames(z) = paste(rownames(z),rep(c('u','d'),each=7))
image(z)


cor(cbind(hex.ups.age03$up$ih.pv[,'brain'],hex.dws.age03$up$ih.pv[,'brain'],hex.ups.age03$dw$ih.pv[,'brain'],hex.dws.age03$dw$ih.pv[,'brain']),m='sp')

t = 'heart'
z = calcAllPairsFT(cbind(hex.ups.age03$up$ih.qv[,t]<0.05,
												 hex.dws.age03$up$ih.qv[,t]<0.05,
												 hex.ups.age03$dw$ih.qv[,t]<0.05,
												 hex.dws.age03$dw$ih.qv[,t]<0.05))
plotFTMotifSimMatrix(z,T)

hex.stat = rbind(apply(hex.ups.age03$up$ih.qv<0.05,2,sum),
								 apply(hex.dws.age03$up$ih.qv<0.05,2,sum),
								 apply(hex.ups.age03$dw$ih.qv<0.05,2,sum),
								 apply(hex.dws.age03$dw$ih.qv<0.05,2,sum))
hex2mot = read.table('output/hex2mot2sf.tab.gz')
table(hex.ups.age03$up$ih.qv[,'brain']<0.05,known=hex2mot$V2!='')

all.hex.qv = cbind(hex.ups.age03$up$ih.qv, hex.dws.age03$up$ih.qv, hex.ups.age03$dw$ih.qv, hex.dws.age03$dw$ih.qv)
colnames(all.hex.qv) = paste0(substr(colnames(all.hex.qv),1,1),rep(c('iu','id','eu','ed'),each=7))
t=apply(all.hex.qv<0.05,1,sum)
hex.ups.age03$up$ih.qv[t==8,]
hex2mot[t==8,]
t=table(t,known=hex2mot$V2!='')
plot(t[,2]/apply(t,1,sum))


i = rep(c(1,8,15,22),times=7)+rep(0:6,each=4)
f = apply(all.hex.qv<0.05,1,sum)

z=calcAllPairsFT(all.hex.qv[,i]<0.05)
plotFTMotifSimMatrix(z,F)
abline(v=0:7*4)
abline(h=0:7*4)

rownames(all.hex.qv)[apply(all.hex.qv[,c('biu','bed','hid','heu')]>0.05,1,sum)==0]
hex2mot['actaac',]
rownames(all.hex.qv)[apply(all.hex.qv[,c('bid','beu','hid','heu')]>0.05,1,sum)==0]
hex2mot['tttgct',]


par(mfrow=c(2,2),tck=-0.01,mgp=c(1.1,0.2,0),mar=c(6,3,1.5,0),oma=c(0,0,2,1))
barplotWithText(hex.stat[c(1,3),],col=rep(params$tissue.col,each=2),las=3,main='upstream',beside = T,den=c(-1,30),srt = 90,adj=c(-0.1,0.5),ylim=c(0,max(hex.stat[c(1,3),])*1.2))
legend('topright',fill='black',den=c(-1,30),legend=c('inclusion','exclusion'),bty = 'n')
barplotWithText(hex.stat[c(2,4),],col=rep(params$tissue.col,each=2),las=3,main='downstream',beside = T,den=c(-1,30),srt = 90,adj=c(-0.1,0.5),ylim=c(0,max(hex.stat[c(2,4),])*1.2))




hex.dws.age03$up$ih.qv['actaac',]

### chicken heart RBFOX ####
orth.fa = readRDS('Rdata/ad.orth.fa.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
orth.seg.ad.tsm = readRDS('Rdata/orth.seg.ad.tsm.Rdata')
orth.seg.ad = readRDS('Rdata/orth.seg.ad.Rdata')
border.stages = readRDS('Rdata/border.stages.Rdata')

age.dpsi = lapply(rownames(species),function(s)getAgeASchanges(orth.seg.ad.tsm,meta.tsm,0.1,border.stages,s,get.dPSI=T))
names(age.dpsi) = rownames(species)
age.dpsi$macaque = cbind(age.dpsi$macaque,ovary=1)

cisbp = read.table('input/cisbp-db/human/RBP_Information.txt',sep='\t',header = T)
dim(cisbp)
length(unique(cisbp$DBID))
table(cisbp$Motif_ID!='.')
length(unique(cisbp$DBID[cisbp$Motif_ID!='.']))
write.table(unique(cisbp$DBID[cisbp$Motif_ID!='.']),quote = F,row.names = F,col.names = F,'~/sf.with.motifs.txt')

hex.dws.age02sgn = readRDS('Rdata/hex.dws.age02sgn.Rdata')
hex.ups.age02sgn = readRDS('Rdata/hex.ups.age02sgn.Rdata')

hex.dws.age02sgn$up$or['tgcatg',,]
hex.dws.age02sgn$up$pv['tgcatg',,]

n.dw= sapply(names(orth.seg.ad),function(s){
	t=orth.fa[[s]][rownames(orth.seg.ad[[s]]$seg)]
	grepl('n',substr(t,nchar(t)-199,nchar(t)),fixed = TRUE)
})

orth.tgcatg.dw = sapply(names(orth.seg.ad),function(s){
	t=orth.fa[[s]][rownames(orth.seg.ad[[s]]$seg)]
	sapply(gregexpr('gcatg',substr(t,nchar(t)-199,nchar(t)),fixed = TRUE),function(x){sum(x!=-1)})
	})

orth.tgcatg.up = sapply(names(orth.seg.ad),function(s){
	t=orth.fa[[s]][rownames(orth.seg.ad[[s]]$seg)]
	sapply(gregexpr('gcatg',substr(t,1,199),fixed = TRUE),function(x){sum(x!=-1)})
})
# paste(sample(strsplit('tgcatg','')[[1]]),collapse='')
# orth.aggctt.dw = sapply(names(orth.seg.ad),function(s){
# 	t=fa[[s]][rownames(orth.seg.ad[[s]]$seg)]
# 	grepl('aggctt',substr(t,nchar(t)-199,nchar(t)))
# })


apply(orth.tgcatg.dw>0,2,sum)
apply(orth.tgcatg.up>0,2,sum)
#apply(orth.aggctt.dw,2,sum)
table(apply(orth.tgcatg.dw>0,1,sum))
table(apply(orth.tgcatg.up>0,1,sum))
#table(apply(orth.aggctt.dw,1,sum))

ss=setNames(rownames(species),rownames(species))
t=sapply(unique(meta$tissue),function(t)
	sapply(ss,function(s){
		setNames(fisher.test(table(
			factor(orth.tgcatg.dw[,s]>0,levels = c(F,T)),
			factor(orth.per.tissue.age.qv[[s]][,t]<0.05 & age.dpsi[[s]][,t]> 0.5,levels = c(F,T))))$estimate,NULL)
	}))
imageWithText(t,names.as.labs=T,col=c('blue','gray','#FF000010','#FF000020','#FF000040','#FF000080','#FF000080','#FF0000'),breaks=c(0,1,1.3,1.6,2,3,4,5,1e20))

getAdultSpec = function(d,t1,ts2,stages,s){
	stages = stages[[s]]
	rx = d[[s]][,paste(s,t1,stages[t1,2])] - apply(d[[s]][,paste(s,ts2,stages[ts2,2])],1,max,na.rm=T)
	rn = d[[s]][,paste(s,t1,stages[t1,2])] - apply(d[[s]][,paste(s,ts2,stages[ts2,2])],1,min,na.rm=T)
	rx[is.infinite(rx)] = NA
	rx[rx <= 0] = 0
	rn[is.infinite(rn)] = NA
	rn[rn >= 0] = 0
	if(sum(is.na(rx) != is.na(rn))>0 |
		 sum(rn !=0 & rx !=0,na.rm=T)>0) stop('something is wrong')
	rx[!is.na(rn) & rn !=0] = rn[!is.na(rn) & rn !=0]
	rx
}

cnt.tis = c('liver','kidney','ovary','testis')

fn = !n.dw[,'mouse'] & !n.dw[,'chicken']

sp = 'mouse'

mh = getAdultSpec(orth.seg.ad.tsm,'heart',cnt.tis,border.stages,sp) 
ch = getAdultSpec(orth.seg.ad.tsm,'heart',cnt.tis,border.stages,'chicken')
rh = getAdultSpec(orth.seg.ad.tsm,'heart',cnt.tis,border.stages,'rat')

mb = getAdultSpec(orth.seg.ad.tsm,'brain',cnt.tis,border.stages,sp) 
cb = getAdultSpec(orth.seg.ad.tsm,'brain',cnt.tis,border.stages,'chicken')


f=orth.per.tissue.age.qv$mouse[,'heart']<0.05 
plot(mh[f],age.dpsi$mouse[f,'heart'],xlim=c(-1,1),ylim=c(-1,1));abline(v=0,lty=3,col='red');abline(h=0,lty=3,col='red')
f=orth.per.tissue.age.qv$chicken[,'heart']<0.05 
f=T
plot(ch[f],age.dpsi$chicken[f,'heart'],xlim=c(-1,1),ylim=c(-1,1));abline(v=0,lty=3,col='red');abline(h=0,lty=3,col='red')
mm = orth.tgcatg.dw[,sp] > 0

m=apply(table(mh[fn]>0.,pmin(3,orth.tgcatg.dw[fn,'mouse'])),2,function(x)my.binom.test(x[2:1]))
c=apply(table(ch[fn]>0.,pmin(3,orth.tgcatg.dw[fn,'chicken'])),2,function(x)my.binom.test(x[2:1]))

m=apply(table(mb[fn]>0.,pmin(3,orth.tgcatg.dw[fn,'mouse'])),2,function(x)my.binom.test(x[2:1]))
c=apply(table(cb[fn]>0.,pmin(3,orth.tgcatg.dw[fn,'chicken'])),2,function(x)my.binom.test(x[2:1]))


plotArea(0:3,t(c),col='blue',new=T,ylim=c(0.18,0.9))
plotArea(0:3,t(m),col='red',new=F)


getConsByMot = function(s1,s2,t,mot,f,use.dev=FALSE,dpsi.thr=0,comp.dir=`>`){
	if(use.dev){
		up1 = age.dpsi[[s1]][f,t]
		up2 = age.dpsi[[s2]][f,t]
	}else{
		up1 = getAdultSpec(orth.seg.ad.tsm,t,cnt.tis,border.stages,s1)[f]
		up2 = getAdultSpec(orth.seg.ad.tsm,t,cnt.tis,border.stages,s2)[f]
	}
	up1f = comp.dir(up1,dpsi.thr)
	up2f = comp.dir(up2,dpsi.thr)
	
	m1 = mot[f,s1]>0
	m2 = mot[f,s2]>0
	ft = function(p1,p2,f){
		f[is.na(f+p1+p2)]  = FALSE
		r = fisher.test(factor(p1[f],levels = c(TRUE,FALSE)),factor(p2[f],levels = c(TRUE,FALSE)))
		c(r$estimate,
			r$conf.int,
			r$p.value,
			total=sum(f),
			up1=sum(p1[f]),
			up2=sum(p2[f]),
			overlap=sum(p1[f] & p2[f]),
			mean.dpsi1=mean(up1[f & p1]),
			mean.dpsi2=mean(up2[f & p2]),
			median.dpsi1=median(up1[f & p1]),
			median.dpsi2=median(up2[f & p2]))
		
	}
	# rbind('m00'=ft(up1f,up2f,!m1 & !m2),
	# 			'm01'=ft(up1f,up2f,!m1 &  m2),
	# 			'm10'=ft(up1f,up2f, m1 & !m2),
	# 			'm11'=ft(up1f,up2f, m1 &  m2))
	m1,p1==p2
}

t=getConsByMot('mouse','chicken','heart',orth.tgcatg.up,fn,use.dev = F,0,`>`)
t

f = ch !=0 | mh!=0
fisher.test(table(orth.tgcatg.up[f,'mouse']>0,sign(ch[f])==sign(mh[f])))


f = ch !=0 | mh!=0
fisher.test(table(orth.tgcatg.up[f,'mouse']>0,sign(ch[f])==sign(mh[f])))

#plot(mh[!mm],ch[!mm])
thr=0.3
table(mh>thr)
fisher.test(table(chi.devAS=ch[mh>thr]>thr,mouse.mot=mm[mh>thr]))

# check tham mammal heart exclusion is less conserved in chicken if have gcatg
sp = 'mouse'
mh = getAdultSpec(orth.seg.ad.tsm,'heart',cnt.tis,border.stages,sp,fun = min) 
ch = getAdultSpec(orth.seg.ad.tsm,'heart',cnt.tis,border.stages,'chicken',fun = min)
mm = orth.tgcatg.up[,sp] > 0
#plot(mh[!mm],ch[!mm])
thr=0.3
table(mh < -thr)
fisher.test(table(chi.devAS=ch[mh < -thr]< -thr,mouse.mot=mm[mh< -thr]))



f = getAdultSpec(orth.seg.ad.tsm,'heart',cnt.tis,border.stages,'mouse') > 0.3 & apply(orth.tgcatg.dw[,-7]>0,1,sum) == 6
table(f)
par(mfrow=c(3,3),tck=-0.01,mgp=c(1.1,0.2,0),mar=c(6,3,1.5,0),oma=c(0,0,2,1))
for(s in names(orth.seg.ad))
	plotTissueAgeProile(apply(orth.seg.ad[[s]]$ir[f,],2,mean,na.rm=T),meta,age.axis = 'rank',ylim=c(0,1),main=s)


r = NULL
b.thr=0.5
h.thr=0.2
f = getAdultSpec(orth.seg.ad.tsm,'brain',cnt.tis,border.stages,'chicken') > b.thr & getAdultSpec(orth.seg.ad.tsm,'heart',cnt.tis,border.stages,'chicken') < 0.1
table(f)
x = lapply(orth.seg.ad.tsm,function(x)x[f,])

for(s in names(orth.seg.ad)[c(-2,-7)]){
	b = getAdultSpec(x,'brain',cnt.tis,border.stages,s)
	h = getAdultSpec(x,'heart',cnt.tis,border.stages,s)
	t = orth.tgcatg.dw[f,'chicken']>0 & orth.tgcatg.dw[f,s]>0
	r = rbind(r,data.frame(species=s,motif='++',brain=sum(t&b>b.thr),brain.heart=sum(t&b>b.thr&h>h.thr)))
	t = orth.tgcatg.dw[f,'chicken']>0 & !orth.tgcatg.dw[f,s]>0
	r = rbind(r,data.frame(species=s,motif='+-',brain=sum(t&b>b.thr),brain.heart=sum(t&b>b.thr&h>h.thr)))
	t = !orth.tgcatg.dw[f,'chicken']>0 & orth.tgcatg.dw[f,s]>0
	r = rbind(r,data.frame(species=s,motif='-+',brain=sum(t&b>b.thr),brain.heart=sum(t&b>b.thr&h>h.thr)))
	t = !orth.tgcatg.dw[f,'chicken']>0 & !orth.tgcatg.dw[f,s]>0
	r = rbind(r,data.frame(species=s,motif='--',brain=sum(t&b>b.thr),brain.heart=sum(t&b>b.thr&h>h.thr)))
}
z=sapply(split(r[,-1:-2],r$motif),function(x)apply(x,2,sum))
z[2,]/z[1,]



a = sapply(orth.seg.ad,function(x)x$seg$type=='ALT')
f = apply(a[,c(-2,-4)], 1, sum)==5
table(f)

orth.dws.hex = array(FALSE,dim=c(sum(f),nrow(hex.dws.age02sgn$up$or),5),dimnames=list(rownames(orth.seg.ad$human$seg)[f],rownames(hex.dws.age02sgn$up$or),names(orth.seg.ad)[-c(2,4)]))
for(s in dimnames(orth.dws.hex)[[3]]){
	print(s)
	t=orth.fa[[s]][rownames(orth.seg.ad[[s]]$seg)[f]]
	t=substr(t,nchar(t)-195,nchar(t))
	for(h in dimnames(orth.dws.hex)[[2]])
		orth.dws.hex[,h,s] = grepl(h,t,fixed = T)
}
	
z=t(apply(orth.dws.hex,2,function(x){
	r = table(factor(apply(x,1,sum),levels = 0:5))
#	r/dbinom(0:5,ncol(x),mean(x))
	}))

hist(z[,1])
hist(z[,5],ylim=c(0,200),100)
table(z[,5]>50)

plot(1,t='n',xlim=c(0,5),ylim=c(1,max(z[,4:5])),log='')
#plot(1,t='n',xlim=c(0,5),ylim=c(1,max(z)),log='y')
for(i in 1:nrow(z))
		lines(0:5,z[i,],col=ifelse(z[i,5]>=z[i,6],'#00000020','red'))
z[order(z[,6]-z[,5],decreasing = T)[1:10],]



boxplot(log10(z/sum(f)))
hist(log10(z[,2]/sum(f)),100)

seg2ens = readRDS('Rdata/seg2ens.Rdata')

rbfox1.sids = names(seg2ens$human)[sapply(seg2ens$human,function(x)'ENSG00000078328' %in% x)]
rbfox1.sids.orth = intersect(rbfox1.sids,rownames(orth.seg.ad$human$seg))


par(mfrow=c(3,3),tck=-0.01,mgp=c(1.1,0.2,0),mar=c(6,3,1.5,0),oma=c(0,0,2,1))
id = which(rownames(orth.seg.ad$human$seg) == rbfox1.sids.orth[7])
for(s in names(orth.seg.ad))
	plotTissueAgeProile(orth.seg.ad[[s]]$ir[id,],meta,ylim=c(0,1),age.axis = 'rank',main=s)
#1 - const
#2 - chicken brain skip
#3 - mammal bit heart skipp
#4 - const
#5 hmm chi brain skipp
#6 const

h = readRDS("Rdata/human.as.u.all.Rdata")
h = h[h$seg$gene_id=='hum.23240',]
rbfox1.sids.h = intersect(rbfox1.sids,rownames(psi.tsm$human))
par(mfrow=c(3,3),tck=-0.01,mgp=c(1.1,0.2,0),mar=c(6,3,1.5,0),oma=c(0,0,2,1))
for(s in rbfox1.sids.h)
	plotTissueAgeProile(psi.tsm$human[s,],meta.tsm,ylim=c(0,1),age.axis = 'rank',main=s)
#hum.23240.s32 - heart skipping
anns$human['hum.23240.s32',]
plotTissueAgeProile(h$ir['hum.23240.s32',],meta,ylim=c(0,1),age.axis = 'rank',main='hum.23240.s32')
plotTissueAgeProile(h$ir['hum.23240.s33',],meta,ylim=c(0,1),age.axis = 'rank',main='hum.23240.s33')
plotTissueAgeProile(psi.tsm$chicken['chi.8331.s6',],meta.tsm,ylim=c(0,1),age.axis = 'rank',main='chi.8331.s6')
plotTissueAgeProile(psi.tsm$chicken['chi.8331.s4',],meta.tsm,ylim=c(0,1),age.axis = 'rank',main='chi.8331.s6') # should be an ortholog of hum.23240.s35
#exon 10, Fox-1 C-terminal domain, it is not 19th exon from Baralle review (19th exon is hum.23240.s35)
#chi.8331.s4   TGTTGTATACCAGGATGGATTTTATGGTGCAGACATTTAT
#hum.23240.s32 TGTTGTTTACCAGGATGGATTTTATGGTGCAGACATTTAT     VVYQDGFYGADIY 303-316
#hum.23240.s33 AGTAGTGTATCAAGAGCCTGTGTATGGCAATAAATTGCTGCAG  VVYQEPVYGNKLLQ 318-332
#chi:          AGCATTATTGCAAGAGCCTGTGTATGGCAATAAGTTACTACAG (14:10884946-10884989)
plot(h$ir['hum.23240.s32',],h$ir['hum.23240.s33',])

hmo.seg.ad = readRDS('Rdata/hmo.seg.ad.Rdata')
which(rownames(hmo.seg.ad$human$seg)=='hum.23240.s32')

pdf('figures/paper.figures/3.2/examples/rbfox1.mammal-mex.psi.pdf',w=9,h=6)
par(mfcol=c(2,3),tck=-0.01,mgp=c(1.1,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
for(s in names(hmo.seg.ad)){
	plotTissueAgeProile(hmo.seg.ad[[s]]$ir[15194,],meta,ylim=c(0,1),age.axis = 'rank',main=paste(s,'ancient exon'),ylab='PSI')
	plotTissueAgeProile(hmo.seg.ad[[s]]$ir[15193,],meta,ylim=c(0,1),age.axis = 'rank',main=paste(s,'mammal exon'),ylab='PSI')
}
dev.off()

sapply(colnames(orth.seg.ad.all.id),function(s)seg2ens[[s]][[orth.seg.ad.all.id['hum.23240.s30',s]]])
orth.ens.genes['ENSG00000078328',]

orth.seg.ad.all.id = readRDS('Rdata/orth.seg.ad.all.id.Rdata')
rbfox1.seg = all.anns$human[all.anns$human$gene_id=='hum.23240'&all.anns$human$sites=='ad',]
rbfox1.seg$has.orth = rownames(rbfox1.seg) %in% orth.seg.ad.all.id[,1]
orth.seg.ad.all.id[rownames(rbfox1.seg)[rbfox1.seg$has.orth],]
all.anns$chicken[all.anns$chicken$gene_id=='chi.8331'&all.anns$chicken$sites=='ad',]

vastdb=read.table('input/VASTDB_SAMPLE_INFO_Gga61_galGal3.tab.gz',sep='\t',header = T)
cvh = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/chi.muscle/',vastdb$SRA_ID[vastdb$VastDB_subgroup=='Heart'],'.bam')
cvm = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/chi.muscle/',vastdb$SRA_ID[vastdb$VastDB_subgroup=='Muscle'],'.bam')

hh = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/',meta$fname[meta$species=='human' & meta$tissue=='heart'],'.bam')
hb = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/human/',meta$fname[meta$species=='human' & meta$tissue=='brain'],'.bam')
oh = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/opossum/',meta$fname[meta$species=='opossum' & meta$tissue=='heart'],'.bam')
ob = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/opossum/',meta$fname[meta$species=='opossum' & meta$tissue=='brain'],'.bam')
ch = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/chicken/',meta$fname[meta$species=='chicken' & meta$tissue=='heart'],'.bam')
cb = paste0('/home/mazin/skoltech/projects/evo.devo/processed/mapping/hisat2.s/chicken/',meta$fname[meta$species=='chicken' & meta$tissue=='brain'],'.bam')



pdf('figures/paper.figures/3.2/examples/rbfox1.mammal-mex.pdf',w=8,h=14)
par(mfrow=c(8,1),tck=-0.01,mgp=c(1.1,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
getReadCoverage(hh,'16',7703817,7726840,-1,T,min.junc.cov = 20,reverse=F,ylab='coverage',main='Human heart',xlab='position (nt)')
abline(v=c(7714931/2 + 7714970/2,7721559/2+7721601/2),col='blue',lty=3)
getReadCoverage(hb,'16',7703817,7726840,-1,T,min.junc.cov = 20,reverse=F,ylab='coverage',main='Human brain',xlab='position (nt)')
abline(v=c(7714931/2 + 7714970/2,7721559/2+7721601/2),col='blue',lty=3)

getReadCoverage(oh,'6',116107221,116134386,1,T,min.junc.cov = 20,reverse=T,ylab='coverage',main='Opossum heart',xlab='position (nt)')
abline(v=c(116112109/2 + 116112151/2,116120835/2+116120874/2),col='blue',lty=3)
getReadCoverage(ob,'6',116107221,116134386,1,T,min.junc.cov = 20,reverse=T,ylab='coverage',main='Opossum brain',xlab='position (nt)')
abline(v=c(116112109/2 + 116112151/2,116120835/2+116120874/2),col='blue',lty=3)

getReadCoverage(ch,'14',10880884,10901507,1,T,min.junc.cov = 5,reverse=T,ylab='coverage',main='Chicken heart',xlab='position (nt)')
abline(v=c(10889537/2 + 10889576/2,10884946/2+10884989/2),col='blue',lty=3)
getReadCoverage(cb,'14',10880884,10901507,1,T,min.junc.cov = 20,reverse=T,ylab='coverage',main='Chicken brain',xlab='position (nt)')
abline(v=c(10889537/2 + 10889576/2,10884946/2+10884989/2),col='blue',lty=3)

getReadCoverage(cvh,'14',10880884,10901507,NA,T,min.junc.cov = 0,reverse=T,ylab='coverage',main='Chicken vastdb heart',xlab='position (nt)')
abline(v=c(10889537/2 + 10889576/2,10884946/2+10884989/2),col='blue',lty=3)
getReadCoverage(cvm,'14',10880884,10901507,NA,T,min.junc.cov = 20,reverse=T,ylab='coverage',main='Chicken vastdb muscle',xlab='position (nt)')
abline(v=c(10889537/2 + 10889576/2,10884946/2+10884989/2),col='blue',lty=3)

dev.off()


all.anns$chicken[all.anns$chicken$gene_id=='chi.8331' & all.anns$chicken$sites=='ad',]

# compare developmental changes and adult state ####
splitTissuesAtStage = function(p,st){
	st = st[!is.na(st)]
	sp = strsplit(colnames(p)[1],' ')[[1]][1]
	p = p[,paste(sp,names(st),st)]
	tis = substr(names(st),1,1)
	do.call(rbind,apply(p,1,function(x){
		o = order(x)
		diff = x[o[-1]] - x[o[-length(o)]]
		split=order(diff,decreasing = T)[1]
		data.frame(low=paste(tis[sort(o[1:split])],collapse=''),high=paste(tis[sort(o[-(1:split)])],collapse=''),diff=diff[split],diff2=sort(diff,decreasing = T)[2])
		}))
}

mts = splitTissuesAtStage(psi.tsm$mouse,border.stages$mouse[,2])
f = anns$mouse$sites=='ad' & age03$mouse[,'heart'] == 'u'/
sort(table(paste0(mts$low,'->',mts$high)[f]))
f = anns$mouse$sites=='ad' & mts$high %in% c('h','bh','ch','bch','ht','bcht','bht','cht','bc') & mts$diff>0.5
table(age03$mouse[f,'heart'],paste0(mts$low,'->',mts$high)[f])



hist(mts$diff2[f]/mts$diff[f],0:100/100)
hist(mts$diff[f],0:100/100)
hist(mts$diff2[f],0:100/100)

mts[anns$mouse$sites=='ad' & grepl('h',mts$high) & mts$diff>0.5 & age03$mouse[,'heart']=='n',]

mts[f & !is.na(mts$diff) & !is.na(mts$diff2) & mts$diff<0.01,]
plotTissueAgeProile(psi.tsm$mouse['mou.46249.s16',],meta.tsm)

sort(table(paste0(mad$low,'->',mad$high)[mad$diff>0.3 & anns$mouse$sites=='ad']),decreasing = T)[1:20]


# look on QKI ####
bu = rownames(age.segs$mouse)[age.segs$mouse[,'brain']=='u' & anns$mouse$sites=='ad']
bu = intersect(bu,names(fa$mouse))
length(bu)
buqki = bu[grep('actaac',substr(fa$mouse[bu],1,200))]
length(buqki)

hu = rownames(age.segs$mouse)[age.segs$mouse[,'heart']=='u' & anns$mouse$sites=='ad']
hu = intersect(hu,names(fa$mouse))
length(hu)
l = nchar(fa$mouse[hu])
huqki = bu[grep('actaac',substr(fa$mouse[hu],l-200+1,l))]
length(huqki)

length(unique(anns$mouse[buqki,'gene_id']))
length(unique(anns$mouse[huqki,'gene_id']))
length(intersect(unique(anns$mouse[buqki,'gene_id']),unique(anns$mouse[huqki,'gene_id'])))
