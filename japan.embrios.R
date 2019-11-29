options(stringsAsFactors = FALSE)
source('code/r.functions/load.all.data.F.R')
library(SAJR)
library(xlsx)
library(doMC)
anns = readRDS('Rdata/anns.Rdata')
species = readRDS('Rdata/species.Rdata')
meta = readRDS('Rdata/main.set.meta.Rdata')
psi.tsm = readRDS('Rdata/psi.tsm.Rdata')
meta.tsm = readRDS('Rdata/meta.tsm.Rdata')
border.stages = readRDS('Rdata/border.stages.Rdata')
chi.je.as = readRDS('Rdata/japan.embrio/chi.je.as.filtered.Rdata')
mou.je.as = readRDS('Rdata/japan.embrio/mou.je.as.filtered.Rdata')
params = list(species.col=c(human='red',macaque='orange',mouse='gray',rat='black',rabbit='brown',opossum='magenta',chicken='violet'),tissue.col=c(brain="#3399CC",cerebellum="#33CCFF",heart="#CC0000",kidney="#CC9900",liver="#339900",ovary="#CC3399",testis="#FF6600"))
params$species.pch=c(human=15,macaque=0,mouse=19,rat=1,rabbit=18,opossum=17,chicken=8)
meta.je = readRDS('Rdata/japan.embrio/meta.je.Rdata')
orth.seg.ad.all.id = readRDS('Rdata/orth.seg.ad.all.id.Rdata')


mou.je.as.tsm = calcMeanCols(mou.je.as$ir,meta.je$stage[meta.je$species=='mouse'])
chi.je.as.tsm = calcMeanCols(chi.je.as$ir,meta.je$stage[meta.je$species=='chicken'])

# meta.je = read.xlsx('input/japan.embrio.TableS1.sample_information.xlsx',sheetIndex = 1)
# meta.je = meta.je[substr(meta.je$RNAseq.data,1,2) %in% c('Mm','Gg'),-1]
# colnames(meta.je) = c('fname','read','lib.prep','start.mat','species','id.old','number.of.embryos.pooled')
# 
# meta.je$species = ifelse(substr(meta.je$fname,1,2)=='Mm','mouse','chicken')
# meta.je$fname = gsub('-','_',meta.je$fname)
# meta.je$fname = gsub('Mm_morulal','Mm_morula',meta.je$fname)
# meta.je$id = gsub('.fastq.gz','',meta.je$fname)
# meta.je$age.rank = 1:nrow(meta.je)
# meta.je$stage = substr(meta.je$id,1,nchar(meta.je$id)-2)
# meta.je$age.rank = rank(sapply(split(meta.je$age.rank,meta.je$stage),min))[meta.je$stage]
# meta.je$age.rank[meta.je$species=='chicken'] = meta.je$age.rank[meta.je$species=='chicken'] - min(meta.je$age.rank[meta.je$species=='chicken'])+1
# meta.je = meta.je[meta.je$id.old!='6-8cell_3',] #absent in Dimas data
# meta.je$days = as.numeric(gsub('Mm_E','',meta.je$stage))
# meta.je$days[meta.je$stage=='Gg_Prim'] = 0.5
# meta.je$days[meta.je$stage=='Gg_HH6'] = 1
# meta.je$days[meta.je$stage=='Gg_HH8'] = 27/24
# meta.je$days[meta.je$stage=='Gg_HH11'] = 43/24
# meta.je$days[meta.je$stage=='Gg_HH14'] = 52/24
# meta.je$days[meta.je$stage=='Gg_HH16'] = 54/24
# meta.je$days[meta.je$stage=='Gg_HH19'] = 3.3
# meta.je$days[meta.je$stage=='Gg_HH21'] = 3.5
# meta.je$days[meta.je$stage=='Gg_HH24'] = 4.5
# meta.je$days[meta.je$stage=='Gg_HH28'] = 5.8
# meta.je$days[meta.je$stage=='Gg_HH32'] = 7.5
# meta.je$days[meta.je$stage=='Gg_HH34'] = 8
# meta.je$days[meta.je$stage=='Gg_HH38'] = 12
# meta.je$days[meta.je$stage=='Mm_2cell'] = 0.5
# meta.je$days[meta.je$stage=='Mm_6_8cell'] = 2
# meta.je$days[meta.je$stage=='Mm_morula'] = 3
# meta.je$days[meta.je$stage=='Mm_blastocyst'] = 4
# meta.je$days  = round(meta.je$days,digits = 1)
# rownames(meta.je) = meta.je$id
#saveRDS(meta.je,'Rdata/japan.embrio/meta.je.Rdata')

# chi.je.as = loadSAData('processed/annotation/all.species/merged/chicken.sajr',paste0('processed/sajr/uniq.japan.embrio/',meta.je$id[meta.je$species=='chicken']),meta.je$id[meta.je$species=='chicken'])
# mou.je.as = loadSAData('processed/annotation/all.species/merged/mouse.sajr'  ,paste0('processed/sajr/uniq.japan.embrio/',meta.je$id[meta.je$species=='mouse'])  ,meta.je$id[meta.je$species=='mouse'])
# 
# saveRDS(chi.je.as,'Rdata/japan.embrio/chi.je.as.Rdata')
# saveRDS(mou.je.as,'Rdata/japan.embrio/mou.je.as.Rdata')

chi.je.as = chi.je.as[rownames(anns$chicken),]
mou.je.as = mou.je.as[rownames(anns$mouse),]

# saveRDS(chi.je.as,'Rdata/japan.embrio/chi.je.as.filtered.Rdata')
# saveRDS(mou.je.as,'Rdata/japan.embrio/mou.je.as.filtered.Rdata')
m = readRDS('Rdata/mouse.as.u.filtered.Rdata')
mpsi = cbind(m$ir,mou.je.as$ir)[m$seg$sites=='ad',]
mcor = cor(mpsi,u='p')


pdf('figures/MDSs/mouse.with.japan.embryo.pdf',w=11,h=5.5)
f = 1:316
#f = c(which(meta[colnames(m$ir),'tissue'] %in% c('testis','brain','heart')),317:(317+18))
mmds = cmdscale(1-mcor[f,f],k=2)
mcol = c(meta[colnames(m$ir),'col'],rep('cyan',ncol(mou.je.as$ir)))
mcex = meta.je$age.rank[meta.je$species=='mouse']
mcex = c(meta[colnames(m$ir),'cex'],0.2+(mcex-1)/(max(mcex)-1)*1.6)*2
par(mfrow=c(1,2),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,0,1))
plot(mmds[,1:2],col=mcol[f],pch=19,cex=mcex[f],xlab='Dim 1',ylab='Dim 2')
f = 1:(355-22)
mmds = cmdscale(1-mcor[f,f],k=2)
mcol = c(meta[colnames(m$ir),'col'],rep('cyan',ncol(mou.je.as$ir)))
mcex = meta.je$age.rank[meta.je$species=='mouse']
mcex = c(meta[colnames(m$ir),'cex'],0.2+(mcex-1)/(max(mcex)-1)*1.6)*2
plot(mmds[,1:2],col=mcol[f],pch=19,cex=mcex[f],xlab='Dim 1',ylab='Dim 2')
legend('topleft',pch=19,col=c(params$tissue.col,'cyan'),legend=c(names(params$tissue.col),'embryo:2cell-9.5dpc'),bty='n')
dev.off()



hist(mou.je.as.tsm[,'Mm_E9.5']-mou.je.as.tsm[,'Mm_2cell'])
which(mou.je.as.tsm[,'Mm_E9.5']-mou.je.as.tsm[,'Mm_2cell'] > 0.9 & m$seg$sites=='ad')

plotTissueAgeProile(m$ir[11032,],meta)

z=apply(mou.je.as.tsm[m$seg$sites=='ad',],2,function(x){
	x=x[!is.na(x)];
	c(total=length(x),psi0=sum(x<0.2),psi05=sum(x>=0.1 & x<=0.9),psi1=sum(x>0.9))
	})
plot(z[3,]/z[1,])


mdpsi = getAgeASchanges(psi.tsm,meta.tsm,0.1,border.stages,'mouse',get.dPSI=T)
f = which(m$seg$sites=='ad')
f=sample(f,1000)
pairs(cbind(mou.je.as.tsm[f,'Mm_E9.5']-mou.je.as.tsm[f,'Mm_2cell'],mdpsi[f,]),pch='.')



m.alt.pr = apply(psi.tsm$mouse,2,function(x){x=x[!is.na(x)];sum(x>0.1 & x<0.9)/length(x)})
pdf('figures/japan.embrio/AS.compl.for.japan.embryo.pdf',w=5,h=5)
par(mfrow=c(1,1),tck=-0.01,mgp=c(1.3,0.2,0),mar=c(3.5,2.5,1.5,0),oma=c(0,0,0,1))
plotTissueAgeProile(m.alt.pr*100,meta.tsm,age.axis = 'rank',xlim=c(-6,14),ylim=c(0.22,0.6)*100,ylab='% of exons in [0.1,0.9]')
points(c(-6:-1,0:10/2),z[3,]/z[1,]*100,pch=19,col='cyan',t='b')
axis(1,-6:0,c('2cells','8cells','morula','blast.','7.5','8.5','9.5'),las=3)
dev.off()


# check correlation of full embrio with per tissue data
table(meta$stage[meta$species=='mouse'])
ed.by.stage = list()

ed.by.stage$all = calcMeanCols(psi.tsm$mouse,meta.tsm[colnames(psi.tsm$mouse),'days'])
ed.by.stage$all = ed.by.stage$all[,order(as.numeric(colnames(ed.by.stage$all)))]

f = !grepl('cerebellum|kidney',colnames(psi.tsm$mouse))
ed.by.stage$bhlot = calcMeanCols(psi.tsm$mouse[,f],meta.tsm[colnames(psi.tsm$mouse)[f],'days'])
ed.by.stage$bhlot = ed.by.stage$bhlot[,order(as.numeric(colnames(ed.by.stage$bhlot)))]

for(t in c('brain','heart','liver','ovary','testis')){
	f = grepl(t,colnames(psi.tsm$mouse))
	ed.by.stage[[t]] = calcMeanCols(psi.tsm$mouse[,f],meta.tsm[colnames(psi.tsm$mouse)[f],'days'])
	ed.by.stage[[t]] = ed.by.stage[[t]][,order(as.numeric(colnames(ed.by.stage[[t]])))]
}

pdf('figures/japan.embrio/ed2je.pdf',w=15,h=15)
f = anns$mouse$sites=='ad' #& (apply(mou.je.as.tsm,1,sd,na.rm=T) > 0.1 | apply(ed.by.stage.bhlot,1,sd,na.rm=T) > 0.1)
par(mfrow=c(3,3),las=2,mar=c(5,5,2,0),oma=c(0,0,1,1),tck=-0.01,mgp=c(1.3,0.2,0))
for(n in names(ed.by.stage)){
	cc = cor(mou.je.as.tsm[f,],ed.by.stage[[n]][f,],u='p')
	imageWithText(cc,digits = 2,ylab='',xlab='',col=getPal(c('blue','white','red'),100),xaxlab = gsub('Mm_','',rownames(cc)),main=n,zlim=c(0.68,0.94))
}
mtext('Pearson correlation',3,outer = T,las=1)
par(mfrow=c(3,3),las=2,mar=c(5,5,2,0),oma=c(0,0,1,1),tck=-0.01,mgp=c(1.3,0.2,0))
for(n in names(ed.by.stage)){
	cc = cor(mou.je.as.tsm[f,],ed.by.stage[[n]][f,],u='p',m='sp')
	imageWithText(cc,digits = 2,ylab='',xlab='',col=getPal(c('blue','white','red'),100),xaxlab = gsub('Mm_','',rownames(cc)),main=n,zlim=c(0.54,0.81))
}
mtext('Spearman correlation',3,outer = T,las=1)
dev.off()
# max=apply(je2ed.bhlot,2,which.max)
# lines(max,1:ncol(je2ed.bhlot))
# lines(1:nrow(je2ed.bhlot),max)


ed.by.stage.chi = calcMeanCols(psi.tsm$chicken,meta.tsm[colnames(psi.tsm$chicken),'days'])
f = anns$chicken$sites=='ad'
je2ed.chi = cor(chi.je.as.tsm[f,],ed.by.stage.chi[f,],u='p')
s2d = setNames(meta.je$days,meta.je$stage)
imageWithText(je2ed.chi,digits = 2,ylab='',xlab='',col=getPal(c('blue','white','red'),100),xaxlab = s2d[rownames(je2ed.chi)])


#embrio mouse to chicken
f = orth.seg.ad.all.id[,'mouse'] %in% rownames(psi.tsm$mouse) & orth.seg.ad.all.id[,'chicken'] %in% rownames(psi.tsm$chicken)
table(f)

mou2chi.je = cor(mou.je.as.tsm[orth.seg.ad.all.id[f,'mouse'],],chi.je.as.tsm[orth.seg.ad.all.id[f,'chicken'],],u='p')
par(las=2,mar=c(7,7,2,1))
imageWithText(mou2chi.je,col=getPal(c('blue','white','red'),100))

#test age
registerDoMC(3)
f = meta.je[colnames(mou.je.as$ir),'days']<11
m = meta.je[colnames(mou.je.as$ir)[f],]
p = fitSAGLM(mou.je.as[,f],formula(x ~ days + I(days^2) + I(days^3)),m,0.05,.parallel = TRUE,return.pv=TRUE)
p = p[,-1]
p[apply(!is.na(mou.je.as$ir[,f]),1,mean) < 0.6,] = NA
p = apply(p,2,p.adjust,m='BH')
p = apply(p,1,min)
mou.ed.dev = data.frame(qv=p,dpsi=mou.je.as.tsm[,"Mm_E10.5"]-mou.je.as.tsm[,"Mm_2cell"])
hist(mou.ed.dev$dpsi[mou.ed.dev$qv>0.05],-50:50/50)
#saveRDS(mou.ed.dev,'Rdata/japan.embrio/devAS.Rdata')
mou.ed.dev = readRDS('Rdata/japan.embrio/devAS.Rdata')
per.tissue.age.qv = readRDS('Rdata/per.tissue.age.qv.Rdata')
DPSI=0.2  # qv <0.05 and dPSI > 0.2
s = 'mouse'
age.segs = lapply(rownames(species),function(s)getAgeASchanges(psi.tsm,meta.tsm,DPSI,border.stages,s)[anns[[s]]$sites=='ad',])
names(age.segs) = rownames(species)
for(s in names(age.segs)) age.segs[[s]][age.segs[[s]] != '-' & (is.na(per.tissue.age.qv[[s]]) | per.tissue.age.qv[[s]]>0.05)[rownames(age.segs[[s]]),colnames(age.segs[[s]])]] = 'n'

p = mou.ed.dev[rownames(age.segs$mouse),]
age.segs$mouse = cbind(embryo='-',age.segs$mouse)
age.segs$mouse[!is.na(p$qv),'embryo'] = 'n'
age.segs$mouse[!is.na(p$qv) & p$qv<0.05 & p$dpsi>  DPSI,'embryo'] = 'u'
age.segs$mouse[!is.na(p$qv) & p$qv<0.05 & p$dpsi< -DPSI,'embryo'] = 'd'

#age.ad.over.=age.ad.over
source('code/r.functions/paper.figures.F.R')
age.ad.over = lapply(age.segs,function(x){
	x[x=='-'] = NA
	x = cbind(x=='u',x=='d')
	#x = x[apply(x,1,sum,na.rm=T)>0,]
	colnames(x) = paste(colnames(x),rep(c('up','dw'),each=ncol(x)/2))
	caclSegOverlap(x)
})

pdf('figures/japan.embrio/devASoverlap.pdf',w=10,h=10)
par(mfrow=c(1,1),mar=c(5,5,2,0),oma=c(0,0,1,1),tck=-0.01,mgp=c(1.3,0.2,0))
plotAgeSegOverlap(age.ad.over$mouse,main=paste('Overlap of devAS exons across',s,'tissues'))
dev.off()

age.segs$mouse[age.segs$mouse[,1]=='u' & age.segs$mouse[,2]=='d',]
plotTissueAgeProile(psi.tsm$mouse['mou.20277.s14',],meta.tsm,age.axis = 'rank',xlim=c(-6,14))
points(c(-6:-1,0:10/2),mou.je.as.tsm['mou.20277.s14',],pch=19,col='cyan',t='b')
seg2ens$mouse['mou.3552.s40']
