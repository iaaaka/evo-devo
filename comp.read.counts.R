options(stringsAsFactors = FALSE)
setwd('~/skoltech/projects/evo.devo/')
species = c(m='mouse',r='rat',b='rabbit',h='human',q='macaque',o='opossum',c='chicken')
rev.species = setNames(names(species),species)
# number of reads in Margarida bams

fn = paste('processed/fq.from.bams/stat/',toupper(substr(species,1,1)),substr(species,2,20),'.rc.from.bam',sep='')
bam.rc = do.call(rbind,lapply(fn,read.table))


bam.rc[,1] = sapply(strsplit(bam.rc[,1],'/',TRUE),function(x){x[length(x)]})
# parse file names
fnames = strsplit(bam.rc[,1],'.',TRUE)
table(sapply(fnames,function(x){x[length(x)-1]}))
table(sapply(fnames,function(x){x[5]}))
table(sapply(fnames,length))

p=lapply(fnames,function(x){
	if(length(x) %in% 5:6){
		r = c(x[1:3],'-','-')
	}else if(length(x) == 7){
		if(x[5] == '5'){
			r = c(x[1:4],'-')
			r[4] = paste(x[4],x[5],sep='.')
		}else
			r = x[1:5]
	}else if(length(x) == 8){
		r = x[c(1:4,6)]
		r[4] = paste(x[4],x[5],sep='.')
	}
	r
})

bam.rc = data.frame(do.call(rbind,p),bam.rc[,2],fname=bam.rc[,1])
colnames(bam.rc) = c('lib.id','species','tissue','age','sex','read.count','fname')
rownames(bam.rc) = bam.rc$lib.id
bam.rc[1:10,]
table(bam.rc$species)

# load number of rows in fq derived from bams
fn = paste('processed/fq.from.bams/stat/',species,'.wc',sep='')
fqn.rc = do.call(rbind,lapply(fn,read.table))
fqn.rc[,1] = sapply(strsplit(fqn.rc[,1],'/',TRUE),function(x){x[length(x)]})
fqn.rc[,1] = rownames(fqn.rc) = sapply(strsplit(fqn.rc[,1],'.',TRUE),'[',1)
colnames(fqn.rc) = c('lib.id','read.count')
fqn.rc$read.count = fqn.rc$read.count/4
dim(fqn.rc)
dim(bam.rc)
length(unique(rownames(fqn.rc),rownames(bam.rc)))

fqn.rc = fqn.rc[rownames(bam.rc),]
table(fqn.rc$read.count == bam.rc$read.count,bam.rc$species) # my fq are correct (compared to original bams), except four cases where I desided to use fq from HD
bam.rc[fqn.rc$read.count != bam.rc$read.count,]
fqn.rc[fqn.rc$read.count != bam.rc$read.count,]
# lib.id species  tissue  age    sex read.count                                    fname
# 1880sTS 1880sTS   Mouse   Brain 2wpb   Male   32000000 1880sTS.Mouse.Brain.2wpb.Male.sorted.bam # not all read were used
# 2684sTS 2684sTS     Rat   Heart   20   Male   30619512     2684sTS.Rat.Heart.20.Male.sorted.bam # some reads were missed in bam, one truncated
# 3607sTS 3607sTS     Rat   Heart   16 Female   34835228   3607sTS.Rat.Heart.16.Female.sorted.bam # some reads were missed in bam
# 4947sTS 4947sTS     Rat Decidua   10 Female   44032473 4947sTS.Rat.Decidua.10.Female.sorted.bam # some reads were missed in bam
# load output of bam2fq
bam.rc$fq.read.count = fqn.rc$read.count

# load numbers from hisat2.f logos
fn = paste('processed/mapping/hisat2.f/',species,'/mapping.stat',sep='')
hs.rc = do.call(rbind,lapply(fn,read.table))
rownames(hs.rc) = hs.rc[,1]
colnames(hs.rc) = c('lib.id','read.count','mapped.uniq','mapped.mult')
table(sort(rownames(bam.rc)) == sort(rownames(hs.rc)))
hs.rc = hs.rc[rownames(bam.rc),]
table(bam.rc[,'fq.read.count']==hs.rc$read.count)
boxplot((hs.rc$mapped.mult/hs.rc$read.count) ~ bam.rc$species) # means nothing because hisat produces wrong statistica
boxplot(((hs.rc$mapped.uniq+hs.rc$mapped.mult)/hs.rc$read.count) ~ bam.rc$species,las=2)

# second mapping
fn = paste('processed/mapping/hisat2.s/',species,'/mapping.stat',sep='')
hss.rc = do.call(rbind,lapply(fn,read.table))
rownames(hss.rc) = hss.rc[,1]
colnames(hss.rc) = c('lib.id','read.count','mapped.uniq','mapped.mult')
hss.rc = hss.rc[rownames(bam.rc),]
table(hss.rc$read.count == bam.rc$fq.read.count)
hist(hss.rc[,4] - hs.rc[,4])
hist(hss.rc[,3] - hs.rc[,3])
hist(hss.rc[,3]+hss.rc[,4] - hs.rc[,3] - hs.rc[,4])


#old mapping
fn = paste('processed/mapping/old.hisat2.s/',species,'/mapping.stat',sep='')
ohss.rc = do.call(rbind,lapply(fn,read.table))
rownames(ohss.rc) = ohss.rc[,1]
colnames(ohss.rc) = c('lib.id','read.count','mapped.uniq','mapped.mult')
ohss.rc = ohss.rc[rownames(bam.rc),]
plot(ohss.rc$mapped.uniq/ohss.rc$read.count,hss.rc$mapped.uniq/hss.rc$read.count,col=factor(bam.rc$species))
plot((ohss.rc$mapped.uniq+ohss.rc$mapped.mult)/ohss.rc$read.count,(hss.rc$mapped.uniq+hss.rc$mapped.mult)/hss.rc$read.count)
abline(a=0,b=1,col='red')

#### make commands for stringtie

bam.rc$stage = bam.rc$age

# make human stages according to Margarida suggestions
bam.rc$stage[bam.rc$species=='Human'] = sapply(bam.rc$age[bam.rc$species=='Human'],function(x){
	r = x
	if(substr(x,nchar(x),nchar(x))=='w')
		r = paste0(x,'pc')
	if(x %in% c('CS13','CS14'))
		r = '4wpc'
	if(x %in% c('CS15','CS16'))
		r = '5wpc'
	if(x %in% c('CS17','CS18'))
		r = '6wpc'
	if(x %in% c('CS19','CS20','CS21'))
		r = '7wpc'
	if(x %in% c('CS22','CS23'))
		r = '8wpc'
	if(grepl('[ymd]pb',x)){
		age = as.numeric(substr(x,1,nchar(x)-3))*switch(substr(x,nchar(x)-2,nchar(x)-2),y=1,m=1/12,d=1/365)
		if(age <= 100/365) 
			r = 'newborn'
		else if(age <= 9/12)
			r = 'infant' 
		else if(age <= 4)
			r = 'toddler' 
		else if(age <= 9)
			r = 'school' 
		else if(age <= 14)
			r = 'youngTeenager' # margarida uses it only for testis
		else if(age <= 19)
			r = 'teenager' # (Margarida uses teenager (no testis samples) and oldTeenager (only testis samples) with overlapping ages)
		else if(age <= 32)
			r = 'youngAdult'
		else if(age <= 41)
			r = 'youngMidAge'
		else if(age <= 54)
			r = 'olderMidAge'
		else if(age > 54)
			r = 'senior'
	}
	r
	})


# make macaque stages according to Margarida suggestions
bam.rc$stage[bam.rc$species=='Macaque'] = sapply(bam.rc$age[bam.rc$species=='Macaque'],function(x){
	r = x
	if(x %in% c('93','109','112','123'))
		r = paste0('e',x)
	if(x %in% c('129','130'))
		r = 'e130'
	if(grepl('[ymd]pb',x)){
		age = as.numeric(substr(x,1,nchar(x)-3))*switch(substr(x,nchar(x)-2,nchar(x)-2),y=1,m=1/12,d=1/365)
		if(age <= 7/365)
			r = 'P0'
		else if(age <= 30/365)
			r = 'P23'
		else if(age <= 0.5)
			r = 'P152'
		else if(age <= 1)
			r = 'P365'
		else if(age <= 3)
			r = 'P1095'
		else if(age <= 9)
			r = 'P3285'
		else if(age <= 15)
			r = 'P5475'
		else if(age > 15)
			r = 'P8030'
	}
	r
})
table(bam.rc$age[bam.rc$species=='Human'],bam.rc$stage[bam.rc$species=='Human'])
table(bam.rc$age[bam.rc$species=='Macaque'],bam.rc$stage[bam.rc$species=='Macaque'])
bam.rc = cbind(bam.rc,hss.rc[,-1])
#write.csv(bam.rc,'input/sample.info.csv')

# Make commands for stringtie
bam.rc = read.csv('input/sample.info.csv',row.names = 1)
t = bam.rc[,c(2,3,7,9)]
t$sta = paste(t$species,t$tissue,t$stage,sep='_')
t = split(t,t$sta)
length(t)
table(sapply(t,nrow))
fnames = sapply(t,function(x){paste(paste("../../mapping/hisat2.s/",tolower(x$species),"/",gsub('.sorted','',x$fname,fixed = TRUE),sep=''),collapse=' ')})
sp = tolower(sapply(t,function(x){x$species[1]}))

sam.cnt = sapply(t,nrow)

tt = paste(ifelse(sam.cnt==1,'cat ','samtools merge - '),fnames," | stringtie - --bam -f 0.1 -p 12 -j 3 -g 10 --rf -l ",substr(sp,1,3)," -o ",sp,"/",names(t),'.gtf',sep='')
#t=paste("samtools merge - /home/mazin/evo.devo/mapping/hisat2.s/",tolower(t[,1]),"/*",t[,1],'.',t[,2],'.',t[,3],"*.bam | stringtie - --bam -f 0.2 -p 12 -j 3 -g 10 -l ",rev.species[t[,1]]," -o ",t[,1],"/",t[,1],'.',t[,2],'.',t[,3],'.gtf',sep='')
write.table(tt,quote = F,col.names = F,row.names = F,'processed/annotation/all.species/stringtie.commands')

