import gzip
import re

mar = 200
chrs = [str(i) for i in range(1,23)]
chrs.append('X')
chrs.append('Y')
chrs.reverse()
workdir='/uge_mnt/home/mazin/'
segs = map(str.split, open(workdir+'projects/evo.devo/processed/gnomad201/human.ad.tab'))

def loadWigFix(f):
	file = gzip.open(f)
	r = [0]
	for l in file:
		l = l.rstrip()
		if l[0] == 'f':
			m = re.search('start=(\d+) ',l)
			if m:
				i = int(m.group(1))
				r.extend([0]*(i-len(r)))
		else:
			r.append(float(l))
	return r

def getProfSeg(p,f,t,rev):
	r = p[f:(t+1)]
	if rev:
		r.reverse()
	return ','.join([str(i) for i in r])
	

def getProfileForChr(chr):
	wig = loadWigFix(workdir+'projects/evo.devo/raw/primate.Phastcons46/chr'+str(chr)+'.phastCons46way.primates.wigFix.gz')
	print('loaded')
	ss = [s[0]+"\t"+getProfSeg(wig,int(s[2])-mar,int(s[3])+mar,s[4]=='-1') for s in segs if s[1] == chr]
	return ss


f = open(workdir+'/projects/evo.devo/processed/ad.phastcons', 'w')

for c in chrs:
	print(c)
	pc=getProfileForChr(c)
	for p in pc:
		f.write(p + "\n")

f.close()
