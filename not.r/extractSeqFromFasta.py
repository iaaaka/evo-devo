#!/usr/bin/python
import os
import sys

#usage:cat *gtf | python  *.py  fasta [output format: empty for tab, 'fa' for fasta] > .fa
path2fa = sys.argv[1]
complement={'a':'t','t':'a','g':'c','c':'g','n':'n','A':'T','T':'A','G':'C','C':'G','N':'N'}

def getSeq(fa,seqname,start,stop,strand):
	start -= 1 #since gtf coordinates are 1-based
	if seqname in fa:
		r = fa[seqname][start:stop]
	else:
		return('-')
	if strand == '-':
		r = ''.join([complement[i] for i in r[::-1]])
	return(r)

def readFasta(f):
	fa = open(f,buffering=1000000)
	r = {}
	name = ''
	lines = []
	i = 0
	for line in fa:
		i += 1
		line = line.rstrip()
		if line[0:1] == '>':
			if name != '':
				r[name] = ''.join(lines)
			name = line[1:]
			lines = []
		else:
			lines.append(line)
	r[name] = ''.join(lines)
	fa.close()
	return(r)

fa = readFasta(path2fa)

for l in sys.stdin:
	l = l.rstrip()
	line = l.split("\t")
	seq =  getSeq(fa,line[0],int(line[3]),int(line[4]),line[6])
	if (len(sys.argv) == 3) & (sys.argv[2] == 'fa'):
		print(">"+line[8]+"\n"+seq)
	else:
		print(line[8]+"\t"+seq)
