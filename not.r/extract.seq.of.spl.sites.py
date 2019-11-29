#!/usr/bin/python
import os
import sys

#usage:cat *splicesites | python  .py  fasta species > .gtf
path2fa = sys.argv[1]
sp = sys.argv[2]

def getSeq(fa,seqname,start,stop):
	if seqname in fa:
		return(fa[seqname][start:stop])
	else:
		return('-')

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
juncs = {}

for l in sys.stdin:
	l = l.rstrip()
	if l in juncs:
		juncs[l] += 1
	else:
		juncs[l]  = 1

for j in juncs:
	line = j.split("\t")
	start = int(line[1])
	start = getSeq(fa,line[0],start+1,start+3)
	stop =  int(line[2])
	stop =  getSeq(fa,line[0],stop-2,stop)
	print(line[0]+"\t.\t.\t"+line[1]+"\t"+line[2]+"\t.\t"+line[3]+"\t.\t"+sp+":"+":".join(line)+":"+start+stop+":"+str(juncs[j]))



