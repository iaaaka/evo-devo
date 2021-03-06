import numpy as np
import re
import collections
import sys
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from intervaltree import Interval, IntervalTree
import cPickle as pickle
import gc

allSpecies=['human','macaque','mouse','rat','rabbit','opossum','chicken']
#allSpecies=['human','macaque','mouse']
#path2o='/home/mazin/evo.devo/annotation/all.species/merged/'
path='/uge_mnt/home/mazin/'
path2o=path+'data/processed/evo.devo/annotation/hqmrboc.subsample/merged/'
path2f=path2o+'liftover/first/'
path2s=path2o+'liftover/second/'

def readSegSites(fname):
	r = {}
	with open(fname) as f:
		for line in f:
			line = line.rstrip()
			l = line.split("\t")
			r[l[0]] = l[1]
	return r

def loadGTFwithSegID(fname,indexBy='segment_id',allowedSegIDs=None): 
	""" read GTF file into dict with key it according to 'indexBy'
	Args:
	fname (str) file name
	indexBy (str) string that should be used for indexing. Either segment_id, or coors, or attr (in this case whole last column will be used)
	Returns:
	dict af strings
	"""
	r = {}
	with open(fname) as f:
		for line in f:
			line = line.rstrip()
			if line[0:1] == '#':
				continue
			l = line.split("\t")
			if l[2] != 'segment' or (indexBy == 'segment_id' and not re.search('position=INTERNAL',l[8])): #to leave only internal in original
				continue
			coor = [l[i] for i in [0,6,3,4]]
			type = re.search('type=([^;]+);', l[8])
			if type is None:
				coor.append(None)
			else:
				coor.append(type.group(1))
			coor[2] = int(coor[2])
			coor[3] = int(coor[3])
			coor = tuple(coor)
			if indexBy == 'segment_id':
				l[8] = re.search('segment_id=([^;]+);',l[8]).group(1)
				if allowedSegIDs == None or l[8] in allowedSegIDs:
					r[l[8]] = coor
			elif indexBy == 'coors':
				r[coor[0:4]] = l[8]
			elif indexBy == 'attr':
				r[l[8]] = coor
			else:
				raise Exception('indexBy have unexpected value!')
	return r

def gc(h,k1,k2):
	"""
	Extract coordinates (four first element for two-dimensial dict
	:param h:
	:param k1:
	:param k2:
	:return: tuple of coordinates
	"""
	r = h.get(k1,None).get(k2,None)
	if r is not None:
		r = r[0:4]
	return r

def getOrthMatrix(spec,sid,species):
	mergedInxs = {spec : mergs[spec][gc(origs,spec,sid)]}
	#get id from merged files
	for s in [x for x in species if x != spec]:
		mergedInxs[s] = mergs[s].get(gc(frsts,spec+'.'+s,sid),None)
	om = {s:{s:None for s in species} for s in species}
	om[spec][spec] = gc(origs,spec,sid)
	for s in species:
		for ss in [x for x in species if x != s]:
			om[s][ss] = gc(scnds,s+'.'+ss,mergedInxs[s])
	return om

def getOrthSMatricesForSpecies(sps):
	orthSegs = {}
	for s in sps:
		print "\t",s
		for sid in origs[s].keys():
			orthSegs[sid] = getOrthMatrix(s,sid,sps)
	return orthSegs

def getIntervalOverlap(a,b):
	if a is None or b is None or a[0] != b[0] or a[1] != b[1]:
		return 0
	return max(0,float(min(a[3],b[3]) - max(a[2],b[2])+1)/(max(a[3],b[3]) - min(a[2],b[2]) +1))

def checkOrthMatrixes(oms,sps):
	"""
	Check that all coordinates in given species are the same
	:param oms: list of orth matrices
	:param sps: list of species
	:return: minimum overlap of coordinates, number of species lost, number of edges lost, and coordinates in each species
	"""
	r = {}
	for k in oms.keys():
		o = oms[k]
		overlap = 1
		lostSp = 0
		lostEdge = -(len(sps)-1)
		unions = {}
		for s2 in sps:
			spEx = False
			unionCoor = overlapCoor = None
			for s1 in sps:
				lostEdge += o[s1][s2] is None
				spEx = spEx or o[s1][s2] is not None
				if o[s1][s2] is not None:
					if unionCoor is None:
						unionCoor = list(o[s1][s2])
						overlapCoor = list(o[s1][s2])
					else:
						if unionCoor[0] != o[s1][s2][0] or unionCoor[1] != o[s1][s2][1] or overlapCoor is None:
							overlap = 0
							overlapCoor = None
						else:
							unionCoor[2] = min(unionCoor[2], o[s1][s2][2])
							unionCoor[3] = max(unionCoor[3], o[s1][s2][3])
							overlapCoor[2] = max(overlapCoor[2], o[s1][s2][2])
							overlapCoor[3] = min(overlapCoor[3], o[s1][s2][3])
			overlap = min(overlap,getIntervalOverlap(unionCoor,overlapCoor))
			unions[s2] = overlapCoor
			if not spEx:
				lostSp += 1
		r[k] = (overlap,lostSp,lostEdge,unions)
	return r

def getCorrectOrthSegs(oms,species):
	ch = checkOrthMatrixes(oms,species)
	r = {}
	for k in ch.keys():
		if ch[k][0] and ch[k][1] == 0 and ch[k][2] == 0:
			t = oms[k]
			os = []
			for s2 in species:
				for s1 in species:
					if t[s1][s2] != None:
						os.append(t[s1][s2])
						break
			r[k] = tuple(os)
	return r

def getOrthSegs(check,minOverlap,sps):
	raise NameError('Deprecated: use getOrthSegsByOverlap instead')
	gos = {k:check[k][3] for k in check.keys() if check[k][0] >= minOverlap and check[k][1] == 0 and check[k][2] == 0}
	orthIds = {}
	for k in gos.keys():
		c = tuple([tuple(gos[k][s]) for s in sps])
		if c in orthIds:
			orthIds[c].append(k)
		else:
			orthIds[c] = [k]
	return orthIds

def p(m):
	""" function to print orth Matrix """
	l = 0
	x = {k:dict(m[k]) for k in m}
	sps = x.keys()
	sps.sort()
	for s in sps:
		for ss in sps:
			if x[s][ss] is None:
				x[s][ss] = 'None'
			else:
				x[s][ss] = ':'.join([str(i) for i in x[s][ss]])
			l = max(l,len(x[s][ss]))
	for s in sps:
		sys.stdout.write(s + (' ' * (l - len(s)) + " "))
	print
	for s in sps:
		for ss in sps:
			sys.stdout.write( x[s][ss]+(' '*(l-len(x[s][ss])))+" ")
		print


def getOrthSegsByOverlap(check, minOverlap, allLiftovers = True):
	sp =  check.values()[1][3].keys()[0]
	check = {k: check[k] for k in check.keys() if check[k][0] > minOverlap and check[k][1] == 0 and (check[k][2] == 0 or not allLiftovers)}
	segsOverlap = {}
	# create interval trees
	segs = {}
	segIds = check.keys()
	segInx = {segIds[i]: i for i in range(0, len(segIds))}
	for id in segIds:
		s = check[id][3][sp]
		chrKey = tuple(s[0:2])
		if chrKey not in segs:
			segs[chrKey] = IntervalTree()
		segs[chrKey].addi(s[2], s[3] + 1, id)  # IntervalTree is not end-inclusive
	gFrom = []
	gTo = []
	for i in range(0, len(segIds), 1):
		seg1 = check[segIds[i]][3]
		segsOverlap[i] = check[segIds[i]][0]
		s = seg1[sp]
		overlapS = segs[tuple(s[0:2])][s[2]:(s[3] + 1)]
		# check that overlap exists in all species
		for os in overlapS:
			if os[2] == segIds[i]:
				continue
			seg2 = check[os[2]][3]
			bySpecOverlap = 1
			for sps in seg1.keys():
				bySpecOverlap = min(bySpecOverlap,getIntervalOverlap(seg1[sps], seg2[sps]))
			if bySpecOverlap > minOverlap:
				segsOverlap[i] = min(bySpecOverlap,segsOverlap[i],check[os[2]][0])
				gFrom.append(i)
				gTo.append(segInx[os[2]])
	# look for linked components
	graph = csr_matrix(([1] * len(gFrom), (gFrom, gTo)), (len(segIds), len(segIds)))
	ccomp = connected_components(graph, directed=False)[1]
	res = {}
	for i in range(0, len(ccomp)):
		if ccomp[i] not in res:
			res[ccomp[i]] = [1,[]]
		res[ccomp[i]][0] = min(res[ccomp[i]][0],segsOverlap[i])
		res[ccomp[i]][1].append(segIds[i])
	res = res.values()
	for r in res:
		unionCoor = {}
		u = check[r[1][0]][3]
		for sid in r[1]:
			for s in u.keys():
				unionCoor[s[0:3]] = [s[0:3],u[s][0],u[s][1],min(u[s][2],check[sid][3][s][2]),max(u[s][3],check[sid][3][s][3])]
		unionCoor = {k[0:3]:unionCoor[k] for k in unionCoor.keys()}
		for k in (set(unionCoor.keys()) - set([k[0:3] for k in r[1]])):
			r[1].append(unionCoor[k])
	return res

def checkSelfOverlap(os,minOverlap):
	"""takes as input an array of arrays of ids of segments. Segments in each array are assumed to be orthologous
	return indexes of arrays that includes more than one segments from the same species"""
	return [o for o in range(1,len(os)) if os[o][0] > minOverlap and len(os[o][1]) != len(set([i[0:3] for i in os[o][1]]))]

def save_obj(obj, name ):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)
# code

#load files
origs = {}
mergs = {}
frsts = {}
scnds = {}

ad = readSegSites(path2o+'seg2sites.tab')
ad = {k:1 for k in ad.keys() if ad[k] == 'ad'}

for s in allSpecies:
	print s
	origs[s] = loadGTFwithSegID(path2o+s+'.sajr',allowedSegIDs=ad)
	mergs[s] = loadGTFwithSegID(path2s+s+'.allSpecies.sajr','coors')
	for ss in [x for x in allSpecies if x != s]:
		print "\t"+ss
		frsts[s+'.'+ss] = loadGTFwithSegID(path2f+s+'To'+ss+'.out')
		scnds[s+'.'+ss] = loadGTFwithSegID(path2s+s+'To'+ss+'.out','attr')


del(ad)

save_obj(origs,path2o+'liftover/origs')
save_obj(mergs,path2o+'liftover/mergs')
save_obj(frsts,path2o+'liftover/frsts')
save_obj(scnds,path2o+'liftover/scnds')


#first, original annotations were liftovered to all other species
#second, all exons in given species (original, and liftovered from other species) were collapsed and new IDs were introduced
#third, all coordinates obtained on previous step were liftovered again to all species

species = {'r':'rat','m':'mouse','b':'rabbit','h':'human','q':'macaque','o':'opossum','c':'chicken'}

origs=load_obj(path2o+'liftover/origs') #hash of hashes of segment coordinates [species][seg_id] from original annotations
#after first liftover union of all exons was calculated for each species (so it includes all exons of the species and all liftovered one, duplicates were collapsed and new ID were introduced (last column in gtf)
mergs=load_obj(path2o+'liftover/mergs') #hash of hashes of segment_id (new IDs that just identify unique coordinates in given species)
frsts=load_obj(path2o+'liftover/frsts') #[species.species][segment_id] - coordinates of exon from first species in the second one. identified by original segment ids
scnds=load_obj(path2o+'liftover/scnds') #same as frsts but for second liftover and new segment IDs are used

# test
# minOverlap = 0.0
# allLiftovers = False
# o = [species[i] for i in 'hqmrboc']
# t1 = getOrthSMatricesForSpecies(o)
# save_obj(t1,path2o+'liftover/hqmrboc.orth.table')
# p(t1['hum_034353.s8'])
# t2 = checkOrthMatrixes(t1, o) 
# t3 = getOrthSegsByOverlap(t2, minOverlap,allLiftovers=allLiftovers) 


#make orth segments

#names = ['hqmrboc','hqmrbo','hqmrb','hmo','hqmo','hmoc','hqmoc','hm','qm','rm','bm','om','cm','hq','hr','hb','ho','hc']
names = ['hqmrboc','hqm']
# make orth tables
minOverlap = 0.0
allLiftovers = False
for minOverlap in [0.0,0.6]:
	for ns in names:
		print(ns)
		o = [species[i] for i in ns]
		t = getOrthSMatricesForSpecies(o) #it is dictionary by seg id. For each original segment it contains coordinates in all species obtained by all possible ways.
		t = checkOrthMatrixes(t, o) #it is collapsed version of 'getOrthSMatricesForSpecies'. There is a single interval (union) for all coordinates for each species. Plus include  minimum overlap of coordinates, number of species lost, number of edges lost	
		t = getOrthSegsByOverlap(t, minOverlap,allLiftovers=allLiftovers) #ids for segments with overlapped liftovered coordinates
		print(collections.Counter([len(x[1]) for x in t]))
		t = [x for x in t if len(x[1]) == len(o)] #filter out all ambiguous cases
		#if len([i for i in t if len(i[1])>len(o)]) > 0: #looks for cases where "orth segments" contains two segments from same species
		#raise NameError('Self-overlapping segments!!')
		f = open(path+'/projects/evo.devo/processed/annotation/hqmrboc.subsample/merged/orth.seg/' + ns + "." + str(minOverlap) + ".orth.ad.segs", 'w')
		f.write('#minOverlap = '+str(minOverlap)+"\n#allLiftovers = "+str(allLiftovers)+"\n")
		f.write('#all cases that include more than one segments from same species were filtered out\n')
		for s in t:
			f.write(str(s[0])+"\t")
			v = s[1]
			for i in range(0,len(v)):
				if isinstance(v[i],list):
					v[i] = '_'.join([str(x) for x in v[i]])
			v = {i[0:3]:i for i in v}
			v = [str(v.get(s[0:3])) for s in o]
			f.write('\t'.join(v) + "\n")
		f.close()





#pdb.set_trace()
#collections.Counter(types)
#for s1 in z:
#	for s2 in z:
#		sys.stdout.write(str(round(min([getIntervalOverlap(t2[s1][3][s],t2[s2][3][s]) for s in o]),2))+'\t')
#	print()


#with open('/home/mazin/c.pkl','wb') as f: pickle.dump(check, f, -1)
#with open('/home/mazin/iitp.disk/Solexa/ma.disk/mazin/c.pkl', 'rb') as f:
#	check = pickle.load(f)



