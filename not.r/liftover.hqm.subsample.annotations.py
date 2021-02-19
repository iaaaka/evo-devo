#!/usr/bin/python
import os
import sys

chainPath = os.environ['SGE_O_HOME'] + '/data/raw/evo.devo/chains/'
run = 'java -jar ' + os.environ['SGE_O_HOME'] +'/bin/liftOver.jar '


original = ['mouse','human','macaque']

inx = int(sys.argv[1]) -1
f = original.pop(inx / (len(original)-1))
t = original[inx % len(original)]

d = f + 'To' + t
#first
#command = run + "../../" + f + '.sajr' + ' ' + chainPath + d + ' ' + d + '.out ' + d + '.mult '+ d + '.unmap  -tryMapByEndsLength=10 -isZeroBased=false'
#second
command = run + f + '.allSpecies.sajr' + ' ' + chainPath + d + ' ' + d + '.out ' + d + '.mult '+ d + '.unmap  -tryMapByEndsLength=10 -isZeroBased=false'
print command
os.system(command)

