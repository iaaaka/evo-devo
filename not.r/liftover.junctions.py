#!/usr/bin/python
import os
import sys

chainPath = os.environ['HOME']+'/projects/evo.devo/raw/chains/'
run = '/biosoftware/conda/bin/java -jar '+os.environ['HOME']+'/bin/liftOver.jar '


original = ['chicken','mouse','rabbit','rat','human','macaque','opossum']

inx = int(sys.argv[1]) - 1
f = original.pop(inx / (len(original)-1))
t = original[inx % len(original)]

d = f + 'To' + t
command = run + f + '.merged.gff' + ' ' + chainPath + d + ' ' + d + '.out ' + d + '.mult '+ d + '.unmap  -tryMapByEndsLength=10 -isZeroBased=true'
print command
os.system(command)

