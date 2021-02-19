# first liftover merged annotations by code/not.r/liftover.annotations.py
# then merge liftovered annotations:
# cat ../human.sajr first/*human.out* | grep -v '#' | cut -f1-5,7-8 | sort | uniq | perl -e '$i=1;while($l=<>){chomp $l;@l = split "\t",$l;print(join("\t",@l[0..4]),"\t.\t",join("\t",@l[5..6]),"\t",$i++,"\n");}' > test.human.allSpecies.sajr
# then run second round liftover

# seems it is old code, but I didn't find newer one
#!/usr/bin/python
import os
import sys

chainPath = '/home/mazin/chains/'
run = 'java -jar /home/mazin/bin/liftOver.jar '


original = ['chicken','mouse','rabbit','rat','human','macaque','opossum']

inx = int(sys.argv[1])
f = original.pop(inx / (len(original)-1))
t = original[inx % len(original)]

d = f + 'To' + t
command = run + "../../" + f + '.sajr' + ' ' + chainPath + d + ' ' + d + '.out ' + d + '.mult '+ d + '.unmap  -tryMapByEndsLength=10 -isZeroBased=false'
print command
os.system(command)

