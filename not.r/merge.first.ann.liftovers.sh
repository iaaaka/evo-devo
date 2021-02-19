cd $SGE_O_HOME/projects/evo.devo/processed/annotation/hqmrboc.subsample/merged/liftover
for s in human macaque mouse rat rabbit opossum chicken
do
	cat ../$s.sajr first/*$s.out* | grep -v '#' | cut -f1-5,7-8 | sort | uniq | perl -e '$i=1;while($l=<>){chomp $l;@l = split "\t",$l;print(join("\t",@l[0..4]),"\t.\t",join("\t",@l[5..6]),"\t",$i++,"\n");}' > second/$s.allSpecies.sajr
done
