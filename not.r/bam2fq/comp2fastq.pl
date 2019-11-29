#!/usr/bin/perl
use strict;

my %a;
my %b;

open A,$ARGV[0] or die "cannot open $ARGV[0]: $!\n";
open B,$ARGV[1] or die "cannot open $ARGV[1]: $!\n";

while(my $l=<A>){
	$l = (split(' ',$l))[0];
	my $ll = $l."\n".<A>.<A>.<A>;
	$a{$ll}=1;
}

while(my $l=<B>){
	$l = (split(' ',$l))[0];
	my $ll = $l."\n".<B>.<B>.<B>;
	$b{$ll}=1;
}


print("a=".scalar(keys(%a))."\n"."b=".scalar(keys(%b))."\n");
my $c=0;

my %absB;
for my $al (keys %a){
	if($b{$al}){
		$c++;
		delete $b{$al};
	}else{
		$absB{$al} = 1;
	}
}
print("common=$c\n");

print("present only in a:");
my @ak = keys %absB;
for(my $i=0;$i<4;$i++){
	print($ak[$i],"\n");
}

print("present only in b:");
my @bk = keys %b;
for(my $i=0;$i<4;$i++){
	print($bk[$i],"\n");
}
close A;
close B;
