#!/usr/bin/perl
use strict;

if(@ARGV != 2){
	print("usage: fa gtf");
	die;
}

open IN,$ARGV[1] or die "cannot open $ARGV[1]: $!\n";
while(my $l=<IN>){
	next if(substr($l,0,1) eq '#');
	chomp $l;
	my @l = split('\t',$l);
	my @r = `samtools faidx $ARGV[0] $l[0]:$l[3]-$l[4]`;
	chomp @r;
	my $seq = join('',@r[1..$#r]);
	$seq = revCompl($seq) if($l[6] eq '-');
	print ">$l[8]\n".$seq."\n";
}

close IN;

sub revCompl{
	my $r = $_[0];
	$r =~ tr/ATGCatgc/TACGtacg/; 
	$r=reverse($r);
	$r
}
