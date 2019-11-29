#!/usr/bin/perl

use strict;

die("usage: $0 path.to.folder.with.logs") if(@ARGV != 1);

opendir my $dir, $ARGV[0] or die "Cannot open directory '$ARGV[0]': $!\n";
my @files = readdir $dir;
closedir $dir;

for my $f (@files){
	if($f =~ /(\d+b?sTS[^.]*)\..+\.log/){
		my $sid = $1;
		my @st = readLog("$ARGV[0]/$f");
		print($sid,"\t",join("\t",@st),"\n");
	}
}

sub readLog{
	open F,$_[0]  or die "Cannot open file '$_[0]': $!\n";
	my @res = (); #total,mapped uniq, mapped mult
	while(my $l = <F>){
		$res[0] = $1 if($l =~ /(\d+) reads; of these:/);
		$res[1] = $1 if($l =~ /(\d+) \(\d+.\d+%\) aligned exactly 1 time/);
		$res[2] = $1 if($l =~ /(\d+) \(\d+.\d+%\) aligned >1 times/);
	}
	close F;
	@res
}
