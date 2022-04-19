#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util;

my @previously_searched;
my $clustalout = $ARGV[0];
my @cluster_seeds;
open(CLUST, "< $clustalout");
while(<CLUST>){
	if($_=~/^#/){
		next;
	}
	elsif($_=~/^C\t.*/){
		my @split = split(/\t/,$_);
		if(grep( /^$split[8]$/, @previously_searched ) ){
			next;
		}
		else{
			print "$_\n";
			push(@cluster_seeds,"$split[8].fasta");
		}
	}
	else{
		next;
	}
}
close CLUST;
