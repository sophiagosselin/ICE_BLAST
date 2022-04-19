#!/usr/bin/perl -w
use strict;
use warnings;
use File::Copy;

my @clusterids=(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9);
my @evals=(1e-3,1e-5,1e-10,1e-20);

foreach my $evalue (@evals){
  foreach my $id (@clusterids){
    my $dir = "$evalue\_$id";
    mkdir($dir);
    copy("intein.fasta","$dir/intein.fasta");
    copy("ice.sh","$dir/ice.sh");
    copy("iceblast_v1.1.pl","$dir/iceblast_v1.1.pl");
  }
}
