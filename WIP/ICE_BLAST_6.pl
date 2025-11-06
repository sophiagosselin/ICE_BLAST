#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util;
use File::Copy;

#globals
my ($hen_input_file, $splice_input_file);
my $help=0;
#getoptions
GetOptions ('hen_in=s' => \$hen_input_file, 'splice_in=s' => \$splice_input_file, 
            'help+' => \$help, 'h+' => \$help);

if($help==1){
    print 
    "This is a help message. God cannot help you now.";
}


#Readin + write should be removed from subroutines and made separate methods.


MAIN();

sub MAIN {

    INITIALIZE();
    my($hen_input_dataref) = READIN_FASTA($hen_input_file);
    my($splice_input_dataref) = READIN_FASTA($splice_input_file);


}

sub MAKE_DATA_STRUCTURE {
    #checks if directories exists. If not, creates it.
	my @dirs=["run_dir"];
    foreach my $directory (){
		unless(-d $directory){
			mkdir($directory);
		}
	}
}

sub READIN_FASTA {
    my $path = shift;
    my %data;
    my $key;
    open(my $fh, "<" $path) || die "Cannot open $path: $@";
    while <$fh> {
        chomp;
        if($_=~/\>/){
            $key=$_;
        }
        else{
            $data{$key}=$_;
        }
    }
    close $fh;
    return(\%data);
}

sub WRITE_FASTA {
    my $path = shift;
    my 
}

sub TBLASTN {
    #inputs - fasta file, databse path
    #outputs - blast hits table
}

sub FILTER_RESULTS {
    #inputs - blast hit table
    #outputs - filtered matches

}

sub FIND_OVERLAP {
    #inputs - filtered match tables from HEN and splice searches
    #outputs - table of matches (using splice bounds) where HEN match was within splice match
}

sub EXTRACT_SEQ_MATCHES {
    #inputs - blast match table
    #outputs - fasta format sequences of matches
}

sub UCLUST {
    #inputs - fasta format file
    #outputs - uclust file paths
}

sub GET_CENTROIDS {
    #inputs - uclust file path, fasta format seq file
    #outputs - fasta format file with centroids only

}
