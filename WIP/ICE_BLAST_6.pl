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

#maybe writes a help message *quelle surprise*
if($help==1){
    print 
    "If god exists then they are a spiteful, uncaring, and mercilless god. 
    They will not help you.";
}

MAIN();

sub MAIN {

    INITIALIZE();

    #readin fasta inputs
    my($hen_indata_ref) = READIN_FASTA($hen_input_file);
    my($splice_indata_ref) = READIN_FASTA($splice_input_file);

    #standardize data
    my($hen_indata_std_ref) = STD_FASTA($hen_indata_ref);
    my($splice_indata_std_ref) = STD_FASTA($splice_indata_ref);


}

sub INITIALIZE {
    #put any one time startup stuff in here!
    
    #make directory structure:
    my @dirs=["run_dir"];
    foreach my $dir (@dirs){
		unless(-d $dir){
			mkdir($dir);
		}
	}

    #
}

sub READIN_FASTA {
    #inputs - fasta file path
    #outputs - hash with acession lines as keys, seqs as data
    my $path = shift;
    my %outdata;
    my $key;
    my(@indata)=READIN_LINES($path);
    #parse fasta data format
    foreach my $line (@indata){
        if($line=~/\>/){
            $key=$line;
        }
        else{
            $data{$key}=$line;
        }
    } 
    return(\%outdata);
}

sub READIN_LINES {
    #readin file, return data line by line in array
    my $path = shift;
    open(my $fh, "<" $path) || die "Cannot open $path: $@";
    chomp (my @outdata = <$fh>);
    close $fh;
    return(@outdata);
}

sub WRITE_LINES {
    #inputs - array of lines to print (no EOL), and a path
    #outputs - prints lines to file.
    my $path = shift;
    my @lines = @_;
    open(my $fh, "+>" $path) || die "Cannot open $path: $@";
    print $fh join ("\n", @lines); #this might not work as intented :)
    close $fh;
}

sub STD_FASTA {
    #inputs - fasta data in hash
    #outputs - fasta data in hash w/o unique chars in annotation line
	my(%fasta_indata) = %{my $ref = shift};
	my %fasta_outdata;

    foreach my $key (keys %fasta_indata){
        my $keyholder = $key;
        $key=~s/[\ \[\]\(\)\:\;\/\.\-\~\`\!\@\#\$\%\^\&\*\=\+\{\}\?\'\"]/\_/g;
        $fasta_outdata{$key}=$fasta_indata{$keyholder};
    }

    return(\%fasta_outdata);
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