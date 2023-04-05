#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util;

#inputs/globals:
my $help = 0;
my $threads = 2;
my $eval = 0.001;
my $iters = 5;
my $clusterid = 0.9;
my $verbosity = 0;
my ($infasta,$database1,$database2,$database3);
GetOptions ('v' =>\$verbosity, 'id' => \$clusterid, 't' => \$threads, 'e' => \$eval, 'i' => \$iters, 'in=s' => \$infasta, 'db1=s' => \$database1, 'db2=s' => \$database2, 'db3=s' => \$database3, 'help+' => \$help, 'h+' => \$help);
if($help == 1){
	die "\nThe following inputs are required: \n
  [in]: Multifasta file you intend to use as seed sequences\n
  [db1]: location of database clustered at lower similarity for intial pssm (i.e. uniref50)\n
  [db2]: location of database clustered at higher similarity for subsequent pssm (i.e. uniref90)\n
	[db3]: location of database you want to pull all matches from\n
	NOTE: Make sure database 3 has been set up with the parse_seqids option.\n\n

The following inputs are optional. Defaults shown.\n
  [t]: number of threads for psiblast. D: 2\n
  [i]: number of psiblast iterations. D: 5\n
  [e]: evalue cutoff for inclusion. D: 0.001\n
	[id]: percent ID for uclust clustering. D: 0.9\n
	[v]: Verbosity toggle. Pass 1 to activate. D: 0\n\n";
}
else{}


#main code
#removes all whitespace and other special characters from acession lines and
#replaces the original file with the new one
if($verbosity==1){print "Removing unique characters in fasta file annotation lines.\n";}
RENAMEFASTA($infasta);
if($verbosity==1){print "Unique characters removed.\n\n";}

#takes an input sequence and splits it into an array of individual fasta files
#each has 1 sequence in it
if($verbosity==1){print "Splitting multifasta input.\n";}
my @searched_seqs = SPLITFASTA($infasta,"");
if($verbosity==1){print "Splitting finished.\n\n";}

#takes the array of fasta files and searches db1 to create a pssm for each one
#this pssm is then used to search the db of actual interest
#a single fasta file with all matches from the db of interest is the output
system("mkdir first_search");
if($verbosity==1){print "Starting PSI-BLAST search of queries.\n";}
my $firstmatch = PSIBLAST($database1,"db1",$database3,"first_search",@searched_seqs);
if($verbosity==1){print "Matches found. Moving all matches to archive.\n\n";}
system("mkdir archive");
system("mv $firstmatch archive");
RENAMEFASTA("archive/$firstmatch");

#nesecarry to make future searches easier.
SPLITFASTA("archive/$firstmatch","archive");

#clusters the output from PSIBLAST
if($verbosity==1){print "Clustering matches from iteration.\n";}
my $clustered = UCLUST("archive/$firstmatch",$clusterid);
if($verbosity==1){print "Clustering finished.\n\n";}

#takes clustered sequences and copys the seed sequence from each cluster out of
#the archive such that it can be easily used for the next search
if($verbosity==1){print "Reseeding query sequences with centroid from clustering.\n";}
my @new_query_seqs = RE_SEED($clustered,@searched_seqs);
if($verbosity==1){print "Performing home directory cleanup.\n";}

#checks the contents of the archive
opendir my $archive, "/archive";
my @old_matches = readdir $archive;
closedir $archive;
push(@searched_seqs,@new_query_seqs);

if($verbosity==1){print "First search and clustering complete. Moving to loop.\n\n";}
system("mkdir iteration_intermediates");
my $loop = 0;
#FFS WHY WON'T THIS INITIALIZE
my $loopnum =1;
#HAH YOU CAN NEVER ESCAPE!
while($loop == 0){
	if($verbosity==1){print "Beginning new iteration.\n\n";}
	my $loop_extension="db2_".$loopnum;
	my @newseedseqs = glob "*.fasta";
	if($verbosity==1){print "Starting PSI-BLAST search of queries.\n";}
	my $matches = PSIBLAST($database2,$loop_extension,$database3,"iteration_intermediates",@newseedseqs);
	if($verbosity==1){print "Matches found. Moving all matches to archive.\n\n";}
	$loopnum++;
	system("mv $matches archive");
	RENAMEFASTA("archive/$matches");
	SPLITFASTA("archive/$matches","archive");
	opendir my $archive, "/archive";
	my @new_matches = readdir $archive;
	closedir $archive;
	if(@new_matches ~~ @old_matches){
		last;
	}
	else{
		if($verbosity==1){print "Clustering matches from iteration.\n";}
		my $newcluster = UCLUST("archive/$matches",$clusterid);
		if($verbosity==1){print "Reseeding query sequences with centroid from clustering.\n";}
		my @new_q_seqs = RE_SEED($newcluster,@searched_seqs);
		@old_matches = @new_matches;
		push(@searched_seqs,@new_q_seqs);
		if($verbosity==1){print "Iteration complete. Starting new iteration of search.\n\n";}
	}
}

#print a single fasta output file with all unique sequences
system("mkdir output");
my %all_sequence_acs;
my @allmatches = glob"/archive/*.fasta";
open(FINAL, "+> output/all_matches.fasta");
foreach my $file (@allmatches){
	my $skip=0;
	open(IN, "< $file");
	while(<IN>){
		if($_=~/\>/){
			$skip=0;
			if(exists $all_sequence_acs{$_}){
				$skip=1;
				next;
			}
			else{
				$all_sequence_acs{$_}=0;
				print FINAL $_;
			}
		}
		else{
			if($skip==1){
				next;
			}
			else{
				print FINAL $_;
			}
		}
	}
	close IN;
}
close FINAL;

sub RENAMEFASTA {
	#removes most unique characters from annotation lines
	#makes later searches and moving of files much easier.
	my $fastafile = shift;
	open(IN, "< $fastafile");
	open(OUT, "+> temp.fasta");
	while(<IN>){
		if($_=~/\>/){
			$_=~s/[\ \[\]\(\)]/\_/g;
			$_=~s/\.//g;
			print OUT $_;
		}
		else{
			print OUT $_;
		}
	}
	close IN;
	close OUT;
	unlink $fastafile;
	rename "temp.fasta", $fastafile;
}

sub SPLITFASTA{
	#splits a multifasta file into several individual fasta files each containing
	#1 sequence
	my $infasta = shift;
	my $directory = shift;
	open(IN, "< $infasta") or die "Error: \n No valid multiple fasta input sequence passed to subroutine!";
	my @fasta_files;
	while(<IN>){
	  if($_=~/\>/){
			close OUT;
			my($fh)=($_=~/\>(.*?)\R/);
			open(OUT, "+> $directory/$fh.fasta");
	    print OUT ">$fh\n";
	    push(@fasta_files,"$directory/$fh.fasta");
	  }
	  else{
	    print OUT $_;
	  }
	}
	close OUT;
	return @fasta_files;
}

sub PSIBLAST{
	#takes an array of fasta files with only 1 sequence in each file
	#runs psiblast on each sequence, extracts matches from the db
	#returns a file with all matches after cleaning up the home directory
	#Inputs need to be: database, file extension, query files. IN THAT ORDER.
	my $database = shift;
	my $extension = shift;
	my $extractdb = shift;
	my $dir = shift;
	my @query_files = @_;
	foreach my $seq (@query_files){
	  system("psiblast -db $database -out temp.txt -query $seq -out_pssm $seq.pssm -out_ascii_pssm $seq.asci.pssm -inclusion_ethresh $eval -outfmt \"6 sseqid\" -num_iterations $iters -num_threads $threads -save_pssm_after_last_round -max_target_seqs 50000");
		system("psiblast -db $extractdb -in_pssm $seq.pssm -out $seq.$extension -inclusion_ethresh $eval -evalue $eval -outfmt \"6 sseqid\" -num_threads $threads");
		open(PSI, "< $seq.$extension");
		open(RANGE, "+> range_file.txt");
		my %matches;
		while(<PSI>){
			chomp;
			if($matches{$_}){
				next;
			}
			else{
				print RANGE "$_\n";
				$matches{$_}=1;
			}
  	}
  	close RANGE;
		close PSI;
  	system("blastdbcmd -db $extractdb -entry_batch range_file.txt -outfmt \"%f\" > all_matches.$database.$extension");
		system("mv $seq.$extension $seq.pssm $seq.asci.pssm range_file.txt $seq temp.txt $dir");
	}
	return("all_matches.$database.$extension");
}

sub UCLUST{
	#INPUTS
	my $allmatches = shift;
	my $id = shift;
	system("uclust --sort $allmatches --output $allmatches.sorted");
	system("uclust --input $allmatches.sorted --uc $allmatches.uc --id $id");
	return("$allmatches.uc");
}


sub RE_SEED{
	my $clustalout = shift;
	my @previously_searched = @_;
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
				push(@cluster_seeds,"$split[8].fasta");
			}
		}
		else{
			next;
		}
	}
	close CLUST;
	foreach my $seedseq (@cluster_seeds){
		system("cp archive/$seedseq $seedseq");
	}
	return(@cluster_seeds);
}
