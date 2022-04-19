#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util;
use Array::Utils qw(:all);

#inputs/globals:
my $threads = 2;
my $eval = 0.001;
my $iters = 5;
my $clusterid = 0.9;
my ($toggleblast,$verbosity,$help) = 0;
my ($infasta,$database1,$database2,$database3);
my(@searched_seqs,@seqs_to_search,$firstmatch,$clustered,@newseedseqs,@seeds_not_searched);
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

if(-e "ice_blast.checkpoint"){
	my(@seq_to_blast,@seq_to_blast2,$readin1,$readin2,$readin3,$readin4,$readin5,$toggle1,$breakpoint,$toggle5) = 0;
	open(CHK, "< ice_blast.checkpoint");
	while(<CHK>){
		chomp;
		if($_=~/\#/){
			if($_=~/#1/){
				$readin1=1;
				$breakpoint++;
			}
			elsif($_=~/#2/){
				$readin2=1;
				$breakpoint++;
			}
			elsif($_=~/#3/){
				$readin3=1;
				$breakpoint++;
			}
			elsif($_=~/#4/){
				$readin4=1;
				$breakpoint++;
			}
			elsif($_=~/#5/){
				$readin5=1;
			}
			elsif($_=~/#6/){
				$breakpoint++;
			}
		}
		elsif($readin1 = 1){
			if($toggle1=0){
				(@seq_to_blast,@searched_seqs) = $_;
				$toggle1++;
			}
			else{
				my @seq_blasted = split(/\t/,$_);
				@seqs_to_search = array_diff(@seq_to_blast,@seq_blasted);
				$readin1 = 0;
			}
		}
		elsif($readin2 = 1){
			$firstmatch = $_;
			$readin2 = 0;
		}
		elsif($readin3 = 1){
			$clustered = $_;
			$readin3 = 0;
		}
		elsif($readin4 = 1){
			@searched_seqs = $_;
			$readin4 = 0;
		}
		elsif($readin5 = 1){
			if($toggle5=0){
				(@seq_to_blast2,@newseedseqs) = $_;
				$toggle5++;
			}
			else{
				my @seq_blasted2 = split(/\t/,$_);
				@seeds_not_searched = array_diff(@seq_to_blast2,@seq_blasted2);
				if(@seq_to_blast2 ~~ @seq_blasted2){
					$toggleblast = 1;
				}
				$readin5 = 0;
			}
		}
		else{
		}
		print "$_ \n This matches no known input marker for this checkpoint file.\nPlease report to the developer :) \n";
	}
	close CHK;
	open(CHK, "> ice_blast.checkpoint");
	if($breakpoint==0){
		goto MAIN;
	}
	elsif($breakpoint==1){
		goto FIRST_PSIBLAST;
	}
	elsif($breakpoint==2){
		goto FIRST_CLUSTERING;
	}
	elsif($breakpoint==3){
		goto FIRST_RESEED;
	}
	elsif($breakpoint==4){
		goto LOOP_START;
	}
	elsif($breakpoint==5){
		goto LOOP_END;
	}
	else{
		print "You somehow have a checkpoint file, but no expected breakpoint value...\nPlease report to the developer :)\n";
	}
}
else{}

MAIN:
open(CHK, "+> ice_blast.checkpoint");
if($verbosity==1){print "Removing unique characters in fasta file annotation lines.\n";}
RENAMEFASTA($infasta);
if($verbosity==1){print "Unique characters removed.\n\n";}
if($verbosity==1){print "Splitting multifasta input.\n";}
(@searched_seqs,@seqs_to_search) = SPLITFASTA($infasta,"");
if($verbosity==1){print "Splitting finished.\n\n";}
print CHK "#1\n@searched_seqs\n";

FIRST_PSIBLAST:
system("mkdir first_search");
if($verbosity==1){print "Starting PSI-BLAST search of queries.\n";}
$firstmatch = PSIBLAST($database1,"db1",$database3,"first_search",@seqs_to_search);
if($verbosity==1){print "Matches found. Moving all matches to archive.\n\n";}
system("mkdir archive");
system("mv $firstmatch archive");
RENAMEFASTA("archive/$firstmatch");
SPLITFASTA("archive/$firstmatch","archive");
print CHK "#2\n$firstmatch\n";

FIRST_CLUSTERING:
if($verbosity==1){print "Clustering matches from iteration.\n";}
$clustered = UCLUST("archive/$firstmatch",$clusterid);
if($verbosity==1){print "Clustering finished.\n\n";}
print CHK "#3\n$clustered\n";

FIRST_RESEED:
if($verbosity==1){print "Reseeding query sequences with centroid from clustering.\n";}
my @new_query_seqs = RE_SEED($clustered,@searched_seqs);
if($verbosity==1){print "Performing home directory cleanup.\n";}
opendir my $archive, "/archive";
my @old_matches = readdir $archive;
closedir $archive;
push(@searched_seqs,@new_query_seqs);
if($verbosity==1){print "First search and clustering complete. Moving to loop.\n\n";}
system("mkdir iteration_intermediates");
print CHK "#4\n@searched_seqs\n";
my ($loop,$toggleuclust,$togglereseed) = 0;
my $loopnum =1;

LOOP_START:
while($loop == 0){
	if($verbosity==1){print "Beginning new iteration.\n\n";}
	my $loop_extension="db2_".$loopnum;
	(@newseedseqs,@seeds_not_searched) = glob "*.fasta";
	print CHK "#5\n@newseedseqs\n";
	if($verbosity==1){print "Starting PSI-BLAST search of queries.\n";}
	if($toggleblast = 0){
		my $matches = PSIBLAST($database2,$loop_extension,$database3,"iteration_intermediates",@newseedseqs);
	}
	else{
		my $matches = PSIBLAST($database2,$loop_extension,$database3,"iteration_intermediates",@seeds_not_searched);
	}
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
print CHK "#6\n";

LOOP_END:
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
close CHK;

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
		print CHK "$seq\t";
	}
	my @blastoutputs = glob"*.$extension";
	open(BESTHITS, "+> besthits.txt");
	my %matches;
	foreach my $psiblastout (@blastoutputs){
		open(PSI, "< $psiblastout");
		while(<PSI>){
			chomp;
			if($matches{$_}){
				next;
			}
			else{
				print BESTHITS "$_\n";
				$matches{$_}=1;
			}
  	}
		close PSI;
	}
	close BESTHITS;
	system("mv *.$extension *.pssm *.fasta temp.txt $dir");
  system("blastdbcmd -db $extractdb -entry_batch besthits.txt -outfmt \"%f\" > all_matches.$database.$extension");
	system("mv besthits.txt $dir");
	print CHK "\n";
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
