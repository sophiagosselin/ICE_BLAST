#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util;
use File::Copy;
no warnings 'experimental';

#directories
unless(-d "intermediates"){
	mkdir("intermediates");
}
unless(-d "output"){
	mkdir("output");
}
unless(-d "archive"){
	mkdir("archive");
}

#inputs/globals:
my $threads = 2;
my $eval = 1e-20;
my $iters = 4;
my $clusterid = 0.7;
my ($verbosity,$help,$loop,$domainspecific,$domain_length,$upperbound,$lowerbound,$iteration) = (0) x 7;
my($infasta,$outdatabase,$psidatabase,@placeholder);

#getoptions
GetOptions ('ds' =>\$domainspecific, 'v' =>\$verbosity, 'id=s' => \$clusterid, 't=s' => \$threads, 'e=s' => \$eval, 'i=s' => \$iters, 'in=s' => \$infasta, 'psidb=s' => \$psidatabase, 'outdb=s' => \$outdatabase, 'help+' => \$help, 'h+' => \$help);

#check for help call
if($help==1){
	die "\nThe following inputs are required: \n
  [in]: Multifasta file you intend to use as seed sequences\n
  [psidb]: location of database intended for PSSM contruction. Suggested that one use a clustered database (i.e. Uniref50) for speed concerns.\n
  [outdb]: location of database you want to pull all matches from\n
	NOTE: Make sure your outdb has been set up with the parse_seqids option.\n\n

The following inputs are optional. Defaults shown.\n
  [t]: Number of threads for psiblast. D: 2\n
  [i]: Number of psiblast iterations. D: 4\n
  [e]: Evalue cutoff for inclusion. D: 1e-20\n
	[id]: Percent ID for uclust clustering. D: 0.7\n
	[v]: Verbosity toggle. Pass 1 to activate. D: 0\n
	[ds]: Domain specific toggle. You will only use the matched region of a subject sequence for future searches. D:0\n\n";
}
else{}

#check_for_pre_run_seqs
if(-d "to_run"){
	RECOVER();
}
else{
	mkdir("to_run");
}

#main
MAIN();

sub MAIN {
	if($domainspecific == 1){
		$domain_length = AVERAGE_QUERY_LENGTH($infasta)
		$upperbound = $domain_length*1.20;
		$lowerbound = $domain_length*.80;
	}
	CORELOOP($infasta,@placeholder);
	OUTPUT();
}

sub RECOVER{
	#subroutine that recovers a previous run based on the condition of
	#the home and to_run directories.
	my @backup_seqs = glob "to_run/*.fasta";

}

sub BACKUP{
	#subroutine that prints to a file, and keeps track of which sequences
	#have been searched with, such that recovering a run is easier
}

sub CORELOOP{
	my $infile = shift;
	my @seqs_to_search = @_;
	STANDARDIZE($infile);
  my(@split_query) = SPLITFASTA($infile);
	push(@seqs_to_search,@split_query);
	my $matches = PSIBLAST(@split_query);
	$clustered = UCLUST($firstmatch);
  my(@new_query_asc) = RE_SEED($clustered,@seqs_to_search);
	if(!@new_query_asc){
		return();
	}
	else{
		my($new_infasta) = EXTRACTFASTA($matches,@new_query_asc);
		CORELOOP($new_infasta,@seqs_to_search);
	}
	return();
}

sub OUTPUT{
	#prints out all unique sequences found during the search 
	my %all_sequence_acs;
  my @allmatches = glob "archive/*";
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
  				print FINAL "$_";
  			}
  		}
  		else{
  			if($skip==1){
  				next;
  			}
  			else{
  				print FINAL "$_";
  			}
  		}
  	}
  	close IN;
  }
  close FINAL;
}

sub AVERAGE_QUERY_LENGTH{
	my $fastafile = shift;
	my $sequencel,@lengths;
	open(IN, "< $fastafile");
	while(<IN>){
		chomp;
		if($_=~/\>/){
			push(@lengths,$sequencel);
			$sequencel = 0;
		}
		else{
			$sequencel += length($_;)
		}
	}
	push(@lengths,$sequencel);
	my ($number,$average) = 0;
	foreach my $length (@lengths){
		$number++;
		$average += $length;
	}
	$average=$average/$number;
	return($average);
}

sub STANDARDIZE {
	#removes most unique characters from annotation lines
	#makes later searches and moving of files much easier.
	my $fastafile = shift;
	open(IN, "< $fastafile");a
	open(OUT, "+> temp.fasta");
	while(<IN>){
		if($_=~/\>/){
			$_=~s/[\ \[\]\(\)\:\;\/\.]/\_/g;
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
	#1 sequence. Use before blasting a file. Otherwise do not use.
	my $infasta = shift;
	open(IN, "< $infasta") or die "Error: \n No valid multiple fasta input sequence passed to subroutine!";
	my @fasta_files;
	while(<IN>){
	  if($_=~/\>/){
			close OUT;
			my($fh)=($_=~/\>(.*?)\R/);
			open(OUT, "+> to_run/$fh.fasta");
			print OUT ">$fh\n";
			push(@fasta_files,"to_run/$fh.fasta");
	  }
	  else{
	    print OUT $_;
	  }
	}
	close OUT;
	move($infasta,"intermediates/$infasta");
	return @fasta_files;
}

sub PSIBLAST{
	#THE OLD DESCRIPTION WAS NO LONGER ACCURATE LOL
	my @query_files = @_;
	my($subjectL,%matches);
	open(BESTHITS, "+> besthits.txt");
	foreach my $infasta (@query_files){
	  system("psiblast -db $psidatabase -out temp.txt -query $infasta -out_pssm $infasta.pssm -inclusion_ethresh $eval -outfmt \"6 sseqid sstart send\" -num_iterations $iters -num_threads $threads -save_pssm_after_last_round -max_target_seqs 50000");
		system("psiblast -db $outdatabase -in_pssm $infasta.pssm -out $infasta.blast6 -inclusion_ethresh $eval -evalue $eval -outfmt \"6 sseqid sstart send\" -num_threads $threads");
		open(PSI, "< $infasta.blast6") or die "No search conducted? \n";
		while(<PSI>){
			chomp;
			if($matches{$_}){
				next;
			}
			else{
				my @tabs = split(/\t/, $_);
				if($domainspecific == 1){
					if(($tabs[1]-$tabs[2]) >=0){
						$subjectL = $tabs[1]-$tabs[2];
						next if($subjectL<$lowerbound);
						next if ($subjectL>$upperbound);
						$strand="minus";
						print BESTHITS "$tabs[0]\ $tabs[2]\-$tabs[1]\ $strand\n";
					}
					else{
						$subjectL = $tabs[2]-$tabs[1];
						next if($subjectL<$lowerbound);
						next if ($subjectL>$upperbound);
						$strand="plus";
						print BESTHITS "$tabs[0]\ $tabs[1]\-$tabs[2]\ $strand\n";
					}
				}
				else{
					print BESTHITS "$tabs[0]\n";
				}
				$matches{$_}=1;
			}
  	}
		close PSI;
		system("mv $infasta.$extension $infasta.pssm $infasta intermediates");
	}
	close BESTHITS;
  system("blastdbcmd -db $outdatabase -entry_batch besthits.txt -outfmt \"%f\" > extracted_matches.$iteration");
	RENAMEFASTA("extracted_matches.$iteration");
	if($domainspecific == 1){
		open(BCMD, "< extracted_matches.$iteration");
		open(OUT, "+> temp.fasta");
		while(<BCMD>){
			chomp;
			if($_=~/\>/){
				my($id,$prot)=($_=~/(\>.*)\_\d+\-\d+(\_+.*?)/);
				print OUT "$id"."$prot\n";
			}
			else{
				print OUT $_;
			}
		}
		close BCMD;
		close OUT;
		unlink "extracted_matches.$iteration";
		rename "temp.fasta", "extracted_matches.$iteration";
	}
	system("mv besthits.txt temp.txt intermediates");
	return("extracted_matches.$iteration");
}

sub UCLUST{
	#INPUTS
	my $allmatches = shift;
	system("uclust --sort $allmatches --output $allmatches.sorted");
	system("uclust --input $allmatches.sorted --uc $allmatches.uc --id $clusterid");
	move("$allmatches.sorted","intermediates/$allmatches.sorted");
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
				push(@cluster_seeds,$split[8]);
			}
		}
		else{
			next;
		}
	}
	close CLUST;
	move("$clustalout","intermediates/$clustalout");
	return(@cluster_seeds);
}

sub EXTRACTFASTA{
	#
	my $infasta = shift;
	my @seqstoget = @_;
	open(IN, "< $infasta") or die "Error: \n No valid multiple fasta input sequence passed to subroutine!";
	my $toggle = 0;
	open(OUT, "+> new_seeds.fasta");
	while(<IN>){
	  if($_=~/\>/){
			$toggle = 0;
			my($fh)=($_=~/\>(.*?)\R/);
			if(grep( /^$fh$/,@seqstoget)){
				print OUT ">$fh\n";
				$toggle = 1;
			}
			else{
				next;
			}
	  }
	  else{
			if($toggle eq 1){
				print OUT $_;
			}
			else{
				next;
			}
	  }
	}
	close OUT;
	move("$infasta","intermediates/$infasta");
	return("new_seeds.fasta");
}
