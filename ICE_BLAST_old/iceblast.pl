#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util;
use File::Copy;
no warnings 'experimental';

#directories
mkdir("intermediates");
mkdir("output");
mkdir("archive");

#inputs/globals:
my $loopnum = 1;
my $threads = 2;
my $eval = 1e-20;
my $iters = 4;
my $clusterid = 0.7;
my ($verbosity,$help,$loop,$toggleuclust,$togglereseed,$matches,$domainspecific) = (0) x 7;
my($infasta,$database1,$database2,$database3,@searched_seqs,@seqs_to_search,$firstmatch,$clustered);

#getoptions
GetOptions ('ds' =>\$domainspecific, 'v' =>\$verbosity, 'id=s' => \$clusterid, 't=s' => \$threads, 'e=s' => \$eval, 'i=s' => \$iters, 'in=s' => \$infasta, 'db1=s' => \$database1, 'db2=s' => \$database2, 'db3=s' => \$database3, 'help+' => \$help, 'h+' => \$help);

#check for help call
if($help==1){
	die "\nThe following inputs are required: \n
  [in]: Multifasta file you intend to use as seed sequences\n
  [db1]: location of database clustered at lower similarity for intial pssm (i.e. uniref50)\n
  [db2]: location of database clustered at higher similarity for subsequent pssm (i.e. uniref90)\n
	[db3]: location of database you want to pull all matches from\n
	NOTE: Make sure database 3 has been set up with the parse_seqids option.\n\n

The following inputs are optional. Defaults shown.\n
  [t]: Number of threads for psiblast. D: 2\n
  [i]: Number of psiblast iterations. D: 4\n
  [e]: Evalue cutoff for inclusion. D: 1e-20\n
	[id]: Percent ID for uclust clustering. D: 0.7\n
	[v]: Verbosity toggle. Pass 1 to activate. D: 0\n
	[ds]: Domain specific toggle. You will only use the matched region of a subject sequence for future searches. D:0\n\n";
}
else{}

#main
MAIN();

sub MAIN {
  RENAMEFASTA($infasta);
  (@seqs_to_search) = SPLITFASTA($infasta);
  move($infasta,"intermediates/$infasta");
  if($verbosity==1){print "Starting PSI-BLAST search of queries.\n";}

  $firstmatch = PSIBLAST($database1,"db1",$database3,"intermediates",@seqs_to_search);
  if($verbosity==1){print "Matches found. Moving all matches to archive.\n\n";}

  if($verbosity==1){print "Clustering matches from iteration.\n";}
	$clustered = UCLUST($firstmatch);
  if($verbosity==1){print "Clustering finished.\n\n";}

  if($verbosity==1){print "Reseeding query sequences with centroid from clustering.\n";}
  my(@new_query_seqs) = RE_SEED($clustered,@seqs_to_search);

  if(@seqs_to_search ~~ @new_query_seqs){
  	die "No new centroid sequences found after first search.\nConsider using a less strict clustering ID cutoff, or a lower E-value cutoff.\n\n";
  }
  push(@searched_seqs,@new_query_seqs);
	my (@newseedseqs) = EXTRACTFASTA($firstmatch,@new_query_seqs);
	my $test = scalar(@newseedseqs);
	print "New Seqs = $test\n";
	move($firstmatch, "archive/$firstmatch");
	move($clustered, "intermediates/$clustered");
  if($verbosity==1){print "First search and clustering complete. Moving to loop.\n\n";}

  while($loop == 0){
  	if($verbosity==1){print "Beginning new iteration.\n\n";}
  	my $loop_extension="db2_".$loopnum;
  	if($verbosity==1){print "Starting PSI-BLAST search of queries.\n";}
  	$matches = PSIBLAST($database2,$loop_extension,$database3,"intermediates",@newseedseqs);
  	if((my $testfilesize = -s $matches) == 0){
  		print "No new matches found moving to output stage.\n";
  		last;
  	}
  	else{
  		$loopnum++;
			if($verbosity==1){print "Clustering matches from iteration.\n";}
			my $newcluster = UCLUST($matches);
			if($verbosity==1){print "Reseeding query sequences with centroid from clustering.\n";}
			my @new_q_seqs = RE_SEED($newcluster,@searched_seqs);
			(@newseedseqs) = EXTRACTFASTA($matches,@new_q_seqs);
  		move($newcluster,"intermediates/$newcluster");
			move($matches,"archive/$matches");
  		push(@searched_seqs,@new_q_seqs);
  		if($verbosity==1){print "Iteration complete. Starting new iteration of search.\n\n";}
  	}
  }
  LOOP_END:
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

sub RENAMEFASTA {
	#removes most unique characters from annotation lines
	#makes later searches and moving of files much easier.
	my $fastafile = shift;
	open(IN, "< $fastafile");
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
			open(OUT, "+> $fh.fasta");
			print OUT ">$fh\n";
			push(@fasta_files,"$fh.fasta");
	  }
	  else{
	    print OUT $_;
	  }
	}
	close OUT;
	return @fasta_files;
}

sub EXTRACTFASTA{
	#splits a multifasta file into several individual fasta files each containing
	#1 sequence. Use before blasting a file. Otherwise do not use.
	my $infasta = shift;
	my @seqstoget = @_;
	open(IN, "< $infasta") or die "Error: \n No valid multiple fasta input sequence passed to subroutine!";
	my @fasta_files;
	my $toggle = 0;
	while(<IN>){
	  if($_=~/\>/){
			close OUT;
			$toggle = 0;
			my($fh)=($_=~/\>(.*?)\R/);
			if(grep( /^$fh$/,@seqstoget)){
				open(OUT, "+> $fh.fasta");
				print OUT ">$fh\n";
				push(@fasta_files,"$fh.fasta");
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
	return(@fasta_files);
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
	my %matches;
	my $smallest = 50;
	open(BESTHITS, "+> besthits.txt");
	foreach my $seq (@query_files){
		if($domainspecific == 1){
			my $qlen = 0;
			open(INPUT, "< $seq");
			while(<INPUT>){
				if($_=~/\>/){
					next;
				}
				else{
						$qlen += length($_);
				}
			}
			close(INPUT);
			$smallest = $qlen*0.8;
		}
	  system("psiblast -db $database -out temp.txt -query $seq -out_pssm $seq.pssm -out_ascii_pssm $seq.asci.pssm -inclusion_ethresh $eval -outfmt \"6 sseqid sstart send\" -num_iterations $iters -num_threads $threads -save_pssm_after_last_round -max_target_seqs 50000");
		system("psiblast -db $extractdb -in_pssm $seq.pssm -out $seq.$extension -inclusion_ethresh $eval -evalue $eval -outfmt \"6 sseqid sstart send\" -num_threads $threads");
		if($verbosity==1){print "$seq search finished\n";}
		open(PSI, "< $seq.$extension") or die "No search conducted? \n";
		while(<PSI>){
			chomp;
			if($matches{$_}){
				next;
			}
			else{
				my @tabs = split(/\t/, $_);
				chomp(@tabs);
				if($domainspecific == 1){
					my $strand;
					if(($tabs[1]-$tabs[2]) >=0){
						my $test1 = $tabs[1]-$tabs[2];
						if($test1<=$smallest){
							next;
						}
						$strand="minus";
						print BESTHITS "$tabs[0]\ $tabs[2]\-$tabs[1]\ $strand\n";
					}
					else{
						my $test2 = $tabs[2]-$tabs[1];
						if($test2 <=$smallest){
							next;
						}
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
    system("mv $seq.$extension $seq.pssm $seq.asci.pssm $seq $dir");
	}
	close BESTHITS;
  system("blastdbcmd -db $extractdb -entry_batch besthits.txt -outfmt \"%f\" > all_matches.$database.$extension");
	RENAMEFASTA("all_matches.$database.$extension");
	if($domainspecific == 1){
		open(BCMD, "< all_matches.$database.$extension");
		open(OUT, "+> temp.fasta");
		while(<BCMD>){
			if($_=~/\>/){
				chomp;
				my($id,$prot)=($_=~/(\>.*)\_\d+\-\d+(\_+.*?)/);
				print OUT "$id"."$prot\n";
			}
			else{
				print OUT $_;
			}
		}
		close BCMD;
		close OUT;
		unlink "all_matches.$database.$extension";
		rename "temp.fasta", "all_matches.$database.$extension";
	}
	system("mv besthits.txt temp.txt $dir");
	return("all_matches.$database.$extension");
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
	return(@cluster_seeds);
}
