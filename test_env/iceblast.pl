#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util;
use File::Copy;
no warnings 'experimental';

#DEV NOTES:
#add a DS output option such that you get two outputs:
#1. that is the domain of interest only
#2. that contains the whole database sequence
#3. check the psi-blast routine for redundant dup checking.
#9. add a description to each subroutine

#inputs/globals:
my $threads = 2;
my $eval = 1e-20;
my $iters = 4;
my $clusterid = 0.7;
my ($verbosity,$help,$loop,$domainspecific,$domain_length,$upperbound,$lowerbound,$iteration) = (0) x 8;
my($infasta,$outdatabase,$psidatabase,@placeholder);

#getoptions
GetOptions ('ds' =>\$domainspecific, 'v' =>\$verbosity, 'id=s' => \$clusterid, 't=s' => \$threads, 'e=s' => \$eval, 'i=s' => \$iters, 'in=s' => \$infasta, 'psidb=s' => \$psidatabase, 'outdb=s' => \$outdatabase, 'help+' => \$help, 'h+' => \$help);

#check for help call
if($help==1){
	die
"ICE-BLAST v1.0.2\n
Comprehensive protein database search for divergent homologs.
Gosselin S. Gogarten J.P. (In Preparation)
Iterative Cluster Expansion BLAST; a tool for comprehensive seuquence extraction.

Usage: perl iceblast.pl -in query -psidb clustered database -outdb database to search.

IMPORTANT: ICE-BLAST has a checkpointing system.
	If your run is interupted simply rerun your original command.

The following inputs are required:
[in]: Multifasta file you intend to use as seed sequence(s).
[psidb]: location of BLAST database intended for PSSM contruction.
	It is suggested that one use a clustered database (i.e. Uniref50).
[outdb]: location of BLAST database you want to pull all matches from.
	Make sure your outdb has been set up with the parse_seqids option.

The following inputs are optional. Defaults shown.

BLAST Options:
[t]: Number of threads for psiblast. D: 2
[i]: Number of psiblast iterations. D: 4
[e]: Evalue cutoff for inclusion. D: 1e-20

UCLUST Options:
	[id]: Percent ID for uclust clustering. D: 0.7

Mode Selection:
[ds]: Domain specific toggle.
	This mode will make ICE-BLAST only use the matched region of a subject sequence for future searches.
	This is useful when your target is a single domain, or a molecular parasite.
	Additionally this option filters out matches that are 20% larger or smaller than your query.
	By default this option is turned off.

Other Options:
[v]: Verbosity toggle. Pass 1 to activate. D: 0
";
}
else{}

VERBOSEPRINT(1, "Initializing");
#check for key inputs
if(!defined $infasta){
	die "No input detected for query sequences (-in)\n";
}
elsif(!defined $psidatabase){
	die "No input detected for clustered database (-psidb)\n";
}
elsif(!defined $outdatabase){
	die "No input detected for search database (-outdb)\n";
}

#check if blast database files exist
if(!-e "$outdatabase.pto"){
	die "BLAST database files to extract sequences from not detected (check your -outdb)\n";
}
if(!-e "$psidatabase.pto"){
	die "BLAST database files for PSSM construction not detected (check your -psidb)\n";
}

#check for previous run outputs, OR check query input
if(-d "to_run"){
	my($backup_ref,$torun_ref)=RECOVER();
	@placeholder = @{$backup_ref};
	my @toblast = @{$torun_ref};
	open(OUT, "+> infasta.fasta");
	foreach my $file (@toblast){
		open(IN, "< $file");
		while(<IN>){
			chomp;
			print OUT "$_\n";
		}
		close IN;
	}
	close OUT;
	$infasta = "infasta.fasta";
}
else{
	mkdir("to_run");
	#check if query input is a fasta file
	open(TEST, "< $infasta") or	die "Input for query sequence cannot be opened ($infasta)\n";
	while(<TEST>){
		if($_=~/\>/){
			last;
		}
		else{
			die "Query file does not seem to be a valid fasta file\n";
		}
	}
	close TEST;
}

#make directories if needed
unless(-d "intermediates"){
	mkdir("intermediates");
}
unless(-d "output"){
	mkdir("output");
}
unless(-d "archive"){
	mkdir("archive");
}

MAIN();

sub MAIN {
	if($domainspecific == 1){
		$domain_length = AVERAGE_QUERY_LENGTH($infasta);
		$upperbound = $domain_length*1.20;
		$lowerbound = $domain_length*.80;
	}
	CORELOOP($infasta,@placeholder);
	OUTPUT();
}

sub RECOVER{
	#subroutine that recovers a previous run based on the condition of
	#the home and to_run directories.
	my @to_run = glob "to_run/*.fasta";
	my @backup;
	open(BACKUP, "< backup.log");
	while(<BACKUP>){
		chomp;
		push(@backup,$_);
	}
	close BACKUP;
	return(\@backup,\@to_run);
}

sub BACKUP{
	#subroutine that prints to a file, and keeps track of which sequences
	#have been searched with, such that recovering a run is easier
	my $searcedseq = shift;
	open(BACKUP, ">> backup.log");
	print BACKUP "$searcedseq\n";
	close BACKUP;
}

sub CORELOOP{
	my $infile = shift;
	my @seqs_to_search = @_;
	STANDARDIZE($infile);
  my(@split_query) = SPLITFASTA($infile);
	push(@seqs_to_search,@split_query);
	my $matches = PSIBLAST(@split_query);
	my $clustered = UCLUST($matches);
  my(@new_query_asc) = RE_SEED($clustered,@seqs_to_search);
	if(!@new_query_asc){
		return();
	}
	else{
		my($new_infasta) = EXTRACTFASTA($matches,@new_query_asc);
		$iteration++;
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
	my @lengths;
	my $sequencel = 0;
	open(IN, "< $fastafile");
	while(<IN>){
		chomp;
		if($_=~/\>/){
			next if($sequencel == 0);
			push(@lengths,$sequencel);
			$sequencel = 0;
		}
		else{
			$sequencel += length($_);
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
	my($subjectL,%matches,$strand);
	open(BESTHITS, ">> besthits.txt");
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
		system("mv $infasta.blast6 $infasta.pssm $infasta intermediates");
		BACKUP($infasta);
	}
	close BESTHITS;
  system("blastdbcmd -db $outdatabase -entry_batch besthits.txt -outfmt \"%f\" > extracted_matches.$iteration");
	STANDARDIZE("extracted_matches.$iteration");
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
				print OUT "$_\n";
			}
		}
		close BCMD;
		close OUT;
		unlink "extracted_matches.$iteration";
		rename "temp.fasta", "extracted_matches.$iteration";
	}
	system("mv besthits.txt intermediates");
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
			my $pattern = "to_run/$split[8].fasta";
			if(grep( /^$pattern$/, @previously_searched ) ){
				next;
			}
			else{
				push(@cluster_seeds,$pattern);
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
	my ($toggle,%toggles);
	open(OUT, "+> new_seeds.fasta");
	while(<IN>){
	  if($_=~/\>/){
			$toggle = 0;
			next if($toggles{$_});
			my($fh)=($_=~/\>(.*?)\R/);
			if(grep( /^to_run\/$fh\.fasta$/,@seqstoget)){
				print OUT ">$fh\n";
				$toggles{$_}=1;
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
	move("$infasta","archive/$infasta");
	return("new_seeds.fasta");
}

sub VERBOSEPRINT{
	(my $verblevel, my $message) =@_;
	if($verblevel <= $verbosity){
		print "$message\n";
	}

}
