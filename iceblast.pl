#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util;
use File::Copy;
no warnings 'experimental';

#inputs/globals:
my $threads = 1;
my $eval = 1e-10;
my $iters = 4;
my $clusterid = 0.7;
my $iterlimit="-1";
my $verbosity = my $help = my $loop = my $domainspecific = my $upperbound = my $lowerbound = my $iteration = 0;
my($input_file,$outdatabase,$psidatabase);

#getoptions
GetOptions ('i_l=s' => \$iterlimit, 'ds=s' =>\$domainspecific, 'v=i' =>\$verbosity, 'id=s' => \$clusterid, 't=s' => \$threads, 'e=s' => \$eval, 'psi_i=s' => \$iters, 'in=s' => \$input_file, 'psidb=s' => \$psidatabase, 'outdb=s' => \$outdatabase, 'help+' => \$help, 'h+' => \$help);

#check for help call
if($help==1){
	die
"ICE-BLAST v1.2.2\n
Comprehensive protein database search for divergent homologs.
Gosselin S. Gogarten J.P. (In Preparation)
Iterative Cluster Expansion BLAST; a tool for comprehensive seuquence extraction.

Usage: perl iceblast.pl -in query -psidb clustered database -outdb database to search.

IMPORTANT: ICE-BLAST has a checkpointing system.
	If your run is interupted simply rerun your original command.

The following inputs are required:
[in]: FASTA file you intend to use as seed sequence(s).
[psidb]: location of BLAST database intended for PSSM contruction.
	It is suggested that one use a clustered database (i.e. Uniref50).
[outdb]: location of BLAST database you want to pull all matches from.
	Make sure your outdb has been set up with the parse_seqids option.

The following inputs are optional. Defaults shown as (D: ###).

BLAST Options:
[t]: Number of threads for psiblast. D: 2
[psi_i]: Number of psiblast iterations. D: 4
[e]: Evalue cutoff for inclusion. D: 1e-20

UCLUST Options:
	[id]: Percent ID for uclust clustering. D: 0.7

Mode Selection:
[ds]: Domain specific toggle. D: Not activated
	This mode will make ICE-BLAST only use the matched region of a subject sequence for future searches.
	This is useful when your target is a single domain, or a molecular parasite.
	Use the flag and specify a cutoff % to filter out matches that are X% larger or smaller than your query.
	Example usage: -ds .25

Other Options:
[v]: Verbosity level. 1 for key checkpoints only. 2 for all messages. D: 0
[i_l]: Iteration limit. Specifically limits the number of iterations (psiBLAST + UCLUST) ICE-BLAST goes through. By default there is no limit.
";
}
else{}

VERBOSEPRINT(1, "Initializing");

#Setup, and seperate completed searches and uncompleted searches if this is a recovered run.
my($ex_queries_refference,$unex_queries_refference)=SETUP();
my @executed_queries = @{$ex_queries_refference};
my @input_queries = @{$unex_queries_refference};

MAIN();

sub MAIN {
	CORELOOP(\@input_queries,\@executed_queries);
	VERBOSEPRINT(1, "All searches completed. No new centroids found. Moving to output.\n");

	#cleanup
	my @seed_sequence_files = glob "new_seeds*";
	foreach my $seed_file (@seed_sequence_files){
		move("$seed_file","intermediates/");
	}
	unlink "temp.txt";

	OUTPUT();
}

sub SETUP{
	#Sets up all nessecary directories and variables.
	#Additionally checks if recovery is nessecary, and if so sends inputs to RECOVER subroutine
	my(@executed_queries,@unexecuted_queries);

	#check if recovery is needed
	if(-d "to_run"){
		my($ex_queries_refference_setup,$unex_queries_refference_setup)=RECOVER($domainspecific);
		@executed_queries = @{$ex_queries_refference_setup};
		@unexecuted_queries = @{$unex_queries_refference_setup};
		#make directories if needed
		DIRECTORY_CHECK("intermediates","output","archive","to_run");
	}

	#if no recovery is needed, prepare input file for searches and get domain length if DS mode is active
	else{
		#check file inputs for privelages and existance
		FILE_I_O_CHECK($input_file);
		FILE_I_O_CHECK($psidatabase);
		FILE_I_O_CHECK($outdatabase);

		#finds upper and lower bounds for BLAST result size filtering if requested
		#also checks if user input is between 0-1 if not convert appropriately
		if($domainspecific != 0){
			if($domainspecific > 1){
				$domainspecific = ($domainspecific/100);
			}
			DOMAIN_SPECIFIC_CUTOFF($input_file);
		}

		#make directories if needed
		DIRECTORY_CHECK("intermediates","output","archive","to_run");
		(@unexecuted_queries) = FASTA_PREP($input_file);
		#create backup file
		BACKUP("");
	}



	VERBOSEPRINT(1, "All inputs checked. Starting run.\n");
	return(\@executed_queries,\@unexecuted_queries);
}

sub CORELOOP{
	VERBOSEPRINT(1, "Beginning iteration $iteration.\n");

	#parse array inputs
	my($unex_queries_refference_core,$ex_queries_refference_core)= @_;
	my @executed_qs = @{$ex_queries_refference_core};
	my @unexecuted_qs = @{$unex_queries_refference_core};

	#takes all unsearched queries, returns all unique matches in the DB from BLAST search.
	VERBOSEPRINT(1, "Starting BLAST searches.\n");
	my $BLAST_matches = PSIBLAST_WORKFLOW(@unexecuted_qs);

	#exits loops if no new matches were found at all
	if(-z $BLAST_matches){
		return()
	}

	#updates array of searched sequences.
	push(@executed_qs,@unexecuted_qs);
	UNIQUE_ARRAY(@executed_qs);

	#clusters matches
	VERBOSEPRINT(1, "Clustering unique matches.\n");
	my $clustered_matches = UCLUST($BLAST_matches);

	#extracts cluster centroids to be used as new queries
	VERBOSEPRINT(1, "Identifying new query sequences.\n");
	my(@new_query_asc) = GET_NEW_QUERIES($clustered_matches,@executed_qs);

	#extracts new query sequences if possible
	my($new_infasta) = EXTRACT_FASTA($BLAST_matches,@new_query_asc);

	#exits loop if no new queries are found
	if(!@new_query_asc){
		return();
	}

	#begins the next iteration
	else{
		my (@new_unexecuted_queries) = FASTA_PREP($new_infasta);

		#checks if iteration limit has been reached
		$iteration++;
		if($iteration==$iterlimit){
			VERBOSEPRINT(1, "Iteration limit reached. Ending run.\n");
		}
		else{
			CORELOOP(\@new_unexecuted_queries,\@executed_qs);
		}
	}
	return();
}

sub OUTPUT{
	#prints out all unique sequences found during the search
	my %all_sequence_acs;
	my @allmatches = glob "archive/*";
	open(my $final, "+> output/all_matches.fasta");
	foreach my $file (@allmatches){
  		my $skip=0;
  		open(my $in, "< $file");
  		while(<$in>){
  			if($_=~/\>/){
  				$skip=0;
  				if(exists $all_sequence_acs{$_}){
  					$skip=1;
  					next;
	  			}
  				else{
  					$all_sequence_acs{$_}=0;
  					print $final "$_";
  				}
  			}
  			else{
  				if($skip==1){
  					next;
  				}
  				else{
  					print $final "$_";
  				}
  			}
  		}
  		close $in;
  	}
	close $final;
	VERBOSEPRINT(1, "ICE-BLAST run completed. Enjoy your sequences!\n");
}

sub PSIBLAST_WORKFLOW{
	#Takes a set of input sequences in FASTA format as input.
	#returns all unique matches
	my @query_files = @_;
	my $filtered_results;

	#uses each input file as a query for PSI-BLAST search then filters results for next step
	foreach my $infasta (@query_files){
		VERBOSEPRINT(2, "Creating PSSM for $infasta.\n");
	  	system("psiblast -db $psidatabase -query $infasta -out temp.txt -out_pssm $infasta.pssm -inclusion_ethresh $eval -outfmt \"6 sseqid sstart send\" -num_iterations $iters -num_threads $threads -save_pssm_after_last_round -max_target_seqs 50000");
		VERBOSEPRINT(2, "Conducting psiBLAST search with PSSM.\n");
		system("psiblast -db $outdatabase -in_pssm $infasta.pssm -out $infasta.blast6 -inclusion_ethresh $eval -evalue $eval -outfmt \"6 sseqid sstart send\" -num_threads $threads");
		VERBOSEPRINT(2, "Filtering matches.\n");
		$filtered_results = FILTER_BLAST("$infasta.blast6");
		move("$infasta.blast6","intermediates/");
		move("$infasta.pssm","intermediates/");
		move("$infasta","intermediates/");
		BACKUP("$infasta\n");
	}

	#extracts BLAST matches from the above searches
  	system("blastdbcmd -db $outdatabase -entry_batch $filtered_results -outfmt \"\%f\" > extracted_matches.$iteration");
	move("$filtered_results","intermediates/$filtered_results");
	if($domainspecific != 0){
		REMOVE_BOUNDS("extracted_matches.$iteration");
	}
	STANDARDIZE_FASTA("extracted_matches.$iteration");
	VERBOSEPRINT(1, "BLAST searches for iteration $iteration completed.\n");
	return("extracted_matches.$iteration");
}

sub FILTER_BLAST{
	#takes a blast output file (format 6) as input.
	#Retrieves the best hits for each match, and makes sure they pass domain specific cutoffs if appropriate.
	my $blast_results = shift;
	my (%best_hits,@entries_for_blastcmd);

	#READIN CURRENT MATCHES FOUND
	if(-e "filtered_matches_iteration_$iteration.txt"){
		open(my $in, "< filtered_matches_iteration_$iteration.txt");
		while(<$in>){
			chomp;
			if($domainspecific != 0){
				my @readin_columns=split(/\ /,$_);
				$best_hits{$readin_columns[0]}=1;
			}
			else{
				$best_hits{$_}=1;
			}
		}
		close $in;
	}
	else{}

	open(my $blast, "< $blast_results") or die VERBOSEPRINT(0, "Check your BLAST software and databases for issues. No BLAST output from search was found.\n");
	open(my $out, ">> filtered_matches_iteration_$iteration.txt");
	while(<$blast>){
		chomp;
		my @output_columns = split(/\t/,$_);
		next if(exists $best_hits{$output_columns[0]});
		if($domainspecific != 0){
			my ($strand,$site_start,$site_end,$match_length);
			if(($output_columns[1]-$output_columns[2]) >=0){
				$strand="minus";
				$site_start = $output_columns[2];
				$site_end = $output_columns[1];
			}
			else{
				$strand="plus";
				$site_start = $output_columns[1];
				$site_end = $output_columns[2];
			}
			$match_length = abs($site_end-$site_start);
			if($match_length<$lowerbound){
				next;
			}
			if($match_length>$upperbound){
				next;
			}
			print $out "$output_columns[0]\ $output_columns[1]\-$output_columns[2]\ $strand\n";
		}
		else{
			print $out "$output_columns[0]\n";
		}
		$best_hits{$output_columns[0]}=1;
	}
	close $blast;
	close $out;
	return("filtered_matches_iteration_$iteration.txt");
}

sub UCLUST{
	#INPUTS
	my $sequences_to_cluster = shift;
	system("usearch -cluster_fast $sequences_to_cluster -sort length -id $clusterid -uc $sequences_to_cluster.uc");
	if(!-e "$sequences_to_cluster.uc"){
		die VERBOSEPRINT(0,"UCLUST error. Check version number, and make sure it can be accessed from the command line via \"usearch\" without an alias.\n");
	}
	move("$sequences_to_cluster.sorted","intermediates/$sequences_to_cluster.sorted");
	return("$sequences_to_cluster.uc");
}


sub GET_NEW_QUERIES{
	#takes a clustal output file, and an array of searched sequences as Inputs
	#returns list of new inputs that have not been searched with yet.
	VERBOSEPRINT(1, "Reseeding for next iteration with centroid sequences.\n");
	my $clustalout = shift;
	my @previously_searched = @_;
	my @cluster_seeds;
	open(my $clust, "< $clustalout");
	while(<$clust>){
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
	close $clust;
	move("$clustalout","intermediates/$clustalout");
	return(@cluster_seeds);
}

sub EXTRACT_FASTA{
	#takes a fasta file, and an array of sequences as input
	#returns a fasta file with only the sequences in the array
	my $infasta = shift;
	my @seqstoget = @_;
	open(my $in, "< $infasta") or die VERBOSEPRINT(0, "No valid multiple fasta input sequence found $infasta. Please try running again, otherwise contact the developer.\n");
	my ($toggle,%toggles);
	open(my $out, "+> new_seeds_iteration_$iteration.fasta");
	while(<$in>){
	  if($_=~/\>/){
			$toggle = 0;
			next if($toggles{$_});
			my($fh)=($_=~/\>(.*?)\R/);
			if(grep( /^to_run\/$fh\.fasta$/,@seqstoget)){
				print $out ">$fh\n";
				$toggles{$_}=1;
				$toggle = 1;
			}
			else{
				next;
			}
	  }
	  else{
			if($toggle eq 1){
				print $out $_;
			}
			else{
				next;
			}
	  }
	}
	close $out;
	close $in;
	move("$infasta","archive/$infasta");
	return("new_seeds_iteration_$iteration.fasta");
}

sub VERBOSEPRINT{
	(my $verblevel, my $message) = @_;
	if($verblevel <= $verbosity){
		$| = 1;
		print "$message\n";
	}
}

sub FASTA_PREP{
	#takes a given multi fasta file as input
	#standardizes the file, then splits it into individual fasta files for each sequence
	my $fasta_file = shift;
	STANDARDIZE_FASTA($fasta_file);
	my(@fasta_files) = SPLIT_FASTA($fasta_file);
	return(@fasta_files);
}

sub FILE_I_O_CHECK{
	#checks a given file path for existance, and R/W privelages
	my ($path) = shift;
	if(!-e $path){
		VERBOSEPRINT(0,"$path does not exist. Check file name for errors.\nAlso make sure to run in a directory that has not previously completed a search.\n");
		die;
	}
	my($read_privelage,$write_privelage) = (-r $path, -w _);
	if(defined($read_privelage & $write_privelage)){
		VERBOSEPRINT(2,"All privelages present for $path.\n");
	}
	else{
		VERBOSEPRINT(0,"Missing privelages for $path. Check read and write privelages. Read $read_privelage, Write $write_privelage.\n");
		die;
	}
}

sub DIRECTORY_CHECK{
	#checks if directories exists. If not, the sub creates it.
	foreach my $directory (@_){
		unless(-d $directory){
			mkdir($directory);
		}
	}
}

sub RECOVER{
	#recovers array values for searched and unsearched queries
	VERBOSEPRINT(1, "Recovering from previous run.\n");
	my $domain_specific_toggle = shift;
	my @unexecuted_queries_backup = glob "to_run/*.fasta";
	my @executed_queries_backup;
	open(my $backup, "< backup.log");
	while(<$backup>){
		chomp;
		if($domain_specific_toggle != 0){
			my @split_recover = split(/\t/,$_);
			$upperbound = $split_recover[0];
			$lowerbound = $split_recover[1];
			$domain_specific_toggle=0;
		}
		else{
			push(@executed_queries_backup,$_);
		}
	}
	close $backup;
	return(\@executed_queries_backup,\@unexecuted_queries_backup);
}

sub BACKUP{
	#subroutine that prints to a file, and keeps track of which sequences
	#have been searched with, such that recovering a run is easier
	my $searcedseq_backup = shift;
	open(my $backup, ">> backup.log");
	print $backup "$searcedseq_backup";
	close $backup;
}

sub SPLIT_FASTA{
	#splits a multifasta file into several individual fasta files each containing
	#1 sequence. Returns list of these files.
	my $infasta_split = shift;
	my @fasta_files_split;
	my $handle_holder = 0;

	open(my $in, "< $infasta_split");
	while(<$in>){
	  if($_=~/\>/){
			chomp;
			if($handle_holder==0){}
			else{
				close $handle_holder;
			}
			my($fh)=($_=~/\>(.*)/);
			open(my $split_handle, "+> to_run/$fh.fasta") or die VERBOSEPRINT(0,"File $infasta_split at annotation line $_ was not properly formatted to create a new fasta file with. Check for unique characters, or windows end lines.\n");
			$handle_holder = $split_handle;
			print $handle_holder ">$fh\n";
			push(@fasta_files_split,"to_run/$fh.fasta");
	  }
	  else{
	    print $handle_holder $_;
	  }
	}
	close $handle_holder;
	close $in;
	move($infasta_split,"intermediates/$infasta_split");
	return @fasta_files_split;
}

sub STANDARDIZE_FASTA {
	#removes most unique characters from annotation lines
	#makes later searches and moving of files much easier.
	my $fastafile = shift;
	open(my $in, "< $fastafile");
	open(my $out, "+> temp.fasta");
	while(<$in>){
		if($_=~/\>/){
			$_=~s/[\ \[\]\(\)\:\;\/\.\-\~\`\!\@\#\$\%\^\&\*\=\+\{\}\?\'\"]/\_/g;
			print $out $_;
		}
		else{
			print $out $_;
		}
	}
	close $in;
	close $out;
	unlink $fastafile;
	rename "temp.fasta", $fastafile;
}

sub DOMAIN_SPECIFIC_CUTOFF{
	#takes a set of input sequences from one FASTA file
	#returns an upper and lower bound for sequence searches when using domain specific mode
	#uses the user specified % for + and -
	my $fastafile = shift;
	my @lengths;
	my $difference = my $sequence_length = 0;
	open(my $in, "< $fastafile");
	while(<$in>){
		chomp;
		if($_=~/\>/){
			next if($sequence_length == 0);
			push(@lengths,$sequence_length);
			$sequence_length = 0;
		}
		else{
			$sequence_length += length($_);
		}
	}
	close $in;
	push(@lengths,$sequence_length);
	my @sorted = (sort {$a <=> $b} @lengths);
	$upperbound = $sorted[-1];
	$lowerbound = $sorted[0];
	$upperbound += $lowerbound*$domainspecific;
	$lowerbound -= $lowerbound*$domainspecific;
	BACKUP("$upperbound\t$lowerbound\n");
}

sub UNIQUE_ARRAY{
	#takes an array as input, returns an array with only unique values.
	my %seen;
  	grep !$seen{$_}++, @_;
}

sub REMOVE_BOUNDS{
	#takes output from blastcmd and removes the bound data for the match
	#If this is not done, you can get 10 matches or more to the same sequence
	#with slightly different bounds.
	my $fasta_file = shift;
	open(my $in, "< $fasta_file");
	open(my $out, "+> temp.fasta");
	while(<$in>){
		chomp;
		if($_=~/\>/){
			my($new_asc_front,$new_asc_rear)=($_=~/(\>.*)\:\d+?\-\d+\ (.*)/);
			print $out "$new_asc_front\ $new_asc_rear\n";
			#print "$_\n";
		}
		else{
			print $out "$_\n";
		}
	}
	close $in;
	close $out;
	unlink $fasta_file;
	rename "temp.fasta", $fasta_file;
}
