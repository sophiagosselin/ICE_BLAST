#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util;
use File::Copy;

#globals and default values
my ($hen_in_fp, #HEN domain nucleotide sequences (fasta format)
    $splice_in_fp, #splicing domain nucleotide sequences (fasta format)
    $nt_blastdb_fp #nucleotide BLAST databse
    );
my $help=0;
my $ovrlp_tolerance=500; #default
my $thd=4; #default
my $blast_eval="1E-20"; #default

#for testing purposes
$hen_in_fp="hen.fna";
$splice_in_fp="splice.fna";
$nt_blastdb_fp="phagesdb/phagesdb.fna";

#getoptions
GetOptions ('hen=s' => \$hen_in_fp, 'splice=s' => \$splice_in_fp, 
            'blastdb=s' => \$nt_blastdb_fp,'help+' => \$help, 'h+' => \$help,
            'thd=i' => \$thd, 'eval=s' => \$blast_eval, 'ovrlp=i' => \$ovrlp_tolerance);

#maybe writes a help message *quelle surprise*
if($help==1){
    print 
    "If god exists then they are a spiteful, uncaring, and mercilless god. 
    They will not help you.";
}

MAIN();

sub MAIN {

    INITIALIZE("run_dir","outputs");

    #readin fasta inputs
    my($hen_indata_ref)=READIN_FASTA($hen_in_fp);
    my($splice_indata_ref)=READIN_FASTA($splice_in_fp);

    #standardize data
    my($hen_indata_std_ref)=STD_FASTA($hen_indata_ref);
    my($splice_indata_std_ref)=STD_FASTA($splice_indata_ref);

    #write std data to file
    my($hen_std_fp)=WRITE_FASTA($hen_indata_std_ref, "run_dir/hen_standardized.fna");
    my($splice_std_fp)=WRITE_FASTA($splice_indata_std_ref, "run_dir/splice_standardized.fna");

    #run tBLASTn 
    my($hen_blast_fp)=TBLASTX($hen_std_fp, $nt_blastdb_fp, "run_dir/hen_blast.txt");
    my($splice_blast_fp)=TBLASTX($splice_std_fp, $nt_blastdb_fp, "run_dir/splice_blast.txt");

    #parse BLAST results
    my($hen_blast_data_ref)=PARSE_BLAST($hen_blast_fp);
    my($splice_blast_data_ref)=PARSE_BLAST($splice_blast_fp);

    #check for overlap
    my($overlap_ref)=FIND_OVERLAP_MATCHES($hen_blast_data_ref,$splice_blast_data_ref);

    #print overlap data to table
    WRITE_OVERLAP($overlap_ref, "outputs/overlap_table.csv");

    #extract matches from overlap data
    EXTRACT_OVERLAP_SEQS($overlap_ref,$hen_blast_data_ref,$splice_blast_data_ref);

    #print fasta of contiguous regions
}

sub INITIALIZE {
    #put any one time startup stuff in here!
    #make directory structure:
    foreach my $dir (@_){
		unless(-d $dir){
			mkdir($dir);
		}
	}
}

sub READIN_FASTA {
    #inputs - fasta file path
    #outputs - hash with acession lines as keys, seqs as data
    my $path=shift;
    my %outdata;
    my $key;
    my(@indata)=READIN_LINES($path);
    #parse fasta data format
    foreach my $line (@indata){
        if($line=~/\>/){
            $key=$line;
        }
        else{
            $outdata{$key}.=$line;
        }
    } 
    return(\%outdata);
}

sub READIN_LINES {
    #readin file, return data line by line in array
    my $path=shift;
    open(my $fh, "< $path") || die "Cannot open $path: $@";
    my @outdata;
    while(<$fh>){
        chomp;
        push(@outdata,$_);
    }
    close $fh;
    return(@outdata);
}

sub WRITE_LINES {
    #inputs - array of lines to print (no EOL), and a path
    #outputs - prints lines to file.
    my $path=shift;
    my @lines=@_;
    open(my $fh, "+> $path") || die "Cannot open $path: $@";
    print $fh join ("\n", @lines); #this might not work as intented :)
    close $fh;
}

sub WRITE_FASTA {
    #inputs - hash of fasta formated data, a file path for print
    #outputs - fasta format file + file path

    my(%fasta_indata)=%{my $ref = shift};
    my $fp = shift;
    my @toprint;

    foreach my $key (keys %fasta_indata){
        push(@toprint,$key);
        push(@toprint,$fasta_indata{$key});
    }

    WRITE_LINES($fp,@toprint);

    return($fp);
}

sub STD_FASTA {
    #inputs - fasta data in hash
    #outputs - fasta data in hash w/o unique chars in annotation line
	my(%fasta_indata)=%{my $ref = shift};
	my %fasta_outdata;

    foreach my $key (keys %fasta_indata){
        my $keyholder=$key;
        $key=~s/[\ \[\]\(\)\:\;\/\.\-\~\`\!\@\#\$\%\^\&\*\=\+\{\}\?\'\"]/\_/g;
        $fasta_outdata{$key}=$fasta_indata{$keyholder};
    }

    return(\%fasta_outdata);
}

sub TBLASTX {
    #inputs - fasta file, databse path
    #outputs - blast hits table

    my $fasta_fp=shift;
    my $db_fp=shift;
    my $blast_results_fp = shift;

    #run tBLASTn
    system("tblastx -db $db_fp -query $fasta_fp -out $blast_results_fp 
            -evalue $blast_eval -outfmt \"6 sseqid sstart send sstrand\" 
            -num_threads $thd");

    return ($blast_results_fp)
}

sub PARSE_BLAST {
    #inputs - blast hit table file path
    #removes redundant overlapping matches
    #outputs - filtered and parsed blast data 

    my $blastin_fp=shift;
    my %blast_data_in;
    my $counter=0; #uinque ID for seqs
    
    my(@blast_data_in)=READIN_LINES($blastin_fp);
    foreach my $line (@blast_data_in){
        #fields are: 0->seqid, 1->match start 2->match end 3->strand
        my(@fields)= split /\t/, $line;
        my $unique_id="$fields[0]\_counter\_$counter";
        $blast_data_in{$unique_id}{"start"}=$fields[1];
        $blast_data_in{$unique_id}{"end"}=$fields[2];
        $blast_data_in{$unique_id}{"strand"}=$fields[3];
        $counter++;
    }

    my($dataref)=RECURSIVE_OVERLAP(\%blast_data_in);
    my(%blast_data_out)=%{$dataref};

    #print data to file for debugging
    my @lines;
    foreach my $seqkey (keys %blast_data_out){
        my $printline = "$seqkey\t$blast_data_out{$seqkey}{'start'}\t$blast_data_out{$seqkey}{'end'}\t$blast_data_out{$seqkey}{'strand'}";
        push(@lines,$printline);
    }
    WRITE_LINES("$blastin_fp.parsed",@lines);

    return(\%blast_data_out);
}

sub RECURSIVE_OVERLAP {
    #recursively look for overlap between matches until only non-overalpping matches remain
    #return final matches as hash
    my(%indata)=%{my $ref = shift};
    my %omatches;
    my %outdata;

    #start
    foreach my $keyq (keys %indata){
        #skip matches that have already been found to have overlap
        next if(exists $omatches{$keyq});
        #get generic key (without counter)
        my($gkeyq)=($keyq=~/(.*?)\_counter\_/);
        my @unsorted;
        
        foreach my $keyr (keys %indata){
            #check if attempting a self v self comparison
            next if($keyr eq $keyq);

            #skip matches that have already been found to have overlap
            next if(exists $omatches{$keyr});

            #check if generic key matches current genome/match
            my($gkeyr)=($keyr=~/(.*?)\_counter\_/);
            next if($gkeyr ne $gkeyq);

            #get a range for the new and old match being compared
            my @rrange;
            my $overlap_toggle=0;
            if(exists $outdata{$keyr}){
                @rrange=($outdata{$keyr}{"start"}..$outdata{$keyr}{"end"});
            }
            else{
                @rrange=($indata{$keyr}{"start"}..$indata{$keyr}{"end"});
            }
            @rrange = sort { $a <=> $b } @rrange;

            my @qrange=($indata{$keyq}{"start"}..$indata{$keyq}{"end"});
            @qrange = sort { $a <=> $b } @qrange; 
            
            #check if there is overlap between the ranges            
            foreach my $i (sort @qrange){
                foreach my $j (sort @rrange){
                    if($i==$j){
                        #debug
                        #print "I:$i J:$j\n";
                        $overlap_toggle=1;
                        $omatches{$keyr}=1;
                        last;
                    }
                    else{
                        next;
                    }
                }
                if($overlap_toggle==1){
                    push(@unsorted,shift(@qrange),pop(@qrange),shift(@rrange),pop(@rrange));
                    last;
                }
                else{
                    next;
                }
            }
            if($overlap_toggle==1){
                last;
            }
        }

        #if overlapping matches were found, get max bounds
        if(@unsorted){
            my @sorted = sort { $a <=> $b } @unsorted;
            my $start=shift(@sorted);
            my $end=pop(@sorted);

            #populate outdata with condensed match    
            $outdata{$keyq}{"start"}=$start;
            $outdata{$keyq}{"end"}=$end;
            $outdata{$keyq}{"strand"}=$indata{$keyq}{"strand"};
        }
        else{
            #need a mechanism to remove this key if you later find a match to it... probably needs to be in a separate hash?
            $outdata{$keyq}{"start"}=$indata{$keyq}{"start"};
            $outdata{$keyq}{"end"}=$indata{$keyq}{"end"};
            $outdata{$keyq}{"strand"}=$indata{$keyq}{"strand"};
        }
    }

    #debug
    my @lines;
    foreach my $seqkey (keys %outdata){
        my $printline = "$seqkey\t$outdata{$seqkey}{'start'}\t$outdata{$seqkey}{'end'}\t$outdata{$seqkey}{'strand'}";
        push(@lines,$printline);
    }
    WRITE_LINES("debug.txt",@lines);

    return(\%outdata);
}


sub FIND_OVERLAP_MATCHES {
    #inputs - filtered match tables from HEN and splice searches
    #outputs - table of matches (using splice bounds) where HEN match was within splice match

    my(%hen_blastin)=%{my $ref1=shift};
    my(%splice_blastin)=%{my $ref2=shift};
    my %overlap_matches;

    foreach my $hen_seq (keys %hen_blastin){
        #check for matches on either side
        my $match_start = my $match_end = "NA";

        foreach my $splice_seq (keys %splice_blastin){

            #skip if matches are from 2 different genomes!
            my($hen_nocounter)=($hen_seq=~/(.*?)\_counter\_.*/);
            my($splice_nocounter)=($splice_seq=~/(.*?)\_counter\_.*/);
            next if($hen_nocounter ne $splice_nocounter);

            #use the strand to determine which ends to add unto
            my($upperbound,$lowerbound);

            if($splice_blastin{$splice_seq}{"strand"} eq "plus"){
                $upperbound = ($splice_blastin{$splice_seq}{"start"})-$ovrlp_tolerance;
                $lowerbound = ($splice_blastin{$splice_seq}{"end"})+$ovrlp_tolerance;

                #then check if HEN is within the tolerance bounded range
                if($upperbound <= $hen_blastin{$hen_seq}{"start"} <= $lowerbound){
                    $match_start=$splice_seq; 
                }
                elsif($upperbound <= $hen_blastin{$hen_seq}{"end"} <= $lowerbound){
                    $match_end=$splice_seq;
                }
            }
            
            elsif($splice_blastin{$splice_seq}{"strand"} eq "minus"){
                $upperbound = ($splice_blastin{$splice_seq}{"start"})+$ovrlp_tolerance;
                $lowerbound = ($splice_blastin{$splice_seq}{"end"})-$ovrlp_tolerance;

                #then check if HEN is within the tolerance bounded range
                if($lowerbound <= $hen_blastin{$hen_seq}{"start"} <= $upperbound){
                    $match_start=$splice_seq;
                }
                elsif($lowerbound <= $hen_blastin{$hen_seq}{"end"} <= $upperbound){
                    $match_end=$splice_seq;
                }
            }
            else{
                #there was no match, skip to next!
                next;
            }
            #if a match was found above, record it!
            #need to check if a match was found to a start region and an end region
            if($match_end ne "NA" && $match_start ne "NA"){
                $overlap_matches{$hen_seq}{"sb"}=$match_start;
                $overlap_matches{$hen_seq}{"eb"}=$match_end;
                last;
            }
            else{
                next;
            }
        }

        #check if there is still an NA match to only 1 of the bounds
        #keep those matches
        if($match_end ne "NA" || $match_start ne "NA"){
            $overlap_matches{$hen_seq}{"sb"}=$match_start;
            $overlap_matches{$hen_seq}{"eb"}=$match_end;
        }
        else{
            next;
        }
    }
    return(\%overlap_matches);
}

sub WRITE_OVERLAP {
    my(%overlap_data)=%{my $hashref = shift};
    my $filepath = shift;
    my @toprint;

    foreach my $hen_key (keys %overlap_data){
        my $line = "$hen_key\t$overlap_data{$hen_key}{'sb'}\t$overlap_data{$hen_key}{'eb'}";
        push(@toprint,$line);
    }

    WRITE_LINES($filepath,@toprint);
}

sub EXTRACT_OVERLAP_SEQS {

    my(%overlap_data)=%{my $ref1 = shift};
    my(%hen_blastdata)=%{my $ref2 = shift};
    my(%splice_blastdata)=%{my $ref3 = shift};

    my %sequences_to_extract;

    foreach my $henkey (keys %overlap_data){
        my($endbound_start,$endbound_end,$startbound_start,$startbound_end);

        #get coordinates for start and end of splice domain matches
        my $splice_start = $overlap_data{$henkey}{"sb"};
        if($splice_start eq "NA"){
            $startbound_start = $hen_blastdata{$henkey}{"start"};
            $startbound_end= $hen_blastdata{$henkey}{"end"};
        }
        else{
            $startbound_start = $splice_blastdata{$splice_start}{"start"};
            $startbound_end = $splice_blastdata{$splice_start}{"end"};
        }
        my $splice_end = $overlap_data{$henkey}{"eb"};
        if($splice_end eq "NA"){
            $endbound_start = $hen_blastdata{$henkey}{"start"};
            $endbound_end= $hen_blastdata{$henkey}{"end"};
        }
        else{
            $endbound_start = $splice_blastdata{$splice_end}{"start"};
            $endbound_end = $splice_blastdata{$splice_end}{"end"};
        }

        #print "Henkey:$henkey\tStartbound_start: $startbound_start\tStartbound_end: $startbound_end\tEndbound_start: $endbound_start\tEndbound_end: $endbound_end\n";
        
        #now with your coordinates, find the bounds to extract
        my @unsorted;
        push(@unsorted, $startbound_end,$startbound_start,$endbound_end,$endbound_start);
        my @sorted = sort { $a <=> $b } @unsorted;
        my $seqstart=shift(@sorted);
        my $seqend=pop(@sorted);

        #print "Seqstart: $seqstart\tSeqend: $seqend\n";

        #assign to hash
        $sequences_to_extract{$henkey}{"start"}=$seqstart;
        $sequences_to_extract{$henkey}{"end"}=$seqend;
    }

    #now print matches to file
    my @lines;
    foreach my $seqs (keys %sequences_to_extract){
        my($match_no_counter)=($seqs=~/(.*?)\_counter\_.*/);
        my $line = "$match_no_counter\ $sequences_to_extract{$seqs}{'start'}\-$sequences_to_extract{$seqs}{'end'}\ $hen_blastdata{$seqs}{'strand'}";
        
        #debug
        #print "$match_no_counter\ $sequences_to_extract{$seqs}{'start'}\-$sequences_to_extract{$seqs}{'end'}\ $hen_blastdata{$seqs}{'strand'}\n";
        
        push(@lines,$line);
    }
    WRITE_LINES("seq_table_for_cmd.txt",@lines);

    #now send to be extracted using blast
    system("blastdbcmd -db $nt_blastdb_fp -entry_batch seq_table_for_cmd.txt -outfmt \"\%f\" > extracted_matches.fna");
}

sub PRINT_TO_SCREEN {
    #stuff
    my $message = shift;
    print "$message";
}
