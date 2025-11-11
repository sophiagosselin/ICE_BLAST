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
my $ovrlp_tolerance=100; #default
my $thd=1; #default
my $blast_eval="1E-20"; #default

#getoptions
GetOptions ('hen_in=s' => \$hen_in_fp, 'splice_in=s' => \$splice_in_fp, 
            'nt_blastdb_fp=s' => \$nt_blastdb_fp,'help+' => \$help, 'h+' => \$help,
            'thd=i' => \$thd, 'eval=s' => \$blast_eval, 'ovrlp=i' => \$ovrlp_tolerance);

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
    my($hen_indata_ref)=READIN_FASTA($hen_in_fp);
    my($splice_indata_ref)=READIN_FASTA($splice_in_fp);

    #standardize data
    my($hen_indata_std_ref)=STD_FASTA($hen_indata_ref);
    my($splice_indata_std_ref)=STD_FASTA($splice_indata_ref);

    #write std data to file
    my($hen_std_fp)=WRITE_FASTA($hen_indata_std_ref, "run_dir/hen_standardized.fna");
    my($splice_std_fp)=WRITE_FASTA($splice_indata_std_ref, "run_dir/splice_standardized.fna");

    #run tBLASTn 
    my($hen_blast_fp)=TBLASTN($hen_std_fp, $nt_blastdb_fp);
    my($splice_blast_fp)=TBLASTN($splice_std_fp, $nt_blastdb_fp);

    #parse BLAST results
    my($hen_blast_data_ref)=PARSE_BLAST($hen_blast_fp);
    my($splice_blast_data_ref)=PARSE_BLAST($splice_blast_fp);

    #check for overlap
    my($overlap_ref)=FIND_OVERLAP($hen_blast_data_ref,$splice_blast_data_ref);

    #print overlap data to table
    WRITE_OVERLAP($overlap_ref, "outputs/overlap_table.csv");

    #extract matches from overlap data
    EXTRACT_OVERLAP_REGION_SEQS($overlap_ref,$hen_blast_data_ref,$splice_blast_data_ref);

    #print fasta of contiguous regions
}

sub INITIALIZE {
    #put any one time startup stuff in here!
    
    #make directory structure:
    my @dirs=["run_dir","outputs"];
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
            $outdata{$key}=$line;
        }
    } 
    return(\%outdata);
}

sub READIN_LINES {
    #readin file, return data line by line in array
    my $path=shift;
    open(my $fh, "< $path") || die "Cannot open $path: $@";
    chomp (my @outdata = <$fh>);
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

sub TBLASTN {
    #inputs - fasta file, databse path
    #outputs - blast hits table

    my $fasta_fp=shift;
    my $db_fp=shift;
    my($blast_results_fp)=($db_fp=~/(.*?)\..*/)."\.blast"; #this might not work 

    #run tBLASTn
    system("tBLASTn -db $db_fp -query $fasta_fp -out $blast_results_fp 
            -inclusion_ethresh $blast_eval -outfmt \"6 sseqid sstart send sstrand\" 
            -thd $thd");

    return ($blast_results_fp)
}

sub PARSE_BLAST {
    #inputs - blast hit table file path
    #removes redundant overlapping matches
    #outputs - filtered and parsed blast data 

    #NOTE: Might need to add filter to remove new matches that contain old matches

    my $blastin_fp=shift;
    my %blast_data_out;
    my $skip_trigger=0;
    my $counter=0; #to ensure every match has a unique ID # in-case of multiple matches to same contig/genome

    my(@blast_data_in)=READIN($blastin_fp);
    foreach my $line (@blast_data_in){
        #fields are: 0->seqid, 1->match start 2->match end 3->strand
        my(@fields)= split /\t/, $line;
        
        #check if match is contained within any current matches
        #for positive strand matches
        if($fields[3] eq "plus"){
            foreach my $seqkey (keys %blast_data_out){
                if($blast_data_out{$seqkey}{"strand"} eq "positive"){
                    if($fields[1]>=$blast_data_out{$fields[0]}{"start"} &&
                       $fields[2]<=$blast_data_out{$fields[0]}{"end"}){
                        
                        #match is contained within an already existing match. Skip.
                        $skip_trigger=1;
                        last; 
                    }
                    else{
                        next;
                    }
                }
                else{
                    next;
                }
            }                    
        }
        #for negative strand matches
        elsif($fields[3] eq "negative"){
            foreach my $seqkey (keys %blast_data_out){
                if($blast_data_out{$seqkey}{"strand"} eq "negative"){
                    if($fields[1]<=$blast_data_out{$fields[0]}{"start"} &&
                       $fields[2]>=$blast_data_out{$fields[0]}{"end"}){
                        
                        #match is contained within an already existing match. Skip.
                        $skip_trigger=1;
                        last; 
                    }
                    else{
                        next;
                    }
                }
                else{
                    next;
                }
            }
        }

        else{
            PRINT_TO_SCREEN("ERROR");
        }

        #check if match was found already
        if($skip_trigger==1){
            $skip_trigger=0;
            next;
        }

        #input to hash
        my $unique_id="$fields[0]\_counter\_$counter";
        $blast_data_out{$fields[0]}{"start"}=$fields[1];
        $blast_data_out{$fields[0]}{"end"}=$fields[2];
        $blast_data_out{$fields[0]}{"strand"}=$fields[3];
        $counter++;
    }

    return(\%blast_data_out);
}

sub FIND_OVERLAP {
    #inputs - filtered match tables from HEN and splice searches
    #outputs - table of matches (using splice bounds) where HEN match was within splice match
    my(%hen_blastin)=%{my $ref1=shift};
    my(%splice_blastin)=%{my $ref2=shift};
    my %overlap_matches;

    foreach my $hen_seq (keys %hen_blastin){
        #check for matches on either side
        my $match_start = my $match_end = "NA";

        foreach my $splice_seq (keys %splice_blastin){

            #use the strand to determine which ends to add unto
            my($upperbound,$lowerbound);

            if($splice_blastin{$splice_seq}{"strand"} eq "positive"){
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
            
            elsif($splice_blastin{$splice_seq}{"strand"} eq "negtive"){
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

sub EXTRACT_OVERLAP_REGION_SEQS {
    #data structure is:
    #primary key -> HEN match
    #secondary key 1 "sb" -> splice match for start of HEN
    #secondary key 2 "eb" -> splice match for end of HEN
    my(%overlap_data)=%{my $ref1 = shift};

    #data structure is:
    #primary key -> HEN OR SPLICE match accession
    #secondary keys "start" and "end" -> start or end of match respectively
    #secondary keys "strand" -> positive or negative
    my(%hen_blastdata)=%{my $ref2 = shift};
    my(%splice_blastdata)=%{my $ref3 = shift};

    



}

sub PRINT_TO_SCREEN {
    #stuff
    my $message = shift;
    print "$message";
}
