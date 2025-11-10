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
    my($hen_std_fp)=WRITE_FASTA($hen_indata_std_ref, "FILE PATH");
    my($splice_std_fp)=WRITE_FASTA($splice_indata_std_ref, "FILE PATH");

    #run tBLASTn 
    my($hen_blast_fp)=TBLASTN($hen_std_fp, $nt_blastdb_fp);
    my($splice_blast_fp)=TBLASTN($splice_std_fp, $nt_blastdb_fp);

    #parse BLAST results
    my($hen_blast_data_ref)=PARSE_BLAST($hen_blast_fp);
    my($splice_blast_data_ref)=PARSE_BLAST($splice_blast_fp);

    #check for overlap
    FIND_OVERLAP($hen_blast_data_ref,$splice_blast_data_ref);

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

    my(@blast_data_in)=READIN($blastin_fp);
    foreach my $line (@blast_data_in){
        #fields are: 0->seqid, 1->match start 2->match end 3->strand
        my(@fields)= split /\t/, $line;
        
        #check if match is contained within any current matches
        #for positive strand matches
        if($fields[3]="plus"){
            foreach my $seqkey (keys %blast_data_out){
                if($blast_data_out{$seqkey}{"strand"} eq "positive"){
                    if($fields[1]>=$blast_data_out{$fields[0]}{"start"} &&
                       $fields[2]<=$blast_data_out{$fields[0]}{"end"}){
                        
                        #match is contained within an already existing match. Skip.
                        $skip_trigger=1;
                        last; 
                    }
                    else{ next; }
                }
                else{ next; }
            }                    
        }
        #for negative strand matches
        elsif($fields[3]="negative"){
            foreach my $seqkey (keys %blast_data_out){
                if($blast_data_out{$seqkey}{"strand"} eq "negative"){
                    if($fields[1]<=$blast_data_out{$fields[0]}{"start"} &&
                       $fields[2]>=$blast_data_out{$fields[0]}{"end"}){
                        
                        #match is contained within an already existing match. Skip.
                        $skip_trigger=1;
                        last; 
                    }
                    else{ next; }
                }
                else{ next; }
            }
        }

        #check if match was found already
        if($skip_trigger==1){
            $skip_trigger=0;
            next;
        }

        #input to hash
        $blast_data_out{$fields[0]}{"start"}=$fields[1];
        $blast_data_out{$fields[0]}{"end"}=$fields[2];
        $blast_data_out{$fields[0]}{"strand"}=$fields[3];
    }

    return(\%blast_data_out);
}

sub FIND_OVERLAP {
    #inputs - filtered match tables from HEN and splice searches
    #outputs - table of matches (using splice bounds) where HEN match was within splice match
    my(%hen_blastin)=%{my $ref1=shift};
    my(%splice_blastin)=%{my $ref2=shift};

    foreach my $hen_seq (keys %hen_blastin){
        #check for matches on either side
        foreach my $splice_seq (keys %splice_blastin){

            if($splice_blastin{$splice_seq}{"start"} <= $hen_blastin{$hen_seq}{"start"} 
                <= $splice_blastin{$splice_seq}{"end"}){

            }
            elsif($splice_blastin{$splice_seq}{"start"} >= $hen_blastin{$hen_seq}{"start"} 
                >= $splice_blastin{$splice_seq}{"end"}){

            }
            #needs to check if HEN start bound contained, or w/in 100 of
            #any splice match start OR end
            
            #then check if HEN end bound countained, or w/in 100 
            #of any splice match start or end
        
            #then ONLY IF both above are true, save the three matches and record them

            #mark whether hen is on sans or antisans strand relative to splice match 


        }
    }
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