#!/usr/bin/perl -w
use strict;
use warnings;

#NOTES:
#probably need to push intein sequence to an array instead of a variable
# will make evolving it easier.
# However!, needs a way to prevent stop codons from being added...
#  Maybe prevents a mutation if it would create a stop codon from the past two nucleotides?????


#create a set of random sequences with fake inteins
#size inputs are in AA, output in Nucleotides

my @codons=(TTT,TTC,TTA,TTG,TCT,TCC,TCA,TCG,TAT,TAC,TGT,TGC,TGG,CTT,CTC,CTA,CTG,CCT,CCC,CCA,CCG,CAT,CAC,CAA,CAG,CGT,CGC,CGA,CGG,ATT,ATC,ATA,ATG,ACT,ACC,ACA,ACG,AAT,AAC,AAA,AAG,AGT,AGC,AGA,AGG,GTT,GTC,GTA,GTG,GCT,GCC,GCA,GCG,GAT,GAC,GAA,GAG,GGT,GGC,GGA,GGG);
my @stop_codons=(TAA,TAG,TGA);

my @number_of_seqs = (1..$ARGV[0]);
my $number_of_inteins = $ARGV[1];
my @size_of_extein = (1..($ARGV[2]-1));
my $intein_size = $ARGV[3];
my @size_of_intein = (1..($intein_size));
my %random_seqs;
my %intein_sequences;

#generate hash of random extein sequences
for my $seq (@number_of_seqs){
  for my $position (@size_of_extein){
    $random_seqs{$seq}.= @codons[int(rand(61))];
  }
  $random_seqs{$seq}.=@stop_codons[int(rand(3))];
}

#generate the original "intein" sequence
my $ancestral_intein_sequence
for my $position (@size_of_intein){
  $ancestral_intein_sequence.=@codons[int(rand(61))];
}


my $og_intein = $ancestral_intein_sequence;
my $number_of_generations = 100000;
my $sub_rate = 0.003;
my @dna=(A,T,G,C);

#mutates the original "intein" sequence based on sub rate, and generations above.
#also prevents nonsense mutations from occuring
#outputs hash of evolved intein sequences
for(my $inteins = 0; $inteins < $number_of_inteins; $inteins++;){
  $ancestral_intein_sequence = $og_intein;
  for(my $gen = 0; $gen < $number_of_generations; $gen++;){
    my @ancestor = split(undef, $ancestral_intein_sequence);
    my $newinteinseq="";
    my $codon_counter = 0;
    for(my $i = 0; $i < $size_of_intein; $i++;){
      if(rand(1) < $sub_rate){
        my $new_nucleotide = @dna[int(rand 4)];
        if($codon_counter==0){
          my $test_codon = $new_nucleotide.@ancestor[($i+1)].@ancestor[($i+2)];
          if(grep( /^$test_codon$/,@stop_codons)){
            $newinteinseq.=$new_nucleotide;
          }
          else{
            $newinteinseq.=@ancestor[$i];
          }
        }
        elsif($codon_counter==1){
          my $test_codon = @ancestor[($i-1)].$new_nucleotide.@ancestor[($i+1)];
          if(grep( /^$test_codon$/,@stop_codons)){
            $newinteinseq.=$new_nucleotide;
          }
          else{
            $newinteinseq.=@ancestor[$i];
          }
        }
        elsif($codon_counter==2){
          my $test_codon = @ancestor[($i-2)].@ancestor[($i-1)].$new_nucleotide;
          if(grep( /^$test_codon$/,@stop_codons)){
            $newinteinseq.=$new_nucleotide;
          }
          else{
            $newinteinseq.=@ancestor[$i];
          }
        }
        else{}
      }
      else{
        #print original nucleotide
        $newinteinseq.=@ancestor[$i];
      }
      if($codon_counter == 2){
        $codon_counter=0;
      }
      else{
        $codon_counter++;
      }
    }
    $ancestral_intein_sequence=$newinteinseq;
  }
  $intein_sequences{$inteins}=$ancestral_intein_sequence;
}

my $insert_site = 200;
#insert the intein sequences into exteins at site determined above
foreach my $extein (values %random_seqs){
  my $exteinl = length($extein);
  m
}
