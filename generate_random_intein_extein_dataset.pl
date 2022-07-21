#!/usr/bin/perl -w
use strict;
use warnings;

#NOTES:
#need to add among site rate variation for the intein evolution...
#Therefore: need a distribution (gamma?), such that:
  #-the c-terminal is under purifying selection
  #-the n-termial is under purifying selection
  #-the HE is under positive selection


#create a set of random sequences with fake inteins
#size inputs are in AA, output in Nucleotides

#dictionaries
my @codons=(TTT,TTC,TTA,TTG,TCT,TCC,TCA,TCG,TAT,TAC,TGT,TGC,TGG,CTT,CTC,CTA,CTG,CCT,CCC,CCA,CCG,CAT,CAC,CAA,CAG,CGT,CGC,CGA,CGG,ATT,ATC,ATA,ATG,ACT,ACC,ACA,ACG,AAT,AAC,AAA,AAG,AGT,AGC,AGA,AGG,GTT,GTC,GTA,GTG,GCT,GCC,GCA,GCG,GAT,GAC,GAA,GAG,GGT,GGC,GGA,GGG);
my @stop_codons=(TAA,TAG,TGA);
my @dna=(A,T,G,C);

#gamma function
sub gamma{
  #take a nucleotide position and decide what it's mutation rate is per gen.
  #F(x;k,0)=(x^(k-1)*e^(-x/0))/(0^k*L(k))
  #where k is the shape, and 0 is the scale, and L(K) is the gamma function
}

#simulation_variables
my $number_of_exteins=1000000;
my @number_of_seqs = (1..$number_of_exteins); #total number of sequences to generate
my $number_of_inteins = 100;
my @size_of_extein = (1..500)); #in AA
my $intein_size = 200; #in AA
my @size_of_intein = (1..($intein_size));
my $insert_site = 200; #where in the extein to insert the fake intein
my $number_of_generations = 100000; #generations for intein evolution
my $sub_rate = 0.0003; #chance of a substitution per site in a generation

#other globals
my (%random_seqs,%intein_sequences);

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

#mutates the original "intein" sequence based on sub rate, and generations above.
#also prevents nonsense mutations from occuring
#outputs hash of evolved intein sequences
for(my $inteins = 0; $inteins < $number_of_inteins; $inteins++;){
  $ancestral_intein_sequence = $og_intein;
  for(my $gen = 0; $gen < $number_of_generations; $gen++;){
    my @ancestor = split(undef, $ancestral_intein_sequence);
    my $newinteinseq="";
    my $codon_checker;
    for(my $i = 0; $i < $size_of_intein; $i++;){
      if(rand(1) < $sub_rate){
        my $new_nucleotide = @dna[int(rand 4)];
        $codon_checker.=$new_nucleotide;
        if(length($codon_checker)>2){
          if(grep( /^$codon_checker$/,@stop_codons)){
            $newinteinseq.=@ancestor[$i];
          }
          else{
            $newinteinseq.=$new_nucleotide;
          }
          $codon_checker="";
        }
        else{
          $newinteinseq.=$new_nucleotide;
        }
      }
      else{
        #print original nucleotide
        $newinteinseq.=@ancestor[$i];
      }
    }
    $ancestral_intein_sequence=$newinteinseq;
  }
  $intein_sequences{$inteins}=$ancestral_intein_sequence;
}

#insert the intein sequences into exteins at site determined above
foreach my $extein (values %random_seqs){
  my $exteinl = length($extein);
  m
}
