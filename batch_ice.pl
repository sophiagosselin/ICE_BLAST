#!/usr/bin/perl -w
use strict;
use warnings;
use File::Copy;

my @clusterids=(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9);
my @evals=(1e-3,1e-5,1e-10,1e-20);

foreach my $evalue (@evals){
  foreach my $id (@clusterids){
    my $dir = "$evalue\_$id";
    mkdir($dir);
    copy("intein.fasta","$dir/intein.fasta");
    copy("iceblast_v1.1.pl","$dir/iceblast_v1.1.pl");
    open(SH, "+> $dir/ice.sh");
    print SH "#!/usr/bin/env bash
#SBATCH --job-name=ice_intein
#SBATCH --nodes=1
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH -t 100:00:00
#SBATCH --chdir=$dir
#SBATCH --mail-type=END
#SBATCH --mail-user=sophia.gosselin\@uconn.edu
#SBATCH -o ice_ds_\%j.out
#SBATCH -e ice_ds_\%j.err

module load blast/2.11.0
module load perl/5.30.1
module load uclust
#put paths for the databases you intend to use here.
export BLASTDB=\$BLASTDB:/labs/Gogarten/ActinophagesSummer21
export BLASTDB=\$BLASTDB:/home/FCAM/sgosselin/uniref50

#refer to the perl script or -H for information on options invoked
perl iceblast_v1.1.pl -in intein.fasta -psidb uniref50.fasta -outdb allPhams.faa -t 8 -id $id -ds 1 -e $evalue";
    close SH;
    system("sbatch $dir/ice.sh");
  }
}
