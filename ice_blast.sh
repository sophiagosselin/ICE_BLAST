#!/usr/bin/env bash
#SBATCH --job-name=iterative_psi
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256gb
#SBATCH -t 100:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=sophia.gosselin@uconn.edu
#SBATCH -o ice_blast_%j.out
#SBATCH -e ice_blast_%j.err

#dependencies
module load blast/2.11.0
module load perl
module load uclust

#put paths for the databases you intend to use here.
export BLASTDB=$BLASTDB:/labs/Gogarten/ActinophagesSummer21
export BLASTDB=$BLASTDB:/home/FCAM/sgosselin/uniref90
export BLASTDB=$BLASTDB:/home/FCAM/sgosselin/uniref50

#refer to the perl script or -H for information on options invoked
perl ice_blast.pl -in peanut_test.fst -db1 uniref50.fasta -db2 uniref90.fasta -db3 allPhams.faa -v 1 -t 10 -id 0.9

seff %j
