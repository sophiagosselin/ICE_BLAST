#!/usr/bin/env bash
#SBATCH --job-name=iterative_psi_xeon
#SBATCH --nodes=1
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH -t 100:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=sophia.gosselin@uconn.edu
#SBATCH -o iterative_psi_%j.out
#SBATCH -e iterative_psi_%j.err


module load blast/2.11.0
module load perl
module load uclust
#put paths for the databases you intend to use here.
export BLASTDB=$BLASTDB:/labs/Gogarten/ActinophagesSummer21
export BLASTDB=$BLASTDB:/home/FCAM/sgosselin/uniref90
export BLASTDB=$BLASTDB:/home/FCAM/sgosselin/uniref50

#refer to the perl script or -H for information on options invoked
perl ice_blast_v2.pl -in peanut_test.fst -db1 uniref50.fasta -db2 uniref50.fasta -db3 allPhams.faa -v 1 -t 32 -id 0.9

seff %j
