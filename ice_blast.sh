#!/usr/bin/env bash
#SBATCH --job-name=ice_GPD
#SBATCH --nodes=1
#SBATCH --qos=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
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
export BLASTDB=$BLASTDB:/home/FCAM/sgosselin/22_08_01_GPD_ice/database
export BLASTDB=$BLASTDB:/home/FCAM/sgosselin/uniref50

#refer to the perl script or -H for information on options invoked
perl iceblast.pl -in intein_centroids.fasta -psidb uniref50.fasta -outdb GPD_proteome.faa -t 16 -id 0.60 -e 1e-20 -ds

