#!/usr/bin/env bash

#This code will run ICE BLAST on as a job on a slurm interface

#This is the SBATCH block.
#SBATCH --job-name=ice_GPD
#SBATCH --nodes=1
#SBATCH --qos=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20gb
#SBATCH -t 100:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=YOUR_EMAIL_ADDRESS@WHATEVER.COM
#SBATCH -o ice_blast_%j.out
#SBATCH -e ice_blast_%j.err

#NOTES: 
#I have found that even large jobs tend to use only ~13-16gb of RAM.
#I put 20 here to be safe, but you are likely fine with 16 for reasonably sized jobs
#CPU count will use most of the threads given PSI-BLAST


#dependencies
module load blast/2.11.0
module load perl
module load usearch

#put paths for the databases you intend to use here.
export BLASTDB=$BLASTDB:/home/FCAM/sgosselin/22_08_01_GPD_ice/database
export BLASTDB=$BLASTDB:/home/FCAM/sgosselin/uniref50
#you most likely don't need to add these exports, but if ICEBLAST can't find your DB
#and you've already given it the full path to the DB try adding that path here!
export PATH=~/home/FCAM/sgosselin/22_08_01_GPD_ice/database:$PATH
export PATH=~/home/FCAM/sgosselin/uniref50:$PATH

#refer to the perl script or -H for information on options invoked
#make sure to include full path to databases, since ICE BLAST checks them for R/W privs
perl iceblast.pl -in intein_centroids.fasta -psidb /home/FCAM/sgosselin/uniref50/uniref50.fasta -outdb /home/FCAM/sgosselin/22_08_01_GPD_ice/database/GPD_proteome.faa -t 16 -id 0.60 -e 1e-20 -ds .25
