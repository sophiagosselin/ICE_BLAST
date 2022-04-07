# ICE-BLAST
ICE-BLAST (Iterative Cluster Expansion BLAST) is designed to capture as many homologs of a given query within your database of choice. 

Dependencies:
- Uclust (1.2.22q)
- BLAST (2.11.0)
- perl (5.30.1)

You will need at least 2 blast databses. 1 to make the PSSM and one to search + extract sequence matches from. The databse you wish to pull sequences from must be formated with the parse_seqids option during database creation. The code will not function otherwise. The database you use to make the PSSM does not need to use the parse_seqids option.

You should be able to just use any up to date version of the dependencies, but versions used in development are provided.

Quick Example:

perl iceblast.pl -in sample.fasta -db1 db_for_pssm.fasta -db3 db_to_search.fasta

You can add options as desired, see -H for more information.

Importantly, if you are interested in just a portion of your query make sure to use the -ds option.
