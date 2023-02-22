# ICE-BLAST
ICE-BLAST (Iterative Cluster Expansion BLAST) is designed to capture as many homologs of a given query within your database of choice. 

Dependencies:
- usearch (v11.0.667)
- BLAST (2.11.0)
- perl (5.30.1)

You will need at least 2 blast databses. 1 to make the PSSM and one to search + extract sequence matches from. The databse you wish to pull sequences from must be formated with the parse_seqids option during database creation. The code will not function otherwise. The database you use to make the PSSM does not need to use the parse_seqids option.

You should be able to just use any up to date version of the dependencies, but versions used in development are provided.

Quick Example:

perl iceblast.pl -in sample.fasta -psidb db_for_pssm.fasta -outdb db_to_search.fasta

You can add options as desired, see -H for more information.

Importantly, if you are interested in just a portion of your query make sure to use the domain specific option (ds). This will ensure that you only caputre the domain or region of interest and will confine your searches to a size window of 80-120% of the average length of your original queries.
