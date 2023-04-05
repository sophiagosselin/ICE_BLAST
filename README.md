# ICE-BLAST
ICE-BLAST (Iterative Cluster Expansion BLAST) is designed to capture large numbers of very divergent homologs of a given query within your database of choice. 

Dependencies:

    -usearch (v11.0.667)
    -BLAST (2.11.0)
    -perl (5.30.1)

You will need at 2 blastable amino acid sequence databses. One will be used to make the PSSM which will then be used to search and extract sequence matches from the second. The databse you wish to pull sequences from **must** be formated with the parse_seqids option during database creation. The code will not function otherwise. Additionally, while the database you use to make the PSSM does not need to use the parse_seqids option, it should be a clustered database (i.e. UniRef) such that search times are not absurdly long. Note that if you build a custom database for this purpose that you will bias your search to sequences present in said database which may limit the scope of your results.

You should be able to just use any up to date version of the dependencies, but versions used in development are provided.

Quick Example:

    perl iceblast.pl -in sample.fasta -psidb db_for_pssm.fasta -outdb db_to_search.fasta

You can add options as desired, see -H for more information.

Importantly, if you are interested in just a portion of your query make sure to use the domain specific option (ds). This will ensure that you only caputre the domain or region of interest and will confine your searches to a size window of 80-120% of the average length of your original queries.
