# EditPredict
EditPredict is a sequence-only, sequencing-independent tool, which is not purported to replace sequencing technology as the definitive method for detecting RNA editing events.
# How EditPredict works
GIREMI calculates the mutual information (MI) of the mismatch pairs identified in the RNA-seq reads to distinguish RNA editing sites and SNPs. It also trains a generalized linear model (GLM) to achieve enhanced predictive power, which makes use of sequence bias information and the difference between the mismatch ratio of the unknown single nucleotide variants (SNVs) and the estimated allelic ratio of the gene.
# Availability of data and materials
EditPredict web application is freely accessible to the public at http://www.innovebioinfo.com/Sequencing_Analysis/RNAediting/RNA1.php
# Instructions for use
1. get sequences: $ python get_seq.py -f [option1] -p [option1] -m [option1] -l [option1]
      * -f, --fasta-ref FILE reference genome sequence file in fasta format (NOTE: the faidx index file generated by samtools should be saved in the same         directory as this fasta file)
      * -p,
      * -m,
      * -l, 


2. 

