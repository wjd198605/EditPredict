# EditPredict
EditPredict is a sequence-only, sequencing-independent tool, which can be used stand-alone to predict novel RNA editing and it can be used to assist in filtering for candidate RNA editing detected from RNA-Seq data
# How EditPredict works
EditPredict calculates the mutual information (MI) of the mismatch pairs identified in the RNA-seq reads to distinguish RNA editing sites and SNPs. It also trains a generalized linear model (GLM) to achieve enhanced predictive power, which makes use of sequence bias information and the difference between the mismatch ratio of the unknown single nucleotide variants (SNVs) and the estimated allelic ratio of the gene.
# Instructions for use
1. ### get RNA sequences: $ python get_seq.py -f [option1] -p [option2] -m [option3] -l [option4]
      * -f, --human reference genome File (hg38)
      * -p, --positions FILE the list of 
      * -m, --genome Directionality the end-to-end chemical orientation of a single strand of nucleic acid at certain locus(b -- bidirectional, l -- upstream, r -- downstream)
      * -l, --length of the RNA sequences 


2. ### predict RNA editing: $ python EditPredict.py -f [option1] -c [option2] -w [option3]
     * -f, --RNA sequences File generated from instruction 1 
     * -c, --model construction File 
     * -w, --model weight File

# Availability of data and materials
EditPredict web application is freely accessible to the public at http://www.innovebioinfo.com/Sequencing_Analysis/RNAediting/RNA1.php
