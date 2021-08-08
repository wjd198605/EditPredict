# EditPredict
EditPredict is a sequence-only, sequencing-independent tool, which can be used stand-alone to predict novel RNA editing and to assist in filtering for candidate RNA editing detected from RNA-Seq data
# How EditPredict works
EditPredict calculates novel pattern identified in flanking sequences reads to distinguish RNA editing sites and SNPs. EditPredict is a Deep Learning solution of RNA editing prediction from genome DNA sequence. Specifically, to utilize CNN to model flanking DNA sequences of known A-to-I RNA editing events, with consideration of alternate flanking directions and variant sequence lengths.
# Instructions for use
1. ### get DNA sequences: $ python get_seq.py -f [option1] -p [option2] -m [option3] -l [option4] -v [options5]
      * -f, --human reference genome File (hg38)
      * -p, --positions FILE the list of chromosome locus (see our paper for details)
      * -m, --genome Directionality the end-to-end chemical orientation of a single strand of nucleic acid at certain locus(b -- bidirectional, l -- upstream, r -- downstream)
      * -l, --length of the DNA sequences 
      * -v, --vcf format File optional for SNPs (Note: not required)


2. ### predict RNA editing: $ python editPredict.py -f [option1] -c [option2] -w [option3]
     * -f, --RNA sequences File generated from instruction 1 
     * -c, --model construction File 
     * -w, --model weight File

# Availability of data and materials
EditPredict web application is freely accessible to the public at http://www.innovebioinfo.com/Sequencing_Analysis/RNAediting/RNA1.php
