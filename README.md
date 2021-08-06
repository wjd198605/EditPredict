# EditPredict
EditPredict is a method that can identify RNA editing sites using one RNA-seq data set without requiring genome sequence data.
# How EditPredict 
GIREMI calculates the mutual information (MI) of the mismatch pairs identified in the RNA-seq reads to distinguish RNA editing sites and SNPs. It also trains a generalized linear model (GLM) to achieve enhanced predictive power, which makes use of sequence bias information and the difference between the mismatch ratio of the unknown single nucleotide variants (SNVs) and the estimated allelic ratio of the gene.
