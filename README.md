# tAI calculator

Goal: Calculate the tRNA adaptation index of any given sequences

### What will be done

1. Module 1: Estimate parameter from scratch
2. Module 2: tAI calculation for the target gene list

### Preparation

- Observed expressional profile:
  - A panda dataframe where the first column is the gene id and the second column is the observed expressional level, such as RPKM of RNA-seq
  - NOTE: the gene id should match with the gene list and the fasta file you provide
- tRNA concentration table: 
  - Format: A dictionary where the key is the anticodon from 5 to 3 in DNA nucleotide and the value is the value representing its concentration. 
  - Concentration can be tRNA gene copy number (tGCN) or TPM based on the tRNA-sequencing
- Gene lists
  - `gene_list_train`: The list of genes for estimating the parameter set (the stability of wobble pairing in calculating adaptiveness S)
  - `gene_list_target`: The list of genes that you want tAI values
- Fasta files
  - Specify a fasta file that saves the DNA sequences of the lists you provided
  - NOTE: Make sure the genes in fasta file, gene list, and genes in the observed expressional profile match each other. 
- Parameter set (optional)
  - If you have already got a list of parameters so that you do not need to estimate them. You can directly provide it use `initial_parameter=<your parameter>` in the argument

### Usage