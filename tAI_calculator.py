from collections import defaultdict
from scipy.stats import gmean
import numpy as np
from scipy.stats import spearmanr
import pandas as pd
from Bio import SeqIO

class tAI_calculator:
    # based on wobble pairing: A-to-I editing
    # I - U
    # I - A
    # I - C
    # G - U

    # From the view of codon 
    # | codon (DNA) | codon (mRNA) | anti (crick) | anti (wobble) | wobble pair | 
    # | --- --- --- | --- --- ---  | --- --- --- -| --- --- --- --| --- --- --- |
    # |     NNA     |     NNA      |     UNN      |      ANN      |    A : I    |    
    # |     NNC     |     NNC      |     GNN      |      ANN      |    C : I    |
    # |     NNG     |     NNG      |     CNN      |      UNN      |    G : U    |
    # |     NNT     |     NNU      |     ANN      |      GNN      |    U : G    | 

    crick_dict = {
        # key: codon ; value : anticodon ; direction : 5 -> 3
        'ATA': 'TAT', 'ATC': 'GAT', 'ATT': 'AAT', 'ATG': 'CAT',
        'ACA': 'TGT', 'ACC': 'GGT', 'ACT': 'AGT', 'ACG': 'CGT',
        'AAA': 'TTT', 'AAC': 'GTT', 'AAT': 'ATT', 'AAG': 'CTT',
        'AGA': 'TCT', 'AGC': 'GCT', 'AGT': 'ACT', 'AGG': 'CCT',
        
        'CTA': 'TAG', 'CTC': 'GAG', 'CTG': 'CAG', 'CTT': 'AAG',
        'CCA': 'TGG', 'CCC': 'GGG', 'CCG': 'CGG', 'CCT': 'AGG',
        'CAA': 'TTG', 'CAC': 'GTG', 'CAG': 'CTG', 'CAT': 'ATG', 
        'CGA': 'TCG', 'CGC': 'GCG', 'CGG': 'CCG', 'CGT': 'ACG',
        
        'GTA': 'TAC', 'GTC': 'GAC', 'GTG': 'CAC', 'GTT': 'AAC',
        'GCA': 'TGC', 'GCC': 'GGC', 'GCG': 'CGC', 'GCT': 'AGC',
        'GAA': 'TTC', 'GAC': 'GTC', 'GAG': 'CTC', 'GAT': 'ATC', 
        'GGA': 'TCC', 'GGC': 'GCC', 'GGG': 'CCC', 'GGT': 'ACC',

        'TTA': 'TAA', 'TTC': 'GAA', 'TTG': 'CAA', 'TTT': 'AAA',  
        'TCA': 'TGA', 'TCC': 'GGA', 'TCG': 'CGA', 'TCT': 'AGA',
        'TAA': 'TTA', 'TAC': 'GTA', 'TAG': 'CTA', 'TAT': 'ATA',
        'TGA': 'TCA', 'TGC': 'GCA', 'TGG': 'CCA', 'TGT': 'ACA'
        # Stop codons: 'TAA': 'TTA', 'TAG': 'CTA','TGA': 'TCA', 
    }

    def __init__(self, gene_list, fasta_file:str , trna_conc, obs_exp_df, initial_parameter=None):
        self.selected_genes = gene_list     # a python list for calculating tAI
        self.fasta_file = fasta_file        # a path to the nucleotide fasta file for the gene list
        self.trna_conc = trna_conc          # a python dictionary saving the tRNA concentration, can be tRNA gene copy number (tGCN) or TPM of tRNA by mapping 
        self.obs_exp_df = obs_exp_df        # a pandas dataframe saving the observed expression level, such as proteomics, Ribo-seq, RNA-seq
        self.parameter_set = None           # a python list saving the parameter_set for estimate the constraint parameter, default is None
        self.crick_dict = self.crick_dict   # a python dictionary saving the pairing rules between codon and anticodon 
        if initial_parameter is not None:
            self.initial_parameter = initial_parameter
        else:
            self.initial_parameter = [0.5, 0.5, 0.5, 0.5]
            
        # the parameter sets means the parameter of constraint / stability of wobble pairing
        # stands for A, C, G, T in the 3rd position of the codon 
        # i.e. the first 0.5 represent the staibility between NNA and ANN is 0.5 
        
    def W_dict(self):
        # this function aims to calculate a dictionary of absolute adaptiveness
        # self.parameter_set = {"A": self.initial_parameter[0], "C": self.initial_parameter[1], "G": self.initial_parameter[2], "T": self.initial_parameter[3]}
        self.W_dict = defaultdict(float)

        for codon in list(self.crick_dict.keys()):
            # if the codon is stop codon, give the adaptiveness to them as 0 
            if codon in ["TAA", "TGA", "TAG"]:
                self.W_dict[codon] = 0

            # else if not stop codon
            else: 
                wobble_base = codon[2]
                crick_pair = self.crick_dict[codon]
                if wobble_base == "A":
                    wobble_pair = "A" + crick_pair[1:3]
                    s = self.initial_parameter[0]
                elif wobble_base == "C":
                    wobble_pair = "A" + crick_pair[1:3]
                    s = self.initial_parameter[1]
                elif wobble_base == "G":
                    wobble_pair = "T" + crick_pair[1:3]
                    s = self.initial_parameter[2]
                elif wobble_base == "T":
                    wobble_pair = "G" + crick_pair[1:3]
                    s = self.initial_parameter[3]
                
                W = self.trna_conc[crick_pair] + (1-s) * self.trna_conc[wobble_pair]
                self.W_dict[codon] = W

        return self.W_dict
    
    def w_dict(self):
        self.w_dict = defaultdict(float)
        # W_max = max(self.W_dict.values()) # Used in dos reis 2004
        # W_max = mean(self.W_dict.values()[0:3])
        W_max = gmean(sorted(self.W_dict.values(), reverse=True)[0:3]) # Take the geometric mean of the top 3 tRNA concentration

        for codon in list(self.crick_dict.keys()):
            # if codon is stop codon, gave them relative adaptiveness of 0
            if codon in ["TAA", "TGA", "TAG"]:
                self.w_dict[codon] = 0
            # if codon is not stop codon and not zero , gave them as W[codon] / W_max
            else:
                if self.W_dict[codon] != 0:
                    self.w_dict[codon] = min(1, self.W_dict[codon] / W_max)

        non_zero_values = [v for v in self.w_dict.values() if v != 0]
        w_gmean = gmean(non_zero_values)
#         w_mean = np.average(list(self.w_dict.values())) # this mean need to be modified 

        for codon in list(self.w_dict.keys()):
            if self.W_dict[codon] == 0:
                self.w_dict[codon] = w_gmean
 
        return self.w_dict
    
    def tAI_calculate(self):
        self.tAI_dict = defaultdict(float)
        # self.error_target_list = []

        target_records_dict = SeqIO.to_dict(SeqIO.parse(self.fasta_file, "fasta"))

        for gene in self.selected_genes:
            sequences = str(target_records_dict[gene].seq)
            product_w = 1
            n_codons = 0

            if len(sequences) % 3 != 0:
                # print(f"The length of {gene} cannot be exactly divided by 3")
                # self.error_target_list.append(gene)
                continue
            else:
                codon_counter = defaultdict(int)
                for i in range(3, len(sequences), 3):
                    codon_code = sequences[i: i+3]
                    codon_counter[codon_code] += 1
                
                for codon in self.w_dict.keys():
                    if codon in ["TAA", "TGA", "TAG"]:
                        continue
                    else:
                        if codon_counter[codon] > 0:
                            product_w *= self.w_dict[codon] ** codon_counter[codon]
                            n_codons += codon_counter[codon]

                tAI = product_w ** (1/n_codons)
                self.tAI_dict[gene] = tAI
        
        # print(f"There are {len(self.error_target_list)} has wrong length in the target gene set")
        # print("The error gene list is:", self.error_target_list)    

    def get_tAI_dict(self):
        return self.tAI_dict
    
    def test_corr(self):
        self.correlation = 0
        self.p_value = 0
        # test Spearman's correlation between tAI of gene and the obs expression
        tAI_df = pd.DataFrame.from_dict(self.tAI_dict, orient="index", columns=[("tAI")]).reset_index()
        tAI_df.columns = ["gene", "tAI"]
        merge_df = pd.merge(self.obs_exp_df, tAI_df, on="gene")
        self.correlation, self.p_value = spearmanr(merge_df.iloc[:,1], merge_df.iloc[:,2])
    
    def get_corr(self):
        return self.correlation
    
    def get_pvalue(self):
        return self.p_value
    
    @classmethod
    def all(cls, gene_list, fasta_file:str , trna_conc, obs_exp_df, initial_parameter=None):
        analysis = cls(gene_list, fasta_file, trna_conc, obs_exp_df, initial_parameter=initial_parameter)
        analysis.W_dict()
        analysis.w_dict()
        analysis.tAI_calculate()
        tAI_dict = analysis.get_tAI_dict()
        analysis.test_corr()        
        spearman_corr = analysis.get_corr() 
        p_value = analysis.get_pvalue()
        return tAI_dict, spearman_corr, p_value