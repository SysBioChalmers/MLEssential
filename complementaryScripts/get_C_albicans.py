#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-12-11

# This python script is to obtain coding sequence and protein sequence for essential genes and non-essential genes

import os
import re
import json
from collections import Counter
from Bio import SeqIO
import pandas as pd


# read the yeast file
def read_xls(filename) :
    data = pd.read_excel("../Data/essential/%s" %(filename))
    df = pd.DataFrame(data)
    genes = df.set_index("Gene id")["Essential genes (E) / Non-essential genes (NE)"].to_dict()
    print("The number of genes collected from experimental data is: %d" %(len(genes)))

    # The number of essential genes and non-essential genes, return results like Counter({'NE': 1714, 'E': 633})
    numbers = Counter(genes.values())

    # output the E and NE gene number for a specific file
    print("The number of essential genes for %s is: %d" %(filename[:-5],numbers["E"]))
    print("The number of non-essential genes for {} is: {}" .format(filename[:-5],numbers["NE"]))
    return genes

# obtain protein sequence by coding sequence.
def translate_DNA(seq) :
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""

    if len(seq)%3 == 0:
       
            for i in range(0, len(seq), 3):
               
                codon = seq[i:i + 3]
                
                protein+= table[codon]

    return protein

# get the coding sequence and protein sequence for each gene id
def get_sequence(genes) :
    geneAll = list()
    # get the gene sequence accoding to gene sequence id
    with open("../Data/geneseq/C_albicans.cds.fasta", "rU") as handleGene :
        geneSeq = dict()
        for record in SeqIO.parse(handleGene, "fasta") :
    # ['__add__', '__bool__', '__class__', '__contains__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__',
    # '__getattribute__', '__getitem__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__iter__', '__le__', '__le___', '__len__', '__lt__', 
    # '__module__', '__ne__', '__new__', '__nonzero__', '__radd__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', 
    # '__subclasshook__', '__weakref__', '_per_letter_annotations', '_seq', '_set_per_letter_annotations', '_set_seq', 'annotations', 'dbxrefs', 'description', 
    # 'features', 'format', 'id', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'translate', 'upper']
            geneSeq[record.id] = str(record.seq)

        print(len(geneSeq))
        print(geneSeq["orf19.5964"])

    # get the protein sequence according to protein sequence id
    # with open("../Data/proteinseq/C_albicans.peptide.fasta", "r") as handleProtein :
    #     proteinSeq = dict()
    #     lines = handleProtein.readlines()
    #     print(len(lines))

    #     for line in lines :
          
    #         proteinannotation = line.split("\t")[:-1]   # Results: ['A0A1D8PHU1', 'A0A1D8PHU1_CANAL', 'orf19.2277', 'CAALFM_C207210CA']
    #         # But Results also exist like this : ['Q5A397', 'Q5A397_CANAL', 'CAALFM_CR08090WA orf19.6367', '']
    #         # print(proteinannotation)
    #         proteinseq= line.split("\t")[-1].strip()

    #         for protein in proteinannotation :
    #             if protein.startswith("orf") :
    #                 proteinSeq[protein] = proteinseq
    #             if protein.split(" ") :
    #                 sub_proteinannotation = protein.split(" ")
    #                 for sub_protein in sub_proteinannotation :
    #                     if sub_protein.startswith("orf") :
    #                         proteinSeq[sub_protein] = proteinseq
    #     print(proteinSeq["orf19.5964"])

    # We have not find over 500 protein sequenes related to gene name or id according to the suggestions as above.

    # proteinSeq = dict()

    # for gene in iter(genes.items()) :
    #     try :
    #         proteinSeq[gene] = translate_DNA(geneSeq[gene])


    # # print(len(proteinSeq))
    # print(proteinSeq["orf19.5964"])

    genewrong = list() # Record all the genes that we can not find sequences
    for gene, essential in iter(genes.items()) :
        gene1 = dict()
        try :
            gene1[gene] = essential
            gene1["gene sequence"] = geneSeq[gene]
            # gene1["protein sequence"] = proteinSeq[gene]
            gene1["protein sequence"] = translate_DNA(gene1["gene sequence"])[:-1]
            geneAll.append(gene1)
        except :
            genewrong.append(gene)

    print(len(geneAll))
    print("The genes that we can not find sequences: %s" % genewrong)  

    return geneAll


def main() :
    # yeasts = ["C. albicans", "P. pastoris", "R. toruloides", "S. cerevisiae", "S. pombe", "Y. lipolytica"]

    # for yeast in yeasts :
    #     read_xls("%s.xlsx" %(yeast))
    genes = read_xls("C. albicans.xlsx")
    geneAll = get_sequence(genes)

    # write geneAll list into processed data folder
    with open("../Data/processed_data/C_albicans.json", "w") as outfile :
        json.dump(geneAll, outfile, indent=4)


if __name__ == "__main__" :
    main()



