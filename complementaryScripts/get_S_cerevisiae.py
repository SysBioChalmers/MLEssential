#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-12-09

# This python script is to obtain coding sequence and protein sequence for essential genes and non-essential genes

import os
import json
from collections import Counter
from Bio import SeqIO
import pandas as pd


# read the yeast file
def read_xls(filename) :
    data = pd.read_excel("../Data/essential/%s" %(filename))
    df = pd.DataFrame(data)
    genes = df.set_index("Gene id")["Essential genes (E) / Non-essential genes (NE)"].to_dict()
    # print(genes)
    print("The number of genes collected from experimental data is: %d" %(len(genes)))

    # The number of essential genes and non-essential genes, return results like Counter({'NE': 1714, 'E': 633})
    numbers = Counter(genes.values())

    # output the E and NE gene number for a specific file
    print("The number of essential genes for %s is: %d" %(filename[:-5],numbers["E"]))
    print("The number of non-essential genes for {} is: {}" .format(filename[:-5],numbers["NE"]))
    return genes


# get the coding sequence and protein sequence for each gene id
def get_sequence(genes) :
    geneAll = list()
    # get the gene sequence accoding to gene sequence id
    with open("../Data/geneseq/S_cerevisiae.cds.fasta", "rU") as handleGene :
        geneSeq = dict()
        for record in SeqIO.parse(handleGene, "fasta") :
    # ['__add__', '__bool__', '__class__', '__contains__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__',
    # '__getattribute__', '__getitem__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__iter__', '__le__', '__le___', '__len__', '__lt__', 
    # '__module__', '__ne__', '__new__', '__nonzero__', '__radd__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', 
    # '__subclasshook__', '__weakref__', '_per_letter_annotations', '_seq', '_set_per_letter_annotations', '_set_seq', 'annotations', 'dbxrefs', 'description', 
    # 'features', 'format', 'id', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'translate', 'upper']
            geneSeq[record.id] = str(record.seq)
        print(len(geneSeq))
        print(geneSeq["YDL217C"])

    # get the protein sequence according to protein sequence id
    with open("../Data/proteinseq/S_cerevisiae.peptide.fasta", "rU") as handleProtein :
        proteinSeq = dict()
        for record in SeqIO.parse(handleProtein, "fasta") :
            proteinSeq[record.id] = str(record.seq)
        print(len(proteinSeq))
        print(proteinSeq["YDL217C"])

    for gene, essential in iter(genes.items()) :
        gene1 = dict()
        try :
            gene1[gene] = essential
            gene1["gene sequence"] = geneSeq[gene]
            gene1["protein sequence"] = proteinSeq[gene]
            geneAll.append(gene1)
        except :
            pass

    print(len(geneAll))
    # print(geneAll[:3])

    return geneAll


def main() :
    # yeasts = ["C. albicans", "P. pastoris", "R. toruloides", "S. cerevisiae", "S. pombe", "Y. lipolytica"]

    # for yeast in yeasts :
    #     read_xls("%s.xlsx" %(yeast))
    genes = read_xls("S. cerevisiae.xlsx")
    geneAll = get_sequence(genes)

    # write geneAll list into processed data folder
    with open("../Data/processed_data/S_cerevisiae.json", "w") as outfile :
        json.dump(geneAll, outfile, indent=4)


if __name__ == "__main__" :
    main()



