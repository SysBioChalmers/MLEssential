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
    # print(genes["YalifMp02"]) # Check why this gene id has no sequence
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
    with open("../Data/geneseq/Y_lipolytica.cds.fasta", "r") as handleGene :
        geneSeq = dict()
        geneprotein = dict() # Establish ID mapping for genes and proteins
        # It can not use Biopython here, what we need to extract is not seq id
        lines = handleGene.readlines()

        for line in lines :
            if line.startswith(">") :
                geneline = re.findall(r"[\[](.*?)[\]]", line)
                # Run an example using re.findall() function:
                # ['locus_tag=YALI0_F31559g',
                #  'db_xref=EnsemblGenomes-Gn:YALI0_F31559g,EnsemblGenomes-Tr:CAG78914,GOA:Q6BZR1,InterPro:IPR021137,UniProtKB/TrEMBL:Q6BZR1,GeneID:2907904',
                #  'protein=YALI0F31559p',
                #  'protein_id=XP_506101.1',
                #  'location=join(3919005..3919157,3919511..3919645)',
                #  'gbkey=CDS']

                if geneline[0].split("=")[0] == "locus_tag" :
                    genename = geneline[0].split("=")[1].replace("_","")
                else :
                    genename = geneline[1].split("=")[1]
                    print(genename)

                geneSeq[genename] = ""

                for gene in geneline :
                    if gene.split("=")[0] == "protein_id" :
                        protein_id = gene.split("=")[1]  # protein_id: XP_506101.1   
                        geneprotein[protein_id] = genename   # genename: YALI0F31559g

            else :
                geneSeq[genename] += line.replace("\n","").strip()

        print(len(geneSeq))
        print(geneSeq["YalifMp02"])

    # get the protein sequence according to protein sequence id
    with open("../Data/proteinseq/Y_lipolytica.peptide.fasta.faa", "r") as handleProtein :
        proteinSeq = dict()
        lines = handleProtein.readlines()

        for line in lines :
            if line.startswith(">") :
                proteiname = line.split(" ")[0][1:]  # proteiname: XP_506101.1
                proteinSeq[geneprotein[proteiname]] = ""
            else :
                proteinSeq[geneprotein[proteiname]] += line.replace("\n","").strip()

        print(len(proteinSeq))
        print(proteinSeq["YalifMp02"])

    genewrong = list() # Record all the genes that we can not find sequences
    for gene, essential in iter(genes.items()) :
        gene1 = dict()
        try :
            gene1[gene] = essential
            gene1["gene sequence"] = geneSeq[gene]
            gene1["protein sequence"] = proteinSeq[gene]
            geneAll.append(gene1)
        except :
            genewrong.append(gene)

    print(len(geneAll))
    print("The genes that we can not find sequences: %s" % genewrong)  # ['YALI0E12442g', 'YALIfMp01', 'YALIfMp02', 'YALIfMp16', 'YALIfMp19', 'YALIfMp20', 'YALIfMp28', 'YALIfMp29']
    print(geneAll[:3])

    return geneAll


def main() :
    # yeasts = ["C. albicans", "P. pastoris", "R. toruloides", "S. cerevisiae", "S. pombe", "Y. lipolytica"]

    # for yeast in yeasts :
    #     read_xls("%s.xlsx" %(yeast))
    genes = read_xls("Y. lipolytica.xlsx")
    geneAll = get_sequence(genes)

    # write geneAll list into processed data folder
    with open("../Data/processed_data/Y_lipolytica.json", "w") as outfile :
        json.dump(geneAll, outfile, indent=4)


if __name__ == "__main__" :
    main()



