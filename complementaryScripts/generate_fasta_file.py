#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-12-22

# This python script is to generate the fasta file used to implement Reciprocal BLAST Hits (RBH)


import json
from Bio import SeqIO



def get_refSeq(strain) :
    # get the protein sequence accoding to protein sequence id
    with open("/Users/leyu/Documents/coding/evolution_code/orthomcl_output/343taxa_proteins.fasta", "rU") as handleGene :
        proteinSeq = dict()
        for record in SeqIO.parse(handleGene, "fasta") :
    # ['__add__', '__bool__', '__class__', '__contains__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__',
    # '__getattribute__', '__getitem__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__iter__', '__le__', '__le___', '__len__', '__lt__', 
    # '__module__', '__ne__', '__new__', '__nonzero__', '__radd__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', 
    # '__subclasshook__', '__weakref__', '_per_letter_annotations', '_seq', '_set_per_letter_annotations', '_set_seq', 'annotations', 'dbxrefs', 'description', 
    # 'features', 'format', 'id', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'translate', 'upper']
            # if record.id.startswith("Candida_albicans") :
            if record.id.startswith(strain) :
                proteinSeq[record.id] = str(record.seq)
        print("The protein number of this strain is:", len(proteinSeq))
        # for key in proteinSeq.keys() :
        #     print(key)
        #     print(proteinSeq[key])
    return proteinSeq

def get_querySeq(organism) :
    with open("../data/processed_data/%s.json" % organism) as f :
    # with open("../data/processed_data/C_albicans.json") as f :
        geneAll = json.load(f)
    
    return geneAll


def main() :
    # proteinSeq = get_refSeq("Candida_albicans")   
    # strains = ["Yarrowia_lipolytica", "Schizosaccharomyces_pombe", "Saccharomyces_cerevisiae", "Komagataella_pastoris"]  # Komagataella pastoris (Pichia pastoris)
    # for strain in strains :
    #     proteinSeq = get_refSeq(strain)
    #     with open("../Data/reciprocal_blast/%s_ref.fasta" % strain, "w") as f :
    #         for proteinId in proteinSeq.keys() :
    #             f.write(">%s\n" % (proteinId))
    #             f.write("%s\n" % (proteinSeq[proteinId]))

    organisms = {"Y_lipolytica" : "Yarrowia_lipolytica", "S_pombe" : "Schizosaccharomyces_pombe", "S_cerevisiae" : "Saccharomyces_cerevisiae", "P_pastoris" : "Komagataella_pastoris"}
    for organism in organisms.keys() :
        geneAll = get_querySeq(organism)
        with open("../Data/reciprocal_blast/%s_query.fasta" % organisms[organism], "w") as f :
            for protein in geneAll :
                f.write(">%s\n" % (list(protein.keys())[0]))
                f.write("%s\n" % (protein["protein sequence"]))

if __name__ == "__main__" :
    main()


