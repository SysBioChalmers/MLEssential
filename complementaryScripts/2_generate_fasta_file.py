#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-12-30

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
        print("The protein number of %s is: %d" % (strain,len(proteinSeq)))
        # for key in proteinSeq.keys() :
        #     print(key)
        #     print(proteinSeq[key])
    return proteinSeq

def get_querySeq(organism) :
    with open("../json/%s.json" % organism) as f :
    # with open("../data/processed_data/C_albicans.json") as f :
       proteinAll = json.load(f)
    
    return proteinAll


def main() :
    # proteinSeq = get_refSeq("Candida_albicans")   
    # strains = ["Candida_albicans", "Candida_glabrata", "Candida_dubliniensis", "Candida_parapsilosis", "Candida_tropicalis", "Yarrowia_lipolytica", "Schizosaccharomyces_pombe", "Saccharomyces_cerevisiae"]
    # for strain in strains :
    #     proteinSeq = get_refSeq(strain)
    #     with open("../reciprocal_blast/%s_ref.fasta" % strain, "w") as f :
    #         for proteinId in proteinSeq.keys() :
    #             f.write(">%s\n" % (proteinId))
    #             f.write("%s\n" % (proteinSeq[proteinId]))

    organisms = ["Candida_albicans", "Candida_glabrata", "Candida_dubliniensis", "Candida_parapsilosis", "Candida_tropicalis", "Yarrowia_lipolytica", "Schizosaccharomyces_pombe", "Saccharomyces_cerevisiae"]
    for organism in organisms :
        proteinAll = get_querySeq(organism)
        with open("../reciprocal_blast/%s_query.fasta" % organism, "w") as f :
            for protein in proteinAll :
                f.write(">%s\n" % (list(protein.keys())[0]))
                f.write("%s\n" % (protein["protein sequence"]))

if __name__ == "__main__" :
    main()

# Results for ref fasta:
# The protein number of Candida_albicans is: 6207
# The protein number of Candida_glabrata is: 5143
# The protein number of Candida_dubliniensis is: 5949
# The protein number of Candida_parapsilosis is: 5280
# The protein number of Candida_tropicalis is: 5975
# The protein number of Yarrowia_lipolytica is: 6433
# The protein number of Schizosaccharomyces_pombe is: 5134
# The protein number of Saccharomyces_cerevisiae is: 5911
# [Finished in 90.5s]


