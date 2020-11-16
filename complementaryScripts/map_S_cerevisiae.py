#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-12-14

# This python script is to do id mapping for the protein id obtained from my processed data and the protein id obtained from that Cell paper


import json
from Bio import SeqIO
from difflib import SequenceMatcher


def get_sequence() :
    # get the gene sequence accoding to gene sequence id
    with open("/Users/leyu/Documents/coding/evolution_code/orthomcl_output/343taxa_proteins.fasta", "rU") as handleGene :
        proteinSeq = dict()
        for record in SeqIO.parse(handleGene, "fasta") :
    # ['__add__', '__bool__', '__class__', '__contains__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__',
    # '__getattribute__', '__getitem__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__iter__', '__le__', '__le___', '__len__', '__lt__', 
    # '__module__', '__ne__', '__new__', '__nonzero__', '__radd__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', 
    # '__subclasshook__', '__weakref__', '_per_letter_annotations', '_seq', '_set_per_letter_annotations', '_set_seq', 'annotations', 'dbxrefs', 'description', 
    # 'features', 'format', 'id', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'translate', 'upper']
            if record.id.startswith("Saccharomyces_cerevisiae") :
                proteinSeq[record.id] = str(record.seq)
        print(len(proteinSeq))
        # for key,value in proteinSeq.items() :
        #     print(key)
        #     print(value)
    return proteinSeq

def get_processdata() :
    with open("../data/processed_data/S_cerevisiae.json") as f :
        geneAll = json.load(f)

    proteinSeq = get_sequence()
    i1 = 0
    gene_mapping = list()
    for gene in geneAll :
        for key_id in list(proteinSeq) :
            if list(gene.keys())[0] == key_id.split("@")[1] :
                i1 += 1
                print('i1 =',i1)
                gene["id"] = str(key_id)

        gene_mapping.append(gene)

    print(len(geneAll))
    print(i1)

    with open("../Data/processed_data/S_cerevisiae_include_id.json", "w") as outfile :
        json.dump(gene_mapping, outfile, indent=4)
                
        # for geneid, protein in proteinSeq.items() :
        #     sequence_alignment(sequence,protein)
            # if sequence == protein :
            #     i += 1
            #     print(i)
            #     print(geneid)  # can find 37 same from 100 sequences

# def sequence_alignment(X,Y) :
#     from Bio import pairwise2
#     from Bio.pairwise2 import format_alignment

#     #A match score = 2, mismatch score = -1, gap opening = -5, gap extension = -2
#     # alignments = pairwise2.align.globalms(X, Y, 2, -1, -5, -2)

#     # for a in alignments:
#     #     print(format_alignment(*a))

#     alignments = pairwise2.align.globalxx(X, Y)
#     print(format_alignment(*alignments[0]))

def sequence_match() : 

    X = "MSNFIEVQRSNIQEIGIIEDSLSRRLLRNPYILPSNLQPRPNILVSKVTKPNSTKRISILQQYELGYFKDRYNQIIKGLNHCLTDEKESFSRTLELITEPGDGFEKFDEMINGISLNDNIVEDPRKLYSSFSSLNKELNPQDVTLKINKKTGKEQEFVKRKHFLSSIASHLKNIEINDILDLHHFHDLYVSNFGPILYIEYLYKFLSFPYTNVNGFYSKYLTELSRFLEATLLKLQPLLDYNALLQNWKKEYDNANEERKSNGNDGDKLFCKACNKLFSKETVYQSHLSGKKHKKNASQQNPDNFVSSLPWLEYFIEKLCQVLAPELEYTRAQVEKLSNLSERELQLDRQVQHDIENEFVAINNEFDDDDLSQNEHGDDDDNDDYLDDSFKNLPLGPDGTPIPFWLYKLQELHKQYNCEICGNISYKGKSVFMKHFSGSKHQYGLKCLGVDEKNMKMFKNITKIDEATELWKVLRKETKLKVTELENSVEVEDKEGNVMLEKDYIDLKKQGLL"
    Y = "MSNFIEVQRSNIQEIGIIEDSLSRRLLRNPYILPSNLQPRPNILVSKVTKPNSTKRISILQQYELGYFKDRYNQIIKGLNHCLTDEKESFSRTLELITEPGDGFEKFDEMINGISSNDNIVEDPRKLYSSFSSLNKELNPQDVTLKINKKTGKEQEFVKRKHFLSSIASHLKNIEINDIMDLHHFHDLYVSNFGPISYIEYLYKFSSFPYANVNGFYSKYLTELSSFLEATLIKLQPLLDYNALLQNWKKEYDNANEERKSNGNDGDKLFCKACNKLFSKETVYQSHLSGKKHKKNASQQNPDNFVSSLPWLEYFIEKLCQVLAPELEYTRAQVEKLSNLSERELQLDRQVQHNIENEFVAINNEFDDDDLSQNEHGDDDDNDNYLDDSFKNLPLGPDGTPIPFWLYKLQGLHKQYNCEICGNISYKGKSVFMKHFSGSKHQYGLKCLGVDEKNMKMFKNITKIDEATELWKVLRKETKLKVTELENSVEVEDKEGNVMSEKDYIDLKKQGLL"

    seq = SequenceMatcher(None, X,Y)
    ratio = seq.ratio()
    print(ratio)

def main() :
    # X = "MSNFIEVQRSNIQEIGIIEDSLSRRLLRNPYILPSNLQPRPNILVSKVTKPNSTKRISILQQYELGYFKDRYNQIIKGLNHCLTDEKESFSRTLELITEPGDGFEKFDEMINGISLNDNIVEDPRKLYSSFSSLNKELNPQDVTLKINKKTGKEQEFVKRKHFLSSIASHLKNIEINDILDLHHFHDLYVSNFGPILYIEYLYKFLSFPYTNVNGFYSKYLTELSRFLEATLLKLQPLLDYNALLQNWKKEYDNANEERKSNGNDGDKLFCKACNKLFSKETVYQSHLSGKKHKKNASQQNPDNFVSSLPWLEYFIEKLCQVLAPELEYTRAQVEKLSNLSERELQLDRQVQHDIENEFVAINNEFDDDDLSQNEHGDDDDNDDYLDDSFKNLPLGPDGTPIPFWLYKLQELHKQYNCEICGNISYKGKSVFMKHFSGSKHQYGLKCLGVDEKNMKMFKNITKIDEATELWKVLRKETKLKVTELENSVEVEDKEGNVMLEKDYIDLKKQGLL"
    # Y = "KVTKPNSTKRISILQQYELGYFKDRYNQIIKGLNHCLTDEKESFSRTLELITEPGDGYLDDSFKNLPLGPDGTPIPFWLYKLQ"
    # sequence_alignment(X,Y)
    get_processdata()
    # sequence_match()
    # with open("../Data/processed_data/C_albicans.json", "w") as outfile :
    #     json.dump(geneAll, outfile, indent=4)


if __name__ == "__main__" :
    main()


