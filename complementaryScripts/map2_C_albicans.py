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
            if record.id.startswith("Candida_albicans") :
                proteinSeq[record.seq] = str(record.id)
        print(len(proteinSeq))
        # for key,value in proteinSeq.items() :
        #     print(key)
        #     print(value)
    return proteinSeq

def get_processdata() :
    with open("../data/processed_data/C_albicans.json") as f :
        geneAll = json.load(f)

    proteinSeq = get_sequence()
    i1 = 0
    i2 = 0
    gene_mapping = list()
    for gene in geneAll[:100] :
        i1 += 1
        print('i1 =',i1)
        # print(gene)
        sequence = gene["protein sequence"]
        i3 = 0
        ratio_all = list()
        for key_seq in list(proteinSeq) :
            seq = SequenceMatcher(None, sequence, key_seq)
            ratio = seq.ratio()

            # For the first ratio calculation, the program will be run as below
            if ratio > 0.4 and i3 == 0: # How about two or three samples ratio > 0.5?
                i3 += 1
                ratio1 = ratio
                ratio_all.append(ratio1)
                i2 += 1
                print('i2 =',i2)
                # print(ratio)
                print(sequence)
                print(key_seq)
                # print(proteinSeq[key_seq])
                gene['id'] = proteinSeq[key_seq]

            # From the second ratio > intinial ratio (0.5), 
            # if new ratio > the first ratio, the dict value will be changed, and i2 doesn't need to be added one
            # Otherwise, the dict value will not be changed,nothing we need do for this.
            # Another idea after testing many times, the new ratio should be compare with the largest ratio.
            if i3 != 0 :
                ratio2 = ratio
                if ratio2 > 0.4 :
                    ratio_all.append(ratio2)
                if ratio2 >= sorted(ratio_all)[-1] and ratio1 < 0.95 :  
                    print('i2 =',i2)
                    # print(ratio2)
                    print(sequence)
                    print(key_seq)
                    # print(proteinSeq[key_seq]) 
                    print(sorted(ratio_all))  
                    gene['id'] = proteinSeq[key_seq]

        gene_mapping.append(gene)

    print(len(geneAll))
    print(i2)
    # print(len(gene_mapping))

    # with open("../Data/processed_data/C_albicans_include_id2.json", "w") as outfile :
    #     json.dump(gene_mapping, outfile, indent=4)
                
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


