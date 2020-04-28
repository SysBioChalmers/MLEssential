#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-12-30

# This python script is to obtain and summarize protein-protein complex data


import json
import re
from collections import defaultdict
from urllib import request


# mapping the protein id and protein sequence from the fasta file downloaded from CGD database
def id_sequence(strain) :
    # get the protein sequence accoding to gene sequence id
    with open("../fasta/%s.fasta" % strain, "r") as handleData :
        lines = handleData.readlines()

    IdSeq = dict()

    for line in lines :
        if line.startswith(">") :
            reg = r".*SGDID:(.*?),.*"
            seqid = re.findall(reg,line)[0]
            # print(seqid)
            IdSeq[seqid] = ""
        else :
            IdSeq[seqid] += line.replace("\n","").strip().replace("*","")

    # print(IdSeq)
    # SeqId = {value:key for key,value in IdSeq.items()}
        
    return IdSeq
                
# Obtain the protein-containing complex ids and the corresponding sequence.
def complex_sequence(strain) :
    with open("../protein/%s.tsv" % strain, "r", encoding="utf-8") as handleData:
        lines = handleData.readlines()

    complexData1 = list()
    for line in lines :
        IdSource = dict()
        data = line.strip().split("\t")
        # ['CGD', 'CAL0000185645', 'ENO1', 'GO:0000015', 'phosphopyruvate hydratase complex', 'CGD_REF:CAL0142013', 'C', '', 
        # 'C1_08500C_A|orf19.8025|IPF26917.1|IPF14429.2|Contig4-3031_0010|orf6.6269|CA3874|ENO|CaO19.395|CaO19.8025|orf19.395|C1_08500C_B|C1_08500C|CAWG_00576', 
        # gene_product', 'NCBITaxon:237561', 'CGD']
        # if data[1] not in IdSource :
        IdSource[data[3]] = data[1],data[0]
        complexData1.append(IdSource)

    complexData2 = defaultdict(list)
    default = [complexData2[k].append(v) for complex in complexData1 for k, v in complex.items()]
    complexData = dict(complexData2)
    for k, v in dict(complexData2).items() :
        complexData[k] = set(v) # delete the data which has the same values
    # print(complexData)
    # Example: {'GO:0000015': {('CAL0000185645', 'CGD')}, 'GO:0000109': {('CAL0000199257', 'CGD'), ('CAL0000195175', 'CGD')}}

    IdSeq = id_sequence(strain)

    complexSeq = dict()
    for k1 in complexData.values() :
        for k2 in k1 :
            try :
                complexSeq[k2[0]] = IdSeq[k2[0]]
            except :
                print("%s can not find from fasta file" % k2[0])


    complexAllData = [{key:"Complex", "protein sequence":value} for key,value in complexSeq.items()]
    Non_complexData = [{key:"Non-complex", "protein sequence":IdSeq[key]} for key in set(IdSeq.keys())-set(complexSeq.keys())]

    complexAll = complexAllData + Non_complexData

    print("The amount of protein-containing complex data for %s is: %d" % (strain,len(complexAll)))
    
    with open("../json/%s.json" % strain.replace(" ", "_"), "w") as outfile :
        json.dump(complexAll, outfile, indent=4)

def main() :
    strains = ["Candida albicans", "Candida glabrata", "Candida dubliniensis", "Candida parapsilosis", "Candida tropicalis", "Yarrowia lipolytica", "Schizosaccharomyces pombe", "Saccharomyces cerevisiae"]
    # uniprot_sequence("A0A1D8PIP5")
    # for strain in strains[] :
    # id_sequence(strains[-1])
    complex_sequence(strains[-1])


if __name__ == "__main__" :
    main()

# Results :
# The amount of protein-containing data for Candida albicans is: 2506 (just include complex no non-complex)
# The amount of protein-containing complex data for Candida glabrata is: 5311
# The amount of protein-containing complex data for Candida dubliniensis is: 5935
# The amount of protein-containing complex data for Candida parapsilosis is: 5863
# The amount of protein-containing complex data for Candida tropicalis is: 6258
# The amount of protein-containing complex data for Yarrowia lipolytica is: 6472
# The amount of protein-containing complex data for Schizosaccharomyces pombe is: 5135
# The amount of protein-containing complex data for Saccharomyces cerevisiae is: 6713