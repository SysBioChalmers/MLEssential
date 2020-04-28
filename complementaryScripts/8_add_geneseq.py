#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-03-11

# This python script is to add gene sequence for the previous processed data

import json
from Bio import SeqIO

with open("/Users/leyu/Documents/coding/evolution_code/orthomcl_output/343taxa_protein_IDs_index.txt","r") as f :
    lines = f.readlines()
    id_mapping = dict()
    for line in lines :
        data = line.strip().split("\t")
        # id_mapping[data[1]] = [data[0].replace(" ","&"), data[2].replace(" ","&")]
        id_mapping[data[1]] = data[0].replace(" ","&")

        # print(data)
    # print(id_mapping)

organisms = ["Candida_glabrata", "Candida_dubliniensis", "Candida_parapsilosis", "Candida_tropicalis", "Yarrowia_lipolytica","Candida_albicans", "Schizosaccharomyces_pombe", "Saccharomyces_cerevisiae"]

for organism in organisms[:5] :
    print("This organism is: %s" % organism.replace("_", " "))

    with open("../complementaryData/generated_fasta/%s.fasta" % organism, "rU") as allGene :
        idSeq = dict()
        for record in SeqIO.parse(allGene, "fasta") :
            # print(type(record.id)) # Type: str
            idSeq[record.id] = str(record.seq).upper()  # allGeneId is all id in fasta file

    with open("../complementaryData/processed_data/%s_include_ortholog.json" % organism, "r") as f :
        data = json.load(f)

    newdata = list()
    for complexData in data :
        modelid = complexData["id"]
        newid = id_mapping[modelid].replace("mRNA-1_1", "mRNA-1")
        # print(newid)
        # print(idSeq[newid])
        complexData["gene sequence"] = idSeq[newid]
        newdata.append(complexData)

    with open("../complementaryData/newjson/%s.json" % organism, "w") as outfile:
        json.dump(newdata, outfile, indent=4)
