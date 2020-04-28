#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-12-30

# This python script is to obtain and summarize protein-protein complex data


import json
import re
import csv
from collections import defaultdict

                
# Write the protein-containing complex data, including complex type, sequence id, index, ortholog into a tsv file 
def write_complex_tsv(strain) :
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
        IdSource[data[3]] = data[1],data[4]
        complexData1.append(IdSource)

    complexData2 = defaultdict(list)
    default = [complexData2[k].append(v) for complex in complexData1 for k, v in complex.items()]
    complexData = dict(complexData2)
    for k, v in dict(complexData2).items() :
        complexData[k] = set(v) # delete the data which has the same values
    # print(complexData)
    # Example: {'GO:0000015': {('CAL0000185645', 'CGD')}, 'GO:0000109': {('CAL0000199257', 'CGD'), ('CAL0000195175', 'CGD')}}


    with open("../processed_data/%s_include_ortholog.json" % strain.replace(" ", "_"), "r") as f:
        allData = json.load(f)
        # print(allData[:2])

    ids = [list(data.keys())[0] for data in allData]
    dataComplex = set([data["id"] for data in allData if list(data.values())[0] == "Complex"])

    allComplex = set()
    for k1 in list(complexData.values()) :
        # print(k1)
        for k2 in k1 :
            if k2[0] in ids :
                index = ids.index(k2[0])
                # print(index)
                # print(k2[1])
                # print(allData[index]["id"])
                # 0
                # nucleotide-excision repair complex
                # Candida_albicans@C1_02810W_A
                # tsv_writer.writerow([k2[1], allData[index]["id"], allData[index]["index"], allData[index]["ortholog"]])
                allComplex.add((k2[1], allData[index]["id"], allData[index]["index"], allData[index]["ortholog"]))
            else :
                print("%s can not find because of no sequence mapping " % k2[0])

    with open("../protein_complex/%s.tsv" % strain, "wt") as outfile :
        tsv_writer = csv.writer(outfile, delimiter="\t")
        tsv_writer.writerow(["complex type", "sequence id", "index", "ortholog"])

        for data in sorted(allComplex) :
            tsv_writer.writerow(list(data))

    print("The amount of protein-containing data without duplication is: %d" % len(dataComplex))
    

def main() :
    strains = ["Candida albicans", "Candida glabrata", "Candida dubliniensis", "Candida parapsilosis", "Candida tropicalis", "Schizosaccharomyces pombe", "Saccharomyces cerevisiae"]
    # uniprot_sequence("A0A1D8PIP5")
    for strain in strains[1:] :
        # id_sequence(strain)
        write_complex_tsv(strain)


if __name__ == "__main__" :
    main()


