#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-12-30

# This python script is to obtain and summarize protein-protein complex data


import json
import re
import csv
from collections import defaultdict
from urllib import request


# mapping the protein id and protein sequence from the fasta file.
def id_sequence(strain) :
    # get the protein sequence accoding to protein sequence id
    with open("../fasta/%s.fasta" % strain, "r") as handleProtein :
        lines = handleProtein.readlines()

    SeqId = dict()
    IdSeq = dict()

    for line in lines :
        if line.startswith(">") :
            proteiname = line.split(" ")[0][1:]  # proteiname: XP_506101.1
            IdSeq[proteiname] = ""
        else :
            IdSeq[proteiname] += line.replace("\n","").strip()

    # print(IdSeq)
    SeqId = {value:key for key,value in IdSeq.items()}
        
    return IdSeq,SeqId

# This function is to obtain the protein sequence according to the protein id from Uniprot API
# https://www.uniprot.org/uniprot/A0A1D8PIP5.fasta 
# https://www.uniprot.org/help/api_idmapping
def uniprot_sequence(id) :
    url = "https://www.uniprot.org/uniprot/%s.fasta" % id
    IdSeq = dict()

    try :
        data = request.urlopen(url)
        respdata = data.read().decode("utf-8").strip()
        IdSeq[id] =  "".join(respdata.split("\n")[1:])
    except :
        print(id, "can not find from uniprot!")
        IdSeq[id] = None
    
    return IdSeq[id]

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
        IdSource[data[3]] = data[1],data[4],data[0]
        complexData1.append(IdSource)

    complexData2 = defaultdict(list)
    default = [complexData2[k].append(v) for complex in complexData1 for k, v in complex.items()]
    complexData = dict(complexData2)
    for k, v in dict(complexData2).items() :
        complexData[k] = set(v) # delete the data which has the same values
    # print(complexData)
    # Example: {'GO:0045259': {('Q36258', 'proton-transporting ATP synthase complex', 'UniProtKB')}}

    # print(complexData)
    IdSeq,SeqId = id_sequence(strain)

    complexSeq = dict()
    for k1 in complexData.values() :
        for k2 in k1 :
            if k2[2] == "UniProtKB" :
                # complexSeq[k2[0]] = uniprot_sequence(k2[0])
                uniprot = uniprot_sequence(k2[0])
                try :
                    new_id = SeqId[uniprot]
                    complexSeq[k2[0]] = new_id  # new_id refers to fasta id, k2[0] refers to uniprot id
                except :
                    print("Uniprot ID mapping not correct: %s" % k2[0])

    with open("../processed_data/%s_include_ortholog.json" % strain.replace(" ", "_"), "r") as f:
        allData = json.load(f)
        # print(allData[:2])

    ids = [list(data.keys())[0] for data in allData]
    # print(ids[:100])
    dataComplex = [data["id"] for data in allData if list(data.values())[0] == "Complex"]

    allComplex = set()
    for k1 in list(complexData.values()) :
        # print(k1)
        for k2 in k1 :
            try : 
                if complexSeq[k2[0]] in ids :  # k2[0] = uniprot id, such as 'Q36258'
                    print(complexSeq[k2[0]])
                    index = ids.index(complexSeq[k2[0]])
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
            except :
                # print("Dict complexSeq mapping not pass!")
                pass

    with open("../protein_complex/%s.tsv" % strain, "wt") as outfile :
        tsv_writer = csv.writer(outfile, delimiter="\t")
        tsv_writer.writerow(["complex type", "sequence id", "index", "ortholog"])

        for data in sorted(allComplex) :
            tsv_writer.writerow(list(data))

    print("The amount of protein-containing data without duplication is: %d" % len(dataComplex))
    

def main() :
    strains = ["Candida albicans", "Candida glabrata", "Candida dubliniensis", "Candida parapsilosis", "Candida tropicalis", "Yarrowia lipolytica", "Schizosaccharomyces pombe", "Saccharomyces cerevisiae"]
    # uniprot_sequence("Q37695")
    for strain in strains[5:6] :
        # id_sequence(strain)
        write_complex_tsv(strain)


if __name__ == "__main__" :
    main()


