#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-03-11

# This python script is to add gene sequence for the previous processed data

import json
from Bio import SeqIO


def pombe(organism) :
    with open("../complementaryData/generated_fasta/%s.fasta" % organism, "rU") as allGene :
        idSeq = dict()
        for record in SeqIO.parse(allGene, "fasta") :
            # print(type(record.id)) # Type: str
            idSeq[record.id] = str(record.seq).upper()  # allGeneId is all id in fasta file

    with open("../complementaryData/processed_data/%s_include_ortholog.json" % organism, "r") as f:
        data =json.load(f)

    newdata = list()
    for complexData in data :
        pombeid = list(complexData.keys())[0]
        print(pombeid)
        complexData["gene sequence"] = idSeq[pombeid]
        newdata.append(complexData)

    with open("../complementaryData/newjson/%s.json" % organism, "w") as outfile:
        json.dump(newdata, outfile, indent=4)

def saccharomyces(organism) :
    with open("../complementaryData/generated_fasta/%s.fasta" % organism, "rU") as allGene :
        idSeq = dict()
        for record in SeqIO.parse(allGene, "fasta") :
            # print(type(record.id)) # Type: str
            idSeq[record.id] = str(record.seq).upper()  # allGeneId is all id in fasta file

    with open("../complementaryData/processed_data/%s_include_ortholog.json" % organism, "r") as f:
        data =json.load(f)

    newdata = list()
    for complexData in data :
        saccharomycesid = complexData["id"].split("@")[1]
        print(saccharomycesid)
        complexData["gene sequence"] = idSeq[saccharomycesid]
        newdata.append(complexData)

    with open("../complementaryData/newjson/%s.json" % organism, "w") as outfile:
        json.dump(newdata, outfile, indent=4)

def albicans(organism) :
    with open("../complementaryData/generated_fasta/%s.fasta" % organism, "rU") as allGene :
        idSeq = dict()
        for record in SeqIO.parse(allGene, "fasta") :
            # print(type(record.id)) # Type: str
            idSeq[record.id] = str(record.seq).upper()  # allGeneId is all id in fasta file

    with open("../complementaryData/processed_data/%s_include_ortholog.json" % organism, "r") as f:
        data =json.load(f)
#   
    complexData2 = [complexData["id"] for complexData in data if list(complexData.values())[0] == "Complex"]
    complexData = set(complexData2)

    allId2 = [complexData["id"] for complexData in data]
    allId = set(allId2) # allId is all the id through reciprocal blast best hits
    non_complexData = allId - complexData

    # print(len(complexData))
    # print(len(non_complexData))

    newdata = list()

    for complex_id in list(complexData) :
        for complexdict in data :
            if list(complexdict.values())[0] == "Complex" and complexdict["id"] == complex_id :
                albicansid = complexdict["id"].split("@")[1]
                complexdict["gene sequence"] = idSeq[albicansid]
                newdata.append(complexdict)
                break

    for non_complex_id in list(non_complexData) :
        for non_complexdict in data :
            if list(non_complexdict.values())[0] == "Non-complex" and non_complexdict["id"] == non_complex_id :
                albicansid = non_complexdict["id"].split("@")[1]
                non_complexdict["gene sequence"] = idSeq[albicansid]
                newdata.append(non_complexdict)
                break

    print(len(newdata))

        # complexdict["gene sequence"] = idSeq[albicansid]
    #     newdata.append(complexdict)

    with open("../complementaryData/newjson/%s.json" % organism, "w") as outfile:
        json.dump(newdata, outfile, indent=4)


def main() :
    # pombe("Schizosaccharomyces_pombe")
    # saccharomyces("Saccharomyces_cerevisiae")
    albicans("Candida_albicans")


if __name__ == "__main__" :
    main()
