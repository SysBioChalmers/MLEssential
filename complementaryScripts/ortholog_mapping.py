#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-12-19

# This python script is to do the id mapping to collect the original gene id with the ortholog groups


import os
import json
import pandas as pd


def getIndex() :
    # get the ortholog accoding to protein sequence id, that means Alloascoidea_hylecoeti@Seq_1 as the key, 0_0 as the value
    with open("/Users/leyu/Documents/coding/evolution_code/orthomcl_output/orthomcl_SeqIDs_index.txt", "r") as indexFile :
        indexs = indexFile.readlines()

    indexSeqId = dict()
    for index in indexs :
        index_Seq = index.strip().split(": ")
        indexSeqId[index_Seq[1]] = index_Seq[0]

    return indexSeqId 

def getOrthologIndex() :
    with open("/Users/leyu/Documents/coding/evolution_code/orthomcl_output/orthomcl_clusters.txt", "r") as orthologFile :
        orthologs = orthologFile.readlines()

    orthologIndex = dict()
    for ortholog in orthologs :
        ortholog_Index = ortholog.strip().split(" ")
        # orthologIndex = {'OG1001': {'328_2397', '189_1696', '279_256',.....}}
        ortholog = ortholog_Index[0][:-1]
        
        for index in ortholog_Index[1:] :
            orthologIndex[index] = ortholog
    # print(orthologIndex)  # {'302_3224': 'OG1000', '317_1502': 'OG1000', '318_1938': 'OG1001', '320_301': 'OG1001', '325_5347': 'OG1001'}

    return orthologIndex

def getOrthologSpecies() :
    file = os.path.join("../Data/processed_data/", "ortholog_occurence_num_all.tsv")

    with open(file, "r") as f :
        data = f.readlines()

    orthologSpecies = dict()

    for species in data :
        # print(species)
        # ['0', 'OG1000', '', '343', '2281', '6.65014577259475', 'with_sce', '16']
        ortholog_species = species.strip().split("\t")

        orthologSpecies[ortholog_species[1]] = ortholog_species[3]
    # {'OG1000': '343', 'OG1001': '343', 'OG1002': '343'}
    # print(orthologSpecies)

    return orthologSpecies


def getId(organism) :
    with open("../data/processed_data/%s_include_id.json" % organism) as f :
        geneAll = json.load(f)

    indexSeqId = getIndex()
    orthologIndex = getOrthologIndex()
    orthologSpecies = getOrthologSpecies()

    for gene in geneAll :
        try :
            if gene["id"] != None:
                gene["index"] = indexSeqId[gene["id"]] # gene["id"] is the protein id, like Alloascoidea_hylecoeti@Seq_1 indexSeqId["Alloascoidea_hylecoeti@Seq_1"] = 0_0
                gene["ortholog"] = orthologIndex[gene["index"]] # gene["index"] is the index, like 302_3224
                gene["species"] = orthologSpecies[gene["ortholog"]]
        except :
            print(gene)
    print("The number of ortholog mapping with protein id: %d" % len(geneAll))

    return geneAll


def main() :
    # getOrthologSpecies()
    # getOrthologIndex()
    organisms = {"Y_lipolytica", "S_pombe", "S_cerevisiae", "P_pastoris", "C_albicans"}
    for organism in organisms :
        geneAll = getId(organism)
        with open("../Data/processed_data/%s_include_ortholog.json" % organism, "w") as outfile :
            json.dump(geneAll, outfile, indent=4)


if __name__ == "__main__" :
    main()


# Results :

# The number of ortholog mapping with protein id: 627
# The number of ortholog mapping with protein id: 5538
# The number of ortholog mapping with protein id: 2307
# The number of ortholog mapping with protein id: 588
# The number of ortholog mapping with protein id: 4566
# [Finished in 19.6s]

