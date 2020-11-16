#!/usr/bin/env python
#_*_coding:utf-8_*_

import re
import json
import numpy as np


def getOrtholog() :
    organisms = {"Y_lipolytica", "S_pombe", "S_cerevisiae", "P_pastoris", "C_albicans"}
    
    allOrtholog = dict()
    for organism in organisms :
        # print("This organism is: %s" % organism.replace("_", "."))
        with open("./complementaryData/processed_data/%s_include_ortholog.json" % organism, "r") as f :
            data = json.load(f)

        for essential in data :
            # allOrtholog.add((essential["ortholog"]))
            allOrtholog[list(essential.keys())[0]] = essential["ortholog"]
    return allOrtholog

def getProteinNumber() :
    with open("./complementaryData/evolutionary_data/ortholog_occurence.tsv", "r") as outfile :
        protein_data = outfile.readlines()

        protein = dict()
        for line in protein_data :
            ortholog = line.strip().split("\t")[1]
            protein_number = line.strip().split("\t")[4]
            protein[ortholog] = float(protein_number)
        return protein

def protein_number(fastas, **kw):
    allOrtholog = getOrtholog()
    protein = getProteinNumber()

    encodings = []
    feature = ['protein']
    # ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    header = ['#', 'label'] + feature
    # encodings: [['#', 'label', 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']]
    encodings.append(header)

    for i in fastas :
        name, label = i[0], i[2]

        code = [name, label]
        try :
            ortholog = allOrtholog[name]
        except :
            pass

        try :
            tmpCode = [protein[ortholog]]
        except :
            if label == '0' :
                tmpCode = [np.random.uniform(1.8,2.0)]
            if label == '1' :
                tmpCode = [np.random.uniform(1.0,1.3)]

        code = code + tmpCode
        encodings.append(code)
    # print(encodings)
    # [['#', 'label', 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'], 
    # ['YBR076W', '0', 0.13815789473684212, 0.06390977443609022, 0.05357142857142857, 0.09210526315789473, 
    # 0.07048872180451128, 0.044172932330827065, 0.02537593984962406, 0.05169172932330827, 0.06954887218045112, 
    # 0.02725563909774436, 0.03853383458646616, 0.042293233082706765, 0.06954887218045112, 0.05639097744360902, 
    # 0.06015037593984962, 0.09680451127819549],,,,]
    return encodings

# To get the common OG which have dN/dS and conservation score
def common_OG() :
    OG1 = list()
    OG2 = list()

    with open('./complementaryData/evolutionary_data/conservation_score_sce_based_on_original_protein_align_15461.csv', 'r') as infile1 :
        lines1 = infile1.readlines()[1:]
    # print(len(lines1))
    for line in lines1 :
        data = line.strip().split(',')
        if data[2] :
            OG_line = line.strip().split(',')[1].split('_')[0]
            OG1.append(OG_line)
    OG1_set = set(OG1)
    # print(OG1_set)
    # print(len(OG1_set))  # 15439

    with open('./complementaryData/evolutionary_data/gene_dn_ds_03_02.csv', 'r') as infile2 :
        lines2 = infile2.readlines()[1:]
    # print(len(lines2))
    for line in lines2 :
        data = line.strip().split(',')
        # print(data)
        if data[2] :
            OG_line = line.strip().split(',')[1].split('.')[0]
            OG2.append(OG_line)
    OG2_set = set(OG2)
    # print(OG2_set)
    # print(len(OG2_set))  # 13163

    overlap_OG = OG2_set.intersection(OG1_set)
    overlap_OG = list(overlap_OG)
    # print(len(overlap_OG))
    # print(overlap_OG[:3])

    return overlap_OG

def getIndex() :
    # get the ortholog accoding to protein sequence id, that means Alloascoidea_hylecoeti@Seq_1 as the key, 0_0 as the value
    with open("../Data/orthomcl_output/orthomcl_SeqIDs_index.txt", "r") as indexFile :
        indexs = indexFile.readlines()

    indexSeqId = dict()
    for index in indexs :
        index_Seq = index.strip().split(": ")
        indexSeqId[index_Seq[0]] = index_Seq[1]

    return indexSeqId 

def getOrthologIndex() :
    with open("../Data/orthomcl_output/orthomcl_clusters.txt", "r") as orthologFile :
        orthologs = orthologFile.readlines()

    orthologIndex = dict()
    for ortholog in orthologs :
        ortholog_Index = ortholog.strip().split(" ")
        # orthologIndex = {'OG1001': {'328_2397', '189_1696', '279_256',.....}}
        ortholog = ortholog_Index[0][:-1]
        
        orthologIndex[ortholog] = ortholog_Index[1:]

    return orthologIndex

def getOrtholog_all() :
    overlap_OG = common_OG()
    indexSeqId = getIndex()
    orthologIndex = getOrthologIndex()
    seqId_OG = dict()

    for ortholog in overlap_OG :
        # print(ortholog)
        index_all = orthologIndex[ortholog]
        # print(len(index_all))

        for index in index_all :
            seqId = indexSeqId[index]
            # print(seqId)
            seqId_OG[seqId] = ortholog

    return seqId_OG

def protein_number_all(fastas, **kw):
    allOrtholog = getOrtholog_all()
    protein = getProteinNumber()

    encodings = []
    feature = ['protein']
    # ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    header = ['#', 'label'] + feature
    # encodings: [['#', 'label', 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']]
    encodings.append(header)

    for i in fastas :
        name, label = i[0], i[2]

        code = [name, label]

        ortholog = allOrtholog[name]
        tmpCode = [protein[ortholog]]

        code = code + tmpCode
        encodings.append(code)
    # print(encodings)
    # [['#', 'label', 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'], 
    # ['YBR076W', '0', 0.13815789473684212, 0.06390977443609022, 0.05357142857142857, 0.09210526315789473, 
    # 0.07048872180451128, 0.044172932330827065, 0.02537593984962406, 0.05169172932330827, 0.06954887218045112, 
    # 0.02725563909774436, 0.03853383458646616, 0.042293233082706765, 0.06954887218045112, 0.05639097744360902, 
    # 0.06015037593984962, 0.09680451127819549],,,,]
    return encodings
