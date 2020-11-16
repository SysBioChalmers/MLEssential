#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import json
import pickle
import numpy as np
from Bio import SeqIO
from sklearn import svm
from pubscripts import *
from descnucleotide import *


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
    print(len(OG1_set))  # 15439

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
    print(len(OG2_set))  # 13163

    overlap_OG = OG2_set.intersection(OG1_set)
    overlap_OG = list(overlap_OG)
    print(len(overlap_OG))
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

def get_refSeq() :
    # get the protein sequence accoding to protein sequence id
    with open("/Users/leyu/Documents/Le/Data/orthomcl_output/343taxa_proteins.fasta", "r") as handleGene :
        proteinSeq = dict()
        for record in SeqIO.parse(handleGene, "fasta") :
    # ['__add__', '__bool__', '__class__', '__contains__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__',
    # '__getattribute__', '__getitem__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__iter__', '__le__', '__le___', '__len__', '__lt__', 
    # '__module__', '__ne__', '__new__', '__nonzero__', '__radd__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', 
    # '__subclasshook__', '__weakref__', '_per_letter_annotations', '_seq', '_set_per_letter_annotations', '_set_seq', 'annotations', 'dbxrefs', 'description', 
    # 'features', 'format', 'id', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'translate', 'upper']
            # if record.id.startswith("Candida_albicans") :
            # if record.id == gene :
            proteinSeq[record.id] = str(record.seq)
        # print("The protein number of %s is: %d" % (gene,len(proteinSeq)))

    return proteinSeq

def get_ML_yeast_data() :
    overlap_OG = common_OG()
    indexSeqId = getIndex()
    orthologIndex = getOrthologIndex()
    gene_seq = get_refSeq()
    seqId_OG = dict()

    with open('./complementaryData/species/yeast_species.txt', 'r') as outfile :
        lines = outfile.readlines()[1:]

    all_species = [line.strip() for line in lines]
    in_species = ['Saccharomyces_cerevisiae', 'Schizosaccharomyces_pombe', 'Yarrowia_lipolytica', 'Candida_albicans', 'Komagataella_pastoris']
    yeast_species = list(set(all_species) - set(in_species))
    print(len(yeast_species))  # 343-5 = 338
    # print(yeast_species[-15:])

    for ortholog in overlap_OG :
        print(ortholog)
        index_all = orthologIndex[ortholog]
        # print(len(index_all))

        for index in index_all :
            seqId = indexSeqId[index]
            # print(seqId)
            seqId_OG[seqId] = ortholog

    for yeast in yeast_species :
        outfile = open("./complementaryData/all_yeast_input/%s.txt" % yeast, 'w')
        for seqId, ortholog in seqId_OG.items() :
            if seqId.split('@')[0] == yeast :
                sequence = gene_seq[seqId]
                outfile.write('>%s|%s|model' % (seqId,2))
                outfile.write('\n')
                outfile.write(sequence)
                outfile.write('\n')
        outfile.close()

def load_model() :
    with open('./model.pickle', 'rb') as pickle_file :
        svm_model = pickle.load(pickle_file)
    return svm_model

def predict_gene() :
    svm_model = load_model()

    parameters = {'Method': "DNDS;Conservation;Occurance;ProteinNumber;DNC;Kmer", 'Kmer_Size': 3}
    dna_cmd_coding = {
        'Kmer': 'Kmer.Kmer(model_data, k=%s, **kw)' % parameters['Kmer_Size'],
        'DNC': 'DNC.DNC(model_data, **kw)',
        'DNDS': 'DNDS.dnds_all(model_data, **kw)',
        'Conservation': 'Conservation.conservation_score_all(model_data, **kw)',
        'Occurance': 'Occurance.occurance_number_all(model_data, **kw)',
        'ProteinNumber': 'ProteinNumber.protein_number_all(model_data, **kw)',
    }
    
    all_files = os.listdir('./complementaryData/all_yeast_input')
    all_files = [file for file in all_files if file.endswith('txt')]
    # print(len(all_files))  # 338 = 343-5
    # print(all_files[:3])

    order = 0
    essential_status = {'1':'Essential', '0': 'Non_essential'}
    for file in all_files :
        fastas = []
        cmd_coding = {}
        model_data = []
        model_code_dict = {}
        features = []
        labels = []

        outfile = open('./complementaryData/all_yeast_output/%s' % file, 'w')
        outfile.write('Gene_id\tGene_phenotype\n')
        order += 1
        print('This is', order, '--------------------------------')
        fastas = read_fasta_sequences.read_nucleotide_sequences('./complementaryData/all_yeast_input/%s' % file)
        cmd_coding = dna_cmd_coding

        for sequence in fastas:
            if sequence[3] == 'model':
                model_data.append(sequence)

        kw = {'nclusters': 3, 'sof': 'sample', 'order': ''}
        method_array = parameters['Method'].split(';')
        for method in method_array :
            if method in ('DNC', 'Kmer'):
                kw['order'] = 'ACGT'
            model_code_dict[method] = eval(cmd_coding[method])

        model_code = np.array(model_code_dict[method_array[0]])

        for i in range(1, len(method_array)):
            # print(model_code)
            # print(type(model_code))
            if model_code_dict[method_array[i]] != 0:
                model_code = np.concatenate((model_code, np.array(model_code_dict[method_array[i]])[:, 2:]), axis=1)

        model_code = model_code.tolist()
        # print(model_code[0])
        # print(model_code[1])
        # print(len(model_code))  # The first list is one explanation for the following lists

        for info in model_code[1:] :
            # print(info)
            # print(info[0])
            features = info[2:]
            prediction = svm_model.predict([features]).tolist()
            # print(prediction[0])
            # print(essential_status[prediction[0]])

            outfile.write('%s\t%s\n' % (info[0], essential_status[prediction[0]]))
        outfile.close()


if __name__ == "__main__" :
    # common_OG()
    # get_ML_yeast_data()
    predict_gene()


