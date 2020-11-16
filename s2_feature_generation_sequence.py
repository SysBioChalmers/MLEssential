#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse
import sys, os, re
import numpy as np
from pubscripts import *
from itertools import combinations
from descnucleotide import *

if __name__ == '__main__':

    parameters = {
        'Sequence_Type': 'DNA',
        'Sequence_File': 'traintest.txt',
        'Kmer_Size': 3,
        'Method': "DNC;Kmer",
        'Output_Format': "tsv",
    }

    # commands for encoding
    dna_cmd_coding = {
        'Kmer': ['Kmer.Kmer(training_data, k=%s, **kw)' % parameters['Kmer_Size'], 'Kmer.Kmer(testing_data, k=%s, **kw)' % parameters['Kmer_Size']],
        'DNC': ['DNC.DNC(training_data, **kw)', 'DNC.DNC(testing_data, **kw)'],
    }

    # Error information
    error_array = []

    # read fasta sequence and specify cmd
    fastas = []
    cmd_coding = {}
    if parameters['Sequence_Type'] in ('DNA', 'RNA'):
        fastas = read_fasta_sequences.read_nucleotide_sequences(parameters['Sequence_File'])
        cmd_coding = dna_cmd_coding
    elif parameters['Sequence_Type'] == 'Protein':
        fastas = read_fasta_sequences.read_protein_sequences(parameters['Sequence_File'])
        cmd_coding = protein_cmd_coding
    else:
        error_array.append('Sequence type can only be selected in "DNA", "RNA" or "Protein".')

    kw = {'nclusters': 3, 'sof': 'sample', 'order': ''}
    kw['order'] = 'ACGT' if parameters['Sequence_Type'] == 'DNA' or parameters['Sequence_Type'] == 'RNA' else 'ACDEFGHIKLMNPQRSTVWY'

    # divide training and testing data
    training_data = []
    testing_data = []
    for sequence in fastas:
        if sequence[3] == 'training':
            training_data.append(sequence)
        else:
            testing_data.append(sequence)
    # print(training_data[:2])  # a list = [genename, geneseq, label(0 or 1), "training"]

    # calculate descriptor for training data
    training_code_dict = {}
    testing_code_dict = {}
    method_array = parameters['Method'].split(';')
    for method in method_array:
        # calculate descriptors
        training_code_dict[method] = eval(cmd_coding[method][0])
        if len(testing_data) > 0:
            testing_code_dict[method] = eval(cmd_coding[method][1])

    training_code = np.array(training_code_dict[method_array[0]])
    testing_code = []
    if len(testing_data) > 0:
        testing_code = np.array(testing_code_dict[method_array[0]])
    print("----------------------------------")
    # print(testing_code[:3])
    # print(testing_code_dict)
    for i in range(1, len(method_array)):
        if training_code_dict[method_array[i]] != 0:
            training_code = np.concatenate((training_code, np.array(training_code_dict[method_array[i]])[:, 2:]), axis=1)
            if len(testing_data) > 0:
                if testing_code_dict[method_array[i]] != 0:
                    testing_code = np.concatenate((testing_code, np.array(testing_code_dict[method_array[i]])[:, 2:]), axis=1)

    if len(testing_data) != 0 and training_code.shape[1] != testing_code.shape[1]:
        error_array.append('Descriptor(s) for testing data calculating failed.')
        testing_data = []

    # print(training_code)
    # print(testing_code)

    training_code = training_code.tolist()
    save_file.save_file(training_code, format=parameters['Output_Format'], file='training_code.txt')
    save_file.save_file(training_code, format='tsv_1', file='training_code_1.tsv')

    if len(testing_data) > 0:
        testing_code = testing_code.tolist()
        save_file.save_file(testing_code, format=parameters['Output_Format'], file='testing_code.txt')
        save_file.save_file(testing_code, format='tsv_1', file='testing_code_1.tsv')

