#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse
import sys, os, re
import numpy as np
from pubscripts import *
from clusters import *
from dimreduction import lda as dimlda
from dimreduction import pca as dimpca
from dimreduction import tsne as dimtsne
from itertools import combinations
from descnucleotide import *

if __name__ == '__main__':
    parameters = {
        'Sequence_Type': 'DNA',
        'Sequence_File': 'traintest.txt',
        'Kmer_Size': 3,
        'Sliding_Window':5,
        'K_Space': 5,
        'Lag_Value': 3,
        'AAindex': "ANDN920101;ARGP820101;ARGP820102;ARGP820103;BEGF750101;BEGF750102;BEGF750103;BHAR880101",
        # 'Method': "DNDS;Conservation;Occurance;ProteinNumber;DNC;Kmer;TACC;TCC;TAC;PseKNC;AAC;CTDC;SOCNumber;GAAC",
        'Method': "DNDS;Conservation;Occurance;ProteinNumber;DNC;Kmer",
        # 'Method': "DNC;Kmer",
        'Output_Format': "tsv",
        'Clustering_Type': "sample",
        'Kmean_Cluster_Number': 2,
        'Clustering_Algorithm': "kmeans",
        'All_Property': False,
        'Tri-DNA-Phychem': "Dnase I;Bendability (DNAse)",
        'Weight_Value': 0.1,
        'Lamada_Value': 2,
        'Di-DNA-Phychem-default6': 'Rise;Roll;Shift;Slide;Tilt;Twist'
    }

    # commands for encoding
    dna_cmd_coding = {
        'Kmer': ['Kmer.Kmer(training_data, k=%s, **kw)' % parameters['Kmer_Size'], 'Kmer.Kmer(testing_data, k=%s, **kw)' % parameters['Kmer_Size']],
        'RCKmer': ['RCKmer.RCKmer(training_data, k=%s, **kw)' % parameters['Kmer_Size'], 'RCKmer.RCKmer(testing_data, k=%s, **kw)' % parameters['Kmer_Size']],
        'NAC': ['NAC.NAC(training_data, **kw)', 'NAC.NAC(testing_data, **kw)'],
        'DNC': ['DNC.DNC(training_data, **kw)', 'DNC.DNC(testing_data, **kw)'],
        'TNC': ['TNC.TNC(training_data, **kw)', 'TNC.TNC(testing_data, **kw)'],
        'ANF': ['ANF.ANF(training_data, **kw)', 'ANF.ANF(testing_data, **kw)'],
        'ENAC': ['ENAC.ENAC(training_data, window=%s, **kw)' % parameters['Sliding_Window'], 'ENAC.ENAC(testing_data, window=%s, **kw)' % parameters['Sliding_Window']],
        'binary': ['binary.binary(training_data, **kw)', 'binary.binary(testing_data, **kw)'],
        'CKSNAP': ['CKSNAP.CKSNAP(training_data, gap=%s, **kw)' % parameters['K_Space'], 'CKSNAP.CKSNAP(testing_data, gap=%s, **kw)' % parameters['K_Space']],
        'NCP': ['NCP.NCP(training_data, **kw)', 'NCP.NCP(testing_data, **kw)'],
        'PSTNPss': ['PSTNPss.PSTNPss(PSTNP_training_data, **kw)', 'PSTNPss.PSTNPss(PSTNP_testing_data, **kw)'],
        'PSTNPds': ['PSTNPds.PSTNPds(PSTNP_training_data, **kw)', 'PSTNPds.PSTNPds(PSTNP_testing_data, **kw)'],
        'EIIP': ['EIIP.EIIP(training_data, **kw)', 'EIIP.EIIP(testing_data, **kw)'],
        'PseEIIP': ['PseEIIP.PseEIIP(training_data, **kw)', 'PseEIIP.PseEIIP(testing_data, **kw)'],
        'DAC': ['ACC.make_ac_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value'], 'ACC.make_ac_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value']],
        'DCC': ['ACC.make_cc_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value'], 'ACC.make_cc_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value']],
        'DACC': ['ACC.make_acc_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value'], 'ACC.make_acc_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value']],
        'TAC': ['ACC.make_ac_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value'], 'ACC.make_ac_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value']],
        'TCC': ['ACC.make_cc_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value'], 'ACC.make_cc_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value']],
        'TACC': ['ACC.make_acc_vector(training_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value'], 'ACC.make_acc_vector(testing_data, my_property_name, my_property_value, %s, my_kmer)' % parameters['Lag_Value']],
        'PseDNC': ['Pse.make_PseDNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight)','Pse.make_PseDNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight)'],
        'PseKNC': ['Pse.make_PseKNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight, %d)' % int(parameters['Kmer_Size']), 'Pse.make_PseKNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight, %d)' % int(parameters['Kmer_Size'])],
        'PCPseDNC': ['Pse.make_PseDNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight)','Pse.make_PseDNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight)'],
        'PCPseTNC': ['Pse.make_PCPseTNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight)','Pse.make_PCPseTNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight)'],
        'SCPseDNC': ['Pse.make_SCPseDNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight)','Pse.make_SCPseDNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight)'],
        'SCPseTNC': ['Pse.make_SCPseTNC_vector(training_data, my_property_name, my_property_value, my_lamada, my_weight)', 'Pse.make_SCPseTNC_vector(testing_data, my_property_name, my_property_value, my_lamada, my_weight)'],
        'DNDS': ['DNDS.dnds(training_data, **kw)', 'DNDS.dnds(testing_data, **kw)'],
        'Conservation': ['Conservation.conservation_score(training_data, **kw)', 'Conservation.conservation_score(testing_data, **kw)'],
        'Occurance': ['Occurance.occurance_number(training_data, **kw)', 'Occurance.occurance_number(testing_data, **kw)'],
        'ProteinNumber': ['ProteinNumber.protein_number(training_data, **kw)', 'ProteinNumber.protein_number(testing_data, **kw)'],
        'AAC': ['AAC.AAC(training_data, **kw)', 'AAC.AAC(testing_data, **kw)'],
        'CTDC': ['CTDC.CTDC(training_data, **kw)', 'CTDC.CTDC(testing_data, **kw)'],
        'SOCNumber': ['SOCNumber.SOCNumber(training_data, nlag=%d, **kw)' % int(parameters['Lag_Value']), 'SOCNumber.SOCNumber(testing_data, nlag=%d, **kw)' % int(parameters['Lag_Value'])],
        'PAAC': ['PAAC.PAAC(training_data, lambdaValue=%d, w=%f, **kw)' % (int(parameters['Lamada_Value']), float(parameters['Weight_Value'])), 'PAAC.PAAC(testing_data, lambdaValue=%d, w=%f, **kw)' % (int(parameters['Lamada_Value']), float(parameters['Weight_Value']))],
        'GAAC': ['GAAC.GAAC(training_data, **kw)', 'GAAC.GAAC(testing_data, **kw)'],
        'GDPC': ['GDPC.GDPC(training_data, **kw)', 'GDPC.GDPC(testing_data, **kw)'],
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
    # kw['order'] = 'ACGT' if parameters['Sequence_Type'] == 'DNA' or parameters['Sequence_Type'] == 'RNA' else 'ACDEFGHIKLMNPQRSTVWY'

    # divide training and testing data
    training_data = []
    testing_data = []
    PSTNP_training_data = []
    PSTNP_testing_data = []
    for sequence in fastas:
        if sequence[3] == 'training':
            training_data.append(sequence)
            PSTNP_training_data.append(sequence)
            PSTNP_training_data.append([sequence[0], sequence[1], sequence[2], 'testing'])
            PSTNP_testing_data.append(sequence)
        else:
            testing_data.append(sequence)
            PSTNP_testing_data.append(sequence)
    # print(training_data[:2])  # a list = [genename, geneseq, label(0 or 1), "training"]
    # get property for AAindex, NMBroto, Geary, Moran
    props = parameters['AAindex'].split(';') if parameters['AAindex'] != '' else ['CIDH920105', 'BHAR880101',
                                                                                  'CHOC760101', 'BIGC670101',
                                                                                  'CHAM810101', 'DAYM780201']


    # get property for ACC descriptors and Pse descriptors
    my_property_name, my_property_value, my_kmer, my_lamada, my_weight = 0, 0, 0, 0, 0

    # calculate descriptor for training data
    training_code_dict = {}
    testing_code_dict = {}
    method_array = parameters['Method'].split(';')
    for method in method_array :
        # before extract these 12 descriptor, the property need be extracted
        if method in ('DAC', 'DCC', 'DACC', 'TAC', 'TCC', 'TACC'):
            my_property_name, my_property_value, my_kmer = check_parameters.check_acc_arguments_pipeline(parameters,
                                                                                                         method)
        if method in ('PseDNC', 'PseKNC', 'PCPseDNC', 'PCPseTNC', 'SCPseDNC', 'SCPseTNC'):
            my_property_name, my_property_value, my_lamada, my_weight = check_parameters.check_Pse_arguments_pipeline(
                parameters, method, fastas)

        if method in ('DNC', 'Kmer', 'TACC', 'TCC', 'TAC', 'PseKNC'):
            kw['order'] = 'ACGT'

        if method in ('AAC', 'CTDC', 'SOCNumber', 'PAAC', 'GAAC', 'GDPC'):
            kw['order'] = 'ACDEFGHIKLMNPQRSTVWY'

        # calculate descriptors
        training_code_dict[method] = eval(cmd_coding[method][0])
        if len(testing_data) > 0:
            testing_code_dict[method] = eval(cmd_coding[method][1])

    training_code = np.array(training_code_dict[method_array[0]])
    testing_code = []
    if len(testing_data) > 0:
        testing_code = np.array(testing_code_dict[method_array[0]])
    # print("----------------------------------")
    # print(testing_code[:3])
    # ----------------------------------
    # [['#' 'label' 'AA' 'AC' 'AG' 'AT' 'CA' 'CC' 'CG' 'CT' 'GA' 'GC' 'GG' 'GT'
    #   'TA' 'TC' 'TG' 'TT']
    #  ['SPAC22A12.15c' '1' '0.09492717227523857' '0.04922149673530889'
    #   '0.06328478151682572' '0.06680060271220492' '0.04018081366147665'
    #   '0.04369663485685585' '0.03867403314917127' '0.07031642390758412'
    #   '0.07734806629834254' '0.04018081366147665' '0.045203415369161226'
    #   '0.06629834254143646' '0.06177800100452034' '0.05976896032144651'
    #   '0.08186840783525866' '0.10045203415369161']
    #  ['SPAC3H8.06' '1' '0.055993690851735015' '0.04810725552050473'
    #   '0.02444794952681388' '0.07807570977917981' '0.057570977917981075'
    #   '0.062302839116719244' '0.031545741324921134' '0.08753943217665615'
    #   '0.033123028391167195' '0.0528391167192429' '0.050473186119873815'
    #   '0.056782334384858045' '0.05993690851735016' '0.07570977917981073'
    #   '0.08675078864353312' '0.138801261829653']]

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

    # print(training_code[:3])
    # print(testing_code)

    training_code = training_code.tolist()
    save_file.save_file(training_code, format=parameters['Output_Format'], file='training_code.txt')
    save_file.save_file(training_code, format='tsv_1', file='training_code_1.tsv')

    if len(testing_data) > 0:
        testing_code = testing_code.tolist()
        save_file.save_file(testing_code, format=parameters['Output_Format'], file='testing_code.txt')
        save_file.save_file(testing_code, format='tsv_1', file='testing_code_1.tsv')

