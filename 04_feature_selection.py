#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse
import sys, os, re
import numpy as np
from pubscripts import *
# from clusters import *
from dimreduction import lda as dimlda
from dimreduction import pca as dimpca
from dimreduction import tsne as dimtsne
from featurenormalization import *
from featureselection import *
# from machinelearning import *
from itertools import combinations

if __name__ == '__main__':
    parameters = {
        'Output_Format': "tsv",
        'Feature_Normalization_Algorithm': "ZScore",
        'Feature_Selection_Algorithm': "CHI2",
        'Dimension_Reduction_Algorithm': "pca",
        'Dimension_Reduction_Number': 3,
        'Selected_Feature_Number': 100
    }

    # prepare data for feature normalization, selection and dimension reduction
    training_code_file = 'training_code.txt'
    testing_code_file = 'testing_code.txt'
    training_code_used, training_labels = read_code.read_code(training_code_file, format=parameters['Output_Format'])
    testing_code_used, testing_labels = [], []
    # if len(testing_data) > 0:
    testing_code_used, testing_labels = read_code.read_code(testing_code_file, format=parameters['Output_Format'])

    # feature normalization
    training_code_normalized = []
    testing_code_normalized = []
    if parameters['Feature_Normalization_Algorithm'] != '' and parameters['Feature_Normalization_Algorithm'] in ('ZScore, MinMax'):
        cmd = parameters['Feature_Normalization_Algorithm'] + '.' + parameters[
            'Feature_Normalization_Algorithm'] + '(training_code_used, training_labels)'
        training_code_normalized, e = eval(cmd)
        save_file.save_file(training_code_normalized, format=parameters['Output_Format'],
                            file='training_code_normalized.txt')
        training_code_file = 'training_code_normalized.txt'

        # if len(testing_data) > 0:
        cmd = parameters['Feature_Normalization_Algorithm'] + '.' + parameters[
            'Feature_Normalization_Algorithm'] + '(testing_code_used, testing_labels)'
        testing_code_normalized, e = eval(cmd)
        save_file.save_file(testing_code_normalized, format=parameters['Output_Format'],
                            file='testing_code_normalized.txt')
        testing_code_file = 'testing_code_normalized.txt'

    # feature selection
    training_code_selected = []
    testing_code_selected = []
    if parameters['Feature_Selection_Algorithm'] != '' and parameters['Feature_Selection_Algorithm'] in (
    'CHI2', 'IG', 'MIC', 'pearsonr', 'Fscore'):
        training_code_used, training_labels = read_code.read_code(training_code_file,
                                                                  format=parameters['Output_Format'])
        # if len(testing_data) > 0:
        testing_code_used, testing_labels = read_code.read_code(testing_code_file,
                                                                format=parameters['Output_Format'])

        cmd = parameters['Feature_Selection_Algorithm'] + '.' + parameters[
            'Feature_Selection_Algorithm'] + '(training_code_used, training_labels)'
        selected_features, e = eval(cmd)
        save_file.save_FS_result(selected_features, e, parameters['Feature_Selection_Algorithm'], 'feaure_rank.txt')
        training_code_selected = select_features.select_features(training_code_used, training_labels,
                                                                 selected_features,
                                                                 parameters['Selected_Feature_Number']).tolist()
        save_file.save_file(training_code_selected, parameters['Output_Format'], 'training_code_selected.txt')
        training_code_file = 'training_code_selected.txt'

        # if len(testing_data) > 0:
        testing_code_selected = select_features.select_features(testing_code_used, testing_labels,
                                                                selected_features,
                                                                parameters['Selected_Feature_Number']).tolist()
        save_file.save_file(testing_code_selected, parameters['Output_Format'], 'testing_code_selected.txt')
        testing_code_file = 'testing_code_selected.txt'

    training_code_file = 'training_code_selected.txt'
    testing_code_file = 'testing_code_selected.txt'    