#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse
import sys, os, re
import numpy as np
from pubscripts import *
from machinelearning import *
from itertools import combinations

if __name__ == '__main__':
    parameters = {
        'Sequence_Type': 'DNA',
        'Sequence_File': 'traintest.txt',
        'Output_Format': "tsv",
        'ML': "SVM;RF",
        'Validation': 5,
        'Kernel': "rbf",
        'Auto_Opterimize': "True",
        'Cost': 1.0,
        'Gamma': '',
        'Tree_Number': 100,
    }

    # Error information
    error_array = []

    # read fasta sequence and specify cmd
    fastas = []
    cmd_coding = {}
    if parameters['Sequence_Type'] in ('DNA', 'RNA'):
        fastas = read_fasta_sequences.read_nucleotide_sequences(parameters['Sequence_File'])
    elif parameters['Sequence_Type'] == 'Protein':
        fastas = read_fasta_sequences.read_protein_sequences(parameters['Sequence_File'])
    else:
        error_array.append('Sequence type can only be selected in "DNA", "RNA" or "Protein".')

    # divide training and testing data
    training_data = []
    testing_data = []

    for sequence in fastas:
        if sequence[3] == 'training':
            training_data.append(sequence)
        else:
            testing_data.append(sequence)
    # print(training_data[2])  # a list = [genename, geneseq, label(0 or 1), "training"]

    # prepare data for machine learning
    # training_code_file = 'training_code_selected.txt'
    # testing_code_file = 'testing_code_selected.txt'
    training_code_file = 'training_code.txt'
    testing_code_file = 'testing_code.txt'

    # machine learning
    ML_array = parameters['ML'].split(';')
    if parameters['ML'] != '' and parameters['Validation'] != '':
        X, y, independent = 0, 0, np.array([])
        X, y = read_code_ml.read_code(training_code_file, format=parameters['Output_Format'])
        classes = sorted(list(set(y)))
        if len(testing_data) > 0:
            ind_X, ind_y = read_code_ml.read_code(testing_code_file, format=parameters['Output_Format'])
            independent = np.zeros((ind_X.shape[0], ind_X.shape[1] + 1))
            independent[:, 0], independent[:, 1:] = ind_y, ind_X

        fold = parameters['Validation'] if parameters['Validation'] != '' else 5
        # svm
        kernel = parameters['Kernel'] if parameters['Kernel'] != '' else 'rbf'
        auto = False
        cost = 1.0
        gamma = 1 / len(training_data)
        if parameters['Auto_Opterimize'] == 'True':
            auto = True
            cost = 1.0
            gamma = 1 / len(training_data)
        else:
            cost = float(parameters['Cost']) if parameters['Cost'] != '' else 1.0
            gamma = float(parameters['Gamma']) if parameters['Gamma'] != '' else 1 / len(training_data)

        # RF
        n_trees = parameters['Tree_Number'] if parameters['Tree_Number'] != '' else 100

        ensemble_train_data = {}
        ensemble_test_data = {}
        # training model
        if len(classes) == 2:
            ML_AUCs_dict = {}
            para_info_dict = {}
            for ML in ML_array:
                para_info, cv_res, ind_res = 0, 0, 0
                if ML == 'SVM':
                    para_info, cv_res, ind_res = SVM.SVM_Classifier(X, y, indep=independent, fold=fold, batch=0.8,
                                                                    auto=auto, kernel=kernel, gamma=gamma, C=cost)
                    con = np.concatenate(tuple((f for f in cv_res)), axis=0)
                    ensemble_train_data[ML] = con
                    ensemble_test_data[ML] = ind_res
                if ML == 'RF':
                    para_info, cv_res, ind_res = RF.RF_Classifier(X, y, indep=independent, fold=fold, n_trees=n_trees)
                    con = np.concatenate(tuple((f for f in cv_res)), axis=0)
                    ensemble_train_data[ML] = con
                    ensemble_test_data[ML] = ind_res
                para_info_dict[ML] = para_info

                save_file.save_CV_result_binary(cv_res, '%s_CV.txt' % ML, para_info)
                ML_AUCs_dict[ML] = draw_plot.plot_roc_cv(cv_res, '%s_ROC_CV.png' % ML, label_column=0, score_column=2)
                mean_auprc = draw_plot.plot_prc_CV(cv_res, '%s_PRC_CV.png' % ML, label_column=0, score_column=2)
                cv_metrics = calculate_prediction_metrics.calculate_metrics_cv(cv_res, label_column=0, score_column=2)
                save_file.save_prediction_metrics_cv(cv_metrics, '%s_metrics_CV.txt' % ML)

                # if len(testing_data) > 0:
                # print(ind_res)
                # print(type(ind_res))
                save_file.save_IND_result_binary(ind_res, '%s_IND.txt' % ML, para_info)
                ind_auc = draw_plot.plot_roc_ind(ind_res, '%s_ROC_IND.png' % ML, label_column=0, score_column=2)
                ind_auprc = draw_plot.plot_prc_ind(ind_res, '%s_PRC_IND.png' % ML, label_column=0, score_column=2)
                ind_metrics = calculate_prediction_metrics.calculate_metrics(ind_res[:, 0], ind_res[:, 2])
                save_file.save_prediction_metrics_ind(ind_metrics, '%s_metrics_IND.txt' % ML)

        if len(classes) >= 3:
            acc_dict = {}
            para_info_dict = {}
            for ML in ML_array:
                para_info, cv_res, ind_res = 0, 0, 0
                if ML == 'SVM':
                    para_info, cv_res, ind_res = SVM.SVM_Classifier(X, y, indep=independent, fold=fold, batch=0.8,
                                                                    auto=auto, kernel=kernel, gamma=gamma, C=cost)
                    con = np.concatenate(tuple((f for f in cv_res)), axis=0)
                    ensemble_train_data[ML] = con
                    ensemble_test_data[ML] = ind_res
                if ML == 'RF':
                    para_info, cv_res, ind_res = RF.RF_Classifier(X, y, indep=independent, fold=fold, n_trees=n_trees)
                    con = np.concatenate(tuple((f for f in cv_res)), axis=0)
                    ensemble_train_data[ML] = con
                    ensemble_test_data[ML] = ind_res
                para_info_dict[ML] = para_info

                save_file.save_CV_result(cv_res, classes, '%s_CV.txt' % ML, para_info)
                cv_metrics = calculate_prediction_metrics.calculate_metrics_cv_muti(cv_res, classes, label_column=0)
                save_file.save_prediction_metrics_cv_muti(cv_metrics, classes, '%s_metrics_CV.txt' % ML)

                # calculate mean acc
                mean_acc = 0
                for res in cv_metrics:
                    for c in classes:
                        mean_acc += res[c]
                acc_dict[ML] = mean_acc

                if len(testing_data) > 0:
                    save_file.save_IND_result(ind_res, classes, '%s_IND.txt' % ML, para_info)
                    ind_metrics = calculate_prediction_metrics.calculate_metrics_ind_muti(ind_res, classes,
                                                                                          label_column=0)
                    save_file.save_prediction_metrics_ind_muti(ind_metrics, classes, '%s_metrics_IND.txt' % ML)

