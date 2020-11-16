#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc


def plot_roc_ind(data1, data2, out, label_column=0, score_column=2):
    fprIndep1, tprIndep1, thresholdsIndep1 = roc_curve(data1[:, label_column], data1[:, score_column])
    fprIndep2, tprIndep2, thresholdsIndep2 = roc_curve(data2[:, label_column], data2[:, score_column])
    ind_auc1 = auc(fprIndep1, tprIndep1)
    ind_auc2 = auc(fprIndep2, tprIndep2)
    # fig = plt.figure(figsize=(3.0,2.1))
    fig = plt.figure(0)
    plt.plot(fprIndep1, tprIndep1, lw=2, alpha=0.7, color='blueviolet',
             label='Seq&Evo (AUC = %0.2f)' % ind_auc1)
    plt.plot(fprIndep2, tprIndep2, lw=2, alpha=0.7, color='cornflowerblue',
             label='Seq (AUC = %0.2f)' % ind_auc2)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=18)
    plt.ylabel('True Positive Rate', fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(loc="lower right", prop={"size":16})

    plt.savefig(out, dpi=400, bbox_inches='tight')
    plt.close(0)
    # return ind_auc1, ind auc2

def svm_model_comparison() :
    ML = "SVM"
    with open('./results_SVM/seq&Evo_oversampling/SVM_IND.txt', 'r') as outfile :
        lines = outfile.readlines()[1:]
    # print(len(lines))
    # print(lines[0])

    ind_res1 = list()
    for line in lines :
        all_info = list()
        data = line.strip().split('\t')
        # print(data[0])
        # print(type(data[0]))
        # print(data[1])

        all_info.append(int(data[0]))
        all_info.append(1-float(data[1]))
        all_info.append(float(data[1]))
        # print(all_info)
        ind_res1.append(list(all_info))

    ind_res1 = np.array(ind_res1)
    # print(ind_res)

    with open('./results_SVM/seq_oversampling/SVM_IND.txt', 'r') as outfile :
        lines = outfile.readlines()[1:]
    # print(len(lines))
    # print(lines[0])

    ind_res2 = list()
    for line in lines :
        all_info = list()
        data = line.strip().split('\t')
        # print(data[0])
        # print(type(data[0]))
        # print(data[1])

        all_info.append(int(data[0]))
        all_info.append(1-float(data[1]))
        all_info.append(float(data[1]))
        # print(all_info)
        ind_res2.append(list(all_info))

    ind_res2 = np.array(ind_res2)
    # print(ind_res)

    plot_roc_ind(ind_res1, ind_res2, './complementaryData/figure/%s_ROC_comparison.pdf' % ML, label_column=0, score_column=2)

def RF_model_comparison() :
    ML = "RF"
    with open('./results_RF/seq&Evo_oversampling/RF_IND.txt', 'r') as outfile :
        lines = outfile.readlines()[1:]
    # print(len(lines))
    # print(lines[0])

    ind_res1 = list()
    for line in lines :
        all_info = list()
        data = line.strip().split('\t')
        # print(data[0])
        # print(type(data[0]))
        # print(data[1])

        all_info.append(int(data[0]))
        all_info.append(1-float(data[1]))
        all_info.append(float(data[1]))
        # print(all_info)
        ind_res1.append(list(all_info))

    ind_res1 = np.array(ind_res1)
    # print(ind_res)

    with open('./results_RF/seq_oversampling/RF_IND.txt', 'r') as outfile :
        lines = outfile.readlines()[1:]
    # print(len(lines))
    # print(lines[0])

    ind_res2 = list()
    for line in lines :
        all_info = list()
        data = line.strip().split('\t')
        # print(data[0])
        # print(type(data[0]))
        # print(data[1])

        all_info.append(int(data[0]))
        all_info.append(1-float(data[1]))
        all_info.append(float(data[1]))
        # print(all_info)
        ind_res2.append(list(all_info))

    ind_res2 = np.array(ind_res2)
    # print(ind_res)

    plot_roc_ind(ind_res1, ind_res2, './complementaryData/figure/%s_ROC_comparison.pdf' % ML, label_column=0, score_column=2)


if __name__ == "__main__" :
    svm_model_comparison()
    RF_model_comparison()





