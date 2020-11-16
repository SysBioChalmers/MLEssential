#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN

# http://cmdlinetips.com/2019/03/how-to-make-grouped-boxplots-in-python-with-seaborn/
# https://github.com/cdanielmachado/cooccurrence/blob/master/notebooks/Figure%205.ipynb


import json
import csv
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ranksums


with open("../complementaryData/evolutionary_data/ortholog_occurence.tsv", "r") as outfile :
    # paralog_data = outfile.readlines()
    paralog_data = csv.reader(outfile,delimiter='\t')
    # print(paralog_data)
    paralog = dict()

    for line in list(paralog_data) :
        ortholog = line[1]
        paralog_number = line[-1]
        # print(type(paralog_number))  # <class 'str'>
        paralog[ortholog] = float(paralog_number)


organisms = ["S_cerevisiae", "S_pombe",  "C_albicans", "Y_lipolytica", "P_pastoris"]

# with open("../Data/complexdata.csv", "w") as outfile :
# outfile = open("../Data/complexdata/complexdata.csv", "w")
# csv_writer = csv.writer(outfile)
# csv_writer.writerow(["type", "organism", "species"])

for organism in organisms :
    # allEssential = set()
    print("This organism is: %s" % organism.replace("_", " "))
    with open("../complementaryData/processed_data/%s_include_ortholog.json" % organism, "r") as f :
        data = json.load(f)


    essential = list()
    non_essential = list()
    for subdata in data :
        ortholog = subdata['ortholog']
        if list(subdata.values())[0] == "E" :
            try :
                essential.append(float(paralog[ortholog]))
            except :
                pass
        if list(subdata.values())[0] == "NE" :
            try :
                non_essential.append(float(paralog[ortholog]))
            except :
                pass

    print("The number of essential genes: %d" % len(essential))
    print("The number of non-essential genes: %d" % len(non_essential))

    print(ranksums(essential,non_essential))


    # # # https://blog.csdn.net/aijiudu/article/details/89387328
    # print(pd.DataFrame(essential).describe())
    # print(pd.DataFrame(non_essential).describe())
    # print("-------------------------------------")


# Results:

# This organism is: S cerevisiae
# The number of essential genes: 1035
# The number of non-essential genes: 4503
# RanksumsResult(statistic=-7.367762667242304, pvalue=1.7351532938677472e-13)
# This organism is: S pombe
# The number of essential genes: 1224
# The number of non-essential genes: 3342
# RanksumsResult(statistic=-3.141419713388134, pvalue=0.0016813089876503247)
# This organism is: C albicans
# The number of essential genes: 625
# The number of non-essential genes: 1682
# RanksumsResult(statistic=-8.861821637073108, pvalue=7.871746884519488e-19)
# This organism is: Y lipolytica
# The number of essential genes: 107
# The number of non-essential genes: 520
# RanksumsResult(statistic=-2.328869232302724, pvalue=0.01986599433739998)
# This organism is: P pastoris
# The number of essential genes: 130
# The number of non-essential genes: 450
# RanksumsResult(statistic=-3.500664623901917, pvalue=0.00046409949475487616)
# [Finished in 1.9s]
