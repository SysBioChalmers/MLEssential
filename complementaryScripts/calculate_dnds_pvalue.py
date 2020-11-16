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


with open("../complementaryData/evolutionary_data/gene_dn_ds.csv", "r") as outfile :
    dnds_data = outfile.readlines()
    # print(dnds_data)
    dnds = dict()

    for line in dnds_data :
        ortholog = line.strip().split(",")[1].split(".")[0]
        dnds_score = line.strip().split(",")[2]
        # print(type(dnds_score))  # <class 'str'>

        if dnds_score :
            dnds[ortholog] = float(dnds_score)

organisms = ["S_cerevisiae", "S_pombe",  "C_albicans", "Y_lipolytica", "P_pastoris"]

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
                essential.append(float(dnds[ortholog]))
            except :
                pass
        if list(subdata.values())[0] == "NE" :
            try :
                non_essential.append(float(dnds[ortholog]))
            except :
                pass

    # print("The number of all collected data: %d" %(len(data)))
    # print("The number of data without duplication: %d" % len(allessential))
    print("The number of essential genes: %d" % len(essential))
    print("The number of non-essential genes: %d" % len(non_essential))

    print(ranksums(essential,non_essential))


# print(essential)
# print(non_essential)


# Results:

# This organism is: S cerevisiae
# The number of essential genes: 980
# The number of non-essential genes: 3612
# RanksumsResult(statistic=-6.859206417906347, pvalue=6.924412471738405e-12)
# This organism is: S pombe
# The number of essential genes: 1065
# The number of non-essential genes: 2216
# RanksumsResult(statistic=-3.9011859936681295, pvalue=9.57225679621583e-05)
# This organism is: C albicans
# The number of essential genes: 577
# The number of non-essential genes: 1379
# RanksumsResult(statistic=-6.467017705352924, pvalue=9.995594198009325e-11)
# This organism is: Y lipolytica
# The number of essential genes: 96
# The number of non-essential genes: 418
# RanksumsResult(statistic=-3.2747574321426938, pvalue=0.0010575273601970865)
# This organism is: P pastoris
# The number of essential genes: 118
# The number of non-essential genes: 374
# RanksumsResult(statistic=-2.1054426774365163, pvalue=0.03525279214434212)
# [Finished in 1.5s]
