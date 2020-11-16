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


with open("../complementaryData/evolutionary_data/conservation_score.txt", "r") as outfile :
    conservation_data = outfile.readlines()

    conservation = dict()
    for line in conservation_data :
        ortholog = line.strip().split(" ")[2][1:-1]
        conservation_score = line.strip().split(" ")[3]
        conservation[ortholog] = float(conservation_score)

    # print(conservation['OG3209'])

# with open("../Data/complexdata.csv", "w") as outfile :
# outfile = open("../Data/complexdata/complexdata.csv", "w")
# csv_writer = csv.writer(outfile)
# csv_writer.writerow(["type", "organism", "species"])

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
                essential.append(float(conservation[ortholog]))
            except :
                pass
        if list(subdata.values())[0] == "NE" :
            try :
                non_essential.append(float(conservation[ortholog]))
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
# The number of essential genes: 1033
# The number of non-essential genes: 4301
# RanksumsResult(statistic=4.84696083548837, pvalue=1.2536716457740872e-06)
# This organism is: S pombe
# The number of essential genes: 1140
# The number of non-essential genes: 2600
# RanksumsResult(statistic=4.304684764366415, pvalue=1.6722367963058152e-05)
# This organism is: C albicans
# The number of essential genes: 625
# The number of non-essential genes: 1665
# RanksumsResult(statistic=10.869464340095213, pvalue=1.6114414777754976e-27)
# This organism is: Y lipolytica
# The number of essential genes: 107
# The number of non-essential genes: 520
# RanksumsResult(statistic=4.8690010107164525, pvalue=1.121638341671566e-06)
# This organism is: P pastoris
# The number of essential genes: 130
# The number of non-essential genes: 450
# RanksumsResult(statistic=3.147419250243309, pvalue=0.0016471859552954343)
# [Finished in 1.5s]
