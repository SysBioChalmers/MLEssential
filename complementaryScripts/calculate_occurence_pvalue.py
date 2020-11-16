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


organisms = ["S_cerevisiae", "S_pombe",  "C_albicans", "Y_lipolytica", "P_pastoris"]

# with open("../Data/essential.csv", "w") as outfile :
# outfile = open("../Data/essential/essential.csv", "w")
# csv_writer = csv.writer(outfile)
# csv_writer.writerow(["type", "organism", "species"])

for organism in organisms :
    # allEssential = set()
    print("This organism is: %s" % organism.replace("_", " "))
    with open("../complementaryData/processed_data/%s_include_ortholog.json" % organism, "r") as f :
        data = json.load(f)


    # for essential in data :
    #     allEssential.add((essential["id"]))
#     # csv_writer.writerow(list(data))
#     for essential in data :
#         csv_writer.writerow([list(essential.values())[0], organism.replace('_', '. '), essential['species']])

# outfile.close()


    essential = [float(essential["species"]) for essential in data if list(essential.values())[0] == "E"]
    non_essential = [float(essential["species"]) for essential in data if list(essential.values())[0] == "NE"]

    # print("The number of all collected data: %d" %(len(data)))
    # print("The number of data without duplication: %d" % len(allessential))
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
# RanksumsResult(statistic=19.78461247730286, pvalue=4.0399131596545685e-87)
# This organism is: S pombe
# The number of essential genes: 1224
# The number of non-essential genes: 3342
# RanksumsResult(statistic=18.611295148994568, pvalue=2.602498980821829e-77)
# This organism is: C albicans
# The number of essential genes: 625
# The number of non-essential genes: 1682
# RanksumsResult(statistic=16.649720694428616, pvalue=3.0401778176405097e-62)
# This organism is: Y lipolytica
# The number of essential genes: 107
# The number of non-essential genes: 520
# RanksumsResult(statistic=5.617063812687873, pvalue=1.942297055931988e-08)
# This organism is: P pastoris
# The number of essential genes: 130
# The number of non-essential genes: 450
# RanksumsResult(statistic=5.925847117741461, pvalue=3.1069154088624106e-09)
# [Finished in 1.5s]
