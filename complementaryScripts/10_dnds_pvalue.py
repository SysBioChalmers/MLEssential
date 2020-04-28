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


with open("../complementaryData/evolution_feature/gene_dn_ds.csv", "r") as outfile :
    dnds_data = outfile.readlines()
    # print(dnds_data)
    dnds = dict()

    for line in dnds_data :
        ortholog = line.strip().split(",")[1].split(".")[0]
        dnds_score = line.strip().split(",")[2]
        # print(type(dnds_score))  # <class 'str'>

        if dnds_score :
            dnds[ortholog] = float(dnds_score)

organisms = ["Candida_glabrata", "Candida_dubliniensis", "Candida_parapsilosis", "Candida_tropicalis", "Candida_albicans", "Yarrowia_lipolytica", "Schizosaccharomyces_pombe", "Saccharomyces_cerevisiae"]

# with open("../Data/complexdata.csv", "w") as outfile :
# outfile = open("../Data/complexdata/complexdata.csv", "w")
# csv_writer = csv.writer(outfile)
# csv_writer.writerow(["type", "organism", "species"])

for organism in organisms :
    # allEssential = set()
    print("This organism is: %s" % organism.replace("_", " "))
    with open("../complementaryData/newjson/%s.json" % organism, "r") as f :
        data = json.load(f)

    complexdata = list()
    non_complexdata = list()
    for subdata in data :
        ortholog = subdata['ortholog']
        if list(subdata.values())[0] == "Complex" :
            try :
                complexdata.append(float(dnds[ortholog]))
            except :
                pass
        if list(subdata.values())[0] == "Non-complex" :
            try :
                non_complexdata.append(float(dnds[ortholog]))
            except :
                pass

    # print("The number of all collected data: %d" %(len(data)))
    # print("The number of data without duplication: %d" % len(allcomplexdata))
    print("The number of complexdata genes: %d" % len(complexdata))
    print("The number of non-complexdata genes: %d" % len(non_complexdata))

    print(ranksums(complexdata,non_complexdata))


    # # # https://blog.csdn.net/aijiudu/article/details/89387328
    # print(pd.DataFrame(complexdata).describe())
    # print(pd.DataFrame(non_complexdata).describe())
    # print("-------------------------------------")

