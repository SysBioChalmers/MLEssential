#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-12-14

# This python script is to plot boxplopt for gene essentiality data and number of yeast species according to their othologs


import json
import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt


organisms = ["Candida_albicans", "Candida_glabrata", "Candida_dubliniensis", "Candida_parapsilosis", "Candida_tropicalis", "Yarrowia_lipolytica", "Schizosaccharomyces_pombe", "Saccharomyces_cerevisiae"]


for organism in organisms :
    print("This organism is: %s" % organism.replace("_", " "))

    allProteinId = set()
    with open("../reciprocal_blast/%s_ref.fasta" % organism, "rU") as allProtein :
        for record in SeqIO.parse(allProtein, "fasta") :
            allProteinId.add(record.id)  # allProteinId is all id in fasta file

    with open("../processed_data/%s_include_ortholog.json" % organism, "r") as f :
        data = json.load(f)

    complexData2 = [complexData["id"] for complexData in data if list(complexData.values())[0] == "Complex"]
    complexData = set(complexData2)

    # complexData = [float(complexData["species"]) for complexData in data if list(complexData.values())[0] == "Complex"]
    # non_complexData = [float(complexData["species"]) for complexData in data if list(complexData.values())[0] == "Non-complex"]

    allId2 = [complexData["id"] for complexData in data]
    allId = set(allId2) # allId is all the id through reciprocal blast best hits
    non_complexData = allId - complexData

    print("The number of all proteins for %s: %d" % (organism, len(allProteinId)))
    print("The number of proteins through reciprocal blast mapping for %s: %d" % (organism, len(allId))) 
    # print("The number of all proteins: %d" % len(complexData))
    # print("The number of complex proteins without duplication: %d" % len(complexData2))
    print("The number of complex proteins for %s: %d" % (organism, len(complexData)))
    print("The number of non-complex proteins: %d" % len(non_complexData))

    # https://blog.csdn.net/aijiudu/article/details/89387328
    # print(pd.DataFrame(complexData).describe())
    # print(pd.DataFrame(non_complexData).describe())

    complexSpecies = list()
    for seqid in complexData :
        for onedict in data :
            if seqid == onedict["id"] :
                complexSpecies.append(float(onedict["species"]))
                break

    non_complexSpecies  = list()
    for seqid in non_complexData :
        for onedict in data :
            if seqid == onedict["id"] :
                non_complexSpecies.append(float(onedict["species"]))
                break

    # print(len(complexSpecies))
    # print(len(non_complexSpecies))

    print(pd.DataFrame(complexSpecies).describe())
    print(pd.DataFrame(non_complexSpecies).describe())

    print("-------------------------------------")

    labels = ["Complex", "Non-complex"]

    figure,axes=plt.subplots(figsize=(6,4))

    # rectangular box plot 
    boxplot = axes.boxplot([np.array(complexSpecies).astype(np.float),np.array(non_complexSpecies).astype(np.float)],
                              vert=True, # vertical box alignment
                              patch_artist=True, # fill with color
                              labels=labels) # will be used to label x-ticks
    axes.set_title(organism.replace("_", " "))

    # fill with colors
    colors = ['lightblue', 'lightgreen']

    for patch, color in zip(boxplot['boxes'], colors):
        patch.set_facecolor(color)

    # adding horizontal grid lines
    axes.yaxis.grid(True)

    axes.set_ylabel("Number of gene occurance")

    plt.savefig("../figure/boxplot/%s_boxplot.png" % organism, dpi=400)


# Results :

# This organism is: Candida albicans
# The number of all proteins for Candida_albicans: 6207
# The number of proteins through reciprocal blast mapping for Candida_albicans: 6184
# The number of complex proteins for Candida_albicans: 1880
# The number of non-complex proteins: 4304
#                  0
# count  1880.000000
# mean    285.186170
# std      80.397631
# min       1.000000
# 25%     273.000000
# 50%     319.000000
# 75%     338.000000
# max     343.000000
#                  0
# count  4304.000000
# mean    197.522537
# std     129.726085
# min       1.000000
# 25%      73.000000
# 50%     246.000000
# 75%     320.000000
# max     343.000000
# -------------------------------------
# This organism is: Candida glabrata
# The number of all proteins for Candida_glabrata: 5143
# The number of proteins through reciprocal blast mapping for Candida_glabrata: 4898
# The number of complex proteins for Candida_glabrata: 1671
# The number of non-complex proteins: 3227
#                  0
# count  1671.000000
# mean    270.296828
# std     103.813175
# min       1.000000
# 25%     269.500000
# 50%     319.000000
# 75%     338.000000
# max     343.000000
#                  0
# count  3227.000000
# mean    233.545398
# std     124.880449
# min       1.000000
# 25%      93.500000
# 50%     301.000000
# 75%     331.000000
# max     343.000000
# -------------------------------------
# This organism is: Candida dubliniensis
# The number of all proteins for Candida_dubliniensis: 5949
# The number of proteins through reciprocal blast mapping for Candida_dubliniensis: 5374
# The number of complex proteins for Candida_dubliniensis: 1629
# The number of non-complex proteins: 3745
#                  0
# count  1629.000000
# mean    294.149171
# std      72.657610
# min       2.000000
# 25%     287.000000
# 50%     323.000000
# 75%     339.000000
# max     343.000000
#                  0
# count  3745.000000
# mean    216.018158
# std     121.798836
# min       1.000000
# 25%      91.000000
# 50%     268.000000
# 75%     325.000000
# max     343.000000
# -------------------------------------
# This organism is: Candida parapsilosis
# The number of all proteins for Candida_parapsilosis: 5280
# The number of proteins through reciprocal blast mapping for Candida_parapsilosis: 4994
# The number of complex proteins for Candida_parapsilosis: 1571
# The number of non-complex proteins: 3423
#                  0
# count  1571.000000
# mean    294.688733
# std      73.967090
# min       1.000000
# 25%     288.500000
# 50%     325.000000
# 75%     339.000000
# max     343.000000
#                  0
# count  3423.000000
# mean    219.730646
# std     122.065204
# min       1.000000
# 25%      93.000000
# 50%     276.000000
# 75%     327.000000
# max     343.000000
# -------------------------------------
# This organism is: Candida tropicalis
# The number of all proteins for Candida_tropicalis: 5975
# The number of proteins through reciprocal blast mapping for Candida_tropicalis: 5437
# The number of complex proteins for Candida_tropicalis: 466
# The number of non-complex proteins: 4971
#                 0
# count  466.000000
# mean   294.362661
# std     75.245715
# min      1.000000
# 25%    285.000000
# 50%    324.000000
# 75%    340.000000
# max    343.000000
#                  0
# count  4971.000000
# mean    229.742909
# std     120.481734
# min       1.000000
# 25%     114.000000
# 50%     288.000000
# 75%     330.000000
# max     343.000000
# -------------------------------------
# This organism is: Yarrowia lipolytica
# The number of all proteins for Yarrowia_lipolytica: 6433
# The number of proteins through reciprocal blast mapping for Yarrowia_lipolytica: 6009
# The number of complex proteins for Yarrowia_lipolytica: 1236
# The number of non-complex proteins: 4773
#                  0
# count  1236.000000
# mean    275.579288
# std      99.997827
# min       1.000000
# 25%     274.000000
# 50%     319.000000
# 75%     337.000000
# max     343.000000
#                  0
# count  4773.000000
# mean    188.885607
# std     142.085271
# min       1.000000
# 25%       6.000000
# 50%     251.000000
# 75%     324.000000
# max     343.000000
# -------------------------------------
# This organism is: Schizosaccharomyces pombe
# The number of all proteins for Schizosaccharomyces_pombe: 5134
# The number of proteins through reciprocal blast mapping for Schizosaccharomyces_pombe: 5016
# The number of complex proteins for Schizosaccharomyces_pombe: 1585
# The number of non-complex proteins: 3431
#                  0
# count  1585.000000
# mean    246.906625
# std     127.910979
# min       1.000000
# 25%     210.000000
# 50%     315.000000
# 75%     336.000000
# max     343.000000
#                  0
# count  3431.000000
# mean    196.136695
# std     147.264894
# min       1.000000
# 25%       5.000000
# 50%     280.000000
# 75%     330.000000
# max     343.000000
# -------------------------------------
# This organism is: Saccharomyces cerevisiae
# The number of all proteins for Saccharomyces_cerevisiae: 5911
# The number of proteins through reciprocal blast mapping for Saccharomyces_cerevisiae: 5730
# The number of complex proteins for Saccharomyces_cerevisiae: 2011
# The number of non-complex proteins: 3719
#                  0
# count  2011.000000
# mean    266.017404
# std     104.363445
# min       1.000000
# 25%     255.000000
# 50%     316.000000
# 75%     337.000000
# max     343.000000
#                  0
# count  3719.000000
# mean    212.595052
# std     133.939336
# min       1.000000
# 25%      62.000000
# 50%     286.000000
# 75%     328.000000
# max     343.000000
# -------------------------------------
# [Finished in 19.7s]
