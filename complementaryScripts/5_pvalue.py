#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-12-14

# This python script is to plot boxplopt for gene essentiality data and number of yeast species according to their othologs


import json
import numpy as np
import pandas as pd
from scipy.stats import ranksums
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

    # print("The number of all proteins for %s: %d" % (organism, len(allProteinId)))
    # print("The number of proteins through reciprocal blast mapping for %s: %d" % (organism, len(allId))) 
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
    print(ranksums(complexSpecies,non_complexSpecies))
    # print(pd.DataFrame(complexSpecies).describe())
    # print(pd.DataFrame(non_complexSpecies).describe())

    print("-------------------------------------")

    # labels = ["Complex", "Non-complex"]

    # figure,axes=plt.subplots(figsize=(6,4))

    # # rectangular box plot 
    # boxplot = axes.boxplot([np.array(complexSpecies).astype(np.float),np.array(non_complexSpecies).astype(np.float)],
    #                           vert=True, # vertical box alignment
    #                           patch_artist=True, # fill with color
    #                           labels=labels) # will be used to label x-ticks
    # axes.set_title(organism.replace("_", " "))

    # # fill with colors
    # colors = ['lightblue', 'lightgreen']

    # for patch, color in zip(boxplot['boxes'], colors):
    #     patch.set_facecolor(color)

    # # adding horizontal grid lines
    # axes.yaxis.grid(True)

    # axes.set_ylabel("Number of yeast species")

    # plt.savefig("../figure/boxplot/%s_boxplot.png" % organism, dpi=400)


# Results :
