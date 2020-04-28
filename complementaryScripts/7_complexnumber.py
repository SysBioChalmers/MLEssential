#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-12-14

# This python script is to plot boxplopt for gene essentiality data and number of yeast species according to their othologs


import json
import numpy as np
import pandas as pd
from scipy.stats import ranksums


organisms = ["Candida_albicans", "Candida_glabrata", "Candida_dubliniensis", "Candida_parapsilosis", "Candida_tropicalis", "Yarrowia_lipolytica", "Schizosaccharomyces_pombe", "Saccharomyces_cerevisiae"]


for organism in organisms :
    print("This organism is: %s" % organism.replace("_", " "))

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

    print("-------------------------------------")

# Results:

# This organism is: Candida albicans
# The number of complex proteins for Candida_albicans: 1880
# The number of non-complex proteins: 4304
# -------------------------------------
# This organism is: Candida glabrata
# The number of complex proteins for Candida_glabrata: 1671
# The number of non-complex proteins: 3227
# -------------------------------------
# This organism is: Candida dubliniensis
# The number of complex proteins for Candida_dubliniensis: 1629
# The number of non-complex proteins: 3745
# -------------------------------------
# This organism is: Candida parapsilosis
# The number of complex proteins for Candida_parapsilosis: 1571
# The number of non-complex proteins: 3423
# -------------------------------------
# This organism is: Candida tropicalis
# The number of complex proteins for Candida_tropicalis: 466
# The number of non-complex proteins: 4971
# -------------------------------------
# This organism is: Yarrowia lipolytica
# The number of complex proteins for Yarrowia_lipolytica: 1236
# The number of non-complex proteins: 4773
# -------------------------------------
# This organism is: Schizosaccharomyces pombe
# The number of complex proteins for Schizosaccharomyces_pombe: 1585
# The number of non-complex proteins: 3431
# -------------------------------------
# This organism is: Saccharomyces cerevisiae
# The number of complex proteins for Saccharomyces_cerevisiae: 2011
# The number of non-complex proteins: 3719
# -------------------------------------
# [Finished in 1.1s]
