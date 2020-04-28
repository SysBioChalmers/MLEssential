#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-02-25

# This python script is to obtain the unique orthologs for all eight yeast species


import json
import csv
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


organisms = ["Candida_albicans", "Candida_glabrata", "Candida_dubliniensis", "Candida_parapsilosis", "Candida_tropicalis", "Yarrowia_lipolytica", "Schizosaccharomyces_pombe", "Saccharomyces_cerevisiae"]

outfile = open("../evolution_feature/ortholog_complex.csv", "w")
csv_writer = csv.writer(outfile)
# csv_writer.writerow(["type", "organism", "species"])

allOrtholog_list = list()
for organism in organisms :
    print("This organism is: %s" % organism.replace("_", " "))
    with open("../processed_data/%s_include_ortholog.json" % organism, "r") as f :
        data = json.load(f)

    for complex in data :
        # allOrtholog.add((complex["ortholog"]))
        allOrtholog_list.append(complex["ortholog"])

allOrtholog = list(set(allOrtholog_list))

    # print("Ortholog number for %s is: %d" %(organism, len(allOrtholog)))

print("The amount of orthologs for all eight yeast species: %d" % len(allOrtholog))
for ortholog in allOrtholog :
    csv_writer.writerow([ortholog])

outfile.close()


# Results :

# The amount of orthologs for all eight yeast species: 11727

