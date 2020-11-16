#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-02-07

# This python script is to plot boxplopt for gene essentiality data and number of yeast species according to their othologs
# http://cmdlinetips.com/2019/03/how-to-make-grouped-boxplots-in-python-with-seaborn/
# https://github.com/cdanielmachado/cooccurrence/blob/master/notebooks/Figure%205.ipynb


import json
import csv
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


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

# organisms = ["Y_lipolytica", "S_pombe", "S_cerevisiae", "P_pastoris", "C_albicans"]
organisms = ["S_cerevisiae", "S_pombe",  "C_albicans", "Y_lipolytica", "P_pastoris",]

# with open("../Data/essential.csv", "w") as outfile :
outfile = open("../complementaryData/boxplot_data/dnds.csv", "w")
csv_writer = csv.writer(outfile)
csv_writer.writerow(["type", "organism", "dnds"])

i = 0
for organism in organisms :
    # allEssential = set()
    print("This organism is: %s" % organism.replace("_", ". "))
    with open("../complementaryData/processed_data/%s_include_ortholog.json" % organism, "r") as f :
        data = json.load(f)

    for essentialdata in data :
        ortholog = essentialdata['ortholog']
        try : 
            csv_writer.writerow([list(essentialdata.values())[0], organism.split('_')[0][0]+'. '+organism.split('_')[1], dnds[ortholog]])
        except :
            # print(ortholog)
            i +=1
print(i)

outfile.close()

alldata = pd.read_csv("../complementaryData/boxplot_data/dnds.csv")

# plt.figure(figsize=(3.,2.4))
plt.figure(figsize=(2.7,2.1))

labels = ["E", "NE"]

    # figure,axes=plt.subplots(figsize=(6,3))

# rectangular box plot 
palette = {"E": '#ed7e17', "NE": '#1ba055'}

for ind in alldata.index:
    alldata.loc[ind,'organism'] = '${0}$'.format(alldata.loc[ind,'organism'])

ax = sns.boxplot(data=alldata, x="organism", y="dnds", hue="type", 
        palette=palette, showfliers=False, linewidth=1)

# ax = sns.stripplot(data=alldata, x="organism", y="dnds", hue="type", palette=palette, 
#           dodge=True, size=3, linewidth=0.5, alpha=0.3)

# https://stackoverflow.com/questions/58476654/how-to-remove-or-hide-x-axis-label-from-seaborn-boxplot
# plt.xlabel(None) will remove the Label, but not the ticks. 
ax.set(xlabel=None)
# plt.xlabel("Organism")

# for tick in ax.get_xticklabels() :
#     tick.set_rotation(30)

plt.ylabel("dN/dS")

plt.ylim(0,0.5)

plt.xticks(rotation=30,ha='right')
plt.yticks([0,0.1,0.2,0.3,0.4,0.5])
# # ax.legend(ax.get_legend_handles_labels()[0], ["E", "NE"])

handles,labels = ax.get_legend_handles_labels()
# # specify just one legend
l = plt.legend(handles[0:2], labels[0:2], loc=0, fontsize=8)

sns.despine(top=True,right=True)

plt.savefig("../complementaryData/figure/dnds_boxplot_italic.pdf", dpi=400, bbox_inches='tight')

