#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-02-14

# This python script is to plot boxplopt for gene essentiality data and number of yeast species according to their othologs
# http://cmdlinetips.com/2019/03/how-to-make-grouped-boxplots-in-python-with-seaborn/
# https://github.com/cdanielmachado/cooccurrence/blob/master/notebooks/Figure%205.ipynb


import json
import csv
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt 


organisms = ["S_cerevisiae", "S_pombe",  "C_albicans", "Y_lipolytica", "P_pastoris"]

# # with open("../Data/essential.csv", "w") as outfile :
outfile = open("../complementaryData/boxplot_data/occurence.csv", "w")
csv_writer = csv.writer(outfile)
csv_writer.writerow(["type", "organism", "species"])

for organism in organisms :
    # allEssential = set()
    print("This organism is: %s" % organism.replace("_", " "))
    with open("../complementaryData/processed_data/%s_include_ortholog.json" % organism, "r") as f :
        data = json.load(f)

    for essentialdata in data :
        csv_writer.writerow([list(essentialdata.values())[0], organism.split('_')[0][0]+'. '+organism.split('_')[1] , essentialdata['species']])

outfile.close()

alldata = pd.read_csv("../complementaryData/boxplot_data/occurence.csv")
print(alldata.head(3))


plt.figure(figsize=(3.,2.4))

# rectangular box plot
palette = {"E": '#ed7e17', "NE": '#1ba055'}

for ind in alldata.index:
    alldata.loc[ind,'organism'] = '${0}$'.format(alldata.loc[ind,'organism'])

ax = sns.boxplot(data=alldata, x="organism", y="species", hue="type",
        palette=palette, showfliers=False, linewidth=1)

# ax = sns.stripplot(data=alldata, x="organism", y="species", hue="type", palette=palette, 
#           dodge=True, size=2, linewidth=0.5, alpha=0.3)

# https://stackoverflow.com/questions/58476654/how-to-remove-or-hide-x-axis-label-from-seaborn-boxplot
# plt.xlabel(None) will remove the Label, but not the ticks. 
ax.set(xlabel=None)
# plt.xlabel("Organism")

# for tick in ax.get_xticklabels() :
#     tick.set_rotation(30)

plt.ylabel("Gene occurence number")
# plt.ylabel("Paralog_num")

# plt.ylim(0,450)
plt.ylim(0,600)

plt.xticks(rotation=30,ha='right')
# plt.yticks([0,150,300,450])
plt.yticks([0,100,200,300,400,500,600])
# # ax.legend(ax.get_legend_handles_labels()[0], ["E", "NE"])

handles,labels = ax.get_legend_handles_labels()
# # specify just one legend
l = plt.legend(handles[0:2], labels[0:2], loc=2)

# https://blog.csdn.net/weixin_38314865/article/details/88633880
plt.savefig("../complementaryData/figure/occurence_boxplot_italic.pdf", dpi=400, bbox_inches = 'tight')


# Results :
