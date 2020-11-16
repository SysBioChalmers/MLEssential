#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-02-14

# http://cmdlinetips.com/2019/03/how-to-make-grouped-boxplots-in-python-with-seaborn/
# https://github.com/cdanielmachado/cooccurrence/blob/master/notebooks/Figure%205.ipynb


import json
import csv
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt 


with open("../complementaryData/evolutionary_data/ortholog_occurence.tsv", "r") as outfile :
    # paralog_data = outfile.readlines()
    paralog_data = csv.reader(outfile,delimiter='\t')
    # print(paralog_data)
    paralog = dict()

    for line in list(paralog_data) :
        ortholog = line[1]
        paralog_number = line[-1]
        # print(type(paralog_number))  # <class 'str'>
        paralog[ortholog] = float(paralog_number)


organisms = ["S_cerevisiae", "S_pombe",  "C_albicans", "Y_lipolytica", "P_pastoris"]


# # with open("../Data/essential.csv", "w") as outfile :
outfile = open("../complementaryData/boxplot_data/paralog.csv", "w")
csv_writer = csv.writer(outfile)
csv_writer.writerow(["type", "organism", "paralog"])

for organism in organisms :
    # allEssential = set()
    print("This organism is: %s" % organism.replace("_", " "))
    with open("../complementaryData/processed_data/%s_include_ortholog.json" % organism, "r") as f :
        data = json.load(f)

    # for essentialdata in data :
    #     allessentialdata.add((essentialdata["id"]))
    # # csv_writer.writerow(list(data))
    for essentialdata in data :
        ortholog = essentialdata['ortholog']
        csv_writer.writerow([list(essentialdata.values())[0], organism.split('_')[0][0]+'. '+organism.split('_')[1] , paralog[ortholog]])

outfile.close()

alldata = pd.read_csv("../complementaryData/boxplot_data/paralog.csv")
print(alldata.head(3))

plt.figure(figsize=(3.,2.4))

# rectangular box plot
palette = {"E": '#ed7e17', "NE": '#1ba055'}

for ind in alldata.index:
    alldata.loc[ind,'organism'] = '${0}$'.format(alldata.loc[ind,'organism'])

ax = sns.boxplot(data=alldata, x="organism", y="paralog", hue="type",
        palette=palette, showfliers=False, linewidth=1)

# ax = sns.stripplot(data=alldata, x="organism", y="paralog", hue="type", palette=palette, 
#           dodge=True, size=2, linewidth=0.5, alpha=0.3)

# https://stackoverflow.com/questions/58476654/how-to-remove-or-hide-x-axis-label-from-seaborn-boxplot
# plt.xlabel(None) will remove the Label, but not the ticks. 
ax.set(xlabel=None)
# plt.xlabel("Organism")

# for tick in ax.get_xticklabels() :
#     tick.set_rotation(30)

plt.ylabel("Number of average paralogs")

plt.ylim(1.0,1.4)

plt.xticks(rotation=30,ha='right')
plt.yticks([1.0,1.1,1.2,1.3,1.4])
# # ax.legend(ax.get_legend_handles_labels()[0], ["E", "NE"])

handles,labels = ax.get_legend_handles_labels()
# # specify just one legend
l = plt.legend(handles[0:2], labels[0:2], loc=0)

# https://blog.csdn.net/weixin_38314865/article/details/88633880
plt.savefig("../complementaryData/figure/paralog_boxplot_italic.pdf", dpi=400, bbox_inches = 'tight')

