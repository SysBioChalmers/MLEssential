#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-12-14

# This python script is to plot boxplopt for gene essentiality data and number of yeast species according to their orthologs


import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


organisms = {"Y_lipolytica", "S_pombe", "S_cerevisiae", "P_pastoris", "C_albicans"}

for organism in organisms :
    print("This organism is: %s" % organism)
    with open("../Data/processed_data/%s_include_ortholog.json" % organism, "r") as f :
        data = json.load(f)

    essential = [float(essential["species"]) for essential in data if list(essential.values())[0] == "E"]
    non_essential = [float(essential["species"]) for essential in data if list(essential.values())[0] == "NE"]

    print("The number of essential genes: %d" % len(essential))
    print("The number of non-essential genes: %d" % len(non_essential))

    # https://blog.csdn.net/aijiudu/article/details/89387328
    print(pd.DataFrame(essential).describe())
    print(pd.DataFrame(non_essential).describe())
    print("-------------------------------------")


    labels = ["Essential", "Non-essential"]

    figure,axes=plt.subplots(figsize=(6,4))

    # rectangular box plot 
    boxplot = axes.boxplot([np.array(essential).astype(np.float),np.array(non_essential).astype(np.float)],
                              vert=True, # vertical box alignment
                              patch_artist=True, # fill with color
                              labels=labels) # will be used to label x-ticks
    axes.set_title(organism)

    # fill with colors
    colors = ['lightblue', 'lightgreen']

    for patch, color in zip(boxplot['boxes'], colors):
        patch.set_facecolor(color)

    # adding horizontal grid lines
    axes.yaxis.grid(True)

    axes.set_ylabel("Number of yeast species")

    plt.savefig("../Data/figure/%s_boxplot.png" % organism, dpi=400)


# Results :

# This organism is: S_cerevisiae
# The number of essential genes: 1035
# The number of non-essential genes: 4503
#                  0
# count  1035.000000
# mean    294.990338
# std      83.666707
# min       2.000000
# 25%     301.000000
# 50%     329.000000
# 75%     340.000000
# max     343.000000
#                  0
# count  4503.000000
# mean    222.506773
# std     127.852372
# min       1.000000
# 25%      71.000000
# 50%     290.000000
# 75%     329.000000
# max     343.000000
# -------------------------------------
# This organism is: P_pastoris
# The number of essential genes: 138
# The number of non-essential genes: 450
#                 0
# count  138.000000
# mean   332.746377
# std     15.763332
# min    264.000000
# 25%    329.000000
# 50%    340.000000
# 75%    342.000000
# max    343.000000
#                 0
# count  450.000000
# mean   320.055556
# std     32.200154
# min    105.000000
# 25%    313.000000
# 50%    331.000000
# 75%    341.000000
# max    343.000000
# -------------------------------------
# This organism is: Y_lipolytica
# The number of essential genes: 107
# The number of non-essential genes: 520
#                 0
# count  107.000000
# mean   324.728972
# std     43.835794
# min      9.000000
# 25%    327.500000
# 50%    339.000000
# 75%    342.000000
# max    343.000000
#                 0
# count  520.000000
# mean   293.844231
# std     74.750148
# min      5.000000
# 25%    284.000000
# 50%    325.000000
# 75%    340.000000
# max    343.000000
# -------------------------------------
# This organism is: C_albicans
# The number of essential genes: 625
# The number of non-essential genes: 1682
#                 0
# count  625.000000
# mean   320.056000
# std     40.112842
# min      4.000000
# 25%    316.000000
# 50%    335.000000
# 75%    342.000000
# max    343.000000
#                  0
# count  1682.000000
# mean    271.024376
# std      86.358259
# min       1.000000
# 25%     244.000000
# 50%     307.000000
# 75%     332.000000
# max     343.000000
# -------------------------------------
# This organism is: S_pombe
# The number of essential genes: 1224
# The number of non-essential genes: 3342
#                  0
# count  1224.000000
# mean    279.892157
# std     106.163174
# min       1.000000
# 25%     291.000000
# 50%     326.000000
# 75%     340.000000
# max     343.000000
#                  0
# count  3342.000000
# mean    197.079294
# std     144.907726
# min       1.000000
# 25%      10.000000
# 50%     277.000000
# 75%     329.000000
# max     343.000000
# -------------------------------------
# [Finished in 4.9s]
#                      