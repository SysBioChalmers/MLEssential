#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-08-07

import csv
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#plt.barh
# feature = ['Kmer-AGA', 'DNC-CT', 'Gene occurence', 'Kmer-AAG', 'dN/dS']

# ChI_value = [748.310, 758.131, 1013.099, 1150.056, 1272.875]

# # plt.figure(figsize=(3.,2.4))
# plt.figure(figsize=(2.7,2.1))

# # https://juejin.im/post/6858230839986421767
# plt.barh(range(len(ChI_value)), ChI_value, tick_label=feature, alpha=0.8, height=0.4, color='pink', edgecolor='r')


# ax = plt.axes()
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)

# plt.xlim(0,1500)

# plt.xticks([0,500,1000,1500])

# plt.xlabel('CHI value')

# plt.savefig("../complementaryData/figure/feature_importance.pdf", dpi=400, bbox_inches='tight')


# plt.bar
# feature = ['dN/dS', 'Kmer-AAG', 'Gene occurence', 'DNC-CT', 'Kmer-AGA']
# feature = ['dN/dS', 'Kmer-AAG', 'Gene occurence', 'DNC-CT', 'Kmer-AGA']
feature = ['Gene occurence', 'dN/dS', 'Conservation score', 'Kmer-TGA', 'Kmer-GAA']

# ChI_value = [1272.875, 1150.056, 1013.099, 758.131, 748.310]
ChI_value = [1765.265, 1649.817, 599.300, 350.908, 340.391]

# feature = ['Kmer-AGA', 'DNC-CT', 'Gene occurence', 'Kmer-AAG', 'dN/dS']

# ChI_value = [748.310, 758.131, 1013.099, 1150.056, 1272.875]

# plt.figure(figsize=(3.,2.4))
plt.figure(figsize=(2.7,2.1))

# https://juejin.im/post/6858230839986421767
plt.bar(range(len(ChI_value)), ChI_value, tick_label=feature, width=0.5, alpha=0.8, color='pink', edgecolor='r')


ax = plt.axes()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.ylim(0,2000)

plt.xticks(rotation=30, ha='right')
plt.yticks([0,500,1000,1500,2000])

# plt.ylabel('CHI value')
plt.ylabel('Feature importance')

plt.savefig("../complementaryData/figure/feature_importance_bar.pdf", dpi=400, bbox_inches='tight')
