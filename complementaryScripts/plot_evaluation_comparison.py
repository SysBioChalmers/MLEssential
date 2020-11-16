#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-02-07

import json
import csv
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

ticks = ['AUC','Accuracy','Precision','F1 score']

seq = [0.68,0.63,0.63,0.68]
seq_evo = [0.91,0.87,0.84,0.88]
# plt.figure(figsize=(3.,2.4))
# plt.figure(figsize=(2.7,2.1))
plt.figure(figsize=(2.7,2.1))
pos = np.arange(len(ticks))*4

plt.bar(pos,seq,label='Seq',width=1.0,color='#1b9e77')
plt.bar(pos+1,seq_evo,label='Seq&Evo',width=1.0,color='#7570b3')

ax = plt.axes()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.ylabel('Percentage')
plt.xticks(pos+0.5,ticks,rotation=30,ha='right')
# plt.ylim([0.5,1.0])
# plt.yticks([0.6,0.7,0.8,0.9,1.0])
plt.ylim([0,1.0])
plt.yticks([0.2,0.4,0.6,0.8,1.0])
# plt.legend(loc='upper right',fontsize=8)
# plt.legend(fontsize=8)
plt.tight_layout()

plt.savefig("../complementaryData/figure/plot_comparison.pdf", dpi=400, bbox_inches='tight')

