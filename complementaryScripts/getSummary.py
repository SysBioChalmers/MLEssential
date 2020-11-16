#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-12-05

# This python script is to obtain coding sequence and protein sequence for essential genes and non-essential genes

import os
from collections import Counter
import pandas as pd


# read the yeast file
def read_xls(filename) :
    data = pd.read_excel("../Data/essential/%s" %(filename))
    df = pd.DataFrame(data)
    genes = df.set_index("Gene id")["Essential genes (E) / Non-essential genes (NE)"].to_dict()
    # print(len(genes))

    # The number of essential genes and non-essential genes, return results like Counter({'NE': 1714, 'E': 633})
    numbers = Counter(genes.values())

    # output the E and NE gene number for a specific file
    print("The number of essential genes for %s is: %d" %(filename[:-5],numbers["E"]))
    print("The number of non-essential genes for {} is: {}" .format(filename[:-5],numbers["NE"]))


def main() :
    yeasts = ["C. albicans", "P. pastoris", "R. toruloides", "S. cerevisiae", "S. pombe", "Y. lipolytica"]

    for yeast in yeasts :
        read_xls("%s.xlsx" %(yeast))


if __name__ == "__main__" :
    main()



