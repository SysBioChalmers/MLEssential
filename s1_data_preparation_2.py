#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-12-12

# This python script is to obtain training data and testing data from preprocessed file for machine learning prediction

import os
import json
import pandas as pd
from sklearn.model_selection import train_test_split


organisms = {"Y_lipolytica", "S_pombe", "S_cerevisiae", "P_pastoris", "C_albicans"}

all_data = list()
essential_data = list()
non_essential_data = list()

for organism in organisms :
    # print("This organism is: %s" % organism)
    with open("./complementaryData/processed_data/%s_include_ortholog.json" % organism, "r") as f :
        data = json.load(f)

    for subdata in data :
        essential_dict = dict()
        non_essential_dict = dict()

        if list(subdata.values())[0] == "E" :
            essential_dict[list(subdata.keys())[0]] = "E"
            essential_dict["gene sequence"] = subdata["gene sequence"]
            essential_data.append(essential_dict)

        if list(subdata.values())[0] == "NE" :
            non_essential_dict[list(subdata.keys())[0]] = "NE"
            non_essential_dict["gene sequence"] = subdata["gene sequence"]
            non_essential_data.append(non_essential_dict)


# print("The number of essential_data genes: %d" % len(essential_data))
# print("The number of non-essential_data genes: %d" % len(non_essential_data))
# print("-----------------------------------------")

all_data = essential_data * 4 + non_essential_data

# print(all_data[:3])
print(len(all_data))

X_train, X_test = train_test_split(all_data, test_size=0.2)

# print(len(X_train))
# print(len(X_test))

with open('traintest.txt', 'w') as f :
    essential = {'E':1, 'NE':0}
    for train in X_train :
        keys = list(train.keys())
        values = list(train.values())
        f.write('>%s|%s|training' % (keys[0],essential[values[0]]))
        f.write('\n')
        f.write(train["gene sequence"])
        f.write('\n')

    for test in X_test :
        keys = list(test.keys())
        values = list(test.values())
        f.write('>%s|%s|testing' % (keys[0],essential[values[0]]))
        f.write('\n')
        f.write(test["gene sequence"])
        f.write('\n')

