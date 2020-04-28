#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-03-20

# This python script is to obtain training data and testing data from preprocessed file for machine learning prediction

import os
import json
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split


def conservation():
    with open("./complementaryData/evolution_feature/conservation_score.txt", "r") as outfile :
        conservation_data = outfile.readlines()

        conservation = dict()
        for line in conservation_data :
            ortholog = line.strip().split(" ")[2][1:-1]
            conservation_score = line.strip().split(" ")[3]
            conservation[ortholog] = float(conservation_score)
        conservation_ortholog = set(conservation.keys())

    return conservation_ortholog

def dnds():
    with open("./complementaryData/evolution_feature/gene_dn_ds.csv", "r") as outfile :
        dnds_data = outfile.readlines()
        # print(dnds_data)
        dnds = dict()

        for line in dnds_data :
            ortholog = line.strip().split(",")[1].split(".")[0]
            dnds_score = line.strip().split(",")[2]
            # print(type(dnds_score))  # <class 'str'>

            if dnds_score :
                dnds[ortholog] = float(dnds_score)

        dnds_ortholog = set(dnds.keys())

    return dnds_ortholog

conservation_ortholog = conservation()
dnds_ortholog = dnds()

interaction_orthlog = conservation_ortholog & dnds_ortholog
# print(len(interaction_orthlog))  6875

def random(data) :
    index = [i for i in range(len(data))]
    np.random.shuffle(index)
    newdata = [data[i] for i in index]
    return newdata

organisms = ["Candida_glabrata", "Candida_dubliniensis", "Candida_parapsilosis", "Candida_tropicalis", "Candida_albicans", "Yarrowia_lipolytica", "Schizosaccharomyces_pombe", "Saccharomyces_cerevisiae"]

training_data = list()
testing_data = list()

for organism in organisms :
    print("This organism is: %s" % organism.replace("_", " "))
    with open("./complementaryData/newjson/%s.json" % organism, "r") as f :
        data = json.load(f)

    if organism in ["Candida_dubliniensis","Schizosaccharomyces_pombe"] :
        for subdata in data :
            essential_dict = dict()
            non_essential_dict = dict()

            if list(subdata.values())[0] == "Complex" and subdata['ortholog'] in interaction_orthlog :
                essential_dict[list(subdata.keys())[0]] = "Complex"
                essential_dict["gene sequence"] = subdata["gene sequence"]
                testing_data.extend([essential_dict]*3)

            if list(subdata.values())[0] == "Non-complex" and subdata['ortholog'] in interaction_orthlog :
                non_essential_dict[list(subdata.keys())[0]] = "Non-complex"
                non_essential_dict["gene sequence"] = subdata["gene sequence"]
                testing_data.append(non_essential_dict)

    else :
        for subdata in data :
            essential_dict = dict()
            non_essential_dict = dict()

            if list(subdata.values())[0] == "Complex" and subdata['ortholog'] in interaction_orthlog :
                essential_dict[list(subdata.keys())[0]] = "Complex"
                essential_dict["gene sequence"] = subdata["gene sequence"]
                training_data.extend([essential_dict]*3)

            if list(subdata.values())[0] == "Non-complex" and subdata['ortholog'] in interaction_orthlog :
                non_essential_dict[list(subdata.keys())[0]] = "Non-complex"
                non_essential_dict["gene sequence"] = subdata["gene sequence"]
                training_data.append(non_essential_dict)



print("The number of training_data genes: %d" % len(training_data))
print("The number of testing_data genes: %d" % len(testing_data))
print("-----------------------------------------")

# This organism is: Y_lipolytica
# This organism is: C_albicans
# This organism is: S_pombe
# This organism is: S_cerevisiae
# This organism is: R_toruloides
# This organism is: P_pastoris
# The number of essential_data genes: 4369
# The number of non-essential_data genes: 16051
# -----------------------------------------
# [Finished in 1.3s]

# all_data = essential_data + non_essential_data

# # print(all_data[:3])
# print(len(all_data))

X_train = random(training_data)
X_test = random(testing_data)

# print(len(X_train))
# print(len(X_test))

with open('traintest.txt', 'w') as f :
    essential = {'Complex':1, 'Non-complex':0}
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

