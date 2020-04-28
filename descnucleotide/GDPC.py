#!/usr/bin/env python
# _*_coding:utf-8_*_

import re
import json


def getProSeq(name) :
    with open("./complementaryData/processed_data/allProtein.json", "r") as file :
        data = json.load(file)
    protein_sequence = data[name]

    return protein_sequence

def GDPC(fastas, **kw):
    group = {
        'alphaticr': 'GAVLMI',
        'aromatic': 'FYW',
        'postivecharger': 'KRH',
        'negativecharger': 'DE',
        'uncharger': 'STCPNQ'
    }

    groupKey = group.keys()
    baseNum = len(groupKey)
    dipeptide = [g1 + '.' + g2 for g1 in groupKey for g2 in groupKey]

    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key

    encodings = []
    header = ['#', 'label'] + dipeptide
    encodings.append(header)

    for i in fastas[:3]:
        name, label = i[0], i[2]
        sequence = getProSeq(name)

        code = [name, label]
        myDict = {}
        for t in dipeptide:
            myDict[t] = 0

        sum = 0
        for j in range(len(sequence) - 2 + 1):
            myDict[index[sequence[j]] + '.' + index[sequence[j + 1]]] = myDict[index[sequence[j]] + '.' + index[
                sequence[j + 1]]] + 1
            sum = sum + 1

        if sum == 0:
            for t in dipeptide:
                code.append(0)
        else:
            for t in dipeptide:
                code.append(myDict[t] / sum)
        encodings.append(code)
        print(encodings)
    return encodings
