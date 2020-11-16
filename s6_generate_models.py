#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import json
import pickle
import numpy as np
from sklearn import svm
from pubscripts import *
from descnucleotide import *

def try_svm() :
    X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    y = np.array([1, 1, 2, 2])
    # y = np.array([[1], [1], [2], [2]])
    model = svm.SVC()
    model.fit(X, y)

    prediction = model.predict([[-0.8, -1], [1, 1], [1, 1]]).tolist()
    print(type(prediction))
    print(prediction)

def random(data) :
    index = [i for i in range(len(data))]
    np.random.shuffle(index)
    newdata = [data[i] for i in index]
    return newdata

def generate_model_data() :
    organisms = ["Y_lipolytica", "S_pombe", "S_cerevisiae", "P_pastoris", "C_albicans"]
    all_data = list()

    for organism in organisms :
        print("This organism is: %s" % organism)
        with open("./complementaryData/processed_data/%s_include_ortholog.json" % organism, "r") as f :
            data = json.load(f)

        for subdata in data :
            essential_dict = dict()
            non_essential_dict = dict()

            if list(subdata.values())[0] == "E" :
                essential_dict[list(subdata.keys())[0]] = "E"
                essential_dict["gene sequence"] = subdata["gene sequence"]
                all_data.extend([essential_dict]*4)
                # all_data.extend([essential_dict])

            if list(subdata.values())[0] == "NE" :
                non_essential_dict[list(subdata.keys())[0]] = "NE"
                non_essential_dict["gene sequence"] = subdata["gene sequence"]
                all_data.append(non_essential_dict)

    print("The number of all_data: %d" % len(all_data))  # 22981
    print("-----------------------------------------")

    model_data = random(all_data)

    with open('model_data.txt', 'w') as f :
        essential = {'E':1, 'NE':0}
        for train in model_data :
            keys = list(train.keys())
            values = list(train.values())
            f.write('>%s|%s|model' % (keys[0],essential[values[0]]))
            f.write('\n')
            f.write(train["gene sequence"])
            f.write('\n')

def pickle_svm_model() :
    fastas = []
    cmd_coding = {}
    model_data = []
    model_code_dict = {}
    features = []
    labels = []

    parameters = {'Method': "DNDS;Conservation;Occurance;ProteinNumber;DNC;Kmer", 'Kmer_Size': 3}
    dna_cmd_coding = {
        'Kmer': 'Kmer.Kmer(model_data, k=%s, **kw)' % parameters['Kmer_Size'],
        'DNC': 'DNC.DNC(model_data, **kw)',
        'DNDS': 'DNDS.dnds(model_data, **kw)',
        'Conservation': 'Conservation.conservation_score(model_data, **kw)',
        'Occurance': 'Occurance.occurance_number(model_data, **kw)',
        'ProteinNumber': 'ProteinNumber.protein_number(model_data, **kw)',
    }
    
    fastas = read_fasta_sequences.read_nucleotide_sequences('model_data.txt')
    cmd_coding = dna_cmd_coding

    for sequence in fastas:
        if sequence[3] == 'model':
            model_data.append(sequence)

    kw = {'nclusters': 3, 'sof': 'sample', 'order': ''}
    method_array = parameters['Method'].split(';')
    for method in method_array :
        if method in ('DNC', 'Kmer'):
            kw['order'] = 'ACGT'
        model_code_dict[method] = eval(cmd_coding[method])

    model_code = np.array(model_code_dict[method_array[0]])

    for i in range(1, len(method_array)):
        # print(model_code)
        # print(type(model_code))
        if model_code_dict[method_array[i]] != 0:
            model_code = np.concatenate((model_code, np.array(model_code_dict[method_array[i]])[:, 2:]), axis=1)

    model_code = model_code.tolist()
    # print(model_code[0])
    # print(model_code[1])
    # print(len(model_code))  # 22982  The first list is one explanation for the following lists

    for info in model_code[1:] :
        features.append(info[2:])
        labels.append(info[1])

    # print(features[:10])
    # print(labels[:10])
    # print(len(features))
    # print(len(labels))

    features = np.array(features)
    labels = np.array(labels)

    svm_model = svm.SVC(C=15, kernel="rbf", degree=3, gamma=8, coef0=0, probability=True, random_state=1)
    svm_model.fit(features,labels)

    file = open('./model.pickle', 'wb')
    pickle.dump(svm_model, file)
    file.close()


if __name__ == "__main__" :
    # try_svm()
    # generate_model_data()
    pickle_svm_model()

