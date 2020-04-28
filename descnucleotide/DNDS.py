#!/usr/bin/env python
#_*_coding:utf-8_*_

import re
import json


def getOrtholog() :
    organisms = ["Candida_glabrata", "Candida_dubliniensis", "Candida_parapsilosis", "Candida_tropicalis", "Candida_albicans", "Yarrowia_lipolytica", "Schizosaccharomyces_pombe", "Saccharomyces_cerevisiae"]
    
    allOrtholog = dict()
    for organism in organisms :
        # print("This organism is: %s" % organism.replace("_", "."))
        with open("./complementaryData/newjson/%s.json" % organism, "r") as f :
            data = json.load(f)

        for essential in data :
            # allOrtholog.add((essential["ortholog"]))
            allOrtholog[list(essential.keys())[0]] = essential["ortholog"]
    return allOrtholog

def getDNDS() :
    with open("./complementaryData/evolution_feature/gene_dn_ds.csv", "r") as outfile :
        dnds_data = outfile.readlines()

        dnds_dict = dict()
        for line in dnds_data :
            ortholog = line.strip().split(",")[1].split(".")[0]
            dnds_score = line.strip().split(",")[2]
            if dnds_score :
                dnds_dict[ortholog] = float(dnds_score)
    return dnds_dict

def dnds(fastas, **kw):
    allOrtholog = getOrtholog()
    dnds_dict = getDNDS()

    encodings = []
    feature = ['dnds']
    # ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    header = ['#', 'label'] + feature
    # encodings: [['#', 'label', 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']]
    encodings.append(header)

    for i in fastas :
        name, label = i[0], i[2]

        code = [name, label]
        try :
            ortholog = allOrtholog[name]
        except :
            pass

        # try :
        #     tmpCode = [dnds_dict[ortholog]]
        # except :
        #     if label == '0' :
        #         tmpCode = [0.12]
        #     if label == '1' :
        #         tmpCode = [0.08]

        tmpCode = [dnds_dict[ortholog]]

        code = code + tmpCode
        encodings.append(code)
    # print(encodings)
    # [['#', 'label', 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'], 
    # ['YBR076W', '0', 0.13815789473684212, 0.06390977443609022, 0.05357142857142857, 0.09210526315789473, 
    # 0.07048872180451128, 0.044172932330827065, 0.02537593984962406, 0.05169172932330827, 0.06954887218045112, 
    # 0.02725563909774436, 0.03853383458646616, 0.042293233082706765, 0.06954887218045112, 0.05639097744360902, 
    # 0.06015037593984962, 0.09680451127819549],,,,]
    return encodings

