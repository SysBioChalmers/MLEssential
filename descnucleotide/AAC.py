#!/usr/bin/env python
#_*_coding:utf-8_*_

import re
import json
from collections import Counter


def getProSeq(name) :
    with open("./complementaryData/processed_data/allProtein.json", "r") as file :
        data = json.load(file)
    protein_sequence = data[name]

    return protein_sequence

def AAC(fastas, **kw):
    AA = kw['order'] if kw['order'] != None else 'ACDEFGHIKLMNPQRSTVWY'
    #AA = 'ARNDCQEGHILKMFPSTWYV'
    encodings = []
    header = ['#', 'label']
    for i in AA:
        header.append(i)
    encodings.append(header)

    for i in fastas:
        # name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        name, label = i[0], i[2]
        sequence = getProSeq(name)
        count = Counter(sequence)
        for key in count:
            count[key] = count[key]/len(sequence)
        code = [name, label]
        for aa in AA:
            code.append(count[aa])
        encodings.append(code)
    return encodings

# print(encodings)
# [['#' 'label' 'A' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'K' 'L' 'M' 'N' 'P' 'Q' 'R'
#   'S' 'T' 'V' 'W' 'Y']
#  ['orf19.4930' '1' '0.057291666666666664' '0' '0.036458333333333336'
#   '0.046875' '0.0625' '0.046875' '0.005208333333333333'
#   '0.06770833333333333' '0.08333333333333333' '0.08333333333333333'
#   '0.005208333333333333' '0.08333333333333333' '0.03125' '0.046875'
#   '0.036458333333333336' '0.08854166666666667' '0.10416666666666667'
#   '0.057291666666666664' '0.026041666666666668' '0.03125']
#  ['RHTO_05824' '0' '0.16192733017377567' '0.004739336492890996'
#   '0.07266982622432859' '0.10821484992101106' '0.011058451816745656'
#   '0.061611374407582936' '0.011848341232227487' '0.018957345971563982'
#   '0.07661927330173776' '0.0947867298578199' '0.008688783570300158'
#   '0.018167456556082148' '0.03712480252764613' '0.05055292259083728'
#   '0.06714060031595577' '0.07345971563981042' '0.04265402843601896'
#   '0.06398104265402843' '0.001579778830963665' '0.014218009478672985']]

