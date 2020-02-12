#!/usr/bin/env python
#_*_coding:utf-8_*_

import re

def DNC(fastas, **kw):
    base = 'ACGT'

    encodings = []
    dinucleotides = [n1 + n2 for n1 in base for n2 in base]
    # ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    header = ['#', 'label'] + dinucleotides
    # encodings: [['#', 'label', 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']]
    encodings.append(header)

    AADict = {}
    for i in range(len(base)):
        AADict[base[i]] = i

    for i in fastas :
        name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        code = [name, label]
        tmpCode = [0] * 16
        for j in range(len(sequence) - 2 + 1):
            tmpCode[AADict[sequence[j]] * 4 + AADict[sequence[j+1]]] = tmpCode[AADict[sequence[j]] * 4 + AADict[sequence[j+1]]] +1
        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]
        code = code + tmpCode
        encodings.append(code)
    # print(encodings)
    # [['#', 'label', 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'], 
    # ['YBR076W', '0', 0.13815789473684212, 0.06390977443609022, 0.05357142857142857, 0.09210526315789473, 
    # 0.07048872180451128, 0.044172932330827065, 0.02537593984962406, 0.05169172932330827, 0.06954887218045112, 
    # 0.02725563909774436, 0.03853383458646616, 0.042293233082706765, 0.06954887218045112, 0.05639097744360902, 
    # 0.06015037593984962, 0.09680451127819549],,,,]
    return encodings
