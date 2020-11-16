#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-03-04

# This python script is to generate only protein id and protein sequence for all yeast species

import json


def main() :
    organisms = ["Y_lipolytica", "S_pombe", "S_cerevisiae", "P_pastoris", "C_albicans", "R_toruloides"]
    id_proseq = dict()
    for organism in organisms :
        with open("../complementaryData/processed_data/%s.json" % organism) as f :
            data = json.load(f)

        for subdata in data :
            geneid = list(subdata.keys())[0]
            id_proseq[geneid] = subdata["protein sequence"]

    with open("../complementaryData/processed_data/allProtein.json", "w") as outfile :
        json.dump(id_proseq, outfile, indent=4)


if __name__ == "__main__" :
    main()


