#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-03-12

import os


# To get new fasta file for 

def uniqueId(inputfile) :
    # for filename in os.listdir("../0_332yeast_genomes/332_genome_annotations/pep"):
    with open("../complementaryData/8_yeast_species/" + inputfile, 'r') as infile :
        lines = infile.readlines()

    with open("../complementaryData/generated_fasta/" + inputfile, 'w') as outfile :
        for line in lines :
            if line.startswith(">") :
                # print(line)             
                sequence_description = line.strip("\n")
                # print(type(sequence_description))
                # print(sequence_description)
                sequence = sequence_description.split(" ")
                # print(sequence)
                # print(sequence[0] + '-'+ sequence[2])
                # sequence[0] = sequence[0] + "-" + sequence[2]  # It cannot work either, because
                # FATAL: Sequence identifiers must be unique. Your fasta file contains two sequences 
                # with the same id (-processed-gene-1.26-mRNA-1_1-CDS=1-675) at pfam_scan.pl line 112.
                sequence[0] = sequence[0] + "&" + sequence[1] + "&" + sequence[2]
                sequence_description = " ".join(sequence)
                new_description = sequence_description + "\n"
                newline = line.replace(line,new_description)
                # print(new_description)
                # print(newline)
                outfile.write(newline)
            else :
                outfile.write(line)

            # if i == 5 :
            #     break

    # print("Running this scripts uses: %ss" %(elapsed))

def main() :
    i = 0
    for inputfile in os.listdir("../complementaryData/8_yeast_species") :
        i += 1
        print(i)
        print(inputfile[:-6])

        organisms = ["Candida_glabrata", "Candida_dubliniensis", "Candida_parapsilosis", "Candida_tropicalis", "Yarrowia_lipolytica"]
        if inputfile[:-6] in organisms :
            uniqueId(inputfile)
    print("Finished------------------------")  


if __name__ == "__main__" :
    main()
