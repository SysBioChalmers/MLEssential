#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2019-12-23

# Plot forward and reverse blast results and reciprocal blast best hits results 
# and then do id mapping for the protein id obtained from my processed data and the protein id obtained from that Cell paper


import os
import json
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# Import Biopython tools for running local BLASTX
from Bio.Blast.Applications import NcbiblastpCommandline
# Colour scale transformation
from matplotlib.colors import LogNorm


def create_blast_command(strain) :
    # strains = ["Yarrowia_lipolytica", "Schizosaccharomyces_pombe", "Saccharomyces_cerevisiae", "Komagataella_pastoris"]  # Komagataella pastoris (Pichia pastoris)
    file1 = os.path.join("../reciprocal_blast/", "%s_query.fasta" % strain)
    file2 = os.path.join("../reciprocal_blast/", "%s_ref.fasta" % strain)

    fwd_out = os.path.join("../reciprocal_blast/", "%s_fwd.tab" % strain)
    rev_out = os.path.join("../reciprocal_blast/", "%s_rev.tab" % strain)

    # file1 = os.path.join("../Data/reciprocal_blast/", "Candida_albicans_query.fasta")
    # file2 = os.path.join("../Data/reciprocal_blast/", "Candida_albicans_ref.fasta")

    # fwd_out = os.path.join("../Data/reciprocal_blast/", "Candida_albicans_fwd.tab")
    # rev_out = os.path.join("../Data/reciprocal_blast/", "Candida_albicans_rev.tab")

    # Create BLAST command-lines for forward and reverse BLAST searches
    # blastp -out ../Data/reciprocal_blast/Candida_albicans_fwd.tab -outfmt "6 qseqid sseqid pident qcovs qlen slen length bitscore evalue" 
    # -query ../Data/reciprocal_blast/Candida_albicans_query.fasta -max_target_seqs 1 -subject ../Data/reciprocal_blast/Candida_albicans_ref.fasta
    fwd_blastp = NcbiblastpCommandline(query=file1, subject=file2, out=fwd_out,
                                       outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
                                       max_target_seqs=3)

    rev_blastp = NcbiblastpCommandline(query=file2, subject=file1, out=rev_out,
                                       outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
                                       max_target_seqs=3)

    # Inspect command-lines
    print("FORWARD: %s" % fwd_blastp)
    print("REVERSE: %s" % rev_blastp)
    os.system(str(fwd_blastp))
    os.system(str(rev_blastp))

    return fwd_blastp, rev_blastp

# https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/05-blast_for_rbh.html#reciprocal
def blast_results(strain) :
    fwd_out = os.path.join("../reciprocal_blast/", "%s_fwd.tab" % strain)
    rev_out = os.path.join("../reciprocal_blast/", "%s_rev.tab" % strain)
    # fwd_out = os.path.join("../Data/reciprocal_blast/", "Candida_albicans_fwd.tab")
    # rev_out = os.path.join("../Data/reciprocal_blast/", "Candida_albicans_rev.tab")

    fwd_results = pd.read_csv(fwd_out, sep="\t", header=None)
    rev_results = pd.read_csv(rev_out, sep="\t", header=None)

    # Add headers to forward and reverse results dataframes
    headers = ["query", "subject", "identity", "coverage",
               "qlength", "slength", "alength",
               "bitscore", "E-value"]

    fwd_results.columns = headers
    rev_results.columns = headers

    # Create a new column in both dataframes: normalised bitscore
    fwd_results['norm_bitscore'] = fwd_results.bitscore/fwd_results.qlength
    rev_results['norm_bitscore'] = rev_results.bitscore/rev_results.qlength

    # Create query and subject coverage columns in both dataframes
    fwd_results['qcov'] = fwd_results.alength/fwd_results.qlength
    rev_results['qcov'] = rev_results.alength/rev_results.qlength
    fwd_results['scov'] = fwd_results.alength/fwd_results.slength
    rev_results['scov'] = rev_results.alength/rev_results.slength

    # Clip maximum coverage values at 1.0
    # fwd_results['qcov'] = fwd_results['qcov'].clip_upper(1)
    # rev_results['qcov'] = rev_results['qcov'].clip_upper(1)
    # fwd_results['scov'] = fwd_results['scov'].clip_upper(1)
    # rev_results['scov'] = rev_results['scov'].clip_upper(1)
    fwd_results['qcov'] = fwd_results['qcov'].clip(upper=1)
    rev_results['qcov'] = rev_results['qcov'].clip(upper=1)
    fwd_results['scov'] = fwd_results['scov'].clip(upper=1)
    rev_results['scov'] = rev_results['scov'].clip(upper=1)

    return fwd_results, rev_results

def plot_oneway(fwd_results, rev_results, strain) :
    # Plot 2D density histograms
    # Calculate 2D density histograms for counts of matches at several coverage levels
    (Hfwd, xedgesf, yedgesf) = np.histogram2d(fwd_results.qcov, fwd_results.scov, bins=20)
    (Hrev, xedgesr, yedgesr) = np.histogram2d(rev_results.qcov, rev_results.scov, bins=20)

    # Create a 1x2 figure array
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharex=True, sharey=True)

    # Plot histogram for forward matches
    im = axes[0].imshow(Hfwd, cmap=plt.cm.Blues, norm=LogNorm(),
                        extent=[xedgesf[0], xedgesf[-1], yedgesf[0], yedgesf[-1]],
                        origin='lower', aspect=1)
    axes[0].set_title("Forward")
    axes[0].set_xlabel("query")
    axes[0].set_ylabel("subject")

    # Plot histogram for reverse matches
    im = axes[1].imshow(Hrev, cmap=plt.cm.Blues, norm=LogNorm(),
                        extent=[xedgesr[0], xedgesr[-1], yedgesr[0], yedgesr[-1]],
                        origin='lower', aspect=1)
    axes[1].set_title("Reverse")
    axes[1].set_xlabel("query")
    axes[1].set_ylabel("subject")

    # Add colourbars
    fig.colorbar(im, ax=axes[0])
    fig.colorbar(im, ax=axes[1])

    plt.savefig("../Data/reciprocal_blast/figure/%s_oneway.png" % strain, dpi=400)

def plot_rbbh(fwd_results, rev_results, strain) : # reciprocal blast best hits
    # Merge forward and reverse results
    rbbh = pd.merge(fwd_results, rev_results[['query', 'subject']],
                    left_on='subject', right_on='query',
                    how='outer')

    # Discard rows that are not RBH
    rbbh = rbbh.loc[rbbh.query_x == rbbh.subject_y]

    # Group duplicate RBH rows, taking the maximum value in each column
    rbbh = rbbh.groupby(['query_x', 'subject_x']).max()

    # Plot 2D density histograms

    # Calculate 2D density histograms for counts of matches at several coverage levels
    # (H, xedges, yedges) = np.histogram2d(rbbh.qcov, rbbh.scov, bins=20)

    # # Create a 1x2 figure array
    # fig, ax = plt.subplots(1, 1, figsize=(6, 6), sharex=True, sharey=True)

    # # Plot histogram for RBBH
    # im = ax.imshow(H, cmap=plt.cm.Blues, norm=LogNorm(),
    #                  extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
    #                  origin='lower', aspect=1)
    # ax.set_title("RBBH")
    # ax.set_xlabel("query")
    # ax.set_ylabel("subject")

    # # Add colourbar
    # fig.colorbar(im, ax=ax)

    # plt.savefig("../figure/rbbh/%s_rbbh.png" % strain, dpi=400)

    return rbbh

def id_mapping(rbbh, strain) :
    rbbh2 = rbbh.loc[rbbh["identity"]>95.0]
    # for col in rbbh2.columns :
    #     print(col)
    gene_protein = rbbh2.set_index('subject_y')['query_y'].to_dict()
    print("Reciprocal blast best hits with Pidentity more than 95%%: %s" % len(gene_protein))

    with open("../json/%s.json" % strain) as f :
        proteinAll = json.load(f)

    print("The number of protein complex data: %s" % len(proteinAll))

    print("Reciprocal blast best hits: %s" % len(rbbh.index))  # including all the reciprocal blast best hits results

    gene_mapping = list()
    gene_unmapping = list()
    for gene in proteinAll :
        try :
            gene["id"] = gene_protein[list(gene.keys())[0]]
            gene_mapping.append(gene)
        except :
            gene_unmapping.append(list(gene.keys())[0])

    print("The number of gene mapping: %s" % len(gene_mapping))  # including all the reciprocal blast best hits results pidentity > 95%
    # print(gene_unmapping)


    with open("../processed_data/%s_include_id.json" % strain, "w") as outfile :
        json.dump(gene_mapping, outfile, indent=4)


    # rbbh.to_excel("../Data/reciprocal_blast/rbbh_Candida_albicans2.xlsx")
    # fwd_results.to_csv("../Data/reciprocal_blast/fwd_Candida_albicans.csv")
    # rev_results.to_csv("../Data/reciprocal_blast/rev_Candida_albicans.csv")
    # len(rbbh.index)

if __name__ == "__main__" :
    strains = ["Candida_albicans", "Candida_glabrata", "Candida_dubliniensis", "Candida_parapsilosis", "Candida_tropicalis", "Yarrowia_lipolytica", "Schizosaccharomyces_pombe", "Saccharomyces_cerevisiae"]
    for strain in strains[:1] : # Run Candida glabrata
        print("This strain is: %s" % strain.replace("_", " "))
        create_blast_command(strain) # get the linux command line for BLAST
        # # create_blast_command()
        # fwd_results, rev_results = blast_results(strain)
        # plot_oneway(fwd_results, rev_results, strain)
        # rbbh = plot_rbbh(fwd_results, rev_results, strain)
        # id_mapping(rbbh, strain)
        # print("-----------------------------------------")


# Results for ref fasta from step 2:
# The protein number of Candida_albicans is: 6207
# The protein number of Candida_glabrata is: 5143
# The protein number of Candida_dubliniensis is: 5949
# The protein number of Candida_parapsilosis is: 5280
# The protein number of Candida_tropicalis is: 5975
# The protein number of Yarrowia_lipolytica is: 6433
# The protein number of Schizosaccharomyces_pombe is: 5134
# The protein number of Saccharomyces_cerevisiae is: 5911
# [Finished in 90.5s]

# Results:

# This strain is: Candida albicans
# Reciprocal blast best hits with Pidentity more than 95%: 6188
# The number of protein complex data: 12421
# Reciprocal blast best hits: 6188
# The number of gene mapping: 6188
# -----------------------------------------
# [Finished in 12.1s]

# This strain is: Candida glabrata
# Reciprocal blast best hits with Pidentity more than 95%: 4903
# The number of protein complex data: 5311
# Reciprocal blast best hits: 4932
# The number of gene mapping: 4903
# -----------------------------------------
# [Finished in 9.8s]

# This strain is: Candida dubliniensis
# Reciprocal blast best hits with Pidentity more than 95%: 5380
# The number of protein complex data: 5935
# Reciprocal blast best hits: 5412
# The number of gene mapping: 5380
# -----------------------------------------
# This strain is: Candida parapsilosis
# Reciprocal blast best hits with Pidentity more than 95%: 5019
# The number of protein complex data: 5863
# Reciprocal blast best hits: 5042
# The number of gene mapping: 5019
# -----------------------------------------
# This strain is: Candida tropicalis
# Reciprocal blast best hits with Pidentity more than 95%: 5442
# The number of protein complex data: 6258
# Reciprocal blast best hits: 5489
# The number of gene mapping: 5442
# -----------------------------------------
# [Finished in 25.9s]

# This strain is: Yarrowia lipolytica
# Reciprocal blast best hits with Pidentity more than 95%: 6013
# The number of protein complex data: 6472
# Reciprocal blast best hits: 6096
# The number of gene mapping: 6013
# -----------------------------------------
# [Finished in 10.9s]

# This strain is: Schizosaccharomyces pombe
# Reciprocal blast best hits with Pidentity more than 95%: 5060
# The number of protein complex data: 5135
# Reciprocal blast best hits: 5063
# The number of gene mapping: 5060
# -----------------------------------------
# [Finished in 9.6s]

# This strain is: Saccharomyces cerevisiae
# Reciprocal blast best hits with Pidentity more than 95%: 5842
# The number of protein complex data: 6713
# Reciprocal blast best hits: 5842
# The number of gene mapping: 5842
# -----------------------------------------
# [Finished in 10.8s]
