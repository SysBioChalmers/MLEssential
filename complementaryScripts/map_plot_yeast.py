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
    file1 = os.path.join("../Data/reciprocal_blast/", "%s_query.fasta" % strain)
    file2 = os.path.join("../Data/reciprocal_blast/", "%s_ref.fasta" % strain)

    fwd_out = os.path.join("../Data/reciprocal_blast/", "%s_fwd.tab" % strain)
    rev_out = os.path.join("../Data/reciprocal_blast/", "%s_rev.tab" % strain)

    # file1 = os.path.join("../Data/reciprocal_blast/", "Candida_albicans_query.fasta")
    # file2 = os.path.join("../Data/reciprocal_blast/", "Candida_albicans_ref.fasta")

    # fwd_out = os.path.join("../Data/reciprocal_blast/", "Candida_albicans_fwd.tab")
    # rev_out = os.path.join("../Data/reciprocal_blast/", "Candida_albicans_rev.tab")

    # Create BLAST command-lines for forward and reverse BLAST searches
    # blastp -out ../Data/reciprocal_blast/Candida_albicans_fwd.tab -outfmt "6 qseqid sseqid pident qcovs qlen slen length bitscore evalue" 
    # -query ../Data/reciprocal_blast/Candida_albicans_query.fasta -max_target_seqs 1 -subject ../Data/reciprocal_blast/Candida_albicans_ref.fasta
    fwd_blastp = NcbiblastpCommandline(query=file1, subject=file2, out=fwd_out,
                                       outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
                                       max_target_seqs=1)

    rev_blastp = NcbiblastpCommandline(query=file2, subject=file1, out=rev_out,
                                       outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
                                       max_target_seqs=1)

    # Inspect command-lines
    print("FORWARD: %s" % fwd_blastp)
    print("REVERSE: %s" % rev_blastp)

    return fwd_blastp, rev_blastp

# https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/05-blast_for_rbh.html#reciprocal
def blast_results(strain) :
    fwd_out = os.path.join("../Data/reciprocal_blast/", "%s_fwd.tab" % strain)
    rev_out = os.path.join("../Data/reciprocal_blast/", "%s_rev.tab" % strain)
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
    (H, xedges, yedges) = np.histogram2d(rbbh.qcov, rbbh.scov, bins=20)

    # Create a 1x2 figure array
    fig, ax = plt.subplots(1, 1, figsize=(6, 6), sharex=True, sharey=True)

    # Plot histogram for RBBH
    im = ax.imshow(H, cmap=plt.cm.Blues, norm=LogNorm(),
                     extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                     origin='lower', aspect=1)
    ax.set_title("RBBH")
    ax.set_xlabel("query")
    ax.set_ylabel("subject")

    # Add colourbar
    fig.colorbar(im, ax=ax)

    plt.savefig("../Data/reciprocal_blast/figure/%s_rbbh.png" % strain, dpi=400)

    return rbbh

def id_mapping(rbbh, organism) :
    rbbh2 = rbbh.loc[rbbh["identity"]>95.0]
    # for col in rbbh2.columns :
    #     print(col)
    gene_protein = rbbh2.set_index('subject_y')['query_y'].to_dict()
    print("Reciprocal blast best hits with Pidentity more than 95%%: %s" % len(gene_protein))

    with open("../data/processed_data/%s.json" % organism) as f :
        geneAll = json.load(f)

    print("The number of gene essentiality data: %s" % len(geneAll))

    print("Reciprocal blast best hits: %s" % len(rbbh.index))

    gene_mapping = list()
    gene_unmapping = list()
    for gene in geneAll :
        try :
            gene["id"] = gene_protein[list(gene.keys())[0]]
            gene_mapping.append(gene)
        except :
            gene_unmapping.append(list(gene.keys())[0])

    print("The number of gene mapping: %s" % len(gene_mapping))
    print(gene_unmapping)


    with open("../Data/processed_data/%s_include_id.json" % organism, "w") as outfile :
        json.dump(gene_mapping, outfile, indent=4)


    # rbbh.to_excel("../Data/reciprocal_blast/rbbh_Candida_albicans2.xlsx")
    # fwd_results.to_csv("../Data/reciprocal_blast/fwd_Candida_albicans.csv")
    # rev_results.to_csv("../Data/reciprocal_blast/rev_Candida_albicans.csv")
    # len(rbbh.index)

if __name__ == "__main__" :
    strains = ["Yarrowia_lipolytica", "Schizosaccharomyces_pombe", "Saccharomyces_cerevisiae", "Komagataella_pastoris", "Candida_albicans"]  # Komagataella pastoris (Pichia pastoris)
    organisms = {"Yarrowia_lipolytica" : "Y_lipolytica", "Schizosaccharomyces_pombe" : "S_pombe", "Saccharomyces_cerevisiae" : "S_cerevisiae", "Komagataella_pastoris" : "P_pastoris", "Candida_albicans" : "C_albicans"}
    for strain in strains :
        print("This strain is: %s" % strain)
    #     create_blast_command(strain)
        # create_blast_command()
        fwd_results, rev_results = blast_results(strain)
        plot_oneway(fwd_results, rev_results, strain)
        rbbh = plot_rbbh(fwd_results, rev_results, strain)
        id_mapping(rbbh, organisms[strain])
        print("-----------------------------------------")


# Results:

# This strain is: Yarrowia_lipolytica
# Reciprocal blast best hits with Pidentity more than 95%: 627
# The number of gene essentiality data: 641
# Reciprocal blast best hits: 637
# The number of gene mapping: 627
# ['YALI0D22396g', 'YALI0B13552g', 'YALI0F19448g', 'YALI0D15598g', 'YalifMp01', 'YalifMp02', 'YalifMp16', 'YalifMp19', 'YalifMp20', 'YalifMp28', 'YalifMp29', 'YalifMp03', 'YalifMp05', 'YalifMp06']
# -----------------------------------------
# This strain is: Schizosaccharomyces_pombe
# Reciprocal blast best hits with Pidentity more than 95%: 4566
# The number of gene essentiality data: 4605
# Reciprocal blast best hits: 4567
# The number of gene mapping: 4566
# ['SPBC1348.06c', 'SPBC18H10.14', 'SPAC11G7.04', 'SPAC1420.03', 'SPAC144.11', 'SPAC1834.03c', 'SPAC1834.04', 'SPAC19B12.04', 'SPAC19B12.12c', 'SPAC1F7.13c', 'SPAC212.01c', 'SPAC212.11', 'SPAC22A12.04c', 'SPAC22F3.11c', 'SPAC23A1.10', 'SPAC23C11.02c', 'SPAC25G10.06', 'SPAC26A3.04', 'SPAC26A3.07c', 'SPAC2F3.12c', 'SPAC3F10.18c', 'SPAC3H5.05c', 'SPAC513.01c', 'SPAC750.02c', 'SPAC977.01', 'SPAC977.07c', 'SPAC977.08', 'SPAC977.09c', 'SPBC1105.11c', 'SPBC1105.12', 'SPBC1348.01', 'SPBC1348.02', 'SPBC16D10.11c', 'SPBC2F12.07c', 'SPBPB10D8.04c', 'SPBPB10D8.05c', 'SPBPB10D8.06c', 'SPCC16C4.13c', 'SPCC622.17']
# -----------------------------------------
# This strain is: Saccharomyces_cerevisiae
# Reciprocal blast best hits with Pidentity more than 95%: 5538
# The number of gene essentiality data: 5577
# Reciprocal blast best hits: 5538
# The number of gene mapping: 5538
# ['YOL120C', 'YAL068C', 'YAR010C', 'YBL027W', 'YBL087C', 'YBR009C', 'YBR010W', 'YBR048W', 'YBR181C', 'YDL083C', 'YDL133C-A', 'YDL136W', 'YDR316W-A', 'YDR365W-A', 'YDR385W', 'YDR418W', 'YDR450W', 'YER074W', 'YER109C', 'YER137C-A', 'YER159C-A', 'YFL002W-A', 'YFL002W-B', 'YFR031C-A', 'YGL135W', 'YGR118W', 'YGR161C-C', 'YHR141C', 'YHR203C', 'YIL148W', 'YJL138C', 'YJR026W', 'YJR094W-A', 'YLR157C-A', 'YLR227W-A', 'YLR287C-A', 'YMR242C', 'YOR031W', 'YOR142W-A']
# -----------------------------------------
# This strain is: Komagataella_pastoris
# Reciprocal blast best hits with Pidentity more than 95%: 588
# The number of gene essentiality data: 608
# Reciprocal blast best hits: 598
# The number of gene mapping: 588
# ['PAS_chr1-4_0423', 'PAS_chr4_0780', 'PAS_chr3_0942', 'PAS_chr3_0676', 'PAS_chr1-4_0303', 'PAS_chr1-3_0211', 'PAS_chr1-1_0389', 'PAS_chr3_0788', 'PAS_chr3_0271', 'PAS_chr3_0597', 'PAS_chr3_0631', 'PAS_chr1-1_0392', 'PAS_chr1-1_0374', 'PAS_chr3_0674', 'PAS_chr4_0422', 'PAS_chr1-4_0447', 'PAS_chr3_1197', 'PAS_chr1-1_0424', 'PAS_chr3_0184', 'PAS_chr2-2_0111']
# -----------------------------------------
# This strain is: Candida_albicans
# Reciprocal blast best hits with Pidentity more than 95%: 2307
# The number of gene essentiality data: 2335
# Reciprocal blast best hits: 2325
# The number of gene mapping: 2307
# ['orf19.2629', 'orf19.5817', 'orf19.3001', 'orf19.2158', 'orf19.2292', 'orf19.2863.1', 'orf19.2666', 'orf19.3586', 'orf19.2124', 'orf19.4984', 'orf19.5221', 'orf19.2943', 'orf19.4769', 'orf19.153', 'orf19.5656', 'orf19.1079', 'orf19.6847', 'orf19.2241', 'orf19.4916', 'orf19.3362', 'orf19.785', 'orf19.7038', 'orf19.645', 'orf19.3659', 'orf19.3937', 'orf19.5895', 'orf19.4195', 'orf19.3365']
# -----------------------------------------
# [Finished in 40.2s]


