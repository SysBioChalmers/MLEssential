#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse
import sys, os, re
import numpy as np
from pubscripts import *
from clusters import *
from dimreduction import lda as dimlda
from dimreduction import pca as dimpca
from dimreduction import tsne as dimtsne
# from featurenormalization import *
# from featureselection import *
# from machinelearning import *
from itertools import combinations
from descnucleotide import *

if __name__ == '__main__':

    parameters = {
        'Sequence_Type': 'DNA',
        # 'Clustering_Type': "sample",
        'Clustering_Type': "feature",
        'Kmean_Cluster_Number': 2,
        'Clustering_Algorithm': "kmeans"
    }

    kw = {'nclusters': 3, 'sof': 'sample', 'order': ''}
    kw['order'] = 'ACGT' if parameters['Sequence_Type'] == 'DNA' or parameters['Sequence_Type'] == 'RNA' else 'ACDEFGHIKLMNPQRSTVWY'
    # clustering for training data
    kw['sof'] = parameters['Clustering_Type'] if parameters['Clustering_Type'] != '' else 'sample'
    kw['nclusters'] = parameters['Kmean_Cluster_Number'] if parameters['Kmean_Cluster_Number'] else 2
    if parameters['Clustering_Algorithm'] in ('kmeans', 'hcluster', 'apc', 'meanshift', 'dbscan'):
        cluster_method = parameters['Clustering_Algorithm'].strip()
        testing_clustering_data = read_code.read_tsv_1('testing_code_1.tsv')
        cmd = cluster_method + '.' + cluster_method + '(testing_clustering_data, **kw)'
        clusters_res, e = eval(cmd)
        save_file.save_cluster_result(clusters_res, e, cluster_method + '.txt')
        draw_plot.plot_clustering_2d(testing_clustering_data, clusters_res, cluster_method, **kw)
        if e != '':
            error_array.append(e)


