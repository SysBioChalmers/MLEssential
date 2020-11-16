#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold


def RF_Classifier(X, y, indep=None, fold=5, n_trees=100, out='RF_output'):
    """
    Parameters:
    ----------
    :param X: 2-D ndarray
    :param y: 1-D ndarray
    :param indep: 2-D ndarray, the first column is labels and the rest are feature values
    :param fold: int, default 5
    :param n_trees: int, number of trees, default: 5
    :param out:
    :return:
        info: str, the model parameters
        cross-validation result: list with element is ndarray
        independent result: ndarray, the first column is labels and the rest are prediction scores.
    """
    classes = sorted(list(set(y)))
    if indep.shape[0] != 0:
        indep_out = np.zeros((indep.shape[0], len(classes) + 1))
        indep_out[:, 0] = indep[:, 0]

    prediction_result_cv = []
    prediction_result_ind = np.array([])
    if indep.shape[0] != 0:
        prediction_result_ind = np.zeros((len(indep), len(classes) + 1))
        prediction_result_ind[:, 0] = indep[:, 0]

    folds = StratifiedKFold(fold).split(X, y)
    for i, (trained, valided) in enumerate(folds):
        train_y, train_X = y[trained], X[trained]
        valid_y, valid_X = y[valided], X[valided]
        model = RandomForestClassifier(n_estimators=n_trees, bootstrap=False)
        rfc = model.fit(train_X, train_y)
        scores = rfc.predict_proba(valid_X)
        tmp_result = np.zeros((len(valid_y), len(classes) + 1))
        tmp_result[:, 0], tmp_result[:, 1:] = valid_y, scores
        prediction_result_cv.append(tmp_result)
        # independent
        if indep.shape[0] != 0:
            prediction_result_ind[:, 1:] += rfc.predict_proba(indep[:, 1:])
    if indep.shape[0] != 0:
        prediction_result_ind[:, 1:] /= fold
    header = 'n_trees: %d' % n_trees
    return header, prediction_result_cv, prediction_result_ind

