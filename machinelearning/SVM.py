#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse
import math
import numpy as np
from sklearn import svm
from sklearn.model_selection import StratifiedKFold, GridSearchCV


def SVM_Classifier(X, y, indep=None, fold=5, batch=None, auto=False, kernel='rbf', degree=3, gamma='auto',
                          coef0=0, C=1.0):
    default_params = {'degree': degree, 'gamma': gamma, 'coef0': coef0, 'C': C}
    if auto:
        data = np.zeros((X.shape[0], X.shape[1] + 1))
        data[:, 0] = y
        data[:, 1:] = X
        np.random.shuffle(data)
        X1 = data[:, 1:]
        y1 = data[:, 0]
        parameters = {'kernel': ['linear'], 'C': [1, 15]} if kernel == 'linear' else {'kernel': [kernel],
                                                                                      'C': [1, 15],
                                                                                      'gamma': 2.0 ** np.arange(-10, 4)}
        optimizer = GridSearchCV(svm.SVC(probability=True), parameters)
        optimizer = optimizer.fit(X1[0:math.ceil(batch * X1.shape[0]), ],
                                  y1[0:math.ceil(batch * y1.shape[0]), ]) if batch else optimizer.fit(X, y)
        params = optimizer.best_params_
        default_params['C'] = params['C']
        if kernel != 'linear':
            default_params['gamma'] = params['gamma']

    classes = sorted(list(set(y)))
    svms = []
    cvs = np.zeros((X.shape[0], len(classes) + 1))
    folds = StratifiedKFold(fold).split(X, y)

    inds = np.array([])
    if indep.shape[0] != 0:
        inds = np.zeros((len(indep), len(classes) + 1))
        inds[:, 0] = indep[:, 0]

    prediction_result_cv = []
    for trained, valided in folds:
        train_y, train_X = y[trained], X[trained]
        valid_y, valid_X = y[valided], X[valided]
        model = svm.SVC(C=default_params['C'], kernel=kernel, degree=default_params['degree'],
                        gamma=default_params['gamma'], coef0=default_params['coef0'], probability=True,
                        random_state=1)
        svc = model.fit(train_X, train_y)
        svms.append(svc)
        proba_ = svc.predict_proba(valid_X)
        cvs[valided, 0], cvs[valided, 1:] = valid_y, proba_
        # save the sample label and prediction result to ndarray
        tmp_result = np.zeros((len(valid_y), len(classes) + 1))
        tmp_result[:, 0], tmp_result[:, 1:] = valid_y, proba_
        prediction_result_cv.append(tmp_result)

        # independent
        if indep.shape[0] != 0:
            inds[:, 1:] += svc.predict_proba(indep[:, 1:])

    header = 'C=%f\tgamma=%s' % (default_params['C'], default_params['gamma'])
    if indep.shape[0] != 0:
        inds[:, 1:] /= fold
    return header, prediction_result_cv, inds

