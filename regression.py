"""Regressions based on linear model"""
import csv
import re
import sys

import numpy as np
import pandas as pd
from scipy import special


def linear_regression(X, y, sample_weight=None, fit_intercept=True):
    """
    Ordinary least squares linear regression.

    If fit_intercept is true, the last elements in stdev, t, and p_t correspond
    to the intercept.
    """
    pd.set_option('precision', 15)
    X = np.array(X, copy=True)
    n = X.shape[0]
    p = X.shape[1]
    y = np.array(y, copy=True).reshape(n, -1)
    n_targets = y.shape[1]
    dof = n - p

    # Preprocess data.
    if fit_intercept:
        X = np.append(X, np.ones((n, 1)), axis=1)
        dof = dof - 1
    if sample_weight is not None:
        sample_weight = np.sqrt(np.array(sample_weight, copy=True))
        try:
            sample_weight = sample_weight.reshape(n, 1)
            X = np.multiply(X, sample_weight)
            y = np.multiply(y, sample_weight)
        except:
            pass

    # Evaluate coefficients and standard deviation.
    XT = X.T
    XTX_inv = np.linalg.inv(np.dot(XT, X))
    beta = np.dot(np.dot(XTX_inv, XT), y)
    intercept = beta[-1] if fit_intercept else np.zeros(y.shape[1])
    coef = beta[0:-1] if fit_intercept else beta
    residual = y - np.dot(X, beta)

    diag = np.diagonal(XTX_inv)
    var_beta = np.outer(diag, np.sum(residual**2, axis=0)) / dof
    stdev = np.sqrt(var_beta).reshape(beta.shape)
    t = np.divide(beta, stdev)
    x = dof / (dof + t**2)
    prob_t = special.betainc(0.5 * dof, 0.5, x)  # Two-sided t test.

    if n_targets == 1:
        coef = coef.reshape(p,)
        intercept = intercept[0]
        residual = residual.reshape(n,)
        dof_regression = n - dof
        stdev = stdev.reshape(dof_regression,)
        t = t.reshape(dof_regression,)
        prob_t = prob_t.reshape(dof_regression,)

    return {
        'coef': coef,
        'intercept': intercept,
        'residual': residual,
        'stdev': stdev,
        't': t,
        'prob_t': prob_t
    }


def quadratic_regression(
        X, y, sample_weight=None, fit_intercept=True, cross_term=True):
    X = np.array(X, copy=True)
    n = X.shape[0]
    p = X.shape[1]
    y = np.array(y, copy=True).reshape(n, -1)

    if cross_term:
        for i in range(p):
            for j in range(i, p):
                new_col = np.multiply(X[:, i], X[:, j]).reshape(n, 1)
                X = np.append(X, new_col, axis=1)
    else:
        for i in range(p):
            new_col = np.multiply(X[:, i], X[:, i]).reshape(n, 1)
            X = np.append(X, new_col, axis=1)

    return linear_regression(X, y, sample_weight, fit_intercept)
