
""" Obtain hci results from output and save as csv."""
import csv
import re
import sys

import numpy as np
import pandas as pd
from scipy import special

import regression

np.set_printoptions(precision=12)

POLYNOMIAL_ORDER = 2

def main():
    """main function"""
    # Check and read res file.
    if len(sys.argv) != 2:
        raise SyntaxError("Usage: extrapolate.py [res_file]")

    parameters = ['n_orbs_var_inv', 'eps_var', 'n_orbs_pt_inv', 'eps_pt']

    # Estimate intercept.
    data = pd.read_csv(sys.argv[1])
    for parameter in ['n_orbs_var', 'n_orbs_pt']:
        data[parameter + '_inv'] = 1.0 / data[parameter]
    if POLYNOMIAL_ORDER == 1:
        initial_model = regression.linear_regression(
            data[parameters].values, data['energy_corr'].values
        )
    else:
        initial_model = regression.quadratic_regression(
            data[parameters].values, data['energy_corr'].values
        )
    initial_intercept = initial_model['intercept']
    print('Intercept estimation: ' + str(initial_intercept))

    # Obtain weight.
    data['uncert'] = 0
    initial_coef = initial_model['coef']
    for i in range(4):
        contribution = data[parameters[i]] * initial_coef[i]
        data['uncert'] = data['uncert'] + contribution**2
    data['weight'] = 1.0 / data['uncert']
    data['weight'] = data['weight'] / np.linalg.norm(data['weight'].values)

    # Weighted regression.
    if POLYNOMIAL_ORDER == 1:
        final_model = regression.linear_regression(
            data[parameters].values,
            data['energy_corr'].values,
            data['weight'].values
        )
    else:
        final_model = regression.quadratic_regression(
            data[parameters].values,
            data['energy_corr'].values,
            data['weight'].values
        )

    final_intercept = final_model['intercept']
    final_stdev = final_model['stdev']

    print('Intercept: ' + str(final_intercept) + ' +- ' + str(final_stdev[-1]))
    print('coef, stdev, prob_t')
    print(np.hstack((
        np.append(final_model['coef'], final_intercept).reshape(-1, 1),
        final_model['stdev'].reshape(-1, 1),
        final_model['prob_t'].reshape(-1, 1)
    )))


if __name__ == '__main__':
    main()
