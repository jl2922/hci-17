""" Obtain hci results from output and save as csv."""
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm

np.set_printoptions(precision=12)


def printCorrelationEnergy(statsResult):
    coefs = statsResult.params.values
    stdevs = statsResult.bse
    print('Correlation Energy: ' + str(coefs[0]) + ' +- ' + str(stdevs[0]))


def BEWRegression(X, y):
    augX = sm.add_constant(X)

    # Backward elimination.
    print()
    print('*' * 80)
    print('Backward elimination:')
    while True:
        results = sm.OLS(y, augX).fit()
        print()
        # print(results.summary())
        printCorrelationEnergy(results)
        intercept = results.params.values[0]
        maxPIndex = np.argmax(results.pvalues)
        maxP = results.pvalues[maxPIndex]

        if maxP < 0.01:
            break
        print('Eliminate: ' + maxPIndex)
        print('P > |t|: ' + str(maxP))
        augX.drop(maxPIndex, axis=1, inplace=True)

    # Weighted OLS
    print('\n[FINAL Weighted OLS]')
    variance = np.square(
        np.dot(augX, np.abs(results.params.values)) + intercept)
    results = sm.WLS(y, augX, weights=1.0 / variance).fit()
    print(results.summary())
    printCorrelationEnergy(results)


def main():
    """main function"""
    # Check and read res file.
    if len(sys.argv) != 2:
        raise SyntaxError("Usage: extrapolate.py [res_file]")

    parameters = ['n_orbs_var_inv', 'eps_var', 'n_orbs_pt_inv', 'eps_pt']

    # Read raw data.
    data = pd.read_csv(sys.argv[1])

    # Add inverse terms.
    for parameter in ['n_orbs_var', 'n_orbs_pt']:
        data[parameter + '_inv'] = 1.0 / data[parameter]

    # Remove parameters not enough for extrapolation.
    for parameter in parameters:
        if data[parameter].value_counts().size < 3:
            parameters.remove(parameter)

    # Add cross terms.
    selectedParameters = parameters[:]
    for i in range(len(parameters)):
        for j in range(i, len(parameters)):
            column = parameters[i] + ' * ' + parameters[j]
            selectedParameters.append(column)
            data[column] = data[parameters[i]] * data[parameters[j]]

    # Estimate intercept.
    X = data[selectedParameters]
    y = data['energy_corr']

    BEWRegression(X, y)

    for i, parameter in enumerate(parameters):
        minValue = X.min()[i]
        keep = X[parameter] != minValue
        X = X[keep]
        y = y[keep]
    BEWRegression(X, y)

    for i, parameter in enumerate(parameters):
        minValue = X.min()[i]
        keep = X[parameter] != minValue
        X = X[keep]
        y = y[keep]
    BEWRegression(X, y)


if __name__ == '__main__':
    main()
