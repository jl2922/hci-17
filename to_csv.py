#!/usr/bin/python
""" Obtain hci results from output and save as csv."""
import csv
import re
import sys


def main():
    """main function"""
    # Check and read res file.
    if len(sys.argv) != 2:
        raise SyntaxError("Usage: to_csv.py [res_file]")
    with open(sys.argv[1], 'r') as res_file:
        res = res_file.read()

    print(','.join([
        'n_orbs_var', 'eps_var', 'n_orbs_pt', 'eps_pt', 'energy_pt', 'energy_corr']))
    pattern = '\n'.join([
        'BEGIN accumulate for eps_pt:.*',
        '.*',
        'n_orbs_var: (.*)',
        'eps_var: (.*)',
        'n_orbs_pt: (.*)',
        'eps_pt: (.*)',
        'Perturbation energy: (.*) Ha',
        'Correlation Energy: (.*) Ha'
    ])
    for match in re.finditer(re.compile(pattern), res):
        print(','.join(match.groups()))


if __name__ == '__main__':
    main()
