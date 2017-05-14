"""Test for regression module"""
import unittest
import regression
import numpy as np
import pandas as pd


class TestRegression(unittest.TestCase):

    def test_linear_regression(self):
        test_data = pd.DataFrame({
            'x': [6, 7, 8],
            'y': [8000, 50000, 116000],
            'weight': [123, 123, 246]})
        coef, intercept, residual, stdev, t, prob_t = regression.linear_regression(
            test_data.as_matrix(['x']),
            test_data['y'].values,
            test_data['weight'].values
        )
        np.testing.assert_array_almost_equal(coef, [55090.90909091], 5)
        np.testing.assert_almost_equal(intercept, -326909.09090909, 5)
        np.testing.assert_array_almost_equal(
            stdev, [6171.11372672, 45032.21987029], 5)
        np.testing.assert_array_almost_equal(prob_t, [0.071016, 0.087147], 5)

    def test_quadratic_regression_with_cross_term(self):
        test_data = pd.DataFrame({
            'x1': [4, 5, 6, 7, 8, 9, 10],
            'x2': [1, 2, 3, 5, 8, 13, 21],
            'y': [3000, 4000, 5000, 8000, 50000, 116000, 200000],
            'weight': [1, 2, 3, 4, 5, 6, 7]})
        coef, intercept, residual, stdev, t, prob_t = regression.quadratic_regression(
            test_data.as_matrix(['x1', 'x2']),
            test_data['y'].values,
            test_data['weight'].values
        )
        np.testing.assert_array_almost_equal(coef, [
            59926.609491270,  # x
            -63583.179647400,  # y
            -7019.859097691,  # x^2
            11165.359300379,  # xy
            -1372.826686406  # y^2
        ], 5)
        np.testing.assert_almost_equal(intercept, -101467.379938516, 5)
        np.testing.assert_array_almost_equal(stdev, [
            90078.711333994,
            108980.816420539,
            8836.035263482,
            11622.378213880,
            950.480040222,
            213425.391893574
        ], 5)
        np.testing.assert_array_almost_equal(prob_t, [
            0.62628, 0.66377, 0.57260, 0.51277, 0.38552, 0.71747
        ], 5)

    def test_quadratic_regression_no_cross_term(self):
        test_data = pd.DataFrame({
            'x1': [4, 5, 6, 7, 8, 9, 10],
            'x2': [1, 2, 3, 5, 8, 13, 21],
            'y': [3000, 4000, 5000, 8000, 50000, 116000, 200000],
            'weight': [1, 2, 3, 4, 5, 6, 7]})
        coef, intercept, residual, stdev, t, prob_t = regression.quadratic_regression(
            test_data.as_matrix(['x1', 'x2']),
            test_data['y'].values,
            test_data['weight'].values,
            cross_term=False
        )
        np.testing.assert_array_almost_equal(coef, [
            1622.948003888,  # x
            38512.336773585,  # y
            -4129.292382104,  # x^2
            -539.871870522  # y^2
        ], 5)
        np.testing.assert_almost_equal(intercept, 26143.556236498, 5)
        np.testing.assert_array_almost_equal(stdev, [
            65269.245149855,
            23666.371976750,
            8146.250299299,
            381.827788241,
            163797.792370775
        ], 5)
        np.testing.assert_array_almost_equal(prob_t, [
            0.98242, 0.24520, 0.66259, 0.29297, 0.88785
        ], 5)


if __name__ == '__main__':
    unittest.main()
