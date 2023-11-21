"""
Piecewise Polynomial Model Evaluation

This script defines functions to create and evaluate a piecewise polynomial function
loaded from a specified file (default: 'PPModel.mat'). The main functionality includes
the 'mkpp' function for creating a piecewise polynomial function and 'get_PPModel'
function for loading a predefined model from a file.

Author: Benjamin Fortuno
Contact: benjaminignacio.fortuno@mail.polimi.it

Dependencies:
- scipy
- numpy
- matplotlib

Usage:
- Run this script to load the piecewise polynomial model and evaluate it at specific points.
"""

import scipy.io as spio
import numpy as np
from matplotlib import pyplot as plt

def mkpp(breaks, coefs):
    """
    Create a piecewise polynomial function using the given breaks and coefficients.

    Parameters:
    - breaks: List or array of break points.
    - coefs: List or array of coefficients for each polynomial segment.

    Returns:
    - A function representing the piecewise polynomial.
    """

    def ppval(x):
        # Find the segment where x belongs
        segment = np.searchsorted(breaks, x, side='right') - 1
        segment = np.clip(segment, 0, len(coefs) - 1)

        # Evaluate the polynomial at x within the selected segment
        result = np.polyval(coefs[segment], x - breaks[segment])
        return result

    return ppval

def get_PPModel(filename='PPModel.mat'):
    ss = spio.loadmat(filename)
    ss = ss['ss'][0][0]
    breaks = ss[1][0]
    coefs = ss[2]

    pp = mkpp(breaks, coefs)
    return pp

if __name__ == '__main__':
    pp = get_PPModel()

    # Evaluate the piecewise polynomial at specific points
    points_to_evaluate = np.linspace(-55, 96, 200)
    y = [pp(x) for x in points_to_evaluate]

    plt.plot(points_to_evaluate, y)
    plt.show()
