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
from os import path


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
    local = path.dirname(path.realpath(__file__))
    filename = path.join(local, filename)
    ss = spio.loadmat(filename)
    ss = ss['ss'][0][0]
    breaks = ss[1][0]
    coefs = ss[2]

    pp = mkpp(breaks, coefs)
    return pp

def hysteresis_model(filename1='antero_posterior_PPModel.mat', filename2='postero_anterior_PPModel.mat'):
    pp1 = get_PPModel(filename1)
    pp2 = get_PPModel(filename2)

    def ppval(x, x_prev):
        if x < x_prev:          # Antero-posterior: if x_prev is greater than x, then the probe is moving backward
            return pp1(x)
        elif x > x_prev:        # Postero-anterior: if x_prev is less than x, then the probe is moving forward
            return pp2(x)
        else:
            raise ValueError('x and x_prev cannot be equal')

    return ppval

if __name__ == '__main__':
    pp = hysteresis_model()

    # Evaluate the piecewise polynomial at specific points
    x = np.linspace(-54, 97, 10000)
    x_min = np.linspace(-55, 96, 10000)
    x_plus = np.linspace(-53, 98, 10000)
    y_pa = [pp(x[i],x_min[i]) for i in range(x.size)]
    y_ap = [pp(x[i],x_plus[i]) for i in range(x.size)]

    # Plotting
    plt.figure(figsize=(8, 6))

    plt.plot(x, y_ap, label='Postero-anterior', linestyle='-', color='blue')
    plt.plot(x, y_pa, label='Antero-posterior', linestyle='-', color='orange')

    plt.title('TEE-Probe Hysteresis Model')
    plt.xlabel('$\\theta$ (degrees)', fontsize=12)
    plt.ylabel('steps', fontsize=12)
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.savefig('hysteresis_plot.png')  # Save the plot as an image file
    plt.show()
