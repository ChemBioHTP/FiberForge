import pandas as pd
import numpy as np
from scipy.interpolate import CubicSpline

def estimate_elastic_modulus(stress, strain):
    """
    Estimate the elastic modulus from a stress-strain curve.
    
    Parameters
    ----------
    stress : array-like
        The stress values.
    strain : array-like
        The strain values.
    
    Returns
    -------
    float
        The elastic modulus.
    """
    window_size = 400
    df = pd.DataFrame({'strain': strain, 'stress': stress})
    moving_avg = df.rolling(window=window_size).mean()
    moving_avg = moving_avg.dropna()
    

    args = np.argsort(moving_avg['strain']).values
    x = moving_avg['strain'].values[args]
    y = moving_avg['stress'].values[args]

    # Make x strictly increasing
    x, idx = np.unique(x, return_index=True)
    y = y[idx]

    # Create a cubic spline
    cs = CubicSpline(x, y)

    # Calculate the first derivative of the spline
    cs_derivative = cs.derivative()

    # Create a fine grid of x values for smooth plotting
    x_fine = np.linspace(x.min(), x.max(), 2000)
    y_fine = cs(x_fine)

    dy2_dx2_fine = cs_derivative.derivative()(x_fine)

    inflection_points = np.isclose(dy2_dx2_fine, 0.0, atol=1e11) == 1

    E = (y_fine[inflection_points][0] - y_fine[0]) / (x_fine[inflection_points][0] - x_fine[0])

    yield_point = x_fine[inflection_points][0]

    return E, yield_point

def calculate_variable_over_time(xvg_file):
    with open(xvg_file, 'r') as f:
        lines = f.readlines()
        time_data = []
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                time_data.append(list(map(float, line.split())))
    return time_data