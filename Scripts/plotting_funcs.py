"""
Some functions for plotting
"""
import matplotlib
import matplotlib.colors as mcolors
import numpy as np

def contour_levels_func(min_contour_level, max_contour_level, levels):
    """Function to define contour levels for contourf"""
    distance_levels = max_contour_level / levels
    contour_levels = np.arange(min_contour_level, max_contour_level, distance_levels)
    return contour_levels