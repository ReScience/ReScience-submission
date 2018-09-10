import matplotlib as mp
from matplotlib.colors import LinearSegmentedColormap

"""
Script to define the custom colormap for the connection experiments.
"""

def myCmap():
    cdict1 = {'red': ((0.0, 0.7, 0.1),
                     (0.4, 0.7, 1.0),
                     (1.0, 0.8, 0.8),
                     (1.0, 0.6, 0.0)),

           'green': ((0.0, 0.7, 0.7),
                     (0.4, 0.2, 1.0),
                     (1.0, 0.0, 0.0),
                     (1.0, 0.0, 0.0)),

           'blue':  ((0.0, 0.2, 0.9),
                     (0.1, 0.0, 0.1),
                     (1.0, 0.0, 0.0),
                     (1.0, 0.0, 0.0))
          }
    blue_red1 = LinearSegmentedColormap('BlueRed1', cdict1)
    return(blue_red1)

