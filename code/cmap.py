import matplotlib as mp
from matplotlib.colors import LinearSegmentedColormap

def myCmap():
    cdict1 = {'red':   ((0.0, 0.7, 0.1),
                     (0.4, 0.7, 1.0),
                     (1.0, 0.8, 0.8)),

           'green': ((0.0, 0.7, 0.1),
                     (0.4, 0.7, 1.0),
                     (1.0, 0.0, 0.0)),

           'blue':  ((0.0, 0.8, 0.8),
                     (0.1, 0.0, 0.1),
                     (1.0, 0.0, 0.0))
          }
    blue_red1 = LinearSegmentedColormap('BlueRed1', cdict1)
    return(blue_red1)
