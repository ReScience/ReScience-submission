from __future__ import print_function
import sys
import platform

print("Platform:", sys.platform)
print("Python:", sys.version)
print("Machine and architecture", platform.machine(), *platform.architecture())

try:
    import numpy as np
    print("NumPy:", np.version.version)
except:
    pass

try:
    import scipy
    print("SciPy:", scipy.version.version)
except:
    pass

try:
    import matplotlib
    print("matplotlib:", matplotlib.__version__)
except:
    pass

try:
    import chaospy
    print("Chaospy:", chaospy.__version__)
except:
    pass

try:
    import uncertainpy
    print("Uncertainpy:", uncertainpy.__version__)
except:
    pass

try:
    import neuron
    print("NEURON:", neuron.__version__)
except:
    pass

try:
    import tqdm
    print("tqdm:", tqdm.__version__)
except:
    pass