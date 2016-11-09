import numpy as np
import matplotlib.pyplot as plot
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import h5py
import sys

# file_name = sys.argv[1]
# event = sys.argv[2]

# with h5py.File(file_name,'r') as hf:
#    d = hf.get(event)
#    data = np.array(data)

length = np.arange(72)+1
x = np.repeat(length,72)
y = np.tile(length,72)
a = np.arange(72**2)+1
   
grid = a.reshape(72, 72)

plot.imshow(grid, vmin=1, vmax=5184, extent=(-x.min(), x.max(), -y.min(), y.max()), 
interpolation='nearest', cmap=cm.winter, aspect="auto")
plot.grid()
plot.show()
