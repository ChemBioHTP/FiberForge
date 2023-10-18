#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

data = np.loadtxt('final_fes.dat')

x = np.unique(data[:, 0])
y = np.unique(data[:, 1])
xy = np.stack([np.tile(x, y.shape[0]), np.repeat(y, x.shape[0])]).T

grid = griddata(data[:, :2], data[:, 2], xy).reshape((x.shape[0], y.shape[0]))

plt.imshow(grid, extent=(x.min(), 1.5, y.max(), y.min()), aspect='auto', interpolation = 'bicubic', cmap='RdBu', vmin=-10, vmax=10)
plt.xlabel('Distance (nm)')
plt.ylabel('Order parameter')
cbar = plt.colorbar()
cbar.set_label('Free energy landscape (kJ/mol)')
#plt.show()

plt.savefig('fes.png')
