from time import sleep
from matplotlib import pyplot, rcParams
import numpy as np
import os, sys
from netCDF4 import Dataset, num2date

rcParams['font.size'] = 8

field = 'hs'
cLevels = np.arange(0.5, 3.601, 0.5)
ncFile = sys.argv[1]


if not os.access(ncFile,(os.F_OK and os.R_OK)):
     raise IOError('Unable to access file \'%s\'.'%ncFile)

tri = None
rect = None
cart = None

with Dataset(ncFile) as nc:
    time = nc.variables['time']
    if 'tri' in nc.variables.keys():
        tri = nc.variables['tri'][:]
        # The offset -1 is required for Python indexing
        tri -= 1
    if 'longitude' in nc.variables.keys():
        y = nc.variables['latitude'][:]
        x = nc.variables['longitude'][:]
        rect = True
    if 'x' in nc.variables.keys():
        y = nc.variables['y'][:]
        x = nc.variables['x'][:]
        x *= 1.E-3
        y *= 1.E-3
        cart = True

    var = nc.variables['hs'][::3]
    T1 = num2date(time[::3], time.units)

ROWS = 3
COLS = 3

fig2, axs = pyplot.subplots(ROWS, COLS, figsize=(2.*COLS, 2.*ROWS),
    layout='constrained',
)

for i, ax in enumerate(axs.flatten()):
  if tri is not None:
      lc = ax.tricontour(
          x, y, tri, var[i,:], cLevels,
          colors="k",
          extend="both",
          zorder=-1,
      )
  else:
      lc = ax.contour(
          x, y, var[i,:], cLevels,
          colors="k",
          extend="both",
          zorder=-1,
      )

  ax.clabel(lc, fontsize=7)
  LABEL = T1[i].strftime("%Y/%m/%d %H:%M")

  if i<6:
      xycoords = (0.5, 0.96)
      _va="top"
  else:
      xycoords = (0.5, 0.04)
      _va="bottom"

  ax.annotate(
      LABEL,
      xycoords,
      xycoords="axes fraction",
      ha="center",
      va=_va,
      backgroundcolor="w",
)

pyplot.show()
