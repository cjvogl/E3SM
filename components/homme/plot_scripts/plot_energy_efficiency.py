from dictionaries import create_energy_error_dict, create_walltime_dict, \
                         methodDict, lineStyleDict, colorDict
import matplotlib.pyplot as pyplot
import matplotlib
from netCDF4 import Dataset
import numpy as np
import os

matplotlib.rcParams.update({'font.size':22})

# Iterate through methods and plot energy efficiency test
f, ax = pyplot.subplots(figsize=(10,10))
f.set_tight_layout(True)
for m,method in enumerate(methodDict.keys()):
  print method
  globstr = 'tsteptype%d_tstep*.out' % methodDict[method]
  relErrorDict = create_energy_error_dict(globstr, 0, np.inf, '', 'restart')
  walltimeDict = create_walltime_dict(globstr, 0, np.inf, '', 'restart')

  # if no solutions were found, goto next method
  dtList = relErrorDict.keys()
  if (len(dtList) == 0):
    continue

  # need dictionary to map strings to floats
  dtDict = {}
  for dt in dtList:
    dtDict[float(dt)] = dt

  # Sort values and add to plot
  dtPlot = np.sort(dtDict.keys())
  relErrorPlot = np.empty(len(dtPlot))
  walltimePlot = np.empty(len(dtPlot))
  for j,dt in enumerate(dtPlot):
    relErrorPlot[j] = abs(relErrorDict[dtDict[dt]])
    walltimePlot[j] = walltimeDict[dtDict[dt]]

  ax.loglog(walltimePlot, relErrorPlot, lineStyleDict[method], label=method, \
            color=colorDict[method], linewidth=2, markersize=12)

ax.set_ylabel('Maximum Relative Energy Error', fontsize='xx-large')
ax.set_xlabel('wall time (s)', fontsize='xx-large')

#f.savefig('efficiencyEnergy.png')

ax.legend(loc='best')
#f.savefig('efficiencyEnergy_legend.png')
pyplot.show()
