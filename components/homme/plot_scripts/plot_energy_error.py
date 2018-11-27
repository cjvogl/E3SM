from dictionaries import create_energy_error_dict, methodDict, lineStyleDict, \
                         colorDict
import matplotlib.pyplot as pyplot
import matplotlib
from netCDF4 import Dataset
import math
import numpy as np
import os

matplotlib.rcParams.update({'font.size':22})

dtTickList = [10, 20, 50, 100, 160, 300]

# Iterate through methods and plot energy test
f1, ax1 = pyplot.subplots(2,1,figsize=(10,10))
f2, ax2 = pyplot.subplots(figsize=(10,10))
f1.set_tight_layout(True)
f2.set_tight_layout(True)
for m,method in enumerate(methodDict.keys()):
  print method
  globstr = 'tsteptype%d_tstep*.out' % methodDict[method]
  relErrorDict = create_energy_error_dict(globstr, 0, np.inf, '', 'restart')

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
  for j,dt in enumerate(dtPlot):
    relErrorPlot[j] = relErrorDict[dtDict[dt]]

  if (relErrorPlot[-1] > 0):
    ind = 0
  else:
    ind = 1
  ax1[ind].plot(dtPlot, relErrorPlot, lineStyleDict[method], label=method, \
                color=colorDict[method], linewidth=3, markersize=12)
  ax2.loglog(dtPlot, abs(relErrorPlot), lineStyleDict[method], label=method, \
                color=colorDict[method], linewidth=3, markersize=12)

for ind in range(2):
  ax1[ind].set_ylabel('Maximum Relative Energy Error', fontsize='xx-large')
  ax1[ind].set_xlabel('dt (s)', fontsize='xx-large')
  ax1[ind].legend(loc='best')

ax2.set_ylabel('Maximum Relative Energy Error', fontsize='xx-large')
ax2.set_xlabel('dt (s)', fontsize='xx-large')
#ax2.axis('equal')

ax2.set_xticks(dtTickList)
ax2.set_xticklabels(dtTickList)

#f.savefig('errorEnergy.png')

ax2.legend(loc='best')
#f.savefig('errorEnergy_legend.png')
pyplot.show()
