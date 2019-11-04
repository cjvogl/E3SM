import glob
import matplotlib
from netCDF4 import Dataset
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
matplotlib.rcParams.update({'font.size': 24})

measureList = [120, 75, 30, 15, 8]
timeList = [3600, 7200, 10800]
markerDict = {3600: 'o',
              7200: 's',
              10800: '^'}
colorDict = {3600: 'blue',
             7200: 'green',
             10800: 'red'}

f, ax = pyplot.subplots(figsize=(10,10))
for time in timeList:
  tmp = np.loadtxt('./data/noSGR_convergence_data_%d.txt' % time)
  dtPlot = tmp[0,:]
  WRMSPlot = tmp[1,:]
  WRMSerror = {}
  for j,dt in enumerate(dtPlot):  
    WRMSerror[int(dt)] = WRMSPlot[j]
 
  measurePlot = np.array([measureList[-1], measureList[0]])
  N = len(measureList)
  sxlogy = 0.0
  slogy = 0.0
  sx2 = 0.0
  sx = 0.0
  for dt in measureList:
    logdt = np.log10(float(dt))
    sxlogy += logdt*np.log10(WRMSerror[dt])
    slogy += np.log10(WRMSerror[dt])
    sx2 += logdt*logdt
    sx += logdt

  order = (N*sxlogy - sx*slogy)/(N*sx2 - sx*sx)
  coeff = np.power(10.0, (sx2*slogy - sx*sxlogy)/(N*sx2 - sx*sx))

  label = "%0.0f hour run (%0.1f)" % (time/3600.0, order)
  ax.loglog(dtPlot, WRMSPlot, markerDict[time], color=colorDict[time], linewidth=3, markersize=12)
  ax.loglog(dtPlot, coeff*np.power(dtPlot,order), ':', color=colorDict[time], linewidth=3, markersize=12)
  ax.loglog(measurePlot, coeff*np.power(measurePlot,order), '-', color=colorDict[time], label=label, linewidth=3, markersize=12)
  print(time, order)

dtPlot = np.array([8, 1800])
ax.loglog(dtPlot, 4e-2/dtPlot[-1]*dtPlot, '--k', linewidth=3)
ax.set_ylabel('WRMS temperature error (K)', fontsize='x-large')
ax.set_xlabel('timestep size (s)', fontsize='x-large')
ax.legend(loc='upper left', fontsize='large')
ax.axis('equal')

f.tight_layout()
f.savefig('convergence_noSGR.pdf')  

  
pyplot.show()
