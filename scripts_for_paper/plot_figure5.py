
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
from netCDF4 import Dataset
import numpy as np

matplotlib.rcParams.update({'font.size': 24})

measureList = [120, 75, 30, 15, 8]
timeList = [10800, 21600, 32400, 43200]
markerDict = {10800: '^',
              21600: 's',
              32400: 'o',
              43200: 'p'}
colorDict = {10800: 'red',
             21600: 'green',
             32400: 'blue',
             43200: 'black'}

fileDict = {'SGR-P[0,0,0]': './data/SGR_P000_RKZ_qlneg_lm4_1069852.txt',
            'SGR-P[1,1,1]': './data/SGR_P111_RKZ_qlneg_lm4_1069852.txt'}
linespecDict = {'SGR-P[0,0,0]': '--g',
                'SGR-P[1,1,1]': '-b'}


f, ax = pyplot.subplots(1, 2, figsize=(20,10))

# plot convergence of SGR-P[1,1,1]
for time in timeList:
  tmp = np.loadtxt('./data/SGR_P111_convergence_data_%d.txt' % time)
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
  ax[0].loglog(dtPlot, WRMSPlot, markerDict[time], color=colorDict[time], linewidth=3, markersize=12)
  ax[0].loglog(dtPlot, coeff*np.power(dtPlot,order), ':', color=colorDict[time], linewidth=3, markersize=12)
  ax[0].loglog(measurePlot, coeff*np.power(measurePlot,order), '-', color=colorDict[time], label=label, linewidth=3, markersize=12)
  print(time, order)

dtPlot = np.array([8, 1800])
ax[0].loglog(dtPlot, 4e-2/dtPlot[-1]*dtPlot, '--k', linewidth=3)
ax[0].set_ylabel('WRMS temperature error (K)', fontsize='x-large')
ax[0].set_xlabel('timestep size (s)', fontsize='x-large')
ax[0].legend(loc='upper left', fontsize='large')
ax[0].axis('equal')

# plot clipping in SGR-P[0,0,0] vs SGR-P[1,1,1]
for sgr in fileDict.keys():
  tmp = np.loadtxt(fileDict[sgr])
  t = tmp[0,:]
  val = tmp[1,:]
  ax[1].plot(t[2:], val[2:], linespecDict[sgr], label=sgr, linewidth=3)

ax[1].set_ylabel('clipped concentration (kg/kg)', fontsize='x-large')
ax[1].set_xlabel('time (s)', fontsize='x-large')
ax[1].legend(loc='best', fontsize='large')

# write out figure
f.tight_layout()
f.savefig('figure5-P111.pdf')


pyplot.show()
