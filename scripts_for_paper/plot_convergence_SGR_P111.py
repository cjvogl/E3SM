import glob
import matplotlib.pyplot as pyplot
import matplotlib
from netCDF4 import Dataset
import numpy as np

matplotlib.rcParams.update({'font.size': 24})

dtList = [1800,450,120,75,30,15,8,4,1]
measureList = [120,75,30,15,8]

timeList = [10800, 21600, 32400, 43200]
markerDict = {10800: '^',
              21600: 's',
              32400: 'o',
              43200: 'p'}
colorDict = {10800: 'red',
             21600: 'green',
             32400: 'blue',
             43200: 'black'}

configuration = 'RKZ_SGR_extrapf_qv1_ql1_Al1_lmt4_adjIC'

f, ax = pyplot.subplots(figsize=(10,10))
for time in timeList:
  solutionDict = {}
  for dt in dtList:
    globstr = configuration + '_ne30_ne30_*/DT%04d_01_dycore_mac/*.cam.h0.0001-01-01-%05d.nc' % (dt, time)
    for fileName in glob.glob(globstr):
      words = fileName.split('_')
      print('Reading solution in ' + fileName)
      data = Dataset(fileName)
      T = data['T'][:]
      T = T[0,:,:]
      p0 = data['P0'][:]
      ps = data['PS'][:]
      ps = ps[0,:]
      area = data['area'][:]
      a = data['hyai'][:]
      b = data['hybi'][:]
      da = a[1:] - a[0:-1]
      db = b[1:] - b[0:-1]
      dp = np.outer(da,p0*np.ones(np.shape(ps))) + np.outer(db,ps)
      solutionDict[dt] = {'T': T, 'p0': p0, 'ps': ps, 'area': area, 'dp': dp}
   
  dtRef = min(dtList)
  dtDict = {}
  WRMSerror = {}
  qRef = solutionDict[dtRef]['T']
  psRef = solutionDict[dtRef]['ps']
  dpRef = solutionDict[dtRef]['dp']
  for dt in dtList:
    if (dt is dtRef):
      continue
    dtDict[float(dt)] = dt
    q = solutionDict[dt]['T']
    area = solutionDict[dt]['area']
    dp = solutionDict[dt]['dp']
    ps = solutionDict[dt]['ps']
    dpw = 0.5*(dpRef + dp)
    psw = 0.5*(psRef + ps)
    dz = dpw
    dxdy = np.outer(np.ones(np.shape(dpw)[0]), area)
    weight = dxdy*dz / np.sum(psw*area)
    WRMSerror[dt] = np.sqrt( np.sum((q-qRef)**2*weight) )
  
  dtPlot = np.sort(dtDict.keys())
  WRMSPlot = np.empty(len(dtPlot))
  measurePlot = np.array([measureList[-1], measureList[0]])
  for j,dt in enumerate(dtPlot):
    WRMSPlot[j] = WRMSerror[dtDict[dt]]

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
  ax.loglog(measurePlot, coeff*np.power(measurePlot,order), '-', color=colorDict[time], label=label, linewidth=3, markersize=12)
  ax.loglog(dtPlot, coeff*np.power(dtPlot,order), ':', color=colorDict[time], linewidth=3, markersize=12)
  print(time, order)

dtPlot = np.array([8, 1800])
ax.loglog(dtPlot, 4e-2/dtPlot[-1]*dtPlot, '--k', linewidth=3)
#ax.loglog(dtPlot, 2e-1/dtPlot[-1]*dtPlot, '--k', linewidth=3)
ax.set_ylabel('WRMS temperature error (K)', fontsize='x-large')
ax.set_xlabel('timestep size (s)', fontsize='x-large')
ax.legend(loc='upper left', fontsize='large')
ax.axis('equal')

f.tight_layout()
f.savefig('convergence_P111.pdf')  

  
pyplot.show()
