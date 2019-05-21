import glob
import matplotlib.pyplot as pyplot
import matplotlib
from netCDF4 import Dataset
import numpy as np

matplotlib.rcParams.update({'font.size': 24})

dtList = [1800,450,120,75,30,15,8,4,1]

timeList = [3600, 7200, 10800]
linespecDict = {3600: '-bo',
                7200: '-gs',
                10800: '-r^'}

configuration = 'RKZ_SGR_extrapf_qv0_ql0_Al0_lmt4_adjIC'

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
  for j,dt in enumerate(dtPlot):
    WRMSPlot[j] = WRMSerror[dtDict[dt]]

  if (len(dtList) > 1):
    orderWRMS = np.log(WRMSPlot[1:]/WRMSPlot[0:-1])/np.log(dtPlot[1:]/dtPlot[0:-1])
  else:
    orderWRMS = (0.0,)

  label = "%0.0f hour run" % (time/3600.0)
  ax.loglog(dtPlot, WRMSPlot, linespecDict[time], label=label, linewidth=3, markersize=12)
  print(np.mean(orderWRMS[0:4]), np.sqrt(np.var(orderWRMS[0:4])))

dtPlot = np.array([4, 1800])
ax.loglog(dtPlot, 4e-2/dtPlot[-1]*dtPlot, '--k', linewidth=3)
#ax.loglog(dtPlot, 2e-1/dtPlot[-1]*dtPlot, '--k', linewidth=3)
ax.set_ylabel('WRMS temperature error (K)', fontsize='x-large')
ax.set_xlabel('timestep size (s)', fontsize='x-large')
ax.legend(loc='best', fontsize='large')
ax.axis('equal')

f.tight_layout()
f.savefig('convergence_P000.pdf')  

  
pyplot.show()
