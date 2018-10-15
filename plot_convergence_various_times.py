import glob
import matplotlib.pyplot as pyplot
import matplotlib
from netCDF4 import Dataset
import numpy as np

matplotlib.rcParams.update({'font.size': 18})

#dtList = [1800,450,120,30,15,8,4,2,1]
dtList = [1800,450,120,75,30,15,8,4,1]
configDict = {}

timeList = [1800, 3600, 5400, 7200, 9000, 10800, 12600]
#timeList = [1800, 3600]

f2, ax2 = pyplot.subplots(figsize=(10,10))
f2.set_tight_layout(True)

for time in timeList:
  for dt in dtList:
    globstr = 'RKZ_P3_adjIC_lmt4_ne30_ne30_*/DT%04d_01_dycore_mac/*.cam.h0.0001-01-01-%05d.nc' % (dt, time)
    for fileName in glob.glob(globstr):
      words = fileName.split('_')
      config = words[0]
      for j in range(1,len(words)):
        if (words[j] == 'ne30'):
          break
        config = config + '_' + words[j]
      if (config not in configDict.keys()):
        solutionDict = {}
        configDict[config] = solutionDict
        print('\nadding ' + config + ' to configuration list\n')
      solutionDict = configDict[config]
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
  for config in configDict.keys(): 
    dtDict = {}
    WRMSerror = {}
    solutionDict = configDict[config]
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

    print(WRMSPlot) 
 
    if (len(dtList) > 1):
      orderWRMS = np.log(WRMSPlot[1:]/WRMSPlot[0:-1])/np.log(dtPlot[1:]/dtPlot[0:-1])
    else:
      orderWRMS = (0.0,)

    label = "%d second run (%3.2f)" % (time, orderWRMS[0])
    ax2.loglog(dtPlot, WRMSPlot, '-o', label=label, linewidth=3, markersize=12)

ax2.set_ylabel('WRMS Error', fontsize='large')
ax2.set_xlabel('dt', fontsize='large')
ax2.axis('equal')

ax2.legend(loc='best')

pyplot.tight_layout()

pyplot.show()






