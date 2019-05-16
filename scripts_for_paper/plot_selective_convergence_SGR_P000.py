from netCDF4 import Dataset
import matplotlib.pyplot as pyplot
import matplotlib
import numpy as np
import glob
from mpl_toolkits.basemap import Basemap

matplotlib.rcParams.update({'font.size': 18})


config = 'RKZ_SGR_extrapf_qv0_ql0_Al0_lmt4_adjIC*'
dtList = [8,4,1]
time = 5400

numberList = [0, 5, 10, 20]

solutionDict = {}
for dt in dtList:

  globstr = '%s/DT%04d_01_dycore_mac/*.cam.h0.0001-01-01-%05d.nc' % (config, dt, time)
  for fileName in glob.glob(globstr):
    print(fileName)
    data = Dataset(fileName)
    lat = data['lat']
    lon = data['lon']
    T = data['T'][0,:,:]
    p0 = data['P0'][:]
    ps = data['PS'][0,:]
    area = data['area'][:]
    a = data['hyai'][:]
    b = data['hybi'][:]
    da = a[1:] - a[0:-1]
    db = b[1:] - b[0:-1]
    dp = np.outer(da,p0*np.ones(np.shape(ps))) + np.outer(db,ps)
    solutionDict[dt] = {'T': T, 'lat': lat, 'lon': lon, 'p0': p0, 'ps': ps, 'area': area, 'dp': dp}

refDict = solutionDict[dtList[2]]
dtDict = solutionDict[dtList[1]]
dtX2dict = solutionDict[dtList[0]]
Tref = refDict['T']
lat = refDict['lat']
lon = refDict['lon']
Tdt = dtDict['T']
TdtX2 = dtX2dict['T']
errDtX2 = abs(TdtX2 - Tref)
errDt = abs(Tdt - Tref)
overage = errDt - 0.5*errDtX2
overage = np.where(overage > 0.0, overage, 0.0)

psRef = refDict['ps']
dpRef = refDict['dp']
levels = np.shape(Tref)[0]
total = levels*np.shape(Tref)[1]
overageSorted = np.sort(overage.flatten())
print("\nWorst overage occurs at index %d\n" % np.argmax(overage))
f, ax = pyplot.subplots() 

for number in numberList:
  overageTol = overageSorted[-(number+1)]

  area = dtDict['area']
  ps = dtDict['ps']
  dp = dtDict['dp']
  psw = 0.5*(psRef + ps)
  dpw = 0.5*(dpRef + dp)
  dz = dpw
  dxdy = np.outer(np.ones(np.shape(dpw)[0]), area)
  weight = dxdy*dz / np.sum(psw*area)
  weight = np.where(overage <= overageTol, weight, 0.0)
  WRMSerrorDt = np.sqrt( np.sum( errDt**2*weight) )
      
  area = dtX2dict['area']
  ps = dtX2dict['ps']
  dp = dtX2dict['dp']
  psw = 0.5*(psRef + ps)
  dpw = 0.5*(dpRef + dp)
  dz = dpw
  dxdy = np.outer(np.ones(np.shape(dpw)[0]), area)
  weight = dxdy*dz / np.sum(psw*area)
  weight = np.where(overage <= overageTol, weight, 0.0)
  WRMSerrorDtX2 = np.sqrt( np.sum( errDtX2**2*weight) )

  order = np.log(WRMSerrorDtX2/WRMSerrorDt)/np.log(2)
  ax.loglog((dtList[1],dtList[0]), (WRMSerrorDt, WRMSerrorDtX2), '-o', \
            label='remove %d of %d cells: %1.2f' % (number, total, order))
ax.legend()
pyplot.show()

