from netCDF4 import Dataset
import matplotlib.pyplot as pyplot
import matplotlib
import numpy as np
import glob
from mpl_toolkits.basemap import Basemap

matplotlib.rcParams.update({'font.size': 18})


config = 'RKZ_SGR_extrapf_qv0_ql0_Al0_lmt4_adjIC*'
dtList = [8,4,1]
#timeList = [1800, 3600, 5400, 7200]
timeList = [7200]


numberList = [0, 1, 10, 20]

#overageDict = {}
#convergenceDict = {}
for time in timeList:
  solutionDict = {}
  f, ax = pyplot.subplots(4,3, figsize=(15,10))
  #f.set_visible(False)
  for dt in dtList:

    globstr = '%s/DT%04d_01_dycore_mac/*.cam.h0.0001-01-01-%05d.nc' % (config, dt, time)
    #globstr = 'RKZ_A1_B1_C2_ql19_lmt4_adjIC_*/DT%04d_01_dycore_mac/*.cam.h0.0001-01-01-%05d.nc' % (dt, time)
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
  #convergenceRate = np.log(errDtX2/errDt)/np.log(2)
  #convergenceDict[time] = convergenceRate
  overage = errDt - 0.5*errDtX2
  overage = np.where(overage > 0.0, overage, 0.0)
  #overageDict[time] = overage

  bmap = Basemap(lon_0=180, ax=ax[0,0])
  x,y = bmap(lon,lat)
  #tmp = bmap.scatter(x, y, c=np.amin(convergenceRate, axis=0), vmin=0, vmax=1.0)
  tmp = bmap.scatter(x, y, c=np.amin(overage, axis=0), cmap='Reds', vmin=0.0)
  f.colorbar(tmp, ax=ax[0,0])
  ax[0,0].set_title('%s (min)' % time)
  ax[0,0].set_ylabel('rate')

  bmap = Basemap(lon_0=180, ax=ax[0,1])
  x,y = bmap(lon,lat)
  #tmp = bmap.scatter(x, y, c=np.mean(convergenceRate, axis=0), vmin=0, vmax=1.0)
  tmp = bmap.scatter(x, y, c=np.mean(overage, axis=0), cmap='Reds', vmin=0.0)
  f.colorbar(tmp, ax=ax[0,1])
  ax[0,1].set_title('%s (mean)' % time)
  ax[0,1].set_ylabel('rate')

  bmap = Basemap(lon_0=180, ax=ax[0,2])
  x,y = bmap(lon,lat)
  #tmp = bmap.scatter(x, y, c=np.amax(convergenceRate, axis=0), vmin=0, vmax=1.0)
  tmp = bmap.scatter(x, y, c=np.amax(overage, axis=0), cmap='Reds', vmin=0.0)
  f.colorbar(tmp, ax=ax[0,2])
  ax[0,2].set_title('%s (max)' % time)
  ax[0,2].set_ylabel('rate')

  bmap = Basemap(lon_0=180, ax=ax[1,0])
  x,y = bmap(lon,lat)
  #indices = np.argmin(convergenceRate, axis=0)
  indices = np.argmin(overage, axis=0)
  val = np.empty(np.shape(x))
  for i,ind in enumerate(indices):
    val[i] = errDt[ind,i]
  tmp = bmap.scatter(x, y, c=val, cmap='Reds', vmin=0.0)
  f.colorbar(tmp, ax=ax[1,0])
  ax[1,0].set_ylabel('error')

  bmap = Basemap(lon_0=180, ax=ax[1,1])
  x,y = bmap(lon,lat)
  tmp = bmap.scatter(x, y, c=np.mean(errDt, axis=0), cmap='Reds', vmin=0.0)
  f.colorbar(tmp, ax=ax[1,1])
  ax[1,1].set_ylabel('error')

  bmap = Basemap(lon_0=180, ax=ax[1,2])
  x,y = bmap(lon,lat)
  #indices = np.argmax(convergenceRate, axis=0)
  indices = np.argmax(overage, axis=0)
  for i,ind in enumerate(indices):
    val[i] = errDt[ind,i]
  tmp = bmap.scatter(x, y, c=val, cmap='Reds', vmin=0.0)
  f.colorbar(tmp, ax=ax[1,2])
  ax[1,2].set_ylabel('error')
  
  levels = np.shape(T)[0]

  #ax[2,0].plot(np.amin(convergenceRate, axis=1), np.arange(levels), '-o')
  ax[2,0].plot(np.amin(overage, axis=1), np.arange(levels), '-o')
  ax[2,0].set_xlabel('convergence order')
  ax[2,0].set_ylabel('level #')

  #ax[2,1].plot(np.mean(convergenceRate, axis=1), np.arange(levels), '-o')
  ax[2,1].plot(np.mean(overage, axis=1), np.arange(levels), '-o')
  ax[2,1].set_xlabel('convergence order')
  ax[2,1].set_ylabel('level #')

  #ax[2,2].plot(np.amax(convergenceRate, axis=1), np.arange(levels), '-o')
  ax[2,2].plot(np.amax(overage, axis=1), np.arange(levels), '-o')
  ax[2,2].set_xlabel('convergence order')
  ax[2,2].set_ylabel('level #')

  #indices = np.argmin(convergenceRate, axis=1)
  indices = np.argmin(overage, axis=1)
  val = np.zeros(levels)
  for i,ind in enumerate(indices):
    val[i] = errDt[i,ind]
  ax[3,0].plot(val, np.arange(levels), '-o')
  ax[3,0].set_xlabel('error')
  ax[3,0].set_ylabel('level #')

  ax[3,1].plot(np.mean(errDt, axis=1), np.arange(levels), '-o')
  ax[3,1].set_xlabel('error')
  ax[3,1].set_ylabel('level #')

  #indices = np.argmax(convergenceRate, axis=1)
  indices = np.argmax(overage, axis=1)
  val = np.zeros(levels)
  for i,ind in enumerate(indices):
    val[i] = errDt[i,ind]
  ax[3,2].plot(val, np.arange(levels), '-o')
  ax[3,2].set_xlabel('error')
  ax[3,2].set_ylabel('level #')
 
  if (time == timeList[-1]):
    psRef = refDict['ps']
    dpRef = refDict['dp']
    total = levels*np.shape(T)[1]
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



#convergenceLast = convergenceDict[timeList[-1]]
#convergenceFirst = convergenceDict[timeList[0]]
#
#f, ax = pyplot.subplots(2)
#bmap = Basemap(lon_0=180, ax=ax[0])
#x,y = bmap(lon,lat)
#tmp = bmap.scatter(x, y, c=np.amin(convergenceLast-convergenceFirst, axis=0))
#f.colorbar(tmp, ax=ax[0])
#ax[0].set_title('degradation from %s to %s' % (timeList[0],timeList[-1]))
#ax[0].set_xlabel('lon')
#ax[0].set_ylabel('lat')
#
#ax[1].plot(np.amin(convergenceLast-convergenceFirst, axis=1), np.arange(levels), '-o')
#ax[1].set_xlabel('convergence order')
#ax[1].set_ylabel('level #')
#
pyplot.show()

