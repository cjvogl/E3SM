from netCDF4 import Dataset
import numpy as np
import glob
import sys

if (len(sys.argv) < 2):
  print "Usage: python determine_problematic_cell_SGR_P000.py <base_output_directory>"
  exit()
outpath = sys.argv[1]

dtList = [8,4,1]

numberList = [0, 1, 10, 20]

solutionDict = {}
for dt in dtList:
  globstr = '%s/RKZ_A1_B1_C2_ql17_lmt4_sgr000*_DT%04d_01_dycore_mac/*.cam.h0.0001-01-01-03600.nc' % (outpath, dt)
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
    solutionDict[dt] = {'T': T, 'lat': lat, 'lon': lon, 'area': area, 'dp': dp}

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

dpRef = refDict['dp']
levels = np.shape(T)[0]
total = levels*np.shape(T)[1]
overageSorted = np.sort(overage.flatten())
print("\nWorst overage occurs at index %d\n" % np.argmax(overage))

for number in numberList:
  overageTol = overageSorted[-(number+1)]

  area = dtDict['area']
  dp = dtDict['dp']
  dpw = 0.5*(dpRef + dp)
  dz = dpw
  dxdy = np.outer(np.ones(np.shape(dpw)[0]), area)
  weight = dxdy*dz / np.sum(dxdy*dz)
  weight = np.where(overage <= overageTol, weight, 0.0)
  WRMSerrorDt = np.sqrt( np.sum( errDt**2*weight) )
      
  area = dtX2dict['area']
  dp = dtX2dict['dp']
  dpw = 0.5*(dpRef + dp)
  dz = dpw
  dxdy = np.outer(np.ones(np.shape(dpw)[0]), area)
  weight = dxdy*dz / np.sum(dxdy*dz)
  weight = np.where(overage <= overageTol, weight, 0.0)
  WRMSerrorDtX2 = np.sqrt( np.sum( errDtX2**2*weight) )
  
  order = np.log(WRMSerrorDtX2/WRMSerrorDt)/np.log(2)
  print('remove %d of %d cells: %1.2f' % (number, total, order))

