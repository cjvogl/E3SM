import glob
from netCDF4 import Dataset
import numpy as np
import sys

if (len(sys.argv) < 2):
  print "Usage: python calc_convergence_SGR_P111.py <base_output_directory>"
  exit()
outpath = sys.argv[1]

dtList = [1800,450,120,75,30,15,8,4,1]
timeList = [10800, 21600, 32400, 43200]

configuration = 'RKZ_SGR_extrapf_qv1_ql1_Al1_lmt4_adjIC'

for time in timeList:
  solutionDict = {}

  for dt in dtList:
    globstr = outpath + '/' + configuration + '_ne30_ne30_*/DT%04d_01_dycore_mac/*.cam.h0.0001-01-01-%05d.nc' % (dt, time)
    fileList = glob.glob(globstr)
    if (len(fileList) > 1):
      print "more than one solution file found..."
      exit()
    fileName = fileList[0]
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
  output = np.vstack((dtPlot, WRMSPlot))
  fileName = './data/SGR_P111_convergence_data_%d.txt' % time
  print('Writing out ' + fileName)
  np.savetxt(fileName, output) 
