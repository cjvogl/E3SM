from netCDF4 import Dataset
import numpy as np
import sys
import glob

if (len(sys.argv) < 2):
  print "Usage: python calc_problematic_cell_SGR_P000.py <base_output_directory>"
  exit()
outpath = sys.argv[1]

configuration = 'RKZ_SGR_extrapf_qv0_ql0_Al0_lmt4_adjIC'
timeList = [0, 1800, 3600, 5400, 7200]
variableList = ['CLDLIQ','RKZ_fac','RKZ_qme','RKZ_term_A','RKZ_term_B','RKZ_term_C']

dt = 4
ind = 1069852

for var in variableList:
  outFile = './data/SGR_P000_' + var + '_%d.txt' % ind
  t = np.empty(0)
  q = np.empty(0)
  for time in timeList:
    globstr = outpath + '/' + configuration + \
      '_ne30_ne30_*/DT%04d_01_dycore_mac/*.cam.h0.0001-01-01-%05d.nc' % (dt, time)
    fileList = glob.glob(globstr)
    if (len(fileList) > 1):
      print "more than one solution file found"
      exit()
    fileName = fileList[0]
    print 'Reading time in from ' + fileName + '...',
    data = Dataset(fileName)
    t_current = data['time'][:]
    print 'done'
    t_current = t_current*3600*24
    numSlices = len(t_current)
    q_current = np.empty(numSlices)
    print ' Reading in ' + var + '...',
    tmp = data[var][:]
    print ' done'
    for i in range(numSlices):
      flat = tmp[i,:,:].flatten()
      q_current[i] = flat[ind]
    t = np.hstack((t, t_current))
    q = np.hstack((q, q_current))
  tmp = np.vstack((t, q))
  print('Writing ' + outFile)
  np.savetxt(outFile, tmp)
