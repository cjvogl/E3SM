#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
from netCDF4 import Dataset
import numpy as np
import os.path

#dtList = [4,8]
dtList = [4]
config = 'cnvg_condensation/run/RKZ_SGR_extrapf_qv1_ql1_Al1_lmt4_adjIC_ne30_ne30_FC5AQUAP_intel_quartz'
timeList = [0, 1800, 3600, 5400, 7200]
variableList = ['T','Q','CLDLIQ','RKZ_RH',\
                'RKZ_AT','RKZ_Al','RKZ_Av','RKZ_fac', \
                'RKZ_dqsatdT', 'RKZ_zqme', 'RKZ_qme',\
                'RKZ_qlneg_lm4', 'RKZ_qvneg_lm4']
variableList = ['CLDLIQ','RKZ_fac','RKZ_qme','RKZ_Al','RKZ_term_A','RKZ_term_B','RKZ_term_C','nonphysical','RKZ_qlneg_lm4', 'RKZ_qvneg_lm4']
plotList = variableList

#ind = 962360
#ind = 963667
#ind = 746876
#ind = 1125452
#ind = 1069850
#ind = 855579
ind = 746876

solutionDict = {}
for var in variableList:
  if (var in plotList):
    f, ax = pyplot.subplots()
  for dt in dtList:
    directory = config + '/DT%04d_01_dycore_mac' % dt
    filePrefix = 'DT%04d_01_dycore_mac.cam.h0.0001-01-01-' % dt
    txtFile = '%s/%s_%d.txt' % (directory,var,ind)
    if (os.path.exists(txtFile)):
      tmp = np.loadtxt(txtFile)
      t = tmp[0,:]
      q = tmp[1,:]
    else:
      t = np.empty(0)
      q = np.empty(0)
      for time in timeList:
        fileName = '%s/%s%05d.nc' % (directory, filePrefix, time)
        print 'Reading time in from ' + fileName + '...',
        data = Dataset(fileName)
        t_current = data['time'][:]
        print 'done'
        t_current = t_current*3600*24
        numSlices = len(t_current)
        q_current = np.empty(numSlices)
        print ' Reading in ' + var + '...',
        if (var is 'nonphysical'):
          ql = data['CLDLIQ'][:]
          f = data['RKZ_fac'][:]
          for i in range(numSlices):
            q_current[i] = np.count_nonzero(np.where(f[i,:,:] < 1.0e-12, ql[i,:,:], 0))
        else:
          tmp = data[var][:]
          print ' done'
          for i in range(numSlices):
            flat = tmp[i,:,:].flatten()
            q_current[i] = flat[ind]
        t = np.hstack((t, t_current))
        q = np.hstack((q, q_current))
      tmp = np.vstack((t, q))
      np.savetxt(txtFile, tmp)
    solutionDict[var] = q
    if (var in plotList):
      ax.plot(t, q, '.', label='dt=%d' % dt)
  if (var in plotList):
    ax.set_title(var + '(%d)' % ind)
    ax.legend()
#    f.savefig(var + '_%d.png' % ind)

fig, ax = pyplot.subplots()
if ('RKZ_qme' in solutionDict.keys()):
  Q = solutionDict['RKZ_qme']
  ax.plot(t, Q, '-r', label='Q')
  if ('RKZ_fac' in solutionDict.keys()):
    f_n = solutionDict['RKZ_fac']
    f_np1 = 2.0*f_n[1:]- f_n[0:-1]
    f_np1 = np.where(f_np1 > 1.0, 1.0, f_np1)
    f_np1 = np.where(f_np1 < 0.0, 0.0, f_np1)
    ax2 = ax.twinx()
    ax2.plot(t[1:], f_n[1:], '-b')
    ax2.tick_params('y', colors='b')
    term_A = None
    term_B = None
    term_C = None
    if ('RKZ_term_A' in solutionDict.keys()):
      term_A = solutionDict['RKZ_term_A']
      ax.plot(t, term_A, '-k', label='in-cloudy condensation')
    if ('RKZ_term_B' in solutionDict.keys()):
      term_B = solutionDict['RKZ_term_B']
      ax.plot(t, term_B, '-m', label='out-cloud evaporation')
    if ('RKZ_term_C' in solutionDict.keys()):
      term_C = solutionDict['RKZ_term_C']
      term_CN = np.where(f_np1 > f_n[1:], term_C[1:], 0.0)
      ax.plot(t[1:], term_CN, '-g', label='cloud nucleation')
      term_CA = np.where(f_np1 < f_n[1:], term_C[1:], 0.0)
      ax.plot(t[1:], term_CA, '-c', label='cloud annihilation')
    if (not (term_A is None or term_B is None or term_C is None)):
      ax.plot(t, term_A + term_B + term_C, '--r')

ax.legend()
ax.set_ylabel('condensation rate', fontsize='xx-large')
ax2.set_ylabel('cloud fraction', fontsize='xx-large', color='b')


pyplot.show()
#fig.savefig('Q_%d.png' % ind)

