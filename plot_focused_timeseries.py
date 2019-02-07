#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
from netCDF4 import Dataset
import numpy as np
import os.path

#dtList = [4,8]
dtList = [4]
config = 'cnvg_condensation/run/RKZ_SGR_qv1_ql1_Al1_lmt4_adjIC_ne30_ne30_FC5AQUAP_intel_quartz'
timeList = [0, 1800, 3600, 5400, 7200]
variableList = ['T','Q','CLDLIQ','RKZ_RH',\
                'RKZ_AT','RKZ_Al','RKZ_Av','RKZ_fac', \
                'RKZ_dqsatdT', 'RKZ_zqme', 'RKZ_qme',\
                'RKZ_qlneg_lm4', 'RKZ_qvneg_lm4']
variableList = ['CLDLIQ','RKZ_fac','RKZ_qme','RKZ_Al','Q', 'T', 'RKZ_Av','RKZ_RH', 'RKZ_qlneg_lm4', 'RKZ_qvneg_lm4']
variableList = ['CLDLIQ','RKZ_fac','RKZ_qme']
plotList = variableList

ind = 962360
#ind = 963667
#ind = 746876

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

#Av = solutionDict['RKZ_Av']
#dqsatdT = solutionDict['RKZ_dqsatdT']
#AT = solutionDict['RKZ_AT']
#gam = solutionDict['gam']
#gamma = 1.0 + gam
#Al = solutionDict['RKZ_Al']
#RH = solutionDict['RKZ_RH']
#dfdRH = solutionDict['RKZ_dfdRH']
#qv = solutionDict['Q']
#ql = solutionDict['CLDLIQ']
#qsat = qv/RH
#zforcing = (Av - RH*dqsatdT*AT)/qsat
#zc3 = (1.0 + gam*RH*dqsatdT)/qsat
#D = Av - dqsatdT*AT + gamma*Al

fig, ax = pyplot.subplots()
Q = solutionDict['RKZ_qme']
ax.plot(t, Q, '-r', label='Q')

#  f = solutionDict['RKZ_f']
#  termA = f*(Av - dqsatdT*AT)/gamma
#  termB = -(1.0 - f)*Al
#  Q = (f*D/gamma - Al + dt*D/gamma*dfdRH*zforcing)/(1.0 + dt*D/gamma*dfdRH*zc3)
#  termC = dt*D/gamma*dfdRH*( (Av - RH*dqsatdT*AT)/qsat - (1.0 + gam*RH*dqsatdT)/qsat*Q)
#  print(np.amax(Q - (termA+termB+termC)))
#  
#  ax.plot(t, termA, '-r', label='term A')
#  ax.plot(t, termB, '-g', label='term B')
#  ax.plot(t, termC, '-b', label='term C')
#  ax.plot(t, Q, '--k', label='Q')
ax.legend()
  
pyplot.show()
