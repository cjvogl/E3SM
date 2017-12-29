import glob
import matplotlib.pyplot as pyplot
from netCDF4 import Dataset
import numpy as np

noHV = False
rtol = '1e-3'
splitting = '1'

methodDict = {#'U35-ref': 5,
                  #'U35': 11,
                  'ARS232-ref': 7,
                  'ARS232': 12,
                  'DBM453': 13,
                  'ARS222': 14,
                  'ARS233': 15,
                  'ARS343': 16,
                  'ARS443': 17,
#                  'ARK324': 18}
                  #'ARK436': 19}
                  'SSP3333b': 20,
                  'SSP3333c': 21}

# Load reference solution
tsteptypeRef = 5
dtRef = 0.000390625
nmaxRef = 768000
if (noHV):
    directory = './output_tsteptype%d_tstep%10.9f_nmax%d_nu0.0/dcmip2012_test31.nc' \
                    % (tsteptypeRef,dtRef,nmaxRef)
else:
    directory = './output_tsteptype%d_tstep%10.9f_nmax%d/dcmip2012_test31.nc' \
                    % (tsteptypeRef,dtRef,nmaxRef)
print 'Reading reference solution from ' + directory
data = Dataset(directory)
Tref = data['T'][:]
Tref = Tref[10,:,:,:]
shape = np.shape(Tref)
numElements = shape[0]*shape[1]*shape[2]

# Iterate through methods and plot convergence test
f, ax1 = pyplot.subplots(figsize=(10,10))
f, ax2 = pyplot.subplots(figsize=(10,10))
for m,method in enumerate(methodDict.keys()):
  dtList = []
  solutionDict = {}
  print method
  # Load timestep and solution data
  if (methodDict[method] < 12):
    globstr = 'tsteptype%d_tstep*.out' % methodDict[method]
  else:
    globstr = 'tsteptype%d_tstep*_rtol%s*_splitting%s*.out' % \
                                         (methodDict[method], rtol, splitting)
  for fileName in glob.glob(globstr):
    if ((noHV and "nu0.0" in fileName) or (not noHV and "nu0.0" not in fileName)):
      words = fileName.split('_')
      dt = words[1].replace('tstep','')
      if (float(dt) > 2.0*dtRef+1e-12 and float(dt) < 10.0):
        directory = './output_'+fileName.replace('.out','')
        print 'Reading solution in ' + directory
        data = Dataset(directory+'/dcmip2012_test31.nc')
        T = data['T'][:]
        if (np.shape(T)[0] == 11):
          dtList.append(float(dt))
          solutionDict[dt] = T[10,:,:,:]
        else:
          print '... skipping due to incomplete results ...'

  # Compute errors from reference solution
  L2error = {}
  LIerror = {}
  for dt in dtList:
    T = solutionDict[str(dt)]
    L2error[str(dt)] = np.sqrt(np.sum((T-Tref)**2)/numElements)
    LIerror[str(dt)] = np.amax(abs(T-Tref))
  # Sort values and add to plot
  dtPlot = np.sort(dtList)
  L2Plot = np.empty(len(dtPlot))
  LIPlot = np.empty(len(dtPlot))
  for j,dt in enumerate(dtPlot):
    L2Plot[j] = L2error[str(dt)]
    LIPlot[j] = LIerror[str(dt)]
  if (m < 7):
    lineStyle = '-o'
  else:
    lineStyle = '-d'
  order = np.log(L2Plot[1:]/L2Plot[0:-1])/np.log(dtPlot[1:]/dtPlot[0:-1])
  ax1.loglog(dtPlot,L2Plot,lineStyle,label='%s (final=%3.2f, max=%3.2f)' % (method,order[0],np.amax(order)))
  order = np.log(LIPlot[1:]/LIPlot[0:-1])/np.log(dtPlot[1:]/dtPlot[0:-1])
  ax2.loglog(dtPlot,LIPlot,lineStyle,label='%s (final=%3.2f, max=%3.2f)' % (method,order[0],np.amax(order)))

#orderList = ('first', 'second', 'third')
#for j, order in enumerate(orderList):
#    coeff = L2Plot[0]/dtPlot[0]**(j+1)
#    ax1.loglog(dtPlot,coeff*dtPlot**(j+1),'--',label=order)
#    coeff = LIPlot[0]/dtPlot[0]**(j+1)
#    ax2.loglog(dtPlot,coeff*dtPlot**(j+1),'--',label=order)

ax1.set_title('L2 Error')
ax1.set_xlabel('dt')
ax1.axis('equal')
ax1.legend(loc='best')
ax2.set_title('LI Error')
ax2.set_xlabel('dt')
ax2.axis('equal')
ax2.legend(loc='best')

pyplot.show()
