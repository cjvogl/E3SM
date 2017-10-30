import glob
import matplotlib.pyplot as pyplot
from netCDF4 import Dataset
import numpy as np

noHV = True

methodDict = {'U35-ref': 5,
                  'U35': 11,
                  'ARS232-ref': 7,
                  'ARS232': 12,
                  'DBM453': 13,
                  'ARS222': 14}
                  #'ARS233': 15,
                  #'ARS343': 16,
                  #'ARS443': 17,
                  #'ARK324': 18}
                  #'ARK436': 19}

# Load reference solution
tsteptypeRef = 5
#dtRef = 0.00078125
#nmaxRef = 768000
dtRef = 0.0015625
nmaxRef = 384000
if (noHV):
    directory = './output_tsteptype%d_tstep%8.7f_nmax%d_nu0.0/dcmip2012_test31.nc' \
                    % (tsteptypeRef,dtRef,nmaxRef)
else:
    directory = './output_tsteptype%d_tstep%8.7f_nmax%d/dcmip2012_test31.nc' \
                    % (tsteptypeRef,dtRef,nmaxRef)
print 'Reading reference solution from ' + directory
data = Dataset(directory)
Tref = data['u'][:]
Tref = Tref[10,:,:,:]
shape = np.shape(Tref)
numElements = shape[0]*shape[1]*shape[2]

# Iterate through methods and plot convergence test
#f, (ax1,ax2) = pyplot.subplots(1,2,figsize=(12,6))
f, ax1 = pyplot.subplots()
f, ax2 = pyplot.subplots()
for method in methodDict.keys():
  dtList = []
  solutionDict = {}
  print method
  # Load timestep and solution data
  for fileName in glob.glob('tsteptype%d_tstep*.out' % methodDict[method]):
    if ((noHV and "nu0.0" in fileName) or (not noHV and "nu0.0" not in fileName)):
      words = fileName.split('_')
      dt = words[1].replace('tstep','')
      if (float(dt) > dtRef+1e-12):
        dtList.append(float(dt))
        directory = './output_'+fileName.replace('.out','')
        print 'Reading solution in ' + directory
        data = Dataset(directory+'/dcmip2012_test31.nc')
        T = data['u'][:]
        solutionDict[dt] = T[10,:,:,:]
  # Compute errors from reference solution
  if (dtRef in dtList):
      dtList.remove(dtRef)
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
  orderA = np.log(L2Plot[1]/L2Plot[0])/np.log(dtPlot[1]/dtPlot[0])
  orderB = np.log(L2Plot[2]/L2Plot[1])/np.log(dtPlot[2]/dtPlot[1])
  orderC = np.log(L2Plot[-1]/L2Plot[-2])/np.log(dtPlot[-1]/dtPlot[-2])
  ax1.loglog(dtPlot,L2Plot,'-o',label='%s (%3.2f, %3.2f, %3.2f)' % (method,orderA,orderB,orderC))
  orderA = np.log(LIPlot[1]/LIPlot[0])/np.log(dtPlot[1]/dtPlot[0])
  orderB = np.log(LIPlot[2]/LIPlot[1])/np.log(dtPlot[2]/dtPlot[1])
  orderC = np.log(LIPlot[-1]/LIPlot[-2])/np.log(dtPlot[-1]/dtPlot[-2])
  ax2.loglog(dtPlot,LIPlot,'-o',label='%s (%3.2f, %3.2f, %3.2f)' % (method,orderA,orderB,orderC))

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
ax1.axis('equal')
ax2.legend(loc='best')

pyplot.show()
