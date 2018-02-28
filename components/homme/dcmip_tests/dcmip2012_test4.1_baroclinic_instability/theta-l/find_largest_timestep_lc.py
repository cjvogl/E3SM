import glob
import matplotlib.pyplot as pyplot
from netCDF4 import Dataset
import numpy as np

splitting = '2'

methodDict = {#'U35-ref': 5,
                  #'U35': 11,
                  'ARS232-ref': 7,
                  'KGS242': 8,
                  'KGS252': 9,
                  'KGS262': 10,
                  'KGS272': 11}
#                  'ARS232': 12,
#                  'DBM453': 13,
#                  'ARS222': 14,
#                  'ARS233': 15,
#                  'ARS343': 16,
#                  'ARS443': 17,
#                  'ARK324': 18,
#                  'ARK436': 19,
#                  'SSP3333b': 20,
#                  'SSP3333c': 21}


# Iterate through methods and find largest step size
maxsizeDict = {}
for m,method in enumerate(methodDict.keys()):
  maxsizeDict[method] = -np.inf
  print method
  # Load timestep and solution data
  if (methodDict[method] < 12):
    globstr = 'tsteptype%d_tstep*.out' % methodDict[method]
  else:
    globstr = 'tsteptype%d_tstep*_splitting%s*.out' % \
                                         (methodDict[method], splitting)
  for fileName in glob.glob(globstr):
    print "reading " + fileName
    words = fileName.split('_')
    dt = words[1].replace('tstep','').replace('.out','')
    f = open(fileName)
    lines = list(f)
    f.close()
    for line in reversed(lines):
      if ("Finished main timestepping loop" in line):
        if (float(dt) > maxsizeDict[method]):
          maxsizeDict[method] = float(dt)
        break

print maxsizeDict
