import glob
import matplotlib.pyplot as pyplot
from netCDF4 import Dataset
import numpy as np

splitting = '2'

methodDict = {#'U35-ref': 5,
                  #'U35': 11,
                  'ARS232-native': 7,
                  'KGS242-native': 8,
                  'KGS252-native': 9,
#                  'KGS262-native': 10,
#                  'KGS272-native': 11,
                  'ARS232': 22,
                  'DBM453': 23,
                  'ARS222': 24,
                  'ARS233': 25,
                  'ARS343': 26,
                  'ARS443': 27,
                  'ARK324': 28,
                  'ARK436': 29,
                  'SSP3333b': 30,
                  'SSP3333c': 31}


# Iterate through methods and find largest step size
maxsizeDict = {}
for m,method in enumerate(methodDict.keys()):
  maxsizeDict[method] = -np.inf
  print method
  # Load timestep and solution data
#  if (methodDict[method] < 12):
  globstr = 'tsteptype%d_tstep*.out' % methodDict[method]
#  else:
#    globstr = 'tsteptype%d_tstep*_splitting%s*.out' % \
#                                         (methodDict[method], splitting)
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

print ""
for key in maxsizeDict.keys():
  print key, methodDict[key], maxsizeDict[key]
