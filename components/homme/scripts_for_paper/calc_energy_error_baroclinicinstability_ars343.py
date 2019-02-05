from dictionaries import methodDict
import numpy as np
import sys

if (len(sys.argv) < 2):
  print "Usage: python calc_error_baroclinicinstability_ars232.py <base_output_directory>"
  exit()
outpath = sys.argv[1]

method = 'ARS343'
tsteps = [10, 20]

outputList = []
for tstep in tsteps:
  fileNames = [outpath + '/tsteptype%d_tstep%d_restart64800_ndays1_nu0.0_rsplit0.out' % (methodDict[method], tstep),
               outpath + '/tsteptype%d_tstep%d_restart64800_ndays1_nu1e15_rsplit0.out' % (methodDict[method], tstep),
               outpath + '/tsteptype%d_tstep%d_restart64800_ndays1_nu0.0_rsplit1.out' % (methodDict[method], tstep),
               outpath + '/tsteptype%d_tstep%d_restart64800_ndays1_nu1e15_rsplit1.out' % (methodDict[method], tstep)]
  for m,fileName in enumerate(fileNames):
    print('Reading in from ' + fileName)
    fileObj = open(fileName)
    lines = list(fileObj)
    fileObj.close()
    time = 0.0
    timeList = [time]
    errorList = [0.0,]
    for line in lines:
      if ('time=' in line):
        time = float(line.split()[3])
        if ('[s]' in line):
          time /= 3600.0
      if ('(E-E0)/E0' in line):
        timeList.append(time)
        errorList.append(float(line.split()[1]))
    output = np.vstack((timeList, errorList))
    words = fileName.split('/')
    filename = words[-1].replace('.out','')
    filename = './data/' + filename + '_energy_error_data.txt'
    print('Writing out ' + filename)
    np.savetxt(filename, output)
