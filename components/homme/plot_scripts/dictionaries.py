from netCDF4 import Dataset
import glob

# This is meant to be an all-inclusive list of possible methods
methodDict = {'KGU35-native': 5,
              'ARS232-native':  7,
              'KGU35': 21,
              'ARS232': 22,
              'DBM453': 23,
              'ARS222': 24,
              'ARS233': 25,
              'ARS343': 26,
              'ARS443': 27,
              'ARK324': 28,
              'ARK436': 29,
              'SSP3333b': 30,
              'SSP3333c': 31,
              'IMEXKG232a': 32,
              'IMEXKG232b': 33,
              'IMEXKG242a': 34,
              'IMEXKG242b': 35,
              'IMEXKG243a': 36,
              'IMEXKG243b': 37,
              'IMEXKG252a': 38,
              'IMEXKG252b': 39,
              'IMEXKG253a': 40,
              'IMEXKG253b': 41,
              'IMEXKG254a': 42,
              'IMEXKG254b': 43,
              'IMEXKG254c': 44,
              'IMEXKG343a': 45,
              'IMEXKG343b': 46}

# Line+marker styles chosen to provide certain information about the method
#   line:
#     solid      = second order
#     dashed     = third order
#     dot-dashed = fourth order
#   marker:
#     circle   = 0 implicit stages
#     x        = 2 implicit stages
#     triangle = 3 implicit stages
#     square   = 4 implicit stages
#     pentagon = 5 implicit stages
lineStyleDict = {'KGU35-native': '--o',
                 'ARS232-native': '-x',
                 'KGU35': '--o',
                 'ARS222': '-x',
                 'ARS232': '-x',
                 'SSP3333b': '--x',
                 'SSP3333c': '--x',
                 'ARS443': '--s',
                 'DBM453': '--s',
                 'ARS233': '--x',
                 'ARS343': '--^',
                 'ARK324': '--^',
                 'ARK436': '-.p',
                 'IMEXKG232a': '-x',
                 'IMEXKG232b': '-x',
                 'IMEXKG242a': '-x',
                 'IMEXKG242b': '-x',
                 'IMEXKG243a': '-^',
                 'IMEXKG243b': '-^',
                 'IMEXKG252a': '-x',
                 'IMEXKG252b': '-x',
                 'IMEXKG253a': '-^',
                 'IMEXKG253b': '-^',
                 'IMEXKG254a': '-s',
                 'IMEXKG254b': '-s',
                 'IMEXKG254c': '-s',
                 'IMEXKG343a': '--^',
                 'IMEXKG343b': '--^'}

colorDict = {'KGU35-native': 'k',
             'KGU35': 'tab:blue',
             'ARS232-native': 'k',
             'ARS222': 'tab:blue',
             'ARS232': 'tab:orange',
             'SSP3333b': 'tab:green',
             'SSP3333c': 'k',
             'ARS233': 'tab:red',
             'IMEXKG232a': 'tab:purple',
             'IMEXKG232b': 'tab:brown',
             'IMEXKG242a': 'tab:pink',
             'IMEXKG242b': 'tab:gray',
             'IMEXKG252a': 'tab:olive',
             'IMEXKG252b': 'tab:cyan',
             'ARS343': 'tab:blue',
             'ARK324': 'tab:orange',
             'IMEXKG243a': 'tab:green',
             'IMEXKG243b': 'tab:red',
             'IMEXKG253a': 'tab:purple',
             'IMEXKG253b': 'tab:brown',
             'IMEXKG343a': 'tab:pink',
             'IMEXKG343b': 'tab:gray',
             'ARS443': 'tab:blue',
             'DBM453': 'tab:orange',
             'IMEXKG254a': 'tab:green',
             'IMEXKG254b': 'tab:red',
             'IMEXKG254c': 'tab:purple',
             'ARK436': 'tab:blue'}

# method to create solution dictionary
def create_solution_dict(globstr, testName, varRef, indRef, tRef, minDt, maxDt, suffix, suffix_omit):
  solutionDict = {}
  for fileName in glob.glob(globstr):
    if (suffix_omit is not None and suffix_omit in fileName):
      continue
    words = fileName.split('/')
    shortName = words[-1]
    words = shortName.split('_')
    dt = words[1].replace('tstep','')
    dt = dt.replace('.out','')
    # If current timestep is appropriate, obtain solution if it exists
    if (float(dt) < maxDt and float(dt) > minDt):
      directory = './output_'+fileName.replace('.out','')
      print 'Reading solution in ' + directory
      data = Dataset(directory+'/'+testName)
      q = data[varRef][:]
      t = data['time'][:]
      if (len(t) > indRef and abs(t[indRef] - tRef) < 1e-10):
        solutionDict[dt] = q[indRef,:,:,:]
      else:
        print '... skipping due to incomplete results ...'
  return solutionDict

# method to create energy error dictionary
def create_energy_error_dict(globstr, minDt, maxDt, suffix, suffix_omit):
  energyErrorDict = {}
  for fileName in glob.glob(globstr):
    if (suffix_omit is not None and suffix_omit in fileName):
      continue
    words = fileName.split('/')
    shortName = words[-1]
    words = shortName.split('_')
    dt = words[1].replace('tstep','')
    dt = dt.replace('.out','')
    if (float(dt) < maxDt and float(dt) > minDt):
      energyErrorMax = 0.0
      print 'Reading diagnostics in ' + fileName
      fileObj = open(fileName)
      lines = list(fileObj)
      fileObj.close()
      flag = False
      for line in reversed(lines):
        if ('Finished main timestepping loop' in line):
          flag = True
        if (flag and '(E-E0)/E0' in line):
          words = line.split()
          currentError = float(words[1])
          if (abs(energyErrorMax) < abs(currentError)):
            energyErrorMax = currentError
          break
      if (flag):
        energyErrorDict[dt] = energyErrorMax
      else:
        print '... skipping due to incomplete results ...'
  return energyErrorDict

def create_walltime_dict(globstr, minDt, maxDt, suffix, suffix_omit):
  walltimeDict = {}
  for fileName in glob.glob(globstr):
    if (suffix_omit is not None and suffix_omit in fileName):
      continue
    words = fileName.split('/')
    shortName = words[-1]
    words = shortName.split('_')
    dt = words[1].replace('tstep','')
    dt = dt.replace('.out','')
    if (float(dt) < maxDt and float(dt) > minDt):
      fileNameTiming = fileName.replace('.out','.err')
      print 'Reading timing info in ' + fileNameTiming
      fileObj = open(fileName)
      lines = list(fileObj)
      fileObj.close()
      flag = False
      for line in reversed(lines):
        if ('Finished main timestepping loop' in line):
          flag = True
          break
      if (not flag):
        print '... skipping due to incomplete results ...'
        continue
      fileObj = open(fileNameTiming)
      lines = list(fileObj)
      fileObj.close()
      words = lines[1].split()
      word = words[1] # ignore the word 'real'
      words = word.split(('m'))
      seconds = int(words[0])*60 + float(words[1].replace('s',''))
      walltimeDict[dt] = seconds
  return walltimeDict

# method to create solution dictionary
def create_solution_dict(globstr, testName, varRef, indRef, tRef, minDt, maxDt, suffix, suffix_omit):
  solutionDict = {}
  for fileName in glob.glob(globstr):
    if (suffix_omit is not None and suffix_omit in fileName):
      continue
    words = fileName.split('/')
    shortName = words[-1]
    words = shortName.split('_')
    dt = words[1].replace('tstep','')
    dt = dt.replace('.out','')
    # If current timestep is appropriate, obtain solution if it exists
    if (float(dt) < maxDt and float(dt) > minDt):
      directory = './output_'+fileName.replace('.out','')
      print 'Reading solution in ' + directory
      data = Dataset(directory+'/'+testName)
      q = data[varRef][:]
      t = data['time'][:]
      if (len(t) > indRef and abs(t[indRef] - tRef) < 1e-10):
        solutionDict[dt] = q[indRef,:,:,:]
      else:
        print '... skipping due to incomplete results ...'
  return solutionDict

# method to create energy error dictionary
def create_energy_error_dict(globstr, minDt, maxDt, suffix, suffix_omit):
  energyErrorDict = {}
  for fileName in glob.glob(globstr):
    if (suffix_omit is not None and suffix_omit in fileName):
      continue
    words = fileName.split('/')
    shortName = words[-1]
    words = shortName.split('_')
    dt = words[1].replace('tstep','')
    dt = dt.replace('.out','')
    if (float(dt) < maxDt and float(dt) > minDt):
      energyErrorMax = 0.0
      print 'Reading diagnostics in ' + fileName
      fileObj = open(fileName)
      lines = list(fileObj)
      fileObj.close()
      flag = False
      for line in reversed(lines):
        if ('Finished main timestepping loop' in line):
          flag = True
        if (flag and '(E-E0)/E0' in line):
          words = line.split()
          currentError = float(words[1])
          if (abs(energyErrorMax) < abs(currentError)):
            energyErrorMax = currentError
          break
      if (flag):
        energyErrorDict[dt] = energyErrorMax
      else:
        print '... skipping due to incomplete results ...'
  return energyErrorDict

def create_walltime_dict(globstr, minDt, maxDt, suffix, suffix_omit):
  walltimeDict = {}
  for fileName in glob.glob(globstr):
    if (suffix_omit is not None and suffix_omit in fileName):
      continue
    words = fileName.split('/')
    shortName = words[-1]
    words = shortName.split('_')
    dt = words[1].replace('tstep','')
    dt = dt.replace('.out','')
    if (float(dt) < maxDt and float(dt) > minDt):
      fileNameTiming = fileName.replace('.out','.err')
      print 'Reading timing info in ' + fileNameTiming
      fileObj = open(fileName)
      lines = list(fileObj)
      fileObj.close()
      flag = False
      for line in reversed(lines):
        if ('Finished main timestepping loop' in line):
          flag = True
          break
      if (not flag):
        print '... skipping due to incomplete results ...'
        continue
      fileObj = open(fileNameTiming)
      lines = list(fileObj)
      fileObj.close()
      words = lines[1].split()
      word = words[1] # ignore the word 'real'
      words = word.split(('m'))
      seconds = int(words[0])*60 + float(words[1].replace('s',''))
      walltimeDict[dt] = seconds
  return walltimeDict
