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
              'IMKG232a': 32,
              'IMKG232b': 33,
              'IMKG242a': 34,
              'IMKG242b': 35,
              'IMKG243a': 36,
              'IMKG243b': 37,
              'IMKG252a': 38,
              'IMKG252b': 39,
              'IMKG253a': 40,
              'IMKG253b': 41,
              'IMKG254a': 42,
              'IMKG254b': 43,
              'IMKG254c': 44,
              'IMKG343a': 45,
              'IMKG343b': 46}

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
                 'IMKG232a': '-x',
                 'IMKG232b': '-x',
                 'IMKG242a': '-x',
                 'IMKG242b': '-x',
                 'IMKG243a': '-^',
                 'IMKG243b': '-^',
                 'IMKG252a': '-x',
                 'IMKG252b': '-x',
                 'IMKG253a': '-^',
                 'IMKG253b': '-^',
                 'IMKG254a': '-s',
                 'IMKG254b': '-s',
                 'IMKG254c': '-s',
                 'IMKG343a': '--^',
                 'IMKG343b': '--^'}

colorDict = {'KGU35-native': 'k',
             'KGU35': 'tab:blue',
             'ARS232-native': 'k',
             'ARS222': 'tab:blue',
             'ARS232': 'tab:orange',
             'SSP3333b': 'tab:green',
             'SSP3333c': 'k',
             'ARS233': 'tab:red',
             'IMKG232a': 'tab:purple',
             'IMKG232b': 'tab:brown',
             'IMKG242a': 'tab:pink',
             'IMKG242b': 'tab:gray',
             'IMKG252a': 'tab:olive',
             'IMKG252b': 'tab:cyan',
             'ARS343': 'tab:blue',
             'ARK324': 'tab:orange',
             'IMKG243a': 'tab:green',
             'IMKG243b': 'tab:red',
             'IMKG253a': 'tab:purple',
             'IMKG253b': 'tab:brown',
             'IMKG343a': 'tab:pink',
             'IMKG343b': 'tab:gray',
             'ARS443': 'tab:blue',
             'DBM453': 'tab:orange',
             'IMKG254a': 'tab:green',
             'IMKG254b': 'tab:red',
             'IMKG254c': 'tab:purple',
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

def create_walltime_dict(globstr, minDt, maxDt, suffix_omit):
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
def create_solution_dict(globstr, testName, varRef, indRef, tRef, minDt, maxDt, suffix_omit):
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
      directory = fileName.replace('.out','').replace('tsteptype','output_tsteptype')
      print 'Reading solution in ' + directory + '/' + testName
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
