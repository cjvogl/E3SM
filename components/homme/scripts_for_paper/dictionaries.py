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
              'IMEX-KG232a': 32,
              'IMEX-KG232b': 33,
              'IMEX-KG242a': 34,
              'IMEX-KG242b': 35,
              'IMEX-KG243a': 36,
              'IMEX-KG243b': 37,
              'IMEX-KG252a': 38,
              'IMEX-KG252b': 39,
              'IMEX-KG253a': 40,
              'IMEX-KG253b': 41,
              'IMEX-KG254a': 42,
              'IMEX-KG254b': 43,
              'IMEX-KG254c': 44,
              'IMEX-KG343a': 45,
              'IMEX-KG343b': 46}

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
                 'IMEX-KG232a': '-x',
                 'IMEX-KG232b': '-x',
                 'IMEX-KG242a': '-x',
                 'IMEX-KG242b': '-x',
                 'IMEX-KG243a': '-^',
                 'IMEX-KG243b': '-^',
                 'IMEX-KG252a': '-x',
                 'IMEX-KG252b': '-x',
                 'IMEX-KG253a': '-^',
                 'IMEX-KG253b': '-^',
                 'IMEX-KG254a': '-s',
                 'IMEX-KG254b': '-s',
                 'IMEX-KG254c': '-s',
                 'IMEX-KG343a': '--^',
                 'IMEX-KG343b': '--^'}

colorDict = {'KGU35-native': 'k',
             'KGU35': 'tab:blue',
             'ARS232-native': 'k',
             'ARS222': 'tab:blue',
             'ARS232': 'tab:orange',
             'SSP3333b': 'tab:green',
             'SSP3333c': 'k',
             'ARS233': 'tab:red',
             'IMEX-KG232a': 'tab:purple',
             'IMEX-KG232b': 'tab:brown',
             'IMEX-KG242a': 'tab:pink',
             'IMEX-KG242b': 'tab:gray',
             'IMEX-KG252a': 'tab:olive',
             'IMEX-KG252b': 'tab:cyan',
             'ARS343': 'tab:blue',
             'ARK324': 'tab:orange',
             'IMEX-KG243a': 'tab:green',
             'IMEX-KG243b': 'tab:red',
             'IMEX-KG253a': 'tab:purple',
             'IMEX-KG253b': 'tab:brown',
             'IMEX-KG343a': 'tab:pink',
             'IMEX-KG343b': 'tab:gray',
             'ARS443': 'tab:blue',
             'DBM453': 'tab:orange',
             'IMEX-KG254a': 'tab:green',
             'IMEX-KG254b': 'tab:red',
             'IMEX-KG254c': 'tab:purple',
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
