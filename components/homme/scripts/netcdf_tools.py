from netCDF4 import Dataset
import numpy as np

class Quantity():
  def __init__(self, name):
    self.name = name
    self.runDictionary = {}
    self.refrunDictionary = {}
    self.errorDictionary = {}

  def getQ(self,run):
    if (run.name not in self.runDictionary.keys()):
      q = run.variables[self.name][:]
      q = np.ma.array(q,mask=run.mask)
      self.runDictionary[run.name] = np.ma.compress_cols(q)
    return self.runDictionary[run.name]

  def getQRef(self,run,refrun):
    if (run.name not in self.refrunDictionary.keys()):
      qRef = refrun.variables[self.name][:]
      qRef = np.ma.array(qRef,mask=run.mask)
      self.refrunDictionary[run.name] = np.ma.compress_cols(qRef)
    return self.refrunDictionary[run.name]

  def computeQError(self,run,refrun=None):
    if (run.name not in self.errorDictionary.keys()):
      if (run.name not in self.runDictionary.keys()):
        q = getQ(run)
      else:
        q = self.runDictionary[run.name]
      if (run.name not in self.refrunDictionary.keys()):
        qRef = getQRef(run,refrun)
      else:
        qRef = self.refrunDictionary[run.name]
      self.errorDictionary[run.name] = np.abs(q-qRef)
    return self.errorDictionary[run.name]

class Run():

    def __init__(self,runName,dt,limiterList=None,quantityList=None):
      self.name = runName
      self.dt = dt
      self.variables = {}
      print('Loading ' + runName)
      dataset = Dataset(runName,mode='r')
      if (quantityList is None):
        for limiter in limiterList:
          if (limiter.q.name not in self.variables.keys()):
            self.variables[limiter.q.name] = dataset.variables[limiter.q.name][:]
          self.variables[limiter.magnitudeName] = dataset.variables[limiter.magnitudeName][:]
        self.mask = self.variables[limiterList[0].magnitudeName].mask
      elif (limiterList is None):
        for quantity in quantityList:
          if (quantity.name not in self.variables.keys()):
            self.variables[quantity.name] = dataset.variables[quantity.name][:]
        # TODO: Make this more dynamic
        self.mask = dataset.variables['qc_conserv'][:].mask

      dataset.close()
