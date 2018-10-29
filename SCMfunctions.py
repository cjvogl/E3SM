from numpy import log10, where

# physical constants
tboil = 373.16 # boiling temperature of water
mwwv = 18.016 # molecular weight of water vapor
mwdair = 28.966 # molecular weight of dry air

# derived values
eps = mwwv/mwdair

def q_sat(p, temp):
  # first compute saturation vapor pressure
  p_sat = 10.0**(-7.90298*(tboil/temp-1.0)+ \
          5.02808*log10(tboil/temp)- \
          1.3816e-7*(10.0**(11.344*(1.0-temp/tboil))-1.0)+ \
          8.1328e-3*(10.0**(-3.49149*(tboil/temp-1.0))-1.0)+ \
          log10(1013.246))*100.0
  # return value that corresponds to vapor pressure
  return where(p_sat < p, eps*p_sat / (p - (1.0-eps)*p_sat), 1.0)

def f(rh, rc, form='original'):
  RH0 = 0.8
  if (form == 'original'):
    return where(rh > RH0, where(rh < 1.0, ((rh-RH0)/(1.0-RH0))**2, 1.0), 0.0)
  if (form == 'alternate'):
    return where(rh > RH0, where(rh < 1.0, ((rh-RH0+rc)/(1.0-RH0+rc))**2, 1.0), (rc/(1.0-RH0+rc))**2)
    
    
      
    
