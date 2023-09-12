import numpy as np 
import scipy.special as sc 
import statistics as st 
import random 
import pandas as pd 
import sys
#from LISAhdf5 import LISAhdf5,ParsUnits 
#%matplotlib inline 


c = 299792.46 # speed of light in Km/sec
G = 6.674*(10.**(-11.)) # Gravitational constant in m^3⋅kg^−1⋅s^−2
sol_mass = 1.988e30 # Value of the Solar Mass in Kg
H_0 = 67.8 # Hubble constant in Km/(s*MPc)
year = 365.25*24*60*60 # Years in second 
f_max = 0.5
max_tc = 10000

#def GetInitialFrequency(m1,m2,coal_T):
#    M = m1 + m2
#   ni = (m1*m2)/(M*M)
#    res = ((256.*ni)/(5.*np.power((c*(10.**3.)),5.)))*np.power((G*M*sol_mass),(5./3.))*coal_T
#   return (np.power(res,(-(3./8.)))/np.pi)
  
#def TimeOutFrqRange(m1,m2,f_max):
#    M = m1 + m2
#    ni = (m1*m2)/(M*M)
#    res = (5.*np.power((c*(10.**3.)),5.))/(256.*ni*np.power((np.pi*f_max),(8./3.))*np.power((G*M*sol_mass),(5./3.)))
#    return res/year  
  

df_key = sys.argv[2]
df_nm = sys.argv[1]

print('We start by loading the dataframe...')

BHCat = pd.read_hdf(df_nm, df_key)

print('We are now redshifting the time and frequency variables !')

#BHCat['CoalTime'] = np.random.rand(len(BHCat.Mass1))*5000./(1. + BHCat.Redshift) # To source frame

#BHCat['InitialFrequency'] = GetInitialFrequency(BHCat.Mass1, BHCat.Mass2, BHCat['CoalTime']*year)*(1./(1. + #BHCat.Redshift))  # Generating if in source frame and then translating into detector frame

#BHCat['InBandTime'] = (BHCat['CoalTime']- TimeOutFrqRange(BHCat['Mass1'], BHCat['Mass2'], f_max*(1. + #BHCat.Redshift))*(1. + BHCat.Redshift))

#BHCat['CoalTime'] *= (1. + BHCat.Redshift) #  Back to detector frame


#print('We now have to redshift the distance !')

#BHCat['Distance'] *= (1. + BHCat.Redshift)

print('and finally the masses !')

BHCat['Mass1'] = BHCat['Mass1']*(1. + BHCat.Redshift)

BHCat['Mass2'] = BHCat['Mass2']*(1. + BHCat.Redshift)

print('Let me save the dataframe and we are done !')

BHCat.to_hdf('DetFrame'+df_nm, df_key, mode='w')
