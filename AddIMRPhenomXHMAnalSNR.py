import numpy as np
import scipy.special as sc
import statistics as st
import random
import os
import pandas as pd
import pickle
import sys
import multiprocessing as mp
import scipy.stats as scst
from tqdm import tqdm
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.integrate import quad, simpson
from pycbc import waveform as wf
from multiprocessing import Pool, Manager, Process, Value, Array, cpu_count
from functools import partial
#from LISAhdf5 import LISAhdf5,ParsUnits

# Parameter defintion

df_nm = sys.argv[1]     # first parameter to be passed on calling, name of the .h5 catalogue file
df_key = sys.argv[2]    # second parameter to be passed on calling, key of the .h5 catalogue file
detector = sys.argv[3]  # third parameter to be passed on calling, detector on which to estimate the ideal SNR (e.g. LISA, aLIGO, aplusLIGO, ET)
SNR_cut = sys.argv[4]   # fifth parameter to be passed on calling, snr cut on events, cut all events from the catalogue with SNR smaller than x

c = 299792.46 # speed of light in Km/sec
G = 6.674*(10.**(-11.)) # Gravitational constant in m^3⋅kg^−1⋅s^−2
sol_mass = 1.988e30 # Value of the Solar Mass in Kg
GPc = 3.086e+25 # GPc to m conversion factor
H_0 = 67.8 # Hubble constant in Km/(s*MPc)
year = 365.25*24*60*60 # Years in second
n_jobs = cpu_count() - 4


if detector == 'LISA':
    f_min = 1.e-4 #Hz
    f_max = 0.5 #Hz
    frq_prec = 1.e-4 # Hz
    T_obs = 4. #Year
    arm_lenght = 2.5e9 # value of the arms lenght
    A = 2 / np.pi**(2 / 3) * np.sqrt(5 / 96) # Constant factor needed to compute LISA anal SNR 

if (detector == 'aLIGO' or detector == 'aplusLIGO'):
    f_min = 5. #Hz
    f_max = 5000. #Hz
    frq_prec = 0.1 # Hz
    T_obs = 1. #Year

if detector == 'ET':
    f_min = 0.1 #Hz
    f_max = 10000. #Hz
    frq_prec = 0.1 # Hz
    T_obs = 1. #Year
    
   


#Additional function definiton



if detector == 'LISA':

    def GetInitialFrequency(m1,m2,coal_T):
        M = m1 + m2
        ni = (m1*m2)/(M*M)
        res = ((256.*ni)/(5.*np.power((c*(10.**3.)),5.)))*np.power((G*M*sol_mass),(5./3.))*coal_T
        freq = (np.power(res,(-(3./8.)))/np.pi)    
        return freq

    def S_oms(freq):
        omega = (2.*np.pi*freq)/(c*1000)
        res = (15*1e-12)*omega*np.sqrt(1. + ((2.*1.e-3)/freq)**4.)
        return res**2.

    def S_acc(freq):
        res = ((3.*1.e-15)/(2.*np.pi*freq*c*1000))*np.sqrt(1. + ((0.4*1.e-3)/(freq))**2.)\
        *np.sqrt(1. + (freq/(8.*1.e-3))**4.)
        return res**2.

    def Snx1p5(freq, arm_lenght):
        omega = (2.*np.pi*freq)/(c*1000)
        res = 16.*(np.sin(omega*arm_lenght)**2.)*(S_oms(freq) + \
              (3. + np.cos(2.*omega*arm_lenght))*S_acc(freq)) 
        return res

    def S_Hx(freq, arm_lenght):
        omega = (2.*np.pi*freq)/(c*1000)
        return (20./3.)*(1. + 0.6*(omega*arm_lenght)**2.)*\
        ((Snx1p5(freq, arm_lenght))/((np.sin(omega*arm_lenght)**2.)*(4.*omega*arm_lenght)**2))

if detector == 'aLIGO':
    aLIGO_Sens = pd.read_csv('ALigoSens.txt', sep = "  ", engine = 'python')
    S_h = interp1d(aLIGO_Sens.Frequency, (aLIGO_Sens.aLIGO_Sh**2.), fill_value="extrapolate")

if detector == 'aplusLIGO':
    aplusLIGO_Sens = pd.read_csv('AplusDesign.txt', sep = "  ", engine = 'python')
    S_h = interp1d(aplusLIGO_Sens.Frequency, (aplusLIGO_Sens.aplusLIGO_Sh**2.), fill_value="extrapolate")

if detector == 'ET':
    ET_Sens = pd.read_csv('ETSens.txt', sep = "   ", engine = 'python')
    S_h = interp1d(ET_Sens.Frequency, (ET_Sens.ETSensD_Sum**2.), fill_value="extrapolate")



print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~~-~-~-~-~-~-~-~-~-~-~-~')
print('We start by loading the dataframe...')

BHCat = pd.read_hdf(df_nm, df_key)
print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~')
print('First of all we have to compute the frequency at the end of the observation time for the various events')


if detector == 'LISA':
    BHCat['CoalTime'] = BHCat['CoalTime']/(1. + BHCat.Redshift)
    BHCat.InitialFrequency[BHCat['InitialFrequency'] < f_min] = f_min
    BHCat['LISA_fend'] = (1./(1. + BHCat.Redshift))*GetInitialFrequency(BHCat.Mass1, BHCat.Mass2, (BHCat.CoalTime - T_obs)*year)
    BHCat.LISA_fend[BHCat['LISA_fend'] > f_max] = f_max


print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~~-~-~-~-~-~-~-~-~-~-~-~')

print('We will now define the function to estimate the analytical SNR using the IMRPhenomD Waveform')

def IMRPhenomD_AnalSNR(Dl, m1, m2):
    
    # Creating the waveform using pycbc
    
    wave = wf.get_fd_waveform(approximant = 'IMRPhenomD', mass1 = m1, mass2 = m2, distance = Dl, delta_f = frq_prec, f_lower = f_min/2., f_final = f_max) # Mass need to be redshifted, distance in megaparsec
    frq_span = wave[0].get_sample_frequencies() # 0 is for the + waveform, 1 is for x
    
    # Getting the mid point of the frequency and waveform array to integrate using trapeze method
    
    mid_frq = np.array(0.5*(frq_span[1::] + frq_span[:-1:]))
    mid_wave = np.array(0.5*(abs(wave[0])[1::] + abs(wave[0])[:-1:]))
    df = np.array(frq_span[1::] - frq_span[:-1:])
    
    # Cutting the waveform for frequency > LIGO
    
    if min(mid_frq) < f_min:
        good_idx = mid_frq > f_min
        mid_frq = mid_frq[good_idx]
        df = df[good_idx]
        mid_wave = mid_wave[good_idx]
    
    if max(mid_frq) > f_max:
        good_idx = mid_frq < f_max
        mid_frq = mid_frq[good_idx]
        df = df[good_idx]
        mid_wave = mid_wave[good_idx]
        
    if (detector == 'aLIGO' or detector == 'aplusLIGO'):   
    	SNR = 4.*(2./5.)*simpson((mid_wave**2.)/S_h(mid_frq), mid_frq)
    if detector == 'ET':
    	SNR = 4.*(2./5.)*(3./4.)*simpson((mid_wave**2.)/S_h(mid_frq), mid_frq)
    	
    return np.sqrt(SNR)
    
def IMRPhenomXHM_AnalSNR(Dl, m1, m2):
    
    # Creating the waveform using pycbc
    
    wave = wf.get_fd_waveform(approximant = 'IMRPhenomXHM', mass1 = m1, mass2 = m2, distance = Dl, delta_f = frq_prec, f_lower = f_min/2., f_final = f_max) # Mass need to be redshifted, distance in megaparsec
    frq_span = wave[0].get_sample_frequencies() # 0 is for the + waveform, 1 is for x
    
    # Getting the mid point of the frequency and waveform array to integrate using trapeze method
    
    mid_frq = np.array(0.5*(frq_span[1::] + frq_span[:-1:]))
    mid_wave = np.array(0.5*(abs(wave[0])[1::] + abs(wave[0])[:-1:]))
    df = np.array(frq_span[1::] - frq_span[:-1:])
    
    # Cutting the waveform for frequency > LIGO
    
    if min(mid_frq) < f_min:
        good_idx = mid_frq > f_min
        mid_frq = mid_frq[good_idx]
        df = df[good_idx]
        mid_wave = mid_wave[good_idx]
    
    if max(mid_frq) > f_max:
        good_idx = mid_frq < f_max
        mid_frq = mid_frq[good_idx]
        df = df[good_idx]
        mid_wave = mid_wave[good_idx]
        
    if (detector == 'aLIGO' or detector == 'aplusLIGO'):   
    	SNR = 4.*(2./5.)*simpson((mid_wave**2.)/S_h(mid_frq), mid_frq)
    if detector == 'ET':
    	SNR = 4.*(2./5.)*(3./4.)*simpson((mid_wave**2.)/S_h(mid_frq), mid_frq)
    	
    return np.sqrt(SNR)

# Function to parallelize the SNR calculation

def AnalSNR_ToPar(i):
    if (((i + 1)*10)%len(BHCat['Mass1']) == 0) :
        print('Percentage of completition : ',((i+1)*100.)/(len(BHCat['Mass1'])), '%', flush=True)

    if (BHCat['Mass1'][i]/BHCat['Mass2'][i])<1000.:
        return IMRPhenomXHM_AnalSNR(BHCat['Distance'][i]*1000., BHCat['Mass1'][i], BHCat['Mass2'][i])
    else:
        return IMRPhenomD_AnalSNR(BHCat['Distance'][i]*1000., BHCat['Mass1'][i], BHCat['Mass2'][i])

print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~~-~-~-~-~-~-~-~-~-~-~-~')

print('We can finally estimate the analytical SNR for the events of the catalogue')

ris = {}

if __name__ == '__main__':                                    
    # start the worker processes equals to n_jobs
    pool = Pool(n_jobs)
    ris = pool.map(AnalSNR_ToPar, range(len(BHCat['Mass1'])))
    pool.close()
    pool.join()

if detector == 'LISA':

    BHCat['AnalSNRIMRPD_LISA'] = np.array(ris)

if detector == 'aLIGO':

    BHCat['AnalSNRIMRPD_aLIGO'] = np.array(ris)

if detector == 'aplusLIGO':

    BHCat['AnalSNRIMRPD_aplusLIGO'] = np.array(ris)

if detector == 'ET':

    BHCat['AnalSNRIMRPD_ET'] = np.array(ris)

print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~~-~-~-~-~-~-~-~-~-~-~-~')

print('Cutting the Dataframe to one with only events having analytical SNR bigger than ', )    

print('Let me save the dataframe and we are done !')

if detector == 'LISA':
    CUTBHCat = BHCat[BHCat['AnalSNRIMRPD_LISA'] >= float(SNR_cut)]
    print('The number of sources with SNR bigger than '+str(SNR_cut)+' is '+str(len(CUTBHCat.Mass1)))
    CUTBHCat.to_hdf('LISAAnSNRIMRPDB'+SNR_cut+df_nm, df_key, mode='w')

if detector == 'aLIGO':
    CUTBHCat = BHCat[BHCat['AnalSNRIMRPD_aLIGO'] >= float(SNR_cut)]
    print('The number of sources with SNR bigger than '+str(SNR_cut)+' is '+str(len(CUTBHCat.Mass1)))
    CUTBHCat.to_hdf('aLIGOAnSNRIMRPDB'+SNR_cut+df_nm, df_key, mode='w')

if detector == 'aplusLIGO':
    CUTBHCat = BHCat[BHCat['AnalSNRIMRPD_aplusLIGO'] >= float(SNR_cut)]
    print('The number of sources with SNR bigger than '+str(SNR_cut)+' is '+str(len(CUTBHCat.Mass1)))
    CUTBHCat.to_hdf('aplusLIGOAnSNRIMRPDB'+SNR_cut+df_nm, df_key, mode='w')

if detector == 'ET':
    CUTBHCat = BHCat[BHCat['AnalSNRIMRPD_ET'] >= float(SNR_cut)]
    print('The number of sources with SNR bigger than '+str(SNR_cut)+' is '+str(len(CUTBHCat.Mass1)))
    CUTBHCat.to_hdf('ETAnSNRIMRPDB'+SNR_cut+df_nm, df_key, mode='w')
