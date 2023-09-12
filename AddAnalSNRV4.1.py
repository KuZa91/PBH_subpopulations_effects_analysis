import numpy as np 
import scipy.special as sc 
import statistics as st 
import random 
import pandas as pd 
import sys
#from LISAhdf5 import LISAhdf5,ParsUnits 
#%matplotlib inline 
import matplotlib.pyplot as plt 
from tqdm import tqdm
from scipy import interpolate
from scipy.interpolate import interp1d
plt.style.use('seaborn-whitegrid') 

# Parameter defintion

df_nm = sys.argv[1]     # first parameter to be passed on calling, name of the .h5 catalogue file
df_key = sys.argv[2]    # second parameter to be passed on calling, key of the .h5 catalogue file
detector = sys.argv[3]  # third parameter to be passed on calling, detector on which to estimate the ideal SNR (e.g. LISA, aLIGO, aplusLIGO)
Inc_mode = sys.argv[4]  # fourth parameter to be passed on calling, define how to deal with events inclination, either max_i to maximize the resulting SNR or avg_i to use an average value 
SNR_cut = sys.argv[5]   # fifth parameter to be passed on calling, snr cut on events, cut all events from the catalogue with SNR smaller than x

c = 299792.46 # speed of light in Km/sec
G = 6.674*(10.**(-11.)) # Gravitational constant in m^3⋅kg^−1⋅s^−2
sol_mass = 1.988e30 # Value of the Solar Mass in Kg
GPc = 3.086e+25 # GPc to m conversion factor
H_0 = 67.8 # Hubble constant in Km/(s*MPc)
year = 365.25*24*60*60 # Years in second

if Inc_mode == 'max_i':
    inc_fac = 16.
if Inc_mode == 'avg_i':
    inc_fac = 32./5.

if detector == 'LISA':
    f_min = 1.e-4 #Hz
    f_max = 0.5 #Hz
    T_obs = 4. #Year
    arm_lenght = 2.5e9 # value of the arms lenght
    A = 2 / np.pi**(2 / 3) * np.sqrt(5 / 96) # Constant factor needed to compute LISA anal SNR 

if (detector == 'aLIGO' or detector == 'aplusLIGO'):
    f_min = 10. #Hz
    f_max = 5000. #Hz
    T_obs = 1. #Year

span_prec = 50000

#Additional function definiton


def ChirpMass(m1,m2): 
   return ((m1*m2)**(3./5.))/((m1+m2)**(1./5.))

if detector == 'LISA':

    def GetInitialFrequency(m1,m2,coal_T):
        M = (m1 + m2).astype('float32')
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

    def GetFisco(m1,m2):
        M = (m1 + m2).astype('float32')
        freq = (1./(6.*np.sqrt(6)*np.pi))*((c*1000)**3.)/(G*M*sol_mass) # Taken from eq 4.39 Maggiore
        return freq

if detector == 'aplusLIGO':
    aplusLIGO_Sens = pd.read_csv('AplusDesign.txt', sep = "  ", engine = 'python')
    S_h = interp1d(aplusLIGO_Sens.Frequency, (aplusLIGO_Sens.aplusLIGO_Sh**2.), fill_value="extrapolate")

    def GetFisco(m1,m2):
        M = (m1 + m2).astype('float32')
        freq = (1./(6.*np.sqrt(6)*np.pi))*((c*1000)**3.)/(G*M*sol_mass) # Taken from eq 4.39 Maggiore
        return freq



print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~~-~-~-~-~-~-~-~-~-~-~-~')
print('We start by loading the dataframe...')

BHCat = pd.read_hdf(df_nm, df_key)
print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~')
print('First of all we have to compute the frequency at the end of the observation time for the various events')

BHCat['Mass1'] = BHCat['Mass1']/(1. + BHCat.Redshift)
BHCat['Mass2'] = BHCat['Mass2']/(1. + BHCat.Redshift)
BHCat['ChM'] = ChirpMass(BHCat.Mass1,BHCat.Mass2)
if detector == 'LISA':
    BHCat['CoalTime'] = BHCat['CoalTime']/(1. + BHCat.Redshift)
    BHCat.InitialFrequency[BHCat['InitialFrequency'] < f_min] = f_min
    BHCat['LISA_fend'] = (1./(1. + BHCat.Redshift))*GetInitialFrequency(BHCat.Mass1, BHCat.Mass2, (BHCat.CoalTime - T_obs)*year)
    BHCat.LISA_fend[BHCat['LISA_fend'] > f_max] = f_max
if (detector == 'aLIGO' or detector == 'aplusLIGO'):
    BHCat['aLIGO_fend'] = GetFisco(BHCat.Mass1, BHCat.Mass2)
    BHCat.aLIGO_fend[BHCat['aLIGO_fend'] > f_max] = f_max


print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~~-~-~-~-~-~-~-~-~-~-~-~')

print('We will now estimate the cumulative distribution for the integral factor of the SNR estimator')

ran_frq = np.linspace(f_min,f_max, span_prec)
ran_cmi = np.linspace(f_min,f_max, span_prec - 1)*0.
ran_mfrq = np.linspace(f_min,f_max, span_prec - 1)*0.

if detector == 'LISA':

    for i in tqdm(range(len(ran_frq) - 1)):

        ran_mfrq[i] = 0.5*(ran_frq[i] + ran_frq[i + 1])
        
        if (i == 0):
            ran_cmi[i] = (ran_frq[i + 1] - ran_frq[i])*\
            ((ran_mfrq[i]**(-7./3.))/(S_Hx(ran_mfrq[i], arm_lenght)))
        else:
            ran_cmi[i] = (ran_frq[i + 1] - ran_frq[i])*\
            ((ran_mfrq[i]**(-7./3.))/(S_Hx(ran_mfrq[i], arm_lenght)))\
            + ran_cmi[i - 1]

if (detector == 'aLIGO' or detector == 'aplusLIGO'):

    for i in tqdm(range(len(ran_frq) - 1)):

        ran_mfrq[i] = 0.5*(ran_frq[i] + ran_frq[i + 1])
            
        if (i == 0):
            ran_cmi[i] = (ran_frq[i + 1] - ran_frq[i])*\
            ((ran_mfrq[i]**(-7./3.))/(S_h(ran_mfrq[i])))
        else:
            ran_cmi[i] = (ran_frq[i + 1] - ran_frq[i])*\
            ((ran_mfrq[i]**(-7./3.))/(S_h(ran_mfrq[i])))\
            + ran_cmi[i - 1]

print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~~-~-~-~-~-~-~-~-~-~-~-~')

print('Given the cumulative distribution, we can just define an interpolator to fastly compute that for any range')

IntFac = interp1d(ran_mfrq, ran_cmi, fill_value="extrapolate")

print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~~-~-~-~-~-~-~-~-~-~-~-~')

print('We can finally estimate the analytical SNR for the events of the catalogue')

if detector == 'LISA':

    BHCat['AnalSNR_LISA'] = np.sqrt(inc_fac*((A*((BHCat.ChM*sol_mass*G)**(5./6.))/(BHCat.Distance*GPc))**2.)*\
        (1./(c*1000)**3.)*(IntFac(BHCat.LISA_fend) - IntFac(BHCat.InitialFrequency)))

if detector == 'aLIGO':

    BHCat['AnalSNR_aLIGO'] = np.sqrt(inc_fac*(IntFac(BHCat.aLIGO_fend) - IntFac(f_min))*(((np.sqrt(5./24.)*((BHCat.ChM*sol_mass*G)**(5./6.)))/(BHCat.Distance*GPc*np.pi**(2./3.)))**2.)*(1./(c*1000)**3.))

if detector == 'aplusLIGO':

    BHCat['AnalSNR_aplusLIGO'] = np.sqrt(inc_fac*(IntFac(BHCat.aLIGO_fend) - IntFac(f_min))*(((np.sqrt(5./24.)*((BHCat.ChM*sol_mass*G)**(5./6.)))/(BHCat.Distance*GPc*np.pi**(2./3.)))**2.)*(1./(c*1000)**3.))

print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~~-~-~-~-~-~-~-~-~-~-~-~')

print('Cutting the Dataframe to one with only events having analytical SNR bigger than ', )    

print('Let me save the dataframe and we are done !')

if detector == 'LISA':
    CUTBHCat = BHCat[BHCat['AnalSNR_LISA'] >= float(SNR_cut)]
    print('The number of sources with SNR bigger than '+str(SNR_cut)+' is '+str(len(CUTBHCat.Mass1)))
    CUTBHCat.to_hdf('LISAAnSNRB'+SNR_cut+Inc_mode+df_nm, df_key, mode='w')

if detector == 'aLIGO':
    CUTBHCat = BHCat[BHCat['AnalSNR_aLIGO'] >= float(SNR_cut)]
    print('The number of sources with SNR bigger than '+str(SNR_cut)+' is '+str(len(CUTBHCat.Mass1)))
    CUTBHCat.to_hdf('aLIGOAnSNRB'+SNR_cut+Inc_mode+df_nm, df_key, mode='w')

if detector == 'aplusLIGO':
    CUTBHCat = BHCat[BHCat['AnalSNR_aplusLIGO'] >= float(SNR_cut)]
    print('The number of sources with SNR bigger than '+str(SNR_cut)+' is '+str(len(CUTBHCat.Mass1)))
    CUTBHCat.to_hdf('aplusLIGOAnSNRB'+SNR_cut+Inc_mode+df_nm, df_key, mode='w')
