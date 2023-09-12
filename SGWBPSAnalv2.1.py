#!/usr/bin/env python
# coding: utf-8

# <h1> SGWB Parameter Space Analysys</h1> 

# In the following, we'll implement a notebook that, given the required Probability Distribution Functions (PDF) describing a Black Hole population(BH), generates the figure of merit for the predicted analytical Stochastic Gravitational Wave Background(SGWB) in function of the amplitude and redshift range of the merging rate.
# First of all, we need to import some modules ! 

# In[1]:


import numpy as np
import scipy.special as sc
import statistics as st
import random
import os
import IPython
import pandas as pd
import pickle
import multiprocessing as mp
import scipy.stats as scst
#from tqdm import tqdm
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.integrate import quad, simpson
from scipy.stats import poisson
from scipy.special import gamma, hyp1f1
from pycbc import waveform as wf
from multiprocessing import Pool, Manager, Process, Value, Array
from functools import partial
#from LISAhdf5 import LISAhdf5,ParsUnits
#%matplotlib inline
#import matplotlib.pyplot as plt
#plt.style.use("seaborn-v0_8-whitegrid")


# <h2> Global Variables of the Simulation </h2>

# The global variables of the simulation will be set to :

# In[2]:


# Flags for the execution modes, initialized to false, check the the FLAG selection section for additional informations and initializing them !

# Flag needed to simulate the standard LIGO SOBBH population

SOBBH = False
SOBBH_Redevol = False
SOBBH_RSpike = False

# Flags for different types of PBH mass distribution

PBH = False
PBH_fRz = False
PBH_fRt = False
PBH_LogNormal = False
PBH_Gaussian = False

# Flags for fast execution

Compute_SNRMat = False

# Merger distribution parameters

T_obs = 1. # Lisa or LIGO estimated years of observation
efficiency = 1. # Lisa effective usefull time percentage for observations
max_tc = 10000. # max years of coalescence time for a BBH mergine event
frq_min = 3.e-5 # Hertz
frq_max = 0.5 # Maximum frequency in hertz to which the LISA detector is sensitive
frq_star = 1.e-2 # Value of the choosen frequency at which we estimate the SGWB to compare with other results
# The total time used to generate the merging events by multipling for the rate of merging will be set to max_tc


#General Constants 

c = 299792.46 # speed of light in Km/sec
G = 6.674*(10.**(-11.)) # Gravitational constant in m^3⋅kg^−1⋅s^−2
sol_mass = 1.988e30 # Value of the Solar Mass in Kg
MPc = 3.08567758149137*1e22 # meters
GPc = MPc*1e3 # meters
h = 0.678
H_0 = 67.8*1e3/MPc # Hubble constant in 1/(s)
Omega_m = 0.3 # Matter density in our universe
Omega_lambda = 0.7 # Cosmological constant density in our universe
Omega_k = 0. # Curvature density in our universe
rho_c = (3.*(H_0**2.))/(8.*np.pi*G) # Critical density in our universe
year = 365.25*24*60*60 # Years in second 
    
# Precision settings for the binned variables

n_jobs = 80
frq_res = 1e-6
frq_prec = int((frq_max - frq_min)/frq_res) + 1


# <h2> FLAG selection section </h2>

# To begin, we have to decide which types of sources we wish to simulate in our SGWB, in order to use the standard _LIGO-Virgo_ fiducial mass function we have to set the _SOBBH_ flag on :

# In[3]:


#SOBBH = True # If true, add the total SGWB strain produced by stellar origin binary black hole merging on the strain as estimated by LIGO


# if we wish to simulate PBH perturbations, instead, we have to choose between the two following different mass functions :

# In[4]:


PBH = True # If true, the FOM will be generated considering a PBH perturbation to the fiducial model


# In[5]:


PBH_LogNormal = True# This will simulate the Log Normal mass distribution for PBH described in ariv 2109.05836
#PBH_Gaussian = True # This will simulate a Gaussian PBH mass distribution, that can be used to generalize a bit the standard monocromatic mass function for PBH


# We may also decide to simulate the catalogue with a redshift evolving merging rate, by setting to true the Red_evol flag:

# In[6]:


SOBBH_Redevol = True # If true, the merging rate will evolve as a function of redshift, if false it will be assumed constant over the volume
#SOBBH_RSpike = True # If true, generate a spike of merging rate in a small redshift region, the population will follow the standard SOBBH one


# In the case of PBH, instead, we are gonna generate the merging rate following a simple power law with $k = 1.4$ as in [V. Atal et al](https://arxiv.org/abs/2201.12218), for what concern the value of $R_0$, we will consider it to be a certain fraction $f$ of the original _SOBBH_ merging rate. This mode can be run by activating the flag _PBH_fRz_:

# In[7]:


#PBH_fRz = True # If true, the merging rate would assumed to be the simple power law evolution of a fixed k, where the value of R0 would be given as a fraction f of the SOBBH one


# alternatively, we can describe the evolution with redshift of the merging rate using the model by [S. S. Bavera et al](https://arxiv.org/pdf/2109.05836.pdf), this would describe its evolution as a power law of the _Hubble Time_ at redshift $z$:

# In[8]:


PBH_fRt = True # If true, the merging rate would assumed to be a powerlaw of the Hubble time at redshift z, where the value of R0 would be given as a fraction f of the SOBBH one 


# The analytical SNR approximation can run using the following two modes for the inclination approximation:

# In[9]:


Inc_mode = 'max_i'  # Maximize the estimated SNR by assuming perfect inclination in respect to the detector
#Inc_mode = 'avg_i' # Average the estimated SNR by integrating over all the possible inclinations in respect to the detector


# Given the sources, we can furthermore decide if we wish to plot their merging rates in function of $z$ :

# In[10]:


#Plot_Rz = True # If true, generate plots at the end of the simulation


# we can also decide if recomputing the SNR matrix in the phase space or just use the one we already have as a file:

# In[11]:


#Compute_SNRMat = True


# <h2> Utility functions </h2>

# In the following, we are going to define some useful generical functions that will be needed to present the results.
# We will start with a function that can be used to convert matplotlib contour line to arrays.

# In[12]:


def get_contour_verts(cn):
    # Given a set of contour line, save them as a dictionary
    contours = []
    # for each contour line
    for cc in cn.collections:
        paths = []
        # for each separate section of the contour line
        for pp in cc.get_paths():
            xy = []
            # for each segment of that section
            for vv in pp.iter_segments():
                xy.append(vv[0])
            paths.append(np.vstack(xy))
        contours.append(paths)

    return contours


# <h2> Standard Cosmological Functions </h2>

# First of all, we'll need a function that allow us to convert from redshift to Gigaparsec :

# In[13]:


# Just a function to convert from Z to GPC using Hubble Law, in order to obtain the comoving distance

z_max = 1.e8
z_prec = 10000

def H(z):
    return np.sqrt((H_0**2.)*(Omega_m*((1. + z)**3.) + Omega_k*((1. + z)**2.) + Omega_lambda))

def Z_to_Gpc(z):
    
    # Remove the commented part to use a linear approximation of the Hubble law for low z 
    
    #if(zmax <= 0.5):
    #    return ((z*c*(10**(-3)))/(H_0)) # only valid for z < 0.5
    #else:
        
        span_z = np.linspace(0.,z,z_prec)
        span_mz = 0.5*(span_z[1::] + span_z[:-1:])
        
        # Beware, would fail if the span z is created in logarithmic scale !
        
        Int_Z = c*(10**(-3))*simpson(1./(H(span_mz)*(MPc/1.e3)), span_mz, axis=0)
    
        return Int_Z
    
def Z_to_HubbleTime(z):
    
    span_z = np.logspace(np.log10(z),np.log10(z_max),z_prec)
    span_mz = 0.5*(span_z[1::] + span_z[:-1:])
        
    # Beware, would fail if the span z is created in logarithmic scale !
        
    
    Int_Z = simpson(1./(H(span_mz)*(1. + span_mz)), span_mz, axis=0)
    
    return Int_Z
    
t_0 = Z_to_HubbleTime(1.e-12) # Can't put 0 as the logarithmic scale would fail        


# we also need a function that estimates the differential comoving volume in function of the redshift :

# In[14]:


#In the following function, the differential comoving volume in function of the redshift will be estimated as a spherical surface, it need to be integrated over dr to obtain the real volume 

def DeVC(z, Delta_z):
    r = dist_func(z)
    z_2 = z + 0.5*Delta_z
    z_1 = z_2 - Delta_z
    Delta_r = dist_func(z_2) - dist_func(z_1)
    return ((4.*np.pi*(r**2.)*Delta_r)/Delta_z)


# Another recurring parameter for inspiralling events is the Chirp Mass, given the mass of the two events involved in the binary merging :

# In[15]:


# Function that return the Chirp Mass of a binary merging event

def ChirpMass(m1,m2): 
   return ((m1*m2)**(3./5.))/((m1+m2)**(1./5.))


# together with the effective spin :

# In[16]:


#Function that given the spin and spin tilt gives the effective spin

def EffectiveSpin(m1, m2, a1, a2, st_a1, st_a2):
    res = (m1*a1*cos(st_a1))/(m1 + m2) + (m2*a2*cos(st_a1))/(m1 + m2) # Hope so, better to double check


# To represent the signal in units of omega, we are gonna need to convert our strain from units of _h_, _hc_, or _Flux_ to units of $\Omega_{gw}$ :

# In[17]:


def h_to_Omega(ran_frq, spectrum):
    # ran_frq and spectrum need to have same shape
    return ((4*((h*np.pi)**2.)*(ran_frq**3.)*spectrum)/(3.*(H_0**2)))


# In[18]:


def hc_to_Omega(ran_frq, spectrum):
    # ran_frq and spectrum need to have same shape
    return ((2*((h*np.pi)**2.)*(ran_frq**2.)*spectrum)/(3.*(H_0**2)))


# In[19]:


def Flux_to_Omega(ran_frq, Flux):
    # Flux need to be a constant expressing the whole integrated flux in function of z and m
    return ((ran_frq**(2./3.))/(rho_c*(c*1e3)**3))*Flux


# To conclude, we may define the energy loss during the inspiral phase, the procedure implemented is described in [P. Ajith et al.](https://arxiv.org/abs/0909.2867), even though in the LISA case we can use the assumption that all the waveforms appearing in detector, are in the pre-merger phase.
# We have :

# In[20]:


def dE_dnu(m1, m2, freq, a1 = -1, a2 =-1, st_a1 = -1, st_a2 = -1):
    # Compute the energy dispersed during an inspiral phase to a certain post-newtonian order in the pre-merger approximation
    # If the 4 parameters describing the spin configuration are not given, would automatically use only the first post-newtonian term
    
    Ch_M = ChirpMass(m1,m2)
    eta = m1*m2/(Ch_M**2)
    nu_prime = (np.pi*Ch_M*sol_mass*G*freq/c**3)**(1./3.)
    
    alpha_2 = -(323/224) + (451/168)*eta
    if(a1 == -1 or a2 == -1 or st_a1 == -1 or st_a2 == -1):
        alpha_3 = 0
    else:
        chi_spin = EffectiveSpin(m1, m2, a1, a2, st_a1, st_a2) 
        alpha_3 = ((27/8) - (11/6)*eta)*chi_spin
    f1 = 1 + alpha_2*nu_prime**2 + alpha_3*nu_prime**3
    
    res = (((G*np.pi)*(Ch_M**(5./3.)))/3)*(freq**(-1./3.))*(f1**2)
    
    return res
    


# while the total spectrum in Omega given by any BH channel expressed in energy spectral density, can be generally described using :

# In[21]:


def SpectralDens_to_OmegaGW(freq, F_nu):
    res = (freq/(rho_c * c**3))*F_nu


# <h2> LISA sensitivity curve </h2>

# In the following we are going to generate the LISA sensitivity curve, in order to compare our result with the properties of the instrument.
# The shape of the sensitivity curve in units of S can be defined using the following function :

# In[22]:


# return the value of the sensitivity curve S_h given the frequency

def get_SciRD(freq):
    S_2 = 3.6*10.**(-41.) #1/Hz
    S_1 = 5.76*(1. + (0.0004/freq)**2.)*10.**(-48.) # 1/(Hz*s^4)
    S_R = 1. + (freq/0.025)**2.
    S_h = (10./3.)*S_R*((S_1/(2.*np.pi*freq)**4.) + S_2)
    return S_h


# <h2> SOBBH LIGO All Channels SGWB </h2>

# In this section, we are going to initialize all the objects needed to compute the Stellar Origin Binary Black Hole merging(SOBBHm) SGWB.
# The probability distribution implemented for the variables of the events, will be taken from [B. P. Abbott T1](https://arxiv.org/abs/1811.12940), [B. P. Abbott T2](https://arxiv.org/abs/2010.14533).

# <h3> SOBBH - Characteristic strain functions </h3>

# The characteristic strain is given by :

# In[23]:


# Function to estimate the characteristic strain for SOBBHm events

SGWB_zmin = 1.e-2
SGWB_zmax = 1.e6
SGWB_zprec = 50000

def SOBBH_hcsqrd(frq, SOBBH_IntFac):
    return ((4.*G**(5./3.))/(3.*(np.pi**(1./3.))*(c*10**3)**2))*(frq**(-4./3.))*SOBBH_IntFac
            


# <h3> SOBBH - Mass distribution functions </h3>

# Let's start by defining the probability distribution in function of the masses.
# 
# We have :

# In[24]:


# Power law + Peak Mass Model of the paper arxiv 2010.14533

    
# Mass Distribution parameters (values taken from the results of arxiv 2111.03634)

SOBBH_m = 5.0 # + 0.86 - 1.7  minimum mass allowed by the popolation inference 
SOBBH_M = 100. # Solar Masses, taken from the prior of the paper as no real higher mass cutoff was estimated !
SOBBH_massprec = 500 # Binning density for the masses
SOBBH_alpha = 3.5 # + 0.6 - 0.56 Power law index
SOBBH_betaq = 1.1 # + 1.7 - 1.3  index for m2 power law in q
SOBBH_deltam = 4.9 #+ 3.4 - 3.2  used for the low mass smoothing function, generate peak at delta_m + m_min
SOBBH_lambdapeak = 0.038 # + 0.058 - 0.026 Intensity of the gaussian peak
SOBBH_mum = 34.0 # + 2.6 - 4.0 Location of the Gaussian peak in Solar Masses
SOBBH_sigmam = 5.69 # +4.28 - 4.34 Solar Masses, taken from arxiv 2010.14533 as no additional claim was made on last paper

# Defining of the smoothing function for m close to the minimimum mass

def SOBBH_MassSmoothing(m, SOBBH_m, SOBBH_deltam):
    if(m < SOBBH_m):
        return 0.
    else:
        if(m >= (SOBBH_m + SOBBH_deltam)):
            return 1.
        else:
            factor = np.exp((SOBBH_deltam/(m - SOBBH_m)) + (SOBBH_deltam/(m - SOBBH_m - SOBBH_deltam)))
            return 1./(factor + 1.)

# Defining a normalized power law distribution function, needed for the final distribution function        

def SOBBH_MassPowLaw(m, SOBBH_m, SOBBH_M, SOBBH_alpha, SOBBH_PLnorm):
    if(SOBBH_m <= m <= SOBBH_M):
        return (1./SOBBH_PLnorm)*(m**(-SOBBH_alpha))
    else:
        return 0.

# Estimating the Phase space of the Power law distribution using trapezoidal integration

def SOBBH_PowerLawPS(SOBBH_ranm1, SOBBH_m, SOBBH_M, SOBBH_alpha):
    
    mid_m = 0.5*(SOBBH_ranm1[1::] + SOBBH_ranm1[:-1:])
    
    
    if min(mid_m) < SOBBH_m:
        good_idx = mid_m >= SOBBH_m
        mid_m = mid_m[good_idx]
    
    if max(mid_m) > SOBBH_M:
        good_idx = mid_m <= SOBBH_M
        mid_m = mid_m[good_idx]
        
    ris = simpson(np.power(mid_m, (-SOBBH_alpha)), mid_m)

    return ris


# Defining a Gaussian distribution of the mass, needed for the final distribution function

def SOBBH_MassGauss(m, SOBBH_mum, SOBBH_sigmam, SOBBH_GSnorm):
    if (SOBBH_m <= m <= SOBBH_M):
        return ((1./(SOBBH_sigmam*np.sqrt(2.*np.pi)))*np.exp(-0.5*((m-SOBBH_mum)/SOBBH_sigmam)**2.))*1./SOBBH_GSnorm
    else:
        return 0.
    
def SOBBH_GaussPS(SOBBH_ranm1, SOBBH_m, SOBBH_M, SOBBH_mum, SOBBH_sigmam):

    mid_m = 0.5*(SOBBH_ranm1[1::] + SOBBH_ranm1[:-1:])
    
    
    if min(mid_m) < SOBBH_m:
        good_idx = mid_m >= SOBBH_m
        mid_m = mid_m[good_idx]
    
    if max(mid_m) > SOBBH_M:
        good_idx = mid_m <= SOBBH_M
        mid_m = mid_m[good_idx]
    
    ris = simpson((1./(SOBBH_sigmam*np.sqrt(2.*np.pi)))*np.exp(-0.5*((mid_m-SOBBH_mum)/SOBBH_sigmam)**2.), mid_m)
    
    return ris

# Defining the normalization constant for the q dependancy of the total mass distribution

def SOBBH_P2PS(SOBBH_ranm2, SOBBH_betaq, SOBBH_m, SOBBH_deltam):

    q_norm = np.linspace(0,1,len(SOBBH_ranm2))
    mid_m = 0.5*(SOBBH_ranm2[1::] + SOBBH_ranm2[:-1:])
    dm = (SOBBH_ranm2[1::] - SOBBH_ranm2[:-1:])

    for i in range(len(SOBBH_ranm2) - 1):

        q_norm[i] = 0.
        
        if (SOBBH_m <= mid_m[i] <= SOBBH_M):

            for j in range(i + 1):

                q_norm[i] += ((mid_m[j]/(mid_m[i]))**(SOBBH_betaq))*dm[j]*SOBBH_MassSmoothing(mid_m[j], SOBBH_m, SOBBH_deltam)
        else:
            q_norm[i] = 1.
        
    q_norm[len(SOBBH_ranm2) - 1] = q_norm[len(SOBBH_ranm2) - 2]

    return q_norm   


# Defining the proper Mass distribution function

def SOBBH_MassDistr(m1, m2, SOBBH_m, SOBBH_M, SOBBH_alpha, SOBBH_betaq, SOBBH_deltam, SOBBH_lambdapeak, SOBBH_mum, SOBBH_sigmam, SOBBH_PLnorm, SOBBH_GSnorm, SOBBH_qnorm, SOBBH_MassPS):

    if(m1 >= m2 and m1 <= SOBBH_M and m2 <= SOBBH_M):
        return ((1. - SOBBH_lambdapeak)*SOBBH_MassPowLaw(m1, SOBBH_m, SOBBH_M, SOBBH_alpha, SOBBH_PLnorm) + \
                SOBBH_lambdapeak*SOBBH_MassGauss(m1, SOBBH_mum, SOBBH_sigmam, SOBBH_GSnorm))*\
                SOBBH_MassSmoothing(m1, SOBBH_m, SOBBH_deltam)*\
                ((m2/m1)**(SOBBH_betaq))*(1./SOBBH_qnorm)*\
                SOBBH_MassSmoothing(m2, SOBBH_m, SOBBH_deltam)*(1./SOBBH_MassPS)
    else:
        return 0.


# Estimating the Phase space for the Model C Mass distribution function using trapezoidal integration

def SOBBH_ModCPS(SOBBH_ranm1, SOBBH_ranm2, SOBBH_m, SOBBH_M, SOBBH_alpha, SOBBH_betaq, SOBBH_deltam, SOBBH_lambdapeak, SOBBH_mum, SOBBH_sigmam, SOBBH_PLnorm, SOBBH_GSnorm, SOBBH_qnorm):

    mid_m = 0.5*(SOBBH_ranm1[1::] + SOBBH_ranm1[:-1:])
    dm = (SOBBH_ranm1[1::] - SOBBH_ranm1[:-1:])
    ris = 0.
    
    if min(mid_m) < SOBBH_m:
        good_idx = mid_m >= SOBBH_m
        mid_m = mid_m[good_idx]
        dm = dm[good_idx]
        SOBBH_qnorm = SOBBH_qnorm[good_idx]
    
    if max(mid_m) > SOBBH_M:
        good_idx = mid_m <= SOBBH_M
        mid_m = mid_m[good_idx]
        SOBBH_qnorm = SOBBH_qnorm[good_idx]

    for i in range(len(mid_m)):
        for j in range(i + 1):
                    q = mid_m[j]/mid_m[i] 
                    ris +=  dm[i]*dm[j]*\
                    ((1. - SOBBH_lambdapeak)*SOBBH_MassPowLaw(mid_m[i], SOBBH_m, SOBBH_M, SOBBH_alpha, SOBBH_PLnorm)\
                    + SOBBH_lambdapeak*SOBBH_MassGauss(mid_m[i], SOBBH_mum, SOBBH_sigmam, SOBBH_GSnorm))\
                    *SOBBH_MassSmoothing(mid_m[i], SOBBH_m, SOBBH_deltam)*(q**(SOBBH_betaq))\
                    *(1./SOBBH_qnorm[i])*SOBBH_MassSmoothing(mid_m[j], SOBBH_m, SOBBH_deltam)

    return ris
    


# <h3> SOBBH - Redshift dependent statistic </h3>

# We may now define, the various implemented merging rates as a function of the redshift _z_ as :

# In[25]:


# Function for the merging rate as described in the paper arxiv 2010.14533, the flag Red_evol will decide if adopting a merging rate the evolve with redshift (true) or not (false)

if SOBBH:
    SOBBH_z = 1.e-4 # to avoid SNR divergence due to extremely close events
    SOBBH_Zlog = 0.1 # max z value generated in log scale
    SOBBH_Zlin = 10.0 # max z value generated in lin scale
    SOBBH_zprec = 300 # Binning density for the redshift

SOBBH_k = 2.9 # + 1.7 - 1.8  VALID FOR REDSHIFT EVOLVING POWER LAW + PEAK MODEL MASS DISTRIBUTION, total agreement with SFR
SOBBH_CorrRz = (((1. + 0.2)**SOBBH_k)/(1. + ((1. + 0.2)/2.9)**(SOBBH_k + 2.9)))**(-1) # Normalization factor estimated at z = 0.2
    
# Defining the value of R0, the 0 index will have the value for redshift evolution merging rate, the 1 index would have the one for constant merging rate

SOBBH_R0 = {}
SOBBH_R0[0] = 28.3/(year*GPc**3.)# +13.9 - 9.1 GPC⁻³ yr^⁻¹ Value of the merging rate fitted at z = 0.2
SOBBH_R0[1] = 23.9/(year*GPc**3.) # +14.9 - 8.6 m^-3 s^-1 Middle value fitted using a Power Law + Peak mass model and a non evolving merging rate

def SOBBH_R(z):
    if(SOBBH_Redevol):
        # This merging rate was interpolated by Angelo Ricciardone and Daniel Figueroa based on arxiv 2010.14533 and arxiv 1907.12562
        return SOBBH_R0[0]*SOBBH_CorrRz*((1. + z)**SOBBH_k)/(1. + ((1. + z)/2.9)**(SOBBH_k + 2.9))
    else:
        return SOBBH_R0[1]

# If we wish to generate just a spike of events at a certain redshift range coming from a merging rate with fixed amplitude, we fix the following        
        
if SOBBH_RSpike:
    # These variables would set the location of the spike in the redshift range
    SOBBH_Rzmin = 2.
    SOBBH_Rzmax = 10.
    SOBBH_zprec = 80
    # These variables, will set the phase space to span in order to obtain the figure of merit 
    SOBBH_SpikeAmplMin = 0.
    SOBBH_SpikeAmplMax = 10000.
    SOBBH_SpikeAmplPrec = 500
    # This will be the value of the amplitude of the merging rate for the constant spike
    SOBBH_SpikeAmpl = 1.
    
    def SOBBH_R(z):
        # Pass the amplitude in units of 1/[yr*GPc], tipically the value is between [1, 200]
        return SOBBH_SpikeAmpl/(year*(GPc**3.))


# <h3> SOBBH - Number density of events</h3>

# We may finally define the distribution function for the number of events,in particular let's start with the function that describes the merging rate dependancy on the reference frame time: 

# In[26]:


def DtrDz(z):
    ris = 1./(H_0*(1. + z)*np.sqrt(Omega_m*((1. + z)**3.) + Omega_k*((1. + z)**2.) + Omega_lambda))
    return ris
        


# we can now integrate the mass and redshift dependant factor in order to get a constant that will multiply the frequency dependance of the characteristic strain function.
# After putting together all the integral dependant factors, we just have to integrate :

# In[27]:


if SOBBH :
    def SOBBH_IntND(i):
        
        ris = 0.
        
        if ((i*10)%len(SOBBH_ranz) == 0) :
            print('Percentage of completition : ',(i*100.)/(len(SOBBH_ranz)), '%')
                
        for j in range(len(SOBBH_ranm1)-1):
            for k in range(j + 1):
                deltas = (SOBBH_ranz[i + 1] - SOBBH_ranz[i])*(SOBBH_ranm1[j + 1] - SOBBH_ranm1[j])*(SOBBH_ranm2[k + 1] - SOBBH_ranm2[k])
                ris += deltas*SOBBH_R(0.5*(SOBBH_ranz[i + 1] + SOBBH_ranz[i]))*\
                            SOBBH_MassDistr(0.5*(SOBBH_ranm1[j + 1] + SOBBH_ranm1[j]), 0.5*(SOBBH_ranm2[k + 1] + SOBBH_ranm2[k]),\
                                           SOBBH_m, SOBBH_M, SOBBH_alpha, SOBBH_betaq, SOBBH_deltam, \
                                           SOBBH_lambdapeak, SOBBH_mum, SOBBH_sigmam, SOBBH_PLnorm, SOBBH_GSnorm, \
                                           SOBBH_qnorm[j], SOBBH_MassPS)*\
                            DtrDz(0.5*(SOBBH_ranz[i + 1] + SOBBH_ranz[i]))*\
                            ((ChirpMass(0.5*(SOBBH_ranm1[j + 1] + SOBBH_ranm1[j]),\
                                        0.5*(SOBBH_ranm2[k + 1] + SOBBH_ranm2[k]))*sol_mass)**(5./3.))\
                            /((1. + 0.5*(SOBBH_ranz[i + 1] + SOBBH_ranz[i]))**(1./3.)) 
                        
        return [0.5*(SOBBH_ranz[i + 1] + SOBBH_ranz[i]),ris]
                               
                        


# The number of events predicted at each $z$ can finally be obtained using:

# In[28]:


def SOBBH_NDistrib(z, m1, m2, Delta_z, SOBBH_qnorm):
    n = SOBBH_R(z)*(year*GPc**3.)*DeVC(z, Delta_z)*(T_obs/(1. + z)) \
        *SOBBH_MassDistr(m1, m2, SOBBH_m, SOBBH_M, SOBBH_alpha, SOBBH_betaq, SOBBH_deltam, SOBBH_lambdapeak, SOBBH_mum, SOBBH_sigmam, SOBBH_PLnorm, SOBBH_GSnorm, SOBBH_qnorm, SOBBH_MassPS)
    return n


# <h2> PBH population functions </h2>

# In this section, we are going to initialize the population functions needed to simulate the _Primordial Black Holes (PBH)_ SGWB.
# 

# <h3> PBH - Merging Rates </h3>

# In order to compute a PBH perturbation analysis, we are gonna define the PBH merging rate as a fraction of the fiducial LIGO merging rate.
# We start by defining the model presented in [V. Atal et al.](https://arxiv.org/abs/2201.12218) evolving with a simple power broken power law having $k = 1.1$ before $z_*$ and $k = 1.4$ after. 
# The model is as follows :

# In[29]:


if PBH_fRz :
    
    PBH_zmin = 0.5 # minimum value of the PBH merging rate
    PBH_zmax = 10.0 # max z value generated in lin scale
    PBH_zprec = 300 # Binning density for the redshift

    # Defining the value of R0, the 0 index will have the value for redshift evolution merging rate, the 1 index would have the one for constant merging rate

    PBH_R0 = 28.3/(year*GPc**3.) # +14.8 - 10.0 GPC⁻³ yr^⁻¹ Value of the merging rate fitted in at z = 0.2 in ligo population inference paper arxiv2111.03634
    PBH_CorrfRz = 1./(1. + 0.2)**SOBBH_k  # normalization factor needed to express the value of the LIGO merging rate in z=0
    
    def PBH_fR(z,f):
        if(z <= 1.):
            PBH_k = 1.1 # Value taken from arxiv 2201.12218, valid for small z !!
            return f*PBH_R0*PBH_CorrfRz*((1. + z)**PBH_k)
        else:
            PBH_k = 1.4 # Value taken from arxiv 2201.12218, valid for high z !!
            PBH_R1_corr = f*PBH_R0*PBH_CorrfRz*(((2.)**PBH_k) - ((2.)**1.1))
            return f*PBH_R0*PBH_CorrfRz*((1. + z)**PBH_k) - PBH_R1_corr
    
    def PBH_fRVec(z,f, tilt_low=1.1, tilt_high=1.4):
    
        to_ret = f*PBH_R0*PBH_CorrfRz
        z_fac  = (1. + z)**tilt_low
        z_fac[z > 1.] = (1. + z[z > 1])**tilt_high - (((2.)**tilt_high) - ((2.)**tilt_low))
        return to_ret*z_fac


# alternatively, we can use the same model described by [S. S. Bavera et al](https://arxiv.org/pdf/2109.05836.pdf) for the redshift evolution of the merging rate. The amplitude of the perturbation can still be parametrized using the $fR$ approach as in the previous model :

# In[30]:


if PBH_fRt:
    PBH_zmin = 0.5 # minimum value of the PBH merging rate
    PBH_zmax = 10.0 # max z value generated in lin scale
    PBH_zprec = 300 # Binning density for the redshift

    PBH_R0 = 28.3/(year*GPc**3.) # +14.8 - 10.0 GPC⁻³ yr^⁻¹ Value of the merging rate fitted in at z = 0.2 in ligo population inference paper arxiv2111.03634
    PBH_CorrfRz = 1./(1. + 0.2)**SOBBH_k  # normalization factor needed to express the value of the LIGO merging rate in z=0
    
    def PBH_fR(z,f):
        return f*PBH_R0*PBH_CorrfRz*((t_z(z)/t_0)**(-34./37.))
    
    def PBH_fRVec(z,f):
        return f*PBH_R0*PBH_CorrfRz*((t_z(z)/t_0)**(-34./37.))


# <h3> PBH - Gaussian Mass Distribution </h3>

# We can define a Gaussian mass distribution for PBH as :

# In[31]:


if PBH_Gaussian:
    PBH_m = 0. # Solar Masses. Minimum value assumed for the PBH mass
    PBH_M = 150. # Solar Masses. Maximum value assumed for the PBH mass
    PBH_massprec = 300 # Binning density for the mass range
    PBH_pdfmspan = np.linspace(0., 100., 150) # this span will be needed to compute the figures of merit
    PBH_mu = 5. # mean of the Gaussian distribution
    PBH_sigmam = 1. # sigma of the variance distribution
    PBH_sigmamspan = [1. ,5. ,10. ,15.] # Values of sigma_m to be spanned by the simulation
    
    # We use the following distribution for the mass, this tend to a monochromatic mass function for small values of sigma, yet it can be used to generalize the result to a wider subset of cases
    def PBH_MassGauss(m, PBH_mu, PBH_sigmam, PBH_GSnorm):
        return ((1./(PBH_sigmam*np.sqrt(2.*np.pi)))*np.exp(-0.5*((m-PBH_mu)/PBH_sigmam)**2.))*1./PBH_GSnorm
    
    # This function is to estimate the normalization constant
    def PBH_GaussPS(PBH_ranm, PBH_mu, PBH_sigmam):

        PBH_midm = 0.5*(PBH_ranm[1::] + PBH_ranm[:-1:])

        ris =  simpson(PBH_MassGauss(PBH_midm, PBH_mu, PBH_sigmam, 1.), PBH_midm)
            
        return ris


# <h3> PBH - Log-Normal Mass Distribution </h3>

# We can define a Log-Normal mass distribution for PBH as described in the paper by [S. S. Bavera et al ](https://arxiv.org/abs/2109.05836):

# In[32]:


if PBH_LogNormal:
    # We use the following distribution for the mass
    PBH_m = 0. # Solar Masses. Minimum value assumed for the PBH mass
    PBH_M = 150. # Solar Masses. Maximum value assumed for the PBH mass
    PBH_massprec = 300 # Binning density for the mass range
    PBH_pdfmspan = np.linspace(0, 100., 150) # this span will be needed to compute the figures of merit
    PBH_Mc = 34.54 # Solar masses, taken from the main paper by Bavera
    PBH_sigmamn = 0.41 # Taken from the main paper by Bavera
    PBH_sigmamnspan = [0.1 ,0.5 ,1. ,2.5] # Values of sigma_m to be spanned by the simulation
    
    def PBH_MassLNorm(m, PBH_Mc, PBH_sigmamn, PBH_LNnorm):
        return (1./(np.sqrt(2*np.pi)*PBH_sigmamn*m))*np.exp(-(np.log(m/PBH_Mc)**2)/(2*PBH_sigmamn**2))*1./PBH_LNnorm
    
    # This function is to estimate the normalization constant
    def PBH_LNnormPS(PBH_ranm, PBH_Mc, PBH_sigmamn):
        
        PBH_midm = 0.5*(PBH_ranm[1::] + PBH_ranm[:-1:])

        ris =  simpson(PBH_MassLNorm(PBH_midm, PBH_Mc, PBH_sigmamn, 1.), PBH_midm)
            
        return ris


# <h3> PBH - Number density of events</h3>

# In order to obtain the number of resolvable events at each z for the fiducial and sub-populations, we need first of all to know the SNR for an event in function of its parameter space $(z, m1, m2)$.
# We will hence define a function that gives the entry at each bin in the parameter space:

# In[33]:


if PBH :
    def AnalSNR_ParamSpace(z, mat):
    
        if (0.5*(PBH_ranz[z + 1] + PBH_ranz[z]) >= 0.5):
            for j in range(len(PBH_ranm1) - 1):
                for k in range(j + 1):
                    if (0.5*(PBH_ranm1[j + 1] + PBH_ranm1[j]))/(0.5*(PBH_ranm2[k + 1] + PBH_ranm2[k])) < 1000.:
                        vals = IMRPhenomXHM_AnalSNR(0.5*(PBH_ranz[z + 1] + PBH_ranz[z]), 0.5*(PBH_ranm1[j + 1] + PBH_ranm1[j]), 0.5*(PBH_ranm2[k + 1] + PBH_ranm2[k]))
                    else:
                        vals = IMRPhenomD_AnalSNR(0.5*(PBH_ranz[z + 1] + PBH_ranz[z]), 0.5*(PBH_ranm1[j + 1] + PBH_ranm1[j]), 0.5*(PBH_ranm2[k + 1] + PBH_ranm2[k]))
                    delta_val = pd.DataFrame([[z, j, k, vals[0], vals[1]],], columns = ['idx_z', 'idx_m1', 'idx_m2', 'aplusLIGO_val', 'ET_val'])
                    mat.append(delta_val)                
        
        if ((z*10)%(len(PBH_ranz) - 1) == 0) :
            print('Percentage of completition : ',(z*100.)/(len(PBH_ranz) - 1), '%', flush=True)
        
                               


# In[34]:


if PBH :
    def PBHvsFid_ResSrc(i, mat):
        #Initializing instance variables
        ris = 0.
        idx_shared = []
        LIGO_idxtofill = [count for count in range(len(PBH_frange))]
        ET_idxtofill = [count for count in range(len(PBH_frange))]
        
        # Estimating the normalization constant for the Gaussian Case 
        
        if PBH_Gaussian:
            PBH_GSnorm = PBH_GaussPS(PBH_ranm1, 0.5*(PBH_pdfmspan[i + 1] + PBH_pdfmspan[i]), PBH_sigmam)
        
        # Estimating the normalization constant for the LogNormal case
        
        if PBH_LogNormal:
            PBH_LNnorm = PBH_LNnormPS(PBH_ranm1, 0.5*(PBH_pdfmspan[i + 1] + PBH_pdfmspan[i]), PBH_sigmamn)
        
        for z in range(len(PBH_ranz) - 1):
            if (0.5*(PBH_ranz[z + 1] + PBH_ranz[z]) >= 0.5):
                if (len(LIGO_idxtofill) > 0 or len(ET_idxtofill) > 0 ):                
                    
                    # Initializing the resolvable sources at distance z for fiducial and subpopulation
                    
                    LIGO_fidres = 0.
                    LIGO_pertres = PBH_frange*0.
                    ET_fidres = 0.
                    ET_pertres = PBH_frange*0.
                    
                    # Spanning over all the mass couples
            
                    for j in range(len(PBH_ranm1) - 1):
                        for k in range(j + 1):
                            deltas = (PBH_ranz[z + 1] - PBH_ranz[z])*(PBH_ranm1[j + 1] - PBH_ranm1[j])*(PBH_ranm2[k + 1] - PBH_ranm2[k])
                            
                            # Due to triangulation of the phase space elements outside the diagonal should be counted twice
                            
                            if j == k:
                                sym_fac = 1.
                            else:
                                sym_fac = 2.
                                
                            if LIGO_SNRs[z][j][k] >= 8.:

                                LIGO_fidres += SOBBH_NDistrib(0.5*(PBH_ranz[z + 1] + PBH_ranz[z]), 0.5*(PBH_ranm1[j + 1] + PBH_ranm1[j]), 0.5*(PBH_ranm2[k + 1] + PBH_ranm2[k]), (PBH_ranz[z + 1] - PBH_ranz[z]), q_norm[j])*deltas

                                if PBH_Gaussian:
                                    LIGO_pertres += sym_fac*PBH_NDistrib(0.5*(PBH_ranz[z + 1] + PBH_ranz[z]), 0.5*(PBH_ranm1[j + 1] + PBH_ranm1[j]), 0.5*(PBH_ranm2[k + 1] + PBH_ranm2[k]), (PBH_ranz[z + 1] - PBH_ranz[z]), PBH_frange, 0.5*(PBH_pdfmspan[i + 1] + PBH_pdfmspan[i]), PBH_sigmam, PBH_GSnorm)*deltas

                                if PBH_LogNormal:
                                    LIGO_pertres += sym_fac*PBH_NDistrib(0.5*(PBH_ranz[z + 1] + PBH_ranz[z]), 0.5*(PBH_ranm1[j + 1] + PBH_ranm1[j]), 0.5*(PBH_ranm2[k + 1] + PBH_ranm2[k]), (PBH_ranz[z + 1] - PBH_ranz[z]), PBH_frange, 0.5*(PBH_pdfmspan[i + 1] + PBH_pdfmspan[i]), PBH_sigmamn, PBH_LNnorm)*deltas

                            if ET_SNRs[z][j][k] >= 8.:

                                ET_fidres += SOBBH_NDistrib(0.5*(PBH_ranz[z + 1] + PBH_ranz[z]), 0.5*(PBH_ranm1[j + 1] + PBH_ranm1[j]), 0.5*(PBH_ranm2[k + 1] + PBH_ranm2[k]), (PBH_ranz[z + 1] - PBH_ranz[z]), q_norm[j])*deltas

                                if PBH_Gaussian:
                                    
                                    ET_pertres += sym_fac*PBH_NDistrib(0.5*(PBH_ranz[z + 1] + PBH_ranz[z]), 0.5*(PBH_ranm1[j + 1] + PBH_ranm1[j]), 0.5*(PBH_ranm2[k + 1] + PBH_ranm2[k]), (PBH_ranz[z + 1] - PBH_ranz[z]), PBH_frange, 0.5*(PBH_pdfmspan[i + 1] + PBH_pdfmspan[i]), PBH_sigmam, PBH_GSnorm)*deltas

                                if PBH_LogNormal:
                                    
                                    ET_pertres += sym_fac*PBH_NDistrib(0.5*(PBH_ranz[z + 1] + PBH_ranz[z]), 0.5*(PBH_ranm1[j + 1] + PBH_ranm1[j]), 0.5*(PBH_ranm2[k + 1] + PBH_ranm2[k]), (PBH_ranz[z + 1] - PBH_ranz[z]), PBH_frange, 0.5*(PBH_pdfmspan[i + 1] + PBH_pdfmspan[i]), PBH_sigmamn, PBH_LNnorm)*deltas
        
                            
                    # Normalizing the results to make bin indipendent
                
                    LIGO_fidres /= (PBH_ranz[z + 1] - PBH_ranz[z])
                    ET_fidres /= (PBH_ranz[z + 1] - PBH_ranz[z])
                    LIGO_pertres /= (PBH_ranz[z + 1] - PBH_ranz[z])
                    ET_pertres /= (PBH_ranz[z + 1] - PBH_ranz[z])
                
                    # Saving the index for the aplusLIGO case
                    
                    idx_mod = [i for i,v in enumerate(LIGO_pertres) if v > 3.*np.sqrt(LIGO_fidres)]
                    idx_shared = list(set(idx_mod) & set(LIGO_idxtofill))
                    if len(idx_shared) > 0:
                        for idx in idx_shared:
                            delta_val = pd.DataFrame([[0, i, idx, 0.5*(PBH_ranz[z + 1] + PBH_ranz[z])],], columns = ['Detector', 'idx_pdfm', 'idx_Rf', 'v'])
                            mat.append(delta_val)
                        LIGO_idxtofill = [val for val in LIGO_idxtofill if val not in idx_shared]
                    
                   # Saving the index for the ET case
                    
                    idx_mod = [i for i,v in enumerate(ET_pertres) if v > 3.*np.sqrt(ET_fidres)]
                    idx_shared = list(set(idx_mod) & set(ET_idxtofill))
                    if len(idx_shared) > 0:
                        for idx in idx_shared:
                            delta_val = pd.DataFrame([[1, i, idx, 0.5*(PBH_ranz[z + 1] + PBH_ranz[z])],], columns = ['Detector', 'idx_pdfm', 'idx_Rf', 'v'])
                            mat.append(delta_val)
                        ET_idxtofill = [val for val in ET_idxtofill if val not in idx_shared]
                else:
                    break
        
        if ((i*10)%len(PBH_pdfmspan) == 0) :
            print('Percentage of completition : ',(i*100.)/(len(PBH_pdfmspan)), '%', flush=True)
        
                               


# The number of predicted events for the _PBH_ subpopulation at each $z$, can again be obtained using:

# In[35]:


def PBH_NDistrib(z, m1, m2, Delta_z, f, PBH_avg, PBH_sig, PBH_norm):
    if PBH_Gaussian:
        n = PBH_fR(z,f)*(year*GPc**3.)*DeVC(z, Delta_z)*(T_obs /(1. + z)) \
            *PBH_MassGauss(m1, PBH_avg, PBH_sig, PBH_norm)*PBH_MassGauss(m2, PBH_avg, PBH_sig, PBH_norm)
    if PBH_LogNormal:
        n = PBH_fR(z,f)*(year*GPc**3.)*DeVC(z, Delta_z)*(T_obs /(1. + z)) \
            *PBH_MassLNorm(m1, PBH_avg, PBH_sig, PBH_norm)*PBH_MassLNorm(m2, PBH_avg, PBH_sig, PBH_norm)
    return n


# In order to estimate the analytical SGWB instead, we can define the integrated mass factor function as :

# In[36]:


def MassFac_func(m1, m2, mu, sigma_m, Norm):
       if PBH_Gaussian:
           ris = PBH_MassGauss(m1, mu, sigma_m, Norm)*PBH_MassGauss(m2, mu, sigma_m, Norm)\
                  *((ChirpMass(m1,m2)*sol_mass)**(5./3.))
       if PBH_LogNormal:
           ris = PBH_MassLNorm(m1, mu, sigma_m, Norm)*PBH_MassLNorm(m2, mu, sigma_m, Norm)\
                  *((ChirpMass(m1,m2)*sol_mass)**(5./3.))
       return ris                       


# In[37]:


if PBH :
    def PBH_AnalSGWB(i):
        #Initializing instance variables
        
        z_fac = simpson(PBH_fRVec(SGWB_ranz,1.)*DtrDz(SGWB_ranz)/((1. + SGWB_ranz)**(1./3.)), SGWB_ranz)
        
        # Estimating the integral of the mass factor 
        
        if PBH_Gaussian:
            # Estimating the normalization constant for the Gaussian Case 
            PBH_GSnorm = PBH_GaussPS(SGWB_ranm1, 0.5*(PBH_pdfmspan[i + 1] + PBH_pdfmspan[i]), PBH_sigmam)
            Part_MF = MassFac_func(SGWB_ranm1[None,: ], SGWB_ranm2, 0.5*(PBH_pdfmspan[i + 1] + PBH_pdfmspan[i]), PBH_sigmam, PBH_GSnorm)
            mass_fac = simpson(simpson(Part_MF, x=SGWB_ranm2, axis=0), SGWB_ranm1 , axis=0)
        
        if PBH_LogNormal:
            # Estimating the normalization constant for the LogNormal case
            PBH_LNnorm = PBH_LNnormPS(SGWB_ranm1, 0.5*(PBH_pdfmspan[i + 1] + PBH_pdfmspan[i]), PBH_sigmamn)    
            Part_MF = MassFac_func(SGWB_ranm1[None,: ], SGWB_ranm2, 0.5*(PBH_pdfmspan[i + 1] + PBH_pdfmspan[i]), PBH_sigmamn, PBH_LNnorm)
            mass_fac = simpson(simpson(Part_MF, x=SGWB_ranm2, axis=0), SGWB_ranm1, axis=0)
        
        return [0.5*(PBH_pdfmspan[i + 1] + PBH_pdfmspan[i]),mass_fac*z_fac];
                               


# <h2> Detector sensitivity curves </h2>

# The list of available detectors is given by:
# 

# In[38]:


Det_names = ['aplusLIGO', 'ET']


# We start by implementing the analytical SNR on the [a+ LIGO](https://dcc.ligo.org/public/0149/T1800042/004/T1800042-v4.pdf) configuration, the detector frequency range for this configuration is given by :

# In[39]:


LIGO_f = 5. #Hz
LIGO_F = 5000. #Hz
LIGO_fprec = 5000

LIGO_fran = np.linspace(LIGO_f,LIGO_F, LIGO_fprec)
LIGO_mfran = 0.5*(LIGO_fran[1::] + LIGO_fran[:-1:])


# and its sensitivity curve can be loaded from the following .txt:

# In[40]:


aplusLIGO_Sens = pd.read_csv('AplusDesign.txt', sep = "  ", engine = 'python')
LIGO_Sh = interp1d(aplusLIGO_Sens.Frequency, (aplusLIGO_Sens.aplusLIGO_Sh**2.), fill_value="extrapolate")

LIGO_cmi = LIGO_mfran*0.

# Estimating the integrated factor of the analytical SNR estimator in order to fit an interpolator

for i in range(len(LIGO_fran) - 1):
            
    if (i == 0):
        LIGO_cmi[i] = (LIGO_fran[i + 1] - LIGO_fran[i])*((LIGO_mfran[i]**(-7./3.))/(LIGO_Sh(LIGO_mfran[i])))
    else:
        LIGO_cmi[i] = (LIGO_fran[i + 1] - LIGO_fran[i])*((LIGO_mfran[i]**(-7./3.))/(LIGO_Sh(LIGO_mfran[i])))\
        + LIGO_cmi[i - 1]

aplusAnSNR_IntFac = interp1d(LIGO_mfran, LIGO_cmi, fill_value="extrapolate")


# we can now implement the analytical SNR on the [ET](https://www.et-gw.eu/index.php) configuration, the detector frequency range in this case is given by:
# 

# In[41]:


ET_f = 0.1 #Hz
ET_F = 10000. #Hz
ET_fprec = 5000

ET_fran = np.linspace(LIGO_f,LIGO_F, LIGO_fprec)
ET_mfran = 0.5*(ET_fran[1::] + ET_fran[:-1:])


# while its sensitivity curve can be loaded from the following .txt:

# In[42]:


ET_Sens = pd.read_csv('ETSens.txt', sep = "   ", engine = 'python')
ET_Sh = interp1d(ET_Sens.Frequency, (ET_Sens.ETSensD_Sum**2.), fill_value="extrapolate")

ET_cmi = ET_mfran*0.

# Estimating the integrated factor of the analytical SNR estimator in order to fit an interpolator

for i in range(len(ET_fran) - 1):
            
    if (i == 0):
        ET_cmi[i] = (ET_fran[i + 1] - ET_fran[i])*((ET_mfran[i]**(-7./3.))/(ET_Sh(ET_mfran[i])))
    else:
        ET_cmi[i] = (ET_fran[i + 1] - ET_fran[i])*((ET_mfran[i]**(-7./3.))/(ET_Sh(ET_mfran[i])))\
        + ET_cmi[i - 1]

ETAnSNR_IntFac = interp1d(ET_mfran, ET_cmi, fill_value="extrapolate")


# <h2> Analytcal SNR estimator </h2>

# In order to understand the number of resolvable events for each population, we are now gonna implement an analytical SNR approximator as presented in [S. Babak et al.](https://arxiv.org/abs/2108.01167). The frequency of coalescence can be approximated by using the $f_{ISCO}$.

# In[43]:


def GetFisco(m1,m2):
    M = m1 + m2 # Masses need to be in the source frame !
    freq = (1./(6.*np.sqrt(6)*np.pi))*((c*1000)**3.)/(G*M*sol_mass) # Taken from eq 4.39 Maggiore
    return freq


# To fasten the code, we will choose two different approximations for the inclination. Depending on the choosen mode, the inclination of the events will be used differently in the waveform approximation:

# In[44]:


if Inc_mode == 'max_i':
    inc_fac = 8. # Signal maximized using best inclination
if Inc_mode == 'avg_i':
    inc_fac = 16./5. # Signal averaged over the possible inclinations


# We can finally define the following function to estimate the analytical SNR:

# In[45]:


# Shape of the waveform used for the SNR calculation
def Stas_WF(f, m1, m2, DL):
    Ch_M = ChirpMass(m1, m2)
    res = (2./(np.pi)**(2./3.))*np.sqrt(5./96.)*(((sol_mass*Ch_M*G)**(5./6.))/(DL*MPc))*(1./(c*1000)**(3./2.))*(f**(-7./6.))
    return res

# Estimating the analytical SNR

def AnalSNR(z, m1, m2, detector):
    #Redshifting the masses
    m1 = m1*(1. + z)
    m2 = m2*(1. + z)
    ChMass = ChirpMass(m1,m2)
    Dl = (1. + z)*dist_func(z)
    if detector == 'aplusLIGO':
        LIGO_fend = min(GetFisco(m1, m2), LIGO_F)
        res = np.sqrt(inc_fac*(aplusAnSNR_IntFac(LIGO_fend) - aplusAnSNR_IntFac(LIGO_f))*(((np.sqrt(5./96.)*((ChMass*sol_mass*G)**(5./6.)))/(Dl*GPc*np.pi**(2./3.)))**2.)*(1./(c*1000)**3.))
    if detector == 'ET':
        ET_fend = min(GetFisco(m1, m2), ET_F)
        res = np.sqrt(inc_fac*(ETAnSNR_IntFac(ET_fend) - ETAnSNR_IntFac(ET_f))*(((np.sqrt(5./96.)*((ChMass*sol_mass*G)**(5./6.)))/(Dl*GPc*np.pi**(2./3.)))**2.)*(1./(c*1000)**3.))
    
    return res 


# alternatively, we can also estimate the analytical SNR using the PyCBC waveforms:

# In[46]:


def IMRPhenomD_AnalSNR(z, m1, m2):
    # Return an array composed of [SNR_aplusLIGO, SNR_ET]
    Dl = (1. + z)*dist_func(z)*(1.e3)
    
    # Creating the waveform using pycbc
    
    wave = wf.get_fd_waveform(approximant = 'IMRPhenomD', mass1 = m1*(1. + z), mass2 = m2*(1. + z), distance = Dl, delta_f = 0.5, f_lower = 0.5) # Mass need to be redshifted, distance in megaparsec
    frq_span = wave[0].get_sample_frequencies() # 0 is for the + waveform, 1 is for x
    
    # Getting the mid point of the frequency and waveform array to integrate using trapeze method
    
    mid_frq = np.array(0.5*(frq_span[1::] + frq_span[:-1:]))
    mid_wave = np.array(0.5*(abs(wave[0])[1::] + abs(wave[0])[:-1:]))
    df = np.array(frq_span[1::] - frq_span[:-1:])
    
    # Cutting the waveform for frequency not in LIGO
    
    if min(mid_frq) < LIGO_f:
        good_idx = mid_frq > LIGO_f
        LIGO_frq = mid_frq[good_idx]
        LIGO_df = df[good_idx]
        LIGO_wave = mid_wave[good_idx]
    else:
        LIGO_frq = mid_frq
        LIGO_df = df
        LIGO_wave = mid_wave
    
    
    if max(LIGO_frq) > LIGO_F:
        good_idx = LIGO_frq < LIGO_F
        LIGO_frq = LIGO_frq[good_idx]
        LIGO_df = LIGO_df[good_idx]
        LIGO_wave = LIGO_wave[good_idx]
        
    # Cutting the waveform for frequency not in LIGO
    
    if min(mid_frq) < ET_f:
        good_idx = mid_frq > ET_f
        ET_frq = mid_frq[good_idx]
        ET_df = df[good_idx]
        ET_wave = mid_wave[good_idx]
    else:
        ET_frq = mid_frq
        ET_df = df
        ET_wave = mid_wave
    
    if max(ET_frq) > ET_F:
        good_idx = ET_frq < ET_F
        ET_frq = ET_frq[good_idx]
        ET_df = ET_df[good_idx]
        ET_wave = ET_wave[good_idx]
    
    #Now estimating the SNRs 
       
    aplusLIGO_SNR = 4.*(2./5.)*simpson((LIGO_wave**2.)/LIGO_Sh(LIGO_frq), LIGO_frq)
    ET_SNR = 4.*(2./5.)*(3./4.)*simpson((ET_wave**2.)/ET_Sh(ET_frq), ET_frq)
    
    return [np.sqrt(aplusLIGO_SNR), np.sqrt(ET_SNR)];


# In[47]:


def IMRPhenomXHM_AnalSNR(z, m1, m2):
    # Return an array composed of [SNR_aplusLIGO, SNR_ET]
    Dl = (1. + z)*dist_func(z)*(1.e3)
    
    # Creating the waveform using pycbc
    
    wave = wf.get_fd_waveform(approximant = 'IMRPhenomXHM', mass1 = m1*(1. + z), mass2 = m2*(1. + z), distance = Dl, delta_f = 0.5, f_lower = 0.5) # Mass need to be redshifted, distance in megaparsec
    frq_span = wave[0].get_sample_frequencies() # 0 is for the + waveform, 1 is for x
    
    # Getting the mid point of the frequency and waveform array to integrate using trapeze method
    
    mid_frq = np.array(0.5*(frq_span[1::] + frq_span[:-1:]))
    mid_wave = np.array(0.5*(abs(wave[0])[1::] + abs(wave[0])[:-1:]))
    df = np.array(frq_span[1::] - frq_span[:-1:])
    
    # Cutting the waveform for frequency not in LIGO
    
    if min(mid_frq) < LIGO_f:
        good_idx = mid_frq > LIGO_f
        LIGO_frq = mid_frq[good_idx]
        LIGO_df = df[good_idx]
        LIGO_wave = mid_wave[good_idx]
    else:
        LIGO_frq = mid_frq
        LIGO_df = df
        LIGO_wave = mid_wave
    
    
    if max(LIGO_frq) > LIGO_F:
        good_idx = LIGO_frq < LIGO_F
        LIGO_frq = LIGO_frq[good_idx]
        LIGO_df = LIGO_df[good_idx]
        LIGO_wave = LIGO_wave[good_idx]
        
    # Cutting the waveform for frequency not in LIGO
    
    if min(mid_frq) < ET_f:
        good_idx = mid_frq > ET_f
        ET_frq = mid_frq[good_idx]
        ET_df = df[good_idx]
        ET_wave = mid_wave[good_idx]
    else:
        ET_frq = mid_frq
        ET_df = df
        ET_wave = mid_wave
    
    if max(ET_frq) > ET_F:
        good_idx = ET_frq < ET_F
        ET_frq = ET_frq[good_idx]
        ET_df = ET_df[good_idx]
        ET_wave = ET_wave[good_idx]
    
    #Now estimating the SNRs 
       
    aplusLIGO_SNR = 4.*(2./5.)*simpson((LIGO_wave**2.)/LIGO_Sh(LIGO_frq), LIGO_frq)
    ET_SNR = 4.*(2./5.)*(3./4.)*simpson((ET_wave**2.)/ET_Sh(ET_frq), ET_frq)
    
    return [np.sqrt(aplusLIGO_SNR), np.sqrt(ET_SNR)];


# In[49]:


PBH_ranz = np.logspace(np.log10(PBH_zmin), np.log10(PBH_zmax), int(PBH_zprec/10))
ran_d = Z_to_Gpc(PBH_ranz)
dist_func = interp1d(PBH_ranz, ran_d, fill_value="extrapolate")
ran_d = 0.


# <h2> Setting of the analyzed phase space </h2>

# The simulation will be spanned over the following range of variables :

# In[ ]:


# Inizialization of the frequency range and spectrum

ran_frq = np.linspace(frq_min, frq_max, frq_prec)
sensitivity = get_SciRD(ran_frq)
spectrum = ran_frq * 0.
t_0 = Z_to_HubbleTime(1.e-12) # Can't put 0 as the logarithmic scale would fail

# Definition of the fiducial level of the SGWB and the various n-sigma values at the frequency of f_star in function of the noise level

SGWB_FidNoise = [
    1.8653859774892988e-12, # Interpolated from the analytical SGWB at frequency equal to f_star
    1.9210318549184913e-12, # Obtained by the 1-sigma confidence ellipses with respect to the fiducial noise level
    1.95702878750183e-12, # Obtained by the 2-sigma confidence ellipses with respect to the fiducial noise level
    1.9937002425570923e-12 # Obtained by the 3-sigma confidence ellipses with respect to the fiducial noise level
]

# Initialization of the redshift phase space

PBH_ranz = np.logspace(np.log10(PBH_zmin), np.log10(PBH_zmax), PBH_zprec)
ran_d = Z_to_Gpc(PBH_ranz)
dist_func = interp1d(PBH_ranz, ran_d, fill_value="extrapolate")
ran_d = 0.

SGWB_ranz = np.logspace(np.log10(SGWB_zmin), np.log10(SGWB_zmax), SGWB_zprec)

# Initialization of the SOBBH phase space
    
# Mass phase space

if SOBBH:
    SOBBH_ranm1 = np.logspace(np.log10(SOBBH_m),np.log10(SOBBH_m + 5. - (SOBBH_M - (SOBBH_m + 5))/SOBBH_massprec), int(SOBBH_massprec/10))
    SOBBH_ranm1 = np.append(SOBBH_ranm1, np.linspace(SOBBH_m + 5., SOBBH_M,SOBBH_massprec))
    SOBBH_ranm2 = SOBBH_ranm1
    SOBBH_PLnorm = SOBBH_PowerLawPS(SOBBH_ranm1, SOBBH_m, SOBBH_M, SOBBH_alpha)
    SOBBH_GSnorm = SOBBH_GaussPS(SOBBH_ranm1, SOBBH_m, SOBBH_M, SOBBH_mum, SOBBH_sigmam) 
    SOBBH_qnorm = SOBBH_P2PS(SOBBH_ranm2, SOBBH_betaq, SOBBH_m, SOBBH_deltam)
    SOBBH_MassPS = SOBBH_ModCPS(SOBBH_ranm1, SOBBH_ranm2, SOBBH_m, SOBBH_M, SOBBH_alpha, SOBBH_betaq, SOBBH_deltam, SOBBH_lambdapeak, SOBBH_mum, SOBBH_sigmam, SOBBH_PLnorm, SOBBH_GSnorm, SOBBH_qnorm)
    SOBBH_ranz = np.logspace(np.log10(SOBBH_z), np.log10(SOBBH_Zlog), SOBBH_zprec*2)
    SOBBH_ranz = np.append(SOBBH_ranz, np.linspace(SOBBH_Zlog + (SOBBH_ranz[(SOBBH_zprec*2) - 1] - SOBBH_ranz[(SOBBH_zprec*2) - 2]), SOBBH_Zlin, SOBBH_zprec*100))
    SOBBH_ranz = np.sort(SOBBH_ranz, kind = 'mergesort')
# Distance phase space 
    
if SOBBH_RSpike:
    SOBBH_ranz = np.linspace(SOBBH_Rzmin, SOBBH_Rzmax, SOBBH_zprec)
    SOBBH_ranampl = np.linspace(SOBBH_SpikeAmplMin, SOBBH_SpikeAmplMax, SOBBH_SpikeAmplPrec)

    

# Initialization of the PBH phase space        

if PBH:
    
    # Mass phase space
    
    PBH_ranm1 = np.linspace(PBH_m, PBH_M, PBH_massprec)
    PBH_ranm2 = PBH_ranm1
    SGWB_ranm1 = 0.5*(PBH_ranm1[1::] + PBH_ranm1[:-1:])
    SGWB_ranm2 = np.linspace(SGWB_ranm1[0]/2., SGWB_ranm1, int(PBH_massprec/2))
    PBH_frange = np.logspace(-3.,0.,400)
    
    # Matching SOBBH normalization precision to PBH parameter space
    
    SOBBH_PLnorm = SOBBH_PowerLawPS(PBH_ranm1, SOBBH_m, SOBBH_M, SOBBH_alpha)
    SOBBH_GSnorm = SOBBH_GaussPS(PBH_ranm1, SOBBH_m, SOBBH_M, SOBBH_mum, SOBBH_sigmam) 
    app = SOBBH_P2PS(PBH_ranm2, SOBBH_betaq, SOBBH_m, SOBBH_deltam)
    q_norm = 0.5*(app[1::] + app[:-1:])
    SOBBH_MassPS = SOBBH_ModCPS(PBH_ranm1, PBH_ranm2, SOBBH_m, SOBBH_M, SOBBH_alpha, SOBBH_betaq, SOBBH_deltam, SOBBH_lambdapeak, SOBBH_mum, SOBBH_sigmam, SOBBH_PLnorm, SOBBH_GSnorm, q_norm)
    
    if PBH_fRt:
        t_span = Z_to_HubbleTime(SGWB_ranz)
        t_z = interpolate.interp1d(SGWB_ranz, t_span)


# <h2> Main body of the simulation </h2>

# We may finally launch the pipeline to generate SGWB spectrum on every frequency bin of the frequency range, as well as the resolvable sources in function of $z$ both for the fiducial and subpopulation. 

# We can start by estimating the analytical SNR in each bin of the phase space by running the following function:

# In[ ]:


# Creating the SNR matrix and defining their names

LIGO_SNRs = np.zeros((len(PBH_ranz) - 1 ,len(PBH_ranm1) - 1, len(PBH_ranm2) - 1))
ET_SNRs = np.zeros((len(PBH_ranz) - 1 ,len(PBH_ranm1) - 1, len(PBH_ranm2) - 1))
LIGOSNRMat_filenm = Det_names[0]+'SNRMatZ'+str(PBH_zprec)+'M'+str(PBH_massprec)+'.pkl'
ETSNRMat_filenm = Det_names[1]+'SNRMatZ'+str(PBH_zprec)+'M'+str(PBH_massprec)+'.pkl'

if PBH and Compute_SNRMat:
    manager = Manager()
    print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~')
    d_par = manager.list()
    print('We are now estimating the SNR matrix over the parameter space')
    
    if __name__ == '__main__':                                    
        # start the worker processes equals to n_jobs
        pool = Pool(n_jobs)
        pool.map(partial(AnalSNR_ParamSpace, mat = d_par), range(len(PBH_ranz)-1))
        pool.close()
        pool.join()
    
    # Saving the SNR matrix over the phase space
    
    for count in range(len(d_par)):
        LIGO_SNRs[int(d_par[count]['idx_z'])][int(d_par[count]['idx_m1'])][int(d_par[count]['idx_m2'])] = float(d_par[count]['aplusLIGO_val'])
        ET_SNRs[int(d_par[count]['idx_z'])][int(d_par[count]['idx_m1'])][int(d_par[count]['idx_m2'])] = float(d_par[count]['ET_val'])

    
    with open(LIGOSNRMat_filenm,'wb') as file:
        pickle.dump(LIGO_SNRs, file)
    with open(ETSNRMat_filenm,'wb') as file:
        pickle.dump(ET_SNRs, file)
    print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~')


else:
    with open(LIGOSNRMat_filenm,'rb') as file:
          LIGO_SNRs = pickle.load(file)
    with open(ETSNRMat_filenm,'rb') as file:
          ET_SNRs = pickle.load(file)


# If we are analyzing a _PBH_ perturbation, the integrated factor in function of the _Mass PDF_ parameter can be estimated as :

# In[ ]:


if PBH and PBH_LogNormal:
    
    #Summing on the PBH background contribution in the case of a LogNormal PDF
    print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~')
    
    manager = Manager()
    LIGO_Zdom = {}
    ET_Zdom = {}
    d_ris = {}
    d_rist = {}
    
    for i in range(len(PBH_sigmamnspan)):
        print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~')
        PBH_sigmamn = PBH_sigmamnspan[i]
        print('Now simulating the Integrated factor for PBH (part ',i + 1,' of ',len(PBH_sigmamnspan),'), this can take some time !')
        d_par = manager.list()
        LIGO_app = np.zeros((len(PBH_pdfmspan) - 1, len(PBH_frange))) + 10.
        ET_app = np.zeros((len(PBH_pdfmspan) - 1, len(PBH_frange))) + 10.
        
        print('First simulating the background levels for the various sub-populations at the given sigma')
        
        if __name__ == '__main__':                                    
            # start the worker processes equals to n_jobs
            pool = Pool(n_jobs)
            d_ris[i] = pool.map(PBH_AnalSGWB, range(len(PBH_pdfmspan)-1))
            pool.close()
            pool.join()
        
        d_rist[i] = np.transpose(d_ris[i])
            
        print('Now estimating at which z the subpopulation will dominate in resolvable sources over the fiducial')
        
        if __name__ == '__main__':                                    
            # start the worker processes equals to n_jobs
            pool = Pool(n_jobs)
            pool.map(partial(PBHvsFid_ResSrc, mat = d_par), range(len(PBH_pdfmspan)-1))
            pool.close()
            pool.join()

        
        for count in range(len(d_par)):
            if Det_names[int(d_par[count]['Detector'])] == 'aplusLIGO':
                LIGO_app[int(d_par[count]['idx_pdfm'])][int(d_par[count]['idx_Rf'])] = float(d_par[count]['v'])
            if Det_names[int(d_par[count]['Detector'])] == 'ET':
                ET_app[int(d_par[count]['idx_pdfm'])][int(d_par[count]['idx_Rf'])] = float(d_par[count]['v'])
        
        LIGO_Zdom[i] = LIGO_app.transpose()
        ET_Zdom[i] = ET_app.transpose()
        print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~')

if PBH and PBH_Gaussian:
    
    #Summing on the PBH background contribution for the case of a Gaussian PDF with several sigma_m
    print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~')
    
    manager = Manager()
    LIGO_Zdom = {}
    ET_Zdom = {}
    d_ris = {}
    d_rist = {}
    
    for i in range(len(PBH_sigmamspan)):
        print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~')
        PBH_sigmam = PBH_sigmamspan[i]
        print('Now simulating the Integrated factor for PBH (part ',i + 1,' of ',len(PBH_sigmamspan),'), this can take some time !')
        d_par = manager.list()
        LIGO_app = np.zeros((len(PBH_pdfmspan) - 1, len(PBH_frange))) + 10.
        ET_app = np.zeros((len(PBH_pdfmspan) - 1, len(PBH_frange))) + 10.
        
        print('First simulating the background levels for the various sub-populations at the given sigma')
        
        if __name__ == '__main__':                                    
            # start the worker processes equals to n_jobs
            pool = Pool(n_jobs)
            d_ris[i] = pool.map(PBH_AnalSGWB, range(len(PBH_pdfmspan)-1))
            pool.close()
            pool.join()
        
        d_rist[i] = np.transpose(d_ris[i])
        print('Now estimating at which z the subpopulation will dominate in resolvable sources over the fiducial')
            
        if __name__ == '__main__':                                    
            # start the worker processes equals to n_jobs
            pool = Pool(n_jobs)
            pool.map(partial(PBHvsFid_ResSrc, mat = d_par), range(len(PBH_pdfmspan)-1))
            pool.close()
            pool.join()
            
        
        for count in range(len(d_par)):
            if Det_names[int(d_par[count]['Detector'])] == 'aplusLIGO':
                LIGO_app[int(d_par[count]['idx_pdfm'])][int(d_par[count]['idx_Rf'])] = float(d_par[count]['v'])
            if Det_names[int(d_par[count]['Detector'])] == 'ET':
                ET_app[int(d_par[count]['idx_pdfm'])][int(d_par[count]['idx_Rf'])] = float(d_par[count]['v'])
        
        LIGO_Zdom[i] = LIGO_app.transpose()
        ET_Zdom[i] = ET_app.transpose()
        print('-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~')


# the resulting dataset, can be reordered as :

# In[ ]:


if PBH and PBH_LogNormal:
    data = {}
    for i in range(len(PBH_sigmamnspan)):
        data[i] = {'PDF_p1' : d_rist[i][0], 'IntFac' : d_rist[i][1]}


# In[ ]:


if PBH and PBH_Gaussian:
    data = {}
    for i in range(len(PBH_sigmamspan)):
        data[i] = {'PDF_p1' : d_rist[i][0], 'IntFac' : d_rist[i][1]}


# <h2> Estimating the figure of merit grid values </h2>

# The values estimated on the previous subsection, now need to be spanned over the phase space in order to plot the figure of merits in function of the parameters. Let's start by initializing the grid, for the PBH subpopulations we have: 

# In[ ]:


if PBH and PBH_LogNormal:

    X = {}
    Y = {}
    Pert_PS = {}
    App = {}
    
    for i in range(len(PBH_sigmamnspan)):
        PBH_IntFacDP = pd.DataFrame(data[i])
        PBH_IntFacDP = PBH_IntFacDP.sort_values(["PDF_p1", "IntFac"], ascending=True)
        App[i] = PBH_IntFacDP['IntFac']
        X[i], Y[i] = np.meshgrid(PBH_IntFacDP.PDF_p1, PBH_frange)
        Pert_PS[i] = np.zeros((len(PBH_frange),len(PBH_IntFacDP.PDF_p1)))
    
if PBH and PBH_Gaussian:
    
    X = {}
    Y = {}
    Pert_PS = {}
    App = {}
    
    for i in range(len(PBH_sigmamspan)):
        PBH_IntFacDP = pd.DataFrame(data[i])
        PBH_IntFacDP = PBH_IntFacDP.sort_values(["PDF_p1", "IntFac"], ascending=True)
        App[i] = PBH_IntFacDP['IntFac']
        X[i], Y[i] = np.meshgrid(PBH_IntFacDP.PDF_p1, PBH_frange)
        Pert_PS[i] = np.zeros((len(PBH_frange),len(PBH_IntFacDP.PDF_p1)))


# In the two different cases, we can now fill the values of the grid as :

# In[ ]:


if PBH and PBH_LogNormal:
    for k in range(len(PBH_sigmamnspan)):
        for i in range(len(PBH_IntFacDP.PDF_p1)):
            for j in range(len(PBH_frange)):
                Pert_PS[k][j][i] = (hc_to_Omega(1e-2,SOBBH_hcsqrd(1e-2, PBH_frange[j]*App[k][i])) + SGWB_FidNoise[0])

if PBH and PBH_Gaussian:
    for k in range(len(PBH_sigmamspan)):
        for i in range(len(PBH_IntFacDP.PDF_p1)):
            for j in range(len(PBH_frange)):
                Pert_PS[k][j][i] = (hc_to_Omega(1e-2,SOBBH_hcsqrd(1e-2, PBH_frange[j]*App[k][i])) + SGWB_FidNoise[0])


# <h2> Saving the dataset </h2>

# We can save the data using :

# In[ ]:


# Saving the integrated factors for the SGWB Perturbation in each bin of parameter space

if PBH and PBH_LogNormal:
    fname = 'IntFacLNPDF.pickle'
if PBH and PBH_Gaussian:
    fname = 'IntFacGSPDF.pickle'

if PBH_fRt :
    fname = 'Rt' + fname
        
if PBH_fRz :
    fname = 'Rz' + fname
    
file_to_write = open(fname, "wb")
pickle.dump(Pert_PS, file_to_write)


# In[ ]:


# Saving the values of z at which the perturbation will overtake the fiducial models

if PBH and PBH_LogNormal:
    fname = 'RedDomLNPDF.pickle'
if PBH and PBH_Gaussian:
    fname = 'RedDomGSPDF.pickle'

if PBH_fRt :
    fname = 'aplusLIGO_Rt' + fname
        
if PBH_fRz :
    fname = 'aplusLIGO_Rz' + fname
    
file_to_write = open(fname, "wb")
pickle.dump(LIGO_Zdom, file_to_write)


# In[ ]:


# Saving the values of z at which the perturbation will overtake the fiducial models

if PBH and PBH_LogNormal:
    fname = 'RedDomLNPDF.pickle'
if PBH and PBH_Gaussian:
    fname = 'RedDomGSPDF.pickle'

if PBH_fRt :
    fname = 'ET_Rt' + fname
        
if PBH_fRz :
    fname = 'ET_Rz' + fname
    
file_to_write = open(fname, "wb")
pickle.dump(ET_Zdom, file_to_write)


# <h3> Setting alarm to inform when simulation is over </h3>

# In[ ]:


#file = 'Alarm-ringtone.mp3'
#os.system("mpg123 "+file)

