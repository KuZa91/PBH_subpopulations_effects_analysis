
# PBH subpopulations effects analysis
**Paolo Marcoccia<sup>1</sup>, Mauro Pieroni<sup>2</sup> and Germano Nardini<sup>1</sup>**

<sub>1. University of Stavanger, Institutt for Matematikk og Fysikk, Kjølv Egelands hus, 5.etg, E-blokk, 4021 Stavanger, Norway </sub>  
<sub>2. Department of Theoretical Phyics, CERN, 1211 Geneva 23, Switzerland</sub>  

## Introduction ##

In this repository, we give all the notebooks required in order to perform the analysis described in [Article to appear](https://www.lifewire.com/thmb/YKCp3LuI-r3vTaaQufVOETpI-CM=/1500x0/filters:no_upscale():max_bytes(150000):strip_icc()/google-404-error-0f9029ad5ea14b2db1cddb65d8188f69.png).
The analysis is based on the results obtained in [S.Babak et al.](https://inspirehep.net/literature/2651157) and [S.S.Bavera et al.](https://arxiv.org/abs/2109.05836), and will require to run the notebook and scripts presented in [Generating-a-BH-Merging-Catalogue](https://github.com/KuZa91/Generating-a-BH-Merging-Catalogue).  
The motivation for this study is to understand the effects of possible _Primordial Black Holes (PBHs)_ sub-populations over the fiducial _Black Holes (BHs)_ population observed by the _LIGO-Virgo-Kagra ([LVK](https://www.ligo.org/))_ collaboration, when the next generation of _Gravitational Wave (GW)_ detectors (such as [LISA](https://www.elisascience.org/) and [ET](https://www.einsteintelescope.nl/en/)) will run.
The fiducial _BH_ population is taken in agreement with what is presented in the latest _LVK_ [inference paper]((https://arxiv.org/abs/2111.03634)), and the results are extended at higher redshifts using the approach described in [S.Babak et al.](https://iopscience.iop.org/article/10.1088/1475-7516/2023/08/034). The _PBH_ sub-populations will instead be generated with a _merger rate_ in agreement with what is presented either in [S.S.Bavera et al.](https://arxiv.org/abs/2109.05836) or [V. Atal et al.](https://arxiv.org/abs/2201.12218), in particular, the _PBH_ abundance in $R(z = 0)$ will be set as a fraction of the _LVK_ fiducial population abundance.
We are going to test sub-population having a mass function that follows either a _Log-Normal_ distribution (as in [S.S.Bavera et al.](https://arxiv.org/abs/2109.05836)) or a _Gaussian_ distribution (see [V. Atal et al.](https://arxiv.org/abs/2201.12218)).

## Analysis details ##

Details of the analysis can be found in our [Article to appear](https://www.lifewire.com/thmb/YKCp3LuI-r3vTaaQufVOETpI-CM=/1500x0/filters:no_upscale():max_bytes(150000):strip_icc()/google-404-error-0f9029ad5ea14b2db1cddb65d8188f69.png).

## Additional material ##

In order to reproduce the results of our paper, an [Anaconda](https://www.anaconda.com/distribution/) environment with _python 3.x_ is needed.
This analysis, was performed using _python 3.8_. In order to estimate the _SNR_ for the generated events, some supplementary materials and libraries are needed.
In particular :

- The notebooks and SNR estimation scripts require the installation of [PyCBC](https://pycbc.org/). This library can be installed using [pip](https://pip.pypa.io/en/stable/) with the command:
```sh
pip install pycbc
```
- The SNR estimation scripts will need the sensitivity curves of the chosen detector in order to run, these files will be available in this directory for the [ALIGO](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/ALigoSens.txt), [A+LIGO](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/AplusDesign.txt) and [ET](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/ETSens.txt) designs.

- Some function will require a _C_, _C++_ compiler in order to run, the [GCC](https://gcc.gnu.org/) compiler may be easily installed by running : 

```sh
sudo apt-get update
sudo apt-get install build-essential
```

- A copy of the notebooks developed in [Generating-a-BH-Merging-Catalogue]([https://github.com/gwastro/gw150914_investigation](https://github.com/KuZa91/Generating-a-BH-Merging-Catalogue/tree/master)) is required in order to generate the catalogs of the population. We recommend using the latest version which allows to generate _PBH_ sub-populations catalogs, and be found in the linked [BH-SynthesisNotebook]([https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/tree/master/Code](https://github.com/KuZa91/Generating-a-BH-Merging-Catalogue/blob/master/BHCatalogV8.0.ipynb)).

## How to run the PBH-subpopulation analysis ##

The steps required in order to reproduce the results presented in [Article to appear](https://www.lifewire.com/thmb/YKCp3LuI-r3vTaaQufVOETpI-CM=/1500x0/filters:no_upscale():max_bytes(150000):strip_icc()/google-404-error-0f9029ad5ea14b2db1cddb65d8188f69.png) are as follows:

1. Use the [BH-SynthesisNotebook]([https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/tree/master/Code](https://github.com/KuZa91/Generating-a-BH-Merging-Catalogue/blob/master/BHCatalogV8.0.ipynb)) to generate $\sim 10$ _LVK_ fiducial population catalogs. The catalogs are constructed in the source frame, and need to be converted into the detector frame for the analysis using the script [to_DetFrameV3.py](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/to_DetFrameV3.py). After converting them, the _SNR_ of the events on the various detectors can be generated by means of [IMRPhenomXHM](https://arxiv.org/abs/2004.06503) waveforms using the [XHM-SNR Script](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/AddIMRPhenomXHMAnalSNR.py);

2. 


## Additional information about the execution of the notebooks

- In order to run the [CreateResiduals.ipynb](https://github.com/gwastro/gw150914_investigation/blob/master/CreateResiduals.ipynb) for events different from _GW150914_, one should first of all move to the choosed event directory using <code> cd _directoryname_ </code>, then check that all the additional data said in points __1)__ and __2)__ of the __How to run__ section was correctly downloaded, and finally launch the two scripts [init_module.py](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/init_module.py) and [GW*event-name*data.py](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/GW151012/GW151012data.py) using the magic command <code> %run _scriptname_ </code>.
Once all these preparation steps are done, the script [CreateResiduals.ipynb] (https://github.com/gwastro/gw150914_investigation/blob/master/CreateResiduals.ipynb) may be easily run for different events by replacing <code> In[3] </code> with the version using auxiliary variables :

  <code> data_files = {'H1': 'H-H1_' + hd_nm, 'L1': 'L-L1_' + hd_nm} </code> 

  and replacing <code> In[9] </code> with line :

  <code> fp = InferenceFile(pst_nm,'r') </code> 
  
- Depending on how old your strain data is, it may be saved either with the acronym _LOSC_ or _GWOSC_, differences in the data sets may results in errors while trying to load the data, so always double check that the name of the files are matching. In particular, the strain arrays inside the data files would have their channel named differently according to what stated at the [GWOSC](https://www.gw-openscience.org/o2_details/), to solve the channel name problem one simply need to replace the '_LOSC-STRAIN_' string passed in <code> In[5] </code> and <code> In[7] </code> at [CreateResiduals.ipynb](https://github.com/gwastro/gw150914_investigation/blob/master/CreateResiduals.ipynb) using '<em>GWOSC-4KHZ_R1_STRAIN</em>'.

- The strain data of event [GW151226](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/tree/master/Code/GW151226) has got _NAN_ values at the start, this will result in numerical errors during the _PSD_ estimation, in particular, you would notice that from a _NAN_ result from <code> Out[20] </code> of [CreateResiduals.ipynb](https://github.com/gwastro/gw150914_investigation/blob/master/CreateResiduals.ipynb). To avoid problems one simply need to add the two labels <code> start_time=psd_start_time-pad_data,
                                     end_time=psd_end_time+pad_data </code> when running reading commands like <code> In[7] </code> in  [CreateResiduals.ipynb](https://github.com/gwastro/gw150914_investigation/blob/master/CreateResiduals.ipynb).
                                     
- If you wish to reproduce the multiple subplots, as done for [AllMaxCorrVsTimeShift.png](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/AllMaxCorrVsTimeShift.png) and [AllBackgrounds.png](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/AllBackgrounds.png), a guide on how to generate multiple subplots may be found in the [Jake VanderPlas D.S Handbook](https://jakevdp.github.io/PythonDataScienceHandbook/04.08-multiple-subplots.html)                                    

