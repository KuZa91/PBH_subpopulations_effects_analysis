
# PBH subpopulations effects analysis
**Paolo Marcoccia<sup>1</sup>, Mauro Pieroni<sup>2</sup> and Germano Nardini<sup>1</sup>**

<sub>1. University of Stavanger, Institutt for Matematikk og Fysikk, Kjølv Egelands hus, 5.etg, E-blokk, 4021 Stavanger, Norway </sub>  
<sub>2. Department of Theoretical Phyics, CERN, 1211 Geneva 23, Switzerland</sub>  

## Introduction ##

In this repository, we give all the notebooks required in order to perform the analysis described in [Article to appear](https://www.lifewire.com/thmb/YKCp3LuI-r3vTaaQufVOETpI-CM=/1500x0/filters:no_upscale():max_bytes(150000):strip_icc()/google-404-error-0f9029ad5ea14b2db1cddb65d8188f69.png).
The analysis is based on the results obtained in [S.Babak et al.](https://inspirehep.net/literature/2651157) and [S.S.Bavera et al.](https://arxiv.org/abs/2109.05836), and will require to run the notebook and scripts presented in [Generating-a-BH-Merging-Catalogue](https://github.com/KuZa91/Generating-a-BH-Merging-Catalogue).  
The motivation for this study is to understand the effects of possible _Primordial Black Holes (PBHs)_ sub-populations over the fiducial _Black Holes (BHs)_ population observed by the _LIGO-Virgo-Kagra ([LVK](https://www.ligo.org/))_ [collaboration](https://www.ligo.org/), when the next generation of _Gravitational Wave (GW)_ detectors (such as [LISA](https://www.elisascience.org/) and [ET](https://www.einsteintelescope.nl/en/)) will run.
The fiducial _BH_ population is taken in agreement with what presented in 

## Analysis Details ##

Details of the analaysis can be found in our [Article to appear](https://www.lifewire.com/thmb/YKCp3LuI-r3vTaaQufVOETpI-CM=/1500x0/filters:no_upscale():max_bytes(150000):strip_icc()/google-404-error-0f9029ad5ea14b2db1cddb65d8188f69.png)

## Additional Material ##

In order to reproduce the results of our paper, an [Anaconda](https://www.anaconda.com/distribution/) environment with _python 2.7_ is needed.
The supplementary libraries needed, may be checked from the file [init_module.py](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/init_module.py).
In particular :

- The notebooks, require the installation of [PyCBC](https://pycbc.org/) v1.12.4 and [LALSuite](https://git.ligo.org/lscsoft/lalsuite) 6.49, which contains version 1.8.0 of the LALSimulation library used to generate the maximum likelihood waveform. Both of these libraries can be installed using [pip](https://pip.pypa.io/en/stable/) with the command:
```sh
pip install 'pycbc==1.12.4' 'lalsuite==6.49'
```
- Some function will require a _C_, _C++_ compiler in order to run, the [GCC](https://gcc.gnu.org/) compiler may be easily installed by running : 

```sh
sudo apt-get update
sudo apt-get install build-essential
```

- A new version of the file [res.py](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/res.py) is available in the code directory. The latter, need to be copied and pasted in the directory  
```sh
/home/*username*/anaconda2/lib/python2.7/site-packages
```

- The standard version of [Astropy](https://www.astropy.org/) will result into a problem while downloading the additional needed data due to a server that went down, a new version of the [iers.py](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/iers.py) was added in the [code](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/tree/master/Code) directory in order to fix the problem.
The new file should be copied in the directory :

```sh
/home/*username*/anaconda2/lib/python2.7/site-packages/astropy/utils/iers
```

- A copy of the notebooks developed by [Nielsen et al.](https://github.com/gwastro/gw150914_investigation) may be found in the [code](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/tree/master/Code) directory.

## How to run the Correlation Analysis ##

The Correlation Analysis, will be run for a single event at a time.
The steps to do in order to generate the results of our paper are :

1. Download the data of the two detectors _H1_,_L1_ for the event you wish to analyze from [GWOSC](https://www.gw-openscience.org/catalog/GWTC-1-confident/) and copy them in their respective directory, we used the _4096 s_, _4 KHz_ _.gwf_ files for our analysis;

2. Download the posteriors for the events of _O1_ from [Biwer et al.](https://github.com/gwastro/pycbc-inference-paper/tree/master/posteriors), the posteriors for _GW170104_ instead may be found [here](https://github.com/gwastro/o2-bbh-pe/tree/master/posteriors);

3. Run the [CreateResiduals.ipynb](https://github.com/gwastro/gw150914_investigation/blob/master/CreateResiduals.ipynb) to create the _residuals.hdf_ file, the previous notebook was built for _GW150914_, how to run that for other events will be briefly described in the **Additional information** section below. The correctness of the results for the _GW150914_ may be checked from the previous notebook, the results for the other events are written in the files [GW151012obtinf.txt](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/GW151012/151012obtinf.txt), [GW151226obtinf.txt](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/GW151226/151226obtinf.txt) and [GW170104obtinf.txt](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/GW170104/170104obtinf.txt) (There could be small variations in the value of _SNR_ and _Phase Shift_);

4. Run the [CorrVsTime.ipynb](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/CorrVsTime.ipynb) to generate the [CorrVsTime.csv](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/GW150914/CorrVsTime.csv) for the desired event.
The previous notebook is setted for [GW150914](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/tree/master/Code/GW150914) but the changing needed to be run for other events are pretty straight-forward !

5. Run the [PltCorrVsGPSTime.ipynb](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/PltCorrVsGPSTime.ipynb) to plot the [GW*event-name*CorrVsTime.png](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/GW150914/GW150914CorrVsTime.png) and find the _GPS_ Time of the max strain correlation.
The latter, will be needed to run the next point and is saved in [GW*event-name*data.py](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/GW150914/GW150914data.py) as <code> cort </code>;

6. Run the [PltCorrVsTimeLag.ipynb](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/PltCorrVsTimeLag.ipynb) to plot the equivalent of [Fig1_Fig2_Correlation.ipynb](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/Fig1_Fig2_Correlation.ipynb) as a multiple plot and generate the [AllMaxCorrVsTimeShift.png](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/AllMaxCorrVsTimeShift.png) figure.
The only correlations needed for our analysis will be generated by <code> In[4] </code> and <code> In[6]</code> of the [Fig1_Fig2_Correlation.ipynb](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/Fig1_Fig2_Correlation.ipynb) notebook, in order to create the [Fig2.png](https://arxiv.org/pdf/1811.04071.pdf) of the paper.

7. Run the [PltBackgrounds.ipynb](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/PltBackgrounds.ipynb) in order to generate the equivalent of [Fig3_Background.ipynb](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/Fig3_Background.ipynb) as a multiple plot. This notebook will describe the significance curve of correlation in function of the whitening band of the various events [AllBackgrounds.png](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/AllBackgrounds.png).
This process, would need the [strain]((https://www.gw-openscience.org/catalog/GWTC-1-confident/)) data of _GW150914_ to be run, and would be really slow !!
In the [Backgrounds](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/tree/master/Code/Backgrounds) directory, you would find the [tevent.py](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/Backgrounds/tevent.py) script, that once runned will automatically set up the psd starting time _5_ minutes after _GW150914_ coalescence time, as well as an already run version of the background for the various frequency bands.
If you're aiming to rerun the [Fig3_Background.ipynb](https://github.com/GravWaves-IMF/Correlation-Method-first-2019-/blob/master/Code/Fig3_Background.ipynb) by yourself, please note that <code> In[5] </code> should be replaced with the whitening function, and you need to generate one background for each whitening frequency range used in your event analysis ! 


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

