
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

- A copy of the notebooks developed in [Generating-a-BH-Merging-Catalogue](https://github.com/KuZa91/Generating-a-BH-Merging-Catalogue/tree/master) is required in order to generate the catalogs of the population. We recommend using the latest version which allows to generate _PBH_ sub-populations catalogs and be found in the linked [BH-SynthesisNotebook](https://github.com/KuZa91/Generating-a-BH-Merging-Catalogue/blob/master/BHCatalogV8.0.ipynb).

## How to run the PBH-subpopulation analysis ##

The steps required in order to reproduce the results presented in [Article to appear](https://www.lifewire.com/thmb/YKCp3LuI-r3vTaaQufVOETpI-CM=/1500x0/filters:no_upscale():max_bytes(150000):strip_icc()/google-404-error-0f9029ad5ea14b2db1cddb65d8188f69.png) are as follows:

1. The comparison among the fiducial and primordial merger rates can be computed using the [MergerRateComparison](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/MergerRateComparison.ipynb) notebook. We get the [RtzMergingRateComparison](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/RtzMergingRateComparison.png) figure.

2. We can compute the detectors horizon by running the [ComputeCosmicHorizon](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/ComputeCosmicHorizon.ipynb) notebook. This will generate the [DetectorHorizonVsz](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/DetectorHorizonVsz.png) figure. 

3. Use the [BH-SynthesisNotebook](https://github.com/KuZa91/Generating-a-BH-Merging-Catalogue/blob/master/BHCatalogV8.0.ipynb) to generate $\sim 10$ _LVK_ fiducial population catalogs. The catalogs are constructed in the source frame, and need to be converted into the detector frame for the analysis using the script [to_DetFrameV3.py](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/to_DetFrameV3.py). After converting them, the _SNR_ of the events on the various detectors can be generated by means of [IMRPhenomXHM](https://arxiv.org/abs/2004.06503) waveforms using the [XHM-SNR Script](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/AddIMRPhenomXHMAnalSNR.py);

4. The analytical _SGWB_ for the fiducial population can be generated using the [AnalyticalSGWB-Notebook](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/AnalyticSGWB.ipynb), this is based on the code presented in the [SGWB_Directory](https://github.com/KuZa91/Constructing-an-Analytic-SGWB). The notebook will present the results both graphically (see [SGWBPicture](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/SGWB%26Sensitivity_curves.png)) and numerically. The value of the SGWB-amplitude reported will be needed in order to perform the next step!

5. The $\sigma_\alpha$ value of the fiducial _SGWB_ amplitude can be estimated by running the [FisherMatrixSGWB-Notebook](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/FisherMatrixSGWB.ipynb). The results will be presented both graphically (see [ConfidenceEllipsesFigure](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/ConfidenceEllipses.png)) and numerically, the latters will be needed as an input for the next step !

6. We can finally run the [SGWBPSAnalysis-Notebook](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/SGWBPSAnalv2.1.ipynb) to generate the datasets required to compute the main plots of [Article to appear](https://www.lifewire.com/thmb/YKCp3LuI-r3vTaaQufVOETpI-CM=/1500x0/filters:no_upscale():max_bytes(150000):strip_icc()/google-404-error-0f9029ad5ea14b2db1cddb65d8188f69.png), the setups for the notebook are described in the latter. The generated datasets can then be plotted using the [PSPlotter](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/PSPlotter.ipynb), this will generate the results both for the [A+LIGOPS](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/1yrLIGORtSGWBPSLNPDF.png) and for the [ETPS](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/1yrETRtSGWBPSLNPDF.png).

7. We can choose some benchmark points from the obtained results and generate a catalog using the same steps as point $1.$  The details of the benchmarks points used in our analysis are presented in the [Article to appear](https://www.lifewire.com/thmb/YKCp3LuI-r3vTaaQufVOETpI-CM=/1500x0/filters:no_upscale():max_bytes(150000):strip_icc()/google-404-error-0f9029ad5ea14b2db1cddb65d8188f69.png).

8. The dataset of the analytical distribution of the _SNR_ for the considered benchmark points can be run using the [ResolvableDistEstimator-Notebook](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/ResolvableDistEstimator.ipynb), the results can be plotted using the [ResDistPlotter-Notebook](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/ResDistPlotter.ipynb). The figures will be automatically be generated both for the [A+LIGOSNRDist](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/apLIGOSNRB8ResDistGs.png) and [ETSNRDist](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/ETLIGOSNRB8ResDistGs.png).

9. By converting into _.pkl_ format using the [h5topickle-script](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/h5toPickle.py) the catalogs of point $1.$ and $5.$, we can run the [ResolvableSourcePlotter-Notebook](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/ResSourcePlotter.ipynb) to plot the resolvable source distribution of the generated catalogs. We will obtain both the statistical distribution of the fiducial sources on [A+LIGOFiducialDistribution](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/apLIGOSNRB8EventDistLogBinFidMod.png) and on [ETFiducialDistribution](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/ETSNRB8EventDistLogBinFidMod.png). Furthermore, the software will also plot the distribution of Fiducial + Sub-Population on [A+LIGOFiducial+SubPopDistribution](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/apLIGOSNRB8EventDistLogBinPlusLNSubPopNum.png) and [ETFiducial+SubPopDistribution](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/ETSNRB8EventDistLogBinPlusLNSubPopNum.png)


## Additional information about the execution of the notebooks

- The multiprocessing part of the code is set in order to run on a server, if you wish to run the codes on a personal computer it is reccomended to change:

  <code> n_jobs = 80  ->  (mp.cpu_count() - 4) </code>

  in the _Global variable of the simulation_ sections !

- The [XHM-SNR Script](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/AddIMRPhenomXHMAnalSNR.py) can be run on postprocessing on any _DetFrame_ catalog and can add to 
  the latter the _SNR_ column for various
  detector. The standard way to call the script is:
   
  <code> python AddIMRPhenomXHMAnalSNR.py _path-to-detframe-catalogue_ _catalogue-key_ _Detector-Name_ _SNR-CutTreshold_ </code>

  e.g

  <code> python AddIMRPhenomXHMAnalSNR.py DetFrameSubPop3.h5 PBH aplusLIGO 8 </code>

- The [BH-SynthesisNotebook](https://github.com/KuZa91/Generating-a-BH-Merging-Catalogue/blob/master/BHCatalogV8.0.ipynb) can be used to generate both catalogs for population coming from the Fiducial model or from a _PBH_ sub-population model. To this extent, read the informations and set the flags properly in the _FLAG selection section_ of said notebook. Once the flags are well defined, the particular parameters of the chosen mode can instead be set in the _Mass distribution functions_ section and in the _Redshift dependent statistic_ section. The spin variables for both populations will be generated in agreement with [S.Babak et al.](https://inspirehep.net/literature/2651157).

- The code for estimating the _SNR_ using the full parameters in the waveform is not available in this directory, for info, ask to _mauro.pieroni@cern.ch_

- On the first run, the [SGWBPSAnalysis-Notebook](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/SGWBPSAnalv2.1.ipynb) requires to compute a support _SNR_ mat and need to be run with the flag:

  <code> Compute_SNRMat = True </code> 

  in <code>In[11]</code>. After the first run, as long as the precision of the variable span is kept, it is possible to set said flag to false and just copy the generated _SNRMat.pkl_ files on the run directory. 
  
- When plotting using the [PSPlotter](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/PSPlotter.ipynb), the flags for the run needs to be in agreement with the generated datasets, this is also true for the set precision !

- In order to run the [ResolvableDistEstimator-Notebook](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/ResolvableDistEstimator.ipynb) it is necessary to use the same _SNRMat.pkl_ generated when using the [SGWBPSAnalysis-Notebook](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/SGWBPSAnalv2.1.ipynb).

- Before running the [ResolvableSourcePlotter-Notebook](https://github.com/KuZa91/PBH_subpopulations_effects_analysis/blob/main/ResSourcePlotter.ipynb) it is highly recommended to check the default paths of the code to understand the hierarchy that the notebook use in order to read the data! It is also necessary to compute the _SNR_ for both the Fiducial and Sub-pop catalogs, and both for $A^+$ _LIGO_ and $ET$ at the same time !

