# MicroBooNE HNL analysis code

This repo contains python3 code used in the search for heavy neutral leptons (HNLs) decaying to electron-positron pairs and neutral pions in MicroBooNE (https://arxiv.org/abs/2310.07660).
The inputs required to run the analysis chain are .root file ntuples that have gone through the searchingfornues (PeLEE) analysis module (https://github.com/ubneutrinos/PELEE/). 
The notebooks here can be used to make various plots and ultimately calculate exclusion limits on the HNL mixing angle $|U_{\mu4}|^2$.
Some of the code here has been adapted from code written by Owen Goodwin in his HNL and Higgs portal-scalar analysis (https://arxiv.org/abs/2207.03840, https://github.com/ogoodwin505).

Most of the functions used in the notebooks are defined in the Functions.py file that is stored inside ./Utilities. Similarly most of the plotting functions are defined in the Plotter.py file. This means that these python files are very large with many functions defined inside. However, it makes the notebooks themselves less cluttered. If I were to rewrite this code I would probably add further .py files and make them smaller and more modular. 
Only a few of the notebooks in the repo are required to run the analysis, to see which refer to **Usage** below.

## Usage

Various python packages are required to run these notebooks, the HNL_env.yml file for the environment I used is included in the repo.
Use conda to install a new virtual environment from this .yml file. Note I am on linux, so there may be some differences for mac users.

To run the analysis from scratch and display a limit on $|U_{\mu4}|^2$ for Majorana HNLs, run the following notebooks:
1. 0\_Make\_pkls.ipynb
1. 2\_Preselection\_script.ipynb
1. 3\_BDT\_training\_optimizing.ipynb
1. 3.5\_BDT\_Results.ipynb
1. 3.6\_Signal\_DetVar\_sys.ipynb
1. 3.7\_Bkg\_DetVar\_sys.ipynb
1. 3.8\_Bkg\_reweight\_sys.ipynb
1. 3.9\_Final\_BDT\_output.ipynb
1. 4\_pyhf\_limit.ipynb
1. 5\_Limit\_plotting.ipynb

There are many other notebooks that can be used to make calculations and plots although they are not required for the basic analysis, that *basic analysis* being a search with NuMI Run1 and Run3 data for HNLs decaying to $e^+e^-$ and $\nu\pi^0$. 
For information on what individual notebooks do, please refer to the markdown text at the start of the individual notebooks. 

## Signal sample production

The signal samples used in this analysis were produced using a HNL generator written by Pawel Guzowski (pawel.guzowski@manchester.ac.uk, https://github.com/pguzowski). This generator has been added into uboonecode.
To produce more signal samples, a series of fhicl files need to be run.

## Issues to be aware of

- The uncertainties displayed on the Pre-selection and BDT input variable plots (as shown in David Marsden's thesis) were incorrect. They did not properly include the systematic uncertainties on the in-cryostat $\nu$ sample. This was due to a plotting bug which has been fixed as of Feb 2024. Therefore, if you are comparing Pre-selection or BDT input variable plots you make to those in the thesis, the uncertainties should be larger on the new plots. The uncertainties on the BDT output have always been correct, so they should be the same.
- The number of input variables used in the BDT models is probably higher than necessary. You can probably get the same separation after dropping several of these (especially the shower-fit ones, which are highly correlated to the track-fit ones). It might also be worth looking into calculating new variables and adding them into the list used for BDT training, for example (Total shower energy)/(Total track energy).
- There may be additional/updated limits that should be plotted in the 5\_Limit\_plotting.ipynb script. The standard ones plotted in the paper are not exhaustive.
- The HNL generator used has been validated in so far as the decay widths match those of the BeamHNL generator (https://github.com/kjplows/Generator) BELOW the threshold of decays to $\mu\pi$. If BeamHNL is used to produce samples from NuMI kaon flux directly then the output could be directly compared to the output of Pawel's generator and this would be a much more complete validation.   
- The limit calculation does not directly account for potential correlations between different types of uncertainties. This might be able to be *hacked in*, however, in its current iteration (version 0.7.1) pyhf does not natively allow for this.

## Authors

David Marsden, david.marsden@manchester.ac.uk OR d.marsden01@gmail.com

