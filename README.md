## Analysis for the publication: "Independent and Cumulative Fission Yield Covariance Matrices for 61 Compound Systems"

### Written by: Eric F. Matthews, Department of Nuclear Engineering, University of California - Berkeley

This repo contains the reproducable workflow that was used to generate the correlation and covariance matrices presented in the impending publication "Independent and Cumulative Fission Yield Covariance Matrices for 61 Compound Systems". This code and set of Jupyter notebooks are presented so that readers may review and scrutinize the dataset and analysis for this publication. 

The method for generating these matrices can be run for an individual compound system (see `yields/systems.txt` for a list) using the Jupyter Notebooks `fit_P_nu_A.ipynb`, `generate.ipynb`, and `ind_to_cml.ipynb`. To access these notebooks in a reproducable environment click the buttons below to go to a Binder.

`fit_P_nu_A.ipynb` - generates fission neutron multiplicity distributions as a function of mass number for the compound system.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/zenodo/10.5281/zenodo.3579556/?filepath=Analysis.ipynb)

`generate.ipynb` - generates independent yield correlation and covariance matrices for the compound system.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/zenodo/10.5281/zenodo.3579556/?filepath=Analysis.ipynb)

`ind_to_cml.ipynb` - converts the independent yield correlation and covariance matrices to cumulative yield matrices for the compound system.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/zenodo/10.5281/zenodo.3579556/?filepath=Analysis.ipynb)

Execute the following steps in order to run `Python` scripts to generate indepenent and cumulative yield correlation and covariance matrices for all 61 systems in the ENDF evaluation (Note: runtimes may be long):

* `cd` into the directory where you have cloned this repo. 

* Enter the command `python fit_P_nu_A.py`. This will generate the necessary fission neutron multiplicity distributions as a function of mass number for all 61 compound systems. 

* Enter the command `python generate.py`. This will generate the independent yield correlation and covariance matrices for all 61 compound systems.

* Enter the command `python ind_to_cml.py`. This will convert the independent yield correlation and covariance matrices to cumulative yield matrices for all 61 compound systems.


Please view `license.txt` to view the licensing and restrictions for this data and code. 
