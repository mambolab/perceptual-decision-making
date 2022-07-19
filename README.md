# Perceptual Decision Making
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6861111.svg)](https://doi.org/10.5281/zenodo.6861111)


This repository contains the code to reproduce the statistical analysis of the paper:


D'Andrea et al. (2022) Magnetoencephalographic spectral fingerprints differentiate evidence accumulation from saccadic motor preparation in perceptual decision-making.

In order to make the code running, you need to download FieldTrip vXXX and the dataset from https://doi.org/10.5281/zenodo.6861111.

## Contents
The folder ```/code/``` contains:
- ```code/alpha_onset-locked.m``` : the script to reproduce the time evolution of alpha power locked to the onset of the stimulus (Figure 4A).
- ```code/alpha_response-locked.m``` : the script to reproduce the time evolution of alpha power locked to the response of the subjects (Figure 4B).
- ```code/beta_response-locked.m``` : the script to reproduce the time evolution of beta power locked to the response of the subjects (Figure 5B).

You should modify the scripts by editing lines related to the path of the FieldTrip toolbox and the downloaded data.


