# Code repository: *The era of the ARG: an empiricist's guide to ancestral recombination graphs*

## Table of Contents
- [Overview](#overview)
- [Author](#author)
- [Contents](#contents)
- [Content details](#content-details)
  - [figures](#figures)
  - [scripts](#project_scripts)

## Author

Name| Contact
:-----:|:-----
Alex Lewanski|[allewanski AT gmail DOT com](mailto:allewanski@gmail.com)
Michael Grundler (add email)
Gideon Bradburd (add email)

## Overview
This `arg_review` repo is associated with the following manuscript: *The era of the ARG: an empiricist's guide to ancestral recombination graphs* (future doi will be posted here). This repo contains the scripts associated with simulations and visualizations included in the paper. Questions about the specific content of this repo and/or the analysis details should be directed to Alex Lewanski (see [Author](#author) section). Higher level questions regarding this project can be directed to any of the authors.

## Contents
```
arg_review
|-- figures
`-- scripts
```


## Content details
### `figures`
The figures in the main text and supplementary materials of the manuscript (this still needs to added)

### `scripts`
The scripts and input files used in simulations and visualizations included in the paper.


## Instructions
All R scripts use the [here](https://here.r-lib.org) package for convenient import and export of files. Hence, to run any of these scripts, you will need to have the here package installed and loaded and make sure that the top level project directory contains the .here file (or use one of the other ways to specify the top level project directory). For further details, visit this [vignette](https://cran.r-project.org/web/packages/here/vignettes/here.html) and links cited within.

Script| Explanation
:-----:|:-----
[`review_sim_arg_intuition.py`](./scripts/review_sim_arg_intuition.py)|A Python script containing the code to complete the [msprime](https://tskit.dev/msprime/docs/stable/intro.html) simulations and perform initial processing of the output. To run the script, you will need the following libraries: msprime, pandas, and numpy.
[`sim_plot_creation.R`](./scripts/sim_plot_creation.R)|An R script that processes and visualizes the msprime simulations. To run this script, in addition the here package, you will need to have the following packages: dplyr, ggplot2, ape, TreeDist, ggforce, egg, ggtree, ggpubr, ggridges, dichromat, cowplot. 
