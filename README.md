# Snakemake workflow: Parity modifies the effect of genetic variants associated with gestational duraiton but not birth weight  

The manuscript can be found [here](https://www.medrxiv.org/content/10.1101/2025.06.17.25329777v1).

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.6.2-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/{{cookiecutter.repo_name}}.svg?branch=master)](https://travis-ci.org/snakemake-workflows/{{cookiecutter.repo_name}})

The scripts in this folder present the Snakemake workflow for a project exploring if parity modifies the maternal genetic effects on gestational duration or birth weight (Gene x Parity interactions). Project uses data from the Norwegian Mother, Father and Child Cohort Study (Moba).  

Each respective folder contains scripts necessary to run the whole pipeline. REGENIE (v3.2.5) and LD score regression (LDSC) (v 1.0.1) used in this pipeline did not exist as a conda environment and need to be installed to be able to run the whole pipeline. 


## Authors

* Karin Ytterberg (@ykarin97})
* Hedvig Sundelin (@hedvigsun)
* Pol Sole-Navais (@psnavais)
* Xiaoping Wu (@Xiaoping-Wu)

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).  
