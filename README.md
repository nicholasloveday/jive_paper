# The Jive Verification System and its Transformative Impact on Weather Forecasting Operations

This repository contains the code to generate figures used in the paper.

[Scores](https://github.com/nci/scores) is used to produce verification results.

The majority of the data used is on [zenodo](https://zenodo.org/records/11015211).

## Reproducing the results

### Python environment

You will need:

* Python >= 3.10
* miniconda or mambaforge
* snakemake >= 5.2.4
* zenodo_get >= 1.3.4

#### Installing miniconda

See [Installing miniconda](https://docs.anaconda.com/free/miniconda/miniconda-install/)

#### Installing mamba (recommended)

```shell
conda install -n base -c conda-forge mamba
conda activate base
mamba init
source ~/.bashrc  # to activate mamba activate functionality
```

#### Creating environment from environment.yml

```shell
mamba env create --file environment.yml
mamba activate jive_paper
```

or if using conda:

```
conda env create --file environment.yml
conda activate jive_paper
```
