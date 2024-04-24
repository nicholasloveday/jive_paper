# The Jive Verification System and its Transformative Impact on Weather Forecasting Operations

This repository contains the code to generate figures used in the paper.

[Scores](https://github.com/nci/scores) is used to produce verification results.

The majority of the data used is on [zenodo](https://zenodo.org/records/11015211). The zip file size is 49.5 MB.

## Reproducing the results

### Prerequisites

You will require:
* [miniconda](https://docs.anaconda.com/free/miniconda/miniconda-install/)
* [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

### Running snakemake

```
snakemake -c 1 all
```

## Directory layout

  - `src/`: The notebooks and python code used to produce the figures and maps in the Jive paper.
  - `data/`: The data used by the notebooks.
  - `results/`: The figures and maps produced by running snakemake to reproduce report contents are placed here.
  - `logs/`: The logs from running snakemake. Only useful if there are issues running the code.
  - `environment.yml`: the conda environment file used by snakemake to run the notebooks.
  - `Snakefile`: the snakemake file to reproduce paper figures and maps.