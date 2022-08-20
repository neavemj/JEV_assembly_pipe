# JEV Assembly Pipeline: Generation of JEV consensus genomes

## Introduction

This is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline to assemble whole Japanese encephalitis virus (JEV) genomes from primer tiling data, such as those designed using [Primal Scheme](https://primalscheme.com/). The pipeline takes raw Illumina reads and trims for quality and adapters using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), aligns to a reference JEV genome with [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), then masks the tiling primers and generates a consensus genome with [iVar](https://andersen-lab.github.io/ivar/html/manualpage.html). 

The results are then summarised by calculating how much of the genome was obtained and coverage over the genome is plotted with [ggplot2](https://ggplot2.tidyverse.org/).

## Installation
This repository first needs to be cloned onto your system. For example:
```
git clone https://github.com/neavemj/JEV_assembly_pipe.git $HOME/software/JEV_assembly_pipe
```

The software required by the pipeline is listed in the environment.yaml file. If you use [Conda](https://conda.io/docs/install/quick.html), you can use the environment.yaml file to install everything you need:

```
conda env create -f $HOME/software/JEV_assembly_pipe/environment.yml
```

Once everything is installed, activate the environment:

```
conda activate JEV_assembly_pipe
```

