# JEV Assembly Pipeline: Generation of JEV consensus genomes

## Introduction

This is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline to assemble whole Japanese encephalitis virus (JEV) genomes from primer tiling data, such as those designed using [Primal Scheme](https://primalscheme.com/). The pipeline takes raw Illumina reads and trims for quality and adapters using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), aligns to a reference JEV genome with [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), then masks the tiling primers and generates a consensus genome with [iVar](https://andersen-lab.github.io/ivar/html/manualpage.html). The genome assembly results are summarised and coverage over the genome is plotted with [ggplot2](https://ggplot2.tidyverse.org/).

## Installation
This repository first needs to be cloned onto your system. For example:
```
git clone https://github.com/neavemj/JEV_assembly_pipe.git $HOME/JEV_assembly_pipe
```

The software required by the pipeline is listed in the environment.yaml file. If you use [Conda](https://anaconda.com), you can use the environment.yaml file to install everything you need:

```
conda env create -f $HOME/JEV_assembly_pipe/environment.yaml
```

Once everything is installed, activate the environment:

```
conda activate JEV_assembly_pipe
```

## Running the pipeline

First create a directory where you want the results to go:

```
mkdir $HOME/JEV_data
```

Then copy the config.yaml file from the github repository into your new data directory: 

```
cp $HOME/JEV_assembly_pipe/config.yaml $HOME/JEV_data/
```

*Note: copying this file is very important and ensures clear separation between the software and the data. Each time you run the pipeline, copy a new config.yaml file to a new data directory* 

Now open the config.yaml file in your favorite text editor and change the installation location, the location of your sequence reads, and any other parameters.

Finally run the snakemake pipeline from your new directory, specifying the number of cores you wish to use. For example, to use 4 cores:

```
cd $HOME/JEV_data/

snakemake -s $HOME/JEV_assembly_pipe/assemble_JEV.snakemake -j 4
```

## Output Files

| Extension | Description |
| --------- | ----------- |
| 04_genomes/*.fa | these are the assembled consensus genomes from iVar |
| 04_genomes/*.qual.txt | the quality scores for each position in the genome |
| 05_report/*.png | coverage depth plots for each sample |
| 05_report/summary_stats.txt | summary of input reads, cleaned reads, and percent genome recovery |



### Example summary output

| Sample_name         | Input_reads | Cleaned_reads | Aligned_reads | Ref_coverage_(%) |
| ------------------- | ----------- | ------------- | ------------- | ---------------- |
| my_sample1              |    140994   |     130440    |      126639   |           98.46  |
| my_sample2              |    254976   |     229112    |      225946   |           39.42  |
| my_unpaired |    976526   |     960010    |      956992   |           68.52  |


### Example coverage figure

![alt text](https://github.com/neavemj/JEV_assembly_pipe/blob/main/config/JEV_sample1_coverage.png)

![alt text](https://github.com/neavemj/JEV_assembly_pipe/blob/main/config/JEV_sample2_coverage.png)


