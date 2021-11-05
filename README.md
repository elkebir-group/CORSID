# CORSID

CORSID is a computational tool to simultaneously identifying TRS sites and gene locations in unannotated coronavirus genomes.
We also provide another tool CORSID-A that only identifies TRS sites given annotated genes.
Given an genome (optionally with their annotation), CORSID(-A) will find the TRS alignment and the core sequence.

The data and results can be found in the repo [CORSID-data](https://github.com/elkebir-group/CORSID-data). The visualized results of our tool applied to 468 coronavirus genomes can be found in [CORSID-viz](https://elkebir-group.github.io/CORSID-viz/).

![Figure](doc/overview.png)

## Contents

  1. [Pre-requisites](#pre-requisites)
  2. [Installation](#install)
      * [Using conda](#conda) (recommended)
      * [Using pip](#pip) (alternative)
  3. [Usage instructions](#usage)
      * [I/O formates](#io)
      * [Conda/compiled usage](#conda-usage)
      * [Docker usage](#docker-usage)

<a name="pre-requisites"></a>

## Pre-requisites
+ python3 (>=3.7)
+ [numpy](https://numpy.org/doc/)
+ [pysam](https://pysam.readthedocs.io/en/latest/)
+ [pandas](https://pandas.pydata.org/pandas-docs/stable/index.html)
+ [pytablewriter](https://pytablewriter.readthedocs.io/en/latest/)
+ (optional for simulation pipeline) [snakemake (>=5.2.0)](https://snakemake.readthedocs.io)

<a name="install"></a>

## Installation

<a name="conda"></a>

### Using conda (recommended)

1. Create a new conda environment named "corsid" and install dependencies:

   ```bash
   conda create -n corsid
   ```

2. Then activate the created environment: `conda activate corsid`.
3. Install the package into current environment "corsid":

    ```bash
    conda install -c bioconda corsid
    ```

<a name="pip"></a>

### Using pip (alternative)

We recommend installing in a virtual environment, as decribed in step 1 and 2 in the previous section.
Use `pip` to install the package:

```bash
pip install corsid
```

<a name="usage"></a>

## Usage instructions

<a name="io"></a>

### I/O formats

CORSID takes a **FASTA file** containing the complete genome as input. Optionally it also takes an **annotation file (GFF format)** to validate the identified genes.

CORSID-A takes a **FASTA file and an annotation file (GFF format)** as input. It will find candidate regions for each gene given the annotation file, and run CORSID-A on candidate regions.

The output is an JSON file containing sorted solutions and auxilary information. This file can be used as the input to the visualization webapp (**link**).
The program also outputs to the standard output, where it shows tables of solutions and visualization of TRS alignment.
