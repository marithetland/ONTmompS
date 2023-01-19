# ONTmomps

**A tool for retrieving mompS2 alleles from closed long-read assemblies for Legionella pneumophila SBT assignment.**

* [Installation](#Installation)
* [Full usage](#Full-usage)
* [Cite](#Cite)


## Quick usage

```
conda activate momps 
python ONTmomps.py -a consensus.fasta
```

## Background



## Installation
Clone the repo:

```
git clone https://github.com/marithetland/ONTmomps.git
```

Install dependencies. We recommend installing in a conda environment like so:

```
mamba create -n ontmomps_env -c bioconda pandas blast samtools emboss
```


## Full usage
Activate the conda environment: 

```
usage: ONTmomps.py [-h] [-v] -r RUN_FOLDER -db DATABASE_FOLDER -a ASSEMBLY
                   [--threads THREADS]

Assign mompS2 allele from a long-read assembly for use in L. pneumophila SBT

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -r, --run_folder RUN_FOLDER
                        Input directory.
  -d, --database_folder DATABASE_FOLDER
                        Provide a path to database location.
  -a, --assembly ASSEMBLY
                        Provide a consensus assembly in fasta format.
  -t, --threads THREADS     Specify number of threads to use. Default: 4
```

## Cite
If you use this tool, please cite:
Soma MA, Hetland MAK, Bjorheim AS, et al. Assign mompS2 allele from a long-read assembly for use in L. pneumophila SBT. https://github.com/marithetland/ONTmomps


## SUS usage - fix this
Go to the folder containing the long-read assembly:

```
conda activate ontmomps_env

cd /path/to/assembly/fasta
```

Now you can run the script:
```
python ~/Scripts/ONT/legionella/ONTmomps/ONTmomps.py -r . -db ~/Scripts/ONT/legionella/ONTmomps -a consensus.fasta
```
