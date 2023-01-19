# ONTmomps

**A tool for retrieving mompS2 alleles from closed long-read assemblies for Legionella pneumophila SBT assignment.**

* [Installation](#Installation)
* [Full usage](#Full-usage)
* [Cite](#Cite)


## Quick usage

```
python ONTmomps.py -a assembly.fasta
```

## Background 
TBD.


## Installation

Clone the repo and install dependencies. We recommend installing in a conda environment:

```
git clone https://github.com/marithetland/ONTmomps.git
mamba create -n ontmomps_env -c bioconda pandas blast samtools emboss
```


## Full usage
Activate the conda environment: 

```
usage: ONTmomps.py [-h] [-v] -a ASSEMBLIES [-o OUTDIR] [-d DATABASE_FOLDER]
                   [-t THREADS] [-k {on,off}]

Assign mompS from long-read assemblies

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

Input options (required):
  -a ASSEMBLIES, --assemblies ASSEMBLIES
                        Provide assemblies in fasta format.

Output options:
  -o OUTDIR, --outdir OUTDIR
                        Output directory for all output files. Default:
                        ./ONTmomps_output/

Optional flags:
  -d DATABASE_FOLDER, --database_folder DATABASE_FOLDER
                        Provide a path to database location if different than
                        that provided by this tool.
  -t THREADS, --threads THREADS
                        Specify number of threads to use. Default: 4
  -k {on,off}, --keep_intermediate_files {on,off}
                        Keep intermediate files. Default=off.
```

## Cite
If you use this tool, please cite:
Krøvel AV et al. Long-read sequencing as a solution to the challenge of calling the mompS-allele for use in L. pneumophila SBT with short-reads. https://github.com/marithetland/ONTmomps
