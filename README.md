# ONTmompS

**A tool for retrieving mompS2 alleles from closed long-read assemblies for Legionella pneumophila SBT assignment.**

* [Quick usage](#Quick-usage)
* [Installation](#Installation)
* [Full usage](#Full-usage)
* [Description](#Description)
* [Update database](#Update-database)

## Quick usage

```
python ONTmompS.py -a assembly.fasta
```

## Description
This tool was built as an _in silico_ approach to identify the Sequence Type (ST) for Legionella pneumophila from long-read assembles. Specifically, it distinguishes the mompS1 and mompS2 alleles from one another.
* A complete ST is reported if allele matches are found to all seven SBT genes in the database (e.g. ST560).
* If there are < 3 inexact matches, the nearest matching ST with the number of locus variants (LVs) is reported, e.g. ST560-1LV.
* For allele matches with <100% sequence identity, the nearest matching allele is noted with "*"
* For incomplete coverage of an allele, the nearest matching allele is noted with "?"
* For loci with no allele matches, sequence identity <90% or sequence coverage <80%, the allele number is reported as "-"

Sequence-based typing (SBT) of _Legionella pneumophila_ is a valuable tool in epidemiological studies and outbreak investigations of Legionnaires’ disease. In the _L. pneumophila_ SBT scheme, _mompS2_ is one of seven genes that determine the ST. The _Legionella_ genome typically contains two copies of _mompS_. When they are non-identical, it can be challenging to determine the mompS2 allele, and subsequently the ST, from Illumina sequencing, due to the short read-length. With long-read sequencing from Oxford Nanopore Technologies (ONT) Kit12 chemistry and R10.4.1 flow cells, together with Trycycler v0.5.3 and Medaka v1.7.2 for long-read assembly and polishing, we were able to identify the mompS2 allele and subsequently the _L. pneumophila_ SBT. 

**Cite**

If you use this tool, please cite: Krøvel AV and Hetland MAK et al. Long-read sequencing as a solution to the challenge of calling the mompS-allele for use in L. pneumophila SBT with short-reads. https://github.com/marithetland/ONTmompS


## Installation

Clone the repo and install dependencies. We recommend installing in a conda environment:

```
git clone https://github.com/marithetland/ONTmompS.git
mamba create -n ontmomps_env -c bioconda -c conda-forge pandas blast emboss parallel
```


## Full usage
Activate the conda environment: 

```
usage: ONTmompS.py [-h] [-v] -a ASSEMBLIES [ASSEMBLIES ...]
                   [-d DATABASE_FOLDER] [--store_mompS_alleles]
                   [--store_novel_alleles] [--store_all_alleles] [--verbose]
                   [--outfilename OUTFILENAME] [-outdir OUTDIR]

In silico SBT of Legionella pneumophila from long-read or hybrid assemblies

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

Input options (required):
  -a ASSEMBLIES [ASSEMBLIES ...], --assemblies ASSEMBLIES [ASSEMBLIES ...]
                        FASTA file(s) for assemblies (*.fasta)

Optional flags:
  -d DATABASE_FOLDER, --database_folder DATABASE_FOLDER
                        Provide a path to database location if different than
                        that provided by this tool.
  --store_mompS_alleles
                        Print mompS alleles to files named
                        {assembly}_{allele}.fna.
  --store_novel_alleles
                        Print novel alleles to files named
                        {assembly}_{allele}.fna.
  --store_all_alleles   Print all alleles (7 genes in SBT scheme + mompS1) to
                        files named {assembly}_{allele}.fna.
  --verbose             Keep intermediate files for debugging.

Output options:
  --outfilename OUTFILENAME
                        Output filename for STs. Default: ./LpST_ONTmompS.tsv
  -outdir OUTDIR        Output directory to store novel alleles in. Default is
                        current working directory
```

## Update database
TBD. Currently from Legsta. Should be updated with UKHSA db for more updated.