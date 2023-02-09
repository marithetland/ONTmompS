# ONTmompS

**ONTmompS is a tool to perform _in silico_ Sequence Based Typing (SBT) of _Legionella pneumophila_ long-read/hybrid assemblies.**

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
This tool was built as an _in silico_ approach to identify the Sequence Type (ST) for Legionella pneumophila from long-read or hybrid assemblies. It first identifies the _mompS1_ and _mompS2_ alleles and then assigns ST.
* A complete ST is reported if allele matches are found to all seven SBT genes in the database (e.g. ST560).
* If there are < 3 inexact matches, the nearest matching ST with the number of locus variants (LVs) is reported, e.g. ST560-1LV.
* For allele matches with <100% sequence identity, the nearest matching allele is noted with "*"
* For incomplete coverage of an allele, the nearest matching allele is noted with "?"
* For loci with no allele matches, sequence identity <90% or sequence coverage <80%, the allele number is reported as "-"

Sequence-based typing (SBT) of _Legionella pneumophila_ is a valuable tool in epidemiological studies and outbreak investigations of Legionnaires’ disease. In the _L. pneumophila_ SBT scheme, _mompS2_ is one of seven genes that determine the ST. The _Legionella_ genome typically contains two copies of _mompS_ (designated _mompS1_ and _mompS2_). When they are non-identical, it can be challenging to determine the _mompS2_ allele, and subsequently the ST, from Illumina sequencing, due to the short read-length. With long-read sequencing from Oxford Nanopore Technologies (ONT) Kit12/Kit10 chemistry and R10.4.1/R9.4 flow cells, together with Trycycler v0.5.3 and Medaka v1.7.2 for long-read assembly and polishing, we were able to identify the _mompS2_ allele and subsequently the _L. pneumophila_ SBT, using this tool. 

**Cite**

If you use this tool, please cite: Krøvel AV and Hetland MAK et al. Long-read sequencing as a solution to the challenge of calling _mompS_ for _L. pneumophila_ SBT with short-reads. https://github.com/marithetland/ONTmompS


## Installation

Clone the repo and install dependencies. We recommend installing in a conda (mamba) environment:

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
  --ST_outfile OUTFILENAME
                        Output filename for STs. Default: ./LpST_ONTmompS.tsv
  --mompS_outfile OUTFILENAME
                        Output filename for mompS copy allele numbers. Default: ./mompS_alleles_ONTmompS.tsv
  -outdir OUTDIR        Output directory to store novel alleles in. Default is
                        current working directory
```

## Update database
The database in this repository is the same version as that in https://github.com/tseemann/legsta (https://github.com/tseemann/legsta/tree/master/db). Please contact UKHSA if you want to obtain a more recent database version. When you have your desired database, you can either specify the path in the command with flag `--db /path/to/db` or simply replace the db directory in this directory (that way you avoid having to specify db each time you run it).