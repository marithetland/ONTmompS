# ONTmompS

**ONTmompS is a tool to perform _in silico_ Sequence Based Typing (SBT) of _Legionella pneumophila_ long-read/hybrid assemblies.**

* [Quick usage](#Quick-usage)
* [Installation](#Installation)
* [Full usage](#Full-usage)
* [Description](#Description)
* [Output](#Output)
* [Update database](#Update-database)

## Quick usage

```
python ONTmompS.py -a assembly.fasta
```

## Installation

Clone the repo and install dependencies. We recommend installing in a conda (mamba) environment:

```
git clone https://github.com/marithetland/ONTmompS.git
mamba create -n ontmomps_env -c bioconda -c conda-forge pandas blast emboss
```

## Description 
This tool was built as an _in silico_ approach to identify the Sequence Type (ST) for _Legionella pneumophila_ from long-read or hybrid assemblies. It first identifies the _mompS1_ and _mompS2_ alleles and then assigns allele numbers and ST. We recommend using long-read or hybrid assemblies that have circular chromosomes when running this tool. It will not work with short-read assemblies.

**Background**

Sequence-based typing (SBT) of _Legionella pneumophila_ is a valuable tool in epidemiological studies and outbreak investigations of Legionnaires’ disease. In the _L. pneumophila_ SBT scheme, _mompS2_ is one of seven genes that determine the ST. The _Legionella_ genome typically contains two copies of _mompS_ (designated _mompS1_ and _mompS2_). When they are non-identical, it can be challenging to determine the _mompS2_ allele, and subsequently the ST, from Illumina sequencing, due to the short read-length. Using long-read sequencing from Oxford Nanopore Technologies (ONT) Kit12/Kit10 chemistry and R10.4.1/R9.4 flow cells, together with Trycycler v0.5.3 and Medaka v1.7.2 for long-read assembly and polishing, we were able to identify the _mompS2_ allele and subsequently the _L. pneumophila_ SBT, when using this tool. 

**Cite**

If you use this tool, please cite: Krøvel AV and Hetland MAK et al. Long-read sequencing as a solution to the challenge of calling _mompS_ for _L. pneumophila_ SBT with short-reads. https://github.com/marithetland/ONTmompS

## Output
With default settings, two output files are created: `LpST_ONTmompS.tsv` reports the ST and SBT alleles. `mompS_alleles_ONTmompS.tsv` reports the mompS alleles. The ST and alleles will be annotated if there are mismatching or missing alleles:

* A complete ST is reported if allele matches are found to all seven SBT genes in the database (e.g. ST560).
* If there are < 3 inexact matches, the nearest matching ST with the number of locus variants (LVs) is reported, e.g. ST560-1LV.
* For allele matches with <100% sequence identity, the nearest matching allele is noted with "*"
* For incomplete coverage of an allele, the nearest matching allele is noted with "?"
* For loci with no allele matches, sequence identity <90% or sequence coverage <80%, the allele number is reported as "-"

It is possible to store the allele sequences to files, using flags: `store_novel_alleles` to store only novel allele sequences, `store_mompS_alleles` to store the mompS allele sequences, or `store_all_alleles` to store all allele sequences. The files will by default be placed in a folder named `ONTmompS_allele_sequences`.




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
The database in this repository is the same version as that in https://github.com/tseemann/legsta/tree/master/db. Please contact the Legionella-SBT team at UKHSA if you want to obtain a more recent database version. When you have your desired database, you need to make sure the files follow the same format as those in the db for this repo (see below). If the sequences are provided in csv format with 'sequence,number', you can run the commands below to make the database compatible with ONTmompS:

```
unzip sbt_schema_10_11_2022_12_59.zip ; cd sbt_schema_10_11_2022_12_59 
cat sbt.csv | sed 's/st/ST/g' | sed 's/,/\t/g' >> lpneumophila.txt ; rm sbt.csv 
cat neuAh.csv >> neuA.csv ; rm neuAh.csv 
for f in $(ls *.csv | sed 's/.csv//g') ; do paste <(cat ${f}.csv | cut -d"," -f2 | sed "s/^/>${f}_/g" ) <(cat ${f}.csv | cut -d"," -f1) |  sed 's/\t/\n/g' | grep -v "number\|sequence" >> ${f}.fna ; done ; rm *.csv
```

Once you have converted the database files, you can either 1) specify the path to your new db with the flag `--db /path/to/db` or 2) move your new db files to this repo's db directory (i.e. `mv *.fna lpneumophila.txt /path/to/ONTmompS/db/`). 
