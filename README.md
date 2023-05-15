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

Clone this repository and install the dependencies. We recommend installing in a conda (mamba) environment:

```
git clone https://github.com/marithetland/ONTmompS.git
mamba create -n ontmomps_env -c bioconda -c conda-forge pandas blast emboss biopython
```

## Description 
This tool was built as an _in silico_ approach to identify the Sequence Type (ST) of _Legionella pneumophila_ genomes from long-read or hybrid assemblies. It first identifies the _mompS1_ and _mompS2_ alleles and then assigns allele numbers and ST. We recommend using long-read or hybrid assemblies that have circular chromosomes when running this tool. It is not intended for short-read assemblies.

**Background**

Sequence-based typing (SBT) of _Legionella pneumophila_ is a valuable tool in epidemiological studies and outbreak investigations of Legionnaires’ disease. In the _L. pneumophila_ SBT scheme, _mompS2_ is one of seven genes that determine the ST. The _Legionella_ genome typically contains two copies of _mompS_ (designated _mompS1_ and _mompS2_). When they are non-identical, it can be challenging to determine the _mompS2_ allele, and subsequently the ST, from Illumina sequences, due to the short read-length. Using long-read sequencing from Oxford Nanopore Technologies (ONT) Kit12/Kit9 chemistry and R10.4/R9.4.1 flow cells, together with Trycycler v0.5.3 and Medaka v1.7.2 for long-read assembly and polishing, we were able to identify the _mompS2_ allele and subsequently the _L. pneumophila_ ST of 81/81 genomes when using this tool.

**Citation**

If you use ONTmompS, please cite the paper: [Krøvel AV, Hetland MAK, Bernhoff E, et al. Long-read sequencing for reliably calling the _mompS_ allele in _Legionella pneumophila_ sequence-based typing. Front. Cell. Infect. Microbiol (2023)](https://doi.org/10.3389/fcimb.2023.1176182)


## Output
With default settings, two output files are created: `LpST_ONTmompS.tsv` reports the ST and allele numbers of the 7 SBT loci. `mompS_alleles_ONTmompS.tsv` reports the mompS alleles. The ST and alleles will be annotated if there are mismatching or missing alleles:

* A complete ST is reported if allele matches are found to all seven SBT genes in the database (e.g. ST560).
* If there are < 3 inexact matches, the nearest matching ST with the number of locus variants (LVs) is reported, e.g. ST560-1LV.
* For allele matches with <100% sequence identity, the nearest matching allele is noted with "*"
* For incomplete coverage of an allele, the nearest matching allele is noted with "?"
* For loci with no allele matches, sequence identity <90% or sequence coverage <80%, the allele number is reported as "-"

It is possible to store the allele sequences to files, using the flags: `store_novel_alleles` to store only novel allele sequences, `store_mompS_alleles` to store the mompS allele sequences, or `store_all_alleles` to store all allele sequences. The files will by default be placed in a folder named `ONTmompS_allele_sequences`.



## Full usage

```
usage: ONTmomps.py [-h] [-v] -a ASSEMBLIES [ASSEMBLIES ...] [--db DB]
                   [--store_mompS_alleles] [--store_novel_alleles]
                   [--store_all_alleles] [--verbose] [-l LOG]
                   [--ST_outfile ST_OUTFILE] [--mompS_outfile MOMPS_OUTFILE]
                   [-o OUTDIR]

In silico SBT of Legionella pneumophila from long-read or hybrid assemblies

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

Input options (required):
  -a ASSEMBLIES [ASSEMBLIES ...], --assemblies ASSEMBLIES [ASSEMBLIES ...]
                        FASTA file(s) for assemblies (*.fasta)

Optional flags:
  --db DB               Provide a path to database location if different than
                        that provided by this tool.
  --store_mompS_alleles
                        Print mompS alleles to files named
                        {assembly}_{allele}.fna.
  --store_novel_alleles
                        Print novel alleles to files named
                        {assembly}_{allele}.fna.
  --store_all_alleles   Print all alleles (7 genes in SBT scheme + mompS1) to
                        files named {assembly}_{allele}.fna.
  --verbose             Log more details and keep intermediate files for
                        debugging.
  -l LOG, --log LOG     Write logging to specified file name instead of stdout

Output options:
  --ST_outfile ST_OUTFILE
                        Output filename for STs. Default: ./LpST_ONTmompS.tsv
  --mompS_outfile MOMPS_OUTFILE
                        Output filename for mompS copy allele numbers.
                        Default: ./mompS_alleles_ONTmompS.tsv
  -o OUTDIR, --outdir OUTDIR
                        Output directory to store novel alleles in. Default is
                        current working directory
```

## Update database
The database in this repository is the same version as that in https://github.com/tseemann/legsta/tree/master/db. Please contact the Legionella-SBT team at UKHSA if you want to obtain a more recent database version. When you have your desired database, you need to make sure the files follow the same format as those in the db for this repo. If the sequences are provided in csv format with 'sequence,number', you can run the commands below to make the database compatible with ONTmompS:

```
unzip sbt_schema_10_11_2022_12_59.zip ; cd sbt_schema_10_11_2022_12_59 
cat sbt.csv | sed 's/st/ST/g' | sed 's/,/\t/g' >> lpneumophila.txt ; rm sbt.csv 
cat neuAh.csv >> neuA.csv ; rm neuAh.csv 
for f in $(ls *.csv | sed 's/.csv//g') ; do paste <(cat ${f}.csv | cut -d"," -f2 | sed "s/^/>${f}_/g" ) <(cat ${f}.csv | cut -d"," -f1) |  sed 's/\t/\n/g' | grep -v "number\|sequence" >> ${f}.fna ; done ; rm *.csv
```

Once you have converted the database files, you can either 1) specify the path to your new db with the flag `--db /path/to/db` or 2) move your new db files to this repo's db directory (i.e. `mv *.fna lpneumophila.txt /path/to/ONTmompS/db/`). 
