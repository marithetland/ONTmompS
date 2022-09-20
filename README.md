# ONTmomps

A script to retrieve mompS1 and mompS2 from a closed long-read assembly.

## Usage
Activate the conda environment: 

```
conda activate momps
```

Go to the folder containing the long-read assembly:

```
cd /path/to/assembly/fasta
```

Now you can run the script:
```
python ~/Scripts/ONT/legionella/ONTmomps/ONTmomps.py -r . -db ~/Scripts/ONT/legionella/ONTmomps -a consensus.fasta
```


## Installation
Install conda environment and install dependencies: 
```
mamba create -n momps python=3.9 
conda activate momps
mamba install -c anaconda pandas
mamba install -c bioconda blast
mamba install -c bioconda samtools
mamba install -c bioconda emboss
```

Clone the github repo and set up:
```
cd ~/Programs/
git clone https://github.com/markus-soma/ONTmomps.git
```
