# MIrROR: Microbial Identification using rRNA Operon Region
<img src="https://user-images.githubusercontent.com/31500750/100863170-817d0800-34d7-11eb-9928-8a19502762d0.png" height="180">

**Platform for metataxonomics with 16S-23S rRNA operon.**</br>
Recently, with the development of Nanopore sequencing technology, there have been studies evaluating taxonomic resolution using 16S-23S rRNA operon region as a taxonomic marker in bacteria. Analyzing 16S-23S rRNA operon has the advantage of having more information than analyzing only 16S rRNA sequence with short-read sequencing, enabling accurate analysis at the species level with compensating for high error rates of Nanopore read. MIrROR provides a curated database for rRNA operon and analysis tools.

## Quick Start

```
# For dependencies
conda install -c bioconda minimap2    # long-read mapper
conda install -c bioconda krona       # Krona plot
pip install pandas                    # Stacked bar plot
pip install matplotlib                # Stacked bar plot

# Getting Started
git clone https://github.com/seoldh/MIrROR.git
cd MIrROR
chmod +x lib/OTUsamples2krona.v0.2.2.sh
MIrROR.py -h

# Download 16S-23S rRNA operon database
mkdir DBDIR && cd DBDIR
wget http://mirror.egnome.co.kr/media/ToolsDatabase/2021_06/MIrROR_DB_r01.mmi
wget http://mirror.egnome.co.kr/media/ToolsDatabase/2021_06/MIrROR_DB_r01.tsv

# Usage examples
MIrROR.py -K -d DBDIR sample.fastq         # for single FASTQ file with Krona plot
MIrROR.py -d DBDIR fastq_list.txt          # for multi sample
MIrROR.py -V -d DBDIR paf_list.txt         # for alignment file with Krona, stacked plot
```

## rRNA operon Database

MIrROR uses a 16S-23S rRNA operon database curated over several steps.</br>
See [MIrROR Website](http://mirror.egnome.co.kr/) for more information.

## Dependencies
MIrROR uses long-read mapper, Minimap2. Please make sure that Minimap2 added to the path.
* [Minimap2](https://github.com/lh3/minimap2) (Tested under 2.17)

### For visualizations
For visualizations, MIrROR has the following dependencies:
* [KronaTools](https://github.com/marbl/Krona)
* Python packages
  + [pandas](https://pandas.pydata.org/)
  + [matplotlib](https://matplotlib.org/) 


## Usage
The program can be executed without installation.
```bash
git clone http://github.com/seoldh/MIrROR.git
cd MIrROR
python MIrROR.py -h
```

```
Usage: MIrROR.py [options] (-d DBDIR) INPUTFILE

MIrROR: a tool for metataxonomics with 16S-23S rRNA operon region
INPUTFILE: FASTA/FASTQ/PAF | SAMPLELIST

  Main options:
 
       -d DBDIR,      --db_dir=DBDIR
                        path to 16S-23S rRNA operon database
       -o OUTDIR,     --output_dir=OUT_DIR
                        write report to FOLDER (default : ./Result)
       -t INT,        --threads=INT
                        number of threads (default : 4)

  Mapping option:
    this option control for Mapping

       -M NUM,        --minibatch=NUM
                        the number of query bases loaded to the memory at once. K/M/G/k/m/g suffix is accepted. (default : 500M)

  Threshold options:
    these options control for Threshold

       -m INT,        --residuemaches=INT
                        number of residue matches in mapping (default : 2500)
       -b INT,        --blocklength=INT
                        alignment block length in mapping (default : 3500)
       -n,            --Normalization
                        normalize by 16S-23S rRNA operon copy number
                        
  Visualization options:
    these options control for Visualization

       -K,            --Krona
                        create a Krona plot. Requires KronaTools to be installed.
       -S,            --Stacked_plot
                        create a stacked plot. Requires Python packages to be installed.
       -V,            --visualized
                        perform all visualization tasks  (Same as: -K -S)
  Help:
       -h,            --help    Show this help message and exit
       -v,            --version Show program's version number and exit

```
### Input

MIrROR requires FASTA/FASTQ/PAF or sample list as input. In case of using PAF file as input, the mapping process is skipped and only abundance table is created.</br>
The sample list can contain one or more groups separated by tabs along with the header. Also, it can contain PAF and FASTA/FASTQ together.
```
more list.txt
Sample         Status      Smoke
BC01.fastq     Healthy     smoker
BC02.fastq     Healthy     non-smoker
BC03.paf       Sick        smoker
BC04.paf       Sick        non-smoker
```
The database path must be specified, and you can download the latest release from [here](http://mirror.egnome.co.kr/). In DBDIR, there are two requirements:
1. MIrROR_DB_r01.mmi (minimap2 index)
2. MIrROR_DB_r01.tsv (taxonomy, `accession<tab>taxonomy`)
```
cd DBDIR
ls
MIrROR_DB_r01.mmi MIrROR_DB_r01.tsv
```


### Output
All output files are located in the directory specified by OUT_DIR.

|Directory|Output file name|Description|
|---|------|---|
|./|||
||RESULT.log|log file|
|./mapping|||
||SAMPLE.paf|result of minimap2|
||SAMPLE.paf.log|log file for mapping|
|./classification||
||SAMPLE.txt|classification result for each sample|
|./OTUtable|||
||OUTPUT_std.txt|abundance table (standard version)|
||OUTPUT_std_type2.txt|abundance table (standard version)|
||OUTPUT_mpa.txt|abundance table (MetaPhlAn version)|
|./graph|||
||stacked_phylum.png|stacked bar plot for phylum level|
||···|···|
||stacked_species.png|stacked bar plot for species level|
|./graph/krona|||
||SAMPLE.html|krona plot for each sample|
||SAMPLE.tsv|percent information for krona|


#### Abundance table
MIrROR provides abundance table in two forms.

```
vi OUTPUT_std.txt
Name                                                                                                                  BC01     BC02        BC03    BC04
status                                                                                                                Healthy  Healthy     Sick    Sick
smoke                                                                                                                 smoker   non-smoker  smoker  non-smoker
p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Butyrivibrio;s__Butyrivibrio_fibrisolvens          51       19          49      13
p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Faecalibacterium;s__Faecalibacterium_prausnitzii   25       24          30      7

vi OUTPUT_mpa.txt
Name                          BC01     BC02        BC03    BC04
status                        Healthy  Healthy     Sick    Sick
smoke                         smoker   non-smoker  smoker  non-smoker
p__Firmicutes;                99       78          60      76
p__Firmicutes;c__Clostridia;  15       8           16      7
```


## Publications

Donghyeok Seol, et al. MIrROR: Microbial Identification using rRNA Operon Region. Submitted.
