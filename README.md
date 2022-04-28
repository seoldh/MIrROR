# MIrROR: Microbial Identification using rRNA Operon Region
<img src="https://user-images.githubusercontent.com/31500750/100863170-817d0800-34d7-11eb-9928-8a19502762d0.png" height="180">

**Analysis tool for metataxonomics using 16S-23S rRNA operon region.**</br>
With the advancement of long-read sequencing technologies, the field of metataxonomics has entered a new phase. Analyzing the 16S-23S rRNA operon region (~4,300 bp) for bacterial community profiling provides more taxonomic information than analyzing only partial 16S rRNA gene sequences using short-read sequencing, allowing for species-level analysis. MIrROR provides a curated database and analysis tool for metaxonomics using 16S-23S rRNA operon region.

## Quick Start

```
# For dependencies

conda install -c bioconda minimap2      # long-read mapper
conda install -c bioconda krona rename  # Krona plot
pip install pandas matplotlib           # Stacked bar plot

# Getting Started
git clone https://github.com/seoldh/MIrROR.git
cd MIrROR
./MIrROR.py -h

# Download 16S-23S rRNA operon database
mkdir DBDIR
wget -P DBDIR http://mirror.egnome.co.kr/media/ToolsDatabase/2021_06/MIrROR_DB_r01.mmi
wget -P DBDIR http://mirror.egnome.co.kr/media/ToolsDatabase/2021_06/MIrROR_DB_r01.tsv

# Usage examples
./MIrROR.py -K -d DBDIR -o result_BC01 ./sample_data/BC01.fastq           # for single FASTQ file with Krona plot
./MIrROR.py -V -d DBDIR -o result_paf ./sample_data/list_PAF.txt          # for alignment files with Krona and stacked bar plot
```

## 16S-23S rRNA operon Database

MIrROR uses a 16S-23S rRNA operon database curated over several steps.</br>
See [MIrROR Website](http://mirror.egnome.co.kr/) for more information.

## Dependencies
MIrROR uses long-read mapper, Minimap2. Please make sure that Minimap2 added to the path.
* [Minimap2](https://github.com/lh3/minimap2)

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
./MIrROR.py -h
```

```
Usage: MIrROR.py [options] (-d DBDIR) INPUTFILE

MIrROR: a tool for metataxonomics with 16S-23S rRNA operon region.
INPUTFILE: FASTA/FASTQ/PAF | SAMPLELIST

  Main options:
 
       -d DBDIR,      --db_dir=DBDIR
                        directory containing MIrROR database [required]
       -o OUTDIR,     --output_dir=OUT_DIR
                        specify directory to output files (default : ./Result)
       -t INT,        --threads=INT
                        number of threads (default : 4)

  Mapping option:
    this option control for Mapping

       -M NUM,        --minibatch=NUM
                        number of query bases loaded to the memory at once. K/M/G/k/m/g suffix is accepted. (default : 500M)

  Threshold options:
    these options control for Threshold

       -m INT,        --residuemaches=INT
                        number of residue matches (default : 2500)
       -b INT,        --blocklength=INT
                        alignment block length (default : 3500)
       -n,            --Normalization
                        normalize by 16S-23S rRNA operon copy number
                        
  Visualization options:
    these options control for Visualization

       -K,            --Krona
                        create a Krona plot. Requires KronaTools to be installed.
       -S,            --Stacked_plot
                        create a stacked bar plot. Requires Python packages to be installed.
       -V,            --visualized
                        perform all visualization tasks  (Same as: -K -S)
  Help:
       -h,            --help    show this help message and exit
       -v,            --version show version number and exit

```
### Input

MIrROR requires a FASTA/FASTQ/PAF file or a sample list as input. The mapping process is skipped if the user only provides PAF file(s).</br>
The sample list can contain one or more group information (tab-delimited). Also, it can contain PAF and FASTA/FASTQ files together.
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
|./visualization|||
||stacked_phylum.png|stacked bar plot for phylum level|
||···|···|
||stacked_species.png|stacked bar plot for species level|
|./visualization/krona|||
||SAMPLE.html|krona plot for each sample|
||SAMPLE.tsv|percent information for krona|


#### Abundance table
MIrROR provides abundance tables in two forms.

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


## Citing MIrROR

Please cite this paper when using MIrROR for your publications.
> **Microbial Identification Using rRNA Operon Region: Database and Tool for Metataxonomics with Long-Read Sequence.** </br>
> *Microbiology Spectrum*, 10(2):e02017-21 (2022). https://doi.org/10.1128/spectrum.02017-21 </br>
> Donghyeok Seol, Jin Soo Lim, Samsun Sung, Young Ho Lee, Misun Jung, Seoae Cho, Woori Kwak, and Heebal Kim.
