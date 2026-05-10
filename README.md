# MIrROR: Microbial Identification using rRNA Operon Region
<img src="https://user-images.githubusercontent.com/31500750/100863170-817d0800-34d7-11eb-9928-8a19502762d0.png" height="180">

**Analysis tool for metataxonomics using 16S-ITS-23S rRNA operon region.**</br>
With the advancement of long-read sequencing technologies, the field of metataxonomics has entered a new phase. Analyzing the 16S-ITS-23S rRNA operon region (~4,300 bp) for bacterial community profiling provides more taxonomic information than analyzing only partial 16S rRNA gene sequences using short-read sequencing, allowing for species-level analysis. MIrROR provides a curated database and analysis tool for metaxonomics using 16S-ITS-23S rRNA operon region.

---

## 📢 Database Update: MIrROR release 02 (2026)

**MIrROR release 02** database has been published in *Scientific Data*:
Key highlights of release 02:
- Built from **1,690,470 genomes** (1,674,514 bacterial + 15,956 archaeal) from NCBI
- Final curated dataset: **476,579 sequences**, **249,907 genomes**, **29,051 species**
- **Archaeal genomes included** for the first time
- Taxonomy reclassified using **GTDB R220**

---

## Quick Start

```bash
# Install dependencies
conda install -c bioconda minimap2      # long-read mapper (required)
conda install -c bioconda krona rename  # Krona visualization (optional)
pip install pandas matplotlib           # stacked bar plots (optional)

# Clone repository
git clone https://github.com/seoldh/MIrROR.git
cd MIrROR
./MIrROR.py -h

# Download MIrROR release 02 database
mkdir DBDIR
wget -P DBDIR https://zenodo.org/records/17639192/files/MIrROR_r02.fa
wget -P DBDIR https://zenodo.org/records/17639192/files/MIrROR_r02.tsv
# After primer-aware curation (see below), build your own index:
# minimap2 -d DBDIR/MIrROR_r02.mmi MIrROR_r02_curated.fa

# Usage examples
MIrROR.py -K -d DBDIR -o result_BC01 ./sample_data/BC01.fastq
MIrROR.py -V -d DBDIR -o result_paf ./sample_data/list_PAF.txt
```

---

## 16S-ITS-23S rRNA operon Database

MIrROR uses a curated 16S-ITS-23S rRNA operon database.</br>
See [MIrROR Website](http://mirror2.egnome.co.kr/) for more information.

### Available releases
| Release | Sequences | Genomes | Species | Taxonomy | Download |
|---------|-----------|---------|---------|----------|----------|
| r01 | 97,781 | 43,653 | 9,485 | GTDB R89 | [legacy](http://mirror.egnome.co.kr/media/ToolsDatabase/2021_06/) |
| **r02** (recommended) | **476,579** | **249,907** | **29,051** | **GTDB R220** | [**Zenodo**](https://doi.org/10.5281/zenodo.17639192) |

### ⚠️ Recommended database preparation workflow
The Zenodo repository for release 02 contains three files:
- `MIrROR_r02.fa` — full-length FASTA sequences (start here for custom curation)
- `MIrROR_r02.mmi` — pre-built minimap2 index (built from full, untrimmed sequences)
- `MIrROR_r02.tsv` — taxonomy + operon copy number table (**required for MIrROR**)

**We recommend building your own minimap2 index from `MIrROR_r02.fa` after trimming sequences to match the primer set used in your experiment.**

A practical workflow:

```bash
# 1. Import FASTA into QIIME2
qiime tools import --type 'FeatureData[Sequence]' \
  --input-path MIrROR_r02.fa \
  --output-path MIrROR_r02.qza
 
# 2. Extract amplicon region with your primer sequences (example: primer set #5)
qiime feature-classifier extract-reads \
  --i-sequences MIrROR_r02.qza \
  --p-f-primer CCTACGGGNBGCWSCAG \
  --p-r-primer ACCRCCCCAGTHRAACT \
  --p-n-jobs 4 \
  --o-reads MIrROR_r02_trimmed.qza
 
# 3. Export and build minimap2 index
qiime tools export --input-path MIrROR_r02_trimmed.qza --output-path trimmed/
minimap2 -d DBDIR/MIrROR_r02.mmi trimmed/dna-sequences.fasta
```

## Dependencies
| Dependency | Purpose | Install |
|------------|---------|---------|
| [Minimap2](https://github.com/lh3/minimap2) | Read mapping (required) | `conda install -c bioconda minimap2` |
| [KronaTools](https://github.com/marbl/Krona) | Krona plots (`-K`) | `conda install -c bioconda krona rename` |
| [pandas](https://pandas.pydata.org/) | Stacked bar plots (`-S`) | `pip install pandas` |
| [matplotlib](https://matplotlib.org/) | Stacked bar plots (`-S`) | `pip install matplotlib` |


## Usage
```
 __  __ ___      ____   ___  ____          _   _
|  \/  |_ _|_ __|  _ \ / _ \|  _ \  __   _/ | / |
| |\/| || || '__| |_) | | | | |_) | \ \ / / | | |
| |  | || || |  |  _ <| |_| |  _ <   \ V /| |_| |
|_|  |_|___|_|  |_| \_\\___/|_| \_\   \_/ |_(_)_|

usage: MIrROR.py [options] (-d DBDIR) INPUTFILE

Input:
  INPUTFILE             FASTA/FASTQ/PAF file(s) or a sample list [required]

Main options:
  -d, --db_dir DBDIR    directory containing MIrROR database [required]
  -o, --output_dir OUTDIR
                        specify directory to output files (default: ./Result)
  -t, --threads INT     number of threads (default: 4)

Mapping options:
  -x, --preset STR      preset options to optimize alignment for different platforms. map-pb/map-hifi/map-ont -
                        CLR/HiFi/Nanopore (default: map-ont)
  -M, --minibatch NUM   number of query bases loaded to memory at once. K/M/G suffix accepted. (default: 500M)

Threshold options:
  -m, --residuematches INT
                        minimum number of residue matches (default: 2500)
  -b, --blocklength INT
                        minimum alignment block length (default: 3500)
  -n, --Normalization   normalize counts by 16S-23S rRNA operon copy number

Visualization options:
  -K, --Krona           create a Krona plot (requires KronaTools)
  -S, --Stackedplot     create stacked bar plots (requires pandas, matplotlib)
  -V, --visualized      perform all visualization tasks (same as -K -S)

Others:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
```

### Input

MIrROR accepts a single FASTA/FASTQ/PAF file or a sample list. Mapping is skipped when only PAF files are provided.</br>
**Sample list format** (tab-delimited; group columns are optional):
```
Sample         Status      Smoke
BC01.fastq     Healthy     smoker
BC02.fastq     Healthy     non-smoker
BC03.paf       Sick        smoker
BC04.paf       Sick        non-smoker
```

---
## Primer selection
 
Primer binding site sequences vary across species and lineages. As a result, no single primer set amplifies all taxa with equal efficiency, and **certain species in your sample may be missed depending on the primer pair used.** We strongly recommend evaluating candidate primer sets against the taxa most relevant to your study before finalizing your experimental design.
 
Use the [MIrROR Primer Checker](https://github.com/jisoll/MirrorPrimerChecker) for *in silico* coverage assessment against the release 02 database.

### Available primer sets
 
| Set | Forward primer | Reverse primer | Notes |
|-----|----------------|----------------|-------|
| **#4** (519F–2428R) | `CAGCMGCCGCGGTAA` | `CCRAMCTGTCTCACGACG` | Recommended for Bacteria + Archaea; best archaeal coverage (75.78% in silico) |
| **#5** (341F–2241R) | `CCTACGGGNBGCWSCAG` | `ACCRCCCCAGTHRAACT` | Recommended for Bacteria; highest bacterial coverage (98.81% in silico) |
| r01 original (27F–2241R) | `AGRGTTYGATYHTGGCTCAG` | `ACCRCCCCAGTHRAACT` | Used in MIrROR release 01 studies; bacteria only |

> Primer sets #4 and #5 are recommended based on *in silico* analysis in Lee et al. (2026).
> The r01 original set remains a valid choice for bacteria-focused studies.
> Regardless of primer choice, always verify coverage for the specific taxa of interest in your sample type before proceeding to wet lab work.

### Output

| Directory | File | Description |
|-----------|------|-------------|
| `./` | `RESULT.log` | run log |
| `./ReadMapping/` | `SAMPLE_minimap.paf` | minimap2 alignment |
| `./Classification/` | `SAMPLE.txt` | per-sample classification |
| `./FeatureTable/` | `OUTPUT_std.txt` | abundance table (standard) |
| `./FeatureTable/` | `OUTPUT_mpa.txt` | abundance table (MetaPhlAn-style) |
| `./FeatureTable/` | `OUTPUT_std_type2.txt` | abundance table (for Krona) |
| `./Visualization/` | `stacked_*.pdf` | stacked bar plots per taxonomic level |
| `./Visualization/krona/` | `SAMPLE.html` | interactive Krona chart |

#### Feature table example

```
#Name                                                                             BC01  BC02
d__Bacteria;p__Bacillota_A;c__Clostridia;...;s__Faecalibacterium_prausnitzii      25    24
```

---

## Citing MIrROR

If you use MIrROR in your research, **please cite both papers**:

**MIrROR tool and release 01 database:**
> Seol D, Lim JS, Sung S, Lee YH, Jung M, Cho S, Kwak W, Kim H.  
> **Microbial Identification Using rRNA Operon Region: Database and Tool for Metataxonomics with Long-Read Sequence.**  
> *Microbiology Spectrum*, 10(2): e02017-21 (2022). https://doi.org/10.1128/spectrum.02017-21
 
**MIrROR release 02 database:**
> Lee J, Hong J, Seol D, Lee W, Lee J, Kim G, Cho S, Kim H.  
> **MIrROR release 02: Expanded and refined 16S-ITS-23S rRNA operon dataset.**  
> *Scientific Data*, 13: 714 (2026). https://doi.org/10.1038/s41597-026-06729-y
