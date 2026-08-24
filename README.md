# miRDeep-P3
A comprehensive pipeline for plant miRNA discovery, annotation and downstream functional analysis from small RNA-seq data.

## Highlights

- **End-to-end miRNA discovery** — from raw sRNA-seq FASTQ to annotated
  mature/precursor miRNAs in a single workflow (identification →
  annotation → downstream analysis).

- **Multi-group comparison** — parallel processing of replicates
  (`-r`/`-p`) with consistent miRNA family naming across groups and runs
  (`--common`/`--consistency`).

- **Rich downstream analysis** — target prediction (psRNATarget-based),
  TFBS/promoter analysis (PlantTFDB), expression, enrichment and network
  modules, all in one tool.

- **Reproducible & publication-ready** — Docker/Conda packaged, multithreaded,
  and every figure exported in SVG/PDF/PNG formats.


## Table of Contents
- [Introduction](#introduction)
- [Workflow](#workflow)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input & Output](#input--output)
- [Configuration](#configuration)
- [Usage](#usage)
- [Examples](#examples)
- [Reproducibility](#reproducibility)
- [Testing](#testing)
- [Troubleshooting / FAQ](#troubleshooting--faq)
- [Citation](#citation)
- [License](#license)
- [Contributing](#contributing)
- [Contact](#contact)
- [Acknowledgements](#acknowledgements)

---

## Introduction

**mirdeep-p3** is an end-to-end pipeline for **plant miRNA discovery and
functional annotation** using **small RNA-seq data and comparative
genomics approaches**. It extends the classic miRDeep framework with
multi-sample comparison, family-level naming consistency, and a complete
suite of downstream analyses.

Typical use cases:
- Identify **novel and conserved miRNAs** from raw sRNA-seq data of any
  plant species, with stem-loop structure validation and mature/star
  annotation
- Perform **downstream analysis** in one place: target
  prediction (psRNATarget-based), promoter/TFBS analysis (PlantTFDB),
  differential expression, GO/KEGG enrichment and regulatory networks
- Generate **publication-ready figures** (SVG/PDF/PNG) and tables for
  manuscripts without extra scripting


If you use this software in academic work, please see [Citation](#citation).

## Workflow
> 建议放一张流程图（`docs/workflow.png`）或用 Mermaid。

```mermaid
flowchart LR
  A[FASTQ] --> B[QC/Trim]
  B --> C[Align/Map]
  C --> D[Quantification]
  D --> E[Downstream analysis]
  E --> F[Reports/Plots]
```

## Requirements

All dependencies are managed via conda. See
[`mirdp3_environment.yml`](mirdp3_environment.yml) for the full list.

**Core tools**: Python ≥3.9, R ≥4.2, Bowtie 1.3, Bowtie2, SAMtools,
BEDTools, SeqKit, ViennaRNA (RNAfold), Trim Galore

**Annotation/homology**: BLAST, HMMER, FASTA3 (ssearch), MMseqs2,
DIAMOND, eggNOG-mapper, MEME (FIMO), Prodigal, GFFread

**R plotting**: ggplot2, dplyr, tidyr, circlize, ComplexHeatmap,
svglite, patchwork, viridis, igraph/ggraph

**Bioconductor**: clusterProfiler, DESeq2, apeglm, enrichplot,
GOSemSim, DOSE, AnnotationForge, ggtree, Biostrings, GenomicRanges

## Installation
### Option A: Conda (recommended)
```bash
conda install mirdeep-p3
```

### Option B: From source
```bash
git clone https://github.com/YangXZ-lab/mirdeep-p3.git
cd mirdeep-p3
conda env create -f mirdp3_environment.yml -n mirdp3
conda activate mirdp3
chmod 755 mirdeep-p3
mirdeep-p3 -h
```

### Option C: Docker
```bash
# c1
## Extract and import
gunzip -c mirdeep-p3-3.1.4c-full.tar.gz | docker load
## or
docker load -i mirdeep-p3-3.1.4c-full.tar.gz

## Make sure the image is loaded.
docker images | grep mirdeep

# c2
## Pull from Docker Hub / GHCR
docker pull <yourname>/mirdeep-p3:3.1.4c-full
## Optional: Rename to short name
docker tag <yourname>/mirdeep-p3:3.1.4c-full mirdeep-p3:3.1.4c-full

# View help
docker run --rm mirdeep-p3:3.1.4c-full -h
# Enter the container interactively (debug/view output)
docker run --rm -it -v $(pwd):/data mirdeep-p3:3.1.4c-full /bin/bash
```

## Usage
### Show help
```bash
mirdeep-p3 --help
```
### miRNA identification from sRNA-Seq
```bash
##1 Single input
mirdeep-p3 identification \
  -i <sample_1.fq.gz> \
  -o <output_dir> \
  -g <genome.fasta> \
  -d <bowtie_index_prefix> \
  -t <threads>

##2 Multiple inputs
mirdeep-p3 identification \
  -i <sample_1.fq.gz,sample_2.fq.gz,sample_3.fq.gz> \
  -o <output_dir> \
  -g <genome.fasta> \
  -d <bowtie_index_prefix> \
  -t <threads> \
  -p <progress>
```
| Parameter | Description | Example |
|---|---|---|
| `-i, --input` | Raw sequencing data (FASTQ/FASTA/compressed files，multiple samples can be separated by commas). | `[e.g., sample_1.fq or sample_1.fq,sample_2.fq,sample_3.fq]` |
| `-o, --output` | Output dir. | `[e.g., output]` |
| `-g, --genome` | Reference genome FASTA file. | `[e.g., genome.fasta]` |
| `-d, --index` | Path prefix of a pre-built bowtie index (optional). If not specified, the index is built automatically under the output directory; provide a prefix here to reuse an existing index and skip building. | `[e.g., bowtie_index_prefix]` |
| `-t, --threads` | Number of threads used in mirdeep-p3 (optional, default: 1). | `[e.g., 14]` |
| `-p, --progress` | Number of samples/files processed in parallel (optional, default: 1). For example, if 3 input files are provided, `-p 3` will process all three simultaneously. | `[e.g., 3]` |

| Path | Description |
|---|---|
| `output/mirdp3-identification-<time>.pipe` | Log for this task. |
| `output/sample_1/<sample>_trimming_report.txt` | trim_galore reslut. |
| `output/sample_1/<sample>_identification.log` | Log for each step. |
| `output/sample_1/<sample>.total_reads` | Total reads count. |
| `output/sample_1/<sample>_filter_P_prediction` | Result of indentification step. |

### Advanced
- **Custom reference build**: pass your own genome with `-g`; if `-d`
  is omitted, the bowtie index is auto-built under the output/index directory.
- **Skip steps**: `--no-reads_clean` skip the reads clean step.

### miRNA annotation after identification
```bash
##1 Single input
mirdeep-p3 annotation \
  -i <output/sample_1> \
  -g <genome.fasta> \
  -d <bowtie_index_prefix> \
  --prefix <output prefix> \
  --prefix_miRNA <miRNA prefix> \
  --species <species> \
  -t <threads> \
  -o <output_dir>

##2 Multiple inputs
mirdeep-p3 annotation \
  -i <output/sample_1,output/sample_2,output/sample_3> \
  -g <genome.fasta> \
  -d <bowtie_index_prefix> \
  --prefix <output prefix> \
  --prefix_miRNA <miRNA prefix> \
  --species <species> \
  -t <threads> \
  -o <output_dir>
```
| Parameter | Description | Example |
|---|---|---|
| `-i, --input` | Output dir of miRNA identification step (multiple samples can be separated by commas). | `[e.g., sample_1 or sample_1,sample_2,sample_3]` |
| `-o, --output` | Output dir. | `[e.g., output]` |
| `-g, --genome` | Reference genome FASTA file. | `[e.g., genome.fasta]` |
| `-d, --index` | Path prefix of a pre-built bowtie index (optional). If not specified, the index is built automatically under the output directory; provide a prefix here to reuse an existing index and skip building. | `[e.g., bowtie_index_prefix]` |
| `-t, --threads` | Number of threads used in mirdeep-p3 (optional, default: 1). | `[e.g., 14]` |
| `--species` | Species name in the output (must be quoted). | `[e.g., "Arabidopsis thaliana"]` |
| `--prefix_miRNA` | miRNA prefix in the output (must be quoted). | `[e.g., "Ath"]` |
| `--prefix` | Prefix of output file (optional). | `[e.g., flower]` |

| Path | Description |
|---|---|
| `output/mirdp3-annotation-<time>.pipe` | Log for this task. |
| `output/<prefix>/<prefix>-basic-info` | Result of annotation step. |
| `output/<prefix>/<prefix>_annotation.log` | Log for each step. |
| `output/<prefix>/<prefix>-basic-info-cluster` | miRNA cluster result. |
| `output/<prefix>/<prefix>-mature.count` | Read count matrix. |
| `output/<prefix>/<prefix>-mature.exp` | Expression matrix. |

### Advanced
- **Custom reference build**: pass your own genome with `-g`; if `-d`
  is omitted, the bowtie index is auto-built under the output/index directory.
- **Skip steps**: `--no-reads_clean` skip the reads clean step.
- **Multi group**: `-r/--replicate` used for multigroup; if -r 3,3: the first three inputs form group A, the last three form group B. Each group is processed and reported independently.
- **prefix**: one prefix per group, comma-separated. The number of prefixes must match the number of groups (2 groups → 2 prefixes).
- **--common**: when more than one group is present, `--common` unifies the naming of non-conserved miRNA families across groups, ensuring family names are consistent between groups for downstream comparison.
- **--consistency**: path to a reference basic-info file produced by a previous `mirdeep-p3 annotation` run (same format). When provided, the current annotation is forced to keep miRNA naming consistent with that reference — i.e. the same miRNA will receive the same name in both datasets. This is useful when re-annotating the same species or when merging annotations across groups.

### miRNA annotation result statistics
```bash
mirdeep-p3 analysis Stat \
  -i <basic-info> \
  -o <output_dir> \
  --rnaplot
```
| Parameter | Description | Example |
|---|---|---|
| `-i, --input` | Result of annotation step (basic info file). | `[e.g., sample-basic-info]` |
| `-o, --output` | Output dir. | `[e.g., output]` |
| `--rnaplot` | Plotting miRNA Stem-loop structure (optional). | `[e.g., --rnaplot]` |

| Path | Description |
|---|---|
| `output/base_dist.svg` | Base composition distribution. |
| `output/first_base_dist.svg` | First base composition distribution. |
| `output/length_dist.svg` | miRNA length distribution. |
| `output/miRNA_family.svg` | miRNA family distribution. |
| `output/Ath-MIR156a_ss.svg` | Stem-loop structure of Ath-MIR156a. |
> **Note**: All figures except miRNA Stem-loop structure are also exported in **PDF and PNG** formats
> (e.g. `base_dist.pdf`, `base_dist.png`) — use whichever suits your.

### miRNA target identification
```bash
mirdeep-p3 analysis Target_finder \
  -i <mature_miRNA.fasta> \
  -c <cds.fasta> \
  -o <output_dir> \
  -t <threads>
```
| Parameter | Description | Example |
|---|---|---|
| `-i, --input` | miRNA mature sequences (fasta). | `[e.g., sample-mature.fasta]` |
| `-o, --output` | Output dir. | `[e.g., output]` |
| `-t, --threads` | Number of threads used in mirdeep-p3 (optional, default: 1). | `[e.g., 20]` |
| `-b, --basic` | Result of annotation step (basic info file, conflict with `-i`). | `[e.g., sample-basic-info]` |
| `-e, --evalue` | E-value threshold for target prediction (optional, default: 2.5). | `[e.g., 30]` |
| `--GUs` | Allowed G:U mismatches in the complementary region (optional, default: 1). | `[e.g., 0.5]` |

| Path | Description |
|---|---|
| `output/target_finder.tsv` | miRNA target identification results. |
> **Notes**:
> 1. The scoring scheme follows **psRNATarget** (https://www.zhaolab.org/psRNATarget/), implemented on top of the scripts from [jtremblay/mirnatarget](https://github.com/jtremblay/mirnatarget).
> 2. Prediction stringency can be tuned with `-e/--evalue` (E-value threshold) and `--GUs` (allowed G:U mismatches).
> 3. When a basic-info file from a previous `annotation` run is available, `-b/--basic` is preferred over `-i/--input`, since it carries family and strand information that improves target prediction.

### miRNA promoter analysis
```bash
mirdeep-p3 analysis TFBS \
  -i <basic-info> \
  --fai <genome.fasta.fai> \
  -g <genome.fasta> \
  -o <output_dir> \
  -p
```
| Parameter | Description | Example |
|---|---|---|
| `-i, --input` | Result of annotation step. | `[e.g., sample-basic-info]` |
| `-o, --output` | Output dir. | `[e.g., output]` |
| `--fai` | fai index of reference genome. | `[e.g., genome.fasta.fai]` |
| `-g, --genome` | Reference genome FASTA file. | `[e.g., genome.fasta]` |
| `-p, --picture` | Generate TFBS report picture. | `[e.g., -p]` |
| `-s, --species` | Species name, quoted (optional, default: Arabidopsis_thaliana). | `[e.g., "Arabidopsis_lyrata"]` |
| `--list` | List available species. | `[e.g., --list]` |
| `-b, --bed` | Bed format input (conflict with `-i`). | `[e.g., sample-basic-info.bed]` |
| `-u, --upstream` | Upstream length to extract (default: 2000). | `[e.g., 3000]` |
| `-e, --evalue` | E-value threshold for FIMO (default: 1e-6). | `[e.g., 1e-6]` |

| Path | Description |
|---|---|
| `output/TFBS_count.svg` | TFBS count per miRNA. |
| `output/TFBS_distribution.svg` | TFBS distribution per miRNA. |
| `output/tfbs.tsv` | Result of miRNA promoter analysis. |
| `output/TF_family_count.svg` | TFBS Count per TF Family. |
| `output/TF_miRNA_network.svg` | TF–miRNA regulatory network. |
| `output/tfbs.log` | Log for each step. |
> **Note**:
> 1. All figures are also exported in **PDF and PNG** formats
>    (e.g. `TFBS_count.pdf`, `TFBS_count.png`) — use whichever suits your needs.
> 2. TFBS prediction follows the workflow of **PlantTFDB** (https://planttfdb.gao-lab.org/); all available species data are also sourced from PlantTFDB.

### miRNA differential expression analysis
```bash
mirdeep-p3 analysis Differential_expression \
  -c <mature.count> \
  -r <mature.exp> \
  --case1 <Columns corresponding to control group> \
  --case2 <Columns corresponding to treatment group> \
  -o <output_dir> \
  --case1name <control> \
  --case2name <treatment> \
  --DEOnly
```
| Parameter | Description | Example |
|---|---|---|
| `-c, --count` | Read count matrix. | `[e.g., sample-mature.count]` |
| `-r, --rpm` | Expression matrix. | `[e.g., sample-mature.exp]` |
| `--case1` | Columns corresponding to samples in the control group (separated by commas). | `[e.g., 3,4,5]` |
| `--case2` | Columns corresponding to samples in the treatment group (separated by commas). | `[e.g., 6,7,8]` |
| `-o, --output` | Output dir. | `[e.g., output]` |
| `--case1name` | Name of control group (optional, default: case1, must be quoted). | `[e.g., "flower"]` |
| `--case2name` | Name of treeatment group (optional, default: case2, must be quoted). | `[e.g., "root"]` |
| `--DEOnly` | Display only the expression patterns of differentially expressed miRNAs (optional). | `[e.g., --DEOnly]` |
| `--miRNA` | Comma-separated list of miRNA names to display (optional). | `[e.g., Ath-miR156a,Ath-miR157a]` |
| `-f, --file` | File containing miRNA names, one per line, equivalent to `--miRNA` (optional, conflict with `--miRNA`). | `[e.g., miRNA_list.txt]` |
| `--min-expr` | Minimum expression in at least one sample to retain miRNA (default: 5.0). | `[e.g., 10]` |

| Path | Description |
|---|---|
| `output/final_miRNA_expression.txt` | miRNA differential expression analysis results. |
| `output/heatmap.svg` | miRNA expression heatmap in all select samples. |
| `output/PCA_scatter.svg` | PCA results. |
| `output/sample_correlation_heatmap.svg` | Sample correlation heatmap. |
| `output/volcano_flower_root.svg` | Volcano plot of differential expressed miRNA. |
| `output/miRNA_Sta-MIR160a_expression.svg` | The expression patterns of a differentially expressed miRNA in different groups. |
> **Note**:
> 1. All figures are also exported in **PDF and PNG** formats
>    (e.g. `PCA_scatter.pdf`, `PCA_scatter.png`) — use whichever suits your needs.
> 2. `--case1` and `--case2` must have the same number of columns (replicates). `--case1` and `--case2` cannot have any intersection.

### miRNA functional analysis
#### download emapperdb form **eggnog** (http://eggnog5.embl.de/download/emapperdb-5.0.2/) (Choose the latest version)
```bash
##1 download and decompress (recommend in /PATH_to_mirdeep-p3/data/eggnog_data_dir/)
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/mmseqs.tar.gz
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/pfam.tar.gz
gunzip eggnog.db.gz
gunzip eggnog_proteins.dmnd.gz
tar xzvf eggnog.taxa.tar.gz
tar xzvf mmseqs.tar.gz
tar xzvf pfam.tar.gz
##2 set EGGNOG_DATA_DIR
conda env config vars set EGGNOG_DATA_DIR=/PATH_to_mirdeep-p3/data/eggnog_data_dir
##3.1 run miRNA functional analysis using protein file
mirdeep-p3 analysis Functional_analysis \
  -p <protein.fasta> \
  -g <target_gene_list> \
  -t <threads> \
  -o <output_dir>
##3.2 run miRNA functional analysis using exsit orgdb
###3.2.1 run emapper
emapper.py --cpu 20 \
  -m diamond \
  --override \
  --dbmem \
  -d euk \
  --tax_scope Viridiplantae \
  -i <protein.fasta> \
  -o <output_dir>
sed '/^##/d' output_dir/*.emapper.annotations | sed 's/#//g'| awk -vFS="\t" -vOFS="\t" '{print $1,$9,$10,$12}' > output_dir/GO.emapper.annotations
###3.2.2 build orgdb with emapper annotations
Rscript bin/build_orgdb.R \
  -i output_dir/GO.emapper.annotations \
  --kojson data/ko00001.json \
  -o <OrgDB_output_dir>
###3.2.3 run miRNA functional analysis
mirdeep-p3.py analysis Functional_analysis \
  --orgdb <OrgDB_output_dir>/org.Morg.eg.db \
  -f <OrgDB_output_dir> \
  -t <threads> \
  --gene <target_gene_list> \
  -o <output_dir>
```
| Parameter | Description | Example |
|---|---|---|
| `-p, --protein` | Protein sequence FASTA file. | `[e.g., sample-protein.fasta]` |
| `-g, --gene` | Target gene list (a gene per line). | `[e.g., target_gene.list]` |
| `-t, --threads` | Number of threads used in mirdeep-p3 (optional, default: 1). | `[e.g., 20]` |
| `-o, --output` | Output dir. | `[e.g., output]` |
| `--orgdb` | Path to OrgDb directory (org.Morg.eg.db) if already built. | `OrgDB_output_dir/org.Morg.eg.db` |
| `-f, --file` | Path to directory containing pathway2gene and pathway2name files. | `[e.g., OrgDB_output_dir]` |
| `--EGGNOG_DATA_DIR` | Path to eggNOG data directory (optional, default: data/eggnog_data_dir). | `[e.g., /PATH/eggnog_data_dir/]` |
| `--kojson` | Path to ko00001.json (otional, default: data/ko00001.json) | `[e.g., ko00001.json]` |
| `--target` | miRNA-target file (at least two columns: miRNA, target, conflict with `-g, --gene`). | `[e.g., output/target_finder.tsv]` |
| `--chord` | Generate chord diagram (optional, only with --target). | `[e.g., --chord]` |

| Path | Description |
|---|---|
| `output/Go.eggnog.emapper.annotations` | Emapper results used for functional analysis. |
| `output/org.Morg.eg.db/` | OrgDb directory. |
| `output/ego_list` | miRNA functional analysis results. |
| `output/GO_bar.svg` | GO bar plot. |
| `output/go.svg` | GO dot plot. |
| `output/kegg.svg` | KEGG dot plot. |
> **Note**:
> 1. All figures are also exported in **PDF and PNG** formats
>    (e.g. `GO_bar.pdf`, `GO_bar.png`) — use whichever suits your needs.
> 2. The `--chord` parameter must be used together with the `--target` parameter.
> 3. Please make sure the output folder is empty before running to avoid OrgDB build errors.

### One step analysis
```bash
mirdeep-p3 analysis Onestep \
  -b <basic-info> \
  -c <mature.count> \
  -r <mature.exp> \
  --fai <genome.fasta.fai> \
  -g <genome.fasta> \
  -t <cds.fasta> \
  -p <protein.fasta> \
  -o <output_dir> \
  --case1 <Columns corresponding to control group> \
  --case2 <Columns corresponding to treatment group> \
  --case1name <control> \
  --case2name <treatment> \
  --threads <threads> \
  --rnaplot \
  --tfbsplot \
  --DEOnly \
  --chord
```
> **Note**:
> 1. All figures are also exported in **PDF and PNG** formats
>    (e.g. `GO_bar.pdf`, `GO_bar.png`) — use whichever suits your needs.
> 2. All available parameters are the same as those used in the preceding steps.
## Examples
See:
- `examples/` for minimal runnable examples
- `docs/` for extended tutorials

## Reproducibility
This repository provides:
- Exact dependency locking via: [Conda env / Docker image / lockfile]
- Version tracking in `results/version_info.txt`
- Deterministic parameters recorded in `results/run_config.yaml`

Recommended practice:
- Run with fixed versions
- Commit `environment.yml` and reference hashes
- Archive outputs with `results/` + `logs/`

## Testing
```bash
pytest -q
# or
bash tests/run_smoke_test.sh
```
CI (GitHub Actions) should validate:
- Lint/format
- Unit tests
- Small end-to-end toy dataset

## Construct new core dataset and evaluation
```bash
python scripts/mirdp3_core_build.py \
  -p data/PmiREN-20260810-isoform.fa \ ##old core dataset
  -i example/PmiREN2.0_basic_info_MIR_unqiue.fasta \ ##non-redundant novel miRNAs dataset
  -s 12178 \ ##family number
  -t 14 \ ##thraeds
  -o test/PmiREN-v2 ##output

#evaluate the new core dataset
##1. build index
makeblastdb -in test/PmiREN-v2/isoform-in-v2.fa \
  -dbtype nucl \
  -out test/index/isoform-in-v2
##2. self-alignment using blastn
blastn -task blastn-short \
  -query test/PmiREN-v2/isoform-in-v2.fa \
  -db test/index/isoform-in-v2 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qseq sseq" \
  -out test/PmiREN-v2/isoform-in-v2.aln \
  -num_threads 14
##3. filtering and scoring
awk '$4>=13{ print $0 }' test/PmiREN-v2/isoform-in-v2.aln > test/PmiREN-v2/isoform-in-v2.aln.filter

python bin/bowtie_scoring_blastn.py \
  -i test/PmiREN-v2/isoform-in-v2.aln.filter \
  -f test/PmiREN-v2/isoform-in-v2.fa \
  -o test/PmiREN-v2/isoform-in-v2.aln.filter.score
##4. Deduplication
python script/dedup_score_pair.py \
  -i test/PmiREN-v2/isoform-in-v2.aln.filter.score \
  -o test/PmiREN-v2/isoform-in-v2.aln.filter.score.best
##5. Statistics and Visualization
Rscript scripts/classify_scores.R \
  -i test/PmiREN-v2/isoform-in-v2.aln.filter.score.best \
  -o test/PmiREN-in-v2-classify
###category1 only aligned to itself
###category2.1 only aligned to members of the same family
###category2.2 only aligned to members of the different family
###category2.3 aligned to members of the same and the different family

##Perform sequence deduplication before processing##
python bin/dedup_family.py \
  -i example/PmiREN2.0_basic_info.fasta \
  -o example/PmiREN2.0_basic_info_MIR.fasta
python bin/dedup_seq.py \
  -i example/PmiREN2.0_basic_info_MIR.fasta \
  --uni example/PmiREN2.0_basic_info_MIR_unique.fasta \
  --dup example/PmiREN2.0_basic_info_MIR_dup.fasta

##Remove the same seq(optional)
python bin/remove_matched_seq.py \
  -i example/PmiREN2.0_basic_info.fasta \
  -r data/PmiREN-20260810-isoform.fa \
  -o example/PmiREN2.0.fasta

##Check the final miRNA family number
cat data/PmiREN-20260810-isoform.fa | grep ">" | sort -V | tail
```
## Reannotation of a miRNA dataset
```bash
python bin/anno_miRNA.py \
  -i example/PmiREN2.0_basic_info_MIR_unqiue.fasta \ ##miRNA dataset
  -p data/PmiREN-20260810-isoform.fa \ ##pmiREN core dataset
  -o test/anno_miRNA \ ##output
  --threads 14 \
  --type MIRN \  ##prefix of new miRNA family
  --prefix "Ath"  ##prefix of miRNA(optional)
###anno.fasta, renamed miRNA
###anno.map, rename map
```

## Troubleshooting / FAQ
### Q1: [Common error]
**A:** [Fix / command / explanation]

### Q2: Low mapping rate?
**A:** Check:
- adapter trimming
- reference build mismatch
- read length filtering

## Citation
If you use **[Project Name]**, please cite:

- **Software:** [Authors]. *[Project Name]* (Version X.Y.Z). GitHub, Year. URL: <repo-url>
- **Paper (if any):** [Authors]. *Title*. Journal Year. DOI

You can also use the `CITATION.cff` file.

## License
This project is licensed under the **[MIT/BSD-3/GPL-3.0/Apache-2.0]** License. See `LICENSE`.

## Contributing
PRs and issues are welcome.
- Please read `CONTRIBUTING.md`
- Run tests before submitting
- Follow code style: `pre-commit run -a`

## Contact
- Maintainer: [Name] ([email])
- Lab/Org: [Lab name]
- Issues: please open a GitHub Issue

## Acknowledgements
- [Funding]
- [Upstream tools]
- [Contributors]
