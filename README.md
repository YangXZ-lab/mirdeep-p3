# [mirdeep-p3]
[一句话说明：这是做什么的。例：A reproducible pipeline for miRNA discovery and quantification from small RNA-seq.]

## Highlights
- ✅ Feature 1 (e.g., end-to-end pipeline from FASTQ to results)
- ✅ Feature 2 (e.g., supports hg38/mm10 and custom references)
- ✅ Feature 3 (e.g., Docker/Conda ready, fully reproducible)
- ✅ Feature 4 (e.g., produces publication-ready plots/tables)

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
**[Project Name]** is a [tool/pipeline/package] for **[biological goal]** using **[method/algorithm]**.

Typical use cases:
- [Use case 1]
- [Use case 2]
- [Use case 3]

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
### Supported Platforms
- Linux (recommended)
- macOS (partial/optional)
- Windows (WSL recommended)

### Dependencies
- [Python >= 3.10] / [R >= 4.3] / [Nextflow/Snakemake] / ...
- [bowtie/bwa/star/samtools/bedtools/...] (按你的项目写)
- Optional: Docker / Singularity

> 建议把完整依赖放到 `environment.yml`、`requirements.txt` 或 `pyproject.toml`，README 只列关键项。

## Installation
### Option A: Conda (recommended)
```bash
conda env create -f environment.yml
conda activate [env-name]
```

### Option B: Pip
```bash
pip install -r requirements.txt
# or
pip install .
```

### Option C: Docker
```bash
docker build -t [image-name] .
docker run --rm -it [image-name] --help
```

## Quick Start
### 1) Prepare input
```bash
mkdir -p data results
# put your FASTQ/FASTA/metadata into data/
```

### 2) Run
```bash
[main_command] \
  --input data/[input] \
  --outdir results \
  --threads 8
```

### 3) Check outputs
```bash
ls -lah results
```

## Input & Output
### Input
| Name | Type | Description |
|---|---|---|
| `--input` | file/dir | [e.g., FASTQ directory] |
| `--reference` | file | [e.g., genome fasta] |
| `--metadata` | table | [e.g., sample sheet CSV/TSV] |

**Example sample sheet (`samples.tsv`)**
```tsv
sample_id	condition	fastq_1	fastq_2
S1	CTRL	data/S1_R1.fastq.gz	data/S1_R2.fastq.gz
S2	TREAT	data/S2_R1.fastq.gz	data/S2_R2.fastq.gz
```

### Output
| Path | Description |
|---|---|
| `results/summary.tsv` | Main summary table |
| `results/logs/` | Logs for each step |
| `results/qc/` | QC reports (FastQC/MultiQC) |
| `results/plots/` | Publication-ready figures |
| `results/version_info.txt` | Tool + reference versions |

## Configuration
- Default config: `config/default.yaml`
- Example config: `config/example.yaml`

Key parameters:
- `threads`: number of threads
- `genome_build`: hg38/mm10/custom
- `adapter_sequence`: [small RNA adapter]
- `min_length`: [e.g., 18]

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
| `-i,--input` | Raw sequencing data (FASTQ/FASTA/compressed files，multiple samples can be separated by commas). | `[e.g., sample_1.fq or sample_1.fq,sample_2.fq,sample_3.fq]` |
| `-o,--output` | Output dir. | `[e.g., output]` |
| `-g,--genome` | Genome file of fasta file. | `[e.g., genome.fasta]` |
| `-d,--index` | Path prefix of a pre-built bowtie index (optional). If not specified, the index is built automatically under the output directory; provide a prefix here to reuse an existing index and skip building. | `[e.g., bowtie_index_prefix]` |
| `-t,--threads` | Number of threads used in mirdeep-p3 (optional, default: 1). | `[e.g., 14]` |
| `-p,--progress` | Number of samples/files processed in parallel (optional, default: 1). For example, if 3 input files are provided, `-p 3` will process all three simultaneously. | `[e.g., 3]` |

| Path | Description |
|---|---|
| `output/mirdp3-identification-<time>.pipe` | Log for this task |
| `output/sample_1/<sample>_trimming_report.txt` | trim_galore reslut |
| `output/sample_1/<sample>_identification.log` | Log for each step |
| `output/sample_1/<sample>.total_reads` | Total reads count |
| `output/sample_1/<sample>_filter_P_prediction` | Result of indentification step |

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
| `-i,--input` | Raw sequencing data (FASTQ/FASTA/compressed files，multiple samples can be separated by commas). | `[e.g., sample_1.fq or sample_1.fq,sample_2.fq,sample_3.fq]` |
| `-o,--output` | Output dir. | `[e.g., output]` |
| `-g,--genome` | Genome file of fasta file. | `[e.g., genome.fasta]` |
| `-d,--index` | Path prefix of a pre-built bowtie index (optional). If not specified, the index is built automatically under the output directory; provide a prefix here to reuse an existing index and skip building. | `[e.g., bowtie_index_prefix]` |
| `-t,--threads` | Number of threads used in mirdeep-p3 (optional, default: 1). | `[e.g., 14]` |
| `--species` | Species name in the output (must be quoted). | `[e.g., "Arabidopsis thaliana"]` |
| `--prefix_miRNA` | miRNA prefix in the output (must be quoted). | `[e.g., "Ath"]` |
| `--prefix` | Prefix of output file (optional). | `[e.g., flower]` |

| Path | Description |
|---|---|
| `output/annotation/mirdp3-annotation-<time>.pipe` | Log for this task |
| `output/annotation/<prefix>/<prefix>-basic-info` | Result of annotation step |
| `output/annotation/<prefix>/<prefix>_annotation.log` | Log for each step |
| `output/annotation/<prefix>/<prefix>-basic-info-cluster` | miRNA cluster result |
| `output/annotation/<prefix>/<prefix>-mature.count | Raed count matrix |
| `output/annotation/<prefix>/<prefix>-mature.exp | Expression matrix |

### Advanced
- **Custom reference build**: pass your own genome with `-g`; if `-d`
  is omitted, the bowtie index is auto-built under the output/index directory.
- **Skip steps**: `--no-reads_clean` skip the reads clean step.
- **Multi group**: `-r/--replicate` used for multigroup; if -r 3,3: the first three inputs form group A, the last three form group B. Each group is processed and reported independently.
- **prefix**: one prefix per group, comma-separated. The number of prefixes must match the number of groups (2 groups → 2 prefixes).
- **--common**: when more than one group is present, `--common` unifies the naming of non-conserved miRNA families across groups, ensuring family names are consistent between groups for downstream comparison.
- **--consistency**: path to a reference basic-info file produced by a previous `mirdeep-p3 annotation` run (same format). When provided, the current annotation is forced to keep miRNA naming consistent with that reference — i.e. the same miRNA will receive the same name in both datasets. This is useful when re-annotating the same species or when merging annotations across groups.

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
CI (GitHub Actions) should validate:
- Lint/format
- Unit tests
- Small end-to-end toy dataset

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
