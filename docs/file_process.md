# File process

> This guide describes how to preprocess files and convert formats.

## Reformat input fasta
```bash
python src/preprocess_files.py \
  -i input.fa \
  -o output.fa ##reformat fasta
```

## Basic-info file to bed format
```bash
python bin/extract_columns.py \
-i /path/to/basic-info \
-o output.bed
```
