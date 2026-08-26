# Construct new core dataset and evaluation

> Moved from README. This guide describes how to build and evaluate
> a new PmiREN core dataset, and how to prepare input data.

## Construct new core dataset and evaluation
```bash
python scripts/mirdp3_core_build.py \
  -p data/PmiREN-20260810-isoform.fa \ ##old core dataset
  -i examples/PmiREN2.0_basic_info_MIR_unqiue.fasta \ ##non-redundant novel miRNAs dataset
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
  -i examples/PmiREN2.0_basic_info.fasta \
  -o examples/PmiREN2.0_basic_info_MIR.fasta
python bin/dedup_seq.py \
  -i examples/PmiREN2.0_basic_info_MIR.fasta \
  --uni examples/PmiREN2.0_basic_info_MIR_unique.fasta \
  --dup examples/PmiREN2.0_basic_info_MIR_dup.fasta

##Remove the same seq(optional)
python bin/remove_matched_seq.py \
  -i examples/PmiREN2.0_basic_info.fasta \
  -r data/PmiREN-20260810-isoform.fa \
  -o examples/PmiREN2.0.fasta

##Check the final miRNA family number
cat data/PmiREN-20260810-isoform.fa | grep ">" | sort -V | tail
```
