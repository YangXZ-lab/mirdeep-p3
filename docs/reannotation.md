# Reannotation of a miRNA dataset

> Moved from README. Re-annotate an existing miRNA dataset against
> the PmiREN core database.

## Reannotation of a miRNA dataset
```bash
python bin/anno_miRNA.py \
  -i examples/PmiREN2.0_basic_info_MIR_unqiue.fasta \ ##miRNA dataset
  -p data/PmiREN-20260810-isoform.fa \ ##pmiREN core dataset
  -o test/anno_miRNA \ ##output
  --threads 14 \
  --type MIRN \  ##prefix of new miRNA family
  --prefix "Ath"  ##prefix of miRNA(optional)
###anno.fasta, renamed miRNA
###anno.map, rename map
```
