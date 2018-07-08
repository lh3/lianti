## Getting Started

```sh
git clone https://github.com/lh3/lianti
cd lianti && make
# preprocessing, mapping and marking duplicates
seqtk mergepe read1.fq.gz read2.fq.gz | ./lianti trim - | bwa mem -Cpt8 ref.fa - \
  | samtools view -uS - | sambamba sort /dev/stdin -o /dev/stdout | ./lianti ldup - > aln.bam
# calling SNVs
./lianti pileup -ycf ref.fa -P20 -L1 bulk.bam lianti.bam > raw.vcf
k8 plp-diff.js raw.vcf > filtered.txt
```

## Introduction

[LIANTI][lianti-paper] is a single-cell whole-genome amplification method.
This repo implements tools to preprocess raw LIANTI sequence data and to
call sequence variations from the alignment. Probably you would like to use the
`trim` command to trim adapters, identify barcodes and merge overlapping read
ends. It is non-trivial to reimplement these tedious functionality on your own.
`ldup` marks PCR duplicates in a barcode-aware manner. It has been superseded
by the `ldup` command in the [adna][adna] repo which is more general. You may
consider to call SNVs with this toolkit, too, but it is not that hard to roll
your own anyway. Calling SVs and CNVs is hard with any callers. This repo does
consider some LIANTI-specific features, but generally you should not expect it
to be the state of art. Good luck.

[adna]: https://github.com/DReichLab/adna
[lianti-paper]: http://science.sciencemag.org/content/356/6334/189
