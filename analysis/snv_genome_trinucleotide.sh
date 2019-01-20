#!/bin/bash

genome_file="/n/home01/hl141/lh3/bwadb/hs37m.fa"
bed_file="/n/home01/hl141/lh3/bwadb/um75-hs37d5-comp.auto.bed.gz"

gunzip -c ${bed_file} | awk -F $'\t' 'BEGIN {OFS = FS;strands[1]="+";strands[2]="-"} {for(i=$2;i<$3;i++){for(j=1;j<=2;j++){print $1,i-1,i+2,".","0",strands[j]}}}' | bedtools getfasta -s -name -bed stdin -fi ${genome_file} |  awk '{if(NR%2==0){print substr($0,2,1) "_" $0}}' | grep '^[ACGT]_[ACGT]*$' | awk '{print} END {nucleotides[1]="A";nucleotides[2]="C";nucleotides[3]="G";nucleotides[4]="T";for(i=1;i<=4;i++){for(k=1;k<=4;k++){for(l=1;l<=4;l++){print nucleotides[i] "_" nucleotides[k] nucleotides[i] nucleotides[l]}}}}' | awk -F $'\t' 'BEGIN {OFS = FS}{counts[$0]++} END {for(x in counts){print x,counts[x]-1}}' | sort -k1,1

