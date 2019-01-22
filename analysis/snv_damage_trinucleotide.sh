#!/bin/bash

input_file="$1"
genome_file="/n/home01/hl141/lh3/bwadb/hs37m.fa"

grep "^DV" ${input_file} | awk -F $'\t' 'BEGIN {OFS = FS; reverse_complement["A"]="T"; reverse_complement["C"]="G"; reverse_complement["G"]="C"; reverse_complement["T"]="A";} {split($8,a,":");split(a[2],b,",");strand="+";if(b[1]==0){strand="-";$4=reverse_complement[$4];$5=reverse_complement[$5]};print $2,$3-2,$3+1,$1 ":" $4 $5,"0",strand}' | bedtools getfasta -s -name -bed stdin -fi ${genome_file} | awk '{if(NR%2==1){split($0,name_array,"[:(]");}else{print name_array[2] "_" $0;}} END {nucleotides[1]="A";nucleotides[2]="C";nucleotides[3]="G";nucleotides[4]="T";for(i=1;i<=4;i++){for(j=1;j<=4;j++){if(i==j)continue;for(k=1;k<=4;k++){for(l=1;l<=4;l++){print nucleotides[i] nucleotides[j] "_" nucleotides[k] nucleotides[i] nucleotides[l]}}}}}' | sort | uniq -c | awk 'BEGIN {OFS = "\t"} {print $2,$1-1;}'



