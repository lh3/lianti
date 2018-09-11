#!/bin/bash

input_file="$1"
genome_file="/n/home01/hl141/lh3/bwadb/hs37m.fa"

grep "^NV" ${input_file} | awk -F $'\t' 'BEGIN {OFS = FS; reverse_complement["A"]="T"; reverse_complement["C"]="G"; reverse_complement["G"]="C"; reverse_complement["T"]="A";} {if ($4=="A"||$4=="C"){print $2,$3-2,$3+1,$1 ":" $4 $5,"0","+";} else {print $2,$3-2,$3+1,$1 ":" reverse_complement[$4] reverse_complement[$5],"0","-";}}' | bedtools getfasta -s -name -bed stdin -fi ${genome_file} | awk '{if(NR%2==1){split($0,name_array,":");}else{print name_array[2] "_" $0;}}' | sort | uniq -c | awk 'BEGIN {OFS = "\t"} {print $2,$1;}'
