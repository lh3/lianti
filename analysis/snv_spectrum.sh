#!/bin/bash

input_file="$1"



awk -F $'\t' 'BEGIN {OFS = FS; reverse_complement["A"]="T"; reverse_complement["C"]="G"; reverse_complement["G"]="C"; reverse_complement["T"]="A";} {if ($1=="A"||$1=="C"){print $1 $2;} else {print reverse_complement[$1] reverse_complement[$2]}}' ${input_file} | awk '{print} END {print "AC"; print "AG"; print "AT"; print "CA"; print "CG"; print "CT"}' | sort | uniq -c | awk 'BEGIN {OFS = "\t"} {print $2,$1-1;}'
