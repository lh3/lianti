#!/bin/bash

input_file="$1"

grep "^DV" ${input_file} | awk -F $'\t' 'BEGIN {OFS = FS; reverse_complement["A"]="T"; reverse_complement["C"]="G"; reverse_complement["G"]="C"; reverse_complement["T"]="A";} {split($8,a,":");split(a[2],b,",");if(b[1]==0){$4=reverse_complement[$4];$5=reverse_complement[$5]};if($4=="G"||$4=="T"){strand="-";$4=reverse_complement[$4];$5=reverse_complement[$5]}else{strand="+"};print $4""$5""strand} END{print "AC+"; print "AG+"; print "AT+"; print "CA+"; print "CG+"; print "CT+"; print "AC-"; print "AG-"; print "AT-"; print "CA-"; print "CG-"; print "CT-"}' | sort | uniq -c | awk 'BEGIN {OFS = "\t"} {print $2,$1-1;}'

