#!/bin/bash

input_file="$1"
bed_file="$2"

# number of SNVs overlapping with the BED file
region_number="$(grep "^NV" ${input_file} | awk -F $'\t' 'BEGIN {OFS = FS} {print $2,$3-1,$3}' | bedtools intersect -u -a stdin -b ${bed_file} | wc -l)"

# total number of SNVs
total_number="$(grep "^NV" ${input_file} | awk -F $'\t' 'BEGIN {OFS = FS} {print $2,$3-1,$3}' | wc -l)"

# fraction of SNVs overlapping with the BED file
percentage="$(bc -l <<< ${region_number}/${total_number}*100 | awk '{printf "%f", $0}')"

# output
printf "#number in BED\ttotal number\tpercentage in BED\n"
printf "%s\t%s\t%s\n" "${region_number}" "${total_number}" "${percentage}"
