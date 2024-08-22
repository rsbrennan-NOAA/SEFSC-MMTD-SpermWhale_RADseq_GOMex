#!/bin/bash

# Check if the mapping file was given
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <mapping_file>"
    exit 1
fi

mapping_file="$1"

# Read the mapping file and make associative array
declare -A name_map
while IFS=$'\t' read -r key value; do
    name_map["$key"]="$value"
done < "$mapping_file"

# Iterate over all fastq.gz files
for file in *_R1_001.fastq.gz; do
    # Extract the barcode, the 2nd entry between _ (CTAGTAC-AAGGAGT)
    barcode=$(echo "$file" | cut -f 2 -d '_')
    
    # Check if the barcode exists
    if [[ -v "name_map[$barcode]" ]]; then
        new_name="${name_map[$barcode]}.R1.fq.gz"
        mv "$file" "$new_name"
        echo "Renamed $file to $new_name"
    else
        echo "No name found for $file"
    fi
done

echo "Rename done"
