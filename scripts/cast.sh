#!/bin/bash

# Step 1: Read all the files and combine into a temporary file
tmp_file="temp_variant_data.tsv"
output_file="casted_variant_report.tsv"
> $tmp_file

for file in variants_report_*; do
    if [[ $file =~ variants_report_(.*)_(.*).tsv ]]; then
        X_value=${BASH_REMATCH[1]}
        perm_value=${BASH_REMATCH[2]}
        
        while IFS=$'\t' read -r line; do
            echo -e "$line\t$X_value\t$perm_value" >> $tmp_file
        done < "$file"
    fi
done

# Step 2: Pivot the data
awk '
BEGIN {
    FS=OFS="\t"
}
NR==1 {
    # Store header and add new columns
    header = $0
    print header, "X0", "perm0", "X1", "perm1", "X2", "perm2", "X3", "perm3", "X4", "perm4", "X5", "perm5", "X6", "perm6", "X7", "perm7", "X8", "perm8", "X9", "perm9"
    next
}
{
    key = $1 "\t" $2 "\t" $3
    if (!(key in data)) {
        data[key] = $0
    }
    X_index = "X" $NF
    perm_index = "perm" $NF
    data[key] = data[key] "\t" $X_index "\t" $perm_index
}
END {
    for (key in data) {
        print data[key]
    }
}
' $tmp_file > $output_file

# Clean up temporary file
rm $tmp_file
