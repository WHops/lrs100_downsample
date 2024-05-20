#!/bin/bash

output_file="merged_variant_report.tsv"
header_processed=false

for file in variants_report_*; do
    if [[ $file =~ variants_report_(.*)_(.*).tsv ]]; then
        X_value=${BASH_REMATCH[1]}
        perm_value=${BASH_REMATCH[2]}
        
        while IFS=$'\t' read -r -a line_array; do
            if [[ "$header_processed" == false ]]; then
                header_processed=true
                header=("${line_array[@]:0:3}" "X" "perm" "${line_array[@]:3}")
                echo -e "$(IFS=$'\t'; echo "${header[*]}")" > $output_file
            else
                # Skip the header line for subsequent files
                if [[ "${line_array[0]}" == "${header[0]}" ]]; then
                    continue
                fi
                line=("${line_array[@]:0:3}" "$X_value" "$perm_value" "${line_array[@]:3}")
                echo -e "$(IFS=$'\t'; echo "${line[*]}")" >> $output_file
            fi
        done < "$file"
    fi
done