#!/bin/bash

# Check if the file is provided as an argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 <file>"
    exit 1
fi

file=$1

# Use awk to process the file
awk '
NR > 1 { 
    for (i = 2; i <= NF; i++) {
        sum[i] += $i
    }
} 
END {
    for (i = 2; i <= NF; i++) {
        printf "Sum of column %d: %s\n", i, sum[i]
    }
}
' "$file"
