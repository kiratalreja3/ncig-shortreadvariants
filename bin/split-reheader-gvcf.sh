#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status

# Usage statement function
usage() {
    echo "Usage: $0 -input <inputgvcf> -output <outputgvcf> -chr <chromosome>"
    echo
    echo "This script extracts a specific chromosome from an input GVCF file,"
    echo "reheaders it to include only relevant contigs, and indexes the output."
    echo
    echo "Arguments:"
    echo "  -input   Path to the input GVCF file (required)"
    echo "  -output  Path to the output reheadered GVCF file (required)"
    echo "  -chr     Chromosome to extract (e.g., chr1, chrX) (required)"
    echo "  -h, --help  Show this help message and exit"
    echo
    exit 1
}

# Check if no arguments are provided
if [[ "$#" -eq 0 ]]; then
    echo "Error: No arguments provided."
    usage
fi

# Initialize variables
inputgvcf=""
outputgvcf=""
chr=""

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -input) shift; inputgvcf="$1" ;;
        -output) shift; outputgvcf="$1" ;;
        -chr) shift; chr="$1" ;;
        -h|--help) usage ;;  # Help flag
        *) 
            echo "Error: Unknown option $1"
            usage
            ;;
    esac
    shift
done

# Validate required arguments
if [[ -z "$inputgvcf" || -z "$outputgvcf" || -z "$chr" ]]; then
    echo "Error: Missing required arguments."
    usage
fi


tmpfile="${PBS_JOBFS}/tmp.${outputgvcf}"

bcftools view -Oz -o "${tmpfile}" -r ${chr} ${inputgvcf}

# Step 3: Extract the header from the subset
bcftools view -h "${tmpfile}" > ${PBS_JOBFS}/${chr}.header.txt

# Step 4: Generate a list of contigs present in the subset
bcftools query -f '%CHROM\n' "${tmpfile}" | sort -u > ${PBS_JOBFS}/${chr}.contigs.list

# Step 5: Clean the header to include only relevant contigs
awk -v contigs_file="${PBS_JOBFS}/${chr}.contigs.list" '
  BEGIN {
    while ((getline line < contigs_file) > 0) contigs[line] = 1;
    close(contigs_file);
  }
  /^##contig/ {
    match($0, /ID=([^,>]+)/, arr);
    if (arr[1] in contigs) print;
  }
  !/^##contig/ {print}
' ${PBS_JOBFS}/${chr}.header.txt > ${PBS_JOBFS}/${chr}.filtered_header.txt

# Step 6: Replace the header with the cleaned version
bcftools reheader -h ${PBS_JOBFS}/${chr}.filtered_header.txt -o "$outputgvcf" "${tmpfile}"

# Step 7: Index the final cleaned GVCF
tabix -p vcf "$outputgvcf"
