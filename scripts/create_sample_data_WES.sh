#!/bin/bash

# Source directory containing raw FASTQ files
SOURCE_DIR="/nasdata/lmannino/Minna_Genomics/X204SC23123384-Z01-F001/01.RawData"

# Destination directory where processed files will be copied
DEST_DIR="/nasdata/lmannino/pipelines/epitope_generation_pipeline/raw_data"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Function to process files
process_files() {
    local src_dir="$1"
    local dest_dir="$2"

    # Loop through all gzipped FASTQ files in the source directory
    find "$src_dir" -type f -name "*.fq.gz" | while read -r file; do
        # Calculate the relative path from the source directory
        relative_path="${file#$src_dir/}"
        
        # Define the destination path, preserving the relative directory structure
        dest_file="$dest_dir/$relative_path"

        # Create the destination directory structure if it doesn't exist
        mkdir -p "$(dirname "$dest_file")"

        echo "Processing $file"

        # Extract the first 100000 reads and save to the destination directory
        zcat "$file" | head -n 400000 | gzip > "$dest_file"

        echo "Saved processed file to $dest_file"
    done
}

# Start processing
process_files "$SOURCE_DIR" "$DEST_DIR"

echo "Processing completed."
