#!/bin/bash

# Source directory containing raw FASTQ files
SOURCE_DIR="/nasdata/Omics/RNA_Seq_projects/Nittykoski_may_2023"

# Destination directory where processed files will be copied
DEST_DIR="/nasdata/lmannino/pipelines/epitope_generation_pipeline/raw_data"

# Function to process files
process_files() {
    local src_dir="$1"
    local dest_dir="$2"

    # Loop through all gzipped FASTQ files in the source directory
    find "$src_dir" -type f -name "*.fastq.gz" | while read -r file; do
        # Extract the base name of the file (without directory path)
        base_name=$(basename "$file")

        # Define the destination file path (in the DEST_DIR)
        dest_file="$dest_dir/$base_name"

        echo "Processing $file"

        # Extract the first 100,000 reads and save to the destination directory
        zcat "$file" | head -n 400000 | gzip > "$dest_file"

        echo "Saved processed file to $dest_file"
    done
}

# Start processing
process_files "$SOURCE_DIR" "$DEST_DIR"

echo "Processing completed."
