#!/usr/bin/env python
import pandas as pd
from gtfparse import read_gtf
import optparse

def format_attributes(row):
    """Format the attributes field, ensuring only relevant attributes are included."""
    
    attributes = []
    
    # Include only relevant fields, excluding columns like seqname, start, end, etc.
    exclude_keys = {"seqname", "source", "feature", "start", "end", "score", "strand", "frame"}

    # Ensure gene_id is formatted correctly
    if pd.notna(row.get("gene_id")):
        gene_id = row["gene_id"].split(".")[0]  # Remove version
        gene_version = row["gene_id"].split(".")[1] if "." in row["gene_id"] else ""
        attributes.append(f'gene_id "{gene_id}";')
        if gene_version:
            attributes.append(f'gene_version "{gene_version}";')

    # Ensure transcript_id is formatted correctly
    if pd.notna(row.get("transcript_id")) and row["transcript_id"] != "":
        transcript_id = row["transcript_id"].split(".")[0]  # Remove version
        transcript_version = row["transcript_id"].split(".")[1] if "." in row["transcript_id"] else ""
        attributes.append(f'transcript_id "{transcript_id}";')
        if transcript_version:
            attributes.append(f'transcript_version "{transcript_version}";')

    # Add other attributes in the correct order
    for col in row.index:
        if col not in exclude_keys and col not in ["gene_id", "transcript_id"] and pd.notna(row[col]) and row[col] != "":
            attributes.append(f'{col} "{row[col]}";')

    return " ".join(attributes)


def convert_gencode_to_ensembl(gtf_input, gtf_output):
    """
    Convert a Gencode GTF file to Ensembl format by removing gene_id and transcript_id versions
    and adding gene_version and transcript_version attributes.

    Args:
        gtf_input (str): Path to the input GTF file.
        gtf_output (str): Path to the output Ensembl-compatible GTF file.
    """

    # Read the GTF file as a pandas DataFrame
    gtf_df = read_gtf(gtf_input, result_type='pandas')  # Ensure it loads as Pandas

    # Format the attributes correctly
    gtf_df["attributes"] = gtf_df.apply(format_attributes, axis=1)

    # Define the required GTF columns
    required_columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]

    # Replace NaN values in 'score' and 'frame' with `.` to match GTF format
    gtf_df["score"] = gtf_df["score"].fillna(".").astype(str)
    gtf_df["frame"] = gtf_df["frame"].replace([0, "0", None], ".").astype(str)

    # Ensure all necessary columns are present
    gtf_df = gtf_df[required_columns]

    # Write the modified GTF back to a file
    gtf_df.to_csv(gtf_output, sep="\t", index=False, header=False, quoting=3)

    print(f"Conversion complete. Output saved to {gtf_output}")


def main():
    parser = optparse.OptionParser(usage="Usage: %prog -i <input.gtf> -o <output.gtf>", version="%prog 1.0")
    
    parser.add_option("-i", "--input", dest="input_gtf", help="Path to the input GTF file", metavar="FILE")
    parser.add_option("-o", "--output", dest="output_gtf", help="Path to the output GTF file", metavar="FILE")

    (options, args) = parser.parse_args()

    if not options.input_gtf or not options.output_gtf:
        parser.error("Both input and output GTF files are required. Use -h for help.")

    convert_gencode_to_ensembl(options.input_gtf, options.output_gtf)


if __name__ == "__main__":
    main()
