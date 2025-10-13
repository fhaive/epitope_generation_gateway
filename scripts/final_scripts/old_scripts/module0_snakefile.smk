import os
import pandas as pd

# Load the CSV file from the command line config option
sample_df = pd.read_csv(config["csvfile"])

# Create a dictionary for sample inputs
sample_dict = {
    (row['sample_name'], row['datatype'], row['lane']): {
        'fastq_R1': row['fastq_R1'],
        'fastq_R2': row['fastq_R2'],
    }
    for idx, row in sample_df.iterrows()
}

# Function to generate all output files for rule all
def generate_all_outputs():
    outputs = []
    for (sample, datatype, lane) in sample_dict.keys():
        outputs.append(f"results/qc/{sample}_{datatype}_{lane}_R1_fastqc.html")
        outputs.append(f"results/qc/{sample}_{datatype}_{lane}_R2_fastqc.html")
        outputs.append(f"results/trimmed_qc/{sample}_{datatype}_{lane}_trimmed_R1_fastqc.html")
        outputs.append(f"results/trimmed_qc/{sample}_{datatype}_{lane}_trimmed_R2_fastqc.html")
        outputs.append(f"results/trimmed_fastq/{sample}_{datatype}_{lane}_trimmed_R1.fastq.gz")
        outputs.append(f"results/trimmed_fastq/{sample}_{datatype}_{lane}_trimmed_R2.fastq.gz")
    outputs.append("results/qc/multiqc_fastqc_before.html")
    outputs.append("results/trimmed_qc/multiqc_fastqc_after.html")
    outputs.append("results/trimmed_fastq/multiqc_fastp.html")
    return outputs

# Rule all to list all output files (No wildcards, only concrete paths)
rule all:
    input:
        generate_all_outputs()

# Rule to run FastQC before trimming
rule fastqc_before:
    input:
        R1=lambda wildcards: sample_dict[(wildcards.sample, wildcards.datatype, wildcards.lane)]['fastq_R1'],
        R2=lambda wildcards: sample_dict[(wildcards.sample, wildcards.datatype, wildcards.lane)]['fastq_R2']
    output:
        R1_report="results/qc/{sample}_{datatype}_{lane}_R1_fastqc.html",
        R1_zip="results/qc/{sample}_{datatype}_{lane}_R1_fastqc.zip",
        R2_report="results/qc/{sample}_{datatype}_{lane}_R2_fastqc.html",
        R2_zip="results/qc/{sample}_{datatype}_{lane}_R2_fastqc.zip"
    conda:
        "../envs/module_0.yaml"
    shell:
        """
        fastqc {input.R1} {input.R2} --outdir=results/qc/
        base_r1=$(basename {input.R1} | sed -E 's/.fq.gz$|.fastq.gz$|.fq$|.fastq$//')
        base_r2=$(basename {input.R2} | sed -E 's/.fq.gz$|.fastq.gz$|.fq$|.fastq$//')
        mv results/qc/${{base_r1}}_fastqc.html {output.R1_report}
        mv results/qc/${{base_r1}}_fastqc.zip {output.R1_zip}
        mv results/qc/${{base_r2}}_fastqc.html {output.R2_report}
        mv results/qc/${{base_r2}}_fastqc.zip {output.R2_zip}
        """

# Rule to run fastp for paired-end data
rule fastp_trim:
    input:
        R1=lambda wildcards: sample_dict[(wildcards.sample, wildcards.datatype, wildcards.lane)]['fastq_R1'],
        R2=lambda wildcards: sample_dict[(wildcards.sample, wildcards.datatype, wildcards.lane)]['fastq_R2']
    output:
        trimmed_R1="results/trimmed_fastq/{sample}_{datatype}_{lane}_trimmed_R1.fastq.gz",
        trimmed_R2="results/trimmed_fastq/{sample}_{datatype}_{lane}_trimmed_R2.fastq.gz",
        json="results/trimmed_fastq/{sample}_{datatype}_{lane}_fastp.json",
        html="results/trimmed_fastq/{sample}_{datatype}_{lane}_fastp.html"
    conda:
        "../envs/module_0.yaml"
    shell:
        """
        fastp -i {input.R1} -I {input.R2} \
              -o {output.trimmed_R1} -O {output.trimmed_R2} \
              -l 50 --detect_adapter_for_pe \
              --json {output.json} --html {output.html} -q 20 -u 30
        """

# Rule to run FastQC after trimming
rule fastqc_after:
    input:
        R1="results/trimmed_fastq/{sample}_{datatype}_{lane}_trimmed_R1.fastq.gz",
        R2="results/trimmed_fastq/{sample}_{datatype}_{lane}_trimmed_R2.fastq.gz"
    output:
        R1_report="results/trimmed_qc/{sample}_{datatype}_{lane}_trimmed_R1_fastqc.html",
        R1_zip="results/trimmed_qc/{sample}_{datatype}_{lane}_trimmed_R1_fastqc.zip",
        R2_report="results/trimmed_qc/{sample}_{datatype}_{lane}_trimmed_R2_fastqc.html",
        R2_zip="results/trimmed_qc/{sample}_{datatype}_{lane}_trimmed_R2_fastqc.zip"
    conda:
        "../envs/module_0.yaml"
    shell:
        """
        fastqc {input.R1} {input.R2} --outdir=results/trimmed_qc/
        """

# Rule to aggregate FastQC results before trimming using MultiQC
rule multiqc_fastqc_before:
    input:
        [f"results/qc/{sample}_{datatype}_{lane}_R1_fastqc.html" for sample, datatype, lane in sample_dict.keys()],
        [f"results/qc/{sample}_{datatype}_{lane}_R2_fastqc.html" for sample, datatype, lane in sample_dict.keys()]
    output:
        "results/qc/multiqc_fastqc_before.html"
    conda:
        "../envs/module_0.yaml"
    params:
        latency_wait=600
    shell:
        """
        multiqc results/qc/ -o results/qc/ --filename multiqc_fastqc_before --force
        """

# Rule to aggregate FastQC results after trimming using MultiQC
rule multiqc_fastqc_after:
    input:
        [f"results/trimmed_qc/{sample}_{datatype}_{lane}_trimmed_R1_fastqc.html" for sample, datatype, lane in sample_dict.keys()],
        [f"results/trimmed_qc/{sample}_{datatype}_{lane}_trimmed_R2_fastqc.html" for sample, datatype, lane in sample_dict.keys()]
    output:
        "results/trimmed_qc/multiqc_fastqc_after.html"
    conda:
        "../envs/module_0.yaml"
    params:
        latency_wait=600
    shell:
        """
        multiqc results/trimmed_qc/ -o results/trimmed_qc/ --filename multiqc_fastqc_after --force
        """

# Rule to aggregate fastp results using MultiQC
rule multiqc_fastp:
    input:
        [f"results/trimmed_fastq/{sample}_{datatype}_{lane}_fastp.json" for sample, datatype, lane in sample_dict.keys()],
        [f"results/trimmed_fastq/{sample}_{datatype}_{lane}_fastp.html" for sample, datatype, lane in sample_dict.keys()]
    output:
        "results/trimmed_fastq/multiqc_fastp.html"
    conda:
        "../envs/module_0.yaml"
    params:
        latency_wait=600
    shell:
        """
        multiqc results/trimmed_fastq/ -o results/trimmed_fastq/ --filename multiqc_fastp --force
        """

# Increase latency wait for file system delays
latency_wait = 600
