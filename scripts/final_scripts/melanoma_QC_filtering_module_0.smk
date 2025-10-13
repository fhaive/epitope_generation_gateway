import os
import pandas as pd

# Load the CSV file from the command line config option
sample_df = pd.read_csv(config["csvfile"])

# Load config file
configfile: "scripts/final_scripts/config/QC_filtering_config_0.yaml"

# Define the results directory, defaulting to "results" if not specified
RESULTS_DIR = config.get("results_dir", "results")

# Ensure all necessary directories exist
os.makedirs(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering", exist_ok=True)
os.makedirs(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering", exist_ok=True)
os.makedirs(f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq", exist_ok=True)

# Create a dictionary for sample inputs using a unified key
sample_dict = {
    f"{row['sample_name']}|{row['datatype']}|{row['lane']}": {
        'fastq_R1': row['fastq_R1'],
        'fastq_R2': row['fastq_R2'],
    }
    for _, row in sample_df.iterrows()
}

# Function to generate all output files for rule all
def generate_all_outputs():
    outputs = []
    for sample_key in sample_dict.keys():
        outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/{sample_key}_R1_fastqc.html")
        outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/{sample_key}_R2_fastqc.html")
        outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/{sample_key}_trimmed_R1_fastqc.html")
        outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/{sample_key}_trimmed_R2_fastqc.html")
        outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{sample_key}_trimmed_R1.fastq.gz")
        outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{sample_key}_trimmed_R2.fastq.gz")
    outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/multiqc_fastqc_pre_filtering.html")
    outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/multiqc_fastqc_post_filtering.html")
    outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/multiqc_fastp.html")
    return outputs

rule all:
    input:
        generate_all_outputs()

rule fastqc_before:
    input:
        R1=lambda wc: sample_dict[wc.sample_key]['fastq_R1'],
        R2=lambda wc: sample_dict[wc.sample_key]['fastq_R2']
    output:
        R1_report=f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/{{sample_key}}_R1_fastqc.html",
        R1_zip=f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/{{sample_key}}_R1_fastqc.zip",
        R2_report=f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/{{sample_key}}_R2_fastqc.html",
        R2_zip=f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/{{sample_key}}_R2_fastqc.zip"
    conda:
        "../envs/module_0.yaml"
    resources:
        fastqc_slots=1
    shell:
        """
        fastqc --threads {threads} {input.R1} --outdir={RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/
        fastqc --threads {threads} {input.R2} --outdir={RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/

        base_r1=$(basename {input.R1} | sed -E 's/.fq.gz$|.fastq.gz$|.fq$|.fastq$//')
        base_r2=$(basename {input.R2} | sed -E 's/.fq.gz$|.fastq.gz$|.fq$|.fastq$//')

        cp {RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/${{base_r1}}_fastqc.html {output.R1_report}
        cp {RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/${{base_r1}}_fastqc.zip {output.R1_zip}
        cp {RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/${{base_r2}}_fastqc.html {output.R2_report}
        cp {RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/${{base_r2}}_fastqc.zip {output.R2_zip}
        """

rule fastp_trim:
    input:
        R1=lambda wc: sample_dict[wc.sample_key]['fastq_R1'],
        R2=lambda wc: sample_dict[wc.sample_key]['fastq_R2']
    output:
        trimmed_R1=f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{{sample_key}}_trimmed_R1.fastq.gz",
        trimmed_R2=f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{{sample_key}}_trimmed_R2.fastq.gz",
        json=f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{{sample_key}}_fastp.json",
        html=f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{{sample_key}}_fastp.html"
    params:
        detect_adapter="--detect_adapter_for_pe" if config.get("fastp", {}).get("detect_adapter", True) else "",
        custom_adapter_R1=f"--adapter_sequence {config.get('fastp', {}).get('custom_adapter_R1', '')}" if config.get("fastp", {}).get("custom_adapter_R1", "") else "",
        custom_adapter_R2=f"--adapter_sequence_r2 {config.get('fastp', {}).get('custom_adapter_R2', '')}" if config.get("fastp", {}).get("custom_adapter_R2", "") else "",
        min_length=config.get("fastp", {}).get("min_length", 50),
        quality=config.get("fastp", {}).get("quality", 20),
        unqualified_base_limit=config.get("fastp", {}).get("unqualified_base_limit", 30)
    conda:
        "../envs/module_0.yaml"
    shell:
        """
        fastp -i {input.R1} -I {input.R2} \
          -o {output.trimmed_R1} -O {output.trimmed_R2} \
          -l {params.min_length} {params.detect_adapter} \
          --json {output.json} --html {output.html} \
          -q {params.quality} -u {params.unqualified_base_limit} \
          {params.custom_adapter_R1} {params.custom_adapter_R2}
        """

rule fastqc_after:
    input:
        R1=f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{{sample_key}}_trimmed_R1.fastq.gz",
        R2=f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{{sample_key}}_trimmed_R2.fastq.gz"
    output:
        R1_report=f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/{{sample_key}}_trimmed_R1_fastqc.html",
        R1_zip=f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/{{sample_key}}_trimmed_R1_fastqc.zip",
        R2_report=f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/{{sample_key}}_trimmed_R2_fastqc.html",
        R2_zip=f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/{{sample_key}}_trimmed_R2_fastqc.zip"
    conda:
        "../envs/module_0.yaml"
    shell:
        """
        fastqc {input.R1} {input.R2} --outdir={RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/
        """

rule multiqc_fastqc_before:
    input:
        expand(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/{{sample_key}}_R1_fastqc.html", sample_key=sample_dict.keys()),
        expand(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/{{sample_key}}_R2_fastqc.html", sample_key=sample_dict.keys())
    output:
        f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/multiqc_fastqc_pre_filtering.html"
    conda:
        "../envs/module_0.yaml"
    shell:
        """
        multiqc {RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/ \
                -o {RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/ \
                --filename multiqc_fastqc_pre_filtering --force
        """

rule multiqc_fastqc_after:
    input:
        expand(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/{{sample_key}}_trimmed_R1_fastqc.html", sample_key=sample_dict.keys()),
        expand(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/{{sample_key}}_trimmed_R2_fastqc.html", sample_key=sample_dict.keys())
    output:
        f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/multiqc_fastqc_post_filtering.html"
    conda:
        "../envs/module_0.yaml"
    shell:
        """
        multiqc {RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/ \
                -o {RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/ \
                --filename multiqc_fastqc_post_filtering --force
        """

rule multiqc_fastp:
    input:
        expand(f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{{sample_key}}_fastp.json", sample_key=sample_dict.keys()),
        expand(f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{{sample_key}}_fastp.html", sample_key=sample_dict.keys())
    output:
        f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/multiqc_fastp.html"
    conda:
        "../envs/module_0.yaml"
    shell:
        """
        multiqc {RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/ \
                -o {RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/ \
                --filename multiqc_fastp --force
        """

latency_wait = 600
