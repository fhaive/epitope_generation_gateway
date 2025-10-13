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
        outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/{sample}_{datatype}_{lane}_R1_fastqc.html")
        outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/{sample}_{datatype}_{lane}_R2_fastqc.html")
        outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/{sample}_{datatype}_{lane}_trimmed_R1_fastqc.html")
        outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/{sample}_{datatype}_{lane}_trimmed_R2_fastqc.html")
        outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{sample}_{datatype}_{lane}_trimmed_R1.fastq.gz")
        outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{sample}_{datatype}_{lane}_trimmed_R2.fastq.gz")
#    outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/multiqc_fastqc_pre_filtering.html")
#    outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/multiqc_fastqc_post_filtering.html")
#    outputs.append(f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/multiqc_fastp.html")
    return outputs

# Rule all to list all output files (No wildcards, only concrete paths)
rule all:
    input:
        generate_all_outputs()

#Rule to run FastQC before trimming
rule fastqc_before:
    input:
        R1=lambda wildcards: sample_dict[(wildcards.sample, wildcards.datatype, wildcards.lane)]['fastq_R1'],
        R2=lambda wildcards: sample_dict[(wildcards.sample, wildcards.datatype, wildcards.lane)]['fastq_R2']
    output:
        R1_report=f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/{{sample}}_{{datatype}}_{{lane}}_R1_fastqc.html",
        R1_zip=f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/{{sample}}_{{datatype}}_{{lane}}_R1_fastqc.zip",
        R2_report=f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/{{sample}}_{{datatype}}_{{lane}}_R2_fastqc.html",
        R2_zip=f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/{{sample}}_{{datatype}}_{{lane}}_R2_fastqc.zip"
    conda:
        "../envs/module_0.yaml"
    resources:
        fastqc_slots=1
    shell:
        """
        fastqc --threads {threads} {input.R1} --outdir={RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/
        fastqc --threads {threads} {input.R2} --outdir={RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/

        # Wait for FastQC output files to appear
        base_r1=$(basename {input.R1} | sed -E 's/.fq.gz$|.fastq.gz$|.fq$|.fastq$//')
        base_r2=$(basename {input.R2} | sed -E 's/.fq.gz$|.fastq.gz$|.fq$|.fastq$//')

        echo ${{base_r2}}
        echo ${{base_r1}}
         # Wait for FastQC output files to appear
#        until [[ -f {RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/${{base_r1}}_fastqc.html ]]; do sleep 1; done
#        until [[ -f {RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/${{base_r2}}_fastqc.html ]]; do sleep 1; done

        cp {RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/${{base_r1}}_fastqc.html {output.R1_report}
        cp {RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/${{base_r1}}_fastqc.zip {output.R1_zip}
        cp {RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/${{base_r2}}_fastqc.html {output.R2_report}
        cp {RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/${{base_r2}}_fastqc.zip {output.R2_zip}

        """


# Rule to run fastp for paired-end data
rule fastp_trim:
    input:
        R1=lambda wildcards: sample_dict[(wildcards.sample, wildcards.datatype, wildcards.lane)]['fastq_R1'],
        R2=lambda wildcards: sample_dict[(wildcards.sample, wildcards.datatype, wildcards.lane)]['fastq_R2']
    output:
        trimmed_R1=f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{{sample}}_{{datatype}}_{{lane}}_trimmed_R1.fastq.gz",
        trimmed_R2=f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{{sample}}_{{datatype}}_{{lane}}_trimmed_R2.fastq.gz",
        json=f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{{sample}}_{{datatype}}_{{lane}}_fastp.json",
        html=f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{{sample}}_{{datatype}}_{{lane}}_fastp.html"
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
        echo "Running fastp command:"
        fastp -i {input.R1} -I {input.R2} \
          -o {output.trimmed_R1} -O {output.trimmed_R2} \
          -l {params.min_length} {params.detect_adapter} \
          --json {output.json} --html {output.html} \
          -q {params.quality} -u {params.unqualified_base_limit} \
          {params.custom_adapter_R1} {params.custom_adapter_R2}
        """

# Rule to run FastQC after trimming
rule fastqc_after:
    input:
        R1=f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{{sample}}_{{datatype}}_{{lane}}_trimmed_R1.fastq.gz",
        R2=f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{{sample}}_{{datatype}}_{{lane}}_trimmed_R2.fastq.gz"
    output:
        R1_report=f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/{{sample}}_{{datatype}}_{{lane}}_trimmed_R1_fastqc.html",
        R1_zip=f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/{{sample}}_{{datatype}}_{{lane}}_trimmed_R1_fastqc.zip",
        R2_report=f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/{{sample}}_{{datatype}}_{{lane}}_trimmed_R2_fastqc.html",
        R2_zip=f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/{{sample}}_{{datatype}}_{{lane}}_trimmed_R2_fastqc.zip"
    conda:
        "../envs/module_0.yaml"
    shell:
        """
        fastqc {input.R1} {input.R2} --outdir={RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/
        """

# Rule to aggregate FastQC results before trimming using MultiQC
rule multiqc_fastqc_before:
    input:
        [f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/{sample}_{datatype}_{lane}_R1_fastqc.html" for sample, datatype, lane in sample_dict.keys()],
        [f"{RESULTS_DIR}/0_Filtering_and_QC/qc_pre_filtering/{sample}_{datatype}_{lane}_R2_fastqc.html" for sample, datatype, lane in sample_dict.keys()] 
        #can't reference pre filtering names
        
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

# Rule to aggregate FastQC results after trimming using MultiQC
rule multiqc_fastqc_after:
    input:
        [f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/{sample}_{datatype}_{lane}_trimmed_R1_fastqc.html" for sample, datatype, lane in sample_dict.keys()],
        [f"{RESULTS_DIR}/0_Filtering_and_QC/qc_post_filtering/{sample}_{datatype}_{lane}_trimmed_R2_fastqc.html" for sample, datatype, lane in sample_dict.keys()]
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

# Rule to aggregate fastp results using MultiQC
rule multiqc_fastp:
    input:
        [f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{sample}_{datatype}_{lane}_fastp.json" for sample, datatype, lane in sample_dict.keys()],
        [f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{sample}_{datatype}_{lane}_fastp.html" for sample, datatype, lane in sample_dict.keys()]
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

# Increase latency wait for file system delays
latency_wait = 600
