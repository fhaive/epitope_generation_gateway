import os
import subprocess
import pandas as pd


# Define the results directory, defaulting to "results" if not specified
RESULTS_DIR = config.get("results_dir", "results")
RNA_ANALYSIS_DIR = f"{RESULTS_DIR}/1B_RNA_fusion_HLA"
ALT_SPLICING_DIR = f"{RESULTS_DIR}/1C_alternative_splicing"
# Load the CSV file from the command line config option
sample_df = pd.read_csv(config["csvfile"])

# Load the appropiate local conda environment
conda_env = "../envs/module_1C.yaml"

# Create a dictionary for sample inputs
sample_dict = {
    (row['sample_name'], row['datatype'], row['lane']): {
        'fastq_R1': row['fastq_R1'],
        'fastq_R2': row['fastq_R2'],
        'datatype': row['datatype']
    }
    for idx, row in sample_df.iterrows()
}

# Rule to generate all output files for rule all
def generate_all_outputs():
    outputs = []
    for (sample, datatype, lane) in sample_dict.keys():
        if datatype == "CancerRNA":  # Only include outputs for CancerRNA
            outputs.append(f"{ALT_SPLICING_DIR}/Spladder/{sample}_CancerRNA/")
            outputs.append(f"{ALT_SPLICING_DIR}/VCF/{sample}_CancerRNA_altsplicing.vcf")
            outputs.append(f"{ALT_SPLICING_DIR}/VCF/{sample}_CancerRNA_altsplicing.vcf.gz")
    return outputs

# Rule all to list all output files
rule all:
    input:
        generate_all_outputs()

def get_lanes(sample):
    # Filter the dataframe for the sample name and the CancerRNA datatype
    lanes = sample_df[
        (sample_df['sample_name'] == sample) &
        (sample_df['datatype'] == 'CancerRNA')
    ]['lane'].tolist()

    return lanes


rule spladder:
    input:
        sorted_bam= lambda wildcards: expand(
            f"{RNA_ANALYSIS_DIR}/sorted_bam/{{sample}}_CancerRNA_{{lane}}_sorted.bam",
            sample=wildcards.sample,
            lane=get_lanes(wildcards.sample)
        ),
        annotation = f"resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf"
    output:
        sample_dir= directory(f"{ALT_SPLICING_DIR}/Spladder/{{sample}}_CancerRNA/")
    conda:
        conda_env
    shell:
        """
        spladder build \
            -o {output.sample_dir} \
            --parallel 20 \
            -b {input.sorted_bam} \
            -a {input.annotation} \
            --filter-overlap-exons \
            --no-primary-only \
            --quantify-graph \
            --confidence 3 \
            --iterations 5 \
            --ase-edge-limit 250 \
            --qmode all
        """



# altsplc2vcf.py script obtained from Scanneo2 here https://github.com/ylab-hi/ScanNeo2/blob/main/workflow/scripts/altsplc2vcf.py
rule conversion_to_vcf:
    input:
        spladder_dir=f"{ALT_SPLICING_DIR}/Spladder/{{sample}}_CancerRNA/"
    output:
        vcf_file=f"{ALT_SPLICING_DIR}/VCF/{{sample}}_CancerRNA_altsplicing.vcf"
    conda:
        conda_env
    shell:
        """
        python scripts/final_scripts/altsplc2vcf.py \
            -i {input.spladder_dir} \
            -r resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa \
            -g {wildcards.sample}_CancerRNA \
            -o {output.vcf_file}
        """


rule sort_vcf:
    input:
        vcf_file=f"{ALT_SPLICING_DIR}/VCF/{{sample}}_CancerRNA_altsplicing.vcf"
    output:
        sorted_vcf=f"{ALT_SPLICING_DIR}/VCF/{{sample}}_CancerRNA_altsplicing.vcf.gz"
    conda:
        conda_env
    shell:
        """
        bcftools sort {input.vcf_file} -o - | bcftools view -O z -o {output.sorted_vcf}
        """


