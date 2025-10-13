import os
import pandas as pd

global_tmpdir = "/mnt/tmp_luca"
conda_env = "../envs/module_1.yaml"
# Load the CSV file
sample_df = pd.read_csv(config["csvfile"])

# Filter out RNA samples and keep only DNA samples (e.g., CancerDNA, NormalDNA)
dna_sample_df = sample_df[sample_df['datatype'].str.contains("DNA")]

# Create a dictionary for sample inputs, pointing to the trimmed fastq files generated in module 0
sample_dict = {
    (row['sample_name'], row['datatype'], row['lane']): {
        'fastq_R1': f"results/trimmed_fastq/{row['sample_name']}_{row['datatype']}_{row['lane']}_trimmed_R1.fastq.gz",
        'fastq_R2': f"results/trimmed_fastq/{row['sample_name']}_{row['datatype']}_{row['lane']}_trimmed_R2.fastq.gz"
    }
    for idx, row in dna_sample_df.iterrows()  # Only iterate through valid DNA samples and lanes
    if os.path.exists(f"results/trimmed_fastq/{row['sample_name']}_{row['datatype']}_{row['lane']}_trimmed_R1.fastq.gz")
    and os.path.exists(f"results/trimmed_fastq/{row['sample_name']}_{row['datatype']}_{row['lane']}_trimmed_R2.fastq.gz")
}

# Create a filtered dataframe for only "NormalDNA" samples
normal_dna_sample_df = dna_sample_df[dna_sample_df['datatype'] == 'NormalDNA']

# Finding the reference genome
ref_genome = "resources/genome_DNA/GRCh38.p13.genome.fa"
dbsnp_vcf = "resources/dbSNP/Homo_sapiens_assembly38.dbsnp138.vcf"
germline_resource = "resources/gnomad/somatic-hg38_af-only-gnomad.hg38.vcf.gz"
panel_of_normals = "resources/PON/somatic-hg38_1000g_pon.hg38.vcf.gz"


# Function to get all lanes for a given sample and datatype
def get_lanes(sample, datatype):
    return [lane for (s, d, lane) in sample_dict.keys() if s == sample and d == datatype]
# Function to get RG string (simplified)
def get_rg_info(sample, lane, datatype):
    rg = f"@RG\\tID:{lane}\\tPL:ILLUMINA\\tSM:{sample}_{datatype}"
    return rg


# Function to generate all output files for rule all
def generate_all_outputs():
    outputs = []
    for (sample, datatype, lane) in sample_dict.keys():
        outputs.append(f"results/DNA_alignment/bwa/{sample}_{datatype}_{lane}.sam")
        outputs.append(f"results/DNA_alignment/bam/{sample}_{datatype}_{lane}.bam")
        outputs.append(f"results/DNA_alignment/bam/{sample}_{datatype}_{lane}_recal_data.table")
        outputs.append(f"results/DNA_alignment/bqsr/{sample}_{datatype}_final_bqsr.bam")
        outputs.append(f"results/DNA_alignment/bqsr/{sample}_{datatype}_final_bqsr.bam.bai")   # Use "final" for all lanes
        outputs.append(f"results/DNA_metrics/{sample}_{datatype}_alignment_metrics.txt")
        outputs.append(f"results/DNA_metrics/{sample}_{datatype}_insert_size_metrics.txt")
        outputs.append(f"results/DNA_metrics/{sample}_{datatype}_insert_size_histogram.pdf")
        
        # Add VCF output for NormalDNA only
        if datatype == "NormalDNA":
            outputs.append(f"results/VCF/{sample}_{datatype}_HC_raw_variants.vcf")
        if datatype == "CancerDNA":
            outputs.append(f"results/VCF/{sample}_somatic.vcf.gz")
    return outputs

# Rule all to list all output files (No wildcards, only concrete paths)
rule all:
    input:
        generate_all_outputs()

# Rule for aligning with bwa mem
rule bwa_mem:
    input:
        R1=lambda wildcards: sample_dict[(wildcards.sample, wildcards.datatype, wildcards.lane)]['fastq_R1'],
        R2=lambda wildcards: sample_dict[(wildcards.sample, wildcards.datatype, wildcards.lane)]['fastq_R2'],
        ref=ref_genome
    output:
        sam="results/DNA_alignment/bwa/{sample}_{datatype}_{lane}.sam"
    params:
       rg=lambda wildcards: get_rg_info(wildcards.sample, wildcards.lane, wildcards.datatype)  # Use RG string with lane and sample
    conda:
        conda_env
    threads: 32
    shell:
        """
        export TMPDIR={global_tmpdir}
        mkdir -p {global_tmpdir}
        mkdir -p results/DNA_alignment/bwa

        bwa mem -t {threads} -R "{params.rg}" {input.ref} {input.R1} {input.R2} > {output.sam}
        """

# Rule for MarkDuplicates
rule gatk_mark_duplicates:
    input:
        sam="results/DNA_alignment/bwa/{sample}_{datatype}_{lane}.sam"
    output:
        bam="results/DNA_alignment/bam/{sample}_{datatype}_{lane}.bam"
    conda:
        conda_env
    resources:
        markdup_slots=1
    threads: 20
    shell:
        """
        mkdir -p {global_tmpdir} results/DNA_alignment/bam
        gatk MarkDuplicatesSpark --java-options "-Xmx64G -XX:CICompilerCount={threads}" --conf spark.local.dir={global_tmpdir} -I {input.sam} -O {output.bam}

        #gatk MarkDuplicatesSpark --java-options "-Xmx256G -XX:CICompilerCount={threads}" -I {input.sam} -O {output.bam}

        #gatk MarkDuplicatesSpark --java-options "-Xmx128G -XX:CICompilerCount={threads}" --conf spark.executor.memoryOverhead=128G -I {input.sam} -O {output.bam}
        """
#check if MarkDuplicatesSpark works on bam
# Rule for BaseRecalibrator
rule gatk_base_recalibrator:
    input:
        bam="results/DNA_alignment/bam/{sample}_{datatype}_{lane}.bam",
        ref=ref_genome,
        known_sites=dbsnp_vcf
    output:
        recal_table="results/DNA_alignment/bam/{sample}_{datatype}_{lane}_recal_data.table"
    
    conda:
        conda_env
    shell:
        """
        gatk BaseRecalibrator -I {input.bam} -R {input.ref} --known-sites {input.known_sites} -O {output.recal_table}
        """

# Rule for ApplyBQSR
rule gatk_apply_bqsr:
    input:
        bam="results/DNA_alignment/bam/{sample}_{datatype}_{lane}.bam",
        recal_table="results/DNA_alignment/bam/{sample}_{datatype}_{lane}_recal_data.table",
        ref=ref_genome
    output:
        bam_bqsr="results/DNA_alignment/bqsr/{sample}_{datatype}_{lane}_bqsr_reads.bam"
    conda:
        conda_env
    shell:
        """
        gatk ApplyBQSR -I {input.bam} -R {input.ref} --bqsr-recal-file {input.recal_table} -O {output.bam_bqsr}
        """

# Rule for merging BAM files

rule merge_bam:
    input:
        lambda wildcards: expand("results/DNA_alignment/bqsr/{sample}_{datatype}_{lane}_bqsr_reads.bam", sample=wildcards.sample, datatype=wildcards.datatype, lane=get_lanes(wildcards.sample, wildcards.datatype))
    output:
        merged_bam="results/DNA_alignment/bqsr/{sample}_{datatype}_final_bqsr.bam"
    conda:
        conda_env
    shell:
        """
        if [ $(echo {input} | wc -w) -gt 1 ]; then
            # Build the command for MergeSamFiles with multiple input files
            inputs=$(for bam in {input}; do echo -n "I=$bam "; done)
            picard MergeSamFiles $inputs O={output.merged_bam}
            
            # Remove the individual unmerged BAMs after merging, will implement in the end
            #rm {input}
        else
            # If there's only one BAM, simply rename it to the final BAM file
            mv {input} {output.merged_bam}
        fi
        """

# Rule for CollectAlignmentSummaryMetrics
rule gatk_collect_alignment_metrics:
    input:
        bam="results/DNA_alignment/bqsr/{sample}_{datatype}_final_bqsr.bam"
    output:
        metrics="results/DNA_metrics/{sample}_{datatype}_alignment_metrics.txt"
    conda:
        conda_env
    shell:
        """
        mkdir -p results/DNA_metrics
        gatk CollectAlignmentSummaryMetrics -R {ref_genome} -I {input.bam} -O {output.metrics}
        """

# Rule for CollectInsertSizeMetrics
rule gatk_collect_insert_size_metrics:
    input:
        bam="results/DNA_alignment/bqsr/{sample}_{datatype}_final_bqsr.bam"
    output:
        metrics="results/DNA_metrics/{sample}_{datatype}_insert_size_metrics.txt",
        histogram="results/DNA_metrics/{sample}_{datatype}_insert_size_histogram.pdf"
    conda:
        conda_env
    shell:
        """
        mkdir -p results/DNA_metrics
        gatk CollectInsertSizeMetrics -I {input.bam} -O {output.metrics} -H {output.histogram}
        """

# Rule for indexing BAM files with samtools
rule samtools_index_bam:
    input:
        bam="results/DNA_alignment/bqsr/{sample}_{datatype}_final_bqsr.bam"
    output:
        index="results/DNA_alignment/bqsr/{sample}_{datatype}_final_bqsr.bam.bai"
    conda:
        conda_env
    shell:
        """
        samtools index {input.bam}
        """

# Rule for HaplotypeCaller on NormalDNA only
rule gatk_haplotypecaller:
    input:
        bam="results/DNA_alignment/bqsr/{sample}_NormalDNA_final_bqsr.bam",
        index="results/DNA_alignment/bqsr/{sample}_NormalDNA_final_bqsr.bam.bai"
    output:
        vcf="results/VCF/{sample}_NormalDNA_HC_raw_variants.vcf"
    conda:
        conda_env
    shell:
        """
        mkdir -p results/VCF
        gatk HaplotypeCaller -R {ref_genome} -I {input.bam} -O {output.vcf}
        """

rule gatk_mutect2:
    input:
        tumor_bam="results/DNA_alignment/bqsr/{sample}_CancerDNA_final_bqsr.bam",
        tumor_bai="results/DNA_alignment/bqsr/{sample}_CancerDNA_final_bqsr.bam.bai",
        normal_bam= "results/DNA_alignment/bqsr/{sample}_NormalDNA_final_bqsr.bam",
        normal_bai="results/DNA_alignment/bqsr/{sample}_NormalDNA_final_bqsr.bam.bai"

    output:
        vcf="results/VCF/{sample}_somatic.vcf.gz"
    params:
        tumor_sample="{sample}_CancerDNA",  # Set as parameter for better control
        normal_sample="{sample}_NormalDNA"  # Normal sample is based on the same sample name
    conda:
        conda_env
    shell:
        """
        gatk Mutect2 \
            -R {ref_genome} \
            -I {input.tumor_bam} \
            -I {input.normal_bam} \
            -normal {params.normal_sample} \
            -tumor {params.tumor_sample} \
            --germline-resource {germline_resource} \
            --panel-of-normals {panel_of_normals} \
            -O {output.vcf}
        """


