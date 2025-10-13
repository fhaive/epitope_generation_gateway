import os
import pandas as pd

# Load the CSV file from the command line config option
sample_df = pd.read_csv(config["csvfile"])

# Define the results directory, defaulting to "results" if not specified
RESULTS_DIR = config.get("results_dir", "results")
MUTATION_ANALYSIS_DIR = f"{RESULTS_DIR}/1A_mutation_analysis"

# Get the temporary directory from config, defaulting to "/tmp" if not specified
TMPDIR = config.get("tmpdir", "/tmp")
os.environ["TMPDIR"] = TMPDIR
global_tmpdir = TMPDIR

# Ensure all necessary directories exist
os.makedirs(MUTATION_ANALYSIS_DIR, exist_ok=True)

available_cores = workflow.cores

# Filter for rows where 'datatype' contains 'DNA'
dna_sample_df = sample_df[sample_df['datatype'].str.contains("DNA")]

sample_dict = {
    (row['sample_name'], row['datatype'], row['lane']): {
        'fastq_R1': f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{row['sample_name']}_{row['datatype']}_{row['lane']}_trimmed_R1.fastq.gz",
        'fastq_R2': f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{row['sample_name']}_{row['datatype']}_{row['lane']}_trimmed_R2.fastq.gz"
    }
    for idx, row in dna_sample_df.iterrows()
    if os.path.exists(f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{row['sample_name']}_{row['datatype']}_{row['lane']}_trimmed_R1.fastq.gz")
    and os.path.exists(f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{row['sample_name']}_{row['datatype']}_{row['lane']}_trimmed_R2.fastq.gz")
}


def fastq_exists_and_nonzero(path):
    return os.path.exists(path) and os.path.getsize(path) > 0

filtered_sample_dict = {}
for (sample, datatype, lane), v in sample_dict.items():
    R1 = v['fastq_R1']
    R2 = v['fastq_R2']
    if fastq_exists_and_nonzero(R1) and fastq_exists_and_nonzero(R2):
        filtered_sample_dict[(sample, datatype, lane)] = v

sample_dict = filtered_sample_dict

# Create a dictionary for sample inputs
#sample_dict = {
#    (row['sample_name'], row['datatype'], row['lane']): {
#        'fastq_R1': f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{row['sample_name']}_{row['datatype']}_{row['lane']}_trimmed_R1.fastq.gz",
#        'fastq_R2': f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{row['sample_name']}_{row['datatype']}_{row['lane']}_trimmed_R2.fastq.gz"
#    }
#    for idx, row in sample_df.iterrows()
#    if os.path.exists(f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{row['sample_name']}_{row['datatype']}_{row['lane']}_trimmed_R1.fastq.gz")
#    and os.path.exists(f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{row['sample_name']}_{row['datatype']}_{row['lane']}_trimmed_R2.fastq.gz")
#}

# Reference genome and resources
ref_genome = "resources/genome_DNA/hg38.fa"
dbsnp_vcf = "resources/dbSNP/Homo_sapiens_assembly38.dbsnp138.vcf"
germline_resource = "resources/gnomad/af-only-gnomad.hg38.vcf.gz"
panel_of_normals = "resources/PON/1000g_pon.hg38.vcf.gz"

# Helper functions
def get_lanes(sample, datatype):
    return [lane for (s, d, lane) in sample_dict.keys() if s == sample and d == datatype]



# Function to generate all output files for rule all
def generate_all_outputs():
    outputs = []
    # Include required genome resources
    outputs.append("resources/genome_DNA/hg38.dict")
    outputs.append("resources/genome_DNA/hg38.fa")
    outputs.append("resources/genome_DNA/hg38.fa.sa")
    outputs.append("resources/dbSNP/Homo_sapiens_assembly38.dbsnp138.vcf")
    outputs.append("resources/gnomad/af-only-gnomad.hg38.vcf.gz")
    outputs.append("resources/PON/1000g_pon.hg38.vcf.gz")
    outputs.append("resources/WGS_intervals/wgs_calling_regions.hg38.interval_list")
    outputs.append("resources/interval_list/exome_calling_regions.v1.1.interval_list")


    # Add sample-specific outputs
    for (sample, datatype, lane) in sample_dict.keys():
        outputs.append(f"{MUTATION_ANALYSIS_DIR}/bwa/{sample}_{datatype}_{lane}_bwa_done.txt")
#        outputs.append(f"{MUTATION_ANALYSIS_DIR}/bam/{sample}_{datatype}_{lane}.bam")
        outputs.append(f"{MUTATION_ANALYSIS_DIR}/merge/{sample}_{datatype}_merged.bam")
        outputs.append(f"{MUTATION_ANALYSIS_DIR}/dedup/{sample}_{datatype}_dedup.bam") 
        outputs.append(f"{MUTATION_ANALYSIS_DIR}/dedup/{sample}_{datatype}_recal_data.table")
        outputs.append(f"{MUTATION_ANALYSIS_DIR}/bqsr/{sample}_{datatype}_final_bqsr.bam")
        outputs.append(f"{MUTATION_ANALYSIS_DIR}/bqsr/{sample}_{datatype}_final_bqsr.bam.bai")
        outputs.append(f"{MUTATION_ANALYSIS_DIR}/metrics/{sample}_{datatype}_alignment_metrics.txt")
        outputs.append(f"{MUTATION_ANALYSIS_DIR}/metrics/{sample}_{datatype}_insert_size_metrics.txt")
        outputs.append(f"{MUTATION_ANALYSIS_DIR}/metrics/{sample}_{datatype}_insert_size_histogram.pdf")
        if datatype == "CancerDNA":
            outputs.append(f"{MUTATION_ANALYSIS_DIR}/VCF/{sample}_somatic.vcf.gz")
        outputs.append(f"{MUTATION_ANALYSIS_DIR}/VCF/{sample}_f1r2.tar.gz")
        outputs.append(f"{MUTATION_ANALYSIS_DIR}/filtering_tables/{sample}_{datatype}_pileups.table")
        outputs.append(f"{MUTATION_ANALYSIS_DIR}/filtering_tables/{sample}_contamination.table")
        outputs.append(f"{MUTATION_ANALYSIS_DIR}/filtering_tables/{sample}_tumor_segmentation.table")
        outputs.append(f"{MUTATION_ANALYSIS_DIR}/VCF/{sample}_read-orientation-model.tar.gz")
        outputs.append(f"{MUTATION_ANALYSIS_DIR}/VCF_filtered/{sample}_somatic_filtered.vcf.gz")
        if datatype == "NormalDNA":
            outputs.append(f"{MUTATION_ANALYSIS_DIR}/VCF_germline/{sample}_NormalDNA_germline_raw_variants.vcf")
            outputs.append(f"{MUTATION_ANALYSIS_DIR}/VCF_germline_filtered/{sample}_NormalDNA_indel_analysis.vcf")
            outputs.append(f"{MUTATION_ANALYSIS_DIR}/VCF_germline_filtered/{sample}_NormalDNA_snp_analysis.vcf")
            outputs.append(f"{MUTATION_ANALYSIS_DIR}/VCF_germline_filtered/{sample}_NormalDNA_germline_filtered.vcf.gz")
    return outputs

# Rule all
rule all:
    input:
        generate_all_outputs()

# Reference genome rules
rule genome_reference_generation:
    output:
        reference_genome_download="resources/genome_DNA/hg38.fa",
        reference_genome="resources/genome_DNA/hg38.dict",
        indexed="resources/genome_DNA/hg38.fa.sa"
    conda:
        "../envs/module_1A.yaml"
    shell:
        """
        mkdir -p resources/genome_DNA
        export TMPDIR={global_tmpdir}
        cd resources/genome_DNA
        wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
        gunzip hg38.fa.gz
        samtools faidx hg38.fa
        gatk CreateSequenceDictionary R=hg38.fa O=hg38.dict
        bwa index hg38.fa
        """

rule genome_reference_support_files:
    output:
        dbSNP="resources/dbSNP/Homo_sapiens_assembly38.dbsnp138.vcf",
        dbSNP_index="resources/dbSNP/Homo_sapiens_assembly38.dbsnp138.vcf.idx",  # Add index
        gnomad="resources/gnomad/af-only-gnomad.hg38.vcf.gz",
        gnomad_index="resources/gnomad/af-only-gnomad.hg38.vcf.gz.tbi",  # Add index
        PON="resources/PON/1000g_pon.hg38.vcf.gz",
        PON_index="resources/PON/1000g_pon.hg38.vcf.gz.tbi"  # Add index
    conda:
        "../envs/module_1A.yaml"
    shell:
        """
        export TMPDIR={global_tmpdir}
        mkdir -p resources/dbSNP
        cd resources/dbSNP
        wget -O Homo_sapiens_assembly38.dbsnp138.vcf https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
        gatk IndexFeatureFile -I Homo_sapiens_assembly38.dbsnp138.vcf
        cd ../..

        mkdir -p resources/gnomad
        cd resources/gnomad
        wget -O af-only-gnomad.hg38.vcf.gz https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
        tabix -p vcf af-only-gnomad.hg38.vcf.gz
        cd ../..

        mkdir -p resources/PON
        cd resources/PON
        wget -O 1000g_pon.hg38.vcf.gz https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
        tabix -p vcf 1000g_pon.hg38.vcf.gz
        """


rule download_exome_calling_regiorns:
    output:
        interval_list="resources/interval_list/exome_calling_regions.v1.1.interval_list"
    shell:
        """
        export TMPDIR={global_tmpdir}
        mkdir -p resources/interval_list
        cd resources/interval_list
        wget -O exome_calling_regions.v1.1.interval_list https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list 
        cd ../..
        """


rule download_wgs_calling_regions:
    output:
        wgs_intervals="resources/WGS_intervals/wgs_calling_regions.hg38.interval_list"
    conda:
        "../envs/module_1A.yaml"
    shell:
        """
        mkdir -p resources/WGS_intervals
        export TMPDIR={global_tmpdir}
        cd resources/WGS_intervals
        wget -O wgs_calling_regions.hg38.interval_list https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list
        """


import gzip


def make_rg(sample, datatype, lane):
    fq1 = sample_dict[(sample, datatype, lane)]["fastq_R1"]
    with gzip.open(fq1, "rt") as fh:
        line = fh.readline().strip()
    parts = line.split()
    header1 = parts[0]
    header2 = parts[1] if len(parts) > 1 else ""

    # Split by colon
    h1 = header1.lstrip("@")
    cols = h1.split(":")
    if len(cols) >= 4:
        flowcell = cols[2]
        lane_id = cols[3]
    else:
        # fallback: use the lane wildcard, and print a warning
        print(f"[WARN] FASTQ header '{line}' does not have enough fields. Falling back to lane={lane}")
        flowcell = "unknownFC"
        lane_id = lane

    # Barcode fallback
    if header2 and ":" in header2:
        barcode = header2.split(":")[-1]
    else:
        barcode = "N"

    rg_id = f"{flowcell}_{lane_id}"
    lb    = f"{flowcell}_{lane_id}"
    pu    = f"{flowcell}.{lane_id}.{barcode}"

    # Debugging info 
    print(f"RG for {sample}, {datatype}, {lane}: {rg_id} {lb} {pu}")

    return "\\t".join([
        "@RG",
        f"ID:{rg_id}",
        f"SM:{sample}_{datatype}",
        f"LB:{lb}",
        "PL:ILLUMINA",
        f"PU:{pu}"
    ])

rule bwa_mem:
    input:
        R1=lambda wc: sample_dict[(wc.sample, wc.datatype, wc.lane)]["fastq_R1"],
        R2=lambda wc: sample_dict[(wc.sample, wc.datatype, wc.lane)]["fastq_R2"],
        ref=ref_genome
    output:
        bam=temp(f"{MUTATION_ANALYSIS_DIR}/bwa/{{sample}}_{{datatype}}_{{lane}}.bam"),
        done=f"{MUTATION_ANALYSIS_DIR}/bwa/{{sample}}_{{datatype}}_{{lane}}_bwa_done.txt"
    params:
        rg=lambda wc: make_rg(wc.sample, wc.datatype, wc.lane)
    threads: 40
    conda: "../envs/module_1A.yaml"
    shell:
        r"""
        export TMPDIR={global_tmpdir}
        mkdir -p {MUTATION_ANALYSIS_DIR}
        mkdir -p {MUTATION_ANALYSIS_DIR}/bwa
        # bwa mem with a proper tab-delimited RG header:
        bwa mem -t {threads} -R "{params.rg}" \
            {input.ref} {input.R1} {input.R2} \
          | samtools sort -@ {threads} -o {output.bam}

        touch {output.done}
        """


rule merge_lane_bams:
    input:
        lambda wc: expand(f"{MUTATION_ANALYSIS_DIR}/bwa/{{sample}}_{{datatype}}_{{lane}}.bam", 
                           sample=wc.sample, datatype=wc.datatype, lane=get_lanes(wc.sample, wc.datatype))
    output:
        merged_bam=f"{MUTATION_ANALYSIS_DIR}/merge/{{sample}}_{{datatype}}_merged.bam"
    conda: "../envs/module_1A.yaml"
    shell:
        """
        export TMPDIR={global_tmpdir}
        mkdir -p {MUTATION_ANALYSIS_DIR}/merge
        if [ $(echo {input} | wc -w) -gt 1 ]; then
          inputs=$(printf 'I=%s ' {input})
          picard MergeSamFiles $inputs O={output.merged_bam} CREATE_INDEX=false
        else
          mv {input} {output.merged_bam}
        fi
        """


rule gatk_mark_duplicates:
    input:
        bam=f"{MUTATION_ANALYSIS_DIR}/merge/{{sample}}_{{datatype}}_merged.bam"
    output:
        bam=f"{MUTATION_ANALYSIS_DIR}/dedup/{{sample}}_{{datatype}}_dedup.bam"
    threads: available_cores
    conda: "../envs/module_1A.yaml"
    resources:
        markdup_slots=1
    shell:
        """
        export TMPDIR={global_tmpdir}
        mkdir -p {MUTATION_ANALYSIS_DIR}/dedup
        gatk MarkDuplicatesSpark --java-options "-Xmx256G -XX:CICompilerCount=10" \
             --conf spark.local.dir={global_tmpdir} --conf spark.driver.bindAddress=127.0.0.1 --conf spark.sql.shuffle.partitions=500 \
             -I {input.bam} -O {output.bam}
        """


rule gatk_base_recalibrator:
    input:
        bam=f"{MUTATION_ANALYSIS_DIR}/dedup/{{sample}}_{{datatype}}_dedup.bam",
        ref=ref_genome,
        known_sites=dbsnp_vcf
    output:
        recal_table=f"{MUTATION_ANALYSIS_DIR}/dedup/{{sample}}_{{datatype}}_recal_data.table"
    conda: "../envs/module_1A.yaml"
    shell:
        """
        export TMPDIR={global_tmpdir}
        gatk BaseRecalibrator -R {input.ref} -I {input.bam} --known-sites {input.known_sites} -O {output.recal_table}
        """

rule gatk_apply_bqsr:
    input:
        bam=f"{MUTATION_ANALYSIS_DIR}/dedup/{{sample}}_{{datatype}}_dedup.bam",
        recal_table=f"{MUTATION_ANALYSIS_DIR}/dedup/{{sample}}_{{datatype}}_recal_data.table",
        ref=ref_genome
    output:
        bam=f"{MUTATION_ANALYSIS_DIR}/bqsr/{{sample}}_{{datatype}}_final_bqsr.bam"
    conda: "../envs/module_1A.yaml"
    shell:
        """
        export TMPDIR={global_tmpdir}
        mkdir -p {MUTATION_ANALYSIS_DIR}/bqsr
        gatk ApplyBQSR -R {input.ref} -I {input.bam} --bqsr-recal-file {input.recal_table} -O {output.bam}
        samtools index {output.bam}
        """



rule gatk_collect_alignment_metrics:
    input:
        bam=f"{MUTATION_ANALYSIS_DIR}/bqsr/{{sample}}_{{datatype}}_final_bqsr.bam"
    output:
        metrics=f"{MUTATION_ANALYSIS_DIR}/metrics/{{sample}}_{{datatype}}_alignment_metrics.txt"
    conda:
        "../envs/module_1A.yaml"
    shell:
        """
        export TMPDIR={global_tmpdir}
        gatk CollectAlignmentSummaryMetrics -R {ref_genome} -I {input.bam} -O {output.metrics}
        """

rule gatk_collect_insert_size_metrics:
    input:
        bam=f"{MUTATION_ANALYSIS_DIR}/bqsr/{{sample}}_{{datatype}}_final_bqsr.bam"
    output:
        metrics=f"{MUTATION_ANALYSIS_DIR}/metrics/{{sample}}_{{datatype}}_insert_size_metrics.txt",
        histogram=f"{MUTATION_ANALYSIS_DIR}/metrics/{{sample}}_{{datatype}}_insert_size_histogram.pdf"
    conda:
        "../envs/module_1A.yaml"
    shell:
        """
        export TMPDIR={global_tmpdir}
        gatk CollectInsertSizeMetrics -I {input.bam} -O {output.metrics} -H {output.histogram}
        """

rule samtools_index_bam:
    input:
        bam=f"{MUTATION_ANALYSIS_DIR}/bqsr/{{sample}}_{{datatype}}_final_bqsr.bam"
    output:
        index=f"{MUTATION_ANALYSIS_DIR}/bqsr/{{sample}}_{{datatype}}_final_bqsr.bam.bai"
    conda:
        "../envs/module_1A.yaml"
    shell:
        """
        export TMPDIR={global_tmpdir}
        samtools index {input.bam}
        """

# Rule for HaplotypeCaller on NormalDNA only
rule gatk_haplotypecaller:
    input:
        bam=f"{MUTATION_ANALYSIS_DIR}/bqsr/{{sample}}_NormalDNA_final_bqsr.bam",
        index=f"{MUTATION_ANALYSIS_DIR}/bqsr/{{sample}}_NormalDNA_final_bqsr.bam.bai"
    output:
        vcf=f"{MUTATION_ANALYSIS_DIR}/VCF_germline/{{sample}}_NormalDNA_germline_raw_variants.vcf"
    conda:
        "../envs/module_1A.yaml"
    shell:
        """
        mkdir -p {MUTATION_ANALYSIS_DIR}/VCF_germline
        export TMPDIR={global_tmpdir}

        gatk HaplotypeCaller -R {ref_genome} -I {input.bam} -O {output.vcf}
        """


rule gatk_hoplotypecaller_filtering:
    input:
        vcf=f"{MUTATION_ANALYSIS_DIR}/VCF_germline/{{sample}}_NormalDNA_germline_raw_variants.vcf",
        ref="resources/genome_DNA/hg38.fa"
    output:
        snp_raw_vcf=f"{MUTATION_ANALYSIS_DIR}/VCF_germline/{{sample}}_NormalDNA_raw_snps.vcf",
        indel_raw_vcf=f"{MUTATION_ANALYSIS_DIR}/VCF_germline/{{sample}}_NormalDNA_raw_indel.vcf",
        snp_filtered_vcf=f"{MUTATION_ANALYSIS_DIR}/VCF_germline_filtered/{{sample}}_NormalDNA_snp_filtered.vcf",
        indel_filtered_vcf=f"{MUTATION_ANALYSIS_DIR}/VCF_germline_filtered/{{sample}}_NormalDNA_indel_filtered.vcf",
        temp_snp_analysis=f"{MUTATION_ANALYSIS_DIR}/VCF_germline_filtered/{{sample}}_NormalDNA_snp_analysis.vcf",
        temp_indel_analysis=f"{MUTATION_ANALYSIS_DIR}/VCF_germline_filtered/{{sample}}_NormalDNA_indel_analysis.vcf"
    conda:
        "../envs/module_1A.yaml"
    shell:
        """
        export TMPDIR={global_tmpdir}

        gatk SelectVariants -R {input.ref} -V {input.vcf} --select-type SNP -O {output.snp_raw_vcf}
        gatk SelectVariants -R {input.ref} -V {input.vcf} --select-type INDEL -O {output.indel_raw_vcf}

        gatk VariantFiltration \
            -R {input.ref} \
            -V {output.snp_raw_vcf} \
            -O {output.snp_filtered_vcf} \
            -filter-name "QD_filter" -filter "QD < 2.0" \
            -filter-name "FS_filter" -filter "FS > 60.0" \
            -filter-name "MQ_filter" -filter "MQ < 40.0" \
            -filter-name "SOR_filter" -filter "SOR > 4.0" \
            -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
            -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
            -genotype-filter-expression "DP < 10" \
            -genotype-filter-name "DP_filter" \
            -genotype-filter-expression "GQ < 10" \
            -genotype-filter-name "GQ_filter"

        gatk VariantFiltration \
            -R {input.ref} \
            -V {output.indel_raw_vcf} \
            -O {output.indel_filtered_vcf} \
            -filter-name "QD_filter" -filter "QD < 2.0" \
            -filter-name "FS_filter" -filter "FS > 200.0" \
            -filter-name "SOR_filter" -filter "SOR > 10.0" \
            -genotype-filter-expression "DP < 10" \
            -genotype-filter-name "DP_filter" \
            -genotype-filter-expression "GQ < 10" \
            -genotype-filter-name "GQ_filter"

        gatk SelectVariants \
            --exclude-filtered \
            -V {output.snp_filtered_vcf} \
            -O {output.temp_snp_analysis}


        gatk SelectVariants \
            --exclude-filtered \
            -V {output.indel_filtered_vcf} \
            -O {output.temp_indel_analysis}

        """

rule merge_filtered_germline:
    input:
        temp_snp_analysis=f"{MUTATION_ANALYSIS_DIR}/VCF_germline_filtered/{{sample}}_NormalDNA_snp_analysis.vcf",
        temp_indel_analysis=f"{MUTATION_ANALYSIS_DIR}/VCF_germline_filtered/{{sample}}_NormalDNA_indel_analysis.vcf"
    output:
        sorted_snp_analysis=f"{MUTATION_ANALYSIS_DIR}/VCF_germline_filtered/{{sample}}_NormalDNA_snp_sorted.vcf",
        sorted_indel_analysis=f"{MUTATION_ANALYSIS_DIR}/VCF_germline_filtered/{{sample}}_NormalDNA_indel_sorted.vcf",
        compressed_snp_analysis=f"{MUTATION_ANALYSIS_DIR}/VCF_germline_filtered/{{sample}}_NormalDNA_snp_sorted.vcf.gz",
        compressed_indel_analysis=f"{MUTATION_ANALYSIS_DIR}/VCF_germline_filtered/{{sample}}_NormalDNA_indel_sorted.vcf.gz",
        merged_filtered_vcf=f"{MUTATION_ANALYSIS_DIR}/VCF_germline_filtered/{{sample}}_NormalDNA_germline_filtered.vcf.gz"
    conda:
        "../envs/module_1A.yaml"
    shell:
        """
        export TMPDIR={global_tmpdir}

        bcftools sort -O v -o {output.sorted_snp_analysis} {input.temp_snp_analysis}

        bcftools sort -O v -o {output.sorted_indel_analysis} {input.temp_indel_analysis}

        bgzip -c {output.sorted_snp_analysis} > {output.compressed_snp_analysis}
        bgzip -c {output.sorted_indel_analysis} > {output.compressed_indel_analysis}

        # Index the compressed files
        tabix -p vcf {output.compressed_snp_analysis}
        tabix -p vcf {output.compressed_indel_analysis}

        bcftools concat -a \
            {output.compressed_snp_analysis} \
            {output.compressed_indel_analysis} \
            -o {output.merged_filtered_vcf} \
            -O z

        """




rule gatk_mutect2:
    input:
        tumor_bam=f"{MUTATION_ANALYSIS_DIR}/bqsr/{{sample}}_CancerDNA_final_bqsr.bam",
        tumor_bai=f"{MUTATION_ANALYSIS_DIR}/bqsr/{{sample}}_CancerDNA_final_bqsr.bam.bai",
        normal_bam=f"{MUTATION_ANALYSIS_DIR}/bqsr/{{sample}}_NormalDNA_final_bqsr.bam",
        normal_bai=f"{MUTATION_ANALYSIS_DIR}/bqsr/{{sample}}_NormalDNA_final_bqsr.bam.bai"
    output:
        vcf=f"{MUTATION_ANALYSIS_DIR}/VCF/{{sample}}_somatic.vcf.gz",
        f1r2=f"{MUTATION_ANALYSIS_DIR}/VCF/{{sample}}_f1r2.tar.gz"
    params:
        tumor_sample="{sample}_CancerDNA",
        normal_sample="{sample}_NormalDNA",
        intervals=(
            "" if config.get("mode", "all_intervals") == "all_intervals" else
            f"-L {config['wgs_intervals']}" if config.get("mode") == "WGS" else
            f"-L {config['bedfile']}" if config.get("mode") == "WES" and "bedfile" in config else
            ""
        )

    conda:
        "../envs/module_1A.yaml"
    shell:
        """
        mkdir -p {MUTATION_ANALYSIS_DIR}/VCF
        export TMPDIR={global_tmpdir}
        gatk Mutect2 -R {ref_genome} \
            -I {input.tumor_bam} \
            -I {input.normal_bam} \
            -normal {params.normal_sample} \
            -tumor {params.tumor_sample} \
            --germline-resource {germline_resource} \
            --panel-of-normals {panel_of_normals} \
            {params.intervals} \
            --f1r2-tar-gz {output.f1r2} \
            -O {output.vcf}
        """

#rule pile_up_summaries:
#    input:
#        bam=f"{MUTATION_ANALYSIS_DIR}/bqsr/{{sample}}_{{datatype}}_final_bqsr.bam",
#        index=f"{MUTATION_ANALYSIS_DIR}/bqsr/{{sample}}_{{datatype}}_final_bqsr.bam.bai",
#        interval_list="resources/interval_list/exome_calling_regions.v1.1.interval_list"
#    output:
#        table=f"{MUTATION_ANALYSIS_DIR}/filtering_tables/{{sample}}_{{datatype}}_pileups.table"
#    conda:
#        "../envs/module_1A.yaml"
#    threads: available_cores
#    shell:
#        """
#        gatk  GetPileupSummaries \
#            --java-options "-Xmx64G" --tmp-dir {global_tmpdir} \
#            -I {input.bam} \
#            -V {germline_resource} \
#            -L {input.interval_list} \
#            -O {output.table}
#
#        """


rule pile_up_summaries:
    input:
        bam   = f"{MUTATION_ANALYSIS_DIR}/bqsr/{{sample}}_{{datatype}}_final_bqsr.bam",
        index = f"{MUTATION_ANALYSIS_DIR}/bqsr/{{sample}}_{{datatype}}_final_bqsr.bam.bai"
    output:
        table = f"{MUTATION_ANALYSIS_DIR}/filtering_tables/{{sample}}_{{datatype}}_pileups.table"
    params:
        intervals = (
            "" if config.get("mode", "all_intervals") == "all_intervals" else
            f"-L {config['wgs_intervals']}" if config.get("mode") == "WGS" else
            f"-L {config['bedfile']}" if (config.get("mode") == "WES" and "bedfile" in config) else
            ""
        )
    conda:
        "../envs/module_1A.yaml"
    threads: available_cores
    shell:
        r"""
        gatk GetPileupSummaries \
          --java-options "-Xmx64G" --tmp-dir {global_tmpdir} \
          -I {input.bam} \
          -V {germline_resource} \
          {params.intervals} \
          -O {output.table}
        """


rule gatk_contamination:
    input:
        tumor_table=f"{MUTATION_ANALYSIS_DIR}/filtering_tables/{{sample}}_CancerDNA_pileups.table",
        normal_table=f"{MUTATION_ANALYSIS_DIR}/filtering_tables/{{sample}}_NormalDNA_pileups.table"
    output:
        contamination_table=f"{MUTATION_ANALYSIS_DIR}/filtering_tables/{{sample}}_contamination.table",
        tumor_segmentation=f"{MUTATION_ANALYSIS_DIR}/filtering_tables/{{sample}}_tumor_segmentation.table"
    conda:
        "../envs/module_1A.yaml"
    shell:
        """
        export TMPDIR={global_tmpdir}
        gatk CalculateContamination \
            -I {input.tumor_table} \
            -matched {input.normal_table} \
            -O {output.contamination_table} \
            --tumor-segmentation {output.tumor_segmentation}
        """
rule read_orientation_estimation:
    input:
        f1r2=f"{MUTATION_ANALYSIS_DIR}/VCF/{{sample}}_f1r2.tar.gz"
    output:
        read_oriantation_model= f"{MUTATION_ANALYSIS_DIR}/VCF/{{sample}}_read-orientation-model.tar.gz"
    conda:
        "../envs/module_1A.yaml"
    shell:
        """
        export TMPDIR={global_tmpdir}
        gatk LearnReadOrientationModel \
            -I {input.f1r2} \
            -O {output.read_oriantation_model}
        """

rule gatk_filter_mutect_calls:
    input:
        vcf=f"{MUTATION_ANALYSIS_DIR}/VCF/{{sample}}_somatic.vcf.gz",
        contamination_table=f"{MUTATION_ANALYSIS_DIR}/filtering_tables/{{sample}}_contamination.table",
        tumor_segmentation=f"{MUTATION_ANALYSIS_DIR}/filtering_tables/{{sample}}_tumor_segmentation.table",
        read_oriantation_model= f"{MUTATION_ANALYSIS_DIR}/VCF/{{sample}}_read-orientation-model.tar.gz"
    output:
        filtered_vcf=f"{MUTATION_ANALYSIS_DIR}/VCF_filtered/{{sample}}_somatic_filtered.vcf.gz"
    params:
        intervals=(
            "" if config.get("mode", "all_intervals") == "all_intervals" else
            f"-L {config['wgs_intervals']}" if config.get("mode") == "WGS" else
            f"-L {config['bedfile']}" if config.get("mode") == "WES" and "bedfile" in config else
            ""
        )
    conda:
        "../envs/module_1A.yaml"
    threads: available_cores
    shell:
        """
        export TMPDIR={global_tmpdir}

        gatk FilterMutectCalls \
            -R {ref_genome} \
            -V {input.vcf} \
            --contamination-table {input.contamination_table} \
            --tumor-segmentation {input.tumor_segmentation} \
            --ob-priors {input.read_oriantation_model} \
            {params.intervals} \
            -O {output.filtered_vcf}
        """
