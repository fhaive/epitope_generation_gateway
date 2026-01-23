import os
import subprocess
import pandas as pd


# Define the results directory, defaulting to "results" if not specified
RESULTS_DIR = config.get("results_dir", "results")
SPLICING_EPITOPES_DIR = f"{RESULTS_DIR}/2C_splicing_epitopes"
RNA_ANALYSIS_DIR = f"{RESULTS_DIR}/1B_RNA_fusion_HLA"
SOMATIC_EPITOPES_DIR = f"{RESULTS_DIR}/2A_somatic_mutation_epitopes"

HLA_ANALYSIS_DIR = f"{RESULTS_DIR}/1B_RNA_fusion_HLA"
MUTATION_ANALYSIS_DIR= f"{RESULTS_DIR}/1A_mutation_analysis"
# Get the temporary directory from config, defaulting to "/tmp" if not specified
TMPDIR = config.get("tmpdir", "/tmp")
os.environ["TMPDIR"] = TMPDIR
global_tmpdir = TMPDIR

# Ensure all necessary directories exist
os.makedirs(SPLICING_EPITOPES_DIR, exist_ok=True)

# Load config file
configfile: "scripts/final_scripts/config/splicing_epitopes_config_2C.yaml"

pvacsplice_options = config.get("pvacsplice_options", {})
pvacsplice_threads = pvacsplice_options.get("threads", 20)  # Get threads from YAML
pvacsplice_depth = pvacsplice_options.get("depth", 1000)
epitope_lengths_1 = ",".join(map(str, pvacsplice_options.get("epitope_lengths_1", [8, 9, 10, 11, 12, 13, 14, 15])))
epitope_lengths_2 = ",".join(map(str, pvacsplice_options.get("epitope_lengths_2", [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30])))
pvacsplice_algorithms =  pvacsplice_options.get("algorithm", "all")
# Function to check if the Docker image exists
def check_docker_image(image_name):
    try:
        # Use subprocess to check if the image exists
        result = subprocess.run(
            ["docker", "inspect", image_name],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        return result.returncode == 0  # If the return code is 0, the image exists
    except Exception as e:
        print(f"Error checking Docker image: {e}")
        return False

# Function to build the Docker image
def build_docker_image(image_name):
    try:
        subprocess.run(
            ["docker", "pull",  image_name],
            check=True,
            env={**os.environ, "DOCKER_TMPDIR": os.path.abspath(global_tmpdir)}
        )
    except subprocess.CalledProcessError as e:
        print(f"Error building Docker image: {e}")
        raise

# check pvactools docker
if not check_docker_image("griffithlab/pvactools"):
    print("Docker image 'griffithlab/pvactools' not found. Building the image...")
    build_docker_image("griffithlab/pvactools")

# Load the CSV file from the command line config option
sample_df = pd.read_csv(config["csvfile"])
conda_env = "../envs/module_2C.yaml"
conda_env_HLA = "../envs/module_2B.yaml"
conda_env_VEP = "../envs/module_2C_VEP.yaml"

# Check if the "HLA-types" column exists
hla_column_exists = "HLA-types" in sample_df.columns

# Filter for rows where 'datatype' contains 'RNA'
RNA_sample_df = sample_df[sample_df['datatype'].str.contains("RNA")]

# Create a dictionary for sample inputs
sample_dict = {
    (row['sample_name'], row['datatype'], row['lane']): {
        'fastq_R1': row['fastq_R1'],
        'fastq_R2': row['fastq_R2'],
        'datatype': row['datatype'],
        'hla_types': row['HLA-types'] if hla_column_exists and not pd.isna(row['HLA-types']) else None
    }
    for idx, row in RNA_sample_df.iterrows()
}


def generate_all_outputs():
    outputs = []
    outputs.append(f"resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf.gz.tbi")
    for (sample, datatype, lane) in sample_dict.keys():
        if datatype == "CancerRNA":  # Only include outputs for CancerRNA
            outputs.append(f"{SPLICING_EPITOPES_DIR}/annotated_somatic_VCF/{sample}_somatic_VEP.vcf")
            outputs.append(f"{SPLICING_EPITOPES_DIR}/regtools_genomic_VCF_genecode/{sample}_splice_effects.tsv")
            outputs.append(f"{SPLICING_EPITOPES_DIR}/pvacSplice/{sample}/")
    return outputs


rule all:
    input:
        generate_all_outputs()

def get_lanes(sample):
    # Function to retrieve all lanes for a given sample
    return RNA_sample_df[sample_df['sample_name'] == sample]['lane'].tolist()



rule bgzip_and_tabix_gtf:
    input:
        gtf="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf"
    output:
        gtf_gz="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf.gz",
        tbi="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf.gz.tbi"
    conda:
        conda_env_VEP
    shell:
        """
        sort -k1,1 -k4,4n {input} > temp.sorted.gtf
        bgzip -c temp.sorted.gtf > {output.gtf_gz}
        tabix -p gff {output.gtf_gz}
        rm temp.sorted.gtf
        """


rule somatic_vep_BAM_matching_fasta:
    input:
        vcf=f"{MUTATION_ANALYSIS_DIR}/VCF_filtered/{{sample}}_somatic_filtered.vcf.gz",
        fasta="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa",
        gtf="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf.gz"
        #gtf="resources/updated_gtf/ref_annot_updated.gtf.gz"
    output:
        vep_output=f"{SPLICING_EPITOPES_DIR}/annotated_somatic_VCF/{{sample}}_somatic_VEP.vcf"
    params:
        #cache_dir="resources/VEP",  # Directory for the VEP cache inside the Docker container
        plugins_dir="resources/VEP/VEP_plugins"  # Directory for VEP plugins inside the Docker container
    conda:
        conda_env_VEP
    shell:
        """
        mkdir -p {SPLICING_EPITOPES_DIR}/annotated_somatic_VCF
        vep \
          --input_file {input.vcf} \
          --output_file {output.vep_output} \
          --format vcf \
          --vcf \
          --symbol \
          --terms SO \
          --tsl \
          --biotype \
          --hgvs \
          --fasta {input.fasta} \
          --gtf {input.gtf} \
          --no_stats \
          --gene_version \
          --transcript_version \
          --pick
         """



rule bgzip_and_index_VEP:
    input:
        vcf=f"{SPLICING_EPITOPES_DIR}/annotated_somatic_VCF/{{sample}}_somatic_VEP.vcf"
    output:
        vcf_gz=f"{SPLICING_EPITOPES_DIR}/annotated_somatic_VCF/{{sample}}_somatic_VEP.vcf.gz",
        vcf_tbi=f"{SPLICING_EPITOPES_DIR}/annotated_somatic_VCF/{{sample}}_somatic_VEP.vcf.gz.tbi"
    conda:
        conda_env_VEP
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf_gz}
        tabix -p vcf {output.vcf_gz}
        """



rule regtools_splice_effects_genomic_VCF_genecode:
    input:
        vcf=f"{SPLICING_EPITOPES_DIR}/annotated_somatic_VCF/{{sample}}_somatic_VEP.vcf",
        bam=lambda wildcards: expand(
            f"{RNA_ANALYSIS_DIR}/sorted_bam/{{sample}}_CancerRNA_{{lane}}_sorted.bam",
            sample=wildcards.sample,
            lane=get_lanes(wildcards.sample)
        ),
        ref="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa",
        gtf="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf"
    output:
        tsv=f"{SPLICING_EPITOPES_DIR}/regtools_genomic_VCF_genecode/{{sample}}_splice_effects.tsv"
    conda:
        conda_env
    shell:
        """
        mkdir -p {SOMATIC_EPITOPES_DIR}/regtools

        regtools cis-splice-effects identify \
            -o {output.tsv} \
            -s XS \
            {input.vcf} \
            {input.bam} \
            {input.ref} \
            {input.gtf}
        """



if not hla_column_exists:
    rule extract_hla_types:
        input:
            f"{RNA_ANALYSIS_DIR}/ArcasHLA/{{sample}}_CancerRNA_{{lane}}_sorted.genotype.json"
        output:
            f"{RNA_ANALYSIS_DIR}/ArcasHLA/{{sample}}_CancerRNA_{{lane}}_HLAs_pvacFormat.txt"
        conda:
            conda_env_HLA
        shell:
            """
            HLA_TYPES=$(jq -r 'to_entries | .[] | "\(.value[])"' {input} | awk '
            BEGIN {{OFS=","}}
            {{
                if ($1 ~ /^D/) {{  # Do not add HLA- for entries starting with D
                    split($1, parts, ":")
                    type = parts[1]
                    if (length(parts) > 1) {{
                        type = type ":" parts[2]
                    }}
                }} else {{  # Add HLA- prefix for other entries
                    split($1, parts, ":")
                    type = parts[1]
                    if (length(parts) > 1) {{
                        type = type ":" parts[2]
                    }}
                    type = "HLA-" type
                }}
                # Append type to list of types
                types = (NR == 1 ? type : types OFS type)
            }}
            END {{print types}}
            ')
            echo $HLA_TYPES > {output}
            """





if not hla_column_exists:
    rule PvacSplice_HLA:
        input:
            main_vcf=f"{SPLICING_EPITOPES_DIR}/annotated_somatic_VCF/{{sample}}_somatic_VEP.vcf.gz",
            hla_types=lambda wildcards: expand(
            f"{HLA_ANALYSIS_DIR}/ArcasHLA/{{sample}}_CancerRNA_{{lane}}_HLAs_pvacFormat.txt",
            sample=wildcards.sample,
            lane=get_lanes(wildcards.sample)
        ),
            regtools_tsv=f"{SPLICING_EPITOPES_DIR}/regtools_genomic_VCF_genecode/{{sample}}_splice_effects.tsv",
            ref="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa",
            gtf="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf.gz"
        output:
            directory(f"{SPLICING_EPITOPES_DIR}/pvacSplice/{{sample}}/")
        threads:
            pvacsplice_threads
        params:
            epitope_lengths_1=epitope_lengths_1,
            epitope_lengths_2=epitope_lengths_2,
            algorithms=pvacsplice_algorithms
        shell:
            """
            docker run -u $(id -u):$(id -g) -v `pwd`:/data -v {global_tmpdir}:{global_tmpdir} -e TMPDIR={global_tmpdir} --rm griffithlab/pvactools \
                /bin/bash -c "
                pvacsplice run \
/data/{input.regtools_tsv} \
{wildcards.sample}_CancerDNA \
\$(cat /data/{input.hla_types}) \
MHCnuggetsI MHCnuggetsII \
/data/{output} \
/data/{input.main_vcf} \
/data/{input.ref} \
/data/{input.gtf} \
--iedb-install-directory /opt/iedb \
-t {threads} \
-e1 {params.epitope_lengths_1} \
-e2 {params.epitope_lengths_2} \
--normal-sample-name  {wildcards.sample}_NormalDNA
"
        """

if hla_column_exists:
    rule PvacSplice_UserHLA:
        input:
            main_vcf=f"{SPLICING_EPITOPES_DIR}/annotated_somatic_VCF/{{sample}}_somatic_VEP.vcf.gz",
            regtools_tsv=f"{SPLICING_EPITOPES_DIR}/regtools_genomic_VCF_genecode/{{sample}}_splice_effects.tsv",
            ref="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa",
            gtf="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf.gz"
        output:
            directory(f"{SPLICING_EPITOPES_DIR}/pvacSplice/{{sample}}/")
        threads:
            pvacsplice_threads
        params:
            hla_types=lambda wildcards: sample_dict[(wildcards.sample, "CancerDNA", wildcards.lane)]["hla_types"],
            epitope_lengths_1=epitope_lengths_1,
            epitope_lengths_2=epitope_lengths_2,
            algorithms=pvacsplice_algorithms
        shell:
            """
            docker run -u $(id -u):$(id -g) -v `pwd`:/data -v {global_tmpdir}:{global_tmpdir} -e TMPDIR={global_tmpdir} --rm griffithlab/pvactools \
                /bin/bash -c "
                pvacsplice run \
/data/{input.regtools_tsv} \
{wildcards.sample}_CancerDNA \
{params.hla_types} \
MHCnuggetsI MHCnuggetsII \
/data/{output} \
/data/{input.main_vcf} \
/data/{input.ref} \
/data/{input.gtf} \
--iedb-install-directory /opt/iedb \
-t {threads} \
-e1 {params.epitope_lengths_1} \
-e2 {params.epitope_lengths_2} \
--normal-sample-name {wildcards.sample}_NormalDNA
"
        """

