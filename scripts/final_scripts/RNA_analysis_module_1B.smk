import os
import subprocess
import pandas as pd


# Define the results directory, defaulting to "results" if not specified
RESULTS_DIR = config.get("results_dir", "results")
RNA_ANALYSIS_DIR = f"{RESULTS_DIR}/1B_RNA_fusion_HLA"

# Check if HLA Typing is Enabled, if --config HLA="no_HLA" then no HLA done if its "yes_HLA" or its missing then its carried out
HLA_ENABLED = config.get("HLA", "yes_HLA") != "no_HLA"

# Optional docker installation directory
docker_host = config.get("dockerdir", "")
docker_host_flag = ["-H", docker_host] if docker_host else []


def check_docker_image(image_name):
    try:
        result = subprocess.run(
            ["docker"] + docker_host_flag + ["inspect", image_name],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        return result.returncode == 0
    except Exception as e:
        print(f"Error checking Docker image: {e}")
        return False

def build_docker_image(image_name):
    try:
        subprocess.run(
            ["docker"] + docker_host_flag + ["pull", image_name],
            check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error pulling Docker image: {e}")
        raise

# Ensure Docker images are available
if not check_docker_image("hubentu/arcas-hla"):
    print("Docker image 'hubentu/arcas-hla' not found. Pulling the image...")
    build_docker_image("hubentu/arcas-hla")

if not check_docker_image("trinityctat/starfusion"):
    print("Docker image 'trinityctat/starfusion' not found. Pulling the image...")
    build_docker_image("trinityctat/starfusion")


# Load the CSV file from the command line config option
sample_df = pd.read_csv(config["csvfile"])

conda_env = "../envs/module_1B.yaml"

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
    outputs.append("resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir")
#    outputs.append("resources/GTF/Homo_sapiens.GRCh38.103.gtf")
    for (sample, datatype, lane) in sample_dict.keys():
        if datatype == "CancerRNA":  # Only include outputs for CancerRNA


            outputs.append(f"{RNA_ANALYSIS_DIR}/StarFusionOut/{sample}_CancerRNA_{lane}/{sample}_CancerRNA_{lane}.bam")
            outputs.append(f"{RNA_ANALYSIS_DIR}/sorted_bam/{sample}_CancerRNA_{lane}_sorted.bam")
            outputs.append(f"{RNA_ANALYSIS_DIR}/RNA_Counts/{sample}_CancerRNA_{lane}_counts.txt")
            
            if HLA_ENABLED:
                outputs.append(f"{RNA_ANALYSIS_DIR}/ArcasHLA/{sample}_CancerRNA_{lane}_sorted.genotype.json")
                outputs.append(f"{RNA_ANALYSIS_DIR}/ArcasHLA/{sample}_CancerRNA_{lane}_HLAs_pvacFormat.txt")
    return outputs

# Rule all to list all output files
rule all:
    input:
        generate_all_outputs()



rule download_fusion_reference:
    output:
#        reference_genome_download="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play",
        genome_lib_dir=directory("resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir")
    shell:
        """
        mkdir -p resources/genome_RNA_fusion
        
        cd resources/genome_RNA_fusion

        wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz
        tar -xvzf GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz
        cd ../..
        """
#rule download_RNA_GTF:
#    output:
#        reference_GTF_download="resources/GTF/Homo_sapiens.GRCh38.103.gtf"
#    shell:
#        """
#        mkdir -p resources/GTF
#        cd resources/GTF
#        wget ftp://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz
#        gunzip -c Homo_sapiens.GRCh38.103.gtf.gz > Homo_sapiens.GRCh38.103.gtf
#        cd ../..
#        """


# Rule to run STAR-Fusion only for CancerRNA
rule star_fusion:
    input:
        R1=lambda wildcards: f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{wildcards.sample}_CancerRNA_{wildcards.lane}_trimmed_R1.fastq.gz",
        R2=lambda wildcards: f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{wildcards.sample}_CancerRNA_{wildcards.lane}_trimmed_R2.fastq.gz",
        genome_lib_dir="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir"
    output:
        bam_renamed=f"{RNA_ANALYSIS_DIR}/StarFusionOut/{{sample}}_CancerRNA_{{lane}}/{{sample}}_CancerRNA_{{lane}}.bam",        
        sample_dir = directory(f"{RNA_ANALYSIS_DIR}/StarFusionOut/{{sample}}_CancerRNA_{{lane}}/")
    params:
        docker_host_flag=lambda wildcards: f"-H {config['dockerdir']}" if config.get("dockerdir") else ""    
    shell:
        """
        docker {params.docker_host_flag} run -u $(id -u):$(id -g) -v `pwd`:/data --rm trinityctat/starfusion \
            /bin/bash -c "
                STAR-Fusion \
                --left_fq /data/{input.R1} \
                --right_fq /data/{input.R2} \
                --genome_lib_dir /data/{input.genome_lib_dir} \
                -O /data/{output.sample_dir}
            "
        mv {output.sample_dir}/Aligned.out.bam {output.bam_renamed}
        """


rule samtools:
    input:
        bam_renamed=f"{RNA_ANALYSIS_DIR}/StarFusionOut/{{sample}}_CancerRNA_{{lane}}/{{sample}}_CancerRNA_{{lane}}.bam"
    output:
        sorted_bam=f"{RNA_ANALYSIS_DIR}/sorted_bam/{{sample}}_CancerRNA_{{lane}}_sorted.bam"
    conda:
        conda_env

    shell:
        """
        samtools sort -o {output.sorted_bam} {input.bam_renamed}
        """

if HLA_ENABLED:
    rule arcasHLA:
        input:
            sorted_bam=f"{RNA_ANALYSIS_DIR}/sorted_bam/{{sample}}_CancerRNA_{{lane}}_sorted.bam"
        output:
            HLA_genotyping = f"{RNA_ANALYSIS_DIR}/ArcasHLA/{{sample}}_CancerRNA_{{lane}}_sorted.genotype.json"
        params:
            docker_host_flag=lambda wildcards: f"-H {config['dockerdir']}" if config.get("dockerdir") else ""
        shell:
            """
            mkdir -p {RNA_ANALYSIS_DIR}/ArcasHLA/  && \
            docker {params.docker_host_flag} run -u $(id -u):$(id -g) -v `pwd`:/data --rm hubentu/arcas-hla \
                /bin/bash -c "
                    arcasHLA extract --unmapped /data/{input.sorted_bam} -o /data/{RNA_ANALYSIS_DIR}/ArcasHLA/ && \
                    arcasHLA genotype /data/{RNA_ANALYSIS_DIR}/ArcasHLA/{wildcards.sample}_CancerRNA_{wildcards.lane}_sorted.extracted.1.fq.gz /data/{RNA_ANALYSIS_DIR}/ArcasHLA/{wildcards.sample}_CancerRNA_{wildcards.lane}_sorted.extracted.2.fq.gz -o /data/{RNA_ANALYSIS_DIR}/ArcasHLA
                "
            """


    rule extract_hla_types:
        input:
            f"{RNA_ANALYSIS_DIR}/ArcasHLA/{{sample}}_CancerRNA_{{lane}}_sorted.genotype.json"
        output:
            f"{RNA_ANALYSIS_DIR}/ArcasHLA/{{sample}}_CancerRNA_{{lane}}_HLAs_pvacFormat.txt"
        conda:
            conda_env
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



rule featureCounts:
    input:
        bam_renamed=f"{RNA_ANALYSIS_DIR}/StarFusionOut/{{sample}}_CancerRNA_{{lane}}/{{sample}}_CancerRNA_{{lane}}.bam"
    output:
        raw_count=f"{RNA_ANALYSIS_DIR}/RNA_Counts/{{sample}}_CancerRNA_{{lane}}_counts.txt"
    conda:
        conda_env

    shell:
        """
        
        featureCounts -p -a resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf -o {output.raw_count} {input.bam_renamed}

        """
