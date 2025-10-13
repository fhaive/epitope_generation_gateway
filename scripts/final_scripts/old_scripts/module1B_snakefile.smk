import os
import subprocess
import pandas as pd


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
            check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error building Docker image: {e}")
        raise

# Check if the Docker image exists, and build if it does not
if not check_docker_image("jfx319/arcashla"):
    print("Docker image 'jfx319/arcashla' not found. Building the image...")
    build_docker_image("jfx319/arcashla")

if not check_docker_image("trinityctat/starfusion"):
    print("Docker image 'trinityctat/starfusion' not found. Building the image...")
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
    for (sample, datatype, lane) in sample_dict.keys():
        if datatype == "CancerRNA":  # Only include outputs for CancerRNA
            outputs.append(f"results/StarFusionOut/{sample}_CancerRNA_{lane}/")
            outputs.append(f"results/StarFusionOut/{sample}_CancerRNA_{lane}/{sample}_CancerRNA_{lane}.bam")
            outputs.append(f"results/sorted_bam/{sample}_CancerRNA_{lane}_sorted.bam")
            outputs.append(f"results/ArcasHLA/{sample}_CancerRNA_{lane}_sorted.partial_genotype.json")
            outputs.append(f"results/RNA_Counts/{sample}_CancerRNA_{lane}_counts.txt")
            #outputs.append(f"results/ArcasHLA/run.genotypes.tsv")
            #outputs.append(f"results/ArcasHLA/run.partial_genotypes.tsv")
            #outputs.append(f"results/ArcasHLA/{sample}_CancerRNA_{lane}_sorted.partial_genotype.json")
    return outputs

# Rule all to list all output files
rule all:
    input:
        generate_all_outputs()

# Rule to run STAR-Fusion only for CancerRNA
rule star_fusion:
    input:
        R1=lambda wildcards: f"results/trimmed_fastq/{wildcards.sample}_CancerRNA_{wildcards.lane}_trimmed_R1.fastq.gz",
        R2=lambda wildcards: f"results/trimmed_fastq/{wildcards.sample}_CancerRNA_{wildcards.lane}_trimmed_R2.fastq.gz",
        genome_lib_dir="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir"
    output:
        directory=directory("results/StarFusionOut/{sample}_CancerRNA_{lane}/"),
        bam_renamed="results/StarFusionOut/{sample}_CancerRNA_{lane}/{sample}_CancerRNA_{lane}.bam"
    shell:
        """
        docker run -u $(id -u):$(id -g) -v `pwd`:/data --rm trinityctat/starfusion \
            /bin/bash -c "
                STAR-Fusion \
                --left_fq /data/{input.R1} \
                --right_fq /data/{input.R2} \
                --genome_lib_dir /data/{input.genome_lib_dir} \
                -O /data/{output.directory}
                #chmod -R 777 /data/{output.directory}
            "
        mv {output.directory}/Aligned.out.bam {output.bam_renamed}
        """


rule samtools:
    input:
        bam_renamed="results/StarFusionOut/{sample}_CancerRNA_{lane}/{sample}_CancerRNA_{lane}.bam"
    output:
        sorted_bam="results/sorted_bam/{sample}_CancerRNA_{lane}_sorted.bam"
    conda:
        conda_env

    shell:
        """
        samtools sort -o {output.sorted_bam} {input.bam_renamed}
        """

rule arcasHLA:
    input:
        sorted_bam="results/sorted_bam/{sample}_CancerRNA_{lane}_sorted.bam"
    output:
        #HLA_typing = "results/ArcasHLA/run.genotypes.tsv",
        #partial_HLA_typing = "results/ArcasHLA/run.partial_genotypes.tsv"
        partial_genotyping = "results/ArcasHLA/{sample}_CancerRNA_{lane}_sorted.partial_genotype.json"
        #merged_genotyping= "results/ArcasHLA/run.genotypes.tsv",
        #merged_partial_genotyping="results/ArcasHLA/run.partial_genotypes.tsv"

    shell:
        """
        mkdir -p results/ArcasHLA && \
        docker run -u $(id -u):$(id -g) -v `pwd`:/data --rm jfx319/arcashla \
            /bin/bash -c "
                chmod -R 777 /data/results/ArcasHLA
                arcasHLA extract /data/{input.sorted_bam} -o /data/results/ArcasHLA && \
                arcasHLA genotype /data/results/ArcasHLA/{wildcards.sample}_CancerRNA_{wildcards.lane}_sorted.extracted.fq.gz -o /data/results/ArcasHLA && \
                arcasHLA partial -G /data/results/ArcasHLA/{wildcards.sample}_CancerRNA_{wildcards.lane}_sorted.genotype.json /data/results/ArcasHLA/{wildcards.sample}_CancerRNA_{wildcards.lane}_sorted.extracted.fq.gz -o /data/results/ArcasHLA 
                #arcasHLA merge --i /data/results/ArcasHLA
                #chmod -R 777 /data/results/ArcasHLA
            "
        """

rule featureCounts:
    input:
        bam_renamed="results/StarFusionOut/{sample}_CancerRNA_{lane}/{sample}_CancerRNA_{lane}.bam"
    output:
        raw_count="results/RNA_Counts/{sample}_CancerRNA_{lane}_counts.txt"
    conda:
        conda_env

    shell:
        """
        gunzip -c resources/GTF/Homo_sapiens.GRCh38.103.gtf.gz > resources/GTF/Homo_sapiens.GRCh38.103.gtf
        featureCounts -p -a resources/GTF/Homo_sapiens.GRCh38.103.gtf -o {output.raw_count} {input.bam_renamed}        

        """
