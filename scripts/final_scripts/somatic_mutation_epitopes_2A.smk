import os
import subprocess
import pandas as pd


# Define the results directory, defaulting to "results" if not specified
RESULTS_DIR = config.get("results_dir", "results")
SOMATIC_EPITOPES_DIR = f"{RESULTS_DIR}/2A_somatic_mutation_epitopes"
HLA_ANALYSIS_DIR = f"{RESULTS_DIR}/1B_RNA_fusion_HLA"
MUTATION_ANALYSIS_DIR = f"{RESULTS_DIR}/1A_mutation_analysis"

# Get the temporary directory from config, defaulting to "/tmp" if not specified
TMPDIR = config.get("tmpdir", "/tmp")
os.environ["TMPDIR"] = TMPDIR
global_tmpdir = TMPDIR

# Ensure all necessary directories exist
os.makedirs(SOMATIC_EPITOPES_DIR, exist_ok=True)

conda_env = "../envs/module_2A.yaml"
conda_env_second = "../envs/module_2A_second_part.yaml"
conda_env_Picard = "../envs/module_2A_Picard.yaml"
conda_env_HLA = "../envs/module_2B.yaml"
# Load config file
configfile: "scripts/final_scripts/config/somatic_mutation_epitopes_config_2A.yaml"

pvacseq_options = config.get("pvacseq_options", {})
pvacseq_threads = pvacseq_options.get("threads", 20)  # Get threads from YAML
pvacseq_depth = pvacseq_options.get("depth", 1000)
epitope_lengths_1 = ",".join(map(str, pvacseq_options.get("epitope_lengths_1", [8, 9, 10, 11, 12, 13, 14, 15])))
epitope_lengths_2 = ",".join(map(str, pvacseq_options.get("epitope_lengths_2", [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30])))
pvacseq_algorithms =  pvacseq_options.get("algorithm", "all")

# Function to check if the Docker image exists
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


# check pvactools docker

if not check_docker_image("griffithlab/pvactools"):
    print("Docker image 'griffithlab/pvactools' not found. Building the image...")
    build_docker_image("griffithlab/pvactools")

# Load the CSV file from the command line config option
sample_df = pd.read_csv(config["csvfile"])
conda_env = "../envs/module_2A.yaml"

# Check if the "HLA-types" column exists
hla_column_exists = "HLA-types" in sample_df.columns

# Create a dictionary for sample inputs
sample_dict = {
    (row['sample_name'], row['datatype'], row['lane']): {
        'fastq_R1': row['fastq_R1'],
        'fastq_R2': row['fastq_R2'],
        'datatype': row['datatype'],
        'hla_types': row['HLA-types'] if hla_column_exists and not pd.isna(row['HLA-types']) else None
    }
    for idx, row in sample_df.iterrows()
}


def generate_all_outputs():
    outputs = []
    outputs.append(f"resources/VEP/setup_done.flag")
    outputs.append(f"resources/VEP/VEP_plugins/vep_plugins_installed.flag")
    outputs.append(f"resources/Kallisto_ref/Homo_sapiens.GRCh38.cdna.all.fa.gz")
    outputs.append(f"resources/Kallisto_ref/Homo_sapiens.GRCh38.cdna.all.idx")
    for (sample, datatype, lane) in sample_dict.keys():
        if datatype == "CancerDNA":
            outputs.append(f"{SOMATIC_EPITOPES_DIR}/annotated_somatic_VCF/{sample}_somatic_VEP.vcf")
        if datatype == "CancerRNA":
            if not hla_column_exists:
                outputs.append(f"{HLA_ANALYSIS_DIR}/ArcasHLA/{sample}_CancerRNA_{lane}_HLAs_pvacFormat.txt")
            outputs.append(f"{SOMATIC_EPITOPES_DIR}/kallisto_quantification/{sample}/abundance.tsv")
            outputs.append(f"{SOMATIC_EPITOPES_DIR}/kallisto_somatic_VCF/{sample}_somatic_quantification.vcf")
            outputs.append(f"{SOMATIC_EPITOPES_DIR}/tumor_only_VCF/{sample}_tumor_only.vcf")
            outputs.append(f"{SOMATIC_EPITOPES_DIR}/tumor_only_VCF/{sample}_tumor_only.vcf.gz.tbi")
        if datatype == "NormalDNA":
            outputs.append(f"{SOMATIC_EPITOPES_DIR}/annotated_germline_VCF/{sample}_germline_VEP.vcf")
            outputs.append(f"{SOMATIC_EPITOPES_DIR}/annotated_germline_VCF_name_updated/{sample}_germline_VEP_name_updated.vcf")
            outputs.append(f"{SOMATIC_EPITOPES_DIR}/combined_VCF/{sample}_combined_somatic_plus_germline.vcf")
            outputs.append(f"{SOMATIC_EPITOPES_DIR}/combined_VCF/{sample}_combined.sorted.vcf")
            outputs.append(f"{SOMATIC_EPITOPES_DIR}/phased_VCF/{sample}_phased.vcf")
            outputs.append(f"{SOMATIC_EPITOPES_DIR}/annotated_phased_VCF/{sample}_phased_VEP.vcf")
            outputs.append(f"{SOMATIC_EPITOPES_DIR}/annotated_phased_VCF/{sample}_phased_VEP.vcf.gz.tbi")
            outputs.append(f"{SOMATIC_EPITOPES_DIR}/kallisto_somatic_VCF/{sample}_somatic_quantification.vcf.gz.tbi")
            outputs.append(f"{SOMATIC_EPITOPES_DIR}/pvacSeq/{sample}/")
    return outputs

rule all:
    input:
        generate_all_outputs()



rule setup_vep_with_conda:
    conda:
        conda_env
    output:
        touch("resources/VEP/setup_done.flag")
    params:
        vep_cache_url="https://ftp.ensembl.org/pub/release-103/variation/indexed_vep_cache/homo_sapiens_vep_103_GRCh38.tar.gz",
        plugins_repo="https://github.com/Ensembl/VEP_plugins.git"
    shell:
        """
        mkdir -p resources/VEP &&
        cd resources/VEP &&
        
        # Download and extract the VEP cache
        #curl -O {params.vep_cache_url} &&
        #tar xzf homo_sapiens_vep_103_GRCh38.tar.gz &&
        
        # Clone the VEP plugins repository
        if [ ! -d VEP_plugins ]; then
            git clone -b release/103 {params.plugins_repo}
        fi &&
        cd ../..
        # Mark setup as done
        touch {output}
        """


rule install_vep_and_plugins_in_docker:
    output:
        touch("resources/VEP/VEP_plugins/vep_plugins_installed.flag")
    params:
        docker_host_flag=lambda wildcards: f"-H {config['dockerdir']}" if config.get("dockerdir") else "",
        vep_plugins_dir="/data/resources/VEP/VEP_plugins"
    shell:
        """
    # Run the command inside the Docker container
        docker {params.docker_host_flag} run -u $(id -u):$(id -g) -v `pwd`:/data -v {global_tmpdir}:{global_tmpdir} -e TMPDIR={global_tmpdir} --rm griffithlab/pvactools \
            /bin/bash -c "
            pvacseq install_vep_plugin {params.vep_plugins_dir}
            "

        # Create the flag file to indicate success
        touch {output}
        """

rule kallisto_prep:
    output:
        cdna="resources/Kallisto_ref/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    shell:
        """
        mkdir -p resources/Kallisto_ref

        wget -O {output.cdna} http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
        """


rule somatic_vep:
    input:
        vcf=f"{MUTATION_ANALYSIS_DIR}/VCF_filtered/{{sample}}_somatic_filtered.vcf.gz",
        fasta="resources/genome_DNA/hg38.fa",
        vep_cache_ready="resources/VEP/setup_done.flag",
        vep_plugins_installed="resources/VEP/VEP_plugins/vep_plugins_installed.flag"
    output:
        vep_output=f"{SOMATIC_EPITOPES_DIR}/annotated_somatic_VCF/{{sample}}_somatic_VEP.vcf",
    params:
        cache_dir="resources/VEP",  # Directory for the VEP cache inside the Docker container
        plugins_dir="resources/VEP/VEP_plugins",  # Directory for VEP plugins inside the Docker container
    conda:
        conda_env
    shell:
        """
        mkdir -p {SOMATIC_EPITOPES_DIR}/annotated_somatic_VCF

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
            --offline \
            --cache \
            --dir_cache {params.cache_dir} \
            --plugin Frameshift \
            --plugin Wildtype \
            --dir_plugins {params.plugins_dir} \
            --pick \
            --transcript_version
         """


rule germline_vep:
    input:
        vcf=f"{MUTATION_ANALYSIS_DIR}/VCF_germline_filtered/{{sample}}_NormalDNA_germline_filtered.vcf.gz",
        fasta="resources/genome_DNA/hg38.fa",
        vep_cache_ready="resources/VEP/setup_done.flag",
        vep_plugins_installed="resources/VEP/VEP_plugins/vep_plugins_installed.flag"
    output:
        vep_output=f"{SOMATIC_EPITOPES_DIR}/annotated_germline_VCF/{{sample}}_germline_VEP.vcf",
    params:
        cache_dir="resources/VEP",  # Directory for the VEP cache inside the Docker container
        plugins_dir="resources/VEP/VEP_plugins",  # Directory for VEP plugins inside the Docker container
    conda:
        conda_env
    shell:
        """
        mkdir -p {SOMATIC_EPITOPES_DIR}/annotated_germline_VCF

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
            --offline \
            --cache \
            --dir_cache {params.cache_dir} \
            --plugin Frameshift \
            --plugin Wildtype \
            --dir_plugins {params.plugins_dir} \
            --pick \
            --transcript_version
         """


rule kallisto_index:
    input:
        transcriptome="resources/Kallisto_ref/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    output:
        index="resources/Kallisto_ref/Homo_sapiens.GRCh38.cdna.all.idx"
    conda:
        conda_env
    shell:
        """
        mkdir -p resources/Kallisto_ref
        kallisto index -i {output.index} {input.transcriptome}
        """

def get_lanes(sample):
    # Filter the dataframe for the sample name and the CancerRNA datatype
    lanes = sample_df[
        (sample_df['sample_name'] == sample) & 
        (sample_df['datatype'] == 'CancerRNA')
    ]['lane'].tolist()
    
    return lanes

rule kallisto_quant:
    input:
        index="resources/Kallisto_ref/Homo_sapiens.GRCh38.cdna.all.idx",

        R1=lambda wildcards: expand(
            f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{{sample}}_CancerRNA_{{lane}}_trimmed_R1.fastq.gz",
            sample=wildcards.sample,
            lane=get_lanes(wildcards.sample)
        ),
        R2=lambda wildcards: expand(
            f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{{sample}}_CancerRNA_{{lane}}_trimmed_R2.fastq.gz",
            sample=wildcards.sample,
            lane=get_lanes(wildcards.sample)
        )
    output:
        quant=f"{SOMATIC_EPITOPES_DIR}/kallisto_quantification/{{sample}}/abundance.tsv"
    conda:
        conda_env
    shell:
        """
        mkdir -p {SOMATIC_EPITOPES_DIR}/kallisto_quantification

        # Aggregate inputs across lanes
        R1_inputs=$(echo {input.R1} | tr ' ' ',')
        R2_inputs=$(echo {input.R2} | tr ' ' ',')

        kallisto quant -i {input.index} -o $(dirname {output.quant}) $R1_inputs $R2_inputs
        """

rule annotate_vcf_with_expression:
    input:
        vcf=f"{SOMATIC_EPITOPES_DIR}/annotated_somatic_VCF/{{sample}}_somatic_VEP.vcf",
        quant=f"{SOMATIC_EPITOPES_DIR}/kallisto_quantification/{{sample}}/abundance.tsv"
    output:
        annotated_vcf=f"{SOMATIC_EPITOPES_DIR}/kallisto_somatic_VCF/{{sample}}_somatic_quantification.vcf"
    params:
        sample=lambda wildcards: f"{wildcards.sample}_CancerDNA"
    conda:
        conda_env
    shell:
        """
        mkdir -p {SOMATIC_EPITOPES_DIR}/kallisto_somatic_VCF
        vcf-expression-annotator {input.vcf} {input.quant} -s {params.sample}  kallisto transcript -o {output.annotated_vcf}
        """




rule tumor_only_vcf:
    input:
        annotated_vcf=f"{SOMATIC_EPITOPES_DIR}/kallisto_somatic_VCF/{{sample}}_somatic_quantification.vcf"
    output:
        tumor_vcf=f"{SOMATIC_EPITOPES_DIR}/tumor_only_VCF/{{sample}}_tumor_only.vcf",
        tumor_vcf_gz=f"{SOMATIC_EPITOPES_DIR}/tumor_only_VCF/{{sample}}_tumor_only.vcf.gz",
        tumor_vcf_tbi=f"{SOMATIC_EPITOPES_DIR}/tumor_only_VCF/{{sample}}_tumor_only.vcf.gz.tbi"
    params:
        sample=lambda wildcards: f"{wildcards.sample}_CancerDNA"
    conda: 
       conda_env_second
    shell:
        """
        mkdir -p {SOMATIC_EPITOPES_DIR}/tumor_only_VCF/
#        gatk_jar=$(find $(conda info --base) -name GenomeAnalysisTK.jar)
        gatk_jar=$(find $CONDA_PREFIX -name GenomeAnalysisTK.jar)
        if [ -z "$gatk_jar" ]; then
            echo "Error: GenomeAnalysisTK.jar not found in the conda environment." >&2
            exit 1
        fi

        # Run GATK with the dynamically resolved JAR path
        java -Xmx32G -jar $gatk_jar \
        -T SelectVariants \
        -R resources/genome_DNA/hg38.fa \
        --variant {input.annotated_vcf} \
        -sn {params.sample} \
        -o {output.tumor_vcf}

        bgzip -c {output.tumor_vcf} > {output.tumor_vcf_gz}
        tabix -p vcf {output.tumor_vcf_gz}

        """

rule update_germline_sample_name:
    input:
        germline_vcf=f"{SOMATIC_EPITOPES_DIR}/annotated_germline_VCF/{{sample}}_germline_VEP.vcf"
    output:
        updated_vcf=f"{SOMATIC_EPITOPES_DIR}/annotated_germline_VCF_name_updated/{{sample}}_germline_VEP_name_updated.vcf"
    params:
        old_name=lambda wildcards: f"{wildcards.sample}_NormalDNA",
        new_name=lambda wildcards: f"{wildcards.sample}_CancerDNA"
    conda:
        conda_env
    shell:
        """
        mkdir -p {SOMATIC_EPITOPES_DIR}//annotated_germline_VCF_name_updated
        bcftools reheader -s <(echo "{params.old_name} {params.new_name}") {input.germline_vcf} -o {output.updated_vcf}
        """

rule combine_variants:
    input:
        germline_vcf_updated=f"{SOMATIC_EPITOPES_DIR}/annotated_germline_VCF_name_updated/{{sample}}_germline_VEP_name_updated.vcf",
        tumor_vcf=f"{SOMATIC_EPITOPES_DIR}/tumor_only_VCF/{{sample}}_tumor_only.vcf.gz"
    output:
        germline_vcf_updated_gz=f"{SOMATIC_EPITOPES_DIR}/annotated_germline_VCF_name_updated/{{sample}}_germline_VEP_name_updated.vcf.gz",
        combined_vcf=f"{SOMATIC_EPITOPES_DIR}/combined_VCF/{{sample}}_combined_somatic_plus_germline.vcf"
    conda:
        conda_env_second
    shell:
        """
        bgzip -c {input.germline_vcf_updated} > {output.germline_vcf_updated_gz}
        tabix -p vcf {output.germline_vcf_updated_gz}

        # Create output directory
        mkdir -p {SOMATIC_EPITOPES_DIR}/combined_VCF/

        # Dynamically find the GATK JAR path
        #gatk_jar=$(find $(conda info --base) -name GenomeAnalysisTK.jar)  #this version works if server is using conda connected to conda base
        gatk_jar=$(find $CONDA_PREFIX -name GenomeAnalysisTK.jar)

        # Check if the JAR file was found
        if [ -z "$gatk_jar" ]; then
            echo "Error: GenomeAnalysisTK.jar not found in the conda environment." >&2
            exit 1
        fi

        # Run GATK CombineVariants
        java -Xmx16G -jar $gatk_jar \
            -T CombineVariants \
            -R resources/genome_DNA/hg38.fa \
            --variant {output.germline_vcf_updated_gz} \
            --variant {input.tumor_vcf} \
            -o {output.combined_vcf} \
            --assumeIdenticalSamples
        """

rule sort_combined_vcf:
    input:
        combined_vcf=f"{SOMATIC_EPITOPES_DIR}/combined_VCF/{{sample}}_combined_somatic_plus_germline.vcf",
        reference_dict="resources/genome_DNA/hg38.dict"
    output:
        sorted_vcf=f"{SOMATIC_EPITOPES_DIR}/combined_VCF/{{sample}}_combined.sorted.vcf"
    conda:
        conda_env_Picard
    shell:
        """
        mkdir -p {SOMATIC_EPITOPES_DIR}/sorted_VCF
        #picard_jar=$(find $(conda info --base) -name picard.jar | head -n 1)
        picard_jar=$(find $CONDA_PREFIX -name picard.jar)

        if [ -z "$picard_jar" ]; then
            echo "Error: picard.jar not found in the conda environment." >&2
            exit 1
        fi

        java -Xmx16G -jar $picard_jar SortVcf \
            I={input.combined_vcf} \
            O={output.sorted_vcf} \
            SEQUENCE_DICTIONARY={input.reference_dict}
        """

rule phase_variants:
    input:
        sorted_vcf=f"{SOMATIC_EPITOPES_DIR}/combined_VCF/{{sample}}_combined.sorted.vcf",
        tumor_bam=f"{MUTATION_ANALYSIS_DIR}/bqsr/{{sample}}_CancerDNA_final_bqsr.bam",
        reference_fasta="resources/genome_DNA/hg38.fa"
    output:
        phased_vcf=f"{SOMATIC_EPITOPES_DIR}/phased_VCF/{{sample}}_phased.vcf"
    conda:
        conda_env_second
    shell:
        """
        mkdir -p {SOMATIC_EPITOPES_DIR}/phased_VCF
        #gatk_jar=$(find $(conda info --base) -name GenomeAnalysisTK.jar)
        gatk_jar=$(find $CONDA_PREFIX -name GenomeAnalysisTK.jar)
        if [ -z "$gatk_jar" ]; then
            echo "Error: GenomeAnalysisTK.jar not found in the conda environment." >&2
            exit 1
        fi

        java -Xmx16G -jar $gatk_jar \
            -T ReadBackedPhasing \
            -R {input.reference_fasta} \
            -I {input.tumor_bam} \
            --variant {input.sorted_vcf} \
            -L {input.sorted_vcf} \
            -o {output.phased_vcf}
        """


rule phased_vep:
    input:
        phased_vcf=f"{SOMATIC_EPITOPES_DIR}/phased_VCF/{{sample}}_phased.vcf",
        fasta="resources/genome_DNA/hg38.fa"
    output:
        vep_output=f"{SOMATIC_EPITOPES_DIR}/annotated_phased_VCF/{{sample}}_phased_VEP.vcf"
    params:
        cache_dir="resources/VEP",
        plugins_dir="resources/VEP/VEP_plugins"
    conda:
        conda_env
    shell:
        """
        mkdir -p {SOMATIC_EPITOPES_DIR}/annotated_phased_VCF

        vep \
            --input_file {input.phased_vcf} \
            --output_file {output.vep_output} \
            --format vcf \
            --vcf \
            --symbol \
            --terms SO \
            --tsl \
            --biotype \
            --hgvs \
            --fasta {input.fasta} \
            --offline \
            --cache \
            --dir_cache {params.cache_dir} \
            --plugin Frameshift \
            --plugin Wildtype \
            --dir_plugins {params.plugins_dir} \
            --pick \
            --transcript_version
         """

rule bgzip_and_index:
    input:
        phased_vep=f"{SOMATIC_EPITOPES_DIR}/annotated_phased_VCF/{{sample}}_phased_VEP.vcf",
        main_vep=f"{SOMATIC_EPITOPES_DIR}/kallisto_somatic_VCF/{{sample}}_somatic_quantification.vcf"
    output:
        phased_vep_gz=f"{SOMATIC_EPITOPES_DIR}/annotated_phased_VCF/{{sample}}_phased_VEP.vcf.gz",
        main_vep_gz=f"{SOMATIC_EPITOPES_DIR}/kallisto_somatic_VCF/{{sample}}_somatic_quantification.vcf.gz",
        phased_vep_gz_tbi=f"{SOMATIC_EPITOPES_DIR}/annotated_phased_VCF/{{sample}}_phased_VEP.vcf.gz.tbi",
        main_vep_gz_tbi=f"{SOMATIC_EPITOPES_DIR}/kallisto_somatic_VCF/{{sample}}_somatic_quantification.vcf.gz.tbi"
    conda:
        conda_env_second
    shell:
        """
        bgzip -c {input.phased_vep} > {output.phased_vep_gz}
        tabix -p vcf {output.phased_vep_gz}

        bgzip -c {input.main_vep} > {output.main_vep_gz}
        tabix -p vcf {output.main_vep_gz}
        """




if not hla_column_exists:
    rule extract_hla_types:
        input:
            f"{HLA_ANALYSIS_DIR}/ArcasHLA/{{sample}}_CancerRNA_{{lane}}_sorted.genotype.json"
        output:
            f"{HLA_ANALYSIS_DIR}/ArcasHLA/{{sample}}_CancerRNA_{{lane}}_HLAs_pvacFormat.txt"
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

BLAST_DB_DIR = config.get("blast_db_dir", "resources/blastdb/refseq_select_prot")

rule setup_blastdb:
    """
    Install BLAST+ via conda and download the RefSeq-Select protein DB.
    Creates a sentinel file so the step runs only once.
    """
    output:
        db_dir = directory(BLAST_DB_DIR),
        ready  = touch(f"{BLAST_DB_DIR}/DB_READY")
    threads: 4
    conda:  "../envs/blast.yaml"
    shell:
        r"""
        # Make sure target directory exists, then move into it
        mkdir -p {output.db_dir}
        cd {output.db_dir}

        # Download DB into this directory
        export BLASTDB=$(pwd)
        update_blastdb.pl --decompress refseq_select_prot

        # Copy blastp and the NCBI shared libs into the DB dir (portable bundle)
        mkdir -p bin lib
        BLAST_BIN=$(dirname $(which blastp))
        BLAST_LIB=$BLAST_BIN/../lib
        cp $BLAST_BIN/blastp bin/
        cp $BLAST_LIB/libncbi* lib/ 2>/dev/null || true

        # Mark everything as complete
        touch DB_READY
        """
 
if not hla_column_exists:
    rule PvacSeq_HLA:
        input:
            main_vcf=f"{SOMATIC_EPITOPES_DIR}/kallisto_somatic_VCF/{{sample}}_somatic_quantification.vcf.gz",
            hla_types=lambda wildcards: expand(
            f"{HLA_ANALYSIS_DIR}/ArcasHLA/{{sample}}_CancerRNA_{{lane}}_HLAs_pvacFormat.txt",
            sample=wildcards.sample,
            lane=get_lanes(wildcards.sample)
        ),
            phased_vcf=f"{SOMATIC_EPITOPES_DIR}/annotated_phased_VCF/{{sample}}_phased_VEP.vcf.gz",
            blastdb    = rules.setup_blastdb.output.ready
        output:
            directory(f"{SOMATIC_EPITOPES_DIR}/pvacSeq/{{sample}}/")
        threads:
            pvacseq_threads
        params:
            docker_host_flag=lambda wildcards: f"-H {config['dockerdir']}" if config.get("dockerdir") else "",
            depth=pvacseq_depth,
            epitope_lengths_1=epitope_lengths_1,
            epitope_lengths_2=epitope_lengths_2,
            algorithms=pvacseq_algorithms,
            blastdb_dir = os.path.abspath(BLAST_DB_DIR)
        shell:
            """
            docker {params.docker_host_flag} run -u $(id -u):$(id -g) -v `pwd`:/data -v {global_tmpdir}:/tmp -v {params.blastdb_dir}:/blastdb -e BLASTDB=/blastdb -e TMPDIR=/tmp --rm griffithlab/pvactools \
                /bin/bash -c "
                pvacseq run \
/data/{input.main_vcf} \
{wildcards.sample}_CancerDNA \
\$(cat /data/{input.hla_types}) \
MHCnuggetsI MHCnuggetsII  \
/data/{output} \
-e1 {params.epitope_lengths_1} \
-e2 {params.epitope_lengths_2} \
--iedb-install-directory /opt/iedb \
-t 10 \
--normal-sample-name {wildcards.sample}_NormalDNA \
-d 1000 \
-s 100 \
--phased-proximal-variants-vcf /data/{input.phased_vcf} "
        """

#--run-reference-proteome-similarity  \
#--blastp-path /blastdb/bin/blastp \
#--blastp-db refseq_select_prot "

#            docker {params.docker_host_flag} run -u $(id -u):$(id -g) -v `pwd`:/data -v {global_tmpdir}:/tmp -e TMPDIR=/tmp --rm \
#  --cpus 80 \
#  --cpuset-cpus 0-79 \
#  -e OMP_NUM_THREADS=4 \
#  -e MKL_NUM_THREADS=4 \
#  -e OPENBLAS_NUM_THREADS=4 \
#  griffithlab/pvactools \
#default for -s is 100 this is the length of tsv files are split into. -s 1000 could make things faster but uses more ram
if hla_column_exists:
    rule PvacSeq_UserHLA:
        input:
            star_fusion_tsv=f"{HLA_ANALYSIS_DIR}/StarFusionOut/{{sample}}_CancerRNA_{{lane}}/star-fusion.fusion_predictions.tsv",
            ag_fusion_directory=f"{FUSION_EPITOPES_DIR}/AGfusion/{{sample}}_CancerRNA_{{lane}}/"
        output:
            directory(f"{FUSION_EPITOPES_DIR}/pvacFuse/{{sample}}_CancerDNA_{{lane}}/")
        threads:
            pvacfuse_threads
        resources:
            pvacslot = 1
        params:
            docker_host_flag=lambda wildcards: f"-H {config['dockerdir']}" if config.get("dockerdir") else "",
            hla_types=lambda wildcards: sample_dict[(wildcards.sample, "CancerDNA", wildcards.lane)]["hla_types"],
            depth=pvacfuse_depth,
            epitope_lengths_1=epitope_lengths_1,
            epitope_lengths_2=epitope_lengths_2,
            algorithms=pvacfuse_algorithms
        shell:
            """
            docker run -u $(id -u):$(id -g) -v `pwd`:/data -v {global_tmpdir}:{global_tmpdir} -e TMPDIR={global_tmpdir} --rm griffithlab/pvactools \
                /bin/bash -c "
                pvacfuse run \
/data/{input.ag_fusion_directory} \
{wildcards.sample}_CancerRNA_{wildcards.lane} \
{params.hla_types} \
{params.algorithms} \
/data/{output} \
-d {params.depth} \
-t {threads} \
--starfusion-file /data/{input.star_fusion_tsv} \
--iedb-install-directory /opt/iedb \
-e1 {params.epitope_lengths_1} \
-e2 {params.epitope_lengths_2} "
        """







