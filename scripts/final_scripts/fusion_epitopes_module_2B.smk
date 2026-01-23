import os
import subprocess
import pandas as pd


# Define the results directory, defaulting to "results" if not specified
RESULTS_DIR = config.get("results_dir", "results")
FUSION_EPITOPES_DIR = f"{RESULTS_DIR}/2B_fusion_epitopes"
RNA_ANALYSIS_DIR = f"{RESULTS_DIR}/1B_RNA_fusion_HLA"


# Get the temporary directory from config, defaulting to "/tmp" if not specified
TMPDIR = config.get("tmpdir", "/tmp")
os.environ["TMPDIR"] = TMPDIR
global_tmpdir = TMPDIR

# Ensure all necessary directories exist
os.makedirs(FUSION_EPITOPES_DIR, exist_ok=True)

# Load config file
configfile: "scripts/final_scripts/config/fusion_epitopes_config_2B.yaml"

pvacfuse_options = config.get("pvacfuse_options", {})
pvacfuse_threads = pvacfuse_options.get("threads", 20)  # Get threads from YAML
pvacfuse_depth = pvacfuse_options.get("depth", 1000)
epitope_lengths_1 = ",".join(map(str, pvacfuse_options.get("epitope_lengths_1", [8, 9, 10, 11, 12, 13, 14, 15])))
epitope_lengths_2 = ",".join(map(str, pvacfuse_options.get("epitope_lengths_2", [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30])))
pvacfuse_algorithms =  pvacfuse_options.get("algorithm", "all")
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
conda_env = "../envs/module_2B.yaml"

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
    for (sample, datatype, lane) in sample_dict.keys():
        if datatype == "CancerRNA":  # Only include outputs for CancerRNA
            outputs.append(f"resources/agfusion/AGFusion_database/agfusion.homo_sapiens.103.db")
            outputs.append(f"{FUSION_EPITOPES_DIR}/AGfusion/{sample}_CancerRNA_{lane}/")
            if not hla_column_exists: 
                outputs.append(f"{RNA_ANALYSIS_DIR}/ArcasHLA/{sample}_CancerRNA_{lane}_HLAs_pvacFormat.txt")
            outputs.append(f"{FUSION_EPITOPES_DIR}/pvacFuse/{sample}_CancerRNA_{lane}/")
    return outputs


rule all:
    input:
        generate_all_outputs()


rule setup_env:
    conda:
        conda_env
    output:
        "resources/agfusion/AGFusion_database/agfusion.homo_sapiens.103.db"
    shell:
        """
        # Run pyensembl installation
        mkdir -p resources/agfusion/AGFusion_database
        cd resources/agfusion/AGFusion_database
        wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
        gunzip Pfam-A.clans.tsv.gz
        pyensembl install --release 103 --species homo_sapiens
        agfusion build -d . -s homo_sapiens -r 103 --pfam Pfam-A.clans.tsv              
        """


rule AGFusion:
    input:
        star_fusion_tsv=f"{RNA_ANALYSIS_DIR}/StarFusionOut/{{sample}}_CancerRNA_{{lane}}/star-fusion.fusion_predictions.tsv",
        AGFusionDatabase="resources/agfusion/AGFusion_database/agfusion.homo_sapiens.103.db"
    output:
        directory(f"{FUSION_EPITOPES_DIR}/AGfusion/{{sample}}_CancerRNA_{{lane}}/")
    conda:
        conda_env
    shell:
        """
        agfusion batch \
            -f {input.star_fusion_tsv} \
            -a starfusion \
            -db {input.AGFusionDatabase} \
            -o {output} \
            --middlestar \
            --noncanonical
        """
if not hla_column_exists:
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



if not hla_column_exists:
    rule PvacFuse_HLA:
        input:
            star_fusion_tsv=f"{RNA_ANALYSIS_DIR}/StarFusionOut/{{sample}}_CancerRNA_{{lane}}/star-fusion.fusion_predictions.tsv",
            hla_types=f"{RNA_ANALYSIS_DIR}/ArcasHLA/{{sample}}_CancerRNA_{{lane}}_HLAs_pvacFormat.txt",
            ag_fusion_directory=f"{FUSION_EPITOPES_DIR}/AGfusion/{{sample}}_CancerRNA_{{lane}}/"
        output:
            directory(f"{FUSION_EPITOPES_DIR}/pvacFuse/{{sample}}_CancerRNA_{{lane}}/")
        threads:
            pvacfuse_threads
        params:
            docker_host_flag=lambda wildcards: f"-H {config['dockerdir']}" if config.get("dockerdir") else "",
            depth=pvacfuse_depth,
            epitope_lengths_1=epitope_lengths_1,
            epitope_lengths_2=epitope_lengths_2,
            algorithms=pvacfuse_algorithms
        shell:
            """
            docker {params.docker_host_flag} run -u $(id -u):$(id -g) -v `pwd`:/data -v {global_tmpdir}:{global_tmpdir} -e TMPDIR={global_tmpdir} --rm griffithlab/pvactools \
                /bin/bash -c "
                pvacfuse run \
/data/{input.ag_fusion_directory} \
{wildcards.sample}_CancerRNA_{wildcards.lane} \
\$(cat /data/{input.hla_types}) \
MHCnuggetsI MHCnuggetsII \
/data/{output} \
-d {params.depth} \
-t {threads} \
--starfusion-file /data/{input.star_fusion_tsv} \
--iedb-install-directory /opt/iedb \
-e1 {params.epitope_lengths_1} \
-e2 {params.epitope_lengths_2} "
        """


if hla_column_exists:
    rule PvacFuse_UserHLA:
        input:
            star_fusion_tsv=f"{RNA_ANALYSIS_DIR}/StarFusionOut/{{sample}}_CancerRNA_{{lane}}/star-fusion.fusion_predictions.tsv",
            ag_fusion_directory=f"{FUSION_EPITOPES_DIR}/AGfusion/{{sample}}_CancerRNA_{{lane}}/"
        output:
            directory(f"{FUSION_EPITOPES_DIR}/pvacFuse/{{sample}}_CancerRNA_{{lane}}/")
        threads:
            pvacfuse_threads
        params:
            docker_host_flag=lambda wildcards: f"-H {config['dockerdir']}" if config.get("dockerdir") else "",
            hla_types=lambda wildcards: sample_dict[(wildcards.sample, "CancerRNA", wildcards.lane)]["hla_types"],
            depth=pvacfuse_depth,
            epitope_lengths_1=epitope_lengths_1,
            epitope_lengths_2=epitope_lengths_2,
            algorithms=pvacfuse_algorithms
        shell:
            """
            docker {params.docker_host_flag} run -u $(id -u):$(id -g) -v `pwd`:/data -v {global_tmpdir}:{global_tmpdir} -e TMPDIR={global_tmpdir} --rm griffithlab/pvactools \
                /bin/bash -c "
                pvacfuse run \
/data/{input.ag_fusion_directory} \
{wildcards.sample}_CancerRNA_{wildcards.lane} \
{params.hla_types} \
MHCnuggetsI MHCnuggetsII \
/data/{output} \
-d {params.depth} \
-t {threads} \
--starfusion-file /data/{input.star_fusion_tsv} \
--iedb-install-directory /opt/iedb \
-e1 {params.epitope_lengths_1} \
-e2 {params.epitope_lengths_2} "
        """

