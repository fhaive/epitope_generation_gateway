import os
import subprocess
import pandas as pd


# Define the results directory, defaulting to "results" if not specified
RESULTS_DIR = config.get("results_dir", "results")
SPLICING_EPITOPES_DIR = f"{RESULTS_DIR}/2C_splicing_epitopes"
RNA_ANALYSIS_DIR = f"{RESULTS_DIR}/1B_RNA_fusion_HLA"
SOMATIC_EPITOPES_DIR = f"{RESULTS_DIR}/2A_somatic_mutation_epitopes"
ALT_SPLICING_DIR = f"{RESULTS_DIR}/1C_alternative_splicing"
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
    outputs.append("resources/updated_gtf/ref_annot_updated.gtf")
    for (sample, datatype, lane) in sample_dict.keys():
        if datatype == "CancerRNA":  # Only include outputs for CancerRNA
            outputs.append(f"{SPLICING_EPITOPES_DIR}/hisat2/{sample}_CancerRNA_{lane}_sorted.bam")
            outputs.append(f"{SPLICING_EPITOPES_DIR}/annotated_somatic_VCF_ensembl/{sample}_somatic_VEP.vcf")
            outputs.append(f"{SPLICING_EPITOPES_DIR}/regtools_genomic_VCF_ensembl/{sample}_splice_effects.tsv")
            outputs.append(f"{SPLICING_EPITOPES_DIR}/regtools/{sample}_splice_effects.tsv")
#            outputs.append(f"{SPLICING_EPITOPES_DIR}/annotated_somatic_VCF/{sample}_somatic_VEP.vcf")
            outputs.append(f"{SPLICING_EPITOPES_DIR}/regtools_genomic_VCF/{sample}_splice_effects.tsv")
            outputs.append(f"{SPLICING_EPITOPES_DIR}/regtools_genomic_VCF_genecode/{sample}_splice_effects.tsv")
#            outputs.append(f"{SPLICING_EPITOPES_DIR}/regtools_genomic_VCF/tvariant_{sample}_splice_effects.tsv")
            outputs.append(f"{SPLICING_EPITOPES_DIR}/pvacSplice/{sample}/")
    return outputs


rule all:
    input:
        generate_all_outputs()

def get_lanes(sample):
    # Function to retrieve all lanes for a given sample
    return RNA_sample_df[sample_df['sample_name'] == sample]['lane'].tolist()



rule convert_gtf:
    input:
        "resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf"
    output:
        "resources/updated_gtf/ref_annot_updated.gtf"
    conda:
        conda_env
    shell:
        """
        python scripts/final_scripts/convert_gtf.py -i {input} -o {output}
        """





rule download_ensembl_reference_files:
    output:
        fasta="resources/FASTA/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    shell:
        """
        mkdir -p resources/FASTA

        wget -O resources/FASTA/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
            ftp://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

        gunzip -c resources/FASTA/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > {output.fasta}
        """

rule extract_splice_sites_exons:
    input:
        gtf="resources/GTF/Homo_sapiens.GRCh38.103.gtf"
    output:
        splice_sites="resources/GTF/ensembl.splice_sites.txt",
        exons="resources/GTF/ensembl.exons.txt"
    conda:
        conda_env  # should include hisat2
    shell:
        """
        hisat2_extract_splice_sites.py {input.gtf} > {output.splice_sites}
        hisat2_extract_exons.py {input.gtf} > {output.exons}
        """



rule build_hisat2_index:
    input:
        fasta="resources/FASTA/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        splice_sites="resources/GTF/ensembl.splice_sites.txt",
        exons="resources/GTF/ensembl.exons.txt"
    output:
        index=expand("resources/HISAT2_index/ensembl_GRCh38.{i}.ht2", i=range(1, 9))
    params:
        prefix="resources/HISAT2_index/ensembl_GRCh38"
    threads: 20
    conda:
        conda_env  # should include hisat2
    shell:
        """
        mkdir -p resources/HISAT2_index
        hisat2-build -p {threads} \
            --ss {input.splice_sites} \
            --exon {input.exons} \
            {input.fasta} {params.prefix}
        """

rule hisat2_align_ensembl:
    input:
        r1 = lambda wc: f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{wc.sample}_CancerRNA_{wc.lane}_trimmed_R1.fastq.gz",
        r2 = lambda wc: f"{RESULTS_DIR}/0_Filtering_and_QC/trimmed_fastq/{wc.sample}_CancerRNA_{wc.lane}_trimmed_R2.fastq.gz",
        index = expand("resources/HISAT2_index/ensembl_GRCh38.{i}.ht2", i=range(1, 9)),
        splice_sites = "resources/GTF/ensembl.splice_sites.txt"
    output:
        bam = f"{SPLICING_EPITOPES_DIR}/hisat2/{{sample}}_CancerRNA_{{lane}}_sorted.bam"
    threads: 20
    conda:
        conda_env  # should include hisat2 + samtools
    shell:
        """
        mkdir -p {RNA_ANALYSIS_DIR}/hisat2
        hisat2 -x resources/HISAT2_index/ensembl_GRCh38 \
               --known-splicesite-infile {input.splice_sites} \
               --rna-strandness RF \
               -p {threads} \
               -1 {input.r1} -2 {input.r2} | \
        samtools sort -@ {threads} -o {output.bam}
        """

#rule regtools_splice_effects:
#    input:
#        vcf=f"{ALT_SPLICING_DIR}/VCF/{{sample}}_CancerRNA_altsplicing.vcf.gz",
#        bam=lambda wildcards: expand(
#            f"{RNA_ANALYSIS_DIR}/sorted_bam/{{sample}}_CancerRNA_{{lane}}_sorted.bam",
#            sample=wildcards.sample,
#            lane=get_lanes(wildcards.sample)
#        ),
#        ref="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa",
#        gtf="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf"
#    output:
#        tsv=f"{SPLICING_EPITOPES_DIR}/regtools/{{sample}}_splice_effects.tsv"
#    conda:
#        conda_env
#    shell:
#        """
#        mkdir -p {SOMATIC_EPITOPES_DIR}/regtools

#        regtools cis-splice-effects identify \
#            -o {output.tsv} \
#            -s XS \
#            {input.vcf} \
#            {input.bam} \
#            {input.ref} \
#            {input.gtf}
#        """


rule index_gtf_for_vep:
    input:
        gtf = "resources/updated_gtf/ref_annot_updated.gtf"
    output:
        gtf_gz = "resources/updated_gtf/ref_annot_updated.gtf.gz",
        gtf_tbi = "resources/updated_gtf/ref_annot_updated.gtf.gz.tbi"
    conda:
        conda_env_VEP  # or use your main env if htslib/tabix is already there
    shell:
        """
        sort -k1,1 -k4,4n {input.gtf} > sorted.gtf
        bgzip -c sorted.gtf > {output.gtf_gz}
        tabix -p gff {output.gtf_gz}
        rm sorted.gtf
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
          --transcript_version \
          --pick
         """
#          --offline \
#         --dir_cache {params.cache_dir} \

rule somatic_vep_BAM_matching_fasta_ensembl:
    input:
        vcf=f"{MUTATION_ANALYSIS_DIR}/VCF_filtered/{{sample}}_somatic_filtered.vcf.gz",
#        fasta="resources/genome_DNA/hg38.fa"
        fasta="resources/FASTA/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        #gtf="resources/updated_gtf/ref_annot_updated.gtf.gz"
    output:
        vep_output=f"{SPLICING_EPITOPES_DIR}/annotated_somatic_VCF_ensembl/{{sample}}_somatic_VEP.vcf"
    params:
        cache_dir="resources/VEP",  # Directory for the VEP cache inside the Docker container
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
          --offline \
          --dir_cache {params.cache_dir} \
          --no_stats \
          --transcript_version \
          --pick
         """





rule regtools_splice_effects_genomic_VCF_ensembl:
    input:
        vcf=f"{SPLICING_EPITOPES_DIR}/annotated_somatic_VCF_ensembl/{{sample}}_somatic_VEP.vcf",
        bam=lambda wildcards: expand(
            f"{SPLICING_EPITOPES_DIR}/hisat2/{{sample}}_CancerRNA_{{lane}}_sorted.bam",
            sample=wildcards.sample,
            lane=get_lanes(wildcards.sample)
        ),
#        ref="resources/genome_DNA/hg38.fa",
        ref="resources/FASTA/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        gtf="resources/GTF/Homo_sapiens.GRCh38.103.gtf"
    output:
        tsv=f"{SPLICING_EPITOPES_DIR}/regtools_genomic_VCF_ensembl/{{sample}}_splice_effects.tsv",
        vcf=f"{SPLICING_EPITOPES_DIR}/regtools_genomic_VCF_ensembl/{{sample}}_splice_effects.vcf"
    conda:
        conda_env
    shell:
        """
        mkdir -p {SOMATIC_EPITOPES_DIR}/regtools
        echo "Running: regtools cis-splice-effects identify -o {output.tsv} -s RF -w 100 -E -I {input.vcf} {input.bam} {input.ref} {input.gtf}"

        regtools cis-splice-effects identify \
            -o {output.tsv} \
            -v {output.vcf} \
            -s RF \
            -I \
            -E \
            -w 100 \
            {input.vcf} \
            {input.bam} \
            {input.ref} \
            {input.gtf}
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
 #       gtf="resources/updated_gtf/ref_annot_updated.gtf"
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
#           -I \
#            -E \


#rule add_transcript_version_column:
#    input:
#        tsv = f"{SPLICING_EPITOPES_DIR}/regtools_genomic_VCF/{{sample}}_splice_effects.tsv"
#    output:
#        patched_tsv = f"{SPLICING_EPITOPES_DIR}/regtools_genomic_VCF/tvariant_{{sample}}_splice_effects.tsv"
#    shell:
#        r"""
#        awk 'BEGIN{{FS=OFS="\t"}} \
#            NR==1{{$(NF+1)="transcript_version"}} \
#            NR>1{{split($(NF-1),a,","); split(a[1],b,"."); $(NF+1)=b[2]}}1' \
#            {input.tsv} > {output.patched_tsv}
#        """



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
            #main_vcf=f"{ALT_SPLICING_DIR}/VCF/{{sample}}_CancerRNA_altsplicing.vcf.gz",
            main_vcf=f"{SPLICING_EPITOPES_DIR}/annotated_somatic_VCF/{{sample}}_somatic_VEP.vcf",
            hla_types=lambda wildcards: expand(
            f"{HLA_ANALYSIS_DIR}/ArcasHLA/{{sample}}_CancerRNA_{{lane}}_HLAs_pvacFormat.txt",
            sample=wildcards.sample,
            lane=get_lanes(wildcards.sample)
        ),
            #regtools_tsv=f"{SPLICING_EPITOPES_DIR}/regtools/{{sample}}_splice_effects.tsv",
            regtools_tsv=f"{SPLICING_EPITOPES_DIR}/regtools_genomic_VCF_genecode/{{sample}}_splice_effects.tsv",
            ref="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa",
            gtf="resources/updated_gtf/ref_annot_updated.gtf"
#            gtf="resources/genome_RNA_fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf"
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
{params.algorithms} \
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
            star_fusion_tsv=f"{HLA_ANALYSIS_DIR}/StarFusionOut/{{sample}}_CancerRNA_{{lane}}/star-fusion.fusion_predictions.tsv",
            ag_fusion_directory=f"{FUSION_EPITOPES_DIR}/AGfusion/{{sample}}_CancerRNA_{{lane}}/"
        output:
            directory(f"{FUSION_EPITOPES_DIR}/pvacFuse/{{sample}}_CancerDNA_{{lane}}/")
        threads:
            pvacfuse_threads
        params:
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

