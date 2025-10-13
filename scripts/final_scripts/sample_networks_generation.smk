import os
import subprocess
import pandas as pd


# Define the results directory, defaulting to "results" if not specified
RESULTS_DIR = config.get("results_dir", "results")
EPITOPE_PRIORITISATION = f"{RESULTS_DIR}/sample_specific_networks"
RNA_ANALYSIS_DIR = f"{RESULTS_DIR}/1B_RNA_fusion_HLA"


# Get the temporary directory from config, defaulting to "/tmp" if not specified
TMPDIR = config.get("tmpdir", "/tmp")
os.environ["TMPDIR"] = TMPDIR
global_tmpdir = TMPDIR

# Ensure all necessary directories exist
os.makedirs(EPITOPE_PRIORITISATION, exist_ok=True)

# Load the CSV file from the command line config option
sample_df = pd.read_csv(config["csvfile"])

conda_env = "../envs/module_3_R.yaml"
conda_env_go = "../envs/conda_env_go.yaml"
conda_env_filtering = "../envs/module_3_filtering.yaml"
# Create a dictionary for sample inputs
sample_dict = {
    (row['sample_name'], row['datatype'], row['lane']): {
        'fastq_R1': row['fastq_R1'],
        'fastq_R2': row['fastq_R2'],
        'datatype': row['datatype']
    }
    for idx, row in sample_df.iterrows()
}




def generate_all_outputs():
    outputs = []
    f"{EPITOPE_PRIORITISATION}/normalised_counts/vst_mat.csv"
#    for (sample, datatype, lane) in sample_dict.keys():
#        if datatype == "CancerRNA":  # Only include outputs for CancerRNA
    return outputs
# Generate LIONESS output file paths for CancerRNA samples
def generate_lioness_outputs():
    outputs = []
    for (sample, datatype, lane), info in sample_dict.items():
        if datatype == "CancerRNA":
            outputs.append(f"{EPITOPE_PRIORITISATION}/Sample_Specific_Networks/{sample}.rds")
    return outputs

# Get list of sample names from CSV where datatype is CancerRNA
samples = list(sample_df[sample_df['datatype'] == "CancerRNA"]['sample_name'].unique())

rule all:
    input:
        f"{EPITOPE_PRIORITISATION}/normalised_counts/vst_mat.csv",
#        generate_lioness_outputs()
        f"{EPITOPE_PRIORITISATION}/gene_lists/membrane_ensembl.rds",
        f"{EPITOPE_PRIORITISATION}/gene_lists/HumanNet_ensembl.rds",
#        expand(os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_PPI_filtered", "filtered_networks_rds", "{sample}_filtered.rds"), sample=samples),
#        expand(os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_PPI_filtered", "filtered_networks_matrix", "{sample}.mtx"), sample=samples),
#        expand(os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_PPI_filtered", "filtered_networks_matrix", "{sample}_genes.txt"), sample=samples)
        # Centrality metrics outputs
        expand(os.path.join(EPITOPE_PRIORITISATION, "Network_Metrics_Betweenness", "{sample}_betweenness.tsv"), sample=samples),
        expand(os.path.join(EPITOPE_PRIORITISATION, "Network_Metrics_Degree", "{sample}_degree.tsv"), sample=samples),
        expand(os.path.join(EPITOPE_PRIORITISATION, "Network_Metrics_Strength", "{sample}_strength.tsv"), sample=samples),
        expand(os.path.join(EPITOPE_PRIORITISATION, "Network_Metrics_LargestComponentImpact", "{sample}_impact.tsv"), sample=samples),
        expand(os.path.join(EPITOPE_PRIORITISATION, "Network_Metrics_Full/Strength", "{sample}_strength.tsv"), sample=samples),
        expand(os.path.join(EPITOPE_PRIORITISATION, "Network_Metrics_Full/WCI", "{sample}_wci.tsv"), sample=samples),
        expand(os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_Distances", "{kind}_dist_matrix.rds"), kind=["filtered", "unfiltered"]),
        expand(os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_Distances", "{kind}_dist_heatmap.png"), kind=["filtered", "unfiltered"]),
        expand(os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_Distances", "{kind}_sim_heatmap.png"), kind=["filtered", "unfiltered"]),
        expand(os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_Distances", "{kind}_dist_matrix.csv"), kind=["filtered", "unfiltered"]),
        expand(os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_PPI_filtered", "qc_plots", "{sample}_qc_plots.pdf"), sample=samples)


# Generate raw count file paths for CancerRNA samples
def generate_raw_counts():
    counts = []
    for (sample, datatype, lane), info in sample_dict.items():
        if datatype == "CancerRNA":
            counts.append(f"{RNA_ANALYSIS_DIR}/RNA_Counts/{sample}_CancerRNA_{lane}_counts.txt")
    return counts



# Rule to run RNA vst normalisation
rule rna_analysis:
    input:
        counts=generate_raw_counts()
    output:
        f"{EPITOPE_PRIORITISATION}/normalised_counts/vst_mat.csv"
    conda:
        conda_env
    shell:
        """
        echo "Input count files: {input.counts}"
        Rscript scripts/final_scripts/R_scripts/RNA_analysis.R {output} {input.counts}
        """



# Rule to run LIONESS analysis
rule lioness_analysis:
    input:
        vst_mat=f"{EPITOPE_PRIORITISATION}/normalised_counts/vst_mat.csv"
    output:
        generate_lioness_outputs()
    conda:
        conda_env
    log:
        "logs/LIONESS_analysis.log"
    shell:
        """
        echo "Input vst_mat: {input.vst_mat}"
        Rscript scripts/final_scripts/R_scripts/LIONESS_analysis.R {input.vst_mat} {EPITOPE_PRIORITISATION}/Sample_Specific_Networks 2>&1 | tee {log}
        """


# Rule to run GO and HumanNet analysis
rule go_humannet_analysis:
    output:
        combined=f"{EPITOPE_PRIORITISATION}/gene_lists/combined_ensembl_ids.tsv",
        not_in_humannet=f"{EPITOPE_PRIORITISATION}/gene_lists/not_in_humannet.tsv",
        HumanNet=f"{EPITOPE_PRIORITISATION}/gene_lists/HumanNet-XN.tsv",
        humannet_rds =f"{EPITOPE_PRIORITISATION}/gene_lists/HumanNet_ensembl.rds",
        membrane_ensembl =f"{EPITOPE_PRIORITISATION}/gene_lists/membrane_ensembl.rds"
    conda:
        conda_env_go
    shell:
        """
        wget -O {output.HumanNet} https://www.inetbio.org/humannetv2/networks/HumanNet-XN.tsv
        echo "Input HumanNet file {output.HumanNet}"
        Rscript scripts/final_scripts/R_scripts/GO_HumanNet_analysis.R {output.HumanNet} {output.combined} {output.not_in_humannet} {output.humannet_rds} {output.membrane_ensembl}
        """


#rule filter_sample_networks:
#    input:
#        rds = os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks", "{sample}.rds"),
#        config = "scripts/final_scripts/config/filter_config.yaml",
#        humannet = f"{EPITOPE_PRIORITISATION}/gene_lists/HumanNet_ensembl.rds",
#        mem_genes = f"{EPITOPE_PRIORITISATION}/gene_lists/membrane_ensembl.rds"
#    output:
#        rds = os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_PPI_filtered", "filtered_networks_rds", "{sample}_filtered.rds"),
#        mm = os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_PPI_filtered", "filtered_networks_matrix", "{sample}.mtx"),
#        genes = os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_PPI_filtered", "filtered_networks_matrix", "{sample}_genes.txt")
#    conda:
#        conda_env_filtering
#    shell:
#        """
#        Rscript scripts/final_scripts/R_scripts/Filter_Sample_Networks.R \
#            {input.rds} \
#            {output.rds} \
#            {output.mm} \
#            {output.genes} \
#            {input.config} \
#            {input.humannet} \
#            {input.mem_genes}
#        """


rule network_metrics_igraph_rds:
    input:
        rds = os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_PPI_filtered", "filtered_networks_rds", "{sample}_filtered.rds")
    output:
        betweenness = os.path.join(EPITOPE_PRIORITISATION, "Network_Metrics_Betweenness", "{sample}_betweenness.tsv"),
        degree      = os.path.join(EPITOPE_PRIORITISATION, "Network_Metrics_Degree", "{sample}_degree.tsv"),
        strength    = os.path.join(EPITOPE_PRIORITISATION, "Network_Metrics_Strength", "{sample}_strength.tsv")
    conda:
        conda_env_filtering
    shell:
        """
        Rscript scripts/final_scripts/R_scripts/compute_network_metrics.R {input.rds} {output.betweenness} {output.degree} {output.strength}
        """ 






rule network_largest_component_impact:
    input:
        rds = os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_PPI_filtered", "filtered_networks_rds", "{sample}_filtered.rds")
    output:
        impact = os.path.join(EPITOPE_PRIORITISATION, "Network_Metrics_LargestComponentImpact", "{sample}_impact.tsv")
    conda:
        conda_env_filtering
    shell:
        """
        Rscript scripts/final_scripts/R_scripts/compute_largest_component_impact.R {input.rds} {output.impact}
        """


rule network_metrics_wci_and_strength:
    input:
        rds = os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks", "{sample}.rds")
    output:
        strength = os.path.join(EPITOPE_PRIORITISATION, "Network_Metrics_Full/Strength", "{sample}_strength.tsv"),
        wci      = os.path.join(EPITOPE_PRIORITISATION, "Network_Metrics_Full/WCI", "{sample}_wci.tsv"),
    conda:
        conda_env_filtering
    shell:
        """
        Rscript scripts/final_scripts/R_scripts/compute_WCI_and_strength.R {input.rds} {output.strength} {output.wci}
        """


rule cosine_distance_filtered:
    input:
        rds = expand(os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_PPI_filtered", "filtered_networks_rds", "{sample}_filtered.rds"), sample=samples)
    output:
        dist_mat = os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_Distances", "filtered_dist_matrix.rds"),
        heatmap = os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_Distances", "filtered_dist_heatmap.png"),
        simmap  = os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_Distances", "filtered_sim_heatmap.png"),
        csv     = os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_Distances", "filtered_dist_matrix.csv")
    conda:
        conda_env_filtering
    shell:
        """
        Rscript scripts/final_scripts/R_scripts/compute_cosine_distance.R \
            {EPITOPE_PRIORITISATION}/Sample_Specific_Networks_Distances/filtered \
            "Filtered PPI networks" \
            {input.rds}
        """

rule cosine_distance_unfiltered:
    input:
        rds = expand(os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks", "{sample}.rds"), sample=samples)
    output:
        dist_mat = os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_Distances", "unfiltered_dist_matrix.rds"),
        heatmap = os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_Distances", "unfiltered_dist_heatmap.png"),
        simmap  = os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_Distances", "unfiltered_sim_heatmap.png"),
        csv     = os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_Distances", "unfiltered_dist_matrix.csv")
    conda:
        conda_env_filtering
    shell:
        """
        Rscript scripts/final_scripts/R_scripts/compute_cosine_distance.R \
            {EPITOPE_PRIORITISATION}/Sample_Specific_Networks_Distances/unfiltered \
            "Unfiltered networks" \
            {input.rds}
        """



rule network_qc_plots:
    input:
        rds = os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_PPI_filtered", "filtered_networks_rds", "{sample}_filtered.rds")
    output:
        pdf = os.path.join(EPITOPE_PRIORITISATION, "Sample_Specific_Networks_PPI_filtered", "qc_plots", "{sample}_qc_plots.pdf")
    conda:
        conda_env_filtering
    shell:
        """
        Rscript scripts/final_scripts/R_scripts/qc_network_plots.R {input.rds} {output.pdf}
        """











