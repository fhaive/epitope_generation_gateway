# metrics_only.smk
import os
import pandas as pd

# ---- Config & paths (same pattern you used) ---------------------------------
RESULTS_DIR = config.get("results_dir", "results")
EPITOPE_PRIORITISATION = f"{RESULTS_DIR}/sample_specific_networks"

# optional TMPDIR handling (harmless if unused by this rule)
TMPDIR = config.get("tmpdir", "/tmp")
os.environ["TMPDIR"] = TMPDIR

# env for this rule
conda_env_filtering = "../envs/module_3_filtering.yaml"

# ---- Samples from your CSV (same columns expected) --------------------------
# expects --config csvfile=path/to/your.csv (just like your original setup)
sample_df = pd.read_csv(config["csvfile"])
# use the same CancerRNA filter you used elsewhere
samples = list(
    sample_df[sample_df["datatype"] == "CancerRNA"]["sample_name"].unique()
)

# ---- Targets: the three metrics per sample ----------------------------------
def _betw_out(s):    return os.path.join(EPITOPE_PRIORITISATION, "Network_Metrics_Betweenness", f"{s}_betweenness.tsv")
def _degree_out(s):  return os.path.join(EPITOPE_PRIORITISATION, "Network_Metrics_Degree",      f"{s}_degree.tsv")
def _strength_out(s):return os.path.join(EPITOPE_PRIORITISATION, "Network_Metrics_Strength",    f"{s}_strength.tsv")
def _impact_out(s):  return os.path.join(EPITOPE_PRIORITISATION, "Network_Metrics_LargestComponentImpact", f"{s}_impact.tsv")

rule all:
    input:
        # Only these outputs are targeted; nothing else will run.
        expand(_betw_out("{sample}"),   sample=samples),
        expand(_degree_out("{sample}"), sample=samples),
        expand(_strength_out("{sample}"), sample=samples),
        expand(_impact_out("{sample}"), sample=samples)

# ---- The only rule: compute metrics from *existing* filtered RDS ------------
rule network_metrics_igraph_rds:
    input:
        rds = os.path.join(
            EPITOPE_PRIORITISATION,
            "Sample_Specific_Networks_PPI_filtered",
            "filtered_networks_rds",
            "{sample}_filtered.rds"
        )
    output:
        betweenness = _betw_out("{sample}"),
        degree      = _degree_out("{sample}"),
        strength    = _strength_out("{sample}")
    conda:
        conda_env_filtering
    shell:
        r"""
        # ensure output dirs exist (no-ops if already there)
        mkdir -p "{EPITOPE_PRIORITISATION}/Network_Metrics_Betweenness" \
                 "{EPITOPE_PRIORITISATION}/Network_Metrics_Degree" \
                 "{EPITOPE_PRIORITISATION}/Network_Metrics_Strength"

        Rscript scripts/final_scripts/R_scripts/compute_network_metrics.R \
            "{input.rds}" "{output.betweenness}" "{output.degree}" "{output.strength}"
        """
rule network_largest_component_impact:
    input:
        rds = os.path.join(
            EPITOPE_PRIORITISATION,
            "Sample_Specific_Networks_PPI_filtered",
            "filtered_networks_rds",
            "{sample}_filtered.rds"
        )
    output:
        impact = os.path.join(
            EPITOPE_PRIORITISATION,
            "Network_Metrics_LargestComponentImpact",
            "{sample}_impact.tsv"
        )
    conda:
        conda_env_filtering
    shell:
        """
        Rscript scripts/final_scripts/R_scripts/compute_largest_component_impact.R {input.rds} {output.impact}
        """
