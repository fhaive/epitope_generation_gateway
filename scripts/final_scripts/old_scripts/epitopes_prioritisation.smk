import os
import subprocess
import pandas as pd

# ---------- paths ----------
RESULTS_DIR = config.get("results_dir", "results")
SOMATIC_EPITOPES_DIR   = f"{RESULTS_DIR}/2A_somatic_mutation_epitopes"
FUSION_EPITOPES_DIR    = f"{RESULTS_DIR}/2B_fusion_epitopes"
SPLICING_EPITOPES_DIR  = f"{RESULTS_DIR}/2C_splicing_epitopes"
SAMPLE_NETWORKS_DIR    = f"{RESULTS_DIR}/sample_specific_networks"
PRIORITISATION_DIR     = f"{RESULTS_DIR}/epitopes_prioritisation"

# tmp
TMPDIR = config.get("tmpdir", "/tmp")
os.environ["TMPDIR"] = TMPDIR
global_tmpdir = TMPDIR

# make dirs you need
os.makedirs(SOMATIC_EPITOPES_DIR, exist_ok=True)
os.makedirs(PRIORITISATION_DIR, exist_ok=True)

# ---------- envs ----------
conda_env = "../envs/module_prioritisation.yaml"
conda_env_prioritisation = "../envs/module_prioritisation2.yaml"

# ---------- unique sample names ----------
# Load the CSV file from the command line config option
sample_df = pd.read_csv(config["csvfile"])


# Create a dictionary for sample inputs
sample_dict = {
    (row['sample_name'], row['datatype'], row['lane']): {
        'fastq_R1': row['fastq_R1'],
        'fastq_R2': row['fastq_R2'],
        'datatype': row['datatype']
    }
    for idx, row in sample_df.iterrows()
}

sample_df = pd.read_csv(config["csvfile"])
unique_samples = sorted(sample_df["sample_name"].unique())

# ---------- depmap config ----------
DEPMAP_GENE_EFFECT_URL = config.get("depmap_gene_effect_url", "")
DEPMAP_FIGSHARE_ID = str(config.get("depmap_figshare_article_id", "27993248"))

CRISPR_GENE_EFFECT = f"{PRIORITISATION_DIR}/CRISPRGeneEffect.csv"
PAN_GENE_SCORES    = f"{PRIORITISATION_DIR}/depmap_pan_cancer_gene_score.csv"

# ---------- intogen config/paths ----------
INTOGEN_ZIP_URL   = "https://www.intogen.org/download?file=IntOGen-Drivers-20230531.zip"
INTOGEN_ZIP       = f"{PRIORITISATION_DIR}/IntOGen-Drivers-20230531.zip"
INTOGEN_COMP_TSV  = f"{PRIORITISATION_DIR}/Compendium_Cancer_Genes.tsv"
INTOGEN_BY_GENE   = f"{PRIORITISATION_DIR}/intogen_compendium_by_gene.csv"

def generate_all_outputs():
    outputs = []
    outputs.append(PAN_GENE_SCORES)
    outputs.append(INTOGEN_BY_GENE)  # ensure IntOGen rule runs
    return outputs

rule all:
    input:
        generate_all_outputs()

# ---------- rules ----------
rule depmap_pan_cancer_gene_scores:
    output:
        raw   = CRISPR_GENE_EFFECT,
        score = PAN_GENE_SCORES
    params:
        url       = DEPMAP_GENE_EFFECT_URL,
        fig_id    = DEPMAP_FIGSHARE_ID,
        outdir    = PRIORITISATION_DIR
    conda:
        conda_env
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"

        URL="{params.url}"
        if [ -z "$URL" ]; then
          echo "[depmap] Resolving CRISPRGeneEffect.csv from Figshare article {params.fig_id}"
          URL="$(curl -s "https://api.figshare.com/v2/articles/{params.fig_id}" \
                 | jq -r '.files[] | select(.name=="CRISPRGeneEffect.csv") | .download_url')"
          if [ -z "$URL" ] || [ "$URL" = "null" ]; then
            echo "ERROR: Could not find CRISPRGeneEffect.csv in Figshare article {params.fig_id}" >&2
            exit 1
          fi
        fi

        echo "[depmap] Downloading CRISPRGeneEffect.csv"
        curl -L -f -S "$URL" -o "{output.raw}"

        test -s "{output.raw}" || (echo "Downloaded file is empty" >&2; exit 1)
        head -c 1024 "{output.raw}" | grep -qi "<html" && {{ echo "Got HTML instead of CSV (likely wrong URL)"; exit 1; }}

        echo "[depmap] Computing pan-cancer median per gene (Chronos)"
        Rscript -e '
          suppressPackageStartupMessages({{ library(data.table) }})
          infile <- "{output.raw}"
          outfile <- "{output.score}"

          x <- fread(infile)
          id_cols <- intersect(names(x), c("ModelID","model_id","DepMap_ID","depmap_id"))
          gene_cols <- setdiff(names(x), id_cols)

          # Split "SYMBOL (ENTREZ)" into symbol & EntrezID
          gene_info <- tstrsplit(gene_cols, " \\(|\\)", keep = c(1, 2))
          gene_symbols <- gene_info[[1]]
          entrez_ids   <- gene_info[[2]]

          M <- as.matrix(x[, ..gene_cols])
          storage.mode(M) <- "double"

          meds <- apply(M, 2, median, na.rm = TRUE)

          out <- data.table(
            Gene     = gene_symbols,
            EntrezID = entrez_ids,
            PanCancer_Median_Chronos = as.numeric(meds)
          )
          fwrite(out, outfile)
        '

        echo "[depmap] Wrote: {output.score}"
        """

rule intogen_cancer_driver_genes_table:
    """
    Download IntOGen driver gene compendium, keep all cancer driver genes
    (oncogenes 'Act' and tumor suppressors 'LoF'), and collapse to one row per gene.
    """
    output:
        table = f"{PRIORITISATION_DIR}/intogen_compendium_by_gene.csv"
    params:
        url    = "https://www.intogen.org/download?file=IntOGen-Drivers-20230531.zip",
        outdir = PRIORITISATION_DIR
    conda:
        conda_env
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.outdir}"

        echo "[intogen] Downloading {params.url}"
        curl -L -f -S "{params.url}" -o "{params.outdir}/IntOGen-Drivers.zip"

        echo "[intogen] Locating Compendium_Cancer_Genes.tsv inside ZIP"
        tsv_path=$(unzip -Z1 "{params.outdir}/IntOGen-Drivers.zip" | grep "Compendium_Cancer_Genes.tsv" | head -n1)

        if [ -z "$tsv_path" ]; then
            echo "ERROR: Could not find Compendium_Cancer_Genes.tsv in ZIP" >&2
            exit 1
        fi

        echo "[intogen] Extracting TSV"
        unzip -p "{params.outdir}/IntOGen-Drivers.zip" "$tsv_path" > "{params.outdir}/Compendium_Cancer_Genes.tsv"

        echo "[intogen] Collapsing to one row per gene (CANCER DRIVER GENES: Act + LoF)"
        Rscript -e '
          suppressPackageStartupMessages({{ library(data.table) }})
          infile  <- "{params.outdir}/Compendium_Cancer_Genes.tsv"
          outfile <- "{output.table}"

          x <- fread(infile, sep = "\t", header = TRUE)

          # Normalize column names
          setnames(x, tolower(names(x)))

          # Keep only rows where ROLE contains Act or LoF (case-insensitive)
          x <- x[ grepl("act", tolower(role), fixed = TRUE) | grepl("lof", tolower(role), fixed = TRUE) ]

          # Aggregate per SYMBOL
          agg <- x[, .(
              ROLE       = paste(sort(unique(role[!is.na(role) & role != ""])), collapse = ";"),
              TRANSCRIPT = paste(sort(unique(transcript[!is.na(transcript) & transcript != ""])), collapse = ";"),
              COHORT     = paste(sort(unique(cohort[!is.na(cohort) & cohort != ""])), collapse = ";"),
              CANCER_TYPE= paste(sort(unique(cancer_type[!is.na(cancer_type) & cancer_type != ""])), collapse = ";")
          ), by = .(SYMBOL = symbol)]

          fwrite(agg, outfile)
        '

        echo "[intogen] Wrote: {output.table}"
        """



# ---------- config-driven columns for epitopes ----------
EPITOPE_COLUMNS_YAML = config.get("epitope_columns_config", "config/epitope_columns.yaml")


# Utility: expected file locations for each sample
def epitope_candidates(sample):
    return {
        "somatic_mhci":  f"{SOMATIC_EPITOPES_DIR}/pvacSeq/{sample}_CancerDNA/MHC_Class_I/{sample}_CancerDNA.filtered.tsv",
        "somatic_mhcii": f"{SOMATIC_EPITOPES_DIR}/pvacSeq/{sample}_CancerDNA/MHC_Class_II/{sample}_CancerDNA.filtered.tsv",
        "fusion_mhci":   f"{FUSION_EPITOPES_DIR}/pvacFuse/{sample}_CancerRNA_{lane}/MHC_Class_I/{sample}_CancerRNA_{lane}.filtered.tsv",
        "fusion_mhcii":  f"{FUSION_EPITOPES_DIR}/pvacFuse/{sample}_CancerRNA_{lane}/MHC_Class_II/{sample}_CancerRNA_{lane}.filtered.tsv",
        "splice_mhci":   f"{SPLICING_EPITOPES_DIR}/pvacSplice/{sample}_CancerDNA/MHC_Class_I/{sample}_CancerDNA.filtered.tsv",
        "splice_mhcii":  f"{SPLICING_EPITOPES_DIR}/pvacSplice/{sample}_CancerDNA/MHC_Class_II//{sample}_CancerDNA.filtered.tsv",
    }

# Collect only files that actually exist so the rule can run even if some are absent
def existing_epitope_inputs(wc):
    cand = epitope_candidates(wc.sample)
    return [p for p in cand.values() if os.path.exists(p)]

# Add per-sample merged outputs to the "all" target
def generate_all_outputs():
    outputs = []
    outputs.append(PAN_GENE_SCORES)
    outputs.append(INTOGEN_BY_GENE)
    # add per-sample merged epitopes
    for s in unique_samples:
        outputs.append(f"{PRIORITISATION_DIR}/{s}_epitopes_merged_ic50.csv")
    return outputs

# ---------- rule: merge epitopes using IC50 (per sample) ----------
rule merge_epitopes_ic50:
    input:
        existing_epitope_inputs
    output:
        merged = lambda wc: f"{PRIORITISATION_DIR}/{wc.sample}_epitopes_merged_ic50.csv"
    params:
        somatic_mhci  = lambda wc: epitope_candidates(wc.sample)["somatic_mhci"],
        somatic_mhcii = lambda wc: epitope_candidates(wc.sample)["somatic_mhcii"],
        fusion_mhci   = lambda wc: epitope_candidates(wc.sample)["fusion_mhci"],
        fusion_mhcii  = lambda wc: epitope_candidates(wc.sample)["fusion_mhcii"],
        splice_mhci   = lambda wc: epitope_candidates(wc.sample)["splice_mhci"],
        splice_mhcii  = lambda wc: epitope_candidates(wc.sample)["splice_mhcii"],
        columns_yaml  = EPITOPE_COLUMNS_YAML,
        out_csv       = lambda wc: f"{PRIORITISATION_DIR}/{wc.sample}_epitopes_merged_ic50.csv",
        scripts_dir   = "scripts"
    conda:
        conda_env_prioritisation
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{PRIORITISATION_DIR}"

        # Ensure script is available
        if [ ! -f "{params.scripts_dir}/merge_epitopes_ic50.R" ]; then
          echo "ERROR: {params.scripts_dir}/merge_epitopes_ic50.R not found." >&2
          exit 1
        fi

        Rscript "{params.scripts_dir}/merge_epitopes_ic50.R" \
          --somatic_mhci  "{params.somatic_mhci}" \
          --somatic_mhcii "{params.somatic_mhcii}" \
          --fusion_mhci   "{params.fusion_mhci}" \
          --fusion_mhcii  "{params.fusion_mhcii}" \
          --splice_mhci   "{params.splice_mhci}" \
          --splice_mhcii  "{params.splice_mhcii}" \
          --columns_yaml  "{params.columns_yaml}" \
          --out_csv       "{params.out_csv}"

        test -s "{output.merged}"
        """











