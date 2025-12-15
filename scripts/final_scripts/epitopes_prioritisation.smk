import os
import pandas as pd
from glob import glob

# ---------- paths ----------
RESULTS_DIR = config.get("results_dir", "results")
SOMATIC_EPITOPES_DIR   = f"{RESULTS_DIR}/2A_somatic_mutation_epitopes"
FUSION_EPITOPES_DIR    = f"{RESULTS_DIR}/2B_fusion_epitopes"
SPLICING_EPITOPES_DIR  = f"{RESULTS_DIR}/2C_splicing_epitopes"
SAMPLE_NETWORKS_DIR    = f"{RESULTS_DIR}/sample_specific_networks"
PRIORITISATION_DIR     = f"{RESULTS_DIR}/epitopes_prioritisation"
COMBINED_EPITOPES_DIR = f"{PRIORITISATION_DIR}/combined_epitopes"
ANNOTATED_EPITOPES_DIR = f"{PRIORITISATION_DIR}/annotated_epitopes"
ANNOTATED_WITH_NETWORK_DIR = f"{PRIORITISATION_DIR}/annotated_epitopes_with_network"
FINAL_EPITOPES_DIR = f"{PRIORITISATION_DIR}/final_epitopes"
ARCASHLA_DIR = f"{RESULTS_DIR}/1B_RNA_fusion_HLA/ArcasHLA"

# --- Borda config/paths ---
BORDA_DIR     = f"{PRIORITISATION_DIR}/Borda_Epiotpes_Prioritisation"
BORDA_CONFIG  = config.get("borda_config", "scripts/final_scripts/config/borda_config.yaml")
COMPUTE_BORDA_R = "scripts/final_scripts/R_scripts/compute_borda_rank.R"

os.makedirs(BORDA_DIR, exist_ok=True)

# ---- AnnotationHub cache ----
ANNOTATIONHUB_CACHE = config.get(
    "annotationhub_cache",
    f"{PRIORITISATION_DIR}/annotationhub_cache"
)
os.environ["ANNOTATIONHUB_CACHE"] = ANNOTATIONHUB_CACHE
os.makedirs(ANNOTATIONHUB_CACHE, exist_ok=True)

envvars:
    "ANNOTATIONHUB_CACHE"

# R script path
MERGE_EPITOPES_R = "scripts/final_scripts/R_scripts/merge_epitopes_ic50.R"
ANNOTATE_EPITOPES_R = "scripts/final_scripts/R_scripts/annotate_epitopes_depmap_intogen.R"
APPEND_NETWORK_R = "scripts/final_scripts/R_scripts/append_network_metrics.R"
GO_HALLMARK_R      = "scripts/final_scripts/R_scripts/annotate_go_hallmark.R"

# tmp
TMPDIR = config.get("tmpdir", "/tmp")
os.environ["TMPDIR"] = TMPDIR
global_tmpdir = TMPDIR

# make dirs needed
os.makedirs(SOMATIC_EPITOPES_DIR, exist_ok=True)
os.makedirs(PRIORITISATION_DIR, exist_ok=True)
os.makedirs(COMBINED_EPITOPES_DIR, exist_ok=True)
os.makedirs(ANNOTATED_EPITOPES_DIR, exist_ok=True)
os.makedirs(ANNOTATED_WITH_NETWORK_DIR, exist_ok=True)
os.makedirs(FINAL_EPITOPES_DIR, exist_ok=True)

# ---------- envs ----------
conda_env = "../envs/module_prioritisation.yaml"
conda_env_prioritisation = "../envs/module_prioritisation2.yaml"
conda_env_go_hallmark = "../envs/go_hallmark_env.yaml"

# ---------- samples ----------
sample_df = pd.read_csv(config["csvfile"])
unique_samples = sorted(sample_df["sample_name"].unique())

# ---------- depmap config ----------
DEPMAP_GENE_EFFECT_URL = config.get("depmap_gene_effect_url", "")
DEPMAP_FIGSHARE_ID     = str(config.get("depmap_figshare_article_id", "27993248"))

CRISPR_GENE_EFFECT = f"{PRIORITISATION_DIR}/CRISPRGeneEffect.csv"
PAN_GENE_SCORES    = f"{PRIORITISATION_DIR}/depmap_pan_cancer_gene_score.csv"

# path to Hallmark GMT 
HALLMARK_GMT = config.get("hallmark_gmt", "")

# ---------- intogen config/paths ----------
INTOGEN_ZIP_URL   = "https://www.intogen.org/download?file=IntOGen-Drivers-20230531.zip"
INTOGEN_ZIP       = f"{PRIORITISATION_DIR}/IntOGen-Drivers-20230531.zip"
INTOGEN_COMP_TSV  = f"{PRIORITISATION_DIR}/Compendium_Cancer_Genes.tsv"
INTOGEN_BY_GENE   = f"{PRIORITISATION_DIR}/intogen_compendium_by_gene.csv"

# ---------- CACHE sentinel paths ----------
RESET_AH_SENTINEL     = f"{ANNOTATIONHUB_CACHE}/_reset_ok"
ENSDB_V103_SENTINEL   = f"{ANNOTATIONHUB_CACHE}/EnsDb.Hsapiens.v103.ready"

# ---------- config-driven columns for epitopes ----------
EPITOPE_COLUMNS_YAML = config.get("epitope_columns_config", "scripts/final_scripts/config/epitope_columns.yaml")

# ---------- helper: candidate path patterns  ----------
from glob import glob


def epitope_patterns(sample):
    return {
        # pVACseq (somatic) — dir is exactly {sample}, file ends with _CancerDNA.filtered.tsv
        "somatic_mhci":  f"{SOMATIC_EPITOPES_DIR}/pvacSeq/{sample}/MHC_Class_I/{sample}_CancerDNA.filtered.tsv",
        "somatic_mhcii": f"{SOMATIC_EPITOPES_DIR}/pvacSeq/{sample}/MHC_Class_II/{sample}_CancerDNA.filtered.tsv",

        # pVACfuse (fusion) — lane wildcard keeps working with full {sample}
        "fusion_mhci":   f"{FUSION_EPITOPES_DIR}/pvacFuse/{sample}_CancerRNA_*/MHC_Class_I/{sample}_CancerRNA_*.filtered.tsv",
        "fusion_mhcii":  f"{FUSION_EPITOPES_DIR}/pvacFuse/{sample}_CancerRNA_*/MHC_Class_II/{sample}_CancerRNA_*.filtered.tsv",

        # pVACsplice (splicing) — dir is exactly {sample}, file ends with _CancerDNA.filtered.tsv
        "splice_mhci":   f"{SPLICING_EPITOPES_DIR}/pvacSplice/{sample}/MHC_Class_I/{sample}_CancerDNA.filtered.tsv",
        "splice_mhcii":  f"{SPLICING_EPITOPES_DIR}/pvacSplice/{sample}/MHC_Class_II/{sample}_CancerDNA.filtered.tsv",
    }

def _pick_first(pattern):
    hits = sorted(glob(pattern))
    return hits[0] if hits else os.devnull

def existing_epitope_inputs_list(wc):
    """Return a LIST (fixed order) to satisfy older Snakemake:
       [som_mhci, som_mhcii, fus_mhci, fus_mhcii, spl_mhci, spl_mhcii]"""
    pats = epitope_patterns(wc.sample)
    return [
        _pick_first(pats["somatic_mhci"]),
        _pick_first(pats["somatic_mhcii"]),
        _pick_first(pats["fusion_mhci"]),
        _pick_first(pats["fusion_mhcii"]),
        _pick_first(pats["splice_mhci"]),
        _pick_first(pats["splice_mhcii"]),
    ]

# ---------- all outputs ----------
def generate_all_outputs():
    outputs = []
    #outputs.append(PAN_GENE_SCORES)
    #outputs.append(INTOGEN_BY_GENE)
    for s in unique_samples:
        outputs.append(f"{COMBINED_EPITOPES_DIR}/{s}_epitopes_merged.csv")
        outputs.append(f"{ANNOTATED_EPITOPES_DIR}/{s}_epitopes_annotated.csv")
        outputs.append(f"{ANNOTATED_WITH_NETWORK_DIR}/{s}_epitopes_annotated_network.csv")
        outputs.append(f"{BORDA_DIR}/{s}_epitopes_borda.csv")
        outputs.append(f"{FINAL_EPITOPES_DIR}/{s}_epitopes_final.csv")
        outputs.append(f"{FINAL_EPITOPES_DIR}/HLA/{s}_epitopes_final.with_hla.csv")
        outputs.append(f"{FINAL_EPITOPES_DIR}/HTML/{s}_epitopes.html")
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
        head -c 1024 "{output.raw}" | grep -qi "<html" && { echo "Got HTML instead of CSV (likely wrong URL)"; exit 1; }

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
    output:
        table = INTOGEN_BY_GENE
    params:
        url    = INTOGEN_ZIP_URL,
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

# ---------- rule: merge epitopes using IC50 (per sample) ----------
# Point to R script
MERGE_EPITOPES_R = "scripts/final_scripts/R_scripts/merge_epitopes_ic50.R"
rule merge_epitopes_ic50:
    input:
        existing_epitope_inputs_list
    output:
        merged = f"{COMBINED_EPITOPES_DIR}/{{sample}}_epitopes_merged.csv"
    params:
        columns_yaml = EPITOPE_COLUMNS_YAML,
        script       = MERGE_EPITOPES_R
    conda:
        conda_env_prioritisation
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.merged}")"

        if [ ! -f "{params.script}" ]; then
          echo "ERROR: {params.script} not found." >&2
          exit 1
        fi

        # input order:
        # 0 somatic_mhci, 1 somatic_mhcii, 2 fusion_mhci,
        # 3 fusion_mhcii, 4 splice_mhci, 5 splice_mhcii
        Rscript "{params.script}" \
          --somatic_mhci  "{input[0]}" \
          --somatic_mhcii "{input[1]}" \
          --fusion_mhci   "{input[2]}" \
          --fusion_mhcii  "{input[3]}" \
          --splice_mhci   "{input[4]}" \
          --splice_mhcii  "{input[5]}" \
          --columns_yaml  "{params.columns_yaml}" \
          --out_csv       "{output.merged}"

        test -s "{output.merged}"
        """


rule annotate_epitopes_depmap_intogen:
    input:
        merged  = f"{COMBINED_EPITOPES_DIR}/{{sample}}_epitopes_merged.csv"
    output:
        annotated = f"{ANNOTATED_EPITOPES_DIR}/{{sample}}_epitopes_annotated.csv"
    params:
        script  = ANNOTATE_EPITOPES_R,
        depmap  = PAN_GENE_SCORES,      
        intogen = INTOGEN_BY_GENE       
    conda:
        conda_env_prioritisation
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.annotated}")"


        Rscript "{params.script}" \
          --merged_in "{input.merged}" \
          --depmap    "{params.depmap}" \
          --intogen   "{params.intogen}" \
          --out_csv   "{output.annotated}"

        test -s "{output.annotated}"
    """
rule reset_annotationhub_cache:
    output:
        RESET_AH_SENTINEL
    run:
        import os, time
        cache = ANNOTATIONHUB_CACHE
        os.makedirs(cache, exist_ok=True)
        # remove only the two hub metadata files if present
        for fname in ("annotationhub.sqlite3", "annotationhub.index.rds"):
            fpath = os.path.join(cache, fname)
            if os.path.exists(fpath):
                try:
                    os.remove(fpath)
                except Exception:
                    pass
        with open(output[0], "w") as f:
            f.write(time.strftime("%Y-%m-%d %H:%M:%S") + " reset\n")

PRECACHE_R = "scripts/final_scripts/R_scripts/precache_ensdb_v103.R"

rule precache_ensdb_v103:
    input:
        RESET_AH_SENTINEL
    output:
        ENSDB_V103_SENTINEL
    params:
        script = PRECACHE_R
    conda:
        conda_env_prioritisation
    shell:
        r"""
        set -euo pipefail
        Rscript "{params.script}" --sentinel "{output}"
        """

rule annotate_epitopes_with_network:
    """
    Append per-gene network centrality metrics (betweenness, degree, impact,
    strength [filtered network], and WCI [full network]) to the already
    annotated epitope table. Also ensures EnsDb v103 is precached.
    """
    input:
        annotated = f"{ANNOTATED_EPITOPES_DIR}/{{sample}}_epitopes_annotated.csv",
        sentinel  = ENSDB_V103_SENTINEL
    output:
        out = f"{ANNOTATED_WITH_NETWORK_DIR}/{{sample}}_epitopes_annotated_network.csv"
    params:
        script     = APPEND_NETWORK_R,
        betweenness = lambda wc: f"{SAMPLE_NETWORKS_DIR}/Network_Metrics_Betweenness/{wc.sample}_betweenness.tsv",
        degree      = lambda wc: f"{SAMPLE_NETWORKS_DIR}/Network_Metrics_Degree/{wc.sample}_degree.tsv",
        impact      = lambda wc: f"{SAMPLE_NETWORKS_DIR}/Network_Metrics_LargestComponentImpact/{wc.sample}_impact.tsv",
        strength    = lambda wc: f"{SAMPLE_NETWORKS_DIR}/Network_Metrics_Strength/{wc.sample}_strength.tsv",
        wci         = lambda wc: f"{SAMPLE_NETWORKS_DIR}/Network_Metrics_Full/WCI/{wc.sample}_wci.tsv"
    conda:
        conda_env_prioritisation
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.out}")"

        Rscript "{params.script}" \
          --annot_in   "{input.annotated}" \
          --betweenness "{params.betweenness}" \
          --degree      "{params.degree}" \
          --impact      "{params.impact}" \
          --strength    "{params.strength}" \
          --wci         "{params.wci}" \
          --out_csv     "{output.out}"

        test -s "{output.out}"
        """


rule compute_borda_rank:
    """
    Compute per-epitope Borda score from configured features.
    Input: annotated + network metrics CSV.
    Output: same table + Borda_Score (and Borda_Rank).
    """
    input:
        annotated_net = f"{ANNOTATED_WITH_NETWORK_DIR}/{{sample}}_epitopes_annotated_network.csv"
    output:
        out = f"{BORDA_DIR}/{{sample}}_epitopes_borda.csv"
    params:
        script = COMPUTE_BORDA_R,
        cfg    = BORDA_CONFIG
    conda:
        conda_env_prioritisation
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.out}")"

        Rscript "{params.script}" \
          --input_csv "{input.annotated_net}" \
          --config    "{params.cfg}" \
          --out_csv   "{output.out}"

        test -s "{output.out}"
        """


rule annotate_go_and_hallmark:
    input:
        borda = f"{BORDA_DIR}/{{sample}}_epitopes_borda.csv"
    output:
        out   = f"{FINAL_EPITOPES_DIR}/{{sample}}_epitopes_final.csv"
    params:
        script   = GO_HALLMARK_R,
        hallmark = config.get("hallmark", "")   # "", "msigdbr", or path to .gmt
    conda:
        conda_env_go_hallmark
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.out}")"

        hallmark_arg=""
        if [ -n "{params.hallmark}" ]; then
          hallmark_arg="--hallmark {params.hallmark}"
        fi

        Rscript "{params.script}" \
          --input_csv "{input.borda}" \
          $hallmark_arg \
          --out_csv   "{output.out}"

        test -s "{output.out}"
        """


conda_env_epitope_html = "../envs/epitope_html_tools.yaml"

rule add_hla_counts:
    input:
        final_csv = f"{FINAL_EPITOPES_DIR}/{{sample}}_epitopes_final.csv",
        genes_json = lambda wc: f"{ARCASHLA_DIR}/{wc.sample}_CancerRNA_LRNA_sorted.genes.json",
        geno_json  = lambda wc: f"{ARCASHLA_DIR}/{wc.sample}_CancerRNA_LRNA_sorted.genotype.json"
    output:
        augmented = f"{FINAL_EPITOPES_DIR}/HLA/{{sample}}_epitopes_final.with_hla.csv"
    conda:
        conda_env_epitope_html
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{FINAL_EPITOPES_DIR}/HLA"
        python "scripts/final_scripts/python/augment_epitopes_with_hla.py" \
          --sample "{wildcards.sample}" \
          --final-in  "{input.final_csv}" \
          --arcas-dir "{ARCASHLA_DIR}" \
          --final-out "{output.augmented}"
        """

rule build_interactive_epitope_html:
    input:
        augmented = f"{FINAL_EPITOPES_DIR}/HLA/{{sample}}_epitopes_final.with_hla.csv",
        script    = f"scripts/final_scripts/python/build_epitope_table.py"
    output:
        html = f"{FINAL_EPITOPES_DIR}/HTML/{{sample}}_epitopes.html"
    conda:
        conda_env_epitope_html
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{FINAL_EPITOPES_DIR}/HTML"
        python "{input.script}" \
          -i "{input.augmented}" \
          -o "{output.html}" \
          --title "Epitope Prioritisation — {wildcards.sample}"
        """
