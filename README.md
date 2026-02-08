# bulk-rnaseq-deg-toolkit

Reusable downstream analysis scripts for bulk RNA-seq differential expression (DEG) result tables. Designed for (1) a clean GitHub portfolio artifact, (2) reuse in a future lab, and (3) immediate use on local DEG CSVs.

## What this repo does

Given DEG result tables (one file per comparison), the scripts generate:

- Volcano plots plus per-plot annotated TSVs (Up, Down, NS calls)
- GO Biological Process over-representation analysis (ORA) using clusterProfiler  
  - ALL terms and immune-only subsets  
  - Within-sex timepoint overlap visualizations (E12.5 vs E17.5)
- Driver inference from DEG signatures  
  - TF activity via DoRothEA (mouse, confidence A–C) plus VIPER  
  - Pathway activity via PROGENy (mouse)  
  - Heatmaps plus ranked tables plus top drivers per list
- GO term–gene network graphs per sex and timepoint, plus overlap networks  
  - Robust to duplicate vertices and duplicate edges

## Repo layout

Recommended structure:

- README.md
- scripts/
  - 01_volcano_brain.R
  - 02_go_bp_ora_brain.R
  - 03_networks_brain_timepoints.R
  - 06_driver_inference_dorothea_progeny.R
- example_inputs/ (optional demo CSVs)

## Input format

Each DEG CSV should contain at minimum:

- symbol (gene symbol)
- logfc (log2 fold change)
- adjpv (BH-adjusted p value / FDR)

Optional columns are allowed and ignored unless a script uses them.

Example file set used in this project:

- E12.5_Benzene_Brain_Female.csv
- E12.5_Benzene_Brain_Male.csv
- E17.5_Benzene_Brain_Female.csv
- E17.5_Benzene_Brain_Male.csv

## How to run

Each script supports a project directory passed as an argument, or a base directory defined by environment variable.

Option A: pass the project directory as an argument (run from repo root)

Rscript scripts/01_volcano_brain.R "C:/path/to/your/project_dir"  
Rscript scripts/02_go_bp_ora_brain.R "C:/path/to/your/project_dir"  
Rscript scripts/06_driver_inference_dorothea_progeny.R "C:/path/to/your/project_dir"  
Rscript scripts/03_networks_brain_timepoints.R "C:/path/to/your/project_dir"

Option B: set RNASEQ_BASE_DIR once, then run scripts without arguments

export RNASEQ_BASE_DIR="C:/path/to/your/project_dir"  
Rscript scripts/01_volcano_brain.R  
Rscript scripts/02_go_bp_ora_brain.R  
Rscript scripts/06_driver_inference_dorothea_progeny.R  
Rscript scripts/03_networks_brain_timepoints.R

Windows PowerShell variant

setx RNASEQ_BASE_DIR "C:\path\to\your\project_dir"  
Rscript scripts/01_volcano_brain.R  
Rscript scripts/02_go_bp_ora_brain.R  
Rscript scripts/06_driver_inference_dorothea_progeny.R  
Rscript scripts/03_networks_brain_timepoints.R

## Dependency install behavior (GitHub-friendly)

By default, scripts do not auto-install packages. If packages are missing, the script stops and lists what to install.

To allow auto-install on a new machine, set:

RNASEQ_TOOLKIT_INSTALL_DEPS=true

Then re-run the script.

## Outputs

All outputs are written inside your project directory.

Volcano outputs

- GO_BP_Results_Brain/volcano_brain/
  - volcano_<label>.png
  - volcano_<label>.tsv
  - volcano_counts_brain.tsv

GO ORA outputs

- GO_BP_Results_Brain/
  - ORA_BP_<label>_ALL.csv and dotplot
  - ORA_BP_<label>_IMMUNE.csv and dotplot
  - overlap plots
  - diagnostics.txt

Driver inference outputs

- Driver_Inference/
  - DoRothEA_TF_activity_NES.csv
  - DoRothEA_TF_activity_top40_heatmap.png
  - PROGENy_pathway_activity_NES.csv
  - PROGENy_pathway_activity_top14_heatmap.png
  - Top10_TFs_per_list.csv
  - Top10_Pathways_per_list.csv
  - sessionInfo.txt

Network outputs

- Networks_Brain/
  - <Sex>_E12.5_GO_network.png and hubs table
  - <Sex>_E17.5_GO_network.png and hubs table
  - <Sex>_Overlap_GO_network.png and hubs table
  - sessionInfo.txt

## Notes on interpretation

VIPER and PROGENy are generally most stable when run on a ranked signature across many genes. If you only have significant-only DEG lists, results can be more sensitive to thresholds. If you later switch to full DE tables (all tested genes with logFC), consider relaxing DEG filtering in the driver inference script.

GO ORA depends strongly on the gene universe and on how the DEG list was defined. Keep DEG generation consistent across comparisons for interpretable differences.

## Customization points (safe edits)

Each script exposes parameters near the top, for example:

- Volcano: adjpv_cut, use_lfc_cut, lfc_cut, label_top_n
- GO ORA: fdr_cut, optional LFC gating, immune-term filtering rules
- Driver inference: deg_fdr_cut, use_lfc_cut
- Networks: TOP_TERMS, JACCARD_THRESH, GO_TERM_FDR_CUT

## Suggested citation

If you use this repo in a manuscript, cite the upstream methods or packages:

- clusterProfiler (GO enrichment / ORA)
- DoRothEA plus VIPER (TF activity inference)
- PROGENy (pathway activity inference)
