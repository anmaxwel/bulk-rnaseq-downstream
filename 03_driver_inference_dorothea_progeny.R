# ============================================================
# 06_driver_inference_dorothea_progeny.R
# Upstream driver inference from DEG lists:
# - TF activity: DoRothEA (mouse, conf A-C) + VIPER
# - Pathway activity: PROGENy (Mouse)
# Inputs: CSVs with columns: symbol, entrez (optional), logfc, adjpv
# Outputs: TF_NES.csv, PROGENy_NES.csv, heatmaps
# ============================================================

# 0) Packages -------------------------------------------------
safe_install <- function(pkgs){
  to_get <- pkgs[!(pkgs %in% rownames(installed.packages()))]
  if (length(to_get)) install.packages(to_get, dependencies = TRUE)
}

safe_bioc_install <- function(pkgs){
  to_get <- pkgs[!(pkgs %in% rownames(installed.packages()))]
  if (length(to_get)) BiocManager::install(to_get, ask = FALSE, update = FALSE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# For GitHub: do not auto-install unless you explicitly enable it
# Set env var: RNASEQ_TOOLKIT_INSTALL_DEPS=true
do_install <- identical(tolower(Sys.getenv("RNASEQ_TOOLKIT_INSTALL_DEPS")), "true")

cran_pkgs <- c("data.table","dplyr","readr","tibble","ggplot2","pheatmap","tidyr")
bioc_pkgs <- c("dorothea","viper","progeny","org.Mm.eg.db","AnnotationDbi")

if (do_install) {
  safe_install(cran_pkgs)
  safe_bioc_install(bioc_pkgs)
} else {
  missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_cran) || length(missing_bioc)) {
    stop(
      "Missing packages:\n",
      if (length(missing_cran)) paste0("CRAN: ", paste(missing_cran, collapse = ", "), "\n") else "",
      if (length(missing_bioc)) paste0("Bioc: ", paste(missing_bioc, collapse = ", "), "\n") else "",
      "Install them (or set RNASEQ_TOOLKIT_INSTALL_DEPS=true)."
    )
  }
}

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(tibble)
  library(ggplot2);    library(pheatmap); library(tidyr)
  library(dorothea);   library(viper);  library(progeny)
})

# 1) Inputs ---------------------------------------------------
# OPTIONAL: allow base_dir override from command line or env var
# Usage:
#   Rscript scripts/06_driver_inference_dorothea_progeny.R "path/to/project_dir"
args <- commandArgs(trailingOnly = TRUE)
base_dir_default <- "C:/Users/antho/OneDrive - Wayne State University/Wayne/Mor Lab/Projects/Fetal Immune System Project/Experiments/17. Benzene Brain/1 - Benzene Brain Initial Analysis"
base_dir <- if (length(args) >= 1) args[1] else Sys.getenv("RNASEQ_BASE_DIR", unset = base_dir_default)

paths <- list(
  Female_E12 = file.path(base_dir, "E12.5_Benzene_Brain_Female.csv"),
  Female_E17 = file.path(base_dir, "E17.5_Benzene_Brain_Female.csv"),
  Male_E12   = file.path(base_dir, "E12.5_Benzene_Brain_Male.csv"),
  Male_E17   = file.path(base_dir, "E17.5_Benzene_Brain_Male.csv")
)

outdir <- file.path(base_dir, "Driver_Inference")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# NEW: DEG gates (filters input rows before making score_mat)
# If your files are already "significant-only", set deg_fdr_cut <- 1 to disable.
deg_fdr_cut <- 0.05

# Optional per-list override
# deg_fdr_cut_by_list <- c(Female_E12=0.05, Female_E17=0.1, Male_E12=0.05, Male_E17=0.1)
deg_fdr_cut_by_list <- NULL

# Optional LFC gate
use_lfc_cut <- FALSE
lfc_cut <- 0.25

# 2) Read and build gene-by-list matrix ----------------------
read_deg <- function(fp, nm, deg_fdr_cut, use_lfc_cut, lfc_cut){
  if (!file.exists(fp)) stop("Missing input file: ", fp)

  x <- data.table::fread(fp, data.table = FALSE, showProgress = FALSE)
  names(x) <- tolower(names(x))
  need <- c("symbol","logfc","adjpv")
  miss <- setdiff(need, names(x))
  if (length(miss)) stop("Missing columns in ", basename(fp), ": ", paste(miss, collapse = ", "))

  x <- x %>%
    mutate(
      symbol = as.character(symbol),
      logfc  = suppressWarnings(as.numeric(logfc)),
      adjpv  = suppressWarnings(as.numeric(adjpv))
    ) %>%
    filter(!is.na(symbol) & nzchar(symbol)) %>%
    distinct(symbol, .keep_all = TRUE)

  # NEW: significance gate
  if (!is.null(deg_fdr_cut) && !is.na(deg_fdr_cut)) {
    x <- x %>% filter(!is.na(adjpv), adjpv <= deg_fdr_cut)
  }

  # NEW: optional LFC gate
  if (isTRUE(use_lfc_cut)) {
    x <- x %>% filter(!is.na(logfc), abs(logfc) >= lfc_cut)
  }

  x
}

diag_lines <- c()

degs <- lapply(names(paths), function(nm){
  this_cut <- deg_fdr_cut
  if (!is.null(deg_fdr_cut_by_list) && nm %in% names(deg_fdr_cut_by_list)) {
    this_cut <- deg_fdr_cut_by_list[[nm]]
  }

  df <- read_deg(paths[[nm]], nm, this_cut, use_lfc_cut, lfc_cut)

  diag_lines <- c(
    diag_lines,
    sprintf("[%s] deg_fdr_cut=%s use_lfc_cut=%s lfc_cut=%s n_genes=%d",
            nm, as.character(this_cut), as.character(use_lfc_cut), as.character(lfc_cut), nrow(df))
  )

  df
})
names(degs) <- names(paths)

# Combine symbols; fill missing with 0 (neutral) so VIPER/PROGENy can run
all_syms <- sort(unique(unlist(lapply(degs, function(d) d$symbol))))
diag_lines <- c(diag_lines, sprintf("[Union] n_unique_symbols=%d", length(all_syms)))

score_mat <- sapply(names(degs), function(nm){
  v <- setNames(degs[[nm]]$logfc, degs[[nm]]$symbol)
  out <- v[all_syms]; out[is.na(out)] <- 0
  out
})
rownames(score_mat) <- all_syms
# score_mat: rows = gene symbols, cols = Female_E12, Female_E17, Male_E12, Male_E17

# 3) DoRothEA + VIPER (mouse regulons, A-C) ------------------
data(dorothea_mm, package = "dorothea")
regulon_mm <- dorothea_mm %>% dplyr::filter(confidence %in% c("A","B","C"))
regulon_viper <- dorothea::df2regulon(regulon_mm)

tf_nes <- viper::viper(
  eset    = score_mat,
  regulon = regulon_viper,
  nes     = TRUE,
  minsize = 5,
  verbose = FALSE
)

tf_tbl <- as.data.frame(tf_nes) %>%
  tibble::rownames_to_column("TF") %>%
  dplyr::arrange(dplyr::desc(abs(rowMeans(dplyr::across(-TF), na.rm = TRUE))))
readr::write_csv(tf_tbl, file.path(outdir, "DoRothEA_TF_activity_NES.csv"))

top_tf <- head(tf_tbl$TF, n = min(40, nrow(tf_tbl)))
tf_mat <- as.matrix(tf_nes[top_tf, , drop = FALSE])
png(file.path(outdir, "DoRothEA_TF_activity_top40_heatmap.png"),
    width = 1200, height = 1000, res = 150)
pheatmap::pheatmap(tf_mat, cluster_rows = TRUE, cluster_cols = TRUE,
                   main = "DoRothEA TF activity (NES) - top 40 by magnitude")
dev.off()

# 4) PROGENy (Mouse) -----------------------------------------
perm_n <- 1000L

path_nes <- tryCatch(
  progeny::progeny(
    score_mat,
    scale     = TRUE,
    organism  = "Mouse",
    top       = 100,
    perm      = perm_n,
    z_scores  = TRUE
  ),
  error = function(e) {
    progeny::progeny(
      score_mat,
      scale       = TRUE,
      organism    = "Mouse",
      top         = 100,
      permutation = perm_n,
      z_scores    = TRUE
    )
  }
)

path_tbl <- as.data.frame(path_nes) %>%
  tibble::rownames_to_column("Pathway") %>%
  dplyr::arrange(dplyr::desc(abs(rowMeans(dplyr::across(-Pathway), na.rm = TRUE))))
readr::write_csv(path_tbl, file.path(outdir, "PROGENy_pathway_activity_NES.csv"))

top_pw <- head(path_tbl$Pathway, n = min(14, nrow(path_tbl)))
pw_mat <- as.matrix(path_nes[top_pw, , drop = FALSE])
png(file.path(outdir, "PROGENy_pathway_activity_top14_heatmap.png"),
    width = 1100, height = 800, res = 150)
pheatmap::pheatmap(
  pw_mat, cluster_rows = TRUE, cluster_cols = TRUE,
  main = "PROGENy pathway activity (NES) - top pathways"
)
dev.off()

# 5) Small summaries -----------------------------------------
summ_top <- function(mat, k = 10){
  as.data.frame(mat) %>%
    tibble::rownames_to_column("feature") %>%
    tidyr::pivot_longer(-feature, names_to = "list", values_to = "NES") %>%
    dplyr::group_by(list) %>%
    dplyr::arrange(dplyr::desc(abs(NES)), .by_group = TRUE) %>%
    dplyr::mutate(.rk = dplyr::row_number()) %>%
    dplyr::filter(.rk <= k) %>%
    dplyr::select(-.rk) %>%
    dplyr::arrange(list, dplyr::desc(abs(NES))) %>%
    dplyr::ungroup()
}

readr::write_csv(summ_top(tf_nes,   10), file.path(outdir, "Top10_TFs_per_list.csv"))
readr::write_csv(summ_top(path_nes, 10), file.path(outdir, "Top10_Pathways_per_list.csv"))

# 6) Diagnostics and session info ----------------------------
writeLines(diag_lines, file.path(outdir, "diagnostics.txt"))
writeLines(c(capture.output(sessionInfo())), file.path(outdir, "sessionInfo.txt"))
message("Done. Outputs in: ", normalizePath(outdir, winslash = "/", mustWork = FALSE))
