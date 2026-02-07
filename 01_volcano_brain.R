# ============================================================
# Volcano plots for Brain Bulk RNA-seq CSVs
# Required columns (aliases OK): symbol, logfc, adjpv
# Outputs: PNG plots + TSV per plot
# ============================================================

# ---------- 0) Base paths ----------
# OPTIONAL: allow base_dir override from command line or env var
# Usage:
#   Rscript scripts/01_volcano_brain.R "path/to/project_dir"
args <- commandArgs(trailingOnly = TRUE)

base_dir_default <- "C:/Users/antho/OneDrive - Wayne State University/Wayne/Mor Lab/Projects/Fetal Immune System Project/Experiments/17. Benzene Brain/1 - Benzene Brain Initial Analysis"
base_dir <- if (length(args) >= 1) args[1] else Sys.getenv("RNASEQ_BASE_DIR", unset = base_dir_default)

# Parent directory where your GO results are going
outdir <- file.path(base_dir, "GO_BP_Results_Brain")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Brain DEG CSVs
paths <- list(
  Female_E12 = file.path(base_dir, "E12.5_Benzene_Brain_Female.csv"),
  Male_E12   = file.path(base_dir, "E12.5_Benzene_Brain_Male.csv"),
  Female_E17 = file.path(base_dir, "E17.5_Benzene_Brain_Female.csv"),
  Male_E17   = file.path(base_dir, "E17.5_Benzene_Brain_Male.csv")
)

# All volcano outputs live here
out_dir <- file.path(outdir, "volcano_brain")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Volcano plots
# ============================================================

# ---------- Packages ----------
safe_install <- function(pkgs){
  to_get <- pkgs[!(pkgs %in% rownames(installed.packages()))]
  if (length(to_get)) install.packages(to_get, dependencies = TRUE)
}

# For GitHub: do not auto-install unless you explicitly enable it
# Set env var: RNASEQ_TOOLKIT_INSTALL_DEPS=true
do_install <- identical(tolower(Sys.getenv("RNASEQ_TOOLKIT_INSTALL_DEPS")), "true")

pkgs_needed <- c("data.table","dplyr","stringr","ggplot2","ggrepel","scales","tidyr")
if (do_install) {
  safe_install(pkgs_needed)
} else {
  missing <- pkgs_needed[!vapply(pkgs_needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         "\nInstall them (or set RNASEQ_TOOLKIT_INSTALL_DEPS=true).")
  }
}

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(stringr)
  library(ggplot2); library(ggrepel); library(scales)
})

# ---------- Parameters ----------
adjpv_cut   <- 0.05   # adjusted p-value threshold
use_lfc_cut <- FALSE  # set TRUE to also require abs(logFC) >= lfc_cut
lfc_cut     <- 0.25
label_top_n <- 0      # number of top significant genes to label

# ---------- Helpers ----------
read_robust <- function(path){
  stopifnot(file.exists(path))
  for (sep in c(",", "\t", ";")) {
    dt <- try(
      fread(path, sep = sep, header = TRUE,
            data.table = FALSE, showProgress = FALSE),
      silent = TRUE
    )
    if (inherits(dt, "try-error")) next
    if (!is.data.frame(dt) || ncol(dt) < 3) next

    nm <- tolower(names(dt))
    nm <- gsub("\\s+","_", nm)
    nm <- gsub("\\.+","_", nm)
    nm <- gsub('"',"", nm, fixed = TRUE)
    names(dt) <- nm

    if (!"symbol" %in% names(dt)) {
      cand_sym <- c("gene","gene_symbol","hgnc_symbol","id","feature","name")
      hit <- cand_sym[cand_sym %in% names(dt)]
      if (length(hit)) dt$symbol <- dt[[hit[1]]]
    }
    if (!"logfc" %in% names(dt)) {
      cand_lfc <- c("log2fc","avg_log2fc","avg_logfc",
                    "log_fold_change","lfc","logfc_rna","logfc_sc")
      hit <- cand_lfc[cand_lfc %in% names(dt)]
      if (length(hit)) dt$logfc <- dt[[hit[1]]]
    }
    if (!"adjpv" %in% names(dt)) {
      cand_p <- c("padj","p_adj","adj_p","adj_p_val",
                  "fdr","qvalue","q_value","bh_adj_p")
      hit <- cand_p[cand_p %in% names(dt)]
      if (length(hit)) dt$adjpv <- dt[[hit[1]]]
    }

    req <- c("symbol","logfc","adjpv")
    if (!all(req %in% names(dt))) next

    dt$symbol <- as.character(dt$symbol)
    dt$logfc  <- suppressWarnings(as.numeric(dt$logfc))
    dt$adjpv  <- suppressWarnings(as.numeric(dt$adjpv))

    # basic clean
    is_char <- vapply(dt, is.character, logical(1))
    if (any(is_char)) {
      dt[is_char] <- lapply(dt[is_char], function(x){
        x <- gsub('^"+|"+$', "", x)
        x <- gsub("\\s+", " ", x)
        trimws(x)
      })
    }

    dt <- dplyr::distinct(dt, symbol, .keep_all = TRUE)
    return(dt[, c("symbol","logfc","adjpv",
                  setdiff(names(dt), c("symbol","logfc","adjpv"))),
              drop = FALSE])
  }
  stop("Failed to read a valid table: ", path)
}

callout_png <- function(file, title, msg, w=8, h=6){
  png(file, width = w*100, height = h*100)
  plot.new()
  text(0.5, 0.65, labels = title, cex = 1.3, font = 2)
  txt <- paste(strwrap(msg, width = 95), collapse = "\n")
  text(0.5, 0.45, labels = txt, cex = 0.95)
  dev.off()
}

volcano_plot <- function(df, label = "Volcano", adjpv_cut = 0.05,
                         use_lfc_cut = FALSE, lfc_cut = 0.25, label_top_n = 20) {
  stopifnot(is.data.frame(df))

  df <- df %>%
    mutate(
      neglog10 = -log10(pmax(adjpv, .Machine$double.xmin)),
      sig_adj  = is.finite(adjpv) & adjpv <= adjpv_cut,
      sig_lfc  = if (use_lfc_cut) is.finite(logfc) & abs(logfc) >= lfc_cut else TRUE,
      sig      = sig_adj & sig_lfc,
      Dir      = case_when(
        sig & logfc >  0 ~ "Up",
        sig & logfc <  0 ~ "Down",
        TRUE             ~ "NS"
      )
    )

  lab_df <- df %>%
    filter(sig) %>%
    arrange(desc(neglog10), desc(abs(logfc))) %>%
    slice_head(n = label_top_n)

  gg <- ggplot(df, aes(x = logfc, y = neglog10)) +
    geom_point(aes(color = Dir), alpha = 0.7, size = 2.6, stroke = 0) +
    scale_color_manual(values = c(Down = "#2c7bb6", NS = "grey75", Up = "#d7191c")) +
    geom_hline(yintercept = -log10(adjpv_cut), linetype = "dashed") +
    labs(
      title = paste0("Volcano: ", label),
      x = "log2 fold change",
      y = "-log10(adj p-value)",
      color = "Direction"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.title.x = element_text(size = 22, face = "bold"),
      axis.title.y = element_text(size = 22, face = "bold"),
      axis.text.x  = element_text(size = 22, face = "bold"),
      axis.text.y  = element_text(size = 22, face = "bold"),
      plot.title   = element_text(size = 22, face = "bold")
    )

  if (use_lfc_cut) {
    gg <- gg + geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed")
  }

  gg <- gg +
    ggrepel::geom_text_repel(
      data = lab_df,
      aes(label = symbol),
      max.overlaps = 1000,
      size = 3,
      box.padding = 0.35,
      point.padding = 0.2,
      min.segment.length = 0
    )

  list(plot = gg, plotted_df = df)
}

save_volcano_safe <- function(df, label, out_png, out_tsv,
                              adjpv_cut, use_lfc_cut, lfc_cut, label_top_n){
  tryCatch({
    res <- volcano_plot(df, label, adjpv_cut, use_lfc_cut, lfc_cut, label_top_n)
    ggsave(out_png, plot = res$plot, width = 7.5, height = 6.5, dpi = 200)
    fwrite(res$plotted_df, out_tsv, sep = "\t")
  }, error = function(e){
    callout_png(out_png, paste0("Volcano: ", label), e$message)
  })
}

# ---------- Run for each Sex x Stage ----------
for (nm in names(paths)) {
  in_file <- paths[[nm]]
  if (!file.exists(in_file)) {
    message("[Skip] Missing input file: ", in_file)
    next
  }

  label   <- nm
  base    <- paste0("volcano_", tolower(nm))
  png_out <- file.path(out_dir, paste0(base, ".png"))
  tsv_out <- file.path(out_dir, paste0(base, ".tsv"))

  df <- read_robust(in_file)

  save_volcano_safe(
    df, label, png_out, tsv_out,
    adjpv_cut   = adjpv_cut,
    use_lfc_cut = use_lfc_cut,
    lfc_cut     = lfc_cut,
    label_top_n = label_top_n
  )
}

message("Brain volcano outputs written to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))

# ============================================================
# Count Up/Down consistent with the volcano plots
# ============================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr)
})

labels <- names(paths)

count_list <- lapply(labels, function(lbl){
  base  <- paste0("volcano_", tolower(lbl))
  tsv   <- file.path(out_dir, paste0(base, ".tsv"))

  if (file.exists(tsv)) {
    df <- data.table::fread(tsv, data.table = FALSE)
  } else {
    if (!file.exists(paths[[lbl]])) return(NULL)
    df <- read_robust(paths[[lbl]])
    df <- df %>%
      dplyr::mutate(
        neglog10 = -log10(pmax(adjpv, .Machine$double.xmin)),
        sig_adj  = is.finite(adjpv) & adjpv <= adjpv_cut,
        sig_lfc  = if (use_lfc_cut) is.finite(logfc) & abs(logfc) >= lfc_cut else TRUE,
        sig      = sig_adj & sig_lfc,
        Dir      = dplyr::case_when(
          sig & logfc >  0 ~ "Up",
          sig & logfc <  0 ~ "Down",
          TRUE             ~ "NS"
        )
      )
  }

  parts <- strsplit(lbl, "_")[[1]]
  sex   <- parts[1]
  stage <- parts[2]

  if ("symbol" %in% names(df)) {
    df <- dplyr::distinct(df, symbol, .keep_all = TRUE)
    total_n <- dplyr::n_distinct(df$symbol)
  } else {
    total_n <- nrow(df)
  }

  df <- df %>%
    dplyr::mutate(Dir = factor(Dir, levels = c("Up","Down","NS")))

  tab <- df %>%
    dplyr::count(Dir, name = "n", .drop = FALSE) %>%
    tidyr::pivot_wider(names_from = Dir, values_from = n, values_fill = 0) %>%
    dplyr::mutate(
      Sex   = sex,
      Stage = stage,
      Total = total_n,
      Sig   = dplyr::coalesce(Up, 0L) + dplyr::coalesce(Down, 0L)
    ) %>%
    dplyr::select(Sex, Stage, dplyr::any_of(c("Up","Down","NS")), Sig, Total)

  # FIXED: use tab columns, not objects
  if (nrow(tab)) {
    sum_udn <- tab$Up + tab$Down + tab$NS
    if (!is.na(sum_udn) && sum_udn != tab$Total) {
      message("[Warn] Totals mismatch for ", sex, " ", stage,
              " -> Up+Down+NS=", sum_udn, " vs Total=", tab$Total)
    }
  }

  tab
})

counts <- dplyr::bind_rows(count_list)
out_file <- file.path(out_dir, "volcano_counts_brain.tsv")
data.table::fwrite(counts, out_file, sep = "\t")
print(counts)
message("Brain counts written to: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
