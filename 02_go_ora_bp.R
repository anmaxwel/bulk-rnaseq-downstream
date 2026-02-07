# ============================================================
# GO:BP ORA from "already significant" brain DEG lists
# Inputs: E12.5_Benzene_Brain_Female/Male, E17.5_Benzene_Brain_Female/Male
# Columns expected: symbol, entrez, logfc, adjpv
# Universe: all mouse Entrez IDs from org.Mm.eg.db (fallback, see note below)
# Outputs: per-file ALL BP tables/plots, immune-only tables/plots,
#          and within-sex overlap plots (E12.5 vs E17.5) for ALL and immune-only
# ============================================================

# ---------- 0) Packages ----------
safe_install <- function(pkgs){
  to_get <- pkgs[!(pkgs %in% rownames(installed.packages()))]
  if (length(to_get)) install.packages(to_get, dependencies = TRUE)
}

safe_bioc_install <- function(pkgs){
  # Minimal change: only install missing Bioc packages
  to_get <- pkgs[!(pkgs %in% rownames(installed.packages()))]
  if (length(to_get)) BiocManager::install(to_get, ask = FALSE, update = FALSE)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

safe_install(c("data.table","dplyr","stringr","ggplot2","forcats","readr"))
safe_bioc_install(c("clusterProfiler","org.Mm.eg.db","GO.db","enrichplot","AnnotationDbi"))

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(stringr)
  library(ggplot2);    library(forcats); library(readr)
  library(clusterProfiler); library(org.Mm.eg.db); library(GO.db)
  library(enrichplot); library(AnnotationDbi)
})

# ---------- 1) Inputs ----------
# NOTE: for GitHub, you will eventually want to replace this with relative paths
base_dir <- "C:/Users/antho/OneDrive - Wayne State University/Wayne/Mor Lab/Projects/Fetal Immune System Project/Experiments/17. Benzene Brain/1 - Benzene Brain Initial Analysis"

paths <- list(
  Female_E12 = file.path(base_dir, "E12.5_Benzene_Brain_Female.csv"),
  Male_E12   = file.path(base_dir, "E12.5_Benzene_Brain_Male.csv"),
  Female_E17 = file.path(base_dir, "E17.5_Benzene_Brain_Female.csv"),
  Male_E17   = file.path(base_dir, "E17.5_Benzene_Brain_Male.csv")
)

# Optional log2FC gates applied to the lists
lfc_cut <- c(Female_E12 = 0,
             Male_E12   = 0,
             Female_E17 = 0,
             Male_E17   = 0)

# NEW: gene-level significance gate (filters your DEG lists before ORA)
# If your DEG files are already "significant only", set this to 1 to effectively disable.
deg_fdr_cut <- 0.05

# Optional: per-list DEG FDR gate (overrides deg_fdr_cut if provided)
# Example:
# deg_fdr_cut_by_list <- c(Female_E12=0.05, Male_E12=0.05, Female_E17=0.1, Male_E17=0.1)
deg_fdr_cut_by_list <- NULL

# FDR threshold used when calling a term "enriched" for overlap plots
fdr_cut <- 0.05

outdir <- file.path(base_dir, "GO_BP_Results_Brain")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---------- 2) Helpers ----------
read_gene_csv <- function(fp){
  df <- suppressWarnings(data.table::fread(fp, data.table = FALSE))
  names(df) <- tolower(names(df))
  need <- c("symbol","logfc","adjpv")
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("Missing columns in ", fp, ": ", paste(miss, collapse = ", "))
  keep <- intersect(c("symbol","entrez","logfc","adjpv"), names(df))
  df[, keep, drop = FALSE]
}

# NEW: apply gene-level significance gate
apply_sig_gate <- function(df, fdr_cut){
  if (is.null(fdr_cut) || is.na(fdr_cut)) return(df)
  df %>%
    mutate(adjpv = as.numeric(adjpv)) %>%
    filter(!is.na(adjpv), adjpv <= fdr_cut)
}

ensure_entrez <- function(df){
  # If 'entrez' missing or mostly NA, map SYMBOL -> ENTREZID (mouse)
  need_map <- (!"entrez" %in% names(df)) || mean(is.na(df$entrez)) > 0.5
  if (need_map) {
    map <- clusterProfiler::bitr(df$symbol,
                                 fromType = "SYMBOL",
                                 toType   = "ENTREZID",
                                 OrgDb    = org.Mm.eg.db)
    df <- df %>% left_join(map, by = c("symbol" = "SYMBOL"))
    df$entrez <- df$ENTREZID
    df$ENTREZID <- NULL
  }
  df %>% mutate(entrez = as.character(entrez)) %>% filter(!is.na(entrez))
}

apply_lfc_gate <- function(df, gate){
  df %>%
    mutate(logfc = as.numeric(logfc)) %>%
    filter(!is.na(logfc), abs(logfc) >= gate)
}

get_immune_bp_ids <- function(){
  # Core immune and inflammation roots
  roots <- c(
    "GO:0002376", # immune system process
    "GO:0006954", # inflammatory response
    "GO:0006952", # defense response
    "GO:0006955", # immune response
    "GO:0045087", # innate immune response
    "GO:0002682"  # regulation of immune system process
  )
  lst <- AnnotationDbi::as.list(GO.db::GOBPOFFSPRING)
  kids <- unique(unlist(lst[roots], use.names = FALSE))
  unique(c(roots, kids))
}

run_enrich_bp <- function(entrez_vec, universe_vec,
                          minGSSize = 5, maxGSSize = 2000){
  if (length(entrez_vec) < minGSSize) return(NULL)
  suppressMessages(
    enrichGO(
      gene          = unique(entrez_vec),
      universe      = unique(universe_vec),
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 1.0,   # do not prefilter at raw p
      qvalueCutoff  = 1.0,
      minGSSize     = minGSSize,
      maxGSSize     = maxGSSize,
      readable      = TRUE
    )
  )
}

er_to_df <- function(er, tag){
  if (is.null(er)) return(tibble())
  res <- as.data.frame(er)
  if (!nrow(res)) return(tibble())
  gr <- suppressWarnings(sapply(strsplit(res$GeneRatio, "/"), function(x){
    as.numeric(x[1]) / as.numeric(x[2])
  }))
  res %>%
    mutate(GeneRatioNum = as.numeric(gr),
           List = tag) %>%
    arrange(p.adjust)
}

safe_dotplot_df <- function(df, title_text, file_out, top_n = 20){
  if (is.null(df) || !nrow(df)) {
    png(file_out, width = 1200, height = 800, res = 150)
    plot.new(); title(main = paste0(title_text, "\n(No enriched terms)"))
    dev.off(); return(invisible(NULL))
  }
  top <- df %>%
    slice_head(n = min(top_n, nrow(.))) %>%
    mutate(Description = forcats::fct_reorder(Description, GeneRatioNum))

  p <- ggplot(top, aes(x = GeneRatioNum, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_size_continuous(range = c(4, 14)) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.05))) +
    scale_color_viridis_c(
      option = "mako",
      begin = 0.25,
      end = 0.9,
      direction = -1,
      guide = guide_colorbar(title = "FDR (p.adjust)")
    ) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
    labs(x = "Gene ratio", y = NULL, title = title_text) +
    theme_classic(base_size = 18) +
    theme(
      axis.text.x  = element_text(face = "bold", size = 30),
      axis.title.x = element_text(face = "bold", size = 30),
      axis.text.y  = element_text(face = "bold", size = 38),
      plot.title   = element_text(face = "bold", size = 28),
      legend.title = element_text(size = 30, face = "bold"),
      legend.text  = element_text(size = 26),
      panel.grid.major.x = element_line(color = "grey80", linewidth = 0.3),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
      panel.grid.minor   = element_blank()
    )
  ggsave(file_out, p, width = 22, height = 14, dpi = 300)
}

filter_to_immune_ids <- function(df, immune_ids){
  if (!nrow(df)) return(df)
  df %>% filter(ID %in% immune_ids)
}

pairwise_overlap_plot <- function(dfA, dfB, labelA, labelB, out_png,
                                  sig_col = c("p.adjust","pvalue"),
                                  sig_thresh = 0.05,
                                  top_n = 30) {

  sig_col <- match.arg(sig_col)

  draw_empty <- function(msg){
    png(out_png, width = 1200, height = 800, res = 150)
    plot.new(); title(main = msg)
    dev.off(); invisible(NULL)
  }

  has_col <- function(x, col) is.data.frame(x) && (col %in% names(x)) && nrow(x) > 0
  if (!has_col(dfA, sig_col) || !has_col(dfB, sig_col))
    return(draw_empty(paste0("Overlap: ", labelA, " vs ", labelB, "\n(No data)")))

  A_sig <- dfA %>% dplyr::filter(!is.na(.data[[sig_col]]), .data[[sig_col]] <= sig_thresh)
  B_sig <- dfB %>% dplyr::filter(!is.na(.data[[sig_col]]), .data[[sig_col]] <= sig_thresh)

  if (!nrow(A_sig) || !nrow(B_sig))
    return(draw_empty(paste0("Overlap: ", labelA, " vs ", labelB,
                             "\n(No significant terms at ", sig_col, " <= ", sig_thresh, ")")))

  inter_ids <- intersect(A_sig$ID, B_sig$ID)
  if (!length(inter_ids))
    return(draw_empty(paste0("Overlap: ", labelA, " vs ", labelB,
                             "\n(No shared significant terms at ", sig_col, " <= ", sig_thresh, ")")))

  A_sub <- A_sig %>% filter(ID %in% inter_ids) %>%
    transmute(ID, Description, List = labelA, neglog10 = -log10(.data[[sig_col]]))
  B_sub <- B_sig %>% filter(ID %in% inter_ids) %>%
    transmute(ID, Description, List = labelB, neglog10 = -log10(.data[[sig_col]]))
  long <- bind_rows(A_sub, B_sub)

  # order by labelB only
  ord_desc <- B_sub %>%
    arrange(desc(neglog10)) %>%
    distinct(Description, .keep_all = TRUE) %>%
    slice_head(n = top_n) %>%
    pull(Description)

  long <- long %>%
    filter(Description %in% ord_desc) %>%
    mutate(Description = factor(Description, levels = rev(ord_desc)))

  p <- ggplot(long, aes(x = List, y = Description)) +
    geom_tile(aes(fill = neglog10),
              width  = 0.8,
              height = 0.8) +
    scale_fill_gradient2(
      low = "purple",
      mid = "red",
      high = "orange",
      midpoint = median(long$neglog10),
      name = paste0("-log10(", sig_col, ")")
    ) +
    labs(
      x = NULL, y = NULL,
      title = paste0("Overlap of enriched BP terms: ", labelA, " vs ", labelB)
    ) +
    theme_minimal(base_size = 20) +
    theme(
      axis.text.x = element_text(size = 35, face = "bold"),
      axis.text.y = element_text(size = 35, face = "bold"),
      plot.title  = element_text(face = "bold")
    )

  ggsave(out_png, p, width = 35, height = 12, dpi = 300)
  invisible(NULL)
}

# ---------- 3) Prepare global universe ----------
# NOTE: Better universe is "all genes tested in DE", but keeping your current behavior.
universe_all <- unique(AnnotationDbi::keys(org.Mm.eg.db, keytype = "ENTREZID"))

# ---------- 4) Run enrichment ----------
immune_ids <- get_immune_bp_ids()

all_tables_all    <- list()
all_tables_immune <- list()
diag_lines <- c()

for (nm in names(paths)) {
  message("Processing: ", nm)

  # NEW: pick per-list DEG FDR cut if provided
  this_deg_fdr <- deg_fdr_cut
  if (!is.null(deg_fdr_cut_by_list) && nm %in% names(deg_fdr_cut_by_list)) {
    this_deg_fdr <- deg_fdr_cut_by_list[[nm]]
  }

  df_raw <- read_gene_csv(paths[[nm]])

  # NEW: apply DEG significance first
  df_sig <- apply_sig_gate(df_raw, this_deg_fdr)

  # then map entrez and apply LFC gate
  df <- df_sig %>% ensure_entrez() %>% apply_lfc_gate(lfc_cut[[nm]])

  genes <- unique(df$entrez)

  diag_lines <- c(diag_lines,
                  sprintf("[%s] deg_fdr_cut = %s", nm, as.character(this_deg_fdr)),
                  sprintf("[%s] n_rows_input = %d", nm, nrow(df_raw)),
                  sprintf("[%s] n_rows_after_deg_fdr = %d", nm, nrow(df_sig)),
                  sprintf("[%s] n_genes_after_entrez_and_lfc = %d", nm, length(genes)))

  er  <- run_enrich_bp(genes, universe_all)
  tbl <- er_to_df(er, nm)

  # save tables
  out_csv_all <- file.path(outdir, paste0("ORA_BP_", nm, "_ALL.csv"))
  if (nrow(tbl)) readr::write_csv(tbl, out_csv_all) else writeLines("", out_csv_all)

  # plot all terms
  safe_dotplot_df(tbl,
                  paste0(nm, " - GO:BP enriched (ALL terms)"),
                  file.path(outdir, paste0("ORA_BP_", nm, "_ALL_dotplot.png")),
                  top_n = 10)

  # immune-only subset
  tbl_imm <- filter_to_immune_ids(tbl, immune_ids)
  out_csv_imm <- file.path(outdir, paste0("ORA_BP_", nm, "_IMMUNE.csv"))
  if (nrow(tbl_imm)) readr::write_csv(tbl_imm, out_csv_imm) else writeLines("", out_csv_imm)

  safe_dotplot_df(tbl_imm,
                  paste0(nm, " - GO:BP enriched (immune-only)"),
                  file.path(outdir, paste0("ORA_BP_", nm, "_IMMUNE_dotplot.png")),
                  top_n = 10)

  all_tables_all[[nm]]    <- tbl
  all_tables_immune[[nm]] <- tbl_imm

  if (is.null(er) || !nrow(tbl)) {
    diag_lines <- c(diag_lines, sprintf("[%s] No enriched terms returned", nm))
  } else {
    sig_n <- sum(tbl$p.adjust <= fdr_cut, na.rm = TRUE)
    diag_lines <- c(diag_lines, sprintf("[%s] Enriched terms (FDR<=%.2f): %d",
                                        nm, fdr_cut, sig_n))
  }
}

# ---------- 5) Overlap plots: E12.5 vs E17.5 within sex ----------
# Here you are using raw pvalue threshold for overlap.
# Change sig_thresh as desired.

# Female - ALL
pairwise_overlap_plot(
  all_tables_all$Female_E12, all_tables_all$Female_E17,
  labelA = "Female E12.5", labelB = "Female E17.5",
  out_png = file.path(outdir, "Overlap_BP_Female_ALL_rawP.png"),
  sig_col = "pvalue", sig_thresh = 0.05, top_n = 10
)

# Female - IMMUNE
pairwise_overlap_plot(
  all_tables_immune$Female_E12, all_tables_immune$Female_E17,
  labelA = "Female E12.5 (immune)", labelB = "Female E17.5 (immune)",
  out_png = file.path(outdir, "Overlap_BP_Female_IMMUNE_rawP.png"),
  sig_col = "pvalue", sig_thresh = 0.05, top_n = 10
)

# Male - ALL
pairwise_overlap_plot(
  all_tables_all$Male_E12, all_tables_all$Male_E17,
  labelA = "Male E12.5", labelB = "Male E17.5",
  out_png = file.path(outdir, "Overlap_BP_Male_ALL_rawP.png"),
  sig_col = "pvalue", sig_thresh = 0.05, top_n = 10
)

# Male - IMMUNE
pairwise_overlap_plot(
  all_tables_immune$Male_E12, all_tables_immune$Male_E17,
  labelA = "Male E12.5 (immune)", labelB = "Male E17.5 (immune)",
  out_png = file.path(outdir, "Overlap_BP_Male_IMMUNE_rawP.png"),
  sig_col = "pvalue", sig_thresh = 0.05, top_n = 10
)

# ---------- 6) Diagnostics ----------
writeLines(diag_lines, file.path(outdir, "diagnostics.txt"))
message("Done. See outputs in: ", normalizePath(outdir))

########################################################################
# Mirrorbar overlap plot (kept, but filenames corrected to match sig_col)
########################################################################

pairwise_overlap_mirrorbar <- function(dfA, dfB, labelA, labelB, out_png,
                                       sig_col   = c("p.adjust","pvalue"),
                                       sig_thresh = 0.05,
                                       top_n      = 20) {
  sig_col <- match.arg(sig_col)

  draw_empty <- function(msg){
    png(out_png, width = 1200, height = 800, res = 150)
    plot.new(); title(main = msg)
    dev.off(); invisible(NULL)
  }

  has_col <- function(x, col) {
    is.data.frame(x) && (col %in% names(x)) && nrow(x) > 0
  }

  if (!has_col(dfA, sig_col) || !has_col(dfB, sig_col)) {
    return(draw_empty(paste0("Overlap: ", labelA, " vs ", labelB, "\n(No data)")))
  }
  if (!("Count" %in% names(dfA)) || !("Count" %in% names(dfB))) {
    return(draw_empty(paste0("Overlap: ", labelA, " vs ", labelB,
                             "\n(Count column not found in enrichment tables)")))
  }

  A_sig <- dfA %>%
    dplyr::filter(!is.na(.data[[sig_col]]), .data[[sig_col]] <= sig_thresh)
  B_sig <- dfB %>%
    dplyr::filter(!is.na(.data[[sig_col]]), .data[[sig_col]] <= sig_thresh)

  if (!nrow(A_sig) || !nrow(B_sig)) {
    return(draw_empty(paste0("Overlap: ", labelA, " vs ", labelB,
                             "\n(No significant terms at ", sig_col, " <= ", sig_thresh, ")")))
  }

  inter_ids <- intersect(A_sig$ID, B_sig$ID)
  if (!length(inter_ids)) {
    return(draw_empty(paste0("Overlap: ", labelA, " vs ", labelB,
                             "\n(No shared significant terms at ", sig_col, " <= ", sig_thresh, ")")))
  }

  A_sub <- A_sig %>%
    dplyr::filter(ID %in% inter_ids) %>%
    dplyr::transmute(
      ID,
      Description,
      value = -log10(.data[[sig_col]]),
      List  = labelA,
      Count = Count
    )

  B_sub <- B_sig %>%
    dplyr::filter(ID %in% inter_ids) %>%
    dplyr::transmute(
      ID,
      Description,
      value = -log10(.data[[sig_col]]),
      List  = labelB,
      Count = Count
    )

  long <- dplyr::bind_rows(A_sub, B_sub)
  if (!nrow(long)) {
    return(draw_empty(paste0("Overlap: ", labelA, " vs ", labelB, "\n(Join produced empty set)")))
  }

  ord_desc <- long %>%
    dplyr::group_by(Description) %>%
    dplyr::summarise(max_val = max(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(max_val)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(Description)

  long <- long %>%
    dplyr::filter(Description %in% ord_desc) %>%
    dplyr::mutate(
      Description = factor(Description, levels = rev(ord_desc)),
      value_signed = dplyr::if_else(List == labelA, -value, value),
      x_label = dplyr::if_else(List == labelA,
                               value_signed - 0.1,
                               value_signed + 0.1),
      hjust_label = dplyr::if_else(List == labelA, 1.0, 0.0)
    )

  abs_lab <- function(x) sprintf("%.1f", abs(x))

  p <- ggplot(long, aes(x = value_signed, y = Description, fill = List)) +
    geom_col(width = 0.7, alpha = 0.9) +
    geom_vline(xintercept = 0, color = "black") +
    geom_text(
      aes(x = x_label, label = Count, hjust = hjust_label),
      size = 3
    ) +
    scale_x_continuous(labels = abs_lab) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
    labs(
      x = paste0("-log10(", sig_col, ")"),
      y = NULL,
      title = paste0("Shared enriched GO:BP terms: ", labelA, " vs ", labelB)
    ) +
    theme_minimal(base_size = 16) +
    theme(
      axis.title.x = element_text(size = 20, face = "bold"),
      axis.text.x  = element_text(size = 20, face = "bold"),
      axis.text.y  = element_text(size = 20, face = "bold"),
      plot.title   = element_text(size = 20, face = "bold"),
      legend.title = element_blank()
    )

  ggsave(out_png, p, width = 20, height = 10, dpi = 300)
  invisible(NULL)
}

# Mirrorbar plots using FDR (p.adjust) for overlap
pairwise_overlap_mirrorbar(
  all_tables_all$Female_E12, all_tables_all$Female_E17,
  labelA = "Female E12.5", labelB = "Female E17.5",
  out_png = file.path(outdir, "Overlap_BP_Female_ALL_FDR_mirrorbar.png"),
  sig_col = "p.adjust", sig_thresh = 0.05, top_n = 10
)

pairwise_overlap_mirrorbar(
  all_tables_immune$Female_E12, all_tables_immune$Female_E17,
  labelA = "Female E12.5 (immune)", labelB = "Female E17.5 (immune)",
  out_png = file.path(outdir, "Overlap_BP_Female_IMMUNE_FDR_mirrorbar.png"),
  sig_col = "p.adjust", sig_thresh = 0.05, top_n = 10
)

pairwise_overlap_mirrorbar(
  all_tables_all$Male_E12, all_tables_all$Male_E17,
  labelA = "Male E12.5", labelB = "Male E17.5",
  out_png = file.path(outdir, "Overlap_BP_Male_ALL_FDR_mirrorbar.png"),
  sig_col = "p.adjust", sig_thresh = 0.05, top_n = 10
)

pairwise_overlap_mirrorbar(
  all_tables_immune$Male_E12, all_tables_immune$Male_E17,
  labelA = "Male E12.5 (immune)", labelB = "Male E17.5 (immune)",
  out_png = file.path(outdir, "Overlap_BP_Male_IMMUNE_FDR_mirrorbar.png"),
  sig_col = "p.adjust", sig_thresh = 0.05, top_n = 10
)
