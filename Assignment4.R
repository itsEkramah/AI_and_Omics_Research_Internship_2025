###############################################################################
# GSE15932_preprocess_all_in_one.R
# Full preprocessing pipeline for Affymetrix CELs (GPL570) — GSE15932
# Produces: QC (before/after) boxplots & PCA PDFs, filtered expression matrix,
#           outlier lists, summary CSV/text, sample group labels.
#
# Edit CONFIG at top, then run the script in R / RStudio.
###############################################################################

## -------------------- CONFIG --------------------
tar_file_path <- "/mnt/E42C87742C874092/AI_and_Omics_Research_Internship_2025/assignment_3/GSE15932_RAW.tar"
work_dir <- file.path(getwd(), "GSE15932_work")
out_dir  <- file.path(work_dir, "outputs")
phenotype_csv <- NULL  # set path to phenotype CSV if you have one (columns: sample, group)
min_sample_fraction_expressed <- 0.25  # fraction of samples that must show expression
expression_threshold_linear <- 100      # linear intensity threshold (pre-log2) used as filter baseline
outlier_pca_z_cutoff <- 3              # >3 SD from center = outlier
save_rds <- TRUE
###############################################################################

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)

## -------------------- helper to install/load --------------------
install_load <- function(pkgs, bioc = FALSE) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install(p, ask = FALSE, update = FALSE)
      } else {
        install.packages(p, repos = "https://cloud.r-project.org")
      }
    }
    library(p, character.only = TRUE)
  }
}
# CRAN
install_load(c("ggplot2", "data.table", "matrixStats"), bioc = FALSE)
# Bioconductor
install_load(c("affy", "affyPLM", "oligo", "limma", "arrayQualityMetrics", "GEOquery", "Biobase"), bioc = TRUE)

## -------------------- Extract tar --------------------
message("Extracting tar: ", tar_file_path)
if (!file.exists(tar_file_path)) stop("Tar file not found: ", tar_file_path)
untar(tar_file_path, exdir = work_dir)
message("Extraction done. Work dir: ", work_dir)

all_files <- list.files(work_dir, recursive = TRUE, full.names = TRUE)
cel_files <- grep("\\.CEL$|\\.cel$", all_files, ignore.case = TRUE, value = TRUE)
if (length(cel_files) == 0) {
  # check for .CEL.gz -> gunzip
  gz <- grep("\\.CEL\\.gz$|\\.cel\\.gz$", all_files, ignore.case = TRUE, value = TRUE)
  if (length(gz) > 0) {
    message("Found .CEL.gz files, decompressing ...")
    for (g in gz) if (file.exists(g)) R.utils::gunzip(g, remove = FALSE, overwrite = FALSE)
    all_files <- list.files(work_dir, recursive = TRUE, full.names = TRUE)
    cel_files <- grep("\\.CEL$|\\.cel$", all_files, ignore.case = TRUE, value = TRUE)
  }
}
if (length(cel_files) == 0) stop("No CEL files found in the tar extraction. Check the tar file.")

cel_dir <- unique(dirname(cel_files))[1]
message("Reading CELs from: ", cel_dir)

## -------------------- Read raw CELs --------------------
# For GPL570, affy::ReadAffy is appropriate
affybatch <- tryCatch({
  affy::ReadAffy(celfile.path = cel_dir)
}, error = function(e) {
  stop("Error reading CELs with affy::ReadAffy(): ", e$message)
})
message("CEL files read. Samples: ", length(sampleNames(affybatch)))

## -------------------- Try to load phenotype (GEO or CSV) --------------------
pheno <- NULL
if (!is.null(phenotype_csv) && file.exists(phenotype_csv)) {
  pheno <- read.csv(phenotype_csv, stringsAsFactors = FALSE, check.names = FALSE)
  message("Phenotype loaded from CSV: ", phenotype_csv)
} else {
  # try GEOquery to fetch metadata (needs internet)
  message("Attempting to fetch phenotype metadata from GEO (GSE15932)...")
  try({
    gse <- GEOquery::getGEO("GSE15932", GSEMatrix = FALSE)
    if (!is.null(gse)) {
      gsml <- GEOquery::GSMList(gse)
      smpl <- names(gsml)
      pheno <- data.frame(sample = smpl, title = sapply(gsml, function(x) Meta(x)$title),
                          characteristics = sapply(gsml, function(x) {
                            ch <- tryCatch(GEOquery::Table(x)$characteristics_ch1, error = function(e) NA)
                            if (is.na(ch)) return(NA)
                            if (length(ch) > 1) paste(ch, collapse=";") else ch
                          }), stringsAsFactors = FALSE)
      message("GEO phenotype metadata retrieved (may take a while).")
    }
  }, silent = TRUE)
}
if (is.null(pheno)) {
  message("No phenotype available. Proceeding with placeholder 'Unknown' labels.")
}

## -------------------- QC BEFORE normalization --------------------
message("QC BEFORE normalization ...")
# Boxplot of raw intensities (log2 scale)
raw_int_matrix <- tryCatch({
  # use affy::exprs on a simple summary? AffyBatch doesn't provide exprs directly until normalization,
  # but we can use affy::intensity to get probe-level intensities; to avoid massive memory use we compute probe medians per array
  pm <- affy::intensity(affybatch)  # may be large; returns matrix probes x arrays
  # Convert to log2
  log2(pm + 1)
}, error = function(e) {
  message("Could not extract PM intensity matrix (will fallback to per-array summaries). Error: ", e$message)
  NULL
})

pdf(file.path(out_dir, "boxplot_raw.pdf"), width = 10, height = 6)
if (!is.null(raw_int_matrix)) {
  boxplot(as.data.frame(raw_int_matrix), main = "Raw probe intensities (log2) - BEFORE normalization", las = 2)
} else {
  # fallback: boxplot of median intensities per array
  medians <- sapply(sampleNames(affybatch), function(s) median(affy::intensity(affybatch)[, s], na.rm = TRUE))
  boxplot(medians, main = "Raw array median intensities - BEFORE normalization")
}
dev.off()

# Density plot (raw)
pdf(file.path(out_dir, "density_raw.pdf"), width = 10, height = 6)
if (!is.null(raw_int_matrix)) {
  den <- apply(raw_int_matrix, 2, function(x) density(x, na.rm = TRUE))
  plot(den[[1]], main = "Raw density plots (log2)", xlab = "log2 intensity")
  for (i in 2:length(den)) lines(den[[i]])
} else {
  plot.new()
  title("Raw densities not available")
}
dev.off()

## -------------------- Outlier detection BEFORE normalization via PCA --------------------
outliers_before <- character(0)
if (!is.null(raw_int_matrix) && nrow(raw_int_matrix) > 50) {
  message("Running PCA on a subset of probes for speed (random 5000 probes) ...")
  set.seed(42)
  idx <- sample(nrow(raw_int_matrix), min(5000, nrow(raw_int_matrix)))
  pca_raw <- prcomp(t(raw_int_matrix[idx, , drop = FALSE]), scale. = TRUE)
  pc1 <- pca_raw$x[,1]; pc2 <- pca_raw$x[,2]
  distc <- sqrt((pc1 - mean(pc1))^2 + (pc2 - mean(pc2))^2)
  z <- scale(distc)
  outliers_before <- colnames(raw_int_matrix)[which(z > outlier_pca_z_cutoff)]
} else {
  message("Raw probe matrix not available or too small for PCA-based outliers; skipping PCA outlier detection BEFORE normalization.")
}

## -------------------- arrayQualityMetrics BEFORE (optional heavy report) --------------------
aqm_before_dir <- file.path(out_dir, "arrayQualityMetrics_before")
dir.create(aqm_before_dir, showWarnings = FALSE, recursive = TRUE)
try({
  # for AQM we need an ExpressionSet - create temporary from affybatch using rma just for diagnostics if necessary
  eset_for_aqm <- affy::rma(affybatch) # this normalizes but AQM will still give diagnostics; NOTE: this changes values but is only for report
  arrayQualityMetrics::arrayQualityMetrics(expressionset = eset_for_aqm, outdir = aqm_before_dir, force = TRUE, do.logtransform = FALSE)
}, silent = TRUE)

## -------------------- NORMALIZATION (RMA) --------------------
message("Normalizing with RMA ...")
eset_rma <- affy::rma(affybatch)   # ExpressionSet (log2 scale)
exprs_norm <- Biobase::exprs(eset_rma)
message("Normalized expression dimensions: ", nrow(exprs_norm), " probes x ", ncol(exprs_norm), " samples")

# Save full normalized matrix head
write.csv(head(exprs_norm), file = file.path(out_dir, "normalized_exprs_head.csv"), row.names = TRUE)

## -------------------- QC AFTER normalization --------------------
message("QC AFTER normalization ...")
# Boxplot normalized
pdf(file.path(out_dir, "boxplot_normalized.pdf"), width = 10, height = 6)
boxplot(as.data.frame(exprs_norm), main = "RMA normalized expression (log2) - AFTER normalization", las = 2)
dev.off()

# Density normalized
pdf(file.path(out_dir, "density_normalized.pdf"), width = 10, height = 6)
den2 <- apply(exprs_norm, 2, function(x) density(x, na.rm = TRUE))
plot(den2[[1]], main = "Normalized density plots (log2)", xlab = "log2 expression")
for (i in 2:length(den2)) lines(den2[[i]])
dev.off()

# PCA AFTER normalization
pca_after <- prcomp(t(exprs_norm), scale. = TRUE)
pca_df <- data.frame(PC1 = pca_after$x[,1], PC2 = pca_after$x[,2], sample = colnames(exprs_norm), stringsAsFactors = FALSE)

# Merge phenotype if available
if (!is.null(pheno) && "sample" %in% colnames(pheno)) {
  # try match by GSM id or by a substring of sample names
  mapping <- match(pca_df$sample, pheno$sample)
  if (all(is.na(mapping))) {
    # try matching GSM ids that may be in sample names
    mapping <- match(pca_df$sample, pheno$sample)
  }
  if (any(!is.na(mapping))) {
    pca_df <- cbind(pca_df, pheno[match(pca_df$sample, pheno$sample), , drop = FALSE])
  }
}

pdf(file.path(out_dir, "PCA_normalized.pdf"), width = 8, height = 6)
plot(pca_df$PC1, pca_df$PC2, pch = 19, xlab = "PC1", ylab = "PC2", main = "PCA - AFTER normalization")
text(pca_df$PC1, pca_df$PC2, labels = pca_df$sample, cex = 0.6, pos = 4)
dev.off()

## -------------------- Outlier detection AFTER normalization --------------------
pc1a <- pca_after$x[,1]; pc2a <- pca_after$x[,2]
dist_after <- sqrt((pc1a - mean(pc1a))^2 + (pc2a - mean(pc2a))^2)
z_a <- scale(dist_after)
outliers_after <- names(which(z_a > outlier_pca_z_cutoff))

# NUSE/RLE using affyPLM (only for AffyBatch)
outliers_nuse <- character(0)
try({
  plm <- affyPLM::fitPLM(affybatch)
  nuse_stats <- affyPLM::NUSE(plm, type = "stats")
  nuse_median <- nuse_stats[,"Median"]
  outliers_nuse <- names(nuse_median)[which(nuse_median > 1.05)]
}, silent = TRUE)

message("Outliers BEFORE normalization (PCA z>", outlier_pca_z_cutoff, "): ", paste(outliers_before, collapse = ", "))
message("Outliers AFTER normalization  (PCA z>", outlier_pca_z_cutoff, "): ", paste(outliers_after, collapse = ", "))
message("Outliers flagged by NUSE (median>1.05) if any: ", paste(outliers_nuse, collapse = ", "))

## -------------------- arrayQualityMetrics AFTER (report) --------------------
aqm_after_dir <- file.path(out_dir, "arrayQualityMetrics_after")
dir.create(aqm_after_dir, showWarnings = FALSE, recursive = TRUE)
try({
  arrayQualityMetrics::arrayQualityMetrics(expressionset = eset_rma, outdir = aqm_after_dir, force = TRUE, do.logtransform = FALSE)
}, silent = TRUE)

## -------------------- FILTERING low-intensity probes --------------------
message("Filtering low-intensity probes ...")
# RMA output is log2 scale: convert expression_threshold_linear to log2 scale for comparison
log2_thresh <- log2(expression_threshold_linear)
n_samples <- ncol(exprs_norm)
min_required <- ceiling(min_sample_fraction_expressed * n_samples)
mask <- rowSums(exprs_norm > log2_thresh, na.rm = TRUE) >= min_required
probes_before <- nrow(exprs_norm)
probes_after  <- sum(mask)
exprs_filtered <- exprs_norm[mask, ]

message("Probes before filtering: ", probes_before)
message("Probes after filtering:  ", probes_after)

write.csv(data.frame(probes_before = probes_before, probes_after = probes_after),
          file = file.path(out_dir, "filter_summary.csv"), row.names = FALSE)

## -------------------- Relabel groups (auto from GEO metadata) --------------------
message("Relabeling target groups ...")
groups_final <- rep("Unknown", n_samples)
names(groups_final) <- colnames(exprs_norm)
if (!is.null(pheno)) {
  # try to detect a characteristic column
  textcols <- sapply(pheno, function(col) any(grepl("cancer|tumor|tumour|diabetes|control|healthy|normal", tolower(as.character(col)), perl = TRUE)))
  if (any(textcols)) {
    colname <- names(which(textcols))[1]
    message("Using phenotype column for labels: ", colname)
    # pheno must have sample column; if not, try to guess sample order
    if ("sample" %in% colnames(pheno)) {
      ordering <- match(colnames(exprs_norm), pheno$sample)
      labs <- tolower(as.character(pheno[[colname]][ordering]))
    } else {
      labs <- tolower(as.character(pheno[[colname]]))
      labs <- labs[1:length(groups_final)]
    }
    # map to canonical Normal / Disease labels
    labs[grepl("normal|control|healthy", labs)] <- "Normal"
    labs[grepl("cancer|tumor|tumour|pancrea", labs)] <- "Disease"
    labs[grepl("diabetes", labs)] <- "Diabetes"
    # fallback keep original if not matched
    groups_final <- ifelse(is.na(labs), "Unknown", labs)
    names(groups_final) <- colnames(exprs_norm)
  } else {
    message("No phenotype text column auto-detected; setting groups to Unknown")
  }
} else {
  message("No phenotype metadata available; groups remain 'Unknown'.")
}

# Save sample->group mapping
write.csv(data.frame(sample = names(groups_final), group = groups_final), file = file.path(out_dir, "sample_groups.csv"), row.names = FALSE)

## -------------------- Save outputs --------------------
if (save_rds) {
  saveRDS(list(affybatch = affybatch, eset_rma = eset_rma, exprs_norm = exprs_norm,
               exprs_filtered = exprs_filtered, groups = groups_final,
               outliers_before = outliers_before, outliers_after = outliers_after,
               outliers_nuse = outliers_nuse),
          file = file.path(out_dir, "preprocessing_results.rds"))
}

# Write short summary file for the assignment form
summary_file <- file.path(out_dir, "preprocessing_summary.txt")
cat("GSE15932 Preprocessing summary\n", file = summary_file)
cat("Dataset accession: GSE15932\n", file = summary_file, append = TRUE)
cat("Samples (normalized): ", ncol(exprs_norm), "\n", file = summary_file, append = TRUE)
cat("Probes before filtering: ", probes_before, "\n", file = summary_file, append = TRUE)
cat("Probes after filtering:  ", probes_after, "\n", file = summary_file, append = TRUE)
cat("Outliers BEFORE normalization (PCA z>", outlier_pca_z_cutoff, "): ", paste(outliers_before, collapse = ", "), "\n", file = summary_file, append = TRUE)
cat("Outliers AFTER normalization  (PCA z>", outlier_pca_z_cutoff, "): ", paste(outliers_after, collapse = ", "), "\n", file = summary_file, append = TRUE)
cat("Outliers flagged by NUSE (if any): ", paste(outliers_nuse, collapse = ", "), "\n", file = summary_file, append = TRUE)
cat("Group counts:\n", file = summary_file, append = TRUE)
capture.output(table(groups_final), file = summary_file, append = TRUE)

message("All outputs written to: ", out_dir)
message("Key files: boxplot_normalized.pdf, PCA_normalized.pdf, filter_summary.csv, preprocessing_summary.txt, sample_groups.csv, preprocessing_results.rds")

###############################################################################
# End of script
###############################################################################
