
setwd("C:/Users/tahaj/Desktop/Assignment 2")

input_dir  <- "Raw Data"
output_dir <- "Results"

if (!dir.exists(output_dir)) dir.create(output_dir)

deg_files <- list.files(input_dir, pattern = "DEG", ignore.case = TRUE)
deg_files <- deg_files[grepl("\\.csv$", deg_files, ignore.case = TRUE)]

if (length(deg_files) == 0) {
  candidates <- c("DEGs_Data_1.csv","DEGs_Data_2.csv","DEGs_data_1.csv","DEGs_data_2.csv")
  deg_files <- candidates[file.exists(file.path(input_dir, candidates))]
}

if (length(deg_files) == 0) {
  stop("No DEGs CSV files found in 'Raw Data'. Please place DEGs_Data_1.csv and DEGs_Data_2.csv inside the 'Raw Data' folder.")
}

find_col <- function(df, patterns) {
  cols <- colnames(df)
  for (p in patterns) {
    idx <- grep(p, cols, ignore.case = TRUE)
    if (length(idx) > 0) return(cols[idx[1]])
  }
  return(NA)
}

classify_gene <- function(logFC, padj) {
  # padj assumed numeric and NA-replaced before calling
  if (is.na(logFC)) {
    return("Not_Significant")
  }
  if (logFC > 1 & padj < 0.05) {
    return("Upregulated")
  } else if (logFC < -1 & padj < 0.05) {
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}


for (file_name in deg_files) {
  cat("\n-------------------------------\n")
  cat("Processing file:", file_name, "\n")
  infile <- file.path(input_dir, file_name)
  df <- read.csv(infile, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  
  
  padj_col  <- find_col(df, c("^padj$", "padj", "adj.p.val", "p.adjust", "adjusted.p.value", "padj"))
  logfc_col <- find_col(df, c("^logfc$", "log2foldchange", "log2fc", "logFC", "log2FoldChange", "log.fold.change"))
  gene_col  <- find_col(df, c("^gene_id$", "gene_id", "^geneid$", "gene", "Gene_Id", "Gene"))
  
  cat("Detected columns -> gene:", ifelse(is.na(gene_col), "<none>", gene_col),
      "; logFC:", ifelse(is.na(logfc_col), "<none>", logfc_col),
      "; padj:", ifelse(is.na(padj_col), "<none>", padj_col), "\n")
  
  
  if (is.na(padj_col)) {
    warning("padj column not found. Creating 'padj' column and filling with 1.")
    df$padj <- 1
    padj_col <- "padj"
  }
  
  
  df[[padj_col]][is.na(df[[padj_col]])] <- 1
  
  df[[padj_col]] <- as.numeric(df[[padj_col]])
  if (is.na(logfc_col)) {
    stop("ERROR: logFC column not found in file '", file_name, "'. Expected 'logFC' or 'log2FoldChange'.")
  }
  df[[logfc_col]] <- as.numeric(df[[logfc_col]])
  
  df$status <- NA_character_
  
  for (i in seq_len(nrow(df))) {
    lf <- df[[logfc_col]][i]
    pv <- df[[padj_col]][i]
    # classify_gene handles NA logFC by returning 'Not_Significant'
    df$status[i] <- classify_gene(lf, pv)
  }
  
  out_name <- paste0(tools::file_path_sans_ext(basename(file_name)), "_classified.csv")
  out_path <- file.path(output_dir, out_name)
  write.csv(df, out_path, row.names = FALSE)
  cat("Saved classified file to:", out_path, "\n")
  
  cat("Summary counts (status):\n")
  print(table(df$status))
}

cat("\nAll files processed. Classified files are in the 'Results' folder.\n")
save.image(file = "C:/Users/tahaj/Desktop/Assignment 2")

