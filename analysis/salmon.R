# Salmon outputs that ship with the airway package (two samples, quantified by the
# package authors against GENCODE v29 with Salmon 0.14.1). Run from the project root:
#   Rscript analysis/salmon.R
# Writes src/data/salmon.json and, if analysis/tx2gene.csv exists, src/data/tximport.json.

suppressPackageStartupMessages({ library(jsonlite); library(tximport); library(readr) })
qdir <- system.file("extdata", "quants", package = "airway")
runs <- list.dirs(qdir, recursive = FALSE, full.names = FALSE)
dir.create("src/data", showWarnings = FALSE, recursive = TRUE)

per_run <- lapply(runs, function(r) {
  d <- file.path(qdir, r)
  meta <- fromJSON(file.path(d, "aux_info", "meta_info.json"))
  lib <- fromJSON(file.path(d, "lib_format_counts.json"))
  cmd <- fromJSON(file.path(d, "cmd_info.json"))
  q <- read_tsv(file.path(d, "quant.sf.gz"), show_col_types = FALSE)
  fld <- as.numeric(strsplit(readLines(file.path(d, "libParams", "flenDist.txt"))[1], "\t")[[1]])
  fld_mean <- sum(seq_along(fld) * fld) / sum(fld)
  top <- q[order(-q$NumReads), ][1:8, ]
  list(
    run = r, salmon_version = meta$salmon_version, index = cmd$index, libType = cmd$libType,
    num_processed = meta$num_processed, num_mapped = meta$num_mapped, percent_mapped = round(meta$percent_mapped, 2),
    num_targets = meta$num_valid_targets, num_eq_classes = meta$num_eq_classes, gc_bias = meta$gc_bias_correct,
    expected_format = lib$expected_format, ISF = lib$ISF, ISR = lib$ISR, strand_bias = round(lib$strand_mapping_bias, 4),
    num_dovetail = meta$num_dovetail_fragments, num_filtered_vm = meta$num_fragments_filtered_vm,
    start = meta$start_time, end = meta$end_time,
    transcripts = nrow(q), transcripts_with_reads = sum(q$NumReads > 0), sum_numreads = round(sum(q$NumReads)), sum_tpm = round(sum(q$TPM)),
    fld_mean = round(fld_mean, 1),
    head = lapply(1:4, function(i) as.list(q[i, ])),
    top = lapply(seq_len(nrow(top)), function(i) list(Name = top$Name[i], Length = top$Length[i], EffectiveLength = round(top$EffectiveLength[i], 1), TPM = round(top$TPM[i], 1), NumReads = round(top$NumReads[i], 1)))
  )
})
names(per_run) <- runs
write_json(per_run, "src/data/salmon.json", auto_unbox = TRUE, digits = NA, pretty = FALSE)
cat("salmon.json:", paste(runs, collapse = ", "), "\n")

# gene-level import, when a transcript-to-gene table is available
t2g_file <- "analysis/tx2gene.csv"
if (file.exists(t2g_file)) {
  tx2gene <- read.csv(t2g_file)
  tx2gene$TXNAME <- sub("[.].*", "", tx2gene$TXNAME)   # drop transcript versions to match quant.sf after ignoreTxVersion
  tx2gene$GENEID <- sub("[.].*", "", tx2gene$GENEID)   # gene versions too, so IDs match org.Hs.eg.db and the airway matrix
  files <- setNames(file.path(qdir, runs, "quant.sf.gz"), runs)
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
  cts <- txi$counts
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  pick <- c("ENSG00000103196", "ENSG00000171819", "ENSG00000152583", "ENSG00000120129", "ENSG00000000003", "ENSG00000000005")
  pick <- pick[pick %in% rownames(cts)]
  sym <- AnnotationDbi::mapIds(org.Hs.eg.db, pick, "SYMBOL", "ENSEMBL", multiVals = "first")
  out <- list(
    runs = runs, n_transcripts_in = nrow(tx2gene),
    n_genes = nrow(cts), genes_with_counts = colSums(cts > 0), total_counts = round(colSums(cts)),
    pick = lapply(pick, function(g) list(gene_id = g, symbol = unname(ifelse(is.na(sym[g]), "", sym[g])),
      counts = round(unname(cts[g, ]), 1), abundance = round(unname(txi$abundance[g, ]), 2), length = round(unname(txi$length[g, ]), 0))),
    countsFromAbundance = txi$countsFromAbundance
  )
  write_json(out, "src/data/tximport.json", auto_unbox = TRUE, digits = NA)
  cat("tximport.json:", nrow(cts), "genes\n")
} else cat("no", t2g_file, "- skipping tximport\n")
