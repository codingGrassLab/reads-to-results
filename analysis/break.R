# Reads to Results: what goes wrong, on the real airway data.
# Every number on the "Break it yourself" and "Mistakes" pages comes from this script.
#   Rscript analysis/break.R        (from the project root; needs the packages of run.R)
# Output: src/data/break.json

suppressPackageStartupMessages({
  library(airway); library(DESeq2); library(apeglm); library(jsonlite)
  library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)
})
set.seed(1)
data(airway)
se <- airway
se$dex <- relevel(se$dex, "untrt")
cts <- assay(se)
cd0 <- data.frame(row.names = colnames(se), cell = factor(se$cell), dex = relevel(factor(se$dex), "untrt"))
sym <- function(ids) { k <- keys(org.Hs.eg.db, "ENSEMBL"); out <- rep(NA_character_, length(ids)); ok <- ids %in% k; if (any(ok)) out[ok] <- unname(mapIds(org.Hs.eg.db, ids[ok], "SYMBOL", "ENSEMBL", multiVals = "first")); out }
markers <- c(FKBP5 = "ENSG00000096060", TSC22D3 = "ENSG00000157514", PER1 = "ENSG00000179094", KLF15 = "ENSG00000163884",
             DUSP1 = "ENSG00000120129", ZBTB16 = "ENSG00000109906", CRISPLD2 = "ENSG00000103196")

# one DESeq2 run on a subset of samples with a given design; the filter follows the vignette (>= 10 counts in the smallest group)
run <- function(samples, design, cd = cd0, counts = cts) {
  cd <- droplevels(cd[samples, , drop = FALSE]); m <- counts[, samples]
  minGroup <- min(table(cd$dex))
  dds <- DESeqDataSetFromMatrix(m, cd, design)
  dds <- dds[rowSums(counts(dds) >= 10) >= minGroup, ]
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("dex", "trt", "untrt"), alpha = 0.05)
  list(dds = dds, res = res, minGroup = minGroup)
}
summ <- function(r, base = NULL) {
  res <- r$res; sig <- which(!is.na(res$padj) & res$padj < 0.05); ids <- rownames(res)[sig]
  out <- list(tested = sum(!is.na(res$padj)), up = sum(res$log2FoldChange[sig] > 0), down = sum(res$log2FoldChange[sig] < 0),
              sig = length(sig), raw_p = sum(!is.na(res$pvalue) & res$pvalue < 0.05), filter_min = r$minGroup,
              markers = lapply(names(markers), function(g) { i <- markers[[g]]; list(symbol = g, padj = if (i %in% rownames(res)) signif(res[i, "padj"], 3) else NA, lfc = if (i %in% rownames(res)) round(res[i, "log2FoldChange"], 2) else NA) }),
              top = { o <- order(res$padj)[1:10]; lapply(o, function(k) list(id = rownames(res)[k], symbol = sym(rownames(res)[k]), lfc = round(res$log2FoldChange[k], 2), padj = signif(res$padj[k], 3))) })
  if (!is.null(base)) { out$kept <- length(intersect(ids, base)); out$lost <- length(setdiff(base, ids)); out$new <- length(setdiff(ids, base)) }
  out$ids <- ids
  out
}
tryMsg <- function(expr) tryCatch({ force(expr); "no error" }, error = function(e) conditionMessage(e))

all8 <- colnames(se)
cellOf <- setNames(as.character(cd0$cell), all8)
# ── 1. The switches: design × samples (8 real runs) ─────────────────────────
sampleSets <- list(
  all = all8,
  drop1 = setdiff(all8, "SRR1039521"),                                   # one treated sample lost (N061011): 3 vs 4
  n2 = all8[cellOf %in% c("N61311", "N052611")],                         # only two donors: 2 vs 2
  swap = all8                                                            # labels of one donor's pair swapped
)
cdSwap <- cd0; cdSwap[c("SRR1039508", "SRR1039509"), "dex"] <- cd0[c("SRR1039509", "SRR1039508"), "dex"]
designs <- list(paired = ~ cell + dex, unpaired = ~ dex)
base <- summ(run(all8, ~ cell + dex))
base_ids <- base$ids
grid <- list()
for (s in names(sampleSets)) for (d in names(designs)) {
  r <- run(sampleSets[[s]], designs[[d]], cd = if (s == "swap") cdSwap else cd0)
  x <- summ(r, base_ids); x$ids <- NULL
  grid[[paste(s, d, sep = "_")]] <- c(list(samples = s, design = d, formula = deparse(designs[[d]]), n = length(sampleSets[[s]])), x)
  cat(s, d, ":", x$sig, "significant (", x$up, "up,", x$down, "down ), kept", x$kept, "lost", x$lost, "new", x$new, "\n")
}

# ── 2. Mistakes ──────────────────────────────────────────────────────────────
mist <- list()
# confounded batch: a 'batch' column identical to the treatment
cdConf <- cd0; cdConf$batch <- factor(ifelse(cd0$dex == "trt", "B2", "B1"))
mist$confounded <- list(error = tryMsg(DESeqDataSetFromMatrix(cts, cdConf, ~ batch + dex)))
# one sample per group
mist$one_vs_one <- list(error = tryMsg({ d1 <- DESeqDataSetFromMatrix(cts[, c("SRR1039508", "SRR1039509")], droplevels(cd0[c("SRR1039508", "SRR1039509"), ]), ~ dex); DESeq(d1, quiet = TRUE) }))
# normalized values instead of raw counts
dds0 <- estimateSizeFactors(DESeqDataSetFromMatrix(cts, cd0, ~ cell + dex))
norm <- counts(dds0, normalized = TRUE)
mist$normalized_error <- list(error = tryMsg(DESeqDataSetFromMatrix(norm, cd0, ~ cell + dex)))
rn <- summ(run(all8, ~ cell + dex, counts = round(norm)), base_ids); rn$ids <- NULL; rn$top <- NULL
mist$normalized_rounded <- rn
tpm_like <- round(sweep(cts, 2, colSums(cts), "/") * 1e6)          # per-million scaling, rounded: what people do with TPM/CPM
rt <- summ(run(all8, ~ cell + dex, counts = tpm_like), base_ids); rt$ids <- NULL; rt$top <- NULL
mist$cpm_rounded <- rt
# fold-change cut after the test vs testing against a fold-change threshold
rb <- run(all8, ~ cell + dex)
r0 <- rb$res; rT <- results(rb$dds, contrast = c("dex", "trt", "untrt"), lfcThreshold = 1, altHypothesis = "greaterAbs", alpha = 0.05)
mist$lfc_threshold <- list(posthoc = sum(r0$padj < 0.05 & abs(r0$log2FoldChange) > 1, na.rm = TRUE), tested = sum(rT$padj < 0.05, na.rm = TRUE))
sh <- lfcShrink(rb$dds, coef = "dex_trt_vs_untrt", type = "apeglm", quiet = TRUE)
# big raw fold changes without evidence: what apeglm shrinkage does to them
noev <- which(!is.na(r0$pvalue) & (is.na(r0$padj) | r0$padj >= 0.05))
big <- noev[abs(r0$log2FoldChange[noev]) > 2]
o <- head(big[order(-abs(r0$log2FoldChange[big]))], 10)
mist$shrinkage <- list(n_big_not_sig = length(big), n_still_big_after_shrink = sum(abs(sh$log2FoldChange[big]) > 2),
                       median_baseMean_big = round(median(r0$baseMean[big]), 1), median_baseMean_all = round(median(r0$baseMean[!is.na(r0$pvalue)]), 1),
                       rows = lapply(o, function(k) list(id = rownames(r0)[k], symbol = sym(rownames(r0)[k]), baseMean = round(r0$baseMean[k], 1), lfc = round(r0$log2FoldChange[k], 2), shrunk = round(sh$log2FoldChange[k], 2), padj = signif(r0$padj[k], 3))))
# enrichment background: tested genes vs the whole genome
sig1 <- rownames(r0)[which(r0$padj < 0.05 & abs(r0$log2FoldChange) > 1)]
universe <- rownames(r0)[!is.na(r0$padj)]
eU <- as.data.frame(enrichGO(sig1, universe = universe, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05))
eG <- as.data.frame(enrichGO(sig1, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05))
onlyG <- setdiff(eG$ID, eU$ID)
mist$background <- list(genes = length(sig1), universe = length(universe), terms_tested_bg = nrow(eU), terms_genome_bg = nrow(eG), shared = length(intersect(eU$ID, eG$ID)), only_genome = length(onlyG),
                        examples = head(eG[eG$ID %in% onlyG, c("Description", "p.adjust")][order(eG$p.adjust[eG$ID %in% onlyG]), ], 5))
# Excel turning gene names into dates, in a real GEO file (GSE311711, 2025)
f <- tempfile(fileext = ".tsv.gz")
ok <- tryCatch({ download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311711/suppl/GSE311711_raw_counts_All_Samples.tsv.gz", f, quiet = TRUE, mode = "wb"); TRUE }, error = function(e) FALSE)
if (ok) {
  g <- read.delim(gzfile(f), check.names = FALSE)
  dates <- g$GeneSymbol[grepl("^[0-9]{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$", g$GeneSymbol)]
  mist$excel <- list(file = "GSE311711_raw_counts_All_Samples.tsv.gz", rows = nrow(g), n = length(dates), examples = head(sort(unique(dates)), 12),
                     ensembl = head(g$ENSG[g$GeneSymbol %in% dates], 12))
}
# p instead of padj on the real comparison
mist$p_vs_padj <- list(raw_p = base$raw_p, padj = base$sig)

write_json(list(base = { b <- base; b$ids <- NULL; b }, grid = grid, mistakes = mist,
                versions = list(deseq2 = as.character(packageVersion("DESeq2")), apeglm = as.character(packageVersion("apeglm")), clusterProfiler = as.character(packageVersion("clusterProfiler")))),
           "src/data/break.json", auto_unbox = TRUE, digits = NA, na = "null", pretty = FALSE)
cat("written src/data/break.json\n")
