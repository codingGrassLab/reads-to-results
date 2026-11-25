# Reads to Results: the whole airway analysis, run once.
# Every number and plot on the site comes from this script. Run from the project root:
#   Rscript analysis/run.R
# Outputs: src/data/*.json (read by the site at build time), public/data/* (downloads), public/plots/*.png

suppressPackageStartupMessages({
  library(airway)
  library(DESeq2)
  library(jsonlite)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(apeglm)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GenomicAlignments)
  library(Rsamtools)
  library(rtracklayer)
  library(EnhancedVolcano)
})

set.seed(1)
out_data <- "src/data"; out_pub <- "public/data"; out_plots <- "public/plots"
for (d in c(out_data, out_pub, out_plots)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
wj <- function(x, name, digits = NA) write_json(x, file.path(out_data, name), digits = digits, auto_unbox = TRUE, pretty = FALSE, na = "null")
png_ <- function(name, w = 1600, h = 1200) png(file.path(out_plots, name), width = w, height = h, res = 200)

# ── 1. The data ──────────────────────────────────────────────────────────────
data(airway)
se <- airway
se$dex <- relevel(se$dex, "untrt")
coldata <- as.data.frame(colData(se))
coldata <- coldata[, c("SampleName", "cell", "dex", "Run", "avgLength", "Experiment", "Sample", "BioSample")]
counts <- assay(se)

samples <- data.frame(
  run = colnames(se), cell = as.character(se$cell), dex = as.character(se$dex),
  libsize = unname(colSums(counts)), detected = unname(colSums(counts > 0)),
  avgLength = unname(se$avgLength), stringsAsFactors = FALSE
)

# ── 2. What a read looks like: three real reads from a subset BAM ─────────────
extdir <- system.file("extdata", package = "airway")
bam1 <- file.path(extdir, "SRR1039508_subset.bam")
reads <- scanBam(BamFile(bam1, yieldSize = 3), param = ScanBamParam(what = c("qname", "flag", "rname", "pos", "mapq", "cigar", "seq", "qual")))[[1]]
read_examples <- lapply(1:3, function(i) list(
  name = reads$qname[i], flag = reads$flag[i], chrom = as.character(reads$rname[i]), pos = reads$pos[i],
  mapq = reads$mapq[i], cigar = reads$cigar[i], seq = as.character(reads$seq[i]), qual = as.character(reads$qual[i])
))

# ── 3. Counting: reads per gene on the subset BAMs (the real counting step, small) ─
gtf <- import(file.path(extdir, "Homo_sapiens.GRCh37.75_subset.gtf"))
exons <- gtf[gtf$type == "exon"]
ebg <- split(exons, exons$gene_id)
bamfiles <- BamFileList(file.path(extdir, paste0(colnames(se), "_subset.bam")), yieldSize = 2e6)
se_sub <- summarizeOverlaps(features = ebg, reads = bamfiles, mode = "Union", singleEnd = FALSE, ignore.strand = TRUE, fragments = TRUE)
sub_counts <- assay(se_sub); colnames(sub_counts) <- colnames(se)
sub_counts <- sub_counts[order(-rowSums(sub_counts)), , drop = FALSE]
sym_sub <- mapIds(org.Hs.eg.db, rownames(sub_counts), "SYMBOL", "ENSEMBL", multiVals = "first")
counting <- list(
  region = list(chrom = as.character(seqnames(gtf))[1], start = min(start(gtf)), end = max(end(gtf))),
  genes = rownames(sub_counts), symbols = unname(ifelse(is.na(sym_sub), "", sym_sub)),
  counts = unname(as.list(as.data.frame(t(sub_counts)))), samples = colnames(se)
)
# full-matrix excerpt for the lesson
head_ids <- c("ENSG00000000003", "ENSG00000000005", "ENSG00000000419", "ENSG00000000457", "ENSG00000000460", "ENSG00000000938")
head_ids <- head_ids[head_ids %in% rownames(counts)]
counts_head <- list(genes = head_ids, counts = unname(as.list(as.data.frame(t(counts[head_ids, ])))), samples = colnames(counts),
                    n_genes = nrow(counts), n_zero = sum(rowSums(counts) == 0))

# ── 4. DESeq2 ────────────────────────────────────────────────────────────────
dds <- DESeqDataSet(se, design = ~ cell + dex)
n_before <- nrow(dds)
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep, ]
n_after <- nrow(dds)
dds <- DESeq(dds)
res <- results(dds, contrast = c("dex", "trt", "untrt"), alpha = 0.05)
resLFC <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm")

sym <- mapIds(org.Hs.eg.db, rownames(res), "SYMBOL", "ENSEMBL", multiVals = "first")
tab <- data.frame(
  gene_id = rownames(res), symbol = unname(ifelse(is.na(sym), "", sym)),
  baseMean = res$baseMean, log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE, stat = res$stat,
  pvalue = res$pvalue, padj = res$padj, shrunk_log2FoldChange = resLFC$log2FoldChange, stringsAsFactors = FALSE
)
tab <- tab[order(tab$padj, tab$pvalue), ]
write.csv(tab, file.path(out_pub, "airway_deseq2_results.csv"), row.names = FALSE)

sig <- !is.na(tab$padj) & tab$padj < 0.05
summary_stats <- list(
  genes_in_matrix = n_before, genes_after_filter = n_after, genes_tested = sum(!is.na(tab$padj)),
  independent_filtered = sum(is.na(tab$padj)),
  sig_005 = sum(sig), up_005 = sum(sig & tab$log2FoldChange > 0), down_005 = sum(sig & tab$log2FoldChange < 0),
  sig_005_lfc1 = sum(sig & abs(tab$log2FoldChange) > 1),
  up_005_lfc1 = sum(sig & tab$log2FoldChange > 1), down_005_lfc1 = sum(sig & tab$log2FoldChange < -1),
  sig_01 = sum(!is.na(tab$padj) & tab$padj < 0.1),
  raw_p_005 = sum(!is.na(tab$pvalue) & tab$pvalue < 0.05),
  size_factors = setNames(as.list(round(sizeFactors(dds), 4)), colnames(dds)),
  design = "~ cell + dex", contrast = "dex: trt vs untrt", alpha = 0.05, filter = "rowSums(counts >= 10) >= 4",
  deseq2 = as.character(packageVersion("DESeq2")), apeglm = as.character(packageVersion("apeglm")),
  airway = as.character(packageVersion("airway")), r = R.version.string
)

# genes the original paper (Himes et al. 2014) discussed, for the validation page
paper_genes <- c("CRISPLD2", "DUSP1", "KLF15", "PER1", "TSC22D3", "FKBP5", "ZBTB16", "SPARCL1", "MAOA", "GLUL", "KCTD12", "C13orf15", "STEAP4", "SAMHD1", "ANGPTL4")
paper_check <- tab[tab$symbol %in% paper_genes, c("gene_id", "symbol", "baseMean", "log2FoldChange", "padj")]
paper_check <- paper_check[order(paper_check$padj), ]; rownames(paper_check) <- NULL

# compact table for the interactive volcano / MA plots: tested genes only
tested <- tab[!is.na(tab$padj), ]
volcano <- list(
  cols = c("id", "symbol", "baseMean", "lfc", "slfc", "pvalue", "padj"),
  rows = unname(Map(function(a, b, c, d, e, f, g) list(a, b, c, d, e, f, g),
    tested$gene_id, tested$symbol, round(tested$baseMean, 1), round(tested$log2FoldChange, 3),
    round(tested$shrunk_log2FoldChange, 3), signif(tested$pvalue, 5), signif(tested$padj, 5)))
)
top <- head(tab[sig, ], 15)
top_genes <- lapply(seq_len(nrow(top)), function(i) list(
  gene_id = top$gene_id[i], symbol = top$symbol[i], baseMean = round(top$baseMean[i], 1),
  log2FoldChange = round(top$log2FoldChange[i], 2), shrunk = round(top$shrunk_log2FoldChange[i], 2),
  pvalue = signif(top$pvalue[i], 2), padj = signif(top$padj[i], 2)
))

# p-value histogram (20 bins)
ph <- hist(tab$pvalue[!is.na(tab$pvalue)], breaks = seq(0, 1, by = 0.05), plot = FALSE)
pvalue_hist <- list(breaks = ph$breaks, counts = ph$counts)

# ── 5. Transformation, PCA, distances, heatmap ──────────────────────────────
vsd <- vst(dds, blind = FALSE)
pca <- plotPCA(vsd, intgroup = c("dex", "cell"), returnData = TRUE)
pv <- round(100 * attr(pca, "percentVar"))
pca_json <- list(
  percentVar = pv,
  points = lapply(seq_len(nrow(pca)), function(i) list(run = rownames(pca)[i], PC1 = round(pca$PC1[i], 3), PC2 = round(pca$PC2[i], 3), dex = as.character(pca$dex[i]), cell = as.character(pca$cell[i])))
)
sampleDists <- dist(t(assay(vsd)))
dist_mat <- as.matrix(sampleDists)
distances <- list(samples = colnames(dist_mat), matrix = unname(as.list(as.data.frame(round(dist_mat, 2)))))

topVar <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 25)
mat <- assay(vsd)[topVar, ]
mat_z <- t(scale(t(mat)))
sym_top <- mapIds(org.Hs.eg.db, rownames(mat), "SYMBOL", "ENSEMBL", multiVals = "first")
heatmap_json <- list(
  genes = rownames(mat), symbols = unname(ifelse(is.na(sym_top), rownames(mat), sym_top)), samples = colnames(mat),
  dex = as.character(colData(vsd)$dex), cell = as.character(colData(vsd)$cell),
  z = unname(lapply(seq_len(nrow(mat_z)), function(i) round(unname(mat_z[i, ]), 3))),
  vst = unname(lapply(seq_len(nrow(mat)), function(i) round(unname(mat[i, ]), 2)))
)

# normalization: per-sample quantiles of log2 counts before and after size-factor scaling
qs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
raw_log <- log2(counts(dds) + 1); norm_log <- log2(counts(dds, normalized = TRUE) + 1)
normalization <- list(
  samples = colnames(dds), quantiles = qs,
  raw = unname(lapply(seq_len(ncol(dds)), function(j) round(unname(quantile(raw_log[, j], qs)), 3))),
  normalized = unname(lapply(seq_len(ncol(dds)), function(j) round(unname(quantile(norm_log[, j], qs)), 3))),
  libsize = unname(colSums(counts(dds))), sizeFactor = unname(round(sizeFactors(dds), 4))
)

# dispersion: a 3000-gene sample of the fitted values for the interactive plot
disp <- data.frame(baseMean = mcols(dds)$baseMean, gene = mcols(dds)$dispGeneEst, fit = mcols(dds)$dispFit, final = dispersions(dds))
disp <- disp[is.finite(disp$gene) & disp$gene > 1e-8, ]
disp_s <- disp[sample(nrow(disp), min(3000, nrow(disp))), ]
dispersion <- list(cols = c("baseMean", "gene", "fit", "final"),
                   rows = unname(Map(function(a, b, c, d) c(a, b, c, d), signif(disp_s$baseMean, 4), signif(disp_s$gene, 4), signif(disp_s$fit, 4), signif(disp_s$final, 4))))

# one gene, all samples: for the "what does log2 fold change mean" page
gene_of_interest <- tab$gene_id[tab$symbol == "CRISPLD2"][1]
cnt <- plotCounts(dds, gene = gene_of_interest, intgroup = c("dex", "cell"), returnData = TRUE)
gene_counts <- list(
  gene_id = gene_of_interest, symbol = "CRISPLD2",
  points = lapply(seq_len(nrow(cnt)), function(i) list(run = rownames(cnt)[i], count = round(cnt$count[i], 1), dex = as.character(cnt$dex[i]), cell = as.character(cnt$cell[i]))),
  result = as.list(tab[tab$gene_id == gene_of_interest, c("baseMean", "log2FoldChange", "shrunk_log2FoldChange", "lfcSE", "pvalue", "padj")])
)

# ── 6. GO enrichment (biological process) ───────────────────────────────────
sig_ids <- tab$gene_id[sig & abs(tab$log2FoldChange) > 1]
universe <- tab$gene_id[!is.na(tab$padj)]
ego <- enrichGO(gene = sig_ids, universe = universe, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
ego_df <- as.data.frame(ego)
go_top <- head(ego_df, 20)
go_json <- list(
  n_genes = length(sig_ids), n_universe = length(universe), n_terms = nrow(ego_df),
  terms = lapply(seq_len(nrow(go_top)), function(i) list(
    id = go_top$ID[i], description = go_top$Description[i], geneRatio = go_top$GeneRatio[i], bgRatio = go_top$BgRatio[i],
    padj = signif(go_top$p.adjust[i], 3), count = go_top$Count[i], genes = strsplit(go_top$geneID[i], "/")[[1]]
  ))
)

# ── 7. Static plots for the lessons ─────────────────────────────────────────
png_("pca.png"); print(ggplot(pca, aes(PC1, PC2, color = dex, shape = cell)) + geom_point(size = 4) +
  xlab(paste0("PC1: ", pv[1], "% variance")) + ylab(paste0("PC2: ", pv[2], "% variance")) + coord_fixed() + theme_bw(base_size = 14)); dev.off()

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
ann <- data.frame(dex = colData(vsd)$dex, cell = colData(vsd)$cell, row.names = colnames(vsd))
png_("sample-distances.png"); pheatmap(dist_mat, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors, annotation_col = ann); dev.off()
png_("heatmap.png", 1600, 1800); pheatmap(mat_z, annotation_col = ann, labels_row = heatmap_json$symbols, fontsize_row = 9); dev.off()
png_("ma.png"); plotMA(res, ylim = c(-5, 5), alpha = 0.05, main = "MA plot (unshrunken)"); dev.off()
png_("ma-shrunk.png"); plotMA(resLFC, ylim = c(-5, 5), alpha = 0.05, main = "MA plot (apeglm shrinkage)"); dev.off()
png_("volcano.png", 1800, 1600); print(EnhancedVolcano(tab, lab = tab$symbol, x = "log2FoldChange", y = "padj", pCutoff = 0.05, FCcutoff = 1,
  title = "Dexamethasone vs untreated", subtitle = "airway, DESeq2", ylab = bquote(~-Log[10]~ "adjusted p"), labSize = 3.5, pointSize = 1.5, legendPosition = "bottom")); dev.off()
png_("gene-counts.png", 1400, 1100); print(ggplot(cnt, aes(dex, count, color = cell, group = cell)) + geom_point(size = 3.5) + geom_line() +
  scale_y_log10() + ggtitle("CRISPLD2, normalized counts") + theme_bw(base_size = 14)); dev.off()
if (nrow(ego_df) > 0) { png_("go-dotplot.png", 1800, 1600); print(dotplot(ego, showCategory = 15) + theme(axis.text.y = element_text(size = 10))); dev.off() }

# ── 8. Write everything ─────────────────────────────────────────────────────
wj(samples, "samples.json"); wj(read_examples, "reads.json")
wj(counting, "counting.json"); wj(counts_head, "counts_head.json")
wj(summary_stats, "summary.json"); wj(top_genes, "top_genes.json")
wj(pvalue_hist, "pvalue_hist.json"); wj(pca_json, "pca.json"); wj(distances, "distances.json"); wj(heatmap_json, "heatmap.json")
wj(normalization, "normalization.json"); wj(gene_counts, "gene_counts.json")
# the two big tables are fetched by the browser, so they live with the downloads
write_json(volcano, file.path(out_pub, "volcano.json"), digits = NA, auto_unbox = TRUE, na = "null")
write_json(dispersion, file.path(out_pub, "dispersion.json"), digits = NA, auto_unbox = TRUE, na = "null")
wj(go_json, "go.json"); wj(paper_check, "paper_check.json")
pk <- c("DESeq2", "apeglm", "airway", "clusterProfiler", "org.Hs.eg.db", "GenomicAlignments", "pheatmap", "ggplot2", "EnhancedVolcano")
wj(list(r = R.version.string, bioc = as.character(BiocManager::version()), packages = setNames(as.list(sapply(pk, function(p) as.character(packageVersion(p)))), pk),
        date = format(Sys.Date())), "versions.json")
cat("done:", summary_stats$sig_005, "genes at padj<0.05;", summary_stats$sig_005_lfc1, "with |LFC|>1;", nrow(ego_df), "GO terms\n")
