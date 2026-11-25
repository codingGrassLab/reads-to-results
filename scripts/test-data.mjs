// Checks that the data files the site is built from are consistent with each other.
// Run with: npm test
import { readFileSync } from "node:fs";
const j = (p) => JSON.parse(readFileSync(new URL(`../${p}`, import.meta.url), "utf8"));
let failed = 0;
const ok = (cond, msg) => { if (!cond) { failed++; console.error("FAIL", msg); } else console.log("ok  ", msg); };
const close = (a, b, tol) => Math.abs(a - b) <= tol;

const summary = j("src/data/summary.json"), samples = j("src/data/samples.json"), pca = j("src/data/pca.json");
const top = j("src/data/top_genes.json"), paper = j("src/data/paper_check.json"), hm = j("src/data/heatmap.json");
const norm = j("src/data/normalization.json"), gene = j("src/data/gene_counts.json"), hist = j("src/data/pvalue_hist.json");
const go = j("src/data/go.json"), counting = j("src/data/counting.json"), head = j("src/data/counts_head.json");
const salmon = j("src/data/salmon.json"), txi = j("src/data/tximport.json");
const volcano = j("public/data/volcano.json"), disp = j("public/data/dispersion.json"), dist = j("src/data/distances.json");

ok(samples.length === 8, "8 samples");
ok(samples.filter((s) => s.dex === "trt").length === 4, "4 treated, 4 untreated");
ok(summary.genes_tested === volcano.rows.length, `volcano has every tested gene (${volcano.rows.length})`);
ok(summary.genes_after_filter === summary.genes_tested + summary.independent_filtered, "filtered + tested adds up");
const sig = volcano.rows.filter((r) => r[6] < 0.05);
ok(sig.length === summary.sig_005, `padj<0.05 count matches (${sig.length})`);
ok(sig.filter((r) => r[3] > 0).length === summary.up_005 && sig.filter((r) => r[3] < 0).length === summary.down_005, "up/down split matches");
// fold changes are stored to 3 decimals, so a gene sitting exactly at |LFC| = 1.000x can fall on either side
ok(Math.abs(sig.filter((r) => Math.abs(r[3]) > 1).length - summary.sig_005_lfc1) <= 5, "padj<0.05 & |LFC|>1 count matches (±5 for rounding)");
ok(volcano.rows.every((r) => r[6] >= r[5] - 1e-12 && r[6] <= 1), "padj ≥ p and ≤ 1 for every gene");
// apeglm can move a fold change *away* from zero when one group has almost no counts (ALOX15B here); it should stay rare
ok(volcano.rows.filter((r) => Math.abs(r[4]) > Math.abs(r[3]) + 0.05).length <= 3, "shrunken LFC larger than raw LFC for at most 3 genes");
ok(hist.counts.reduce((a, b) => a + b, 0) >= summary.genes_tested, "p-value histogram covers the tested genes");
ok(hist.counts[0] === summary.raw_p_005, "first histogram bin equals raw p<0.05 count");
ok(top.every((g) => g.padj < 0.05), "top genes are all significant");
ok(top.every((g, i) => i === 0 || g.padj >= top[i - 1].padj), "top genes sorted by padj");
ok(paper.length >= 12 && paper.every((g) => g.padj < 0.05), `paper genes all significant (${paper.length})`);
const crispld2 = paper.find((g) => g.symbol === "CRISPLD2");
ok(crispld2 && crispld2.log2FoldChange > 2, "CRISPLD2 up by more than 4-fold, as in Himes et al.");
ok(pca.percentVar[0] + pca.percentVar[1] <= 100 && pca.points.length === 8, "PCA variance and points sane");
ok(pca.points.filter((p) => p.dex === "trt").every((p) => p.PC1 > 0) && pca.points.filter((p) => p.dex === "untrt").every((p) => p.PC1 < 0), "PC1 separates treated from untreated");
ok(hm.z.every((row) => close(row.reduce((a, b) => a + b, 0), 0, 0.02)), "heatmap z-score rows average to 0");
ok(hm.genes.length === 25 && hm.samples.length === 8, "heatmap is 25 × 8");
ok(norm.sizeFactor.every((f, i) => close(f, summary.size_factors[norm.samples[i]], 1e-3)), "size factors agree between files");
ok(norm.normalized.every((q) => close(q[2], norm.normalized[0][2], 0.1)), "medians line up after normalization");
ok(gene.points.length === 8 && gene.result.padj < 0.05, "gene example has 8 points and is significant");
ok(dist.samples.length === 8 && dist.matrix.every((row, i) => row[i] === 0), "distance matrix has zero diagonal");
ok(disp.rows.length > 1000 && disp.rows.every((r) => r.every((v) => v > 0)), "dispersion sample is positive");
ok(go.terms.length > 0 && go.terms.every((t) => t.padj < 0.05), "GO terms are significant");
ok(counting.genes.length === counting.counts.length && counting.counts.every((c) => c.length === 8), "subset counting table shape");
ok(head.n_genes === summary.genes_in_matrix, "count matrix size matches summary");

ok(Object.values(salmon).every((r) => r.percent_mapped > 85 && r.expected_format === "IU"), "Salmon runs: mapping rate > 85%, unstranded");
ok(txi.runs.every((r, i) => txi.total_counts[i] === salmon[r].num_mapped), "tximport column sums equal Salmon mapped fragments");
ok(txi.pick.find((g) => g.symbol === "ANGPTL7").counts.join() === counting.counts[counting.symbols.indexOf("ANGPTL7")].slice(0, 2).join(), "ANGPTL7 counts agree between Salmon/tximport and the STAR-based subset count");
ok(txi.pick.find((g) => g.symbol === "TNMD").counts.every((c) => c === 0), "TNMD is zero in both Salmon samples");

// "Break it yourself" and "Mistakes" (analysis/break.R) must agree with the course analysis (analysis/run.R)
const br = j("src/data/break.json"), bg = br.grid, bm = br.mistakes;
ok(br.base.sig === summary.sig_005 && br.base.up === summary.up_005 && br.base.down === summary.down_005 && br.base.tested === summary.genes_tested, `break.R reproduces the course result (${br.base.sig} genes, ${br.base.up} up, ${br.base.down} down)`);
ok(bg.all_paired.sig === br.base.sig && bg.all_paired.lost === 0 && bg.all_paired.new === 0, "the default switch setting is the course analysis");
ok(Object.values(bg).every((x) => x.kept + x.lost === br.base.sig && x.kept + x.new === x.sig && x.up + x.down === x.sig), "every scenario: kept + lost = reference, kept + new = its own count, up + down = significant");
ok(Object.keys(bg).length === 8 && Object.values(bg).every((x) => x.markers.length === 7 && x.top.length === 10), "8 scenarios, each with the 7 marker genes and a top 10");
ok(bg.all_paired.markers.every((x) => x.padj < 0.05), "all 7 marker genes are significant in the correct analysis");
ok(bm.p_vs_padj.raw_p === summary.raw_p_005 && bm.p_vs_padj.padj === summary.sig_005, "p vs padj numbers match the course analysis");
ok(bm.lfc_threshold.posthoc === summary.sig_005_lfc1 && bm.lfc_threshold.tested < bm.lfc_threshold.posthoc, `fold-change cut: ${bm.lfc_threshold.posthoc} post hoc = course number; lfcThreshold test gives fewer (${bm.lfc_threshold.tested})`);
ok(/full rank/.test(bm.confounded.error) && /same number of samples and coefficients/.test(bm.one_vs_one.error) && /not integers/.test(bm.normalized_error.error), "the three DESeq2 error messages were captured");
ok(bm.background.genes === go.n_genes && bm.background.universe === go.n_universe && bm.background.terms_tested_bg === go.n_terms, `enrichment with the tested-gene background equals the course GO result (${go.n_terms} terms)`);
ok(bm.shrinkage.rows.every((r) => Math.abs(r.shrunk) < Math.abs(r.lfc)), "shrinkage pulls every listed fold change toward zero");

console.log(failed ? `\n${failed} check(s) failed` : "\nall checks passed");
process.exit(failed ? 1 : 0);
