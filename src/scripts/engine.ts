// DESeq2 in the browser: webR (R compiled to WebAssembly) with DESeq2 from the Bioconductor
// r-universe WebAssembly builds. Everything runs inside the visitor's browser; the count
// matrix never leaves it. Loaded only when the visitor asks for it (tens of MB, cached).

export const WEBR_URL = "https://webr.r-wasm.org/v0.6.0/webr.mjs";
export const REPOS = {
  bioc: "https://bioc.r-universe.dev",      // DESeq2, apeglm, … built for webR 0.6 / R 4.6
  cran: "https://repo.r-wasm.org/",         // CRAN packages built by the webR project
  locfit: "/wasm/",                          // locfit (DESeq2 dependency without a public wasm build), hosted with the site
};

export type Progress = (msg: string, done?: boolean) => void;
let webR: any = null, ready: Promise<any> | null = null;

/** Start R and install DESeq2 (+ apeglm). Cached by the browser after the first time. */
export function loadEngine(onProgress: Progress): Promise<any> {
  if (ready) return ready;
  ready = (async () => {
    const t0 = performance.now();
    const step = (m: string) => onProgress(`${((performance.now() - t0) / 1000).toFixed(0)} s · ${m}`);
    step("downloading R (webR 0.6)…");
    const mod = await import(/* @vite-ignore */ WEBR_URL);
    const R = new mod.WebR({ repoUrl: REPOS.cran, channelType: mod.ChannelType.PostMessage });
    await R.init();
    (window as any).__R = R;   // for debugging in the browser console
    step(`R ${await R.evalRString('paste(R.version$major, R.version$minor, sep = ".")')} is running; installing locfit…`);
    await R.installPackages(["locfit"], { repos: [location.origin + REPOS.locfit], mount: true });
    step("installing DESeq2 and its dependencies (this is the big download)…");
    await R.installPackages(["DESeq2"], { repos: [REPOS.bioc, REPOS.cran], mount: true });
    step("installing apeglm…");
    try { await R.installPackages(["apeglm"], { repos: [REPOS.bioc, REPOS.cran], mount: true }); } catch { /* optional */ }
    step("loading DESeq2…");
    await R.evalRVoid("suppressPackageStartupMessages(library(DESeq2))");
    const v = await R.evalRString('as.character(packageVersion("DESeq2"))');
    const hasApeglm = await R.evalRBoolean('requireNamespace("apeglm", quietly = TRUE)');
    onProgress(`ready: DESeq2 ${v}${hasApeglm ? " with apeglm" : ""} in R ${await R.evalRString('paste(R.version$major, R.version$minor, sep = ".")')} (webR)`, true);
    webR = R; (R as any).__version = v; (R as any).__apeglm = hasApeglm;
    return R;
  })();
  ready.catch(() => { ready = null; });
  return ready;
}

export interface RunConfig { countsCSV: string; samplesCSV: string; idCol?: string; sampleCol: string; conditionCol: string; refLevel: string; trtLevel: string; batchCol: string; minCounts: number; minSamples: number; alpha: number; shrink: boolean }
export interface RunResult { resultsCSV: string; log: string; nGenes: number; nAfterFilter: number; nTested: number; nUp: number; nDown: number; sizeFactors: Record<string, number>; versions: { deseq2: string; r: string; apeglm: string } }

const rstr = (s: string) => `"${s.replace(/\\/g, "\\\\").replace(/"/g, '\\"')}"`;

/** The R code that runs (reads /counts.csv and /samples.csv, writes /results.csv). */
export function deseq2Code(c: Omit<RunConfig, "countsCSV" | "samplesCSV">): string {
  if (c.idCol === undefined) c = { ...c, idCol: "" };
  const design = c.batchCol ? `~ ${c.batchCol} + ${c.conditionCol}` : `~ ${c.conditionCol}`;
  const code = `
    log <- character(0); note <- function(...) log <<- c(log, paste0(...))
    raw <- read.csv("/counts.csv", check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
    idc <- ${c.idCol ? `match(${rstr(c.idCol)}, names(raw))` : "1L"}; if (is.na(idc)) idc <- 1L
    ids <- trimws(raw[[idc]]); raw <- raw[-idc]
    note("Gene IDs from the column ", ${rstr(c.idCol || "(first column)")})
    samples <- read.csv("/samples.csv", stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
    samples <- samples[!is.na(samples[[${rstr(c.sampleCol)}]]) & samples[[${rstr(c.sampleCol)}]] != "", , drop = FALSE]
    cv <- samples[[${rstr(c.conditionCol)}]]
    drop <- is.na(cv) | trimws(cv) == ""
    if (any(drop)) { note("Left out ", sum(drop), " samples with no ${c.conditionCol.replace(/[`"\$]/g, "")} value"); samples <- samples[!drop, , drop = FALSE]; cv <- cv[!drop] }
    if (!all(c(${rstr(c.refLevel)}, ${rstr(c.trtLevel)}) %in% cv)) stop("Condition levels found: ", paste(unique(cv), collapse = ", "), ". Expected ${c.refLevel.replace(/[`"\$]/g, "")} and ${c.trtLevel.replace(/[`"\$]/g, "")}.")
    # internal, syntactic names: the condition and batch columns become condition / batch and their levels make.names() tokens
    lv <- unique(cv); safe <- make.names(lv, unique = TRUE); names(safe) <- lv
    ref <- safe[[${rstr(c.refLevel)}]]; trt <- safe[[${rstr(c.trtLevel)}]]
    cd <- data.frame(condition = relevel(factor(unname(safe[cv])), ref = ref), row.names = samples[[${rstr(c.sampleCol)}]])
    ${c.batchCol ? `bv <- samples[[${rstr(c.batchCol)}]]; bv[is.na(bv) | bv == ""] <- "none"; cd$batch <- factor(make.names(bv))` : ""}
    missing <- setdiff(rownames(cd), colnames(raw))
    if (length(missing)) stop("Samples in the sheet but not in the count matrix: ", paste(missing, collapse = ", "))
    extra <- setdiff(colnames(raw), rownames(cd))
    if (length(extra)) note("Columns of the matrix not used: ", paste(head(extra, 8), collapse = ", "), if (length(extra) > 8) ", …" else "")
    counts <- vapply(rownames(cd), function(s) suppressWarnings(as.numeric(gsub(",", ".", raw[[s]]))), numeric(length(ids)))
    if (is.null(dim(counts))) counts <- matrix(counts, ncol = nrow(cd), dimnames = list(NULL, rownames(cd)))
    ok <- !is.na(ids) & ids != ""
    counts <- counts[ok, , drop = FALSE]; ids <- ids[ok]
    if (any(is.na(counts))) { note(sum(is.na(counts)), " empty or non-numeric count cells set to 0"); counts[is.na(counts)] <- 0 }
    if (any(counts < 0)) stop("The matrix has negative values: these are not raw counts (log-transformed or normalized data cannot go into DESeq2).")
    if (any(counts %% 1 != 0)) { fr <- mean(counts %% 1 != 0); if (fr > 0.5 && max(counts) < 50) stop("Most values are small decimals: this looks like normalized or log-scale data (TPM, FPKM, log2CPM), not raw counts."); note("Non-integer values (", signif(100 * fr, 2), "% of cells, e.g. Salmon estimated counts) rounded to integers."); counts <- round(counts) }
    if (anyDuplicated(ids)) { note(sum(duplicated(ids)), " duplicated gene IDs: their counts were summed"); counts <- rowsum(counts, ids, reorder = FALSE) } else rownames(counts) <- ids
    storage.mode(counts) <- "integer"
    note("Loaded ", nrow(counts), " genes x ", ncol(counts), " samples; design ", ${rstr(design)}, "; ", ${rstr(c.trtLevel)}, " (n = ", sum(cv == ${rstr(c.trtLevel)}), ") vs ", ${rstr(c.refLevel)}, " (n = ", sum(cv == ${rstr(c.refLevel)}), ")")
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = cd, design = ${c.batchCol ? "~ batch + condition" : "~ condition"})
    keep <- rowSums(counts(dds) >= ${c.minCounts}) >= ${c.minSamples}
    nAfter <- sum(keep); dds <- dds[keep, ]
    note("Kept ", nAfter, " genes with >= ${c.minCounts} counts in >= ${c.minSamples} samples")
    dds <- DESeq(dds, parallel = FALSE, quiet = TRUE)
    res <- results(dds, contrast = c("condition", trt, ref), alpha = ${c.alpha}, parallel = FALSE)
    shr <- NA_character_
    ${c.shrink ? `if (requireNamespace("apeglm", quietly = TRUE)) { r2 <- lfcShrink(dds, coef = paste0("condition_", trt, "_vs_", ref), type = "apeglm", quiet = TRUE); res$shrunk_log2FoldChange <- r2$log2FoldChange; shr <- as.character(packageVersion("apeglm")) } else note("apeglm not available: no shrunken fold changes")` : ""}
    df <- as.data.frame(res); df <- cbind(gene_id = rownames(df), df); df <- df[order(df$padj, df$pvalue), ]
    write.csv(df, "/results.csv", row.names = FALSE)
    sig <- !is.na(df$padj) & df$padj < ${c.alpha}
    c(paste0("nGenes=", nrow(counts)), paste0("nAfter=", nAfter), paste0("nTested=", sum(!is.na(df$padj))), paste0("nUp=", sum(sig & df$log2FoldChange > 0)), paste0("nDown=", sum(sig & df$log2FoldChange < 0)),
      paste0("sf=", paste(colnames(dds), signif(sizeFactors(dds), 5), sep = ":", collapse = ";")), paste0("shr=", shr), paste0("log=", paste(log, collapse = " | ")))
  `;
  return code;
}

/** Run DESeq2 on a count matrix and sample sheet (both as CSV text). */
export async function runDESeq2(c: RunConfig, onProgress: Progress): Promise<RunResult> {
  const R = await loadEngine(onProgress);
  const enc = new TextEncoder();
  await R.FS.writeFile("/counts.csv", enc.encode(c.countsCSV));
  await R.FS.writeFile("/samples.csv", enc.encode(c.samplesCSV));
  const code = deseq2Code(c);
  onProgress("running DESeq2…");
  const shelter = await new R.Shelter();
  try {
    const out = await shelter.captureR(code, { withAutoprint: false, captureStreams: true, captureConditions: false });
    const lines = (await out.result.toArray()) as string[];
    const kv: Record<string, string> = {}; for (const l of lines) { const i = l.indexOf("="); kv[l.slice(0, i)] = l.slice(i + 1); }
    const resultsCSV = new TextDecoder().decode(await R.FS.readFile("/results.csv"));
    const sizeFactors = Object.fromEntries((kv.sf || "").split(";").filter(Boolean).map((x) => { const [n, v] = x.split(":"); return [n, Number(v)]; }));
    return { resultsCSV, log: (kv.log || "").split(" | ").join("\n"), nGenes: +kv.nGenes, nAfterFilter: +kv.nAfter, nTested: +kv.nTested, nUp: +kv.nUp, nDown: +kv.nDown,
      sizeFactors, versions: { deseq2: (R as any).__version, r: await R.evalRString('paste(R.version$major, R.version$minor, sep = ".")'), apeglm: kv.shr === "NA" ? "" : kv.shr } };
  } catch (e: any) {
    const msg = e?.message ?? String(e);
    throw new Error(msg.replace(/^Error in [^:]*: /, "").trim());
  } finally { shelter.purge(); }
}
