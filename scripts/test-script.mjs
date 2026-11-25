// Generates a DESeq2 script for the airway data with the script generator and runs it in R.
// The result must match the site's own analysis (4,081 genes at padj < 0.05, 951 with |LFC| > 1).
// Run with: node scripts/test-script.mjs   (needs R with the packages; takes ~1 minute)
import { buildScript, defaults } from "../src/scripts/rscript.mjs";
import { mkdirSync, writeFileSync, readFileSync, existsSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { gunzipSync } from "node:zlib";

const dir = new URL("../analysis/cache/script-test/", import.meta.url);
mkdirSync(dir, { recursive: true });
const counts = gunzipSync(readFileSync(new URL("../public/data/airway_counts.csv.gz", import.meta.url))).toString("utf8");
writeFileSync(new URL("counts.csv", dir), counts.replace(/^""/, '"gene_id"'));
const samples = readFileSync(new URL("../public/data/airway_samples.csv", import.meta.url), "utf8").split(/\r?\n/).filter(Boolean);
// keep Run, cell, dex
const head = samples[0].split(",").map((s) => s.replace(/"/g, ""));
const iRun = head.indexOf("Run"), iCell = head.indexOf("cell"), iDex = head.indexOf("dex");
writeFileSync(new URL("samples.csv", dir), ["sample,cell,dex", ...samples.slice(1).map((l) => { const f = l.split(",").map((s) => s.replace(/"/g, "")); return [f[iRun], f[iCell], f[iDex]].join(","); })].join("\n") + "\n");

const cfg = { ...defaults, input: "counts", countsFile: "counts.csv", samplesFile: "samples.csv", outDir: "out",
  sampleCol: "sample", conditionCol: "dex", refLevel: "untrt", trtLevel: "trt", batchCol: "cell", minCounts: 10, minSamples: 4, alpha: 0.05, lfc: 1,
  plots: { pca: true, ma: true, volcano: true, heatmap: true, counts: true, phist: true }, countGenes: ["CRISPLD2", "FKBP5"], enrichment: false };
const script = buildScript(cfg);
writeFileSync(new URL("analysis.R", dir), script);
console.log(`script: ${script.split("\n").length} lines → ${dir.pathname}analysis.R`);

const rscript = ["C:/Program Files/R/R-4.5.2/bin/Rscript.exe", "Rscript"].find((p) => p === "Rscript" || existsSync(p));
let log;
try { log = execFileSync(rscript, ["analysis.R"], { cwd: dir, encoding: "utf8", stdio: ["ignore", "pipe", "pipe"] }); }
catch (e) { console.error(e.stdout, e.stderr); process.exit(1); }
const sigLine = log.split("\n").find((l) => l.startsWith("Significant"));
console.log(log.split("\n").filter((l) => /^(Genes before|Loaded|Significant|Done)/.test(l)).join("\n"));
const m = /Significant.*: (\d+)\s+up: (\d+)\s+down: (\d+)/.exec(sigLine ?? "");
const summary = JSON.parse(readFileSync(new URL("../src/data/summary.json", import.meta.url), "utf8"));
const ok = m && +m[1] === summary.sig_005_lfc1 && +m[2] === summary.up_005_lfc1 && +m[3] === summary.down_005_lfc1;
const outputs = ["deseq2_results.csv", "significant_genes_padj0.05_lfc1.csv", "pca.png", "volcano.png", "ma_plot.png", "heatmap_top30.png", "gene_counts.pdf", "sessionInfo.txt"].filter((f) => !existsSync(new URL("out/" + f, dir)));
console.log(ok ? `ok   generated script reproduces the site's numbers (${m[1]} significant, ${m[2]} up, ${m[3]} down)` : `FAIL expected ${summary.sig_005_lfc1}/${summary.up_005_lfc1}/${summary.down_005_lfc1}, got ${sigLine}`);
console.log(outputs.length ? `FAIL missing outputs: ${outputs.join(", ")}` : "ok   all expected output files exist");
process.exit(ok && !outputs.length ? 0 : 1);
