// Human Ensembl gene → symbol table used by the results and enrichment tools.
// HGNC (the naming authority) first; genes HGNC does not name keep the org.Hs.eg.db symbol
// already in the table (written by analysis/symbols.R).
//   node scripts/build-symbols.mjs hgnc_complete_set.txt
// HGNC file: https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
import { readFileSync, writeFileSync } from "node:fs";

const file = process.argv[2];
if (!file) { console.error("usage: node scripts/build-symbols.mjs hgnc_complete_set.txt"); process.exit(1); }
const out = "public/data/symbols_human.json";
const old = JSON.parse(readFileSync(out, "utf8"));
const lines = readFileSync(file, "utf8").split(/\r?\n/).filter(Boolean);
const head = lines[0].split("\t"), iS = head.indexOf("symbol"), iE = head.indexOf("ensembl_gene_id"), iSt = head.indexOf("status");
const hgnc = {};
for (const l of lines.slice(1)) { const f = l.split("\t"); if (f[iE] && f[iSt] === "Approved") hgnc[f[iE]] = f[iS]; }
const merged = { ...old, ...hgnc };
let changed = 0, added = 0; for (const [k, v] of Object.entries(hgnc)) { if (!(k in old)) added++; else if (old[k] !== v) changed++; }
writeFileSync(out, JSON.stringify(Object.fromEntries(Object.entries(merged).sort())));
console.log(`HGNC ${Object.keys(hgnc).length} Ensembl IDs; ${changed} symbols corrected, ${added} added; ${Object.keys(merged).length} in total.`);
