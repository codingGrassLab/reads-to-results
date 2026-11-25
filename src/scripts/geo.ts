// GEO series → count matrix + sample sheet, read straight from NCBI's FTP mirror (it allows
// cross-origin reads). The SOFT "family" file carries every sample's title, characteristics and
// description; the description or supplementary file name usually names the column the submitter
// used in the count matrix, which is how matrix columns are matched to GSM samples.

export const geoDir = (gse: string) => `https://ftp.ncbi.nlm.nih.gov/geo/series/${gse.slice(0, -3)}nnn/${gse}/`;

export interface GeoSample { gsm: string; title: string; source: string; organism: string; chars: Record<string, string>; desc: string[]; files: string[] }
export interface CountFile { name: string; score: number; note: string }

async function textOf(url: string): Promise<string> {
  const r = await fetch(url);
  if (!r.ok) throw new Error(r.status === 404 ? "not found on NCBI" : `NCBI returned ${r.status}`);
  const b = await r.blob();
  const head = new Uint8Array(await b.slice(0, 2).arrayBuffer());
  if (head[0] === 0x1f && head[1] === 0x8b) return new Response(b.stream().pipeThrough(new DecompressionStream("gzip"))).text();
  return b.text();
}
export const fetchGeoText = textOf;

/** Supplementary files of a series, ranked by how much they look like a raw count matrix. */
export async function listCountFiles(gse: string): Promise<CountFile[]> {
  const html = await fetch(geoDir(gse) + "suppl/").then((r) => { if (!r.ok) throw new Error(r.status === 404 ? `${gse} has no supplementary files` : `NCBI returned ${r.status}`); return r.text(); });
  const names = [...html.matchAll(/href="([^"\/?]+)"/g)].map((m) => decodeURIComponent(m[1])).filter((f) => f !== "filelist.txt");
  return names.map((name) => {
    let score = 0, note = "";
    if (!/\.(txt|csv|tsv|tab)(\.gz)?$/i.test(name)) return { name, score: -100, note: /\.xlsx?(\.gz)?$/i.test(name) ? (/fpkm|rpkm|tpm|cpm|norm/i.test(name) ? "Excel, normalized values: not for DESeq2" : "Excel: save as CSV and drop it in") : /_RAW\.tar$/.test(name) ? "per-sample files in one archive" : "not a text table" };
    if (/count/i.test(name)) score += 6;
    if (/raw|reads|featurecounts|htseq|readcount|umi/i.test(name)) score += 3;
    if (/gene/i.test(name)) score += 1;
    if (/fpkm|rpkm|tpm|cpm|norm|log|vst|rlog|abundance|expression_?values/i.test(name)) { score -= 8; note = "looks normalized: DESeq2 needs raw counts"; }
    if (/transcript|isoform/i.test(name)) { score -= 2; note ||= "transcript level"; }
    return { name, score, note };
  }).sort((a, b) => b.score - a.score);
}

/** Samples of a series from its SOFT family file (small: metadata only for sequencing series). */
export async function fetchSamples(gse: string): Promise<GeoSample[]> {
  const txt = await textOf(`${geoDir(gse)}soft/${gse}_family.soft.gz`);
  const out: GeoSample[] = []; let cur: GeoSample | null = null;
  for (const line of txt.split(/\r?\n/)) {
    if (line.startsWith("^SAMPLE")) { cur = { gsm: line.split("=")[1].trim(), title: "", source: "", organism: "", chars: {}, desc: [], files: [] }; out.push(cur); continue; }
    if (line.startsWith("^")) { cur = null; continue; }
    if (!cur || !line.startsWith("!Sample_")) continue;
    const i = line.indexOf(" = "); if (i < 0) continue;
    const k = line.slice(8, i), v = line.slice(i + 3).trim();
    if (k === "title") cur.title = v;
    else if (k === "source_name_ch1") cur.source = v;
    else if (k === "organism_ch1") cur.organism = v;
    else if (k === "characteristics_ch1") { const j = v.indexOf(": "); if (j > 0) cur.chars[v.slice(0, j).trim()] = v.slice(j + 2).trim(); else cur.chars[`characteristic ${Object.keys(cur.chars).length + 1}`] = v; }
    else if (k === "description") cur.desc.push(v);
    else if (/^supplementary_file/.test(k) && v !== "NONE") cur.files.push(v.split("/").pop() ?? "");
  }
  return out;
}

// lower case, runs of non-alphanumerics → "_", and "01" → "1" (titles often pad numbers the headers do not)
const tok = (s: string) => `_${s.toLowerCase().replace(/[^a-z0-9]+/g, "_").replace(/^_|_$/g, "")}_`.replace(/_0+(\d)/g, "_$1");
/** Names a sample could carry as a matrix column header. */
function aliases(s: GeoSample): string[] {
  const a = [s.gsm, s.title, s.source];
  for (const d of s.desc) { a.push(d); const m = /^(?:sample[ _]?(?:name|id)|library(?: name)?|name)\s*[:=]\s*(.+)$/i.exec(d); if (m) a.push(m[1]); }
  for (const f of s.files) { const stem = f.replace(/\.(gz|bz2|zip)$/i, "").replace(/\.[a-z0-9]+$/i, "").replace(/^GSM\d+[_.-]?/, ""); if (stem) a.push(stem); }
  return [...new Set(a.map((x) => x.trim()).filter((x) => x.length >= 2))];
}

/** Match matrix column headers to samples. Exact name matches first, then the longest alias that
 *  appears in the header as whole tokens (so "S1" never matches "S10"). Each sample is used once. */
export function matchColumns(cols: string[], samples: GeoSample[]): Map<string, GeoSample> {
  const pairs: { c: string; s: GeoSample; score: number }[] = [];
  for (const c of cols) {
    const tc = tok(c);
    for (const s of samples) {
      let best = 0;
      for (const al of aliases(s)) {
        const ta = tok(al);
        if (ta.length <= 3) continue;
        if (ta === tc) best = Math.max(best, 1000 + ta.length);
        else if (tc.includes(ta) || ta.includes(tc)) best = Math.max(best, Math.min(ta.length, tc.length));
      }
      if (best) pairs.push({ c, s, score: best });
    }
  }
  pairs.sort((a, b) => b.score - a.score);
  const map = new Map<string, GeoSample>(), used = new Set<string>();
  for (const p of pairs) {
    if (map.has(p.c) || used.has(p.s.gsm)) continue;
    // ambiguous: another sample scores the same for this column
    if (pairs.some((q) => q.c === p.c && q.s !== p.s && q.score === p.score && !used.has(q.s.gsm))) continue;
    map.set(p.c, p.s); used.add(p.s.gsm);
  }
  return map;
}

const csvCell = (x: string) => (/[",\n]/.test(x) ? `"${x.replace(/"/g, '""')}"` : x);
/** Sample sheet: one row per matched column, with the characteristics that differ between samples. */
export function sampleSheet(map: Map<string, GeoSample>, extra: Record<string, Record<string, string>> = {}): { csv: string; columns: string[] } {
  const rows = [...map.entries()];
  const keys = [...new Set(rows.flatMap(([, s]) => Object.keys(s.chars)))];
  const varies = (f: (s: GeoSample) => string) => new Set(rows.map(([, s]) => f(s))).size > 1;
  const charCols = keys.filter((k) => varies((s) => s.chars[k] ?? "")).map((k) => k.replace(/,/g, " "));
  const charKeys = keys.filter((k) => varies((s) => s.chars[k] ?? ""));
  const extraCols = Object.keys(extra);
  const head = ["sample", "gsm", "title", ...extraCols, ...charCols, ...(varies((s) => s.source) ? ["source"] : [])];
  const lines = rows.map(([c, s]) => [c, s.gsm, s.title, ...extraCols.map((e) => extra[e][s.gsm] ?? ""), ...charKeys.map((k) => s.chars[k] ?? ""), ...(head.includes("source") ? [s.source] : [])].map(csvCell).join(","));
  return { csv: [head.map(csvCell).join(","), ...lines].join("\n") + "\n", columns: head };
}

/** g:Profiler organism id from a scientific name ("Mus musculus" → "mmusculus"). */
export const gprofilerOrganism = (sci: string) => { const [g, s] = sci.trim().toLowerCase().split(/\s+/); return g && s ? g[0] + s : ""; };

export interface Availability { level: "ok" | "warn" | "none"; counts: string; meta: string; metaOk: boolean; best?: string }
/** What a series offers for an analysis here: a raw count table? sample metadata beyond the titles? */
export async function availability(gse: string): Promise<Availability> {
  const [files, samples] = await Promise.all([listCountFiles(gse).catch(() => [] as CountFile[]), fetchSamples(gse).catch(() => [] as GeoSample[])]);
  const text = files.filter((f) => f.score > -100);
  const raw = text.filter((f) => !/normalized/.test(f.note));
  let level: Availability["level"] = "ok", counts = "";
  if (raw.length && raw[0].score > 0) counts = `raw count file: ${raw[0].name}`;
  else if (raw.length) { level = "warn"; counts = `text table found (${raw[0].name}), but its name does not say counts: check it`; }
  else if (text.length) { level = "warn"; counts = "only normalized values (FPKM/TPM/CPM): DESeq2 needs raw counts"; }
  else if (files.some((f) => f.note === "Excel: save as CSV and drop it in")) { level = "warn"; counts = "table is an Excel file: download it, save as CSV, drop it in"; }
  else if (files.some((f) => /Excel, normalized/.test(f.note))) { level = "warn"; counts = "only an Excel file of normalized values: DESeq2 needs raw counts"; }
  else if (files.some((f) => /archive/.test(f.note))) { level = "warn"; counts = "only per-sample files in a RAW.tar archive"; }
  else { level = "none"; counts = files.length ? "no count table among the supplementary files" : "no supplementary files"; }
  const keys = [...new Set(samples.flatMap((s) => Object.keys(s.chars)))];
  const varying = keys.filter((k) => new Set(samples.map((s) => s.chars[k] ?? "")).size > 1);
  const metaOk = varying.length > 0;
  const meta = !samples.length ? "sample metadata could not be read" : metaOk ? `sample metadata: ${varying.slice(0, 3).join(", ")}${varying.length > 3 ? "…" : ""}` : "conditions only in the sample titles: check the sample sheet by hand";
  if (!metaOk && level === "ok") level = "warn";
  return { level, counts, meta, metaOk, best: raw[0]?.name };
}
