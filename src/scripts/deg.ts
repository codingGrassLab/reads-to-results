// Parsing of differential-expression result tables (DESeq2, edgeR, limma) and the
// publication figures built from them. No dependencies; used by the results tool page.

export interface Row { id: string; symbol: string; mean: number; lfc: number; p: number; padj: number }
export type Cls = "up" | "down" | "ns";
export interface Thresholds { padj: number; lfc: number; minMean: number }
export interface Parsed { rows: Row[]; columns: string[]; mapping: Mapping; format: string; nInput: number; nDropped: number; delimiter: string }
export interface Mapping { id: string; symbol: string; mean: string; lfc: string; p: string; padj: string }

const known: { format: string; mean: string[]; lfc: string[]; p: string[]; padj: string[] }[] = [
  { format: "DESeq2", mean: ["baseMean"], lfc: ["log2FoldChange"], p: ["pvalue"], padj: ["padj"] },
  { format: "edgeR", mean: ["logCPM"], lfc: ["logFC"], p: ["PValue"], padj: ["FDR"] },
  { format: "limma", mean: ["AveExpr"], lfc: ["logFC"], p: ["P.Value"], padj: ["adj.P.Val"] },
];
const idNames = ["gene_id", "gene", "geneid", "ensembl", "ensembl_id", "id", "row", "rowname", "rownames", "feature", "", "x", "genes"];
const symbolNames = ["symbol", "gene_symbol", "gene_name", "genename", "hgnc_symbol", "name"];

export function detectDelimiter(text: string): string {
  const head = text.slice(0, 5000).split(/\r?\n/).slice(0, 5);
  const score = (d: string) => head.map((l) => l.split(d).length).reduce((a, b) => a + b, 0);
  // tabs and semicolons are unambiguous; commas can also be decimal separators, so they are the fallback
  if (score("\t") > head.length * 1.5) return "\t";
  if (score(";") > head.length * 1.5) return ";";
  return ",";
}

/** Minimal CSV parsing with quoted fields. */
export function parseTable(text: string, delimiter?: string): { columns: string[]; data: string[][]; delimiter: string } {
  const d = delimiter ?? detectDelimiter(text);
  const lines = text.split(/\r?\n/).filter((l) => l.trim().length);
  const split = (line: string) => {
    const out: string[] = []; let cur = "", q = false;
    for (let i = 0; i < line.length; i++) {
      const c = line[i];
      if (q) { if (c === '"') { if (line[i + 1] === '"') { cur += '"'; i++; } else q = false; } else cur += c; }
      else if (c === '"') q = true;
      else if (c === d) { out.push(cur); cur = ""; }
      else cur += c;
    }
    out.push(cur);
    return out;
  };
  const columns = split(lines[0]).map((c) => c.trim());
  const data = lines.slice(1).map(split);
  // R's write.csv leaves the row-name column unnamed: header has one field fewer than the rows
  if (data.length && data[0].length === columns.length + 1) columns.unshift("");
  return { columns, data, delimiter: d };
}

export function autoMap(columns: string[]): { mapping: Mapping; format: string } | null {
  const lower = columns.map((c) => c.toLowerCase());
  const find = (names: string[]) => { for (const n of names) { const i = lower.indexOf(n.toLowerCase()); if (i >= 0) return columns[i]; } return ""; };
  for (const k of known) {
    const mean = find(k.mean), lfc = find(k.lfc), p = find(k.p), padj = find(k.padj);
    if (lfc && padj) {
      const id = find(idNames) || columns[0];
      return { mapping: { id, symbol: find(symbolNames), mean, lfc, p, padj }, format: k.format };
    }
  }
  return null;
}

export function toRows(columns: string[], data: string[][], m: Mapping): { rows: Row[]; nDropped: number } {
  const ix = (c: string) => columns.indexOf(c);
  const iId = ix(m.id), iSym = m.symbol ? ix(m.symbol) : -1, iMean = m.mean ? ix(m.mean) : -1, iLfc = ix(m.lfc), iP = m.p ? ix(m.p) : -1, iPadj = ix(m.padj);
  // decimal commas (write.csv2, European Excel): "495,567" with no dot is read as 495.567
  const num = (s: string | undefined) => { if (s === undefined || s === "" || s === "NA" || s === "NaN" || s === "null") return NaN; let t = s.trim(); if (/^-?\d+,\d+(e[+-]?\d+)?$/i.test(t)) t = t.replace(",", "."); const v = Number(t); return Number.isFinite(v) ? v : NaN; };
  const rows: Row[] = []; let nDropped = 0;
  for (const r of data) {
    const lfc = num(r[iLfc]), padj = num(r[iPadj]);
    if (!Number.isFinite(lfc)) { nDropped++; continue; }
    rows.push({ id: (r[iId] ?? "").replace(/^"|"$/g, ""), symbol: iSym >= 0 ? r[iSym] ?? "" : "", mean: iMean >= 0 ? num(r[iMean]) : NaN, lfc, p: iP >= 0 ? num(r[iP]) : NaN, padj });
  }
  return { rows, nDropped };
}

export function classify(r: Row, t: Thresholds): Cls {
  if (!Number.isFinite(r.padj) || r.padj >= t.padj || Math.abs(r.lfc) < t.lfc) return "ns";
  if (Number.isFinite(r.mean) && t.minMean > 0 && r.mean < t.minMean) return "ns";
  return r.lfc > 0 ? "up" : "down";
}

export function label(r: Row): string { return r.symbol || r.id; }

export function significantCSV(rows: Row[], t: Thresholds, meanName = "baseMean"): string {
  const sig = rows.map((r) => ({ r, c: classify(r, t) })).filter((x) => x.c !== "ns").sort((a, b) => (a.r.padj - b.r.padj) || (Math.abs(b.r.lfc) - Math.abs(a.r.lfc)));
  const esc = (s: string) => (/[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s);
  const head = ["gene_id", "symbol", meanName, "log2FoldChange", "pvalue", "padj", "direction"].join(",");
  const lines = sig.map(({ r, c }) => [esc(r.id), esc(r.symbol), fmtNum(r.mean), fmtNum(r.lfc), fmtNum(r.p), fmtNum(r.padj), c].join(","));
  return [head, ...lines].join("\n") + "\n";
}

export function fmtNum(v: number): string { if (!Number.isFinite(v)) return "NA"; if (v !== 0 && Math.abs(v) < 1e-3) return v.toExponential(3); return String(+v.toPrecision(6)); }

// ── figures ────────────────────────────────────────────────────────────────
// Okabe–Ito colours: distinguishable with the common forms of colour blindness.
export const palette = { up: "#D55E00", down: "#0072B2", ns: "#BBBBBB", ink: "#000000", grid: "#DDDDDD" };

export interface FigureOpts {
  kind: "volcano" | "ma";
  widthMm: number;          // 89 (single column) or 183 (double)
  heightMm: number;
  fontPt: number;           // 7 or 8
  labels: string[];         // ids or symbols to label
  labelTop: number;         // label the N most significant genes
  title?: string;
  thresholds: Thresholds;
  meanLabel: string;        // "mean of normalized counts" / "logCPM" / "average expression"
  meanIsLog: boolean;
  showThresholdLines: boolean;
  yCap?: number;            // cap for -log10(padj)
}

const MM = 72 / 25.4; // pt per mm

/** A publication figure as an SVG string, with physical size in mm and text in pt. */
export function figureSVG(rows: Row[], o: FigureOpts): string {
  const W = o.widthMm * MM, H = o.heightMm * MM, f = o.fontPt;
  const m = { l: f * 4.2, r: f * 0.8, t: (o.title ? f * 1.6 : 0) + f * 1.9, b: f * 3.4 };
  const t = o.thresholds;
  const cls = rows.map((r) => classify(r, t));
  const yVal = (r: Row) => (o.kind === "volcano" ? -Math.log10(Math.max(r.padj, 1e-300)) : r.lfc);
  const xVal = (r: Row) => (o.kind === "volcano" ? r.lfc : o.meanIsLog ? r.mean : Math.log10(Math.max(r.mean, 0.01)));
  const xs = rows.map(xVal).filter(Number.isFinite), ysRaw = rows.map(yVal).filter(Number.isFinite);
  const yCap = o.yCap ?? Infinity;
  const ys = ysRaw.map((v) => Math.min(v, yCap));
  let x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys);
  if (o.kind === "volcano") { const a = Math.max(Math.abs(x0), Math.abs(x1)); x0 = -a; x1 = a; y0 = 0; }
  const padX = (x1 - x0) * 0.04, padY = (y1 - y0) * 0.05; x0 -= padX; x1 += padX; y1 += padY; if (o.kind === "ma") y0 -= padY;
  const sx = (v: number) => m.l + ((v - x0) / (x1 - x0)) * (W - m.l - m.r);
  const sy = (v: number) => H - m.b - ((v - y0) / (y1 - y0)) * (H - m.b - m.t);
  const fontAt = (size: number) => `font-family="Arial, Helvetica, sans-serif" font-size="${+size.toFixed(2)}"`;
  const font = fontAt(f);
  const parts: string[] = [];
  parts.push(`<svg xmlns="http://www.w3.org/2000/svg" width="${o.widthMm}mm" height="${o.heightMm}mm" viewBox="0 0 ${W.toFixed(2)} ${H.toFixed(2)}">`);
  parts.push(`<rect width="100%" height="100%" fill="#fff"/>`);
  // points: grey first, then coloured
  const r0 = f * 0.13, r1 = f * 0.17;
  for (const pass of ["ns", "down", "up"] as Cls[]) {
    const g: string[] = [];
    rows.forEach((r, i) => { if (cls[i] !== pass) return; const xv = xVal(r), yv = yVal(r); if (!Number.isFinite(xv) || !Number.isFinite(yv)) return; g.push(`<circle cx="${sx(xv).toFixed(1)}" cy="${sy(Math.min(yv, yCap)).toFixed(1)}" r="${pass === "ns" ? r0 : r1}"/>`); });
    parts.push(`<g fill="${palette[pass]}" fill-opacity="${pass === "ns" ? 0.5 : 0.85}">${g.join("")}</g>`);
  }
  // threshold guides
  if (o.showThresholdLines) {
    const dash = `stroke="#777" stroke-width="0.5" stroke-dasharray="2,2"`;
    if (o.kind === "volcano") {
      parts.push(`<line x1="${m.l}" x2="${W - m.r}" y1="${sy(-Math.log10(t.padj))}" y2="${sy(-Math.log10(t.padj))}" ${dash}/>`);
      if (t.lfc > 0) for (const s of [-1, 1]) parts.push(`<line x1="${sx(s * t.lfc)}" x2="${sx(s * t.lfc)}" y1="${m.t}" y2="${H - m.b}" ${dash}/>`);
    } else {
      parts.push(`<line x1="${m.l}" x2="${W - m.r}" y1="${sy(0)}" y2="${sy(0)}" stroke="#333" stroke-width="0.5"/>`);
      if (t.lfc > 0) for (const s of [-1, 1]) parts.push(`<line x1="${m.l}" x2="${W - m.r}" y1="${sy(s * t.lfc)}" y2="${sy(s * t.lfc)}" ${dash}/>`);
    }
  }
  // axes
  parts.push(`<g stroke="#000" stroke-width="0.6" fill="none"><path d="M${m.l},${m.t}V${H - m.b}H${W - m.r}"/></g>`);
  const tx = niceTicks(x0, x1, 6), ty = niceTicks(y0, y1, 6);
  for (const v of tx) parts.push(`<line x1="${sx(v)}" x2="${sx(v)}" y1="${H - m.b}" y2="${H - m.b + 2.5}" stroke="#000" stroke-width="0.6"/><text x="${sx(v)}" y="${H - m.b + 2.5 + f}" text-anchor="middle" ${font}>${o.kind === "ma" && !o.meanIsLog ? fmtPow(v) : fmtTick(v)}</text>`);
  for (const v of ty) parts.push(`<line x1="${m.l - 2.5}" x2="${m.l}" y1="${sy(v)}" y2="${sy(v)}" stroke="#000" stroke-width="0.6"/><text x="${m.l - 3.5}" y="${sy(v) + f * 0.35}" text-anchor="end" ${font}>${fmtTick(v)}</text>`);
  const xl = o.kind === "volcano" ? "log<tspan baseline-shift=\"sub\" font-size=\"" + f * 0.7 + "\">2</tspan> fold change" : o.meanLabel;
  const yl = o.kind === "volcano" ? "−log10(adjusted p)" : "log2 fold change";
  parts.push(`<text x="${(m.l + W - m.r) / 2}" y="${H - f * 0.6}" text-anchor="middle" ${font}>${xl}</text>`);
  parts.push(`<text transform="translate(${f * 1.1},${(m.t + H - m.b) / 2}) rotate(-90)" text-anchor="middle" ${font}>${yl.replace("log10", 'log<tspan baseline-shift="sub" font-size="' + f * 0.7 + '">10</tspan>').replace("log2", 'log<tspan baseline-shift="sub" font-size="' + f * 0.7 + '">2</tspan>')}</text>`);
  if (o.title) parts.push(`<text x="${m.l}" y="${f * 1.2}" ${font} font-weight="bold">${esc(o.title)}</text>`);
  // legend: one row above the plot area, so it can never cover data
  const nUp = cls.filter((c) => c === "up").length, nDown = cls.filter((c) => c === "down").length;
  const leg = [[`up (n = ${nUp})`, palette.up], [`down (n = ${nDown})`, palette.down], [`not significant`, palette.ns]] as const;
  let lx = m.l; const ly = m.t - f * 0.7;
  for (const [txt, col] of leg) { parts.push(`<circle cx="${lx + r1}" cy="${ly - f * 0.3}" r="${r1}" fill="${col}"/><text x="${lx + r1 * 2 + f * 0.3}" y="${ly}" ${fontAt(f * 0.9)}>${txt}</text>`); lx += r1 * 2 + f * 0.3 + txt.length * f * 0.5 + f * 1.2; }
  // labels, with simple collision avoidance
  const want = new Set(o.labels.map((s) => s.toLowerCase()));
  const order = rows.map((r, i) => ({ r, i })).filter(({ r, i }) => cls[i] !== "ns" || want.has(r.id.toLowerCase()) || want.has(r.symbol.toLowerCase())).sort((a, b) => a.r.padj - b.r.padj);
  const chosen = order.filter(({ r }, k) => k < o.labelTop || want.has(r.id.toLowerCase()) || want.has(r.symbol.toLowerCase()));
  const placed: { x0: number; y0: number; x1: number; y1: number }[] = [];
  const hit = (b: { x0: number; y0: number; x1: number; y1: number }) => placed.some((p) => b.x0 < p.x1 + 0.5 && b.x1 > p.x0 - 0.5 && b.y0 < p.y1 + 0.5 && b.y1 > p.y0 - 0.5) || b.x0 < m.l || b.x1 > W - m.r || b.y0 < m.t || b.y1 > H - m.b;
  for (const { r } of chosen) {
    const xv = xVal(r), yv = Math.min(yVal(r), yCap); if (!Number.isFinite(xv) || !Number.isFinite(yv)) continue;
    const text = label(r), w = text.length * f * 0.5, h = f * 0.9, px = sx(xv), py = sy(yv);
    // candidate anchor positions around the point, nearest first
    const cands: [number, number][] = [];
    for (const d of [1, 2, 3.2, 4.6, 6.2, 8]) for (const [ax, ay] of [[1, -1], [-1, -1], [1, 1], [-1, 1], [0, -1], [0, 1], [1, 0], [-1, 0]]) cands.push([px + ax * f * 0.5 * d + (ax < 0 ? -w : ax === 0 ? -w / 2 : 0), py + ay * f * 0.55 * d + (ay > 0 ? h : ay === 0 ? h / 2 : 0)]);
    let box = null as null | { x0: number; y0: number; x1: number; y1: number };
    for (const [x, y] of cands) { const b = { x0: x, y0: y - h, x1: x + w, y1: y }; if (!hit(b)) { box = b; break; } }
    if (!box) continue;
    placed.push(box);
    const ax = Math.abs(box.x0 - px) < Math.abs(box.x1 - px) ? box.x0 : box.x1, ay = (box.y0 + box.y1) / 2;
    if (Math.hypot(ax - px, ay - py) > f * 0.8) parts.push(`<line x1="${px.toFixed(1)}" y1="${py.toFixed(1)}" x2="${ax.toFixed(1)}" y2="${ay.toFixed(1)}" stroke="#555" stroke-width="0.3"/>`);
    parts.push(`<text x="${box.x0.toFixed(1)}" y="${(box.y1 - h * 0.2).toFixed(1)}" ${fontAt(f * 0.85)} font-style="italic">${esc(text)}</text>`);
  }
  parts.push("</svg>");
  return parts.join("");
}

export function figureCaption(rows: Row[], o: FigureOpts, source: string): string {
  const t = o.thresholds, cls = rows.map((r) => classify(r, t));
  const nUp = cls.filter((c) => c === "up").length, nDown = cls.filter((c) => c === "down").length;
  const crit = `adjusted p < ${t.padj}${t.lfc > 0 ? ` and |log2 fold change| > ${t.lfc}` : ""}${t.minMean > 0 ? ` and ${o.meanLabel} ≥ ${t.minMean}` : ""}`;
  const nTested = rows.filter((r) => Number.isFinite(r.padj)).length;
  if (o.kind === "volcano") return `Volcano plot of ${nTested.toLocaleString("en-US")} genes tested for differential expression (${source}). Each point is one gene; the x axis shows the log2 fold change and the y axis the −log10 Benjamini–Hochberg-adjusted p-value. Genes meeting ${crit} are coloured orange (up-regulated, n = ${nUp}) or blue (down-regulated, n = ${nDown}); dashed lines mark the thresholds.${o.labelTop > 0 ? ` The ${o.labelTop} genes with the smallest adjusted p-values are labelled.` : ""}`;
  return `MA plot of ${rows.length.toLocaleString("en-US")} genes (${source}; ${nTested.toLocaleString("en-US")} with an adjusted p-value). Each point is one gene; the x axis shows the ${o.meanLabel} and the y axis the log2 fold change. Genes meeting ${crit} are coloured orange (up-regulated, n = ${nUp}) or blue (down-regulated, n = ${nDown}).`;
}

function niceTicks(a: number, b: number, n: number): number[] {
  const span = b - a || 1, step0 = Math.pow(10, Math.floor(Math.log10(span / n))), err = span / n / step0;
  const step = step0 * (err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1);
  const out: number[] = []; for (let v = Math.ceil(a / step) * step; v <= b + 1e-9; v += step) out.push(+v.toFixed(10)); return out;
}
const fmtTick = (v: number) => (Math.abs(v) >= 1000 ? v.toLocaleString("en-US") : String(+v.toFixed(3)));
const fmtPow = (v: number) => (Number.isInteger(v) ? (v >= 0 && v <= 6 ? String(10 ** v) : `10^${v}`) : "");
const esc = (s: string) => s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");

/** Rasterise an SVG string to a PNG blob at the given dpi. */
export async function svgToPng(svg: string, widthMm: number, heightMm: number, dpi = 300): Promise<Blob> {
  const img = new Image();
  const url = URL.createObjectURL(new Blob([svg], { type: "image/svg+xml;charset=utf-8" }));
  await new Promise<void>((res, rej) => { img.onload = () => res(); img.onerror = () => rej(new Error("svg render failed")); img.src = url; });
  const w = Math.round((widthMm / 25.4) * dpi), h = Math.round((heightMm / 25.4) * dpi);
  const c = document.createElement("canvas"); c.width = w; c.height = h;
  const ctx = c.getContext("2d")!; ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, w, h); ctx.drawImage(img, 0, 0, w, h);
  URL.revokeObjectURL(url);
  return new Promise((res) => c.toBlob((b) => res(b!), "image/png"));
}

export function download(name: string, data: Blob | string, type = "text/plain") {
  const blob = data instanceof Blob ? data : new Blob([data], { type });
  const a = document.createElement("a"); a.href = URL.createObjectURL(blob); a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 2000);
}
