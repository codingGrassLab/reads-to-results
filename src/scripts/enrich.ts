// Enrichment through the g:Profiler g:GOSt API, and the publication dot plot drawn from it.

export interface Term { source: string; id: string; name: string; p: number; termSize: number; querySize: number; overlap: number; domain: number; fold: number; genes: string[]; query: string }
export interface EnrichOpts { organism: string; sources: string[]; method: "g_SCS" | "fdr" | "bonferroni"; threshold: number; background: string[] | null; noIEA: boolean }

export const organisms: [string, string][] = [["hsapiens", "Homo sapiens (human)"], ["mmusculus", "Mus musculus (mouse)"], ["rnorvegicus", "Rattus norvegicus (rat)"], ["drerio", "Danio rerio (zebrafish)"], ["dmelanogaster", "Drosophila melanogaster"], ["celegans", "Caenorhabditis elegans"], ["scerevisiae", "Saccharomyces cerevisiae"], ["athaliana", "Arabidopsis thaliana"], ["ggallus", "Gallus gallus (chicken)"], ["sscrofa", "Sus scrofa (pig)"], ["btaurus", "Bos taurus (cattle)"]];
export const sourceNames: Record<string, string> = { "GO:BP": "GO biological process", "GO:MF": "GO molecular function", "GO:CC": "GO cellular component", KEGG: "KEGG", REAC: "Reactome", WP: "WikiPathways", TF: "TRANSFAC motifs", MIRNA: "miRTarBase", HPA: "Human Protein Atlas", CORUM: "CORUM complexes", HP: "Human Phenotype Ontology" };

/** Run g:GOSt for one or more gene lists. Returns terms and the API's version string. */
export async function gprofiler(queries: Record<string, string[]>, o: EnrichOpts): Promise<{ terms: Term[]; version: string; failed: Record<string, string[]>; mapped: Record<string, number> }> {
  const body: Record<string, unknown> = { organism: o.organism, query: queries, sources: o.sources, user_threshold: o.threshold, significance_threshold_method: o.method, no_evidences: false, no_iea: o.noIEA, ordered: false };
  if (o.background && o.background.length) { body.domain_scope = "custom_annotated"; body.background = o.background; }
  const r = await fetch("https://biit.cs.ut.ee/gprofiler/api/gost/profile/", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
  if (!r.ok) throw new Error(`g:Profiler returned ${r.status}`);
  const j = await r.json();
  const qmeta = j.meta?.query_metadata?.queries ?? {}, gmeta = j.meta?.genes_metadata?.query ?? {};
  // map the API's internal gene order back to the user's identifiers
  const userIds: Record<string, string[]> = {};
  for (const q of Object.keys(queries)) {
    const ens: string[] = qmeta[q] ?? [];
    const mapping: Record<string, string[]> = gmeta[q]?.mapping ?? {};
    const back = new Map<string, string>();
    for (const [uid, targets] of Object.entries(mapping)) for (const t of targets) if (!back.has(t)) back.set(t, uid);
    userIds[q] = ens.map((e) => back.get(e) ?? e);
  }
  const terms: Term[] = (j.result ?? []).map((x: any) => {
    const ids = userIds[x.query] ?? [];
    const genes = (x.intersections ?? []).map((ev: string[], i: number) => (ev && ev.length ? ids[i] : null)).filter(Boolean) as string[];
    return { source: x.source, id: x.native, name: x.name, p: x.p_value, termSize: x.term_size, querySize: x.query_size, overlap: x.intersection_size, domain: x.effective_domain_size, fold: (x.intersection_size / x.query_size) / (x.term_size / x.effective_domain_size), genes, query: x.query };
  });
  const failed: Record<string, string[]> = {}, mapped: Record<string, number> = {};
  for (const q of Object.keys(queries)) { failed[q] = gmeta[q]?.failed ?? []; mapped[q] = (qmeta[q] ?? []).length; }
  return { terms, version: j.meta?.version ?? "", failed, mapped };
}

export function termsCSV(terms: Term[]): string {
  const esc = (s: string) => (/[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s);
  return ["query,source,term_id,term_name,adjusted_p,overlap,term_size,query_size,background_size,fold_enrichment,genes", ...terms.map((t) => [t.query, t.source, t.id, esc(t.name), t.p.toExponential(3), t.overlap, t.termSize, t.querySize, t.domain, t.fold.toFixed(2), esc(t.genes.join(";"))].join(","))].join("\n") + "\n";
}

export interface DotOpts { widthMm: number; heightMm: number; fontPt: number; topN: number; title?: string; xAxis: "fold" | "p" }
const MM = 72 / 25.4;
const esc = (s: string) => s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");

/** Dot plot of the top terms: x = fold enrichment (or −log10 adjusted p), dot size = genes in term, colour = adjusted p. */
export function dotPlotSVG(terms: Term[], o: DotOpts): string {
  const top = terms.slice().sort((a, b) => a.p - b.p).slice(0, o.topN).reverse();
  const W = o.widthMm * MM, H = o.heightMm * MM, f = o.fontPt;
  const labelW = Math.min(W * 0.55, Math.max(...top.map((t) => t.name.length), 10) * f * 0.5 + f);
  const m = { l: labelW + f * 0.8, r: f * 9, t: (o.title ? f * 1.6 : 0) + f * 0.8, b: f * 3.2 };
  const xv = (t: Term) => (o.xAxis === "fold" ? t.fold : -Math.log10(Math.max(t.p, 1e-300)));
  const xs = top.map(xv); const x1 = Math.max(...xs) * 1.08, x0 = 0;
  const sx = (v: number) => m.l + ((v - x0) / (x1 - x0 || 1)) * (W - m.l - m.r);
  const rowH = (H - m.t - m.b) / Math.max(top.length, 1);
  const sy = (i: number) => H - m.b - (i + 0.5) * rowH;
  const maxOv = Math.max(...top.map((t) => t.overlap), 1), rMax = Math.min(rowH * 0.42, f * 0.9), rMin = f * 0.2;
  const rad = (n: number) => rMin + (rMax - rMin) * Math.sqrt(n / maxOv);
  const lp = top.map((t) => -Math.log10(Math.max(t.p, 1e-300))); const lpMin = Math.min(...lp), lpMax = Math.max(...lp);
  const col = (v: number) => { const k = lpMax === lpMin ? 1 : (v - lpMin) / (lpMax - lpMin); const a = [198, 219, 239], b = [8, 48, 107]; return `rgb(${a.map((c, i) => Math.round(c + (b[i] - c) * k)).join(",")})`; };
  const font = (s: number) => `font-family="Arial, Helvetica, sans-serif" font-size="${+s.toFixed(2)}"`;
  const P: string[] = [`<svg xmlns="http://www.w3.org/2000/svg" width="${o.widthMm}mm" height="${o.heightMm}mm" viewBox="0 0 ${W.toFixed(2)} ${H.toFixed(2)}"><rect width="100%" height="100%" fill="#fff"/>`];
  top.forEach((t, i) => { P.push(`<line x1="${m.l}" x2="${W - m.r}" y1="${sy(i)}" y2="${sy(i)}" stroke="#eee" stroke-width="0.4"/>`); });
  P.push(`<g stroke="#000" stroke-width="0.6" fill="none"><path d="M${m.l},${m.t}V${H - m.b}H${W - m.r}"/></g>`);
  const ticks = niceTicks(x0, x1, 5);
  for (const v of ticks) P.push(`<line x1="${sx(v)}" x2="${sx(v)}" y1="${H - m.b}" y2="${H - m.b + 2.5}" stroke="#000" stroke-width="0.6"/><text x="${sx(v)}" y="${H - m.b + 2.5 + f}" text-anchor="middle" ${font(f)}>${+v.toFixed(2)}</text>`);
  top.forEach((t, i) => {
    const name = t.name.length > 48 ? t.name.slice(0, 46) + "…" : t.name;
    P.push(`<text x="${m.l - f * 0.5}" y="${sy(i) + f * 0.35}" text-anchor="end" ${font(f * 0.9)}>${esc(name)}</text>`);
    P.push(`<circle cx="${sx(xv(t))}" cy="${sy(i)}" r="${rad(t.overlap).toFixed(2)}" fill="${col(lp[i])}" stroke="#333" stroke-width="0.3"/>`);
  });
  const xl = o.xAxis === "fold" ? "fold enrichment" : "−log<tspan baseline-shift=\"sub\" font-size=\"" + f * 0.7 + "\">10</tspan>(adjusted p)";
  P.push(`<text x="${(m.l + W - m.r) / 2}" y="${H - f * 0.5}" text-anchor="middle" ${font(f)}>${xl}</text>`);
  if (o.title) P.push(`<text x="${m.l - labelW}" y="${f * 1.2}" ${font(f)} font-weight="bold">${esc(o.title)}</text>`);
  // legends: size and colour
  const lx = W - m.r + f * 0.8; let ly = m.t + f;
  P.push(`<text x="${lx}" y="${ly}" ${font(f * 0.85)}>genes</text>`); ly += f * 0.6;
  for (const n of [...new Set([Math.max(1, Math.round(maxOv / 4)), Math.round(maxOv / 2), maxOv])]) { ly += rMax * 1.1 + f * 0.3; P.push(`<circle cx="${lx + rMax}" cy="${ly}" r="${rad(n).toFixed(2)}" fill="#999" stroke="#333" stroke-width="0.3"/><text x="${lx + rMax * 2 + f * 0.3}" y="${ly + f * 0.3}" ${font(f * 0.8)}>${n}</text>`); }
  ly += f * 1.6; P.push(`<text x="${lx}" y="${ly}" ${font(f * 0.85)}>adj. p</text>`); ly += f * 0.4;
  const gh = f * 4; P.push(`<defs><linearGradient id="g" x1="0" x2="0" y1="1" y2="0"><stop offset="0" stop-color="${col(lpMin)}"/><stop offset="1" stop-color="${col(lpMax)}"/></linearGradient></defs><rect x="${lx}" y="${ly}" width="${f * 0.8}" height="${gh}" fill="url(#g)" stroke="#333" stroke-width="0.3"/>`);
  P.push(`<text x="${lx + f * 1.1}" y="${ly + f * 0.5}" ${font(f * 0.75)}>${fmtP(Math.pow(10, -lpMax))}</text><text x="${lx + f * 1.1}" y="${ly + gh}" ${font(f * 0.75)}>${fmtP(Math.pow(10, -lpMin))}</text>`);
  P.push("</svg>");
  return P.join("");
}

export function enrichCaption(terms: Term[], queryName: string, o: DotOpts, eo: EnrichOpts, version: string, nQuery: number): string {
  const shown = Math.min(o.topN, terms.length);
  const method = eo.method === "g_SCS" ? "the g:SCS method" : eo.method === "fdr" ? "the Benjamini–Hochberg procedure" : "Bonferroni correction";
  return `Functional enrichment of the ${nQuery.toLocaleString("en-US")} ${queryName} genes, computed with g:Profiler (g:GOSt, data version ${version}) against ${eo.sources.map((s) => sourceNames[s] ?? s).join(", ")}${eo.background ? ", using all genes tested for differential expression as the statistical background" : ""}. The ${shown} terms with the smallest adjusted p-values (${method}, threshold ${eo.threshold}) are shown; the x axis gives the ${o.xAxis === "fold" ? "fold enrichment (observed over expected overlap)" : "−log10 adjusted p-value"}, dot size the number of query genes in the term, and colour the adjusted p-value.`;
}

function niceTicks(a: number, b: number, n: number): number[] { const span = b - a || 1, step0 = Math.pow(10, Math.floor(Math.log10(span / n))), err = span / n / step0; const step = step0 * (err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1); const out: number[] = []; for (let v = Math.ceil(a / step) * step; v <= b + 1e-9; v += step) out.push(+v.toFixed(10)); return out; }
const fmtP = (v: number) => (v >= 0.01 ? v.toFixed(2) : v.toExponential(0));
