// Set overlaps as publication figures: Venn (2 or 3 sets) and an UpSet-style plot (any number).
// Plain SVG strings, physical size in mm, Arial, no dependencies.

const MM = 72 / 25.4;
const esc = (s: string) => s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
const font = (pt: number) => `font-family="Arial, Helvetica, sans-serif" font-size="${+pt.toFixed(2)}"`;
export const setColours = ["#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"];

export interface SetInput { name: string; items: string[] }

/** Every region of the overlap: key is a bitmask of set membership, value the items. */
export function regions(sets: SetInput[]): Map<number, string[]> {
  const where = new Map<string, number>();
  sets.forEach((s, i) => s.items.forEach((x) => where.set(x, (where.get(x) ?? 0) | (1 << i))));
  const out = new Map<number, string[]>();
  for (const [x, m] of where) { if (!out.has(m)) out.set(m, []); out.get(m)!.push(x); }
  for (const v of out.values()) v.sort();
  return out;
}

export function vennSVG(sets: SetInput[], o: { widthMm: number; fontPt: number; title?: string }): string {
  const n = sets.length; if (n < 2 || n > 3) throw new Error("Venn needs 2 or 3 sets");
  const W = o.widthMm * MM, H = W * (n === 2 ? 0.62 : 0.92), f = o.fontPt;
  const reg = regions(sets);
  const count = (mask: number) => (reg.get(mask) ?? []).length;
  const r = n === 2 ? W * 0.26 : W * 0.24;
  const centres: [number, number][] = n === 2 ? [[W * 0.37, H * 0.52], [W * 0.63, H * 0.52]] : [[W * 0.37, H * 0.4], [W * 0.63, H * 0.4], [W * 0.5, H * 0.64]];
  const P: string[] = [`<svg xmlns="http://www.w3.org/2000/svg" width="${o.widthMm}mm" height="${(H / MM).toFixed(1)}mm" viewBox="0 0 ${W.toFixed(1)} ${H.toFixed(1)}"><rect width="100%" height="100%" fill="#fff"/>`];
  if (o.title) P.push(`<text x="${f}" y="${f * 1.3}" ${font(f)} font-weight="bold">${esc(o.title)}</text>`);
  centres.forEach(([cx, cy], i) => P.push(`<circle cx="${cx.toFixed(1)}" cy="${cy.toFixed(1)}" r="${r.toFixed(1)}" fill="${setColours[i]}" fill-opacity="0.35" stroke="${setColours[i]}" stroke-width="0.8"/>`));
  // region label positions
  const lab = (mask: number, x: number, y: number, big = false) => P.push(`<text x="${x.toFixed(1)}" y="${y.toFixed(1)}" text-anchor="middle" ${font(big ? f * 1.3 : f)} font-weight="${big ? "bold" : "normal"}">${count(mask)}</text>`);
  if (n === 2) {
    const [a, b] = centres; lab(1, a[0] - r * 0.45, a[1] + f * 0.4); lab(2, b[0] + r * 0.45, b[1] + f * 0.4); lab(3, (a[0] + b[0]) / 2, a[1] + f * 0.4, true);
  } else {
    const [a, b, c] = centres;
    lab(1, a[0] - r * 0.45, a[1] - r * 0.2); lab(2, b[0] + r * 0.45, b[1] - r * 0.2); lab(4, c[0], c[1] + r * 0.55);
    lab(3, (a[0] + b[0]) / 2, a[1] - r * 0.35); lab(5, (a[0] + c[0]) / 2 - r * 0.2, (a[1] + c[1]) / 2 + r * 0.2); lab(6, (b[0] + c[0]) / 2 + r * 0.2, (b[1] + c[1]) / 2 + r * 0.2);
    lab(7, (a[0] + b[0] + c[0]) / 3, (a[1] + b[1] + c[1]) / 3 + f * 0.4, true);
  }
  // set names with sizes
  const namePos: [number, number, string][] = n === 2 ? [[centres[0][0] - r * 0.3, centres[0][1] - r - f * 0.6, "middle"], [centres[1][0] + r * 0.3, centres[1][1] - r - f * 0.6, "middle"]] : [[centres[0][0] - r * 0.4, centres[0][1] - r - f * 0.6, "middle"], [centres[1][0] + r * 0.4, centres[1][1] - r - f * 0.6, "middle"], [centres[2][0], centres[2][1] + r + f * 1.4, "middle"]];
  sets.forEach((s, i) => P.push(`<text x="${namePos[i][0].toFixed(1)}" y="${namePos[i][1].toFixed(1)}" text-anchor="${namePos[i][2]}" ${font(f)} fill="${setColours[i]}" font-weight="bold">${esc(s.name)} (${s.items.length})</text>`));
  P.push("</svg>");
  return P.join("");
}

/** UpSet-style plot: bars for each intersection (exclusive regions), dots for set membership, set sizes at left. */
export function upsetSVG(sets: SetInput[], o: { widthMm: number; fontPt: number; title?: string; maxBars?: number }): string {
  const n = sets.length, f = o.fontPt, W = o.widthMm * MM;
  const reg = [...regions(sets).entries()].map(([mask, items]) => ({ mask, size: items.length, k: popcount(mask) })).sort((a, b) => b.k - a.k || b.size - a.size).slice(0, o.maxBars ?? 20);
  const leftW = f * 9, barH = f * 7, rowH = f * 1.4, top = (o.title ? f * 1.6 : 0) + f * 1.2, colW = Math.max(f * 1.4, (W - leftW - f) / Math.max(reg.length, 1));
  const H = top + barH + f * 1.2 + n * rowH + f * 1.2;
  const maxSize = Math.max(...reg.map((r) => r.size), 1);
  const P: string[] = [`<svg xmlns="http://www.w3.org/2000/svg" width="${o.widthMm}mm" height="${(H / MM).toFixed(1)}mm" viewBox="0 0 ${W.toFixed(1)} ${H.toFixed(1)}"><rect width="100%" height="100%" fill="#fff"/>`];
  if (o.title) P.push(`<text x="${f * 0.5}" y="${f * 1.2}" ${font(f)} font-weight="bold">${esc(o.title)}</text>`);
  const y0 = top + barH, x0 = leftW;
  P.push(`<line x1="${x0}" x2="${x0}" y1="${top}" y2="${y0}" stroke="#000" stroke-width="0.6"/>`);
  for (const t of [0, Math.round(maxSize / 2), maxSize]) { const y = y0 - (t / maxSize) * barH; P.push(`<line x1="${x0 - 2.5}" x2="${x0}" y1="${y}" y2="${y}" stroke="#000" stroke-width="0.6"/><text x="${x0 - 3.5}" y="${y + f * 0.35}" text-anchor="end" ${font(f * 0.85)}>${t}</text>`); }
  P.push(`<text transform="translate(${f * 0.9},${(top + y0) / 2}) rotate(-90)" text-anchor="middle" ${font(f * 0.9)}>genes in intersection</text>`);
  reg.forEach((r, j) => {
    const cx = x0 + colW * (j + 0.5), h = (r.size / maxSize) * barH;
    P.push(`<rect x="${cx - colW * 0.3}" y="${y0 - h}" width="${colW * 0.6}" height="${h}" fill="${r.k === n ? "#D55E00" : "#333"}"/><text x="${cx}" y="${y0 - h - f * 0.3}" text-anchor="middle" ${font(f * 0.85)}>${r.size}</text>`);
    const ys = sets.map((_, i) => y0 + f * 1.2 + i * rowH + rowH / 2).filter((_, i) => r.mask & (1 << i));
    if (ys.length > 1) P.push(`<line x1="${cx}" x2="${cx}" y1="${Math.min(...ys)}" y2="${Math.max(...ys)}" stroke="#333" stroke-width="1"/>`);
    sets.forEach((_, i) => P.push(`<circle cx="${cx}" cy="${y0 + f * 1.2 + i * rowH + rowH / 2}" r="${f * 0.42}" fill="${r.mask & (1 << i) ? "#333" : "#ddd"}"/>`));
  });
  sets.forEach((s, i) => P.push(`<text x="${x0 - f * 0.5}" y="${y0 + f * 1.2 + i * rowH + rowH / 2 + f * 0.35}" text-anchor="end" ${font(f * 0.9)}>${esc(s.name)} (${s.items.length})</text>`));
  P.push("</svg>");
  return P.join("");
}

function popcount(m: number) { let c = 0; while (m) { c += m & 1; m >>= 1; } return c; }
