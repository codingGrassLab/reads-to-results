// Small plotting helpers shared by the interactive figures. No dependencies.

export type Scale = ((v: number) => number) & { domain: [number, number]; range: [number, number]; invert: (p: number) => number };

export function linear(domain: [number, number], range: [number, number]): Scale {
  const [d0, d1] = domain, [r0, r1] = range;
  const k = (r1 - r0) / (d1 - d0 || 1);
  const f = ((v: number) => r0 + (v - d0) * k) as Scale;
  f.domain = domain; f.range = range;
  f.invert = (p: number) => d0 + (p - r0) / k;
  return f;
}

export function ticks(d0: number, d1: number, n = 6): number[] {
  const span = d1 - d0 || 1;
  const step0 = Math.pow(10, Math.floor(Math.log10(span / n)));
  const err = (span / n) / step0;
  const step = step0 * (err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1);
  const out: number[] = [];
  for (let v = Math.ceil(d0 / step) * step; v <= d1 + 1e-9; v += step) out.push(+v.toFixed(10));
  return out;
}

export function fmt(v: number, digits = 2): string {
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e4 || a < 1e-3) return v.toExponential(digits);
  return (+v.toFixed(digits)).toString();
}

export function el<K extends keyof SVGElementTagNameMap>(tag: K, attrs: Record<string, string | number> = {}, text?: string): SVGElementTagNameMap[K] {
  const e = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, String(v));
  if (text !== undefined) e.textContent = text;
  return e;
}

/** Draw x and y axes with ticks into `g`. */
export function axes(g: SVGGElement, x: Scale, y: Scale, opts: { xLabel?: string; yLabel?: string; xTicks?: number[]; yTicks?: number[]; xFmt?: (v: number) => string; yFmt?: (v: number) => string } = {}) {
  const ax = el("g", { class: "axis" });
  const [x0, x1] = x.range, [y0, y1] = y.range;
  ax.appendChild(el("path", { d: `M${x0},${y0}H${x1}`, fill: "none" }));
  ax.appendChild(el("path", { d: `M${x0},${y0}V${y1}`, fill: "none" }));
  for (const t of opts.xTicks ?? ticks(x.domain[0], x.domain[1])) {
    if (t < Math.min(...x.domain) - 1e-9 || t > Math.max(...x.domain) + 1e-9) continue;
    ax.appendChild(el("line", { x1: x(t), x2: x(t), y1: y0, y2: y0 + 5 }));
    ax.appendChild(el("text", { x: x(t), y: y0 + 18, "text-anchor": "middle" }, (opts.xFmt ?? fmt)(t)));
  }
  for (const t of opts.yTicks ?? ticks(y.domain[0], y.domain[1])) {
    if (t < Math.min(...y.domain) - 1e-9 || t > Math.max(...y.domain) + 1e-9) continue;
    ax.appendChild(el("line", { x1: x0 - 5, x2: x0, y1: y(t), y2: y(t) }));
    ax.appendChild(el("text", { x: x0 - 8, y: y(t) + 4, "text-anchor": "end" }, (opts.yFmt ?? fmt)(t)));
  }
  g.appendChild(ax);
  if (opts.xLabel) g.appendChild(el("text", { class: "label", x: (x0 + x1) / 2, y: y0 + 38, "text-anchor": "middle" }, opts.xLabel));
  if (opts.yLabel) g.appendChild(el("text", { class: "label", x: x0 - 44, y: (y0 + y1) / 2, "text-anchor": "middle", transform: `rotate(-90 ${x0 - 44} ${(y0 + y1) / 2})` }, opts.yLabel));
}

/** Tooltip helper bound to a `.plot` container. */
export function tooltip(container: HTMLElement) {
  const tip = document.createElement("div");
  tip.className = "tip";
  container.appendChild(tip);
  return {
    show(html: string, clientX: number, clientY: number) {
      const r = container.getBoundingClientRect();
      tip.innerHTML = html;
      tip.style.left = `${clientX - r.left}px`;
      tip.style.top = `${clientY - r.top}px`;
      tip.classList.add("show");
    },
    hide() { tip.classList.remove("show"); },
  };
}

export function readData<T>(root: HTMLElement): T {
  return JSON.parse(root.querySelector('script[type="application/json"]')!.textContent!) as T;
}

export const colors = { trt: "#0f6e6a", untrt: "#b86b1f", up: "#c0392b", down: "#2c5fa8", ns: "#b7bdbb", ink: "#15201f", muted: "#5b6764" };
export const cellShape: Record<string, string> = { N61311: "circle", N052611: "square", N080611: "triangle", N061011: "diamond" };

export function shape(kind: string, cx: number, cy: number, r: number, attrs: Record<string, string | number>) {
  switch (kind) {
    case "square": return el("rect", { x: cx - r, y: cy - r, width: 2 * r, height: 2 * r, ...attrs });
    case "triangle": return el("path", { d: `M${cx},${cy - r * 1.2}L${cx + r * 1.1},${cy + r * 0.8}L${cx - r * 1.1},${cy + r * 0.8}Z`, ...attrs });
    case "diamond": return el("path", { d: `M${cx},${cy - r * 1.25}L${cx + r * 1.25},${cy}L${cx},${cy + r * 1.25}L${cx - r * 1.25},${cy}Z`, ...attrs });
    default: return el("circle", { cx, cy, r, ...attrs });
  }
}

/** Setup a canvas sized to its container with device-pixel scaling. */
export function canvas(container: HTMLElement, width: number, height: number) {
  const c = document.createElement("canvas");
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  c.width = width * dpr; c.height = height * dpr;
  c.style.width = "100%"; c.style.aspectRatio = `${width} / ${height}`;
  container.appendChild(c);
  const ctx = c.getContext("2d")!;
  ctx.scale(dpr, dpr);
  return { c, ctx, dpr };
}

export function canvasAxes(ctx: CanvasRenderingContext2D, x: Scale, y: Scale, xLabel: string, yLabel: string, xFmt: (v: number) => string = fmt, yFmt: (v: number) => string = fmt) {
  const [x0, x1] = x.range, [y0, y1] = y.range;
  ctx.save();
  ctx.strokeStyle = "#c9cec8"; ctx.fillStyle = "#5b6764"; ctx.lineWidth = 1;
  ctx.font = "11px Inter, system-ui, sans-serif";
  ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y0); ctx.moveTo(x0, y0); ctx.lineTo(x0, y1); ctx.stroke();
  ctx.textAlign = "center"; ctx.textBaseline = "top";
  for (const t of ticks(x.domain[0], x.domain[1])) { if (t < x.domain[0] || t > x.domain[1]) continue; ctx.beginPath(); ctx.moveTo(x(t), y0); ctx.lineTo(x(t), y0 + 5); ctx.stroke(); ctx.fillText(xFmt(t), x(t), y0 + 8); }
  ctx.textAlign = "right"; ctx.textBaseline = "middle";
  for (const t of ticks(y.domain[0], y.domain[1])) { if (t < y.domain[0] || t > y.domain[1]) continue; ctx.beginPath(); ctx.moveTo(x0 - 5, y(t)); ctx.lineTo(x0, y(t)); ctx.stroke(); ctx.fillText(yFmt(t), x0 - 8, y(t)); }
  ctx.font = "12px Inter, system-ui, sans-serif"; ctx.textAlign = "center"; ctx.textBaseline = "alphabetic";
  ctx.fillText(xLabel, (x0 + x1) / 2, y0 + 38);
  ctx.save(); ctx.translate(x0 - 44, (y0 + y1) / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(yLabel, 0, 0); ctx.restore();
  ctx.restore();
}
