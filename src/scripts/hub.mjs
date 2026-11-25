// Hub-gene ranking on an undirected protein-protein interaction network, with the
// twelve node scores of cytoHubba (Chin et al. 2014). Plain JavaScript: used by the
// hub page in the browser and by scripts/test-hub.mjs in Node against networkx.

/** Parse an edge table: STRING's TSV export (#node1, node2, …, combined_score), or any
 *  two-column list with an optional numeric score column. Returns undirected, de-duplicated edges. */
export function parseNetwork(text, opts = {}) {
  const minScore = opts.minScore ?? 0;
  const lines = text.split(/\r?\n/).filter((l) => l.trim().length);
  if (!lines.length) throw new Error("The file is empty.");
  const delim = (lines[0].match(/\t/g) || []).length >= 1 ? "\t" : (lines[0].match(/;/g) || []).length >= 1 ? ";" : ",";
  const split = (l) => l.split(delim).map((x) => x.trim().replace(/^"|"$/g, ""));
  let header = split(lines[0]);
  let start = 1, iA = 0, iB = 1, iS = -1, format = "two columns";
  const low = header.map((h) => h.toLowerCase().replace(/^#/, ""));
  const fa = low.findIndex((h) => ["node1", "preferredname_a", "protein1", "source", "from", "gene1", "a"].includes(h));
  const fb = low.findIndex((h) => ["node2", "preferredname_b", "protein2", "target", "to", "gene2", "b"].includes(h));
  const fs = low.findIndex((h) => ["combined_score", "score", "weight", "confidence"].includes(h));
  if (fa >= 0 && fb >= 0) { iA = fa; iB = fb; iS = fs; format = low.includes("combined_score") || low.includes("preferredname_a") ? "STRING export" : "edge list with header"; }
  else if (!isNaN(Number(header[0])) || header.length < 2 || header.every((h) => !/[a-z]/i.test(h))) { start = 1; }
  else { start = 0; header = []; }
  if (fa < 0 && header.length >= 3 && header.every((h) => /[A-Za-z]/.test(h))) { start = 1; }
  const seen = new Set(), edges = [], nodes = new Set();
  let dropped = 0;
  for (let i = start; i < lines.length; i++) {
    const f = split(lines[i]); const a = f[iA], b = f[iB];
    if (!a || !b || a === b) continue;
    const score = iS >= 0 ? Number(f[iS]) : NaN;
    if (iS >= 0 && Number.isFinite(score)) { const sc = score > 1 ? score / 1000 : score; if (sc < minScore) { dropped++; continue; } }
    const key = a < b ? `${a}\u0001${b}` : `${b}\u0001${a}`;
    if (seen.has(key)) continue;
    seen.add(key); edges.push([a, b, Number.isFinite(score) ? (score > 1 ? score / 1000 : score) : null]); nodes.add(a); nodes.add(b);
  }
  if (!edges.length) throw new Error("No edges could be read. Expected two name columns (e.g. STRING's #node1 and node2), one interaction per line.");
  return { nodes: [...nodes].sort(), edges, format, hasScore: iS >= 0, dropped };
}

export function buildGraph(nodes, edges) {
  const idx = new Map(nodes.map((n, i) => [n, i]));
  const adj = nodes.map(() => new Set());
  for (const [a, b] of edges) { adj[idx.get(a)].add(idx.get(b)); adj[idx.get(b)].add(idx.get(a)); }
  return { nodes, idx, adj: adj.map((s) => [...s]) };
}

/** BFS distances from s; returns dist (−1 unreachable), sigma (shortest-path counts), order, predecessors. */
function bfs(g, s) {
  const n = g.nodes.length, dist = new Array(n).fill(-1), sigma = new Array(n).fill(0), pred = g.nodes.map(() => []), order = [];
  dist[s] = 0; sigma[s] = 1; const q = [s];
  for (let h = 0; h < q.length; h++) { const v = q[h]; order.push(v); for (const w of g.adj[v]) { if (dist[w] < 0) { dist[w] = dist[v] + 1; q.push(w); } if (dist[w] === dist[v] + 1) { sigma[w] += sigma[v]; pred[w].push(v); } } }
  return { dist, sigma, pred, order };
}

export const algorithms = ["Degree", "MCC", "DMNC", "MNC", "EPC", "BottleNeck", "EcCentricity", "Closeness", "Radiality", "Betweenness", "Stress", "ClusteringCoefficient"];

/** All twelve scores for every node. EPC is stochastic; `seed` makes it reproducible. */
export function scoreAll(g, opts = {}) {
  const n = g.nodes.length, adj = g.adj;
  const deg = adj.map((a) => a.length);
  const betw = new Array(n).fill(0), stress = new Array(n).fill(0), close = new Array(n).fill(0), rad = new Array(n).fill(0), ecc = new Array(n).fill(0), bottleneck = new Array(n).fill(0);
  const bfsAll = [];
  let diam = 0;
  for (let s = 0; s < n; s++) { const r = bfs(g, s); bfsAll.push(r); for (const d of r.dist) if (d > diam) diam = d; }
  for (let s = 0; s < n; s++) {
    const { dist, sigma, pred, order } = bfsAll[s];
    // Brandes accumulation for betweenness and stress
    const delta = new Array(n).fill(0), sdelta = new Array(n).fill(0);
    for (let i = order.length - 1; i >= 0; i--) { const w = order[i]; for (const v of pred[w]) { delta[v] += (sigma[v] / sigma[w]) * (1 + delta[w]); sdelta[v] += sigma[v] * (1 + sdelta[w] / sigma[w]); } if (w !== s) { betw[w] += delta[w]; stress[w] += sdelta[w]; } }
    // closeness (sum of reciprocal distances, as cytoHubba), radiality, eccentricity
    let sum = 0, radSum = 0, maxd = 0;
    for (let t = 0; t < n; t++) { if (t === s || dist[t] < 0) continue; sum += 1 / dist[t]; radSum += diam + 1 - dist[t]; if (dist[t] > maxd) maxd = dist[t]; }
    close[s] = sum; rad[s] = n > 1 ? radSum / (n - 1) : 0; ecc[s] = maxd ? 1 / maxd : 0;
    // BottleNeck: in a shortest-path tree rooted at s, nodes whose subtree holds more than n/4 nodes
    const parent = new Array(n).fill(-1); for (const w of order) if (w !== s) parent[w] = pred[w][0];
    const sub = new Array(n).fill(1);
    for (let i = order.length - 1; i > 0; i--) { const w = order[i]; sub[parent[w]] += sub[w]; }
    for (const w of order) if (w !== s && sub[w] > n / 4) bottleneck[w] += 1;
  }
  for (let v = 0; v < n; v++) { betw[v] /= 2; stress[v] /= 2; }
  // clustering coefficient
  const cc = adj.map((nb, v) => { const k = nb.length; if (k < 2) return 0; let t = 0; const set = new Set(nb); for (const a of nb) for (const b of adj[a]) if (b > a && set.has(b)) t++; return (2 * t) / (k * (k - 1)); });
  // MNC and DMNC: largest connected component of the neighbourhood subgraph
  const mnc = new Array(n).fill(0), dmnc = new Array(n).fill(0);
  for (let v = 0; v < n; v++) {
    const nb = adj[v]; if (!nb.length) continue; const set = new Set(nb); const seen = new Set(); let best = 0, bestE = 0;
    for (const s0 of nb) { if (seen.has(s0)) continue; const comp = [s0]; seen.add(s0); for (let h = 0; h < comp.length; h++) for (const w of adj[comp[h]]) if (set.has(w) && !seen.has(w)) { seen.add(w); comp.push(w); } let e = 0; const cs = new Set(comp); for (const a of comp) for (const b of adj[a]) if (b > a && cs.has(b)) e++; if (comp.length > best) { best = comp.length; bestE = e; } }
    mnc[v] = best; dmnc[v] = best > 1 ? bestE / Math.pow(best, 1.7) : 0;
  }
  // MCC: sum over maximal cliques containing v of (|C|-1)!; a node whose cliques are all edges scores its degree
  const mcc = new Array(n).fill(0);
  const fact = (k) => { let f = 1; for (let i = 2; i <= k; i++) f *= i; return f; };
  const adjSets = adj.map((a) => new Set(a));
  let cliqueCount = 0;
  const report = (R) => { cliqueCount++; const w = fact(R.length - 1); for (const v of R) mcc[v] += w; };
  const bk = (R, P, X) => { if (cliqueCount > 200000) return; if (!P.length && !X.length) { if (R.length >= 2) report(R); return; } let u = -1, best = -1; for (const c of [...P, ...X]) { let cnt = 0; for (const p of P) if (adjSets[c].has(p)) cnt++; if (cnt > best) { best = cnt; u = c; } } for (const v of P.filter((p) => !adjSets[u].has(p))) { bk([...R, v], P.filter((p) => adjSets[v].has(p)), X.filter((x) => adjSets[v].has(x))); P = P.filter((p) => p !== v); X.push(v); } };
  bk([], [...Array(n).keys()], []);
  // EPC: edge percolated component, mean size of the component containing v over `iters` random edge removals (p = 0.5)
  const iters = opts.epcIters ?? 1000; let seed = opts.seed ?? 42;
  const rnd = () => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296; };
  const edgeList = []; for (let v = 0; v < n; v++) for (const w of adj[v]) if (w > v) edgeList.push([v, w]);
  const epc = new Array(n).fill(0);
  for (let it = 0; it < iters; it++) {
    const keep = adj.map(() => []);
    for (const [a, b] of edgeList) if (rnd() < 0.5) { keep[a].push(b); keep[b].push(a); }
    const comp = new Array(n).fill(-1); let c = 0; const sizes = [];
    for (let v = 0; v < n; v++) { if (comp[v] >= 0) continue; const st = [v]; comp[v] = c; let sz = 0; while (st.length) { const x = st.pop(); sz++; for (const y of keep[x]) if (comp[y] < 0) { comp[y] = c; st.push(y); } } sizes.push(sz); c++; }
    for (let v = 0; v < n; v++) epc[v] += sizes[comp[v]];
  }
  for (let v = 0; v < n; v++) epc[v] /= iters;
  return { Degree: deg, MCC: mcc, DMNC: dmnc, MNC: mnc, EPC: epc, BottleNeck: bottleneck, EcCentricity: ecc, Closeness: close, Radiality: rad, Betweenness: betw, Stress: stress, ClusteringCoefficient: cc, diameter: diam };
}

/** Rank nodes by one score (higher = more central), ties by degree then name. */
export function rank(g, scores, algo) {
  const s = scores[algo];
  return [...Array(g.nodes.length).keys()].sort((a, b) => (s[b] - s[a]) || (scores.Degree[b] - scores.Degree[a]) || g.nodes[a].localeCompare(g.nodes[b]));
}
