// Checks the hub-gene scores against networkx on random graphs. Run: node scripts/test-hub.mjs
import { buildGraph, scoreAll } from "../src/scripts/hub.mjs";
import { execFileSync } from "node:child_process";
import { writeFileSync, mkdtempSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";

let seed = 7; const rnd = () => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296; };
let failed = 0;
const ok = (c, m) => { if (!c) { failed++; console.error("FAIL", m); } else console.log("ok  ", m); };
const close = (a, b, tol = 1e-6) => Math.abs(a - b) <= tol * Math.max(1, Math.abs(b));

for (const [n, p] of [[12, 0.3], [40, 0.12], [80, 0.06]]) {
  const nodes = [...Array(n).keys()].map((i) => `g${i}`), edges = [];
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) if (rnd() < p) edges.push([nodes[i], nodes[j], null]);
  const g = buildGraph(nodes, edges), s = scoreAll(g, { epcIters: 200 });
  const dir = mkdtempSync(join(tmpdir(), "hub-"));
  writeFileSync(join(dir, "edges.txt"), edges.map((e) => `${e[0]} ${e[1]}`).join("\n") + "\n");
  const py = `
import networkx as nx, json, sys, math
G = nx.read_edgelist(sys.argv[1]); G.add_nodes_from(${JSON.stringify(nodes)})
diam = max(max(d.values()) for _, d in nx.all_pairs_shortest_path_length(G))
out = {}
bet = nx.betweenness_centrality(G, normalized=False)
clu = nx.clustering(G)
for v in G.nodes:
    d = nx.single_source_shortest_path_length(G, v)
    out[v] = dict(Degree=G.degree(v), Betweenness=bet[v], ClusteringCoefficient=clu[v],
                  Closeness=sum(1/x for u,x in d.items() if u!=v), Radiality=sum(diam+1-x for u,x in d.items() if u!=v)/(G.number_of_nodes()-1),
                  EcCentricity=(1/max(d.values())) if max(d.values())>0 else 0)
    # MCC: sum over maximal cliques containing v of (|C|-1)!
    out[v]["MCC"] = sum(math.factorial(len(c)-1) for c in nx.find_cliques(G) if v in c and len(c)>=2)
    # MNC / DMNC
    nb = list(G.neighbors(v)); H = G.subgraph(nb)
    comps = list(nx.connected_components(H)) if nb else []
    big = max(comps, key=len) if comps else set()
    out[v]["MNC"] = len(big); out[v]["DMNC"] = (H.subgraph(big).number_of_edges()/len(big)**1.7) if len(big)>1 else 0
print(json.dumps(out))`;
  writeFileSync(join(dir, "ref.py"), py);
  const ref = JSON.parse(execFileSync("python", ["-I", join(dir, "ref.py"), join(dir, "edges.txt")], { encoding: "utf8" }));
  for (const algo of ["Degree", "Betweenness", "ClusteringCoefficient", "Closeness", "Radiality", "EcCentricity", "MCC", "MNC", "DMNC"]) {
    const bad = nodes.filter((v, i) => !close(s[algo][i], ref[v][algo]));
    ok(!bad.length, `n=${n}: ${algo} matches networkx${bad.length ? ` (${bad.length} differ, e.g. ${bad[0]}: ${s[algo][g.idx.get(bad[0])]} vs ${ref[bad[0]][algo]})` : ""}`);
  }
  // Stress: brute-force count of shortest paths through each node (small graphs)
  if (n <= 40) {
    const py2 = `
import networkx as nx, json, sys
G = nx.read_edgelist(sys.argv[1]); G.add_nodes_from(${JSON.stringify(nodes)})
st = {v:0 for v in G.nodes}
nodes = list(G.nodes)
for i,a in enumerate(nodes):
    for b in nodes[i+1:]:
        if nx.has_path(G,a,b):
            for path in nx.all_shortest_paths(G,a,b):
                for v in path[1:-1]: st[v]+=1
print(json.dumps(st))`;
    writeFileSync(join(dir, "ref2.py"), py2);
    const ref2 = JSON.parse(execFileSync("python", ["-I", join(dir, "ref2.py"), join(dir, "edges.txt")], { encoding: "utf8" }));
    const bad = nodes.filter((v, i) => !close(s.Stress[i], ref2[v]));
    ok(!bad.length, `n=${n}: Stress matches brute-force shortest-path counting`);
  }
  ok(s.EPC.every((x) => x >= 1 && x <= n), `n=${n}: EPC within [1, n]`);
  ok(s.BottleNeck.every((x) => Number.isInteger(x) && x >= 0 && x <= n), `n=${n}: BottleNeck is a count in [0, n]`);
}
console.log(failed ? `\n${failed} check(s) failed` : "\nall hub checks passed");
process.exit(failed ? 1 : 0);
