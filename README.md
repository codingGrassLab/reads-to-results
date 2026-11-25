# Reads to Results

**Live site: https://codinggrasslab.github.io/reads-to-results/**

A free course in bulk RNA-seq analysis. One real experiment (airway) is followed from raw reads to differentially expressed genes and pathways with Salmon and DESeq2, every command shown with its actual output. Static site built with Astro.

- **Course** (11 lessons): raw reads → QC → trimming → how Salmon works → salmon quant → tximport → DESeq2 → plots → enrichment → your own project.
- **Break it yourself**: switches (design, lost sample, two donors, swapped label) over eight real DESeq2 runs on airway.
- **Ten mistakes**: confounding, too few replicates, swapped labels, normalized input, p vs padj, post-hoc fold-change cuts, unshrunken fold changes, enrichment background, Excel gene names, each made on purpose on airway.
- **Reading results** (8 pages): PCA, volcano, MA plot, p-values, log2 fold change, normalization, dispersion, heatmaps, each with an interactive figure drawn from the real data.

Dataset: `airway` (Bioconductor), from Himes et al. 2014, PLoS ONE 9(6): e99625; GEO GSE52778.

## Run it

```
npm install
npm run dev        # http://localhost:4321/reads-to-results/
npm test           # consistency checks on the data files
npm run build      # static site in dist/
```

Every push to `main` is tested, built and published to GitHub Pages by `.github/workflows/deploy.yml`.

## Regenerate the data

Only needed to recompute the numbers; editing text and pages needs no R. Needs R ≥ 4.4 with `airway`, `DESeq2`, `apeglm`, `clusterProfiler`, `org.Hs.eg.db`, `GenomicAlignments`, `rtracklayer`, `pheatmap`, `EnhancedVolcano`, `jsonlite`.

```
Rscript analysis/run.R       # the course analysis
Rscript analysis/salmon.R    # Salmon statistics and tximport
Rscript analysis/break.R     # Break it yourself and Ten mistakes
npm test
```

These rewrite `src/data/*.json` (read at build time), `public/data/*` (downloads and the tables the browser fetches) and `public/plots/*.png`. `npm test` checks that the files agree with each other (counts of significant genes, size factors, PCA separation, the paper's marker genes, and that every Break-it scenario adds up to the course result).

## Layout

```
analysis/run.R            the course's DESeq2 analysis
analysis/salmon.R         statistics from the Salmon runs that ship with airway, plus tximport
analysis/break.R          the eight Break-it runs and the ten mistakes
analysis/tx2gene.csv      GENCODE v29 transcript → gene table (built from the GTF, see salmon.R)
src/content/learn/en      course lessons (MDX; numbers are imported from src/data)
src/content/read/en       reading-results pages
src/pages                 home, Break it, Mistakes, About and the page routes
src/components            interactive plots (plain SVG/canvas, no dependencies)
src/i18n/ui.ts            interface text and link helpers
public/data               downloads, volcano.json, dispersion.json
```

## Licence

- Text and figures: [CC BY 4.0](LICENSE-text.md). Reuse freely with credit to Amarta Chowdhury.
- Code: [MIT](LICENSE).
- The `airway` data are not covered by these licences; they come from the Bioconductor package `airway` (Himes et al. 2014, LGPL).
