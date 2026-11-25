# Reads to Results

Analysis tools for RNA-seq differential-expression results (publication-ready figures, tables and methods text, all in the browser), with a complete worked analysis of one real dataset behind them. Static site built with Astro; every number and plot comes from one R script.

- **Run DESeq2 in your browser**: real DESeq2 via webR (R → WebAssembly; packages from the Bioconductor r-universe wasm builds, locfit namespace stub hosted in `public/wasm/`). Count matrix + sample sheet → results table, handed to the results tool.
- **Explore your results**: upload a DESeq2/edgeR/limma table → thresholds → significant genes (CSV/XLSX), volcano and MA plots (SVG, 300-dpi PNG at journal widths), figure legends.
- **Methods writer**: tools, versions and settings → methods paragraph with verified citations and BibTeX.
- **Course** (11 lessons): raw reads → QC → trimming → how Salmon works → salmon quant → tximport → DESeq2 → plots → enrichment → your own project.
- **Reading results** (8 pages): PCA, volcano, MA plot, p-values, log2 fold change, normalization, dispersion, heatmaps, each with an interactive figure drawn from the real data.
- **Data and methods**: inputs, versions, downloads, a check against the original paper, and a Colab notebook.

Dataset: `airway` (Bioconductor), from Himes et al. 2014, PLoS ONE 9(6): e99625; GEO GSE52778.

## Run it

```
npm install
npm run dev        # http://localhost:4321
npm test           # consistency checks on the data files
npm run build      # static site in dist/
```

## Regenerate the data

Needs R ≥ 4.4 with `airway`, `DESeq2`, `apeglm`, `clusterProfiler`, `org.Hs.eg.db`, `GenomicAlignments`, `rtracklayer`, `pheatmap`, `EnhancedVolcano`, `jsonlite`.

```
Rscript analysis/run.R
npm test
```

This rewrites `src/data/*.json` (read at build time), `public/data/*` (downloads and the two tables the browser fetches) and `public/plots/*.png`. The `test` script checks that the files agree with each other (counts of significant genes, size factors, PCA separation, the paper's marker genes).

## Layout

```
analysis/run.R            the whole DESeq2 analysis
analysis/salmon.R         statistics from the Salmon runs that ship with airway, plus tximport
analysis/tx2gene.csv      GENCODE v29 transcript → gene table (built from the GTF, see salmon.R)
src/content/learn/en     course lessons (MDX; numbers are imported from src/data)
src/content/read/en      reading-results pages
src/components            interactive plots (plain SVG/canvas, no dependencies)
src/i18n/ui.ts            interface text
public/data               downloads, volcano.json, dispersion.json
public/notebooks          Colab notebook (R runtime)
```

## Licence

Text CC BY 4.0, code MIT. Data from the airway package, authors Himes et al.
