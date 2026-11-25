# Reads to Results

A free course in bulk RNA-seq analysis. One real experiment (airway) is followed from raw reads to differentially expressed genes and pathways with Salmon and DESeq2: every command with its actual output, and every number computed by one R script. Static site built with Astro.

- **Course** (11 lessons): raw reads → QC → trimming → how Salmon works → salmon quant → tximport → DESeq2 → plots → enrichment → your own project.
- **Break it yourself**: switches (design, lost sample, two donors, swapped label) over eight real DESeq2 runs on airway.
- **Ten mistakes**: confounding, too few replicates, swapped labels, normalized input, p vs padj, post-hoc fold-change cuts, unshrunken fold changes, enrichment background, Excel gene names, each made on purpose on airway. Numbers from `analysis/break.R` → `src/data/break.json`.
- **Reading results** (8 pages): PCA, volcano, MA plot, p-values, log2 fold change, normalization, dispersion, heatmaps, each with an interactive figure drawn from the real data.
- **Data and methods**: inputs, versions, downloads, a check against the original paper, and a Colab notebook.

(The analysis tools that were built for a while are kept on the `analysis-tools` branch.)

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
