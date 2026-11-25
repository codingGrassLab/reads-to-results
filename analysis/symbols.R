# Ensembl gene ID -> symbol tables for the results tool (public/data/symbols_<species>.json).
#   Rscript analysis/symbols.R
suppressPackageStartupMessages({ library(jsonlite); library(AnnotationDbi) })
make <- function(pkg, out) {
  if (!requireNamespace(pkg, quietly = TRUE)) { cat("skip", pkg, "\n"); return(invisible()) }
  db <- get(pkg, asNamespace(pkg))
  ens <- keys(db, keytype = "ENSEMBL")
  sym <- mapIds(db, ens, "SYMBOL", "ENSEMBL", multiVals = "first")
  sym <- sym[!is.na(sym)]
  write_json(as.list(sym), out, auto_unbox = TRUE)
  cat(out, length(sym), "ids;", pkg, as.character(packageVersion(pkg)), "\n")
}
make("org.Hs.eg.db", "public/data/symbols_human.json")
make("org.Mm.eg.db", "public/data/symbols_mouse.json")
