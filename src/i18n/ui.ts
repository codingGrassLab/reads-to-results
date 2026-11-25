export const languages = { en: "English" } as const;
export type Lang = keyof typeof languages;
export const langs = Object.keys(languages) as Lang[];
export const defaultLang: Lang = "en";

/** Text that exists in both languages. */
export type T = Record<Lang, string>;

const ui = {
  en: {
    "site.name": "Reads to Results",
    "site.tagline": "Publication-ready figures, tables and methods text from your RNA-seq results",
    "site.description":
      "Upload a DESeq2, edgeR or limma results table and download significant genes, journal-size volcano and MA plots, figure legends and methods text, all in your browser. With a complete worked analysis of a real dataset for anyone who wants to understand the steps.",
    "nav.home": "Home",
    "nav.dataset": "Find dataset",
    "nav.results": "Results",
    "nav.enrich": "Enrichment",
    "nav.hub": "Hub genes",
    "nav.learn": "Course",
    "nav.read": "Reading plots",
    "nav.methods": "Methods",
    "nav.data": "Data",
    "nav.about": "About",
    "nav.menu": "Menu",
    "nav.language": "Language",
    "common.skip": "Skip to content",
    "common.next": "Next",
    "common.prev": "Previous",
    "common.part": "Part",
    "common.lesson": "Lesson",
    "common.minutes": "min read",
    "common.toc": "In this lesson",
    "common.download": "Download",
    "common.source": "Source",
    "common.figure": "Figure",
    "common.hover": "Hover or tap a point for the gene name.",
    "common.loading": "Loading data…",
    "common.copy": "Copy",
    "common.copied": "Copied",
    "home.eyebrow": "RNA-seq · DESeq2 · real data",
    "home.title": "From raw reads to a result you can defend",
    "home.lead":
      "One real experiment, followed from FASTQ files to a list of genes and pathways. Each step shows the command, the actual output, and why it is done. Nothing here is made up: every number was computed by the code on the data page.",
    "home.ctaStart": "Start the course",
    "home.ctaRead": "How to read the plots",
    "home.learnTitle": "The course",
    "home.learnLead": "Eleven lessons that follow one dataset all the way through, with Salmon and DESeq2.",
    "home.readTitle": "Reading results",
    "home.readLead": "The plots every RNA-seq paper shows, each one interactive, drawn from the same dataset.",
    "home.dataTitle": "The dataset",
    "home.dataLead":
      "Human airway smooth muscle cells from four donors, treated with the steroid dexamethasone or left untreated (Himes et al. 2014). Eight samples, a paired design, and a strong, well-understood signal.",
    "home.honest": "What this site does not do",
    "home.honestText":
      "It does not run analyses in your browser, and it is not a replacement for the DESeq2 documentation. It shows one analysis done carefully, so you can do yours.",
    "learn.title": "The course",
    "learn.lead": "Follow one dataset from raw reads to pathways. Read in order the first time.",
    "read.title": "Reading results",
    "read.lead": "Each page takes one kind of plot or number, shows it on the airway data, and explains what it can and cannot tell you.",
    "data.title": "Data and methods",
    "data.lead": "Everything on this site was produced by one R script. Here are the inputs, the outputs, the versions and the checks.",
    "about.title": "About",
    "footer.built": "Static site. No tracking, no cookies, no server of our own.",
    "footer.license": "Text under CC BY 4.0; code under MIT. Data from the airway Bioconductor package (Himes et al. 2014).",
    "404.title": "Page not found",
    "404.text": "That page does not exist. Start from the home page.",
  },
} as const satisfies Record<Lang, Record<string, string>>;

export type UIKey = keyof (typeof ui)["en"];

export function t(lang: Lang, key: UIKey): string {
  return ui[lang][key] ?? ui.en[key];
}

/** `/en/`, `/en/learn/`, … (always with a trailing slash). */
export function path(lang: Lang, ...parts: string[]): string {
  const p = parts.filter(Boolean).join("/");
  return `/${lang}/${p}${p ? "/" : ""}`;
}

/** The same page in another language. */
export function switchLang(pathname: string, to: Lang): string {
  const segs = pathname.split("/").filter(Boolean);
  if (segs.length && (langs as string[]).includes(segs[0])) segs[0] = to;
  else segs.unshift(to);
  return "/" + segs.join("/") + "/";
}

/** Numbers shown in running text. */
export function num(_lang: Lang, n: number | string): string {
  return typeof n === "number" ? n.toLocaleString("en-US") : n;
}
