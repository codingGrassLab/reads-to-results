import { defineConfig } from "astro/config";
import sitemap from "@astrojs/sitemap";
import mdx from "@astrojs/mdx";

// GitHub Pages: https://codinggrasslab.github.io/reads-to-results/
export default defineConfig({
  site: "https://codinggrasslab.github.io",
  base: "/reads-to-results",
  trailingSlash: "always",
  integrations: [
    mdx(),
    sitemap({
      i18n: { defaultLocale: "en", locales: { en: "en" } },
      filter: (page) => page !== "https://codinggrasslab.github.io/reads-to-results/",
    }),
  ],
  markdown: {
    shikiConfig: { theme: "github-light", wrap: true },
  },
});
