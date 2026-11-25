import { defineConfig } from "astro/config";
import sitemap from "@astrojs/sitemap";
import mdx from "@astrojs/mdx";

// `site` is a placeholder until the GitHub Pages address (or a domain) is decided.
export default defineConfig({
  site: "https://example.github.io",
  trailingSlash: "always",
  integrations: [
    mdx(),
    sitemap({
      i18n: { defaultLocale: "en", locales: { en: "en" } },
      filter: (page) => page !== "https://example.github.io/",
    }),
  ],
  markdown: {
    shikiConfig: { theme: "github-light", wrap: true },
  },
});
