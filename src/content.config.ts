import { defineCollection, z } from "astro:content";
import { glob } from "astro/loaders";

const lesson = z.object({
  title: z.string(),
  description: z.string(),
  order: z.number(),
  minutes: z.number().optional(),
});

// Entry ids look like "en/01-start-here"; the slug is the part after the language.
export const collections = {
  learn: defineCollection({ loader: glob({ pattern: "**/*.mdx", base: "./src/content/learn" }), schema: lesson }),
  read: defineCollection({ loader: glob({ pattern: "**/*.mdx", base: "./src/content/read" }), schema: lesson }),
};
