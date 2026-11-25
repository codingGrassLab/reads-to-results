import { getCollection, type CollectionEntry } from "astro:content";
import type { Lang } from "../i18n/ui";

export type Coll = "learn" | "read";
export type Entry = CollectionEntry<Coll>;

export function split(id: string): { lang: Lang; slug: string } {
  const [lang, ...rest] = id.split("/");
  return { lang: lang as Lang, slug: rest.join("/") };
}

/** Lessons of one collection in one language, in order. */
export async function lessons(coll: Coll, lang: Lang): Promise<Entry[]> {
  const all = await getCollection(coll);
  return all.filter((e) => split(e.id).lang === lang).sort((a, b) => a.data.order - b.data.order);
}

export async function allPaths(coll: Coll) {
  const all = await getCollection(coll);
  return all.map((e) => {
    const { lang, slug } = split(e.id);
    return { params: { lang, slug }, props: { entry: e } };
  });
}
