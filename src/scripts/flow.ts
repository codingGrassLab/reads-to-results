// Step state for the workflow bar. Each tool page calls markDone(step) when it has produced its
// output, and passes the page element that carries the output forward (an <a> whose click handler
// writes the hand-off to sessionStorage) so the bar's "Next" button does the same thing.
const KEY = "rr-done";
const done = (): number[] => { try { return JSON.parse(sessionStorage.getItem(KEY) || "[]"); } catch { return []; } };

function paint() {
  const bar = document.querySelector<HTMLElement>("[data-wf]"); if (!bar) return;
  const d = new Set(done()), cur = +bar.dataset.wf!;
  bar.querySelectorAll<HTMLElement>("[data-wf-step]").forEach((li) => li.classList.toggle("done", d.has(+li.dataset.wfStep!)));
  bar.querySelector("[data-wf-next]")?.classList.toggle("ready", d.has(cur));
}

let forward: HTMLAnchorElement | null = null;
/** Mark a step finished (its output exists) and, optionally, name the element that hands it on. */
export function markDone(step: number, handoff?: HTMLAnchorElement | null) {
  const d = new Set(done()); d.add(step);
  try { sessionStorage.setItem(KEY, JSON.stringify([...d])); } catch {}
  if (handoff) forward = handoff;
  paint();
}
/** The step's output is gone (a new input was loaded): later steps no longer follow from it. */
export function resetFrom(step: number) {
  try { sessionStorage.setItem(KEY, JSON.stringify(done().filter((s) => s < step))); } catch {}
  paint();
}
export function setHandoff(a: HTMLAnchorElement | null) { forward = a; }

export function initFlow() {
  paint();
  const next = document.querySelector<HTMLAnchorElement>("[data-wf-next]");
  next?.addEventListener("click", (e) => {
    if (!forward) return;               // nothing to carry: plain link to the next tool
    e.preventDefault();
    forward.click();                     // runs the page's hand-off and follows its href
  });
}
