# kycg/tmp — dogfooding drop zone

This folder is the durable collection point for **internal use of kycg**. It is
git-tracked (not scratch), so anything left here survives for the maintainers to
review. The paper is published (Goldberg *et al.*, *Sci Adv* 2025) and the goal
is adoption — we want rough edges found here, by us, before external users hit
them.

## When to drop something here

While doing any real methylation functional-enrichment task, use `kycg test` /
`kycg annotate` (or the R `knowYourCG` for arrays) instead of a bespoke
hypergeometric/Fisher script. Then, if you either:

- **hit an issue** — a bug, a confusing interface, a missing feature, or a
  wrong/suspicious number; or
- **make a showcase** — a clean enrichment result, a nice volcano/dot/bar plot,
  a compelling biological use case,

leave it here.

## How to drop it

- Name it date-first and descriptive: `YYYYMMDD_<slug>.<ext>`
  (e.g. `20260724_chromhmm_volcano.png`, `20260724_fdr_off_by_one.md`).
- Always add a sibling `YYYYMMDD_<slug>.md` note with:
  1. **issue** or **showcase**;
  2. the exact `kycg` / `knowYourCG` command that produced it (must reproduce);
  3. for an issue: expected vs actual. For a showcase: one line of biological
     interpretation.
- Keep it small: a PNG/PDF plus a short TSV and the note. Reference large inputs
  by path — do not copy them in.

## Example

```
20260724_cgi_enrichment.png        # the figure
20260724_cgi_enrichment.md         # showcase; command + interpretation
20260724_cgi_enrichment.tsv        # the small kycg test output behind it
```
