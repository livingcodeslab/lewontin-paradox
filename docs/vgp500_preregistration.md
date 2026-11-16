# VGP-500 Heterozygosity Survey — Pre-registered Analysis Plan

**Project:** Lewontin's Paradox / gBGC-saturation test — cross-species empirical arm
**Status:** pre-registration / work plan (v1, 2026-11-16)
**Authors:** Erik Garrison (PI), Shuo (heterozygosity sweep)
**Artifact of record for the hypothesis:** `manuscript.typ` (mechanism + formal refutation of the gene-conversion alternative). This document pre-registers the *empirical* arm that the manuscript currently lacks.

---

## 1. Why this exists

The theory is settled in this repo. The original mechanism (unbiased gene conversion
homogenization) is **refuted** — it is a martingale and leaves `E[π] = 4Nₑμ` unchanged
(Wiuf & Hein; Lean-formalized in `lean/src/`, zero `sorry`/`axiom`). The live hypothesis:

> **gBGC**, always-on and selection-like, plus mutation–selection–drift balance at
> weak/strong (W/S) sites, produces an **`Nₑ`-independent diversity floor `2μ/(uδ)`**
> of the right magnitude (~0.01 for `u_eff≈10⁻⁵`, `δ≈0.18`, `μ≈10⁻⁸`).

What we have so far is consistent but under-powered:

- **Tier 1** (Buffalo 172-taxon compiled catalog): rising-then-plateau; low-`Nc`
  slope +0.118 (p=0.002), high-`Nc` +0.018 (p=0.58). But AIC **ties** a saturating
  power law with a plain power law (ΔAIC≈0) — the compiled data cannot discriminate.
- **Tier 3** (composition proxy, 135 assemblies): null/uncertain (GC3 class-fixed
  slope includes zero; predeclared concavity diagnostic came out convex).
- **VGP pilot** (this repo, `analysis/vgp_*`): a *modality* — same-individual
  phased-assembly (hap1 × hap2) divergence as a direct π estimator — validated on a
  bounded 3-pair core, with one hard-won lesson: **collapsed paralogy inflates π by
  ~25–30%** (P07 *Spinachia*: 0.00117 → ≈0.00084 after triage;
  `analysis/vgp_symmetric_read_test_results_v1.md`).

**The pre-registered claim of this study:** at appropriate power — a uniform,
primary-data π survey across the full VGP phased-diploid catalog (~500 species) —
the predicted saturation of π with census size becomes visible. We have not yet had
the appropriate power. This document fixes the predictions, the data, the pipeline,
and the statistics *before* the data are collected.

## 2. Pre-registered predictions

Let `π_i` = within-individual autosomal π for species *i* (1:1-filtered
hap1-vs-hap2 assembly divergence, exclusion beds applied). Let `Nc_i` = census-size
proxy (from life-history databases, see §5.6).

- **P1 (plateau).** `log π` vs `log Nc` is better described by a saturating
  (rising-then-plateau) model than by a single linear log-log or power-law fit, and
  the high-`Nc` segment slope is statistically indistinguishable from zero.
- **P2 (floor magnitude).** The plateau level is `2μ/(uδ)` — order 0.001–0.02, i.e.
  the observed asymptote should not exceed ~2× the maximum plausible floor.
- **P3 (composition signature).** Species in the plateau regime show elevated
  GC3 relative to the low-`Nc` regime, in the gBGC direction, controlling for
  phylogeny. (Direction predeclared; magnitude not claimed — Tier 3 showed this
  signal is noisy and class-confounded.)
- **P4 (demographic noise).** Residual variation around both segments is partly
  explained by PSMC-derived historical `Nₑ` (recent bottlenecks push π below the
  segment fit; the floor bounds it from below).

**Falsifiers.** A single log-log linear fit wins decisively (ΔAIC > 10 in favor);
the high-`Nc` slope is positive and significant; the asymptote is far above any
plausible floor (e.g. >0.1); or no composition signature in either regime. Any of
these counts *against* the hypothesis, and we commit to reporting it as such.

**Honest boundaries (carried over from the pilot's independent review).** One phased
individual per species ⇒ this design measures *individual* heterozygosity, not
population π. That is adequate for the cross-species plateau test (the π range
across vertebrates is ~100×; individual sampling noise is small against it) but
makes **no within-species population claims**. Also: census-size proxies are
life-history surrogates, not censuses — the x-axis is ordinal-ish, which we handle
in §7.

## 3. Species set and inclusion criteria

Frame: **all VGP phased diploid assemblies** (hap1 + hap2 primary contigs),
prioritizing:

1. the Freeze-1 mirror already on disk (`/moosefs/erikg/vgp/freeze1/objects/`,
   catalogued in `analysis/vgp_freeze1_mirror_manifest.tsv`), and
2. later-phase VGP releases as they are mirrored.

Inclusion gate (pre-registered, before π estimation):

- G1. hap1 and hap2 both present, non-degenerate (≥90% of estimated genome in
  contigs ≥ some N50 floor; record, don't guess).
- G2. HiFi reads available (needed for k-mer profile, merqury, and the symmetric
  read triage).
- G3. BUSCO completeness ≥ 85% (lineage-appropriate) **and** duplication reported.
- G4. Not blacklistable for ploidy ambiguity (SmudgePlot/polyploid suspects).
- Exclusions and the reason for each go to a **failure ledger**
  (`failure_ledger.tsv`, same discipline as `analysis/tier3c_failure_ledger.tsv`).
  Nothing is silently dropped.

Target: every species that passes. **The authoritative species roster is
Shuo's VGP catalog** — he maintains the exact genome list and counts. For
reference only, the Freeze-1 mirror on disk
(`/moosefs/erikg/vgp/freeze1/objects/`, `analysis/vgp_freeze1_mirror_manifest.tsv`)
holds 581 assembly accessions, 483 with assembly reports (264 principal
pseudohaplotype, 71 hap1, 25 hap2, 13 maternal/paternal, 2 alternate, 108 plain
haploid); phased pairs arrive as separate h1/h2 accessions (e.g. P01:
`mCamDro1.hap1` = GCA_036321535.1, `mCamDro1.hap2` = GCA_036321565.1). The
working species list is whatever Shuo's catalog says passes G1–G4.

## 4. Pipeline (what to run — "collect everything now, decide later")

Backend decision (operator, 2026-11-16): **sweepga with the fastga aligner,
sensitivity-first** (`--overlap 0`, super-sensitive seeding), then **exact 1:1
multiplicity filtering** downstream — the filter does the cleanup, not the aligner.
wfmash is retained as the audit/spot-check backend (see §6). Known failure mode
from the pilot wave: fastga's pair-file sweep IO failed deterministically above
~2.1 Gbp query (P09×2, P06, P10, P08); for those genomes fall back to wfmash
full-node (`mapping_stage_wfmash.sh`) or chunked query (parked
`vgp_chunked_mapping.py`). Record the backend per pair in the execution record —
and pre-register here that **backend is a covariate we check, not a filter we hide**.

Per species, collect ALL of the following (each item idempotent, digest-recorded,
restartable — this is the "running now" sweep):

| # | Artifact | Tool / in-repo machinery | Purpose |
|---|---|---|---|
| 5.1 | hap1×hap2 native PAF, 1:1 | `sweepga --aligner fastga --num-mappings 1:1 --scaffold-jump 0 --overlap 0 --scoring log-length-ani --min-aln-length 25k` (template: `analysis/slurm/vgp_10_pilot/mapping_stage_*`) | raw divergence evidence |
| 5.2 | exact-1:1 enforced PAF + multiplicity audit | `python3 -m analysis.vgp_10_pilot enforce-paf` / `audit-paf` | the filter that cleans up fastga's sensitivity |
| 5.3 | exclusion beds (1:1, not-1:1, N/gap, low-complexity, h2-projected) | in-script derivation in `mapping_stage_wfmash.sh` (reuse verbatim) | callable-space definition |
| 5.4 | variant calls from the 1:1 PAF | `minimap2 -x asm5 --cs` re-alignment of the 1:1 blocks + `paftools.js call`-style cs-tag parsing (add to `analysis/vgp_10_pilot`) — call SNVs/indels, never estimate π from unfiltered PAF | π numerator |
| 5.5 | k-mer profile | `meryl` k=21,31 on HiFi reads; `genomescope2` (VGP pipeline stage 1) | independent heterozygosity + repeat fraction + error model |
| 5.6 | merqury hap comparison | `merqury`/`merquryfk` reads-vs-hap1/hap2 k-mer concordance | assembly QC + collapse screen |
| 5.7 | BUSCO | lineage-appropriate | completeness (G3) + duplication as collapse proxy |
| 5.8 | GC / GC3 | `analysis/tier3c_ncbi_gc.py` machinery re-pointed at VGP assemblies | P3 composition covariate |
| 5.9 | PSMC(μ-scaled) | per-pair PSMC machinery already in the repo (`analysis/slurm/vgp_10_pilot/`) | P4 demographic covariate |
| 5.10 | assembly stats | contiguity, size, N50 (BUSCO-adjacent script) | QC + covariates |
| 5.11 | life-history / census proxy row | PanTHERIA-style databases + expert coding; body mass, range, trophic level, IUCN estimate | `Nc_i` (§2) — **non-genome half, start early, it is slow** |

Storage layout mirrors the pilot: per-pair `mapping/`, `kmers/`, `busco/`, `psmc/`
under a versioned wave root (e.g. `vgp500-wave-v1`) with `.complete.json` packets,
SHA-256 manifests, and command records — same contract as
`mapping_stage_wfmash.sh` produces. Idempotency rule: any stage re-run must
short-circuit on an existing digest-complete packet.

## 5. Variant calling & π estimation (beyond the alignment)

1. Re-align each 1:1 syntenic block with `minimap2 -x asm5 --cs` (or parse cs
   directly from sweepga output if provided) and call SNVs/indels per block.
2. Apply the exclusion beds (5.3): callable space = 1:1 ∩ non-N ∩ non-low-complexity.
3. Filter calls: block length ≥ 25k, identity ≥ 95% *for the variant-call step*
   (stricter than the mapping identity — this is where `-p 95` belongs), per-block
   depth/mappability sanity.
4. `π = observed differences / callable bases` per species, with a block-level
   bootstrap for CIs. Pre-register the estimator: **π is computed only on
   callable 1:1 space, post-exclusions** — never on raw PAF.
5. Cross-estimator triage (the scaled-up paralog screen — this replaces the
   pilot's per-pair manual adjudication):
   - `π_align` (5.4) vs `π_kmer` (GenomeScope2 heterozygosity) disagreement
     > 2× or > 0.005 absolute ⇒ flag.
   - merqury hap-concordance anomalies ⇒ flag.
   - BUSCO duplication > 2× median-of-class ⇒ flag.
   - Flagged species get the symmetric read test (pilot machinery) or are
     reported in a sensitivity analysis with and without.
   Pre-registered rule: **the primary analysis excludes flagged species; the
   sensitivity analysis includes them; both are reported.**

## 6. Validation subset

On ~30 species spanning the π range and genome-size range (include the pilot's
P01–P10, the Tier-3A trio, both backends):

- run **both** fastga and wfmash full pipelines; report π agreement per pair
  (the pilot found 0.25 variant Jaccard across backends *pre*-filtering — this
  subset measures how much of that survives exact 1:1 + variant filtering);
- run the symmetric read test to calibrate the §5.5 triage flags against ground
  truth (P07 is a known-positive).

## 7. Analysis plan (pre-declared)

- Model comparison on `log π` vs `log Nc`: (a) single power law, (b) saturating
  power law / rising-then-plateau (two-segment or Michaelis–Menten-like),
  (c) linear. AIC weights + explicit ΔAIC thresholds declared in §2.
- **Phylogenetic non-independence is mandatory:** PGLS (or phylogenetic paired
  contrasts) under a published vertebrate timetree. Report both raw and
  phylogenetically corrected.
- Segment slopes with CIs; test high-`Nc` slope = 0 (P1). Asymptote vs
  plausible `2μ/(uδ)` range (P2).
- Covariate models: GC3 (P3), PSMC-derived recent `Nₑ` (P4), repeat fraction,
  genome size, sequencing depth, backend — each entered as a covariate, with
  pre-declared direction.
- `Nc` proxy uncertainty: rank-based / ordinal treatment and a proxy-subset
  sensitivity analysis (e.g., IUCN-estimate subset only).
- Every number in the eventual manuscript must trace to a digest-complete
  artifact; the repo's audit discipline (headline recomputation at 1e-12,
  `results/tier3/headline_audit.tsv` style) carries over.

## 8. Work plan

| When | What | Who |
|---|---|---|
| now | Start the collection sweep: k-mer profiles + BUSCO + assembly stats for every candidate assembly (5.5–5.10) — cheapest first, no dependencies | Shuo |
| now | hap1×hap2 alignments + 1:1 filtering for the pilot-known pairs first (they validate the scaled pipeline), then the full set in size order | Shuo |
| now | Life-history / `Nc` proxy coding (5.11) — slowest external dependency, start immediately | Erik + Shua |
| +1 wk | QC gate pass/fail list + failure ledger for the species set (§3) | Shuo |
| +2 wk | Variant calling + π estimator + triage flags on the first ~50 species | Shuo + Erik |
| +3–4 wk | Validation subset cross-backend + symmetric read test calibration (§6) | Shuo |
| +5–6 wk | Full π table; freeze; analysis per §7; manuscript section drafted into `manuscript.typ` | Erik |

## 9. Provenance discipline (non-negotiable, carried from the pilot)

- Every stage: command record + SHA-256 of inputs/outputs + `.complete.json`;
  re-runs idempotent; failures to a ledger, never silent.
- Backend recorded per pair; no mixing without a covariate check.
- No π number is reported that isn't reproducible from a digest-complete packet.
- The pilot's evidence-correction history
  (`analysis/vgp_evidence_correction_v1/`) is the precedent: when a correction
  is found, it supersedes via a dated packet, and the paper quotes only the
  latest supersession.

## 10. What this study is not

It is not a within-species population-genetics study (one individual per
species). It is not a claim that gBGC is established by association. It is the
adequately-powered version of the test that Tier 1 and Tier 3 could only gesture
at — with predictions, exclusions, estimators, and statistics fixed in advance,
here, before the data exist.
