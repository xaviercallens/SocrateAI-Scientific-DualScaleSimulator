# TopoDB — a queryable record of every persistent-homology computation this project trusts

**Built 2026-09-20.** Branch `loop/topodb`. The database file lives outside git at
`/mnt/disks/disk-socrateai-local-1/topodb/topodb.sqlite` (SQLite, WAL), with dated backups beside it.

## Why it exists

Four separate campaigns had produced persistence results — cosmology, K3 × T², quantum fluids, genetics —
in incompatible JSON layouts, each with its own conventions for what counts as a control and what a
p-value means. TopoDB makes them one comparable corpus, and enforces in code the rules those campaigns
learned the hard way.

## The rules the schema enforces

- **A p-value is refused unless `null_model` and `n_null` are given.** Across 557 runs written by four
  independent agents, **0 rows carry a p-value without its null**. Where a source had no null, the
  statistic is stored alone.
- **`matches` exists only where an expectation was pre-stated.** 901 Betti rows match a pre-stated
  expectation, 55 mismatch (and are kept), 398 are ungated.
- **Tier is mandatory** and restricted to A/B/L/C/X; method is a closed list.
- Every run records its script, exact command, seed, git commit and branch, wall time and peak memory.

## Contents (2026-09-20)

| | |
|---|---|
| datasets | 237 — biology 98, mathematics 60, synthetic 42, quantum_fluid 22, astro 15 |
| runs | 557 — alpha 176, chain_complex 153, rips 148, lower_star_graph 58, cubical 18, witness 4 |
| tiers | B (exact arithmetic) 159, X (numerics) 398 |
| bars | 7 344 · Betti rows 1 354 · statistics 3 800 (212 with a p-value) |
| controls | 635 — 502 passed, 86 failed, 47 not run |
| findings | 236 — 159 recovered, 28 null, 26 failed, 18 inconclusive, 5 artefact |

**The failed controls and the mismatches are the point.** A corpus where everything passes would mean the
controls were not doing any work.

## Use

```bash
export TOPODB_PATH=/mnt/disks/disk-socrateai-local-1/topodb/topodb.sqlite
python -m topodb.cli summary
python -m topodb.cli search "vortex"
python -m topodb.cli datasets --domain biology
python -m topodb.cli runs --tier B
python -m topodb.cli show 331 --bars 10
python -m topodb.cli nearest 331 --k 5 --dim 1      # topologically similar runs
python -m topodb.cli distance 12 44 --dim 1 --metric bottleneck
python -m topodb.cli export --json /tmp/topodb.json
```

Ingestion from Python is `topodb.api.TopoDB`; see its module docstring.

### Limits of `nearest` / `distance`, stated by the tool itself
Bars are truncated to the top 50 per dimension on ingest (and several sources truncated further before
committing). Essential bars are dropped by default, since the schema stores `death = NULL`. Filtration
units differ across methods, so runs of different methods are not mixed by default; 321 of the
backfilled runs have no bars at all (chain-level computations have no filtration, hence no barcode).

## What is inside, by origin

- **Backfill of committed results** (`topodb/ingest/`, 406 runs): the known-answer suite, the fixed-library
  validation, quantum fluids, genetics, crystallography, K3 × T² (both rounds and the chain-level study),
  cosmic vorticity, and the round-2 CMB run. 102 skips are listed with reasons in `BACKFILL_NOTES.json` —
  the dominant reason being that the source saved Betti *curves* or summary statistics, never birth/death
  pairs, so no barcode exists to ingest.
- **New runs** (`topodb_runs/{astro,quantum_fluid,biology}/`): see each directory's `report.json`.

## Discrepancies the corpus deliberately keeps

- **K3's second Betti number is 22 by exact chain-level computation (tier B) and 27 by point sampling
  (tier X, `expected=22`, `matches=0`).** Both are stored. The sampled value is wrong, and recording it is
  what makes the limits of sampled TDA legible.
- **The original sky complex gives b₂ = 49 147 on a sphere; the fixed one gives (1, 0, 1).** Both runs are
  stored, the former marked superseded.
- The Rips reference engine disagrees with the alpha pipeline on five known-answer cases; both are stored.
- E. coli Hi-C passes on the alpha path and fails on the Rips cross-check of the same matrix.

## Known defects in the corpus, recorded rather than hidden

- Two house null designs (the Hi-C "linear" null and the scRNA gene permutation) return their floor
  p-value for *every* case, controls included: they destroy all structure, not the structure under test.
  The affected p-values were withdrawn.
- A pre-registered 4.5 Å band in the protein statistic does no work: the unbanded statistic gives the
  same AUC. The statistic measures H₁ classes per residue, not a selected scale.
- The influenza site-bootstrap null is ill-posed for reassortment (resampling columns leaves each genome
  intact). The failed run is kept; a pooled null, declared post-hoc, reproduces the published result.

## Operational note

The database is shared. During the first build three agents wrote concurrently and one wholesale
replacement destroyed another agent's rows; run-id reuse then caused a correction script to touch rows it
did not own. Nothing was lost, because every run is reconstructible from a committed result JSON.
**Conventions now in force:** never replace the file; scope any reset to rows you own (the backfill uses
`_backfill_source` in `params_json` plus a sidecar manifest); never resolve a run id by SQL lookup in a
shared database — remember the id you just created.
