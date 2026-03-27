# CL-KG × CellxGene Census Integration

This repo contains the proof-of-concept and roadmap for connecting the **CL-KG knowledge graph** to the **CellxGene Census** via roaring bitmaps, enabling fast, provenance-rich cell-set queries at runtime.

---

## The core idea

The CL-KG graph maps author cell type annotations (`Cell_cluster` nodes) to CL terms via `composed_primarily_of` edges, anchored to `Dataset` nodes. Each `Cell_cluster` represents a named cluster from a source study.

The goal is to extend these nodes with **`soma_joinid` bitmaps** — compact representations of the exact cells in each cluster, indexed against a given Census build. This unlocks:

> *"Given a Census obs filter (disease, tissue, assay, …), which source datasets and author labels are relevant? Return matching cells with full provenance (cluster label, DOI, CL term)."*

This is the primary use case, demonstrated end-to-end in this repo.

---

## How it works

### Graph build (one-time per dataset)
1. Download source H5AD via `cellxgene_census.download_source_h5ad(dataset_id)`
2. Extract `observation_joinid → author_label` pairs from H5AD obs
3. Query Census obs for the same `dataset_id` → get `(observation_joinid, soma_joinid)` pairs
4. Inner join → resolve `soma_joinid` set per `Cell_cluster`
5. Encode as roaring bitmap; store on node keyed by node IRI

### Cache refresh (per Census release)
For each `Cell_cluster` with `observation_joinids` stored:
1. Join `observation_joinid` → `soma_joinid` against the new Census build
2. Re-encode bitmap; record Census version
3. Log join rate — a drop >5% flags potential `observation_joinid` instability

### Query time
```
CL term → graph → Cell_cluster subtree → soma_joinid bitmaps
  → decompress → coords array
  → Census axis_query(coords=soma_joinid_array, var_query=gene_set)
  → expression matrix with full provenance
```

Using `coords=` (integer index lookup) rather than `value_filter` keeps queries fast regardless of Census size.

---

## PoC results

The proof-of-concept ran against 8 macrophage/lung clusters from 2 COVID datasets, queried from CL-KG (`export.csv`).

**Phase 1 — Bitmap build** (`pilot_bitmap_build.py`, Census `2025-11-08`):
- 8/8 clusters built, **100% join rate** across all clusters
- Node IRI sufficient as bitmap key — no separate UUID needed
- `obs_meta` fallback correctly resolved curated-merge case to `cell_type_fine`

**Phase 2 — Scenario B query demo** (`demo_scenario_b.py`, filter: COVID-19):

| Cluster label | Cells | Passing filter | Fraction |
|---|---|---|---|
| Macro (Ren et al.) | 21,471 | 21,445 | 99.9% |
| Monocyte-derived macrophages (Liao et al.) | 9,534 | 7,504 | **78.7%** |
| Macro_c1–c6 (Ren et al.) | 21,471 | ~21,444 | >99% |
| **Total unique** | **31,005** | **28,949** | **93.4%** |

The 78.7% figure for Liao et al. reflects non-COVID donors in a broad lung atlas — a real biological difference surfaced automatically by the filter.

Each output row carries: `soma_joinid`, `cell_type`, `disease`, `tissue`, `cluster_label`, `doi`.

**Phase 2b** (`demo_scenario_b.py`): union bitmap → Census `get_anndata()` with gene filter → annotated AnnData saved as `macro_lung_annotated.h5ad` (31,005 cells, `author_label` populated from bitmap lookup alone — no H5AD re-download at query time).

---

## Key files

| File | Purpose |
|---|---|
| `test_joinid.py` | Step-by-step validation of the `observation_joinid` join workflow |
| `pilot_bitmap_build.py` | Bitmap build for 8 macrophage/lung clusters from `export.csv` |
| `demo_scenario_b.py` | End-to-end Scenario B query — change `OBS_FILTER` to try other filters |
| `pilot_output.txt` | Bitmap build summary output |
| `export.csv` | Source graph query results (macrophage/lung, 8 clusters, 2 datasets) |
| `macro_lung_annotated.h5ad` | Annotated AnnData output from Phase 2b |
| `ROADMAP.md` | Detailed architecture, schema, and implementation plan |
| `census-joinid-pipeline-test.md` | Session notes and architecture decisions from initial join test |

---

## Running the demo

```bash
uv sync
# Try different Census obs filters by editing OBS_FILTER at the top of the script
python demo_scenario_b.py
```

Requires network access to CellxGene Census (public S3). No Sanger network access needed for the demo — only Phase 3 (graph write-back) requires Bolt access to `cl-kg-neo4j-db`.

---

See `ROADMAP.md` for the full architecture, proposed schema extensions, and implementation phases.
