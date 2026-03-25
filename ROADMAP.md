---
date: 2026-03-19
status: draft
tags:
  - cl-kg
  - census
  - cellxgene
  - roadmap
  - pilot
---

# CL-KG × Census Integration — Roadmap

## Background

The CL-KG graph maps author cell type annotations (`Cell_cluster`) to CL terms (`Cell`) via `composed_primarily_of` relationships, anchored to `Dataset` nodes. The goal of this work is to extend each `Cell_cluster` node with a set of `observation_joinid` values drawn from the source H5AD, then use these to resolve `soma_joinid` bitmaps against each Census build — enabling fast, coords-based expression retrieval and cell-set queries at runtime.

See `census-joinid-pipeline-test.md` for the validated proof-of-concept join workflow.

---

## Current graph schema (CellMark / cl-kg-neo4j-db)

```cypher
(Cell_cluster {label, iri})
  -[:composed_primarily_of]-> (Cell {curie, label})
  -[:subcluster_of*0..]-> (Cell_cluster)        // cluster hierarchy, variable depth
  -[:tissue]-> (Class {label})                  // anatomy association (OWL Class nodes)
  -[]-> (Dataset {title})
```

Example query — macrophage subclusters in lung:
```cypher
MATCH P=(c:Cell { label: 'macrophage'})<-[:composed_primarily_of]-(clus_up:Cell_cluster)
        <-[:subcluster_of*0..]-(clus_down:Cell_cluster)
MATCH (clus_down)-[:tissue]-(anat:Class { label: 'lung'})
RETURN P
```

The `subcluster_of` hierarchy means a Scenario B query returns a **subtree** of `Cell_cluster` nodes. The bitmap union across all nodes in that subtree is the full cell set for the query — e.g. "all macrophage subclusters in lung tissue".

**Connectivity:** Bolt on port `7687`, Sanger network only (all external ports blocked).
Public host `cl-kg-neo4j-db.cellgeni.sanger.ac.uk:443` serves the Neo4j browser UI only.
Access from outside Sanger requires VPN or SSH tunnel.

---

## Dataset node — what the graph already provides

Confirmed from `export.csv` (macrophage/lung query results):

| Dataset property | Example value | Notes |
|---|---|---|
| `download_link` | `https://datasets.cellxgene.cziscience.com/75c059c8-...h5ad` | The UUID is the Census `dataset_id` |
| `publication` | `https://doi.org/10.1038/s41586-021-03569-1` | DOI for provenance |
| `obs_meta` | `[[{"field_name":"cell_type_fine","field_type":"author_cell_type_label"}],...]` | Lists all author annotation obs columns and their hierarchy levels |
| `title` | `A molecular single-cell lung atlas of lethal COVID-19` | Human-readable title |

### Recovering the source obs column per Cell_cluster

The obs column that sourced a given cluster label is stored on the `Cell_cluster` node as a **non-standard extra property key**:
- Rows with `celltype` key → label came from `obs["celltype"]`
- Rows with `majorType` key → label came from `obs["majorType"]` (coarser hierarchy level)
- **No extra key present** → cluster was curated/merged; winning label came from one of the `obs_meta` fields but is not stamped on the node

The no-extra-key case requires trying each `obs_meta` field name until one matches — or storing the winning field name in the graph at build time. **This is a recommended schema improvement.**

### Pilot datasets (from export.csv)

| dataset_id | title | clusters in query |
|---|---|---|
| `75c059c8-8fb7-4e6e-a618-a3e01ac42060` | A molecular single-cell lung atlas of lethal COVID-19 (Liao et al. Nature 2021) | "Monocyte-derived macrophages" (obs col: `cell_type_fine`, curated) |
| `891c6eef-a243-41a7-8b47-c8c1a224a193` | Large-scale single-cell analysis... COVID-19 patients (Ren et al. Cell 2021) | Macro_c1–c6 (obs col: `celltype`), "Macro" (obs col: `majorType`) |

---

## Proposed schema extension

```cypher
(Cell_cluster {label, iri, obs_column: str, observation_joinids: [str]})
  -[:composed_primarily_of]-> (Cell {curie, label})
  -[:subcluster_of*0..]-> (Cell_cluster)
  -[:tissue]-> (Class {label})
  -[:has_source]-> (Dataset {title, census_dataset_id: str, download_link: str,
                              publication: str, obs_meta: str,
                              census_version_cached: str})
```

| New property | Node | Notes |
|---|---|---|
| `obs_column` | `Cell_cluster` | Explicit obs column name used to source this cluster's label — resolves the curated-merge ambiguity |
| `observation_joinids` | `Cell_cluster` | Array of stable cross-build IDs from source H5AD obs |
| `census_dataset_id` | `Dataset` | Resolved Census `dataset_id` (differs from `dataset_version_id` in `download_link`) |
| `census_version_cached` | `Dataset` | Census build against which current bitmaps were generated |

Bitmap store keyed by `(node_iri, census_version)` — no separate UUID needed.

---

## Architecture

### Graph build (per new dataset)
1. Download source H5AD via `cellxgene_census.download_source_h5ad(dataset_id)`
2. Read `obs` only (backed mode — no expression matrix loaded)
3. Extract `observation_joinid → author_label` pairs
4. Write `observation_joinids` array and `census_bitmap_uuid` to each `Cell_cluster` node
5. Record `dataset_id` and `collection_doi` on `Dataset` node

### Post-build cache job (per Census release)
For each `Cell_cluster` node with `observation_joinids`:
1. Query Census obs for `dataset_id` → get `(observation_joinid, soma_joinid)` pairs
2. Inner join on `observation_joinid` → resolve `soma_joinid` array
3. Encode as roaring bitmap
4. Store bitmap keyed by `census_bitmap_uuid` (+ Census version)
5. Update `census_version_cached` on `Dataset` node
6. Log join rate — flag nodes where rate drops >5% vs prior build (signals `observation_joinid` instability)

Cache job is obs-only, no expression data. Parallelisable per dataset.

### Query time
```
CL term → graph → Cell_cluster nodes → census_bitmap_uuid
  → decompress bitmap → soma_joinid array
  → Census axis_query(coords=soma_joinid_array, var_query=gene_set)
  → expression matrix
```

---

## Use case scenarios

### Scenario A — annotation update via Census filter
> *"Which author annotations are subsumed by this CL term / obs filter?"*

1. Apply Census obs filter → get `soma_joinid` set S
2. For each `Cell_cluster` bitmap: intersect with S (inner join)
3. Nodes where intersection is large → annotation is consistent with filter
4. Use to validate or update `composed_primarily_of` edges

### Scenario B — author annotation discovery (more imminent)
> *"Given a Census obs filter, which source datasets and author labels are relevant? Return cells with provenance."*

1. Apply obs filter on graph → find `Cell_cluster` nodes with partial overlap
2. Decompress bitmaps for those nodes → `soma_joinid` arrays
3. Query Census obs (no expression) with coords + obs filter → get intersection
4. Return result rows enriched with `collection_doi` from `Dataset` node
5. Optionally: download H5AD for those cells using intersection `soma_joinid` set

> **Join type:** inner join throughout. Left-outer (graph as left) is useful for QA only
> (surfaces cells in author annotation absent from Census — confirms QC exclusions).

---

## Pilot plan

### Phase 0 — Infrastructure (now)
- [x] Validate `observation_joinid` join workflow (`test_joinid.py`, dataset `9bb9596d`)
- [x] Confirm `download_source_h5ad()` works
- [x] Add `neo4j` and `cxg-query-enhancer` to project dependencies
- [x] Real graph query results in hand (`export.csv` — macrophage/lung, 8 clusters, 2 datasets)
- [ ] Confirm Bolt port and access route for `cl-kg-neo4j-db` (7687, internal only — need VPN or SSH tunnel)
- [ ] Verify `observation_joinid` stability across two Census builds

### Phase 1 — Offline bitmap build (no DB connection needed)
Using `export.csv` as a stand-in for a live graph query:
- [x] Run `pilot_bitmap_build.py` against the two pilot datasets
- [x] Download H5ADs, extract `observation_joinid` sets per cluster label + obs column
- [x] Join to Census → `soma_joinid` bitmaps
- [x] All 8/8 clusters built, **100% join rate** across all clusters
- [x] IRIs extracted and used as bitmap keys — no separate UUID needed

**Pilot results (Census stable 2025-11-08):**

| Cluster label | obs_column | Cells | Join rate | Bitmap (bytes) | DOI |
|---|---|---|---|---|---|
| Monocyte-derived macrophages | `cell_type_fine` (fallback) | 9,534 | 100% | 16,810 | 10.1038/s41586-021-03569-1 |
| Macro_c1-C1QC | `celltype` | 11,221 | 100% | 10,088 | 10.1016/j.cell.2021.01.053 |
| Macro_c2-CCL3L1 | `celltype` | 6,727 | 100% | 9,816 | 10.1016/j.cell.2021.01.053 |
| Macro_c3-EREG | `celltype` | 1,030 | 100% | 2,100 | 10.1016/j.cell.2021.01.053 |
| Macro_c4-DNAJB1 | `celltype` | 1,282 | 100% | 2,636 | 10.1016/j.cell.2021.01.053 |
| Macro_c5-WDR74 | `celltype` | 671 | 100% | 1,470 | 10.1016/j.cell.2021.01.053 |
| Macro_c6-VCAN | `celltype` | 540 | 100% | 1,120 | 10.1016/j.cell.2021.01.053 |
| Macro | `majorType` | 21,471 | 100% | 14,654 | 10.1016/j.cell.2021.01.053 |
| **Union** | — | **31,005** | — | — | — |

Note: `dataset_version_id` (UUID in graph `download_link`) ≠ Census `dataset_id`. Resolution required: try `dataset_version_id` direct match in Census datasets table, fall back to DOI match. Both resolved cleanly here.

Also confirmed: `obs_meta` fallback correctly resolved the curated-merge case (no extra key on node) to `cell_type_fine`.

**Schema simplification confirmed:** node `iri` is sufficient as bitmap key — no separate `census_bitmap_uuid` needed if graph build and bitmap build are co-scheduled.

### Phase 2 — Scenario B query demo
- [x] Take union bitmap (31,005 soma_joinids)
- [x] Query Census obs with coords + obs filter → intersection
- [x] Return result DataFrame with per-row DOI and cluster label

**Demo results — filter: `disease_ontology_term_id == 'MONDO:0100096'` (COVID-19)**

Script: `demo_scenario_b.py`

| Cluster label | obs_column | Cells in cluster | Passing filter | Fraction | DOI |
|---|---|---|---|---|---|
| Macro | `majorType` | 21,471 | 21,445 | 99.9% | 10.1016/j.cell.2021.01.053 |
| Macro_c1-C1QC | `celltype` | 11,221 | 11,216 | 99.9% | 10.1016/j.cell.2021.01.053 |
| Monocyte-derived macrophages | `cell_type_fine` | 9,534 | 7,504 | **78.7%** | 10.1038/s41586-021-03569-1 |
| Macro_c2-CCL3L1 | `celltype` | 6,727 | 6,708 | 99.7% | 10.1016/j.cell.2021.01.053 |
| Macro_c4-DNAJB1 | `celltype` | 1,282 | 1,282 | 100% | 10.1016/j.cell.2021.01.053 |
| Macro_c3-EREG | `celltype` | 1,030 | 1,030 | 100% | 10.1016/j.cell.2021.01.053 |
| Macro_c5-WDR74 | `celltype` | 671 | 670 | 99.8% | 10.1016/j.cell.2021.01.053 |
| Macro_c6-VCAN | `celltype` | 540 | 539 | 99.8% | 10.1016/j.cell.2021.01.053 |
| **Total unique** | — | **31,005** | **28,949** | **93.4%** | — |

**Notable result:** "Monocyte-derived macrophages" (Liao et al.) passes at only 78.7% vs >99% for all Ren et al. clusters. Liao et al. is a broad lung atlas including non-COVID donors; Ren et al. is COVID-only. The ~21% gap reflects non-COVID cells in the Liao cluster — a real biological/study-design difference surfaced automatically by the filter.

Each output row carries: `soma_joinid`, `cell_type` (CL-mapped), `disease`, `tissue`, `cluster_label` (author annotation), `doi`. Full provenance from graph node → individual cell → source paper.

**Files to inspect results:**
- `demo_scenario_b.py` — full demo script; change `OBS_FILTER` constant at top to try other Census obs filters (disease, tissue, assay, organism, etc.)
- `pilot_output.txt` — bitmap build summary from Phase 1
- `export.csv` — source graph query results (macrophage/lung, 8 clusters, 2 datasets)
- Cached H5ADs in repo root:
  - `d8da613f-e681-4c69-b463-e94f5e66847f.h5ad` — Liao et al. COVID lung atlas (Census id for `75c059c8`)
  - `9dbab10c-118d-496b-966a-67f1763a6b7d.h5ad` — Ren et al. COVID immune cells (Census id for `891c6eef`)
  - `9bb9596d-f23f-4558-912f-d4dc7d52721b.h5ad` — CGE interneurons MOp (original join test)

### Phase 2b — Annotated AnnData via Census stream
- [x] Union bitmap → Census `get_anndata()` with 3-gene var filter (MARCO, C1QC, CCL3L1)
- [x] `author_label` column added to obs via `soma_joinid` bitmap map — no source H5AD needed at query time
- [x] Result saved to `macro_lung_annotated.h5ad` (2.2MB, 31,005 cells × 2 genes measured)

```
AnnData object with n_obs × n_vars = 31005 × 2
    obs: 'soma_joinid', 'cell_type', 'disease', 'tissue', 'dataset_id', 'author_label'
    var: 'soma_joinid', 'feature_id', 'feature_name', ...
```

Note: only 2 genes returned (MARCO, C1QC) — CCL3L1 not measured in either dataset. `author_label` populated for all 31,005 cells from bitmap lookup alone.

**Inspect:** open `macro_lung_annotated.h5ad` in any AnnData-compatible tool (scanpy, cellxgene Explorer, etc.)

---

### Phase 3 — Graph write (requires Sanger network access)
- [ ] Write `obs_column` and `observation_joinids` to `Cell_cluster` nodes
- [ ] Write `census_dataset_id` and `census_version_cached` to `Dataset` nodes
- [ ] Schema improvement: ensure `obs_column` is stored for all nodes, including curated-merge cases

### Phase 2 — Cache job prototype
- [ ] Script: for a single dataset, join `observation_joinid` → `soma_joinid` against Census stable
- [ ] Encode result as `pyroaring.BitMap`
- [ ] Store bitmap (file or Redis) keyed by `census_bitmap_uuid`
- [ ] Measure join rate, confirm >95%

### Phase 3 — Query demo (Scenario B)
- [ ] Given a CL term: find `Cell_cluster` nodes via graph
- [ ] Decompress bitmap → `soma_joinid` array
- [ ] Query Census obs with coords + an example obs filter (e.g. tissue, disease)
- [ ] Return result DataFrame with `collection_doi` column

### Phase 4 — Scenario A demo
- [ ] Run Census obs filter → `soma_joinid` set
- [ ] Intersect with bitmaps from graph → identify subsumed `Cell_cluster` nodes
- [ ] Compare against existing `composed_primarily_of` edges — report consistency

---

## Open questions

| Question | Notes |
|---|---|
| Is `observation_joinid` stable across Census builds? | Critical assumption — verify before Phase 2 |
| Where does the bitmap cache live? | File store vs Redis vs property on node — depends on query latency needs |
| Does `Cell_cluster` always map 1:1 to a single `dataset_id`? | If a cluster spans datasets, `observation_joinids` join logic needs adjustment |
| Licensing for `download_source_h5ad()` at scale? | CxG ToS — check for bulk download constraints |
| Run cache job from inside Sanger network or expose Bolt externally? | Current firewall means cache job must run on Sanger infra |

---

## Key files

| File | Purpose |
|---|---|
| `test_joinid.py` | Validated 4-step join workflow |
| `census-joinid-pipeline-test.md` | Session notes and architecture decisions |
| `ROADMAP.md` | This file |
| [CellMark neo4j_client.py](https://github.com/Cellular-Semantics/CellMark/blob/main/src/scripts/neo4j_client.py) | Existing graph query patterns |
