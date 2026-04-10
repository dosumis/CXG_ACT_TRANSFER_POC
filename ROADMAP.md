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
(Cell_cluster {
  label,                           // existing — curated/winning display label
  iri,                             // existing
  author_label_column: str,        // NEW — obs column name that sourced the winning label
  author_synonym_columns: [str],   // NEW — other obs column names with identical cell sets
                                   //       labels recoverable from existing node properties
                                   //       keyed by those column names
  // author obs column values already stored as ad-hoc properties e.g. celltype: "Macro_c2-CCL3L1"
})
  -[:composed_primarily_of]-> (Cell {curie, label})
  -[:subcluster_of*0..]-> (Cell_cluster)
  -[:tissue]-> (Class {label})
  -[:has_source]-> (Dataset {title, census_dataset_id: str, download_link: str,
                              publication: str, obs_meta: str,
                              census_version_cached: str})
```

| New property | Node | Notes |
|---|---|---|
| `author_label_column` | `Cell_cluster` | Obs column name that sourced the winning/curated label. Fixes the merged-node gap where no column key is currently stored. Used by bitmap build job to filter H5AD obs. |
| `author_synonym_columns` | `Cell_cluster` | List of other obs column names whose cell sets were identical and merged into this node. Labels recoverable from existing ad-hoc column-keyed properties on the node. Neo4j stores as `[str]` — nested maps not supported as property values. |
| `census_dataset_id` | `Dataset` | Resolved Census `dataset_id` (differs from `dataset_version_id` in `download_link`) |
| `census_version_cached` | `Dataset` | Census build against which current bitmaps were generated |

**Note on `observation_joinid`:** this is an H5AD/Census concept only — it never touches the KG. It is used transiently during the bitmap build job as the join key between H5AD obs and Census obs, then discarded. The KG only needs `author_label_column` to know how to extract the right cells from an H5AD.

**Bitmap store:** keyed by `(node_iri, census_version)` — no separate UUID needed.

### Reconstructing the full author annotation picture from KG alone

Given a `Cell_cluster` node:
1. `label` — the curated display label
2. `author_label_column` → look up `node[author_label_column]` → original author label string (if stored)
3. `author_synonym_columns` → for each column name, look up `node[column_name]` → synonym author labels

This gives the complete set of author annotation strings across all merged obs columns, recoverable without touching the H5AD.

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

## Primary use case

> *"Given a Census obs filter, which source datasets and author labels are relevant? Return matching cells with provenance."*

1. Apply obs filter on graph → find `Cell_cluster` nodes (via `subcluster_of` subtree)
2. Decompress bitmaps for those nodes → `soma_joinid` arrays
3. Query Census obs (no expression) with coords + obs filter → get intersection
4. Return result rows enriched with `collection_doi` from `Dataset` node
5. Optionally: pass intersection `soma_joinid` set to `axis_query` with gene filter → expression slice

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

---

## Implementation plan — Pandasaurus_cxG tickets

The PoC (Phases 0–2b above) is complete and validates the end-to-end workflow. The next stage is integrating this into the Pandasaurus_cxG build pipeline. Three tickets, in dependency order.

---

### Ticket 1 — Schema: add `author_label_column`, `author_synonym_columns`, `census_dataset_id`

**Goal:** Update the `Cell_cluster` and `Dataset` node schemas so each node explicitly records which obs column sourced its label, which other columns were merged, and the resolved Census dataset identifier.

**Changes to `Cell_cluster` nodes:**
- `author_label_column: str` — the obs column name that sourced this node's winning label (e.g. `"cell_type_fine"`)
- `author_synonym_columns: [str]` — other obs column names whose identical cell sets were merged into this node. Labels recoverable from existing ad-hoc properties keyed by those column names.

**Changes to `Dataset` nodes:**
- `census_dataset_id: str` — resolved Census `dataset_id` (differs from `dataset_version_id` in `download_link`; resolution via version_id or DOI match, both confirmed working in PoC)
- `census_version_cached: str` — Census build against which bitmaps were last generated

**Where:** `GraphGenerator.generate_rdf_graph` in Pandasaurus_cxG — the method already iterates `(field_name, value)` pairs per `Cell_cluster`. `author_label_column` is just the `field_name` that won curation. `author_synonym_columns` requires detecting identical cell sets during the merge step (`AnndataAnalyzer`) and recording which fields produced them.

**Pre-requisites:** None — this is a schema-only change with no external dependencies.

---

### Ticket 2 — Bitmap build integrated into graph build pipeline

**Goal:** At graph build time, while the H5AD is still open, build `soma_joinid` roaring bitmaps for each `Cell_cluster` node.

**What it does:** For each `Cell_cluster` node during build:
1. Filter `obs[author_label_column] == label` → extract `observation_joinid` set (H5AD already open)
2. Join to Census obs for the dataset's `census_dataset_id` → resolve `soma_joinid` array
3. Encode as roaring bitmap, persist keyed by `(node_iri, census_version)`

`observation_joinid` is transient — used only as the join key, never stored in the KG.

**Storage:** File-based initially (`{iri_uuid}__{census_version}.bitmap`), co-located with the query server. Redis or node-property storage if latency requirements tighten.

**Why co-schedule with graph build:** The H5AD is already downloaded and open. Extracting `observation_joinid` sets is a cheap obs-only operation. Doing it at build time avoids a second H5AD pass entirely.

**Prototype:** `pilot_bitmap_build.py` — generalise to operate on live graph rather than `export.csv`.

**Pre-requisites:** Ticket 1 (nodes have `author_label_column` and `census_dataset_id`).

---

### Ticket 3 — Bitmap cache refresh job (per Census release)

**Goal:** When a new Census build is released (but the graph hasn't changed), rebuild bitmaps without re-downloading H5ADs — using `author_label_column` + cached H5ADs or re-downloading only obs.

**What it does:**
1. For each `Cell_cluster` node: re-run the `observation_joinid` → `soma_joinid` join against the new Census build
2. Overwrite bitmap files with new Census version key
3. Update `census_version_cached` on `Dataset` nodes
4. Log join rate per node; flag drops >5% vs prior build (signals `observation_joinid` instability)

**Cadence:** Runs once per Census release, independent of graph build cadence.

**Note:** This does need the H5ADs (for `observation_joinid` extraction). Options:
- Keep H5AD cache from last graph build (simplest — just don't delete them)
- Re-download via `download_source_h5ad()` (safe fallback but slow)
- Store `observation_joinid` sets in a lightweight sidecar file at Ticket 2 build time, avoiding H5AD re-reads entirely

**Pre-requisites:** Ticket 2 complete (bitmap build proven). Must run on Sanger infra or via VPN (Census access is public; Bolt access is not).

---

### Server architecture — query endpoint

The bitmap query workflow (graph → bitmaps → Census → results) needs a runtime home. Options:

| Option | Notes |
|---|---|
| **New endpoint on same server** | Lightweight FastAPI service alongside Neo4j. Accepts a CL term or Cell_cluster IRI, looks up cached bitmaps, queries Census, returns obs/expression slice. |
| Extend `cxg_query_enhance` | If it already has a Census query interface, add a CL-KG-aware routing layer that resolves CL terms → bitmaps before querying. |
| Embedded in client | For PoC/dev: client-side script (no server needed). Not suitable for production. |

**Recommended:** New FastAPI endpoint on the same host as `cl-kg-neo4j-db`. Supplements the existing CL-KG Neo4j interface without changing it. Extensions to `cxg_query_enhance` come after the endpoint is stable — that library becomes the client-side interface.

---

## Open questions

| Question | Status | Notes |
|---|---|---|
| Is `observation_joinid` stable across Census builds? | **Open** | Critical — verify across two Census builds before Ticket B goes to production |
| Where does the bitmap cache live? | **Decided (draft)** | File store initially, keyed by `(iri_hash, census_version)`; see server architecture above |
| Does `Cell_cluster` always map 1:1 to a single `dataset_id`? | **Open** | If a cluster spans datasets, observation_joinid join logic needs adjustment |
| Licensing for `download_source_h5ad()` at scale? | **Open** | CxG ToS — check for bulk download constraints before running across full graph |
| Run cache job from inside Sanger network or expose Bolt externally? | **Confirmed** | Must run on Sanger infra (Bolt is internal-only); Census access is public |

---

## Key files

| File | Purpose |
|---|---|
| `test_joinid.py` | Validated 4-step join workflow |
| `pilot_bitmap_build.py` | Bitmap build for pilot clusters from `export.csv` |
| `demo_scenario_b.py` | End-to-end query demo — edit `OBS_FILTER` to try other filters |
| `census-joinid-pipeline-test.md` | Session notes and architecture decisions |
| `ROADMAP.md` | This file |
| [CellMark neo4j_client.py](https://github.com/Cellular-Semantics/CellMark/blob/main/src/scripts/neo4j_client.py) | Existing graph query patterns |
