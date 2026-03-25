---
date: 2026-03-18
updated: 2026-03-19
status: test complete — architecture validated
tags:
  - cl-kg
  - census
  - cellxgene
  - pipeline
  - observation-joinid
  - task
---

# CL-KG Pipeline: Census Join ID Test

**Goal:** Validate the `observation_joinid` join workflow on a small real dataset before integrating into the CL-KG graph construction pipeline.

**Test dataset:**
- **Dataset ID:** `8f98c236-43f0-4dc4-985b-c304499f7b44`
- **Collection ID:** `20a1dadf-a3a7-4783-b311-fcff3c457763`
- Small patch-seq dataset — downloaded locally as H5AD

---

## Background: Architecture Decision

### The problem being solved

The CL-KG pipeline needs to store cell-level evidence linking:
- **CL-level annotations** (`cell_type_ontology_term_id`) — from Census obs
- **Author-level annotations** (`author_cell_type`) — from source H5AD only (not in Census)
- **Expression data** — accessible via Census soma_joinid coordinates

The key question was: what cell identifier should graph nodes store, and how do we join author annotations to Census data?

### Cell identifiers in Census

Two relevant IDs in Census obs:

| Field | Description | Stability |
|---|---|---|
| `soma_joinid` | SOMA integer key; indexes X matrix; NOT positionally indexed | **Reassigned each Census build** |
| `observation_joinid` | Inherited from CELLxGENE dataset schema; also present in source H5AD | Likely stable across builds (to verify) |

> [!note] Key caveat
> `observation_joinid` stability across Census builds is not explicitly documented. If it's stable, it's the linchpin identifier. If not, pinning to an LTS build and using `soma_joinid` is the safe fallback.

### Architecture: Efficient Census access via coords

TileDB-SOMA is optimised for integer coordinate lookups. You can query Census directly by `soma_joinid` without a value_filter scan:

```python
import tiledbsoma as soma

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    exp = census["census_data"]["homo_sapiens"]
    with exp.axis_query(
        measurement_name="RNA",
        obs_query=soma.AxisQuery(coords=(my_joinid_array,)),
    ) as query:
        adata = query.to_anndata(X_name="raw")
```

- `coords` lookup → fast (direct sparse array indexing)
- `value_filter` on `observation_joinid` → slow (string scan over millions of rows)

This means: store `soma_joinid` sets as roaring bitmaps on graph nodes → at query time decompress bitmap → pass directly to Census as coords. Efficient for both graph construction and runtime queries.

---

## Graph Construction Workflow

### One-time per dataset (graph build)

1. **Download source H5AD** for dataset of interest
2. **Extract** `observation_joinid → author_cell_type` mapping from H5AD obs
3. **Query Census obs** for same `dataset_id`, pulling only:
   - `soma_joinid`
   - `observation_joinid`
   - `cell_type_ontology_term_id`
4. **Join** on `observation_joinid` → get triple: `(soma_joinid, CL_term, author_cell_type)`
5. **Build graph nodes** at both CL-level and author-annotation-level
6. **Store** `soma_joinid` sets as roaring bitmaps on each node

### Census obs query (no expression data needed):

```python
import cellxgene_census

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    obs_df = cellxgene_census.get_obs(
        census,
        "homo_sapiens",
        value_filter="dataset_id == '8f98c236-43f0-4dc4-985b-c304499f7b44'",
        column_names=["soma_joinid", "observation_joinid", "cell_type_ontology_term_id"],
    )
```

### Runtime (after graph built)

- Decompress bitmap from graph node → array of `soma_joinid` values
- Pass as `coords` to Census → get back expression matrix for exactly those cells
- No H5AD download needed at runtime

---

## Test Plan

### Step 1: Inspect the H5AD

```python
import anndata
adata = anndata.read_h5ad("path/to/patchseq.h5ad")
print(adata.obs.columns.tolist())
print(adata.obs[["observation_joinid", "author_cell_type"]].head(10))
```

Check: does `observation_joinid` exist in obs? What does it look like (string? int?)?

### Step 2: Query Census obs for dataset

```python
import cellxgene_census

with cellxgene_census.open_soma(census_version="stable") as census:
    obs_df = cellxgene_census.get_obs(
        census,
        "homo_sapiens",
        value_filter="dataset_id == '8f98c236-43f0-4dc4-985b-c304499f7b44'",
        column_names=["soma_joinid", "observation_joinid", "cell_type_ontology_term_id"],
    )
print(obs_df.head())
print(f"Rows: {len(obs_df)}")
```

### Step 3: Verify join

```python
# Join H5AD author annotations to Census obs
h5ad_obs = adata.obs[["observation_joinid", "author_cell_type"]].reset_index()
merged = obs_df.merge(h5ad_obs, on="observation_joinid", how="inner")

print(f"Census rows: {len(obs_df)}")
print(f"H5AD rows: {len(h5ad_obs)}")
print(f"Merged rows: {len(merged)}")
print(merged[["soma_joinid", "observation_joinid", "cell_type_ontology_term_id", "author_cell_type"]].head(20))
```

Key checks:
- [ ] Join rate close to 100% (all cells match)
- [ ] `observation_joinid` format is consistent between sources
- [ ] `author_cell_type` populated correctly on merged rows

### Step 4: Test coords query

```python
import tiledbsoma as soma
import numpy as np

joinid_array = np.array(merged["soma_joinid"].values)

with cellxgene_census.open_soma(census_version="stable") as census:
    exp = census["census_data"]["homo_sapiens"]
    with exp.axis_query(
        measurement_name="RNA",
        obs_query=soma.AxisQuery(coords=(joinid_array,)),
    ) as query:
        mini_adata = query.to_anndata(X_name="raw")

print(mini_adata)
print(f"Cells retrieved: {mini_adata.n_obs}")
```

Check: cell count matches expectation, expression matrix present.

---

## Dependencies

```bash
uv add cellxgene-census tiledbsoma anndata pandas pyroaring
```

---

## Open Questions

1. **`observation_joinid` stability** — is it stable across Census builds? Check `#cellxgene-census-users` Slack or test across two builds.
2. **`author_cell_type` field name** — may vary by dataset. Inspect H5AD obs columns to confirm exact field name.
3. **LTS version to pin to** — `2025-11-08` is current LTS candidate; confirm latest stable via `cellxgene_census.get_census_version_description("stable")`.

---

## Next Steps

- [ ] Run test on patch-seq dataset (Steps 1–4 above)
- [ ] Confirm join rate and `observation_joinid` consistency
- [ ] If successful: design integration into CL-KG graph construction pipeline
- [ ] Consult `#cellxgene-census-users` on `observation_joinid` stability guarantee
- [ ] Consider whether to store `observation_joinid` or just `soma_joinid` + LTS pin on graph nodes

**Feeds into:** [[bican-scellector-rebuild-planning]], [[cl-knowledge-base-white-paper]]

---

## Session Summary 2026-03-19

### Test outcome

Ran all 4 steps successfully against Census stable (2025-11-08).

**Test dataset used:** `9bb9596d-f23f-4558-912f-d4dc7d52721b`
CGE-derived interneurons integrated with 10X sequencing MOp data (*Mus musculus*, MOp)
Downloaded via `cellxgene_census.download_source_h5ad()`.

> [!note] Original dataset_id in this doc (`8f98c236...`) does not exist in Census.
> The collection (`20a1dadf...`) contains 3 Census datasets (CGE interneurons, MGE interneurons, excitatory neurons); the patch-seq H5AD (`89f305e5...`) is mouse-only and was not ingested into Census.

**Results:**
- Join rate: **98.5%** (15,511 / 15,755 cells matched) — gap is expected QC exclusions
- `observation_joinid` format: short random strings (e.g. `DGKA}N^0-J`), consistent between H5AD and Census
- Author annotation field for this dataset: **`BICCN_subclass_label`** (not `author_cell_type`)
- Coords query (Step 4) returned correct cell count with expression matrix (53,384 genes)
- Unique key is `(dataset_id, observation_joinid)` — `observation_joinid` is per-dataset, not globally unique

**Test script:** `test_joinid.py`

---

### Architecture decisions confirmed

#### What graph nodes store
- **`observation_joinid` sets** — stable across Census builds (assumed; to verify)
- **`BICCN_subclass_label` / author annotation** — stored on node as label
- **CL term** (`cell_type_ontology_term_id`) — stored on node
- **DOI / citation** — from `census['census_info']['datasets']` (`collection_doi`), captured at graph-build time

#### Post-build cache job (per Census release)
For each annotation node:
1. Decompress `observation_joinid` set
2. Query Census obs only (no expression) → get `soma_joinid` mapping
3. Store roaring bitmap keyed by `(node_uuid, census_version)`

This decouples three cadences:
- **Graph build** — driven by new datasets/annotations
- **Census release** — drives cache refresh job
- **Query time** — always hits cached `soma_joinid` bitmaps → fast coords lookup

Cache job is also a build-over-build consistency check: a drop in join rate for a node flags `observation_joinid` instability.

#### Query-time workflow
- Decompress cached soma_joinid bitmap → array
- Pass as `coords=` to Census `axis_query` (fast sparse index lookup)
- Optionally add `var_query` to restrict genes — avoids pulling full 53k gene matrix

#### Use case scenarios

**Scenario A — annotation update via Census filter:**
- Query Census obs filter → get soma_joinid set
- Intersect (inner join) with graph node bitmaps → identify which author annotation nodes are subsumed
- Update CL annotations on graph

**Scenario B — author annotation discovery (more imminent):**
- Apply obs filter on graph → find annotation nodes with partial cell overlap
- Retrieve soma_joinid bitmaps for those nodes
- Query Census obs (not expression) with coords + obs filter → get exact intersection
- Use intersection soma_joinids for expression retrieval or H5AD join
- Each result row carries DOI for traceability to source paper

> [!note] Both scenarios use inner join (not outer). Outer join (graph as left) is useful for QA — surfaces cells in author annotation absent from Census — but is not the primary use case.

#### Step 4 performance note
Retrieving expression for 15k cells × 53k genes is slow due to data volume, not indexing. In practice: filter vars to a gene set of interest. The coords path is correct; ceiling is network bandwidth.

---

### Next steps

- [ ] Verify `observation_joinid` stability across two Census builds (query `#cellxgene-census-users`)
- [ ] Pilot: build minimal graph with 2–3 annotation nodes storing `observation_joinid` sets
- [ ] Pilot: implement cache job — join `observation_joinid` → `soma_joinid` per Census build
- [ ] Pilot: demonstrate end-to-end query (CL node → bitmap → Census coords → expression slice)
- [ ] Confirm latest LTS Census version to pin to
