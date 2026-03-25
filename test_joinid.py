"""
Census observation_joinid join workflow test.

Dataset:  9bb9596d-f23f-4558-912f-d4dc7d52721b
          CGE-derived interneurons integrated with 10X sequencing MOp data
          Mus musculus, 15511 cells
H5AD:     downloaded via cellxgene_census.download_source_h5ad()
"""

import anndata
import cellxgene_census
import numpy as np
import pandas as pd
import tiledbsoma as soma

H5AD_PATH = "9bb9596d-f23f-4558-912f-d4dc7d52721b.h5ad"
DATASET_ID = "9bb9596d-f23f-4558-912f-d4dc7d52721b"
ORGANISM = "mus_musculus"
CENSUS_VERSION = "stable"

# ---------------------------------------------------------------------------
# Step 1: Inspect the H5AD
# ---------------------------------------------------------------------------
print("=== Step 1: H5AD inspection ===")
adata = anndata.read_h5ad(H5AD_PATH)
print("obs columns:", adata.obs.columns.tolist())

# Confirm observation_joinid exists and inspect format
assert "observation_joinid" in adata.obs.columns, "observation_joinid missing from H5AD obs!"

# author annotation field name may vary — show all candidate columns
author_cols = [c for c in adata.obs.columns if any(t in c.lower() for t in ("author", "subclass", "biccn", "cluster_label"))]
print("Candidate author annotation columns:", author_cols)

# Known field name for this BICCN dataset; falls back to first candidate
AUTHOR_COL = "BICCN_subclass_label" if "BICCN_subclass_label" in adata.obs.columns else author_cols[0]
print(f"Using author annotation column: {AUTHOR_COL!r}")
print(adata.obs[["observation_joinid", AUTHOR_COL]].head(10))
print(f"H5AD cells: {adata.n_obs}")

# Resolve dataset_id: prefer value embedded in obs, fall back to constant
if "dataset_id" in adata.obs.columns:
    h5ad_dataset_ids = adata.obs["dataset_id"].unique().tolist()
    print(f"dataset_id(s) found in H5AD obs: {h5ad_dataset_ids}")
    DATASET_ID = h5ad_dataset_ids[0]
else:
    print(f"No dataset_id column in H5AD obs — using constant: {DATASET_ID}")
print(f"Querying Census with dataset_id: {DATASET_ID}")

# ---------------------------------------------------------------------------
# Step 2: Query Census obs for this dataset
# ---------------------------------------------------------------------------
print("\n=== Step 2: Census obs query ===")
with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    obs_df = cellxgene_census.get_obs(
        census,
        ORGANISM,
        value_filter=f"dataset_id == '{DATASET_ID}'",
        column_names=["soma_joinid", "observation_joinid", "cell_type_ontology_term_id", "dataset_id"],
    )

if len(obs_df) == 0:
    print(f"WARNING: dataset_id '{DATASET_ID}' returned 0 rows from Census.")
    print("This dataset may not be in Census (e.g. mouse-only, or not yet ingested).")
    print("Check: is this a human dataset? Is it in the CELLxGENE corpus?")
    raise SystemExit(1)

print(obs_df.head())
print(f"Census rows: {len(obs_df)}")
print("observation_joinid dtype:", obs_df["observation_joinid"].dtype)
print("observation_joinid sample:", obs_df["observation_joinid"].iloc[:3].tolist())

# ---------------------------------------------------------------------------
# Step 3: Join H5AD author annotations to Census obs
# ---------------------------------------------------------------------------
print("\n=== Step 3: Join verification ===")
h5ad_obs = adata.obs[["observation_joinid", AUTHOR_COL]].reset_index(drop=True)

# Ensure matching dtypes before join
obs_df["observation_joinid"] = obs_df["observation_joinid"].astype(str)
h5ad_obs["observation_joinid"] = h5ad_obs["observation_joinid"].astype(str)

merged = obs_df.merge(h5ad_obs, on="observation_joinid", how="inner")

print(f"Census rows:  {len(obs_df)}")
print(f"H5AD rows:    {len(h5ad_obs)}")
print(f"Merged rows:  {len(merged)}")
join_rate = len(merged) / max(len(obs_df), len(h5ad_obs))
print(f"Join rate:    {join_rate:.1%}")
assert join_rate > 0.95, f"Join rate too low: {join_rate:.1%}"

print(merged[["soma_joinid", "observation_joinid", "cell_type_ontology_term_id", AUTHOR_COL]].head(20))

# ---------------------------------------------------------------------------
# Step 4: Test coords-based Census expression query
# ---------------------------------------------------------------------------
print("\n=== Step 4: Coords expression query ===")
joinid_array = np.array(merged["soma_joinid"].values)
print(f"Querying {len(joinid_array)} cells by soma_joinid coords...")

with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    exp = census["census_data"][ORGANISM]
    with exp.axis_query(
        measurement_name="RNA",
        obs_query=soma.AxisQuery(coords=(joinid_array,)),
    ) as query:
        mini_adata = query.to_anndata(X_name="raw")

print(mini_adata)
print(f"Cells retrieved: {mini_adata.n_obs}")
assert mini_adata.n_obs == len(merged), (
    f"Cell count mismatch: got {mini_adata.n_obs}, expected {len(merged)}"
)

print("\n=== ALL CHECKS PASSED ===")
