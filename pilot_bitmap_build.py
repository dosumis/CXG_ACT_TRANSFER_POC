"""
Pilot bitmap build from export.csv graph query results.

For each Cell_cluster in the export:
  1. Identify dataset_id (from Dataset.download_link) and obs column (from cluster extra key)
  2. Download source H5AD if not cached (obs only, backed mode)
  3. Filter obs by cluster label → get observation_joinid set
  4. Join to Census obs → soma_joinid array
  5. Encode as pyroaring BitMap
  6. Report join rate and bitmap stats per cluster

Usage:
    uv run python pilot_bitmap_build.py
"""

import re
import json
import os
import anndata
import pandas as pd
import cellxgene_census
from pyroaring import BitMap

CENSUS_VERSION = "stable"
H5AD_CACHE_DIR = "."  # reuse any already-downloaded H5ADs
EXPORT_CSV = "export.csv"

# Standard Cell_cluster property keys — anything else is the obs column name
STANDARD_KEYS = {
    "sl", "qsl", "short_form", "cell_count", "label_rdfs", "label",
    "cell_type", "cell_type_fine", "curie", "uniqueFacets", "iri",
}


# ---------------------------------------------------------------------------
# Parsing helpers (Neo4j export is not valid JSON — values unquoted)
# ---------------------------------------------------------------------------

def extract_str(blob: str, key: str) -> str | None:
    m = re.search(rf'"{key}":([^\[,{{}}]+?)(?:,\"|,\}}|\}}$|$)', blob)
    if m:
        return m.group(1).strip().strip('"')
    # fallback: unquoted value up to next comma or closing brace
    m = re.search(rf'"{key}":(http[^\s,}}\"]+)', blob)
    if m:
        return m.group(1).strip()
    return None

def extract_list_first(blob: str, key: str) -> str | None:
    m = re.search(rf'"{key}":\[([^\]]+)\]', blob)
    if m:
        return m.group(1).strip().strip('"')
    return None

def extract_obs_meta(blob: str) -> list[str]:
    """Return all field_names with field_type == author_cell_type_label."""
    m = re.search(r'"obs_meta":\[(\[\[.*?\]\])\]', blob)
    if not m:
        return []
    try:
        nested = json.loads(m.group(1))
        names = []
        for group in nested:
            for entry in group:
                if entry.get("field_type") == "author_cell_type_label":
                    names.append(entry["field_name"])
        return names
    except Exception:
        return []

def extract_dataset_id(blob: str) -> str | None:
    m = re.search(r'datasets\.cellxgene\.cziscience\.com/([0-9a-f-]{36})\.h5ad', blob)
    return m.group(1) if m else None

def extract_extra_keys(blob: str) -> list[str]:
    all_keys = re.findall(r'"(\w+)":', blob)
    return [k for k in all_keys if k not in STANDARD_KEYS]


# ---------------------------------------------------------------------------
# Census dataset lookup (graph stores dataset_version_id, not dataset_id)
# ---------------------------------------------------------------------------

_census_datasets_cache: pd.DataFrame | None = None

def get_census_datasets() -> pd.DataFrame:
    global _census_datasets_cache
    if _census_datasets_cache is None:
        with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
            _census_datasets_cache = census["census_info"]["datasets"].read().concat().to_pandas()
    return _census_datasets_cache

def resolve_census_dataset_id(dataset_version_id: str, doi: str | None) -> str | None:
    """
    The graph stores dataset_version_id (UUID in download_link filename).
    Census uses a stable dataset_id. Resolve via dataset_version_id match first,
    then fall back to DOI match.
    """
    df = get_census_datasets()
    # Try dataset_version_id direct match
    match = df[df["dataset_version_id"] == dataset_version_id]
    if len(match):
        return match.iloc[0]["dataset_id"]
    # Fall back to DOI
    if doi:
        doi_clean = doi.rstrip("/")
        match = df[df["collection_doi"].str.contains(doi_clean.split("/")[-1], na=False)]
        if len(match):
            print(f"  Resolved via DOI ({doi_clean.split('/')[-1]}): {len(match)} dataset(s) — using first")
            return match.iloc[0]["dataset_id"]
    return None


# ---------------------------------------------------------------------------
# H5AD download / load
# ---------------------------------------------------------------------------

def get_h5ad_obs(census_dataset_id: str, download_url: str | None = None) -> pd.DataFrame:
    path = os.path.join(H5AD_CACHE_DIR, f"{census_dataset_id}.h5ad")
    if not os.path.exists(path):
        print(f"  Downloading {census_dataset_id}.h5ad via Census ...")
        cellxgene_census.download_source_h5ad(
            census_dataset_id, to_path=path, census_version=CENSUS_VERSION
        )
    else:
        print(f"  Using cached {census_dataset_id}.h5ad")
    adata = anndata.read_h5ad(path, backed="r")
    return adata.obs.copy()


# ---------------------------------------------------------------------------
# Census join
# ---------------------------------------------------------------------------

def get_census_obs(census_dataset_id: str) -> pd.DataFrame:
    print(f"  Querying Census obs for {census_dataset_id} ...")
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        for organism in ("homo_sapiens", "mus_musculus"):
            df = cellxgene_census.get_obs(
                census, organism,
                value_filter=f"dataset_id == '{census_dataset_id}'",
                column_names=["soma_joinid", "observation_joinid", "dataset_id"],
            )
            if len(df) > 0:
                print(f"  Found {len(df)} cells in Census ({organism})")
                return df
    print(f"  WARNING: {census_dataset_id} not found in Census")
    return pd.DataFrame(columns=["soma_joinid", "observation_joinid", "dataset_id"])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    df = pd.read_csv(EXPORT_CSV)
    df.columns = ["clus_props", "tissue_pct", "ds_props"]

    # Group by dataset to avoid re-downloading H5ADs
    datasets: dict[str, dict] = {}  # dataset_id -> {obs, census_obs, obs_columns}

    # Parse all rows first
    clusters = []
    for _, row in df.iterrows():
        cp, dp = row["clus_props"], row["ds_props"]

        label = extract_str(cp, "label")
        iri = extract_str(cp, "iri")
        cell_count = extract_list_first(cp, "cell_count")
        extra_keys = extract_extra_keys(cp)
        obs_column = extra_keys[0] if extra_keys else None  # None = curated merge

        dataset_id = extract_dataset_id(dp)
        obs_meta_fields = extract_obs_meta(dp)
        doi = extract_list_first(dp, "publication")

        clusters.append({
            "label": label,
            "iri": iri,
            "cell_count": int(cell_count) if cell_count else None,
            "obs_column": obs_column,
            "obs_meta_fields": obs_meta_fields,
            "dataset_id": dataset_id,
            "doi": doi,
        })

        if dataset_id and dataset_id not in datasets:
            datasets[dataset_id] = {
                "obs": None, "census_obs": None,
                "census_dataset_id": None,
                "doi": doi,
            }

    # Resolve Census dataset_ids and load data
    print(f"\n{'='*60}")
    print(f"Loading data for {len(datasets)} dataset(s)...")
    for dataset_version_id in datasets:
        doi = datasets[dataset_version_id]["doi"]
        print(f"\n--- Dataset version {dataset_version_id} ---")
        census_id = resolve_census_dataset_id(dataset_version_id, doi)
        if census_id is None:
            print(f"  WARNING: could not resolve Census dataset_id — skipping")
            continue
        print(f"  Census dataset_id: {census_id}")
        datasets[dataset_version_id]["census_dataset_id"] = census_id
        datasets[dataset_version_id]["obs"] = get_h5ad_obs(census_id)
        datasets[dataset_version_id]["census_obs"] = get_census_obs(census_id)

        # Build observation_joinid lookup: census obs indexed by observation_joinid
        census_obs = datasets[dataset_version_id]["census_obs"]
        if len(census_obs):
            datasets[dataset_version_id]["joinid_map"] = census_obs.set_index(
                census_obs["observation_joinid"].astype(str)
            )["soma_joinid"]
        else:
            datasets[dataset_version_id]["joinid_map"] = pd.Series(dtype="int64")

    # Build bitmaps per cluster
    print(f"\n{'='*60}")
    print("Building bitmaps...\n")
    results = []
    for clus in clusters:
        dataset_id = clus["dataset_id"]
        label = clus["label"]
        print(f"Cluster: {label!r}  (dataset {dataset_id})")

        if not dataset_id or dataset_id not in datasets or datasets[dataset_id]["census_dataset_id"] is None:
            print("  SKIP: dataset not resolved in Census\n")
            continue

        h5ad_obs = datasets[dataset_id]["obs"]
        joinid_map = datasets[dataset_id].get("joinid_map", pd.Series(dtype="int64"))

        # Determine obs column
        obs_col = clus["obs_column"]
        if obs_col is None:
            # Curated merge: try each obs_meta field
            for field in clus["obs_meta_fields"]:
                if field in h5ad_obs.columns and label in h5ad_obs[field].values:
                    obs_col = field
                    print(f"  obs_column resolved via obs_meta fallback: {obs_col!r}")
                    break
        if obs_col is None or obs_col not in h5ad_obs.columns:
            print(f"  WARNING: could not resolve obs column. obs_meta fields: {clus['obs_meta_fields']}\n")
            continue

        # Filter H5AD obs to this cluster
        mask = h5ad_obs[obs_col].astype(str) == str(label)
        cluster_obs = h5ad_obs[mask]
        h5ad_joinids = cluster_obs["observation_joinid"].astype(str)

        # Join to Census
        matched = h5ad_joinids[h5ad_joinids.isin(joinid_map.index)]
        soma_joinids = joinid_map.loc[matched].values.tolist()

        join_rate = len(soma_joinids) / len(h5ad_joinids) if len(h5ad_joinids) else 0
        print(f"  obs_column   : {obs_col}")
        print(f"  H5AD cells   : {len(h5ad_joinids)}")
        print(f"  Matched      : {len(soma_joinids)}")
        print(f"  Join rate    : {join_rate:.1%}")

        if len(soma_joinids) == 0:
            print("  WARNING: empty soma_joinid set — cluster not in Census?\n")
            continue

        bitmap = BitMap(soma_joinids)
        print(f"  Bitmap size  : {len(bitmap)} cells, {len(bitmap.serialize())} bytes serialised")
        print(f"  IRI          : {clus['iri']}")
        print(f"  DOI          : {clus['doi']}")
        print()

        results.append({
            "label": label,
            "iri": clus["iri"],
            "obs_column": obs_col,
            "dataset_id": dataset_id,
            "doi": clus["doi"],
            "h5ad_cells": len(h5ad_joinids),
            "matched_cells": len(soma_joinids),
            "join_rate": join_rate,
            "bitmap": bitmap,
        })

    # Summary
    print(f"{'='*60}")
    print(f"Summary: {len(results)}/{len(clusters)} clusters built successfully\n")
    summary = pd.DataFrame([{k: v for k, v in r.items() if k != "bitmap"} for r in results])
    print(summary.to_string(index=False))

    # Persist bitmaps keyed by IRI UUID + census version
    bitmap_dir = os.path.join(H5AD_CACHE_DIR, "bitmaps")
    os.makedirs(bitmap_dir, exist_ok=True)
    print(f"\nPersisting bitmaps to {bitmap_dir}/")
    for r in results:
        iri = r["iri"] or ""
        uuid = iri.split("/")[-1]  # extract UUID from http://example.org/{uuid}
        filename = f"{uuid}__{CENSUS_VERSION}.bitmap"
        path = os.path.join(bitmap_dir, filename)
        with open(path, "wb") as f:
            f.write(r["bitmap"].serialize())
        print(f"  {filename}  ({r['matched_cells']} cells, {os.path.getsize(path)} bytes)  [{r['label']}]")

    # Demo: bitmap union across all clusters (simulates a graph query result)
    if results:
        union_bm = BitMap()
        for r in results:
            union_bm |= r["bitmap"]
        print(f"\nUnion bitmap (all {len(results)} clusters): {len(union_bm)} unique cells")


if __name__ == "__main__":
    main()
