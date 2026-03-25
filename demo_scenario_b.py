"""
Scenario B demo: given a graph query result (export.csv), use bitmaps to
retrieve the intersection of author-annotated cells with a Census obs filter,
returning each row with provenance (DOI, cluster label, obs_column).

Pipeline:
  1. Build bitmaps from export.csv (reuses pilot_bitmap_build logic)
  2. Union all bitmaps → full candidate soma_joinid set
  3. Query Census obs with coords=union_bitmap + obs value_filter
  4. Join result back to per-cluster bitmaps → per-cluster cell counts
  5. Print result table with DOI

Usage:
    uv run python demo_scenario_b.py
"""

import re
import json
import os
import numpy as np
import anndata
import pandas as pd
import cellxgene_census
import tiledbsoma as soma
from pyroaring import BitMap

# ── Config ────────────────────────────────────────────────────────────────────

CENSUS_VERSION = "stable"
EXPORT_CSV     = "export.csv"
H5AD_CACHE_DIR = "."

# The Census obs filter that simulates a user query.
# Change this to explore different slices. Examples:
#   "disease_ontology_term_id == 'MONDO:0100096'"          # COVID-19
#   "tissue_ontology_term_id == 'UBERON:0002048'"           # lung
#   "assay_ontology_term_id == '10x 3\\' v3'"               # 10x Chromium
#   "sex_ontology_term_id == 'PATO:0000384'"                # male
OBS_FILTER = "disease_ontology_term_id == 'MONDO:0100096'"  # COVID-19

# Results are printed to stdout:
#   - Per-cluster table: cells_in_cluster, cells_passing_filter, fraction, doi
#   - Sample rows: soma_joinid, cell_type, disease, tissue, cluster_label, doi
# Re-run with a different OBS_FILTER to explore other Census obs dimensions.

STANDARD_KEYS = {
    "sl","qsl","short_form","cell_count","label_rdfs","label",
    "cell_type","cell_type_fine","curie","uniqueFacets","iri",
}

# ── Parsing (Neo4j export is not valid JSON) ───────────────────────────────────

def extract_str(blob, key):
    m = re.search(rf'"{key}":([^\[,{{}}]+?)(?:,\"|,\}}|\}}$|$)', blob)
    if m:
        return m.group(1).strip().strip('"')
    m = re.search(rf'"{key}":(http[^\s,}}\"]+)', blob)
    return m.group(1).strip() if m else None

def extract_list_first(blob, key):
    m = re.search(rf'"{key}":\[([^\]]+)\]', blob)
    return m.group(1).strip().strip('"') if m else None

def extract_obs_meta(blob):
    m = re.search(r'"obs_meta":\[(\[\[.*?\]\])\]', blob)
    if not m:
        return []
    try:
        nested = json.loads(m.group(1))
        return [e["field_name"] for group in nested for e in group
                if e.get("field_type") == "author_cell_type_label"]
    except Exception:
        return []

def extract_dataset_version_id(blob):
    m = re.search(r'datasets\.cellxgene\.cziscience\.com/([0-9a-f-]{36})\.h5ad', blob)
    return m.group(1) if m else None

def extract_extra_keys(blob):
    return [k for k in re.findall(r'"(\w+)":', blob) if k not in STANDARD_KEYS]

# ── Census helpers ─────────────────────────────────────────────────────────────

_census_datasets: pd.DataFrame | None = None

def census_datasets():
    global _census_datasets
    if _census_datasets is None:
        with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
            _census_datasets = census["census_info"]["datasets"].read().concat().to_pandas()
    return _census_datasets

def resolve_census_id(version_id, doi):
    df = census_datasets()
    m = df[df["dataset_version_id"] == version_id]
    if len(m):
        return m.iloc[0]["dataset_id"]
    if doi:
        m = df[df["collection_doi"].str.contains(doi.split("/")[-1], na=False)]
        if len(m):
            return m.iloc[0]["dataset_id"]
    return None

def get_h5ad_obs(census_id):
    path = os.path.join(H5AD_CACHE_DIR, f"{census_id}.h5ad")
    if not os.path.exists(path):
        print(f"  Downloading {census_id}.h5ad ...")
        cellxgene_census.download_source_h5ad(census_id, to_path=path,
                                               census_version=CENSUS_VERSION)
    return anndata.read_h5ad(path, backed="r").obs.copy()

def get_census_joinid_map(census_id):
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        for organism in ("homo_sapiens", "mus_musculus"):
            df = cellxgene_census.get_obs(
                census, organism,
                value_filter=f"dataset_id == '{census_id}'",
                column_names=["soma_joinid", "observation_joinid"],
            )
            if len(df):
                return df.set_index(df["observation_joinid"].astype(str))["soma_joinid"]
    return pd.Series(dtype="int64")

# ── Bitmap build (condensed from pilot) ───────────────────────────────────────

def build_bitmaps(clusters, datasets):
    bitmaps = {}
    for clus in clusters:
        vid = clus["dataset_version_id"]
        if vid not in datasets or datasets[vid].get("census_id") is None:
            continue
        h5ad_obs  = datasets[vid]["obs"]
        joinid_map = datasets[vid]["joinid_map"]

        obs_col = clus["obs_column"]
        if obs_col is None:
            for field in clus["obs_meta_fields"]:
                if field in h5ad_obs.columns and clus["label"] in h5ad_obs[field].values:
                    obs_col = field
                    break
        if obs_col is None or obs_col not in h5ad_obs.columns:
            continue

        mask = h5ad_obs[obs_col].astype(str) == str(clus["label"])
        h5ad_joinids = h5ad_obs[mask]["observation_joinid"].astype(str)
        matched = h5ad_joinids[h5ad_joinids.isin(joinid_map.index)]
        soma_ids = joinid_map.loc[matched].values.tolist()

        if soma_ids:
            bitmaps[clus["iri"]] = {
                "bitmap":     BitMap(soma_ids),
                "label":      clus["label"],
                "obs_column": obs_col,
                "doi":        clus["doi"],
                "dataset_version_id": vid,
            }
    return bitmaps

# ── Scenario B query ───────────────────────────────────────────────────────────

def scenario_b_query(bitmaps, obs_filter):
    # Union all bitmaps
    union_bm = BitMap()
    for v in bitmaps.values():
        union_bm |= v["bitmap"]
    coords = np.array(sorted(union_bm), dtype=np.int64)
    print(f"Union bitmap : {len(union_bm):,} cells")
    print(f"Obs filter   : {obs_filter}")
    print(f"Querying Census obs (coords + filter)...\n")

    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        exp = census["census_data"]["homo_sapiens"]
        obs = exp.axis_query(
            measurement_name="RNA",
            obs_query=soma.AxisQuery(
                coords=(coords,),
                value_filter=obs_filter,
            ),
        ).obs(column_names=["soma_joinid", "observation_joinid",
                             "cell_type_ontology_term_id", "cell_type",
                             "disease", "tissue", "dataset_id"]).concat().to_pandas()

    print(f"Cells passing filter: {len(obs):,}\n")
    if len(obs) == 0:
        print("No cells matched the filter.")
        return

    # Per-cluster breakdown
    filter_set = BitMap(obs["soma_joinid"].values.tolist())
    rows = []
    for iri, v in bitmaps.items():
        intersection = v["bitmap"] & filter_set
        if len(intersection) == 0:
            continue
        rows.append({
            "cluster_label": v["label"],
            "obs_column":    v["obs_column"],
            "cells_in_cluster": len(v["bitmap"]),
            "cells_passing_filter": len(intersection),
            "fraction": len(intersection) / len(v["bitmap"]),
            "doi": v["doi"],
            "iri": iri,
        })

    result = pd.DataFrame(rows).sort_values("cells_passing_filter", ascending=False)

    print("=== Per-cluster breakdown ===")
    display_cols = ["cluster_label","obs_column","cells_in_cluster",
                    "cells_passing_filter","fraction","doi"]
    print(result[display_cols].to_string(index=False))

    print(f"\nTotal unique cells after filter : {len(filter_set & union_bm):,}")
    print(f"Clusters with ≥1 passing cell   : {len(result)}/{len(bitmaps)}")

    # Sample rows with full provenance
    print("\n=== Sample cells (first 10, with provenance) ===")
    sample_ids = list(sorted(filter_set & union_bm))[:10]
    sample_obs = obs[obs["soma_joinid"].isin(sample_ids)].copy()

    # attach cluster label + doi by soma_joinid lookup
    soma_to_meta = {}
    for iri, v in bitmaps.items():
        for sid in v["bitmap"]:
            if sid in filter_set:
                soma_to_meta[sid] = (v["label"], v["doi"])

    sample_obs["cluster_label"] = sample_obs["soma_joinid"].map(
        lambda s: soma_to_meta.get(s, ("?", "?"))[0])
    sample_obs["doi"] = sample_obs["soma_joinid"].map(
        lambda s: soma_to_meta.get(s, ("?", "?"))[1])

    print(sample_obs[["soma_joinid","cell_type","disease","tissue",
                       "cluster_label","doi"]].to_string(index=False))


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    df = pd.read_csv(EXPORT_CSV)
    df.columns = ["clus_props", "tissue_pct", "ds_props"]

    clusters = []
    datasets = {}
    for _, row in df.iterrows():
        cp, dp = row["clus_props"], row["ds_props"]
        vid  = extract_dataset_version_id(dp)
        doi  = extract_list_first(dp, "publication")
        extra = extract_extra_keys(cp)
        clusters.append({
            "label":           extract_str(cp, "label"),
            "iri":             extract_str(cp, "iri"),
            "obs_column":      extra[0] if extra else None,
            "obs_meta_fields": extract_obs_meta(dp),
            "dataset_version_id": vid,
            "doi":             doi,
        })
        if vid and vid not in datasets:
            datasets[vid] = {"census_id": None, "obs": None, "joinid_map": None, "doi": doi}

    print("Resolving Census dataset IDs and loading H5AD obs...\n")
    for vid, d in datasets.items():
        census_id = resolve_census_id(vid, d["doi"])
        if census_id is None:
            print(f"  SKIP {vid}: not found in Census")
            continue
        print(f"  {vid} → {census_id}")
        d["census_id"] = census_id
        d["obs"]        = get_h5ad_obs(census_id)
        d["joinid_map"] = get_census_joinid_map(census_id)

    print("\nBuilding bitmaps...")
    bitmaps = build_bitmaps(clusters, datasets)
    print(f"Built {len(bitmaps)} bitmaps.\n")

    print("=" * 60)
    scenario_b_query(bitmaps, OBS_FILTER)


if __name__ == "__main__":
    main()
