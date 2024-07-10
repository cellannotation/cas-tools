import json
from typing import Optional

import anndata

from cas.file_utils import read_anndata_file, read_json_file


def populate_cell_ids(cas_json_path: str, anndata_path: str, labelsets: list = None):
    """
    Add/update CellIDs to CAS from matching AnnData file.

    Parameters:
        cas_json_path: path to CAS JSON file
        anndata_path: path to anndata file
        labelsets: List of labelsets to update with IDs from AnnData. If value is null, rank '0' labelset is used.
    """
    ad = read_anndata_file(anndata_path)
    if ad is not None:
        cas = read_json_file(cas_json_path)
        cas = add_cell_ids(cas, ad, labelsets)
        if cas:
            with open(cas_json_path, "w") as json_file:
                json.dump(cas, json_file, indent=2)
    else:
        raise Exception("Anndata read operation failed: {}".format(anndata_path))


def add_cell_ids(cas: dict, ad: Optional[anndata.AnnData], labelsets: list = None):
    """
    Add/update CellIDs to CAS from matching AnnData file.

    Parameters:
        cas: CAS JSON object
        ad: anndata object
        labelsets: List of labelsets to update with IDs from AnnData. If value is null, rank '0' labelset is used. The
        labelsets should be provided in order, starting from rank 0 (leaf nodes) and ascending to higher ranks.
    """
    rank_zero_labelset = [lbl_set["name"] for lbl_set in cas["labelsets"]
                          if isinstance(cas["labelsets"][0].get("rank"), int) and lbl_set.get("rank") == 0
                          or lbl_set.get("rank") == "0"][0]
    if not labelsets:
        labelsets = rank_zero_labelset

    cluster_identifier_column = get_obs_cluster_identifier_column(ad, labelsets, rank_zero_labelset)

    if cluster_identifier_column:
        cid_lookup = {}
        for anno in cas["annotations"]:
            if anno["labelset"] == rank_zero_labelset and anno["labelset"] in labelsets:
                cell_ids = []
                if cluster_identifier_column.endswith("_id") or cluster_identifier_column.lower().endswith(" id"):
                    # cluster column value is integer cluster id
                    cluster_id = anno.get("author_annotation_fields", {}).get(cluster_identifier_column)
                    if not cluster_id:
                        raise ValueError("AnnData cluster identifier column ({}) couldn't be find found in"
                                         "CAS author_annotation_fields.".format(cluster_identifier_column))
                    cell_ids = list(
                        ad.obs.loc[
                            ad.obs[cluster_identifier_column] == int(cluster_id),
                            cluster_identifier_column,
                        ].index
                    )
                else:
                    # cluster column value is cell label
                    cluster_label = anno["cell_label"]
                    cell_ids = list(
                        ad.obs.loc[
                            ad.obs[cluster_identifier_column] == cluster_label
                            ].index
                    )
                anno["cell_ids"] = cell_ids
                if "parent_cell_set_name" in anno:
                    lookup_key = anno["parent_cell_set_name"]
                    if lookup_key in cid_lookup:
                        cid_lookup[lookup_key].update(cell_ids)
                    else:
                        cid_lookup[lookup_key] = set(cell_ids)

        for anno in cas["annotations"]:
            if anno["labelset"] in labelsets and anno["cell_label"] in cid_lookup:
                cell_ids = list(cid_lookup.get(anno["cell_label"], []))
                anno["cell_ids"] = cell_ids

        return cas
    else:
        print(
            "WARN: Cluster identifier column couldn't be identified in OBS. Populate cell ids operation is aborted."
        )

    return None


def get_obs_cluster_identifier_column(ad, labelsets: list = None, rank_zero_labelset: str = None):
    """
    Anndata files may use different column names to uniquely identify Clusters. Get the cluster identifier column name for the current file.
    Args:
        ad: anndata object
        labelsets: List of labelsets to update with IDs from AnnData. The labelsets should be provided in order,
        starting from rank 0 (leaf nodes) and ascending to higher ranks.
        rank_zero_labelset: rank 0 labelset name
    Returns:
        cluster identifier column name
    """
    obs_keys = ad.obs_keys()
    cluster_identifier_column = ""
    if labelsets and labelsets[0] in obs_keys:
        cluster_identifier_column = labelsets[0]
    else:
        # check id variants
        if rank_zero_labelset.title() + "_id" in obs_keys:
            cluster_identifier_column = rank_zero_labelset.title() + "_id"
        elif rank_zero_labelset.lower() + "_id" in obs_keys:
            cluster_identifier_column = rank_zero_labelset.lower() + "_id"
        elif rank_zero_labelset.title() + " id" in obs_keys:
            cluster_identifier_column = rank_zero_labelset.title() + " id"
        elif rank_zero_labelset.lower() + " id" in obs_keys:
            cluster_identifier_column = rank_zero_labelset.lower() + " id"
    return cluster_identifier_column


# def populate_cell_ids(cas_json_path: str, anndata_path: str, labelsets: list = None):
#     """
#     Add/update CellIDs to CAS from matching AnnData file.
#
#     Parameters:
#         cas_json_path: path to CAS JSON file
#         anndata_path: path to anndata file
#         labelsets: List of labelsets to update with IDs from AnnData. If value is null, rank '0' labelset is used.
#     """
#     ad = read_anndata_file(anndata_path)
#     cas = read_json_file(cas_json_path)
#
#     if not labelsets:
#         labelsets = [lbl_set["name"] for lbl_set in cas["labelsets"] if lbl_set["rank"] == "0"]
#
#     print(labelsets)
#     out = {}
#     for k, v in ad.obs.iterrows():
#         for r in labelsets:
#             val = v[r]
#             if val:
#                 if not r in out.keys(): out[r] = {}
#                 if not val in out[r].keys(): out[r][val] = []
#                 out[r][v[r]].append(k)
#
#     for a in cas['annotations']:
#         if not (a['labelset'] in out.keys()):
#             print("Unknown labelset %s", (a['labelset']))
#         elif not (a['cell_label']) in out[a['labelset']].keys():
#             print("Unknown value %s in labelset %s", (a['cell_label'], a['labelset']))
#         else:
#             a['cell_ids'] = out[a['labelset']][a['cell_label']]
#
#     with open(cas_json_path, 'w') as f:
#         f.write(json.dumps(cas))
#
#     ad.file.close()
