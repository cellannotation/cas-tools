import json

from cas.file_utils import read_json_file, read_anndata_file


def populate_cell_ids(cas_json_path: str, anndata_path: str, labelsets: list = None):
    """
    Add/update CellIDs to CAS from matching AnnData file.

    Parameters:
        cas_json_path: path to CAS JSON file
        anndata_path: path to anndata file
        labelsets: List of labelsets to update with IDs from AnnData. If value is null, rank '0' labelset is used.
    """
    ad = read_anndata_file(anndata_path)
    cas = read_json_file(cas_json_path)

    if not labelsets:
        labelsets = [lbl_set["name"] for lbl_set in cas["labelsets"] if lbl_set["rank"] == "0"]

    cid_lookup = {}
    for anno in cas["annotations"]:
        if "user_annotations" in anno and "parent_cell_set_name" in anno and anno["labelset"] in labelsets:
            cluster_id = anno["user_annotations"][0]["cell_label"]
            cell_ids = list(ad.obs.loc[ad.obs["cluster_id"] == int(cluster_id), "cluster_id"].index)
            anno["cell_ids"] = cell_ids
            lookup_key = anno["parent_cell_set_name"]
            if lookup_key in cid_lookup:
                cid_lookup[lookup_key].update(cell_ids)
            else:
                cid_lookup[lookup_key] = set(cell_ids)

    for anno in cas["annotations"]:
        if anno["labelset"] in labelsets and anno["cell_label"] in cid_lookup:
            cell_ids = list(cid_lookup.get(anno["cell_label"], []))
            anno["cell_ids"] = cell_ids

    with open(cas_json_path, "w") as json_file:
        json.dump(cas, json_file, indent=2)


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
