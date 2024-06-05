"""
Script Description:
-------------------
This script integrates cell annotations from a CAS (Cell Annotation Schema) JSON file into an AnnData object.
It performs validation checks to ensure data consistency between the CAS file and the AnnData file.
The AnnData file location should ideally be specified as a resolvable path in the CAS file.

Validation Checks:
------------------
1. Check if all barcodes (cell IDs) in the CAS file exist in the AnnData file.
   - If not, a warning is issued with an option to terminate.

2. Check if all barcodes (cell IDs) in the AnnData file exist in the CAS file.
   - If not, a warning is issued with an option to terminate.

3. Check if any obs keys in the AnnData file match labelset names in the CAS.
   - If matches are found:
     - Verify if the cell sets (sets of cell IDs) associated with each annotation for the labelsets in CAS match the AnnData.
     - If they do, check if the cell_labels are identical.
       - If yes, no action is taken.
       - If no, a warning is issued with options to update values to cell_labels in CAS or terminate.
     - If no matches are found:
       - Option to flush and replace obs labelsets with those defined in CAS (labelset:cell_label pairs) or terminate.
         - For example, if a labelset in the AnnData file annotates different sets of cells than those in CAS,
           users have the option to delete the labelset in AnnData and add back the labelset from CAS.

Additional Notes:
-----------------
- JSON data is stored in AnnData.uns (excluding barcodes).
"""

import json
import sys
from typing import Optional

import anndata

from cas.file_utils import read_anndata_file, read_json_file

LABELSET_NAME = "name"

LABELSET = "labelset"

LABELSETS = "labelsets"

ANNOTATIONS = "annotations"

CELL_IDS = "cell_ids"

CELL_LABEL = "cell_label"


def merge(cas_path: str, anndata_path: str, validate: bool, output_file_name: str):
    """
    Tests if CAS json and AnnData are compatible and merges CAS into AnnData if possible.

    Args:
        cas_path: The path to the CAS json file.
        anndata_path: The path to the AnnData file.
        validate: Boolean to determine if validation checks will be performed before writing to the output AnnData file.
        output_file_name: Output AnnData file name.

    """
    input_json = read_json_file(cas_path)
    input_anndata = read_anndata_file(anndata_path)

    merge_cas_object(input_json, input_anndata, validate, output_file_name)


def merge_cas_object(input_json: dict, input_anndata: Optional[anndata.AnnData], validate: bool, output_file_name: str):
    """
    Tests if CAS json and AnnData are compatible and merges CAS into AnnData if possible.

    Args:
        input_json: The CAS json object.
        input_anndata: The AnnData object.
        validate: Boolean to determine if validation checks will be performed before writing to the output AnnData file.
        output_file_name: Output AnnData file name.

    """
    test_compatibility(input_anndata, input_json, validate)

    save_cas_to_uns(input_anndata, input_json)
    write_anndata(input_anndata, output_file_name)


def test_compatibility(input_anndata, input_json, validate):
    """
    Tests if CAS and AnnData can be merged.

     Args:
        input_anndata: The AnnData object.
        input_json: The CAS data json object.
        validate: Boolean to determine if validation checks will be performed before writing to the output AnnData file.
    """
    annotations = get_cas_annotations(input_json)
    validate_cell_ids(input_anndata, annotations, validate)

    matching_obs_keys = get_matching_obs_keys(input_anndata, input_json)
    check_labelsets(input_json, input_anndata, matching_obs_keys, validate)


def check_labelsets(cas_json, input_anndata, matching_obs_keys, validate):
    annotations = get_cas_annotations(cas_json)
    derived_cell_ids = get_derived_cell_ids(cas_json)

    for ann in annotations:
        if ann[LABELSET] in matching_obs_keys:
            anndata_labelset_cell_ids = (
                input_anndata.obs.groupby(ann[LABELSET], observed=False)
                .apply(lambda group: set(group.index), include_groups=False)
                .to_dict()
            )
            for cell_label, cell_list in anndata_labelset_cell_ids.items():
                if cell_list == derived_cell_ids.get(
                    str(ann["cell_set_accession"]), set()
                ):
                    handle_matching_labelset(ann, cell_label, input_anndata, validate)
                elif cell_label == ann[CELL_LABEL]:
                    handle_non_matching_labelset(ann, input_anndata, validate, derived_cell_ids)


def get_cas_annotations(input_json):
    return input_json[ANNOTATIONS]


def get_matching_obs_keys(input_anndata, input_json):
    cas_labelset_names = [item[LABELSET_NAME] for item in input_json[LABELSETS]]
    obs_keys = input_anndata.obs_keys()
    matching_obs_keys = list(set(obs_keys).intersection(cas_labelset_names))
    return matching_obs_keys


def handle_matching_labelset(ann, cell_label, input_anndata, validate):
    # Used for label changes
    if cell_label != ann[CELL_LABEL]:
        print(
            f"{ann[CELL_LABEL]} cell ids from CAS match with the cell ids in {cell_label} from anndata. "
            "But they have different cell label."
        )
        if validate:
            sys.exit()
        # add new category to labelset column
        input_anndata.obs[ann[LABELSET]] = input_anndata.obs[
            ann[LABELSET]
        ].cat.add_categories(ann[CELL_LABEL])
        # Overwrite the labelset value with CAS labelset
        input_anndata.obs.loc[ann[CELL_IDS], ann[LABELSET]] = input_anndata.obs.loc[
            ann[CELL_IDS], ann[LABELSET]
        ].map({cell_label: ann[CELL_LABEL]})


def handle_non_matching_labelset(ann, input_anndata, validate, derived_cell_ids):
    # Used for hierarchy changes
    print(
        f"{ann[CELL_LABEL]} cell ids from CAS do not match with the cell ids from anndata. "
        "Please update your CAS json."
    )
    if validate:
        sys.exit()
    # Flush the labelset from anndata
    # input_anndata.obs.loc[list(cell_list), cell_label] = ""
    # Add labelset from CAS to anndata
    cell_ids = derived_cell_ids.get(str(ann["cell_set_accession"]), set())
    input_anndata.obs.loc[list(cell_ids), ann[LABELSET]] = str(ann[CELL_LABEL])


def save_cas_to_uns(input_anndata, input_json):
    # drop cell_ids
    json_without_cell_ids = {
        "author_name": input_json["author_name"],
        "labelset": input_json[LABELSETS],
        "annotations": [
            {key: value for key, value in annotation.items() if key != "cell_ids"}
            for annotation in input_json["annotations"]
        ],
    }
    input_anndata.uns.update({"cas": json.dumps(json_without_cell_ids)})


def validate_cell_ids(input_anndata, annotations, validate):
    # check cell ids
    cas_cell_ids = set()
    for ann in annotations:
        cas_cell_ids.update(ann.get(CELL_IDS, []))
    anndata_cell_ids = set(input_anndata.obs.index)
    # cas -> anndata
    if not cas_cell_ids.issubset(anndata_cell_ids):
        print("Not all members of cell ids from cas exist in anndata.")
        if validate:
            sys.exit()
    # anndata -> cas
    if not anndata_cell_ids.issubset(cas_cell_ids):
        print("Not all members of cell ids from anndata exist in cas.")
        if validate:
            sys.exit()
    return annotations


def write_anndata(input_anndata, output_file_path):
    # Close the AnnData file to prevent blocking
    input_anndata.file.close()
    input_anndata.write(output_file_path)


def get_derived_cell_ids(cas):
    """
    Using the cluster hierarchy derives cell ids for all nodes.
    Args:
        cas: cas json object

    Returns:
        dictionary of cell_set_accession - set of derived cell ids
    """
    derived_cell_ids = dict()

    labelsets = sorted(cas[LABELSETS], key=lambda x: int(x["rank"]))
    for labelset in labelsets:
        ls_annotations = [
            ann for ann in cas[ANNOTATIONS] if ann["labelset"] == labelset["name"]
        ]

        for ann in ls_annotations:
            if "parent_cell_set_accession" in ann:
                cell_ids = set()
                if CELL_IDS in ann and ann[CELL_IDS]:
                    cell_ids = set(ann[CELL_IDS])
                    derived_cell_ids[ann["cell_set_accession"]] = cell_ids
                elif (
                    "cell_set_accession" in ann
                    and ann["cell_set_accession"] in derived_cell_ids
                ):
                    cell_ids = derived_cell_ids[str(ann["cell_set_accession"])]

                if ann["parent_cell_set_accession"] in derived_cell_ids:
                    derived_cell_ids[ann["parent_cell_set_accession"]].update(cell_ids)
                else:
                    derived_cell_ids[ann["parent_cell_set_accession"]] = set(cell_ids)

    return derived_cell_ids
