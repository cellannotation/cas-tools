#!/usr/bin/env python3

"""
Script Description:
This script processes and integrates information from a JSON file and an AnnData (Annotated Data) file,
creating a new AnnData object that incorporates the metadata. The resulting AnnData object is then saved to a new file.

Key Features:
1. Parses command-line arguments for input JSON file, input AnnData file, and output file.
2. Reads and processes the input JSON file and AnnData file.
3. Updates the AnnData object with information from the JSON annotations and root keys.
4. Writes the modified AnnData object to a specified output file.
"""

from cas.file_utils import read_json_file, read_anndata_file
from cas.anndata_conversion import test_compatibility

LABELSET_NAME = "name"

LABELSET = "labelset"

LABELSETS = "labelsets"

ANNOTATIONS = "annotations"

CELL_IDS = "cell_ids"

CELL_LABEL = "cell_label"


def is_list_of_strings(var):
    """
    Check if a value is a list of strings.

    Parameters:
        var (list or any): The value to be checked.

    Returns:
        bool: True if the value is a list containing only string elements,
              False otherwise.

    """
    return isinstance(var, list) and all(isinstance(item, str) for item in var)


def flatten(json_file_path, anndata_file_path, validate, output_file_path):
    """
     Processes and integrates information from a JSON file and an AnnData (Annotated Data) file, creating a new AnnData
     object that incorporates the metadata. The resulting AnnData object is then saved to a new file.

    Args:
        json_file_path: The path to the CAS json file.
        anndata_file_path: The path to the AnnData file.
        validate: Boolean to determine if validation checks will be performed before writing to the output AnnData file.
        output_file_path: Output AnnData file name.
    """
    input_json = read_json_file(json_file_path)
    input_anndata = read_anndata_file(anndata_file_path)

    if validate:
        test_compatibility(input_anndata, input_json, validate)
    # obs
    annotations = input_json[ANNOTATIONS]

    parent_cell_ids = collect_parent_cell_ids(input_json)

    for ann in annotations:
        cell_ids = []
        if CELL_IDS in ann and ann[CELL_IDS]:
            cell_ids = ann[CELL_IDS]
        elif "cell_set_accession" in ann and ann["cell_set_accession"] in parent_cell_ids:
            cell_ids = list(parent_cell_ids[ann["cell_set_accession"]])

        for k, v in ann.items():
            if k == CELL_IDS or k == LABELSET:
                continue
            if k == CELL_LABEL:
                key = ann[LABELSET]
            else:
                key = f"{ann[LABELSET]}--{k}"

            value = v
            if isinstance(v, list):
                non_dict_v = [value for value in v if not isinstance(value, dict)]
                value = ", ".join(sorted(non_dict_v))
                if len(v) > len(non_dict_v):
                    print("WARN: dict values are excluded on field '{}'".format(key))

            for index_to_insert in cell_ids:
                input_anndata.obs.at[index_to_insert, key] = value
    # uns
    uns_json = {}
    root_keys = list(input_json.keys())
    root_keys.remove(ANNOTATIONS)
    for key in root_keys:
        value = input_json[key]
        if is_list_of_strings(value):
            uns_json[key] = ", ".join(sorted(value))
        elif isinstance(value, str):
            uns_json[key] = value
        else:
            for labelset in value:
                for k, v in labelset.items():
                    if k == LABELSET_NAME:
                        continue
                    new_key = f"{labelset.get(LABELSET_NAME, '')}--{k}"
                    uns_json.update({new_key: v})
    input_anndata.uns.update(uns_json)
    # Close the AnnData file to prevent blocking
    input_anndata.file.close()
    input_anndata.write(output_file_path)


def collect_parent_cell_ids(cas):
    parent_cell_ids = dict()

    labelsets = sorted(cas[LABELSETS], key=lambda x: int(x["rank"]))
    for labelset in labelsets:
        ls_annotations = [ann for ann in cas[ANNOTATIONS] if ann["labelset"] == labelset["name"]]

        for ann in ls_annotations:
            if "parent_cell_set_accession" in ann:
                cell_ids = set()
                if CELL_IDS in ann and ann[CELL_IDS]:
                    cell_ids = set(ann[CELL_IDS])
                elif "cell_set_accession" in ann and ann["cell_set_accession"] in parent_cell_ids:
                    cell_ids = parent_cell_ids[ann["cell_set_accession"]]

                if ann["parent_cell_set_accession"] in parent_cell_ids:
                    parent_cell_ids[ann["parent_cell_set_accession"]].update(cell_ids)
                else:
                    parent_cell_ids[ann["parent_cell_set_accession"]] = set(cell_ids)

    return parent_cell_ids
