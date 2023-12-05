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

LABELSET_NAME = "name"

LABELSET = "labelset"

ANNOTATIONS = "annotations"

CELL_IDS = "cell_ids"


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


def flatten(json_file_path, anndata_file_path, output_file_path):
    input_json = read_json_file(json_file_path)
    input_anndata = read_anndata_file(anndata_file_path)
    # obs
    annotations = input_json[ANNOTATIONS]
    for ann in annotations:
        cell_ids = ann.get(CELL_IDS, [])

        for k, v in ann.items():
            if k == CELL_IDS:
                continue
            key = f"{ann[LABELSET]}--{k}"

            non_dict_v = [value for value in v if not isinstance(value, dict)]
            if len(v) > len(non_dict_v):
                print("WARN: dict values are excluded on field '{}'".format(key))
            value = non_dict_v if not isinstance(non_dict_v, list) else ", ".join(sorted(non_dict_v))

            input_anndata.obs[key] = ""

            for index_to_insert in ann[CELL_IDS]:
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
