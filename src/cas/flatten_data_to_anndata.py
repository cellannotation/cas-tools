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
import shutil

import h5py
import numpy as np
import pandas as pd

from cas.file_utils import read_json_file, update_obs_dataset, write_json_to_hdf5

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


def flatten(json_file_path, anndata_file_path, output_file_path):
    """
     Processes and integrates information from a JSON file and an AnnData (Annotated Data) file, creating a new AnnData
     object that incorporates the metadata. The resulting AnnData object is then saved to a new file.

    Args:
        json_file_path: The path to the CAS json file.
        anndata_file_path: The path to the AnnData file.
        output_file_path: Output AnnData file name.
    """
    input_json = read_json_file(json_file_path)

    flatten_cas_object(input_json, anndata_file_path, output_file_path)


# @profile
def flatten_cas_object(input_json, anndata_file_path, output_file_path):
    """
     Processes and integrates information from a JSON file and an AnnData (Annotated Data) file, creating a new AnnData
     object that incorporates the metadata. The resulting AnnData object is then saved to a new file.

    Args:
        input_json: CAS json file.
        anndata_file_path: The path to the AnnData file.
        output_file_path: Output AnnData file name.
    """
    if output_file_path:
        shutil.copy(anndata_file_path, output_file_path)
        anndata_file_path = output_file_path

    annotations = input_json[ANNOTATIONS]
    parent_cell_ids = collect_parent_cell_ids(input_json)

    with h5py.File(anndata_file_path, "r+") as f:
        obs_dataset = f["obs"]
        obs_index = np.array(obs_dataset["CellID"], dtype=str)

        # obs
        flatten_data = process_annotations(annotations, obs_index, parent_cell_ids)
        update_obs_dataset(obs_dataset, flatten_data)

        # uns
        uns_json = generate_uns_json(input_json)
        uns_dataset = f["uns"]
        write_json_to_hdf5(uns_dataset, uns_json)


def process_annotations(annotations, obs_index, parent_cell_ids):
    """
    Processes annotations and generates flattened data for obs dataset.

    Args:
        annotations (list): List of annotations.
        obs_index (np.ndarray): Array representing the index of the obs dataset.
        parent_cell_ids (dict): Dictionary containing parent cell ids.

    Returns:
        dict: Dictionary containing flattened data.
    """
    flatten_data = {}
    for ann in annotations:
        cell_ids = ann.get(
            CELL_IDS, parent_cell_ids.get(ann.get("cell_set_accession", []))
        )
        # Convert cell_ids to a list if it's not already for np.isin
        if not isinstance(cell_ids, list):
            cell_ids = list(cell_ids)
        mask = np.isin(obs_index, cell_ids)

        for k, v in ann.items():
            if k in [CELL_IDS, LABELSET]:
                continue

            key = f"{ann[LABELSET]}--{k}"
            value = ", ".join(
                sorted([str(value) for value in v] if isinstance(v, list) else [str(v)])
            )

            if key not in flatten_data:
                flatten_data[key] = pd.Series("", index=obs_index)
            new_array = flatten_data[key]
            new_array[mask] = value

    return flatten_data


def generate_uns_json(input_json):
    """
    Generates a dictionary representing the uns (unstructured) field in an AnnData object from a given JSON input.

    This function processes information from a JSON input and generates a dictionary that represents the uns (unstructured)
    field in an AnnData object. The resulting dictionary can be used to populate the uns field in the AnnData object.

    Args:
        input_json (dict): A dictionary representing the input JSON data containing annotations.

    Returns:
        dict: A dictionary representing the uns (unstructured) field in an AnnData object, ready to be used as input
              for writing to an AnnData file.

    """
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
    return uns_json


def collect_parent_cell_ids(cas):
    """
    Collects parent cell IDs from the given CAS (Cluster Annotation Service) data.

    This function iterates through labelsets in the CAS data and collects parent cell IDs
    associated with each labelset annotation. It populates and returns a dictionary
    mapping parent cell set accessions to sets of corresponding cell IDs.

    Args:
        cas (dict): The Cluster Annotation Service data containing labelsets and annotations.

    Returns:
        dict: A dictionary mapping parent cell set accessions to sets of corresponding cell IDs.
    """
    parent_cell_ids = dict()

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
                elif (
                    "cell_set_accession" in ann
                    and ann["cell_set_accession"] in parent_cell_ids
                ):
                    cell_ids = parent_cell_ids[ann["cell_set_accession"]]

                if ann["parent_cell_set_accession"] in parent_cell_ids:
                    parent_cell_ids[ann["parent_cell_set_accession"]].update(cell_ids)
                else:
                    parent_cell_ids[ann["parent_cell_set_accession"]] = set(cell_ids)

    return parent_cell_ids
