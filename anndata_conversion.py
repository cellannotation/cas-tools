#!/usr/bin/env python3

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

Command-line Arguments:
-----------------------
--json      : Path to the CAS JSON schema file.
--anndata   : Path to the AnnData file. Ideally, the location will be specified by a resolvable path in the CAS file.
--validate  : Perform validation checks before writing to the output AnnData file.
--output    : Output AnnData file name (default: output.h5ad).

Usage Example:
--------------
python script.py --json path/to/CAS_schema.json --anndata path/to/input_anndata.h5ad --validate --output path/to/output.h5ad
"""

import argparse
import json
import sys
from typing import Optional

import anndata

LABELSET_NAME = "name"

LABELSET = "labelset"

ANNOTATIONS = "annotations"

CELL_IDS = "cell_ids"

CELL_LABEL = "cell_label"


def read_json_file(file_path):
    """
    Reads and parses a JSON file into a Python dictionary.

    Args:
        file_path (str): The path to the JSON file.

    Returns:
        dict: The JSON data as a Python dictionary.

    Returns None if the file does not exist or if there is an issue
    parsing the JSON content.

    Example:
        json_data = read_json_file('path/to/your/file.json')
        if json_data is not None:
            # Use the parsed JSON data as a dictionary
            print(json_data)
    """
    try:
        with open(file_path, "r") as file:
            data = json.load(file)
            return data
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error reading JSON file: {e}")
        return None


def read_anndata_file(file_path: str) -> Optional[anndata.AnnData]:
    """Load anndata object from a file.

    Args:
        file_path: The path to the file containing the anndata object.

    Returns:
        The loaded anndata object if successful, else None.
    """
    try:
        anndata_obj = anndata.read_h5ad(file_path, backed="r")
        return anndata_obj
    except Exception as e:
        print(f"An error occurred while loading the file: {e}")
        return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", required=True, help="Input JSON file path")
    parser.add_argument("--anndata", required=True, help="Input AnnData file path")
    parser.add_argument(
        "-v",
        "--validate",
        help="Perform validation checks before writing to the output AnnData file. "
        "Checks include:\n"
        "1. Ensure all barcodes (cell IDs) in the CAS file exist in the AnnData file.\n"
        "2. Ensure all barcodes (cell IDs) in the AnnData file exist in the CAS file.\n"
        "3. Check for obs keys in the AnnData file matching labelset names in CAS.\n"
        "   - If matches are found, a warning will be issued with options to overwrite or terminate.",
        default=False,
    )

    parser.add_argument(
        "--output",
        help="Output AnnData file name (default: output.h5ad)",
        default="output.h5ad",
    )

    args = parser.parse_args()
    json_file_path = args.json
    anndata_file_path = args.anndata
    output_file_path = args.output
    validate = args.validate

    if anndata_file_path == output_file_path:
        raise ValueError("--anndata and --output cannot be the same")

    input_json = read_json_file(json_file_path)
    input_anndata = read_anndata_file(anndata_file_path)

    # check cell ids
    annotations = input_json[ANNOTATIONS]
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

    # check obs keys
    cas_labelset_names = [item[LABELSET_NAME] for item in input_json[LABELSET]]
    obs_keys = input_anndata.obs_keys()
    matching_obs_keys = list(set(obs_keys).intersection(cas_labelset_names))

    for ann in annotations:
        if ann[LABELSET] in matching_obs_keys:
            anndata_labelset_cellIDs = (
                input_anndata.obs.groupby(ann[LABELSET], observed=False)
                .apply(lambda group: set(group.index))
                .to_dict()
            )
            for cell_label, cell_list in anndata_labelset_cellIDs.items():
                if cell_list == set(ann[CELL_IDS]):
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
                        input_anndata.obs.loc[
                            ann[CELL_IDS], ann[LABELSET]
                        ] = input_anndata.obs.loc[ann[CELL_IDS], ann[LABELSET]].map(
                            {cell_label: ann[CELL_LABEL]}
                        )
                elif cell_label == ann[CELL_LABEL]:
                    print(
                        f"{ann[CELL_LABEL]} cell ids from CAS do not match with the cell ids from anndata. "
                        "Please update your CAS json."
                    )
                    if validate:
                        sys.exit()
                    # Flush the labelset from anndata
                    input_anndata.obs.loc[cell_list, cell_label] = None
                    # Add labelset from CAS to anndata
                    input_anndata.obs.loc[ann[CELL_IDS], ann[LABELSET]] = ann[
                        CELL_LABEL
                    ]

    # drop cell_ids
    json_without_cell_ids = {
        "author_name": input_json["author_name"],
        "labelset": input_json["labelset"],
        "annotations": [
            {key: value for key, value in annotation.items() if key != "cell_ids"}
            for annotation in input_json["annotations"]
        ],
    }
    input_anndata.uns.update({"cas": json.dumps(json_without_cell_ids)})
    # Close the AnnData file to prevent blocking
    input_anndata.file.close()
    input_anndata.write(output_file_path)
