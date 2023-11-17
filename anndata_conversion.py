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
from typing import Optional

import anndata

from anndata_conversion_utils import (
    check_labelsets,
    get_cas_annotations,
    get_matching_obs_keys,
    save_cas_to_uns,
    validate_cell_ids,
    write_anndata,
)


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--json", required=True, help="Path to the CAS JSON schema file."
    )
    parser.add_argument("--anndata", required=True, help="Path to the AnnData file.")
    # TODO find a better argument name and Help message.
    parser.add_argument(
        "-v",
        "--validate",
        action="store_true",
        help="Perform validation checks before writing to the output AnnData file.",
    )
    parser.add_argument(
        "--output",
        default="output.h5ad",
        help="Output AnnData file name (default: output.h5ad).",
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

    annotations = get_cas_annotations(input_json)
    validate_cell_ids(input_anndata, annotations, validate)

    matching_obs_keys = get_matching_obs_keys(input_anndata, input_json)
    check_labelsets(annotations, input_anndata, matching_obs_keys, validate)

    save_cas_to_uns(input_anndata, input_json)
    write_anndata(input_anndata, output_file_path)


if __name__ == "__main__":
    main()
