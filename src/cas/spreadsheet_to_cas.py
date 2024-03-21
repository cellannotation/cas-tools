import json
import logging
import re
from collections import OrderedDict
from importlib import resources
from typing import Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd
from cas_schema import schemas

from cas.cxg_utils import download_dataset_with_id
from cas.file_utils import read_anndata_file
from cas.validate import SCHEMA_FILE_MAPPING

logging.basicConfig(level=logging.INFO)

LABELSET_COLUMN = "CELL LABELSET NAME"
CELL_LABEL_COLUMN = "CELL TYPE TERM"
CL_TERM_COLUMN = "CL TERM"
EVIDENCE_COLUMN = "EVIDENCE"
MARKER_GENES_COLUMN = "MARKER GENES"
SYNONYM_COLUMN = "SYNONYMS"
CATEGORIES_COLUMN = "CATEGORIES"


def resolve_ref(schema, ref):
    parts = ref[1:].split('/')
    definition = schema
    for part in parts[1:]:
        definition = definition[part]
    return definition


def read_spreadsheet(file_path: str, sheet_name: Optional[str], schema: dict):
    """
    Read the specific sheet from the Excel file into a pandas DataFrame.

    Args:
        file_path (str): Path to the Excel file.
        sheet_name (str, optional): Target sheet name. If not provided, reads the first sheet.
        schema: Cell annotation schema

    Returns:
        tuple: Tuple containing metadata (dict), column names (list), and raw data (pd.DataFrame).
    """
    if sheet_name:
        spreadsheet_df = pd.read_excel(file_path, sheet_name=sheet_name, header=None)
    else:
        spreadsheet_df = pd.read_excel(file_path, header=None)

    meta_data = {}
    column_names = []
    raw_data = pd.DataFrame()

    header_row_index = None
    metadata_properties = {k: v for k, v in schema['properties'].items() if k not in ["labelsets", "annotations"]}
    for index, row in spreadsheet_df.iterrows():
        first_cell = str(row[0])
        if first_cell.startswith("#"):
            key = first_cell[1:].strip()
            value = row[1] if pd.notnull(row[1]) else ""
            if key in metadata_properties:
                meta_data[key] = value
                if metadata_properties[key]["type"] == "array":
                    meta_data[key] = re.split("[,|]", value)
        else:
            header_row_index = index
            break

    if header_row_index is not None:
        column_names = spreadsheet_df.iloc[header_row_index, :].tolist()
        raw_data = spreadsheet_df.iloc[header_row_index + 1:, :]
        raw_data.columns = column_names
        # Iterate through each column and filter out rows with 'cell_type' in any column
        for column in raw_data.columns:
            raw_data = raw_data[raw_data[column] != 'cell_type']
        # resolve_ref(schema, schema["properties"]["annotations"]["items"]["$ref"])["properties"]["labelset"]
        raw_data = raw_data.where(pd.notnull(raw_data), "")
    else:
        raise ValueError("Header row not found in the spreadsheet.")

    return meta_data, column_names, raw_data


def retrieve_schema(schema_name):
    schema_name = str(schema_name).strip().lower()
    if schema_name not in SCHEMA_FILE_MAPPING:
        raise Exception("Schema name should be one of 'base', 'bican' or 'cap'")
    schema_file = resources.files(schemas) / SCHEMA_FILE_MAPPING[schema_name]
    with schema_file.open("rt") as f:
        schema = json.loads(f.read())
    return schema


def get_cell_ids(dataset: ad.AnnData, labelset: str, cell_label: str) -> List[str]:
    """
    Get cell IDs from an AnnData dataset based on a specified labelset and cell label.

    Args:
        dataset (ad.AnnData): AnnData dataset.
        labelset (str): Labelset to filter.
        cell_label (str): Cell label to match.

    Returns:
        List[str]: List of cell IDs.
    """
    cell_label_lower = str(cell_label).lower()
    return dataset.obs.index[
        dataset.obs[labelset].astype(str).str.lower() == cell_label_lower
        ].tolist()


def calculate_labelset_rank(input_list: List[str]) -> Dict[str, int]:
    """
    Assign ranks to items in a list.

    Args:
        input_list (List[str]): The list of items.

    Returns:
        Dict[str, int]: A dictionary where keys are items from the input list and
        values are their corresponding ranks (0-based).

    """
    return {item: rank for rank, item in enumerate(input_list)}


def custom_lowercase_transform(s):
    """
    Transforms the given string to lowercase except for words that are acronyms or specific cell type names
    which are three characters or less.

    Args:
        s (str): The input string.

    Returns:
        str: The transformed string.
    """

    # Define a function to decide whether a word should stay uppercase
    def transform_word(match):
        word = match.group()
        # If a word is three characters or less, it stays uppercase
        if len(word) <= 3:
            return word
        # Otherwise, convert the word to lowercase
        else:
            return word.lower()

    # Use regex to find words and apply the transformation logic
    transformed_string = re.sub(r"\b[A-Za-z]+\b", transform_word, s)

    return transformed_string


def spreadsheet2cas(
        spreadsheet_file_path: str,
        sheet_name: Optional[str],
        anndata_file_path: Optional[str],
        labelset_list: Optional[List[str]],
        schema_name: Optional[str],
        output_file_path: str,
):
    """
    Convert a spreadsheet to Cell Annotation Schema (CAS) JSON.

    Args:
        spreadsheet_file_path (str): Path to the spreadsheet file.
        sheet_name (Optional[str]): Target sheet name in the spreadsheet. Can be a string or None.
        anndata_file_path: The path to the AnnData file.
        labelset_list (Optional[List[str]]): List of names of observation (obs) fields used to record author cell
        type names, which determine the rank of labelsets in a spreadsheet.
        schema_name (Optional[str]): Name of the CAS schema, can be one of 'base', 'bican' or 'cap'.
        output_file_path (str): Output CAS file name.
    """
    cell_annotation_schema = retrieve_schema(schema_name if schema_name else "cap")
    meta_data_result, column_names_result, raw_data_result = read_spreadsheet(
        spreadsheet_file_path, sheet_name, cell_annotation_schema
    )

    dataset_anndata, matrix_file_id = load_or_fetch_anndata(
        anndata_file_path, meta_data_result
    )

    cas = initialize_cas_structure(matrix_file_id, meta_data_result)

    labelsets = process_and_add_annotations(cas, dataset_anndata, raw_data_result,
                                            column_names_result, cell_annotation_schema)

    process_and_add_labelsets(cas, labelset_list, labelsets)

    # Write the JSON data to the file
    with open(output_file_path, "w") as json_file:
        json.dump(cas, json_file, indent=2)


def process_and_add_labelsets(
        cas: dict, labelset_list: Optional[List[str]], labelsets: OrderedDict
):
    """
    Adds labelsets with ranks to the CAS structure, using either a provided list or keys from `labelsets`.

    Args:
        cas (dict): The CAS structure to update.
        labelset_list (Optional[List[str]]): Custom list to define labelset ranks. If None, uses `labelsets` keys.
        labelsets (OrderedDict): Processed labelsets, used if no custom list is provided.

    Uses `calculate_labelset_rank` to determine ranks, appending each labelset to CAS with its rank.
    """
    labelset_rank_dict = calculate_labelset_rank(
        labelset_list if labelset_list else list(labelsets.keys())
    )
    for labelset, rank in labelset_rank_dict.items():
        cas.get("labelsets").append(
            {"name": labelset, "description": None, "rank": str(rank)}
        )


def process_and_add_annotations(cas, dataset_anndata, raw_data_result, columns, schema):
    """
    Adds processed annotations from raw data to the CAS structure and tracks labelsets.
    Assumes certain external definitions for column names and transformation functions.

    Args:
        cas (dict): The CAS structure to update with annotations.
        dataset_anndata (AnnData): Dataset for annotation context.
        raw_data_result (DataFrame): Raw annotation data.
        columns (list): Column names of raw data to process.
        schema (dict): Cell annotation schema

    Returns:
        OrderedDict: Tracks labelsets encountered, initialized to None.

    Note:
        Requires `custom_lowercase_transform`, `get_cell_ids`, and column constants to be defined.
    """
    stripped_data_result = raw_data_result.map(
        lambda x: x.strip() if isinstance(x, str) else x
    )
    labelsets = OrderedDict()
    for index, row in stripped_data_result.iterrows():
        annotation_properties = schema["definitions"]["Annotation"]["properties"]

        labelsets[row["labelset"]] = None
        anno = {}
        user_annotations = {}
        for column_name in columns:
            if column_name in annotation_properties:
                anno[column_name] = row[column_name]
                if annotation_properties[column_name]["type"] == "array":
                    anno[column_name] = re.split("[,|]", row[column_name])
            elif row[column_name]:
                user_annotations[column_name] = row[column_name]

        anno.update({"cell_ids": get_cell_ids(dataset_anndata, anno["labelset"], anno["cell_label"])})
        if user_annotations:
            anno["user_annotations"] = user_annotations
        cas.get("annotations").append(anno)
    return labelsets


def initialize_cas_structure(matrix_file_id: str, meta_data_result: dict):
    """
    Initializes the Cell Annotation Schema (CAS) structure with basic information and placeholders
    for annotations and labelsets. Fields initialized with None values are omitted in the final output.

    Args:
        matrix_file_id (str): The ID of the matrix file, used within the CAS for identification.
        meta_data_result (dict): Metadata containing at least the 'matrix_file_id' for the CAS URL.

    Returns:
        dict: The initial CAS structure with the matrix file ID, annotation URL, and placeholders
              for future data. Excludes fields that remain None.
    """
    cas_init = {k: v for k, v in meta_data_result.items()}
    cas_init.update({
        "matrix_file_id": f"cxg_dataset:{matrix_file_id}",
        "cellannotation_url": meta_data_result["matrix_file_id"],
        "annotations": [],
        "labelsets": [],
    })
    cas = {k: v for k, v in cas_init.items() if v is not None}
    return cas


def load_or_fetch_anndata(anndata_file_path: str, meta_data_result: dict):
    """
    Loads or fetches an AnnData file, based on a local path or a matrix file ID from metadata.

    Args:
        anndata_file_path (str): Path to an AnnData file, or None to fetch using metadata.
        meta_data_result (dict): Metadata with 'matrix_file_id' for fetching the dataset.

    Returns:
        tuple: (AnnData object, matrix file ID), ready for use.

    Raises:
        ValueError: If 'matrix_file_id' is missing from metadata.
    """
    matrix_file_id = (
        meta_data_result["matrix_file_id"].rstrip("/").split("/")[-1].split(".")[0]
    )
    if not anndata_file_path:
        anndata_file_path = download_dataset_with_id(matrix_file_id)
    dataset_anndata = read_anndata_file(anndata_file_path)
    return dataset_anndata, matrix_file_id
