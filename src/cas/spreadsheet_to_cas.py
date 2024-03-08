from collections import OrderedDict
import json
import os
import logging
from typing import Dict, List, Optional

import anndata as ad
import cellxgene_census
import pandas as pd

from cas.cxg_utils import download_dataset_with_id
from cas.file_utils import read_anndata_file


logging.basicConfig(level=logging.INFO)

LABELSET_COLUMN = "CELL LABELSET NAME"
CELL_LABEL_COLUMN = "CELL TYPE TERM"
CL_TERM_COLUMN = "CL TERM"
EVIDENCE_COLUMN = "EVIDENCE"
MARKER_GENES_COLUMN = "MARKER GENES"
SYNONYM_COLUMN = "SYNONYMS"
CATEGORIES_COLUMN = "CATEGORIES"


def read_spreadsheet(file_path: str, sheet_name: Optional[str]):
    """
    Read the specific sheet from the Excel file into a pandas DataFrame.

    Args:
        file_path (str): Path to the Excel file.
        sheet_name (str, optional): Target sheet name. If not provided, reads the first sheet.

    Returns:
        tuple: Tuple containing metadata (dict), column names (list), and raw data (pd.DataFrame).
    """
    # Read the specific sheet or the first sheet if sheet_name is not provided
    if sheet_name:
        spreadsheet_df = pd.read_excel(file_path, sheet_name=sheet_name, header=None)
    else:
        spreadsheet_df = pd.read_excel(file_path, header=None)
    # Extract meta data (first 8 rows)
    meta_data = dict(spreadsheet_df.iloc[:8, 0:2].to_records(index=False))
    # Extract column names (9th row)
    column_names = spreadsheet_df.iloc[8, :].tolist()
    # Extract raw data (from 10th row onwards)
    raw_data = spreadsheet_df.iloc[9:, :]
    raw_data.columns = column_names

    return meta_data, column_names, raw_data


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


def spreadsheet2cas(
    spreadsheet_file_path: str,
    sheet_name: Optional[str],
    anndata_file_path: str,
    labelset_list: Optional[List[str]],
    output_file_path: str,
):
    """
    Convert a spreadsheet to Cell Annotation Schema (CAS) JSON.

    Args:
        spreadsheet_file_path (str): Path to the spreadsheet file.
        sheet_name (Optional[str]): Target sheet name in the spreadsheet. Can be a string or None.
        anndata_file_path: The path to the AnnData file.
        labelset_list (Optional[List[str]]): List to determine the rank of labelsets in spreadsheet.
        output_file_path (str): Output CAS file name.
    """
    meta_data_result, column_names_result, raw_data_result = read_spreadsheet(
        spreadsheet_file_path, sheet_name
    )

    matrix_file_id = (
        meta_data_result["CxG LINK"].rstrip("/").split("/")[-1].split(".")[0]
    )
    if not anndata_file_path:
        download_dataset_with_id(matrix_file_id)
    dataset_anndata = read_anndata_file(anndata_file_path)

    labelsets = OrderedDict()

    # metadata
    cas = {
        "matrix_file_id": matrix_file_id,
        "cellannotation_schema_version": "TBA",
        "cellannotation_timestamp": "TBA",
        "cellannotation_version": "TBA",
        "cellannotation_url": meta_data_result["CxG LINK"],
        "author_name": "TBA",
        "author_contact": "TBA",
        "orcid": "TBA",
        "annotations": [],
        "labelsets": [],
    }

    # annotations
    stripped_data_result = raw_data_result.map(
        lambda x: x.strip() if isinstance(x, str) else x
    )
    for index, row in stripped_data_result.iterrows():
        labelset = row[LABELSET_COLUMN]
        cell_label = row[CELL_LABEL_COLUMN]
        cell_ontology_term_id = row[CL_TERM_COLUMN]
        rationale = row[EVIDENCE_COLUMN]
        marker_gene_evidence = row[MARKER_GENES_COLUMN]
        synonyms = row[SYNONYM_COLUMN]
        category_fullname = row[CATEGORIES_COLUMN]

        # labelsets.add(labelset)
        labelsets[labelset] = None

        anno = {
            "labelset": labelset,
            "cell_label": cell_label,
            "cell_fullname": cell_label,
            "cell_ontology_term_id": cell_ontology_term_id,
            "cell_ontology_term": cell_label,
            "cell_ids": get_cell_ids(dataset_anndata, labelset, cell_label),
            "rationale": rationale,
            "rationale_dois": meta_data_result["PAPER DOI"],
            "marker_gene_evidence": marker_gene_evidence,
            "synonyms": synonyms,
            "category_fullname": category_fullname,
            "category_cell_ontology_exists": "TBA",
            "category_cell_ontology_term_id": "TBA",
            "category_cell_ontology_term": "TBA",
        }
        cas.get("annotations").append(anno)

    # labelsets
    labelset_rank_dict = calculate_labelset_rank(labelset_list if labelset_list else list(labelsets.keys()))
    for labelset, rank in labelset_rank_dict.items():
        cas.get("labelsets").append({"name": labelset, "description": "TBA", "rank": str(rank)})

    # Write the JSON data to the file
    with open(output_file_path, "w") as json_file:
        json.dump(cas, json_file, indent=2)
