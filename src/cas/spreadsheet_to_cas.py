import json
import os
import logging
from typing import List, Optional

import anndata as ad
import cellxgene_census
import pandas as pd

from cas.file_utils import read_anndata_file


logging.basicConfig(level=logging.INFO)


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
    return dataset.obs.index[
        dataset.obs[labelset].str.lower() == cell_label.lower()
    ].tolist()


def download_and_read_dataset_with_id(dataset_id: str) -> ad.AnnData:
    """
    Download and read an AnnData dataset with a specified ID.

    Args:
        dataset_id (str): ID of the dataset.

    Returns:
        ad.AnnData: AnnData object.
    """
    anndata_file_path = f"{dataset_id}.h5ad"
    # Check if the file already exists
    if os.path.exists(anndata_file_path):
        print(f"File '{anndata_file_path}' already exists. Skipping download.")
    else:
        logging.info(f"Downloading dataset with ID '{dataset_id}'...")
        cellxgene_census.download_source_h5ad(dataset_id, to_path=anndata_file_path)
        logging.info(f"Download complete. File saved at '{anndata_file_path}'.")
    anndata = read_anndata_file(anndata_file_path)
    return anndata


def spreadsheet2cas(spreadsheet_file_path: str, sheet_name: Optional[str], output_file_path: str):
    """
    Convert a spreadsheet to Cell Annotation Schema (CAS) JSON.

    Args:
        spreadsheet_file_path (str): Path to the spreadsheet file.
        sheet_name (Optional[str]): Target sheet name in the spreadsheet. Can be a string or None.
        output_file_path (str): Output CAS file name.
    """
    meta_data_result, column_names_result, raw_data_result = read_spreadsheet(
        spreadsheet_file_path, sheet_name
    )

    matrix_file_id = (
        meta_data_result["CxG LINK"].rstrip("/").split("/")[-1].split(".")[0]
    )
    dataset_anndata = download_and_read_dataset_with_id(matrix_file_id)
    labelsets = set()

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
        "labelsets": []
    }

    # annotations
    stripped_data_result = raw_data_result.map(
        lambda x: x.strip() if isinstance(x, str) else x
    )
    for index, row in stripped_data_result.iterrows():
        labelsets.add(row[0])
        anno = {
            "labelset": row[0],
            "cell_label": row[1],
            "cell_fullname": row[1],
            "cell_ontology_term_id": row[2],
            "cell_ontology_term": row[1],
            "cell_ids": get_cell_ids(dataset_anndata, row[0], row[1]),
            "rationale": row[8],
            "rationale_dois": meta_data_result["PAPER DOI"],
            "marker_gene_evidence": row[3],
            "synonyms": row[6],
            "category_fullname": row[7],
            "category_cell_ontology_exists": "TBA",
            "category_cell_ontology_term_id": "TBA",
            "category_cell_ontology_term": "TBA",
        }
        cas.get("annotations").append(anno)

    # labelsets
    for labelset in labelsets:
        labelsets_dict = {
            "name": labelset,
            "description": "TBA",
            "rank": "TBA"
        }
        cas.get("labelsets").append(labelsets_dict)

    # Write the JSON data to the file
    with open(output_file_path, "w") as json_file:
        json.dump(cas, json_file, indent=2)
