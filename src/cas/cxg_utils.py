"""
cxg_utils.py

This module provides utility functions for working with AnnData datasets in the context of the CellxGene Census library.
"""
import os
import logging

import cellxgene_census


def download_dataset_with_id(dataset_id: str):
    """
    Download an AnnData dataset with a specified ID.

    Args:
        dataset_id (str): ID of the dataset.

    """
    anndata_file_path = f"{dataset_id}.h5ad"
    # Check if the file already exists
    # Currently using dataset id as file names for the downloads
    if os.path.exists(anndata_file_path):
        print(f"File '{anndata_file_path}' already exists. Skipping download.")
    else:
        logging.info(f"Downloading dataset with ID '{dataset_id}'...")
        cellxgene_census.download_source_h5ad(dataset_id, to_path=anndata_file_path)
        logging.info(f"Download complete. File saved at '{anndata_file_path}'.")
