import logging
import os
from typing import Optional

import cellxgene_census

from cas.dataset_retrieval.dataset_retriever import DatasetRetriever


class CxGDownloader(DatasetRetriever):
    def download_data(self, file_path: Optional[str] = None) -> str:
        """
        Download an AnnData dataset with the specified ID.

        Args:
            file_path (Optional[str], optional): The file path to save the downloaded AnnData. If not provided,
                the dataset will be saved in the current working directory with the dataset_id as the file name.
                Supports both absolute and relative paths.

        Returns:
            str: The path to the downloaded AnnData dataset
        """
        default_file_name = f"{self.matrix_id}.h5ad"
        anndata_file_path = default_file_name if file_path is None else file_path

        anndata_file_path = os.path.abspath(anndata_file_path)

        # Check if the file already exists
        if os.path.exists(anndata_file_path):
            logging.info(
                f"File '{anndata_file_path}' already exists. Skipping download."
            )
            return anndata_file_path

        # Ensure the directory exists
        directory = os.path.dirname(anndata_file_path)
        if directory and not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)

        logging.info(f"Downloading dataset with ID '{self.matrix_id}'...")
        cellxgene_census.download_source_h5ad(self.matrix_id, to_path=anndata_file_path)
        logging.info(f"Download complete. File saved at '{anndata_file_path}'.")
        return anndata_file_path
