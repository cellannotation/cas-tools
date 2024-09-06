import logging
import os
from typing import Optional

import cellxgene_census

from cas.dataset_retrieval.dataset_retriever import (
    DatasetRetriever,
    check_file_exists,
    construct_full_download_path,
    create_directory_if_missing,
)


class CxGDownloader(DatasetRetriever):
    def download_data(
        self, file_name: Optional[str] = None, download_dir: Optional[str] = None
    ) -> str:
        """
        Download an AnnData dataset with the specified ID.

        Args:
            file_name: The name of the file to save the downloaded AnnData. If not provided, the dataset
                will be saved with the dataset_id as the file name. The file_name parameter only represents
                the name of the file and does not support absolute or relative paths.
                Use download_dir to specify the directory.
            download_dir: The directory where the AnnData file will be downloaded. If not provided, the
                current working directory will be used. The full path is constructed by combining this directory with `file_name`.

        Returns:
            str: The full path to the downloaded AnnData dataset.
        """
        default_file_name = f"{self.matrix_id}.h5ad"
        full_download_path = construct_full_download_path(
            file_name, download_dir, default_file_name
        )
        create_directory_if_missing(full_download_path)
        file_exists = check_file_exists(full_download_path)

        if not file_exists:
            logging.info(f"Downloading dataset with ID '{self.matrix_id}'...")
            cellxgene_census.download_source_h5ad(
                self.matrix_id, to_path=full_download_path
            )
            logging.info(f"Download complete. File saved at '{full_download_path}'.")
        return full_download_path
