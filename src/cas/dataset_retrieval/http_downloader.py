import logging
import os
from typing import Optional

import requests
from tqdm import tqdm

from cas.dataset_retrieval.dataset_retriever import (
    DatasetRetriever,
    check_file_exists,
    construct_full_download_path,
    create_directory_if_missing,
)

logging.basicConfig(level=logging.WARNING)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class HTTPDownloader(DatasetRetriever):
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
        url = self.matrix_id
        raw_matrix_id = self.matrix_id.split("/")[-1].split(".")[0]
        default_file_name = f"{raw_matrix_id}.h5ad"
        full_download_path = construct_full_download_path(
            file_name, download_dir, default_file_name
        )
        create_directory_if_missing(full_download_path)
        file_exists = check_file_exists(full_download_path)

        if not file_exists:
            logging.info(f"Downloading dataset with ID '{raw_matrix_id}'...")
            response = requests.get(url, stream=True)
            HTTPDownloader._download_file(response, full_download_path)
            logging.info(f"Download complete. File saved at '{full_download_path}'.")
        return full_download_path

    @staticmethod
    def _download_file(response, full_download_path):
        # Get the total file size from the response headers (if available)
        total_size = int(response.headers.get("content-length", 0))
        # Initialize the tqdm progress bar
        with tqdm(
            total=total_size,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            desc="Downloading",
        ) as progress_bar, open(full_download_path, "wb") as file:
            HTTPDownloader._write_chunks_to_file(response, file, progress_bar)

    @staticmethod
    def _write_chunks_to_file(response, file, progress_bar, chunk_size=8192):
        for chunk in response.iter_content(chunk_size=chunk_size):
            if chunk:  # filter out keep-alive new chunks
                file.write(chunk)
                # Update the progress bar with the size of the chunk
                progress_bar.update(len(chunk))