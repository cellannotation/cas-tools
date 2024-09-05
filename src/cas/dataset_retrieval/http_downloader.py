import logging
import os
from typing import Optional

import requests
from tqdm import tqdm

from cas.dataset_retrieval.dataset_retriever import (
    DatasetRetriever,
    check_file_exists,
    create_directory_if_missing,
)

logging.basicConfig(level=logging.WARNING)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class HTTPDownloader(DatasetRetriever):
    def download_data(self, file_name: Optional[str] = None) -> str:
        url = self.matrix_id
        raw_matrix_id = self.matrix_id.split("/")[-1].split(".")[0]
        default_file_name = f"{raw_matrix_id}.h5ad"
        anndata_file_path = default_file_name if file_name is None else file_name

        anndata_file_path = check_file_exists(anndata_file_path)

        create_directory_if_missing(anndata_file_path)

        response = requests.get(url, stream=True)

        # Get the total file size from the response headers (if available)
        total_size = int(response.headers.get("content-length", 0))

        logging.info(f"Downloading dataset with ID '{raw_matrix_id}'...")
        # Initialize the tqdm progress bar
        with tqdm(
            total=total_size,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            desc="Downloading",
        ) as progress_bar:
            with open(anndata_file_path, "wb") as file:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:  # filter out keep-alive new chunks
                        file.write(chunk)
                        # Update the progress bar with the size of the chunk
                        progress_bar.update(len(chunk))

        logging.info(f"Download complete. File saved at '{anndata_file_path}'.")
        return anndata_file_path
