import logging
import os
from abc import ABC, abstractmethod
from typing import Optional


class DatasetRetriever(ABC):
    def __init__(self, matrix_id, **kwargs):
        self.matrix_id = matrix_id
        self.kwargs = kwargs

    @abstractmethod
    def download_data(
        self, file_name: Optional[str] = None, download_dir: Optional[str] = None
    ) -> str:
        """
        Download data based on the matrix_id.

        Args:
            file_name: The name of the file where the dataset will be saved. Defaults to None.
            download_dir: The directory of the file where the dataset will be saved. Defaults to None

        Returns:
            str: The path to the downloaded dataset file.
        """
        pass

    @classmethod
    def create(cls, matrix_id: str, **kwargs):
        """
        Factory method to create a dataset downloader based on the matrix_id format.

        Args:
            matrix_id: Identifier or URL for the dataset.
            **kwargs: Additional arguments for downloader initialization.

        Returns:
            DatasetRetriever: An instance of a subclass of DatasetRetriever based on matrix_id type.

        Raises:
            ValueError: If the matrix_id format is unsupported.
        """
        if matrix_id.startswith("cxg_dataset:"):
            from cas.dataset_retrieval.cxg_downloader import CxGDownloader

            return CxGDownloader(matrix_id.split(":")[-1], **kwargs)
        elif matrix_id.startswith(
            "https://datasets.cellxgene.cziscience.com/"
        ) and matrix_id.endswith(".h5ad"):
            from cas.dataset_retrieval.http_downloader import HTTPDownloader

            return HTTPDownloader(matrix_id, **kwargs)
        elif matrix_id.startswith("s3"):
            from cas.dataset_retrieval.s3_downloader import S3Downloader

            return S3Downloader(matrix_id, **kwargs)
        else:
            raise ValueError(f"Unsupported matrix_id: {matrix_id}")


def check_file_exists(file_name: str) -> bool:
    """
    Check if a file exists at the specified path.

    Args:
        file_name: Path to the file.

    Returns:
        True if the file exists, otherwise False.
    """
    anndata_file_path = os.path.abspath(file_name)

    if os.path.exists(anndata_file_path):
        logging.info(f"File '{anndata_file_path}' already exists. Skipping download.")
        return True

    return False


def create_directory_if_missing(file_name: str):
    """
    Ensure the directory for the specified file exists. If it doesn't, create it.

    Args:
        file_name: Path to the file whose directory needs to be checked or created.
    """
    directory = os.path.dirname(file_name)
    if directory and not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)


def construct_full_download_path(
    file_name: Optional[str], download_dir: Optional[str], default_file_name: str
) -> str:
    """
    Constructs the full file path for the AnnData file based on the provided file name and download directory.

    Args:
        file_name: The name of the file. If None, `default_file_name` will be used.
        download_dir: The directory where the file should be saved. If None, the current directory will be used.
        default_file_name: The default file name to use if `file_name` is not provided.

    Returns:
        str: The full file path combining the directory and file name.
    """
    anndata_file_path = default_file_name if file_name is None else file_name
    full_download_path = os.path.join(download_dir or "", anndata_file_path)

    return full_download_path
