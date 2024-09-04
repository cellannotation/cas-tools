from abc import ABC, abstractmethod


class DatasetRetriever(ABC):
    def __init__(self, matrix_id, **kwargs):
        self.matrix_id = matrix_id
        self.kwargs = kwargs

    @abstractmethod
    def download_data(self, file_name=None) -> str:
        """
        Download data based on the matrix_id.

        Args:
            file_name (str, optional): The name of the file where the dataset will be saved. Defaults to None.

        Returns:
            str: The path to the downloaded dataset file.
        """
        pass

    @classmethod
    def create(cls, matrix_id, **kwargs):
        if matrix_id.startswith("cxg_dataset:"):
            from cas.dataset_retrieval.cxg_downloader import CxGDownloader
            return CxGDownloader(matrix_id.split(":")[-1], **kwargs)
        elif matrix_id.startswith("https://datasets.cellxgene.cziscience.com/") and matrix_id.endswith(".h5ad"):
            from cas.dataset_retrieval.http_downloader import HTTPDownloader
            return HTTPDownloader(matrix_id, **kwargs)
        elif matrix_id.startswith("s3"):
            from cas.dataset_retrieval.s3_downloader import S3Downloader
            return S3Downloader(matrix_id, **kwargs)
        else:
            raise ValueError(f"Unsupported matrix_id: {matrix_id}")
