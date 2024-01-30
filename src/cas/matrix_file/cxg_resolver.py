import os
from os.path import isfile, join
from typing import Optional

import anndata

from cas.file_utils import read_anndata_file
from cas.matrix_file.base_resolver import BaseMatrixFileResolver

H5AD_FILE_EXTENSION = ".h5ad"


class CxGDatasetResolver(BaseMatrixFileResolver):
    """
    CellxGene dataset resolver.
    """

    def __init__(self, cache_folder_path):
        """
        Initializer.
        Params:
            cache_folder_path: matrix file cache folder path
        """
        self.cache_folder_path = cache_folder_path
        self.cached_datasets = list_cached_datasets(cache_folder_path)

    def resolve_matrix_file(self, dataset_id) -> Optional[anndata.AnnData]:
        return read_anndata_file(self.resolve_matrix_file_path(dataset_id))

    def resolve_matrix_file_path(self, dataset_id) -> str:
        dataset_file_name = dataset_id + H5AD_FILE_EXTENSION
        if dataset_file_name in self.cached_datasets:
            return os.path.join(self.cache_folder_path, dataset_file_name)

        # TODO implement dataset download operation
        raise Exception(
            "Dataset could not be found at location: '{}'".format(
                os.path.join(self.cache_folder_path, dataset_file_name)
            )
        )


def list_cached_datasets(cache_folder_path):
    if os.path.isdir(cache_folder_path):
        return [
            f
            for f in os.listdir(cache_folder_path)
            if os.path.isfile(os.path.join(cache_folder_path, f))
        ]
    else:
        raise Exception(
            "CellxGene dataset cache folder is not accessible: '{}'".format(
                cache_folder_path
            )
        )
