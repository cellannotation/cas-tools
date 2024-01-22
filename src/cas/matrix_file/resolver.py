from typing import Optional

import anndata

from cas.matrix_file.cxg_resolver import CxGDatasetResolver

CXG_PREFIX = "CellXGene_dataset"


def resolve_matrix_file(
    matrix_file_id: str, cache_folder_path: str = None
) -> Optional[anndata.AnnData]:
    """
    Resolves matrix file identified by the given matrix_file_id.

    Parameters:
        matrix_file_id: dataset identifier
        cache_folder_path: (Optional) matrix file cache folder path
    Returns:
        AnnData object
    """
    protocol = matrix_file_id.split(":")[0]
    dataset_id = matrix_file_id.split(":")[1]

    if protocol == CXG_PREFIX:
        resolver = CxGDatasetResolver(cache_folder_path)
    else:
        raise Exception("Unrecognised matrix file protocol: '{}'".format(protocol))

    return resolver.resolve_matrix_file(dataset_id)


def resolve_matrix_file_path(matrix_file_id: str, cache_folder_path: str = None) -> str:
    """
    Resolves matrix file identified by the given matrix_file_id.

    Parameters:
        matrix_file_id: dataset identifier
        cache_folder_path: (Optional) matrix file cache folder path
    Returns:
        AnnData file path
    """
    protocol = matrix_file_id.split(":")[0]
    dataset_id = matrix_file_id.split(":")[1]

    if protocol == CXG_PREFIX:
        resolver = CxGDatasetResolver(cache_folder_path)
    else:
        raise Exception("Unrecognised matrix file protocol: '{}'".format(protocol))

    return resolver.resolve_matrix_file_path(dataset_id)
