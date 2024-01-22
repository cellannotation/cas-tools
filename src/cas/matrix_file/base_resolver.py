import abc
from typing import Optional

import anndata


class BaseMatrixFileResolver(metaclass=abc.ABCMeta):
    """
    Base abstract matrix file Resolver
    """

    @abc.abstractmethod
    def resolve_matrix_file(self, dataset_id) -> Optional[anndata.AnnData]:
        """
        Resolves matrix file identified by the given dataset_id.
        Parameters:
            dataset_id: dataset identifier
        Returns:
            AnnData object
        """
        pass

    @abc.abstractmethod
    def resolve_matrix_file_path(self, dataset_id) -> str:
        """
        Resolves matrix file identified by the given dataset_id and returns its path.
        Parameters:
            dataset_id: dataset identifier
        Returns:
            AnnData file path
        """
        pass
