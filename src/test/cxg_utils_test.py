import unittest
from unittest.mock import patch

from cas.cxg_utils import download_dataset_with_id


class CxGUtilsTests(unittest.TestCase):

    @patch("cellxgene_census.download_source_h5ad")
    def test_download_and_read_dataset_with_id(
            self, mock_download_source_h5ad
    ):
        dataset_id = "dataset_id"

        download_dataset_with_id(dataset_id)

        mock_download_source_h5ad.assert_called_once_with(dataset_id, to_path=f"{dataset_id}.h5ad")

