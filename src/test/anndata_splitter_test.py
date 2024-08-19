import json
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from cas.anndata_splitter import split_anndata, split_anndata_to_file


class TestSplitAnndata(unittest.TestCase):
    def setUp(self):
        self.current_directory = Path.cwd()
        self.test_data_path = Path("test_data/anndata_splitter")
        self.test_data_path.mkdir(parents=True, exist_ok=True)

        # Sample AnnData with aligned values
        obs_data = {
            "cell_type": [f"type_{i % 3}" for i in range(100)],
        }
        obs = pd.DataFrame(data=obs_data, index=[f"cell{i}" for i in range(100)])

        var_data = {
            "gene_name": [f"gene_{i}" for i in range(1000)],
        }
        var = pd.DataFrame(data=var_data, index=[f"gene{i}" for i in range(1000)])

        X = sp.csr_matrix(np.random.randn(100, 1000))
        self.adata = AnnData(X=X, obs=obs, var=var)

        self.anndata_file_path = self.test_data_path / "sample_anndata.h5ad"
        self.adata.write(self.anndata_file_path)

        # Create 1/3 of the cell{i} from obs and split among two CAS JSON files
        cell_ids_part1 = [f"cell{i}" for i in range(0, 33)]
        cell_ids_part2 = [f"cell{i}" for i in range(33, 67)]

        # Sample CAS JSON data
        self.cas_json_data = [
            {"annotations": [{"cell_ids": cell_ids_part1}]},
            {"annotations": [{"cell_ids": cell_ids_part2}]},
        ]

        # Write the CAS JSON data to files
        self.cas_json_paths = []
        for i, cas_data in enumerate(self.cas_json_data):
            cas_json_path = self.test_data_path / f"cas_data_{i}.json"
            self.cas_json_paths.append(cas_json_path)
            with open(cas_json_path, "w") as f:
                json.dump(cas_data, f)

    def tearDown(self):
        # Clean up the test data directory
        for file_path in self.test_data_path.iterdir():
            if file_path.is_file():
                file_path.unlink()
        self.test_data_path.rmdir()
        # Clean up split h5ad files
        file1 = self.current_directory / "split_cas_data_0.h5ad"
        file2 = self.current_directory / "split_cas_data_1.h5ad"
        if file1.exists():
            file1.unlink()
        if file2.exists():
            file2.unlink()

    @patch("cas.anndata_splitter.read_json_file")
    @patch("cas.anndata_splitter.read_anndata_file")
    def test_split_anndata_to_file(self, mock_read_anndata_file, mock_read_json_file):
        mock_read_anndata_file.return_value = self.adata
        mock_read_json_file.side_effect = [self.cas_json_data[0], self.cas_json_data[1]]

        split_anndata_to_file(
            str(self.anndata_file_path),
            [str(path) for path in self.cas_json_paths],
            multiple_outputs=True,
        )

        # Check that the read functions were called
        mock_read_anndata_file.assert_called_once_with(str(self.anndata_file_path))
        self.assertEqual(mock_read_json_file.call_count, 2)

        # Check that the output files were created
        self.assertTrue((self.current_directory / "split_cas_data_0.h5ad").exists())
        self.assertTrue((self.current_directory / "split_cas_data_1.h5ad").exists())

    def test_split_anndata_multiple_outputs(self):
        result = split_anndata(
            self.adata,
            {
                "cas_data_0.json": self.cas_json_data[0],
                "cas_data_1.json": self.cas_json_data[1],
            },
            multiple_outputs=True,
        )

        # Check that the result is a list of AnnData objects
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)
        self.assertIsInstance(result[0], AnnData)
        self.assertIsInstance(result[1], AnnData)

        # Verify that the cells match the expected splits
        self.assertListEqual(
            result[0].obs.index.tolist(), [f"cell{i}" for i in range(0, 33)]
        )
        self.assertListEqual(
            result[1].obs.index.tolist(), [f"cell{i}" for i in range(33, 67)]
        )

    def test_split_anndata_single_output(self):
        combined_cas_json_data = {
            "annotations": [{"cell_ids": [f"cell{i}" for i in range(0, 67)]}]
        }

        result = split_anndata(
            self.adata,
            {"cas_data_combined.json": combined_cas_json_data},
            multiple_outputs=False,
        )

        # Check that the result is a single AnnData object in a list
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertIsInstance(result[0], AnnData)

        # Verify that the cells match the expected combined set
        self.assertListEqual(
            result[0].obs.index.tolist(), [f"cell{i}" for i in range(0, 67)]
        )
