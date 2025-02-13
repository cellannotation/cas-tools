import unittest
from unittest.mock import MagicMock, patch

import pandas as pd

from cas.utils.validation_utils import (
    validate_labelset_markers,
    validate_markers,
    compare_labelsets_cas_obs,
    validate_labelset_values,
    infer_cas_cell_hierarchy,
    infer_obs_cell_hierarchy,
)


class TestMarkerValidation(unittest.TestCase):
    def setUp(self):
        # Setup a sample dataframe and CAS dictionary
        self.marker_column = "gene_markers"
        self.data = {"gene_markers": ["gene1", "gene2", "gene3"]}
        self.df = pd.DataFrame(self.data)
        self.df_var = pd.DataFrame(self.data)
        self.adata = MagicMock()
        self.adata.var = self.df_var

        self.cas = {
            "annotations": [
                {
                    "labelset": "set1",
                    "cell_label": "A",
                    "marker_gene_evidence": ["gene1", "gene2"],
                },
                {
                    "labelset": "set2",
                    "cell_label": "B",
                    "marker_gene_evidence": ["gene4"],
                }
                # this should trigger a warning
            ]
        }

    def test_validate_markers(self):
        # Test that all markers are found and no exception is raised
        result = validate_markers(self.cas, self.adata, self.marker_column)
        self.assertTrue(result)

        # Test for KeyError when column does not exist
        with self.assertRaises(KeyError):
            validate_markers(self.cas, self.adata, "nonexistent_column")

    @patch("cas.utils.validation_utils.logger.warning")
    def test_validate_labelset_markers(self, mock_logger_warning):
        # Test with missing marker
        annotation = {
            "labelset": "set2",
            "cell_label": "B",
            "marker_gene_evidence": ["gene4"],
        }
        validate_labelset_markers(
            annotation, self.adata.var[self.marker_column].tolist()
        )
        mock_logger_warning.assert_called_once_with(
            "Not all marker genes from set2-B pair exist in anndata's var section. Missing markers: gene4"
        )

        # Test with all markers present
        mock_logger_warning.reset_mock()
        annotation = {
            "labelset": "set1",
            "cell_label": "A",
            "marker_gene_evidence": ["gene1", "gene2"],
        }
        validate_labelset_markers(
            annotation, self.adata.var[self.marker_column].tolist()
        )
        mock_logger_warning.assert_not_called()

    def test_compare_labelsets_cas_obs(self):
        obs = pd.DataFrame({"set1": ["A", "B"], "set2": ["X", "Y"]})

        cas = {"labelsets": [{"name": "set1"}, {"name": "set2"}]}

        # All labelsets exist in obs (should return True)
        result = compare_labelsets_cas_obs(cas, obs)
        self.assertTrue(result)

        # One labelset is missing in obs (should return False)
        modified_obs = obs.drop(columns=["set1"])

        result = compare_labelsets_cas_obs(cas, modified_obs)
        self.assertFalse(result)

        # Empty obs DataFrame (should return False)
        empty_obs = pd.DataFrame()

        result = compare_labelsets_cas_obs(cas, empty_obs)
        self.assertFalse(result)

        # Empty CAS JSON (should return True)
        empty_cas = {"labelsets": []}

        result = compare_labelsets_cas_obs(empty_cas, obs)
        self.assertTrue(result)

    def test_validate_labelset_values(self):
        obs = pd.DataFrame({"set1": ["A", "B", "C"], "set2": ["X", "Y", "Z"]})

        cas = {
            "annotations": [
                {"labelset": "set1", "cell_label": "A"},
                {"labelset": "set1", "cell_label": "B"},
                {"labelset": "set1", "cell_label": "C"},
                {"labelset": "set2", "cell_label": "X"},
                {"labelset": "set2", "cell_label": "Y"},
                {"labelset": "set2", "cell_label": "Z"},
            ]
        }

        # All labelset members exist (should return True)
        valid_obs = pd.DataFrame({"set1": ["A", "B", "C"], "set2": ["X", "Y", "Z"]})
        self.assertTrue(validate_labelset_values(cas, valid_obs))

        # Missing labelset member in obs (should return False)
        modified_obs = pd.DataFrame(
            {"set1": ["A", "B", "C"], "set2": ["X", "Y", "z"]}  # "Z" is missing
        )
        self.assertFalse(validate_labelset_values(cas, modified_obs))

        # Missing entire labelset column in obs (should return False)
        missing_labelset_obs = pd.DataFrame(
            {"set1": ["A", "B", "C"]}  # "set2" is missing
        )
        self.assertFalse(validate_labelset_values(cas, missing_labelset_obs))

        # Empty CAS annotations (should return True)
        empty_cas = {"annotations": []}
        self.assertTrue(validate_labelset_values(empty_cas, obs))

        # Empty obs DataFrame (should return False)
        empty_obs = pd.DataFrame()
        self.assertFalse(validate_labelset_values(cas, empty_obs))

    def test_infer_cas_cell_hierarchy(self):
        cas = {
            "annotations": [
                {"cell_label": "Neuron", "cell_set_accession": "CS001"},
                {
                    "cell_label": "Excitatory",
                    "cell_set_accession": "CS002",
                    "parent_cell_set_accession": "CS001",
                },
                {
                    "cell_label": "Inhibitory",
                    "cell_set_accession": "CS003",
                    "parent_cell_set_accession": "CS001",
                },
                {
                    "cell_label": "Layer5",
                    "cell_set_accession": "CS004",
                    "parent_cell_set_accession": "CS002",
                },
            ]
        }
        # Case 1: Correct hierarchy inference
        expected_hierarchy = {
            "Neuron": None,
            "Excitatory": "Neuron",
            "Inhibitory": "Neuron",
            "Layer5": "Excitatory",
        }
        self.assertEqual(infer_cas_cell_hierarchy(cas), expected_hierarchy)

        # Case 2: No parent relationships (All should be None)
        cas_no_parents = {
            "annotations": [
                {"cell_label": "Neuron", "cell_set_accession": "CS001"},
                {"cell_label": "Excitatory", "cell_set_accession": "CS002"},
                {"cell_label": "Inhibitory", "cell_set_accession": "CS003"},
            ]
        }
        expected_no_parents = {
            "Neuron": None,
            "Excitatory": None,
            "Inhibitory": None,
        }
        self.assertEqual(infer_cas_cell_hierarchy(cas_no_parents), expected_no_parents)

        # Case 3: Empty CAS Annotations (Should return an empty dictionary)
        empty_cas = {"annotations": []}
        self.assertEqual(infer_cas_cell_hierarchy(empty_cas), {})

        # Case 4: Invalid Parent Accession (Should raise KeyError)
        cas_invalid_parent = {
            "annotations": [
                {"cell_label": "Neuron", "cell_set_accession": "CS001"},
                {
                    "cell_label": "Excitatory",
                    "cell_set_accession": "CS002",
                    "parent_cell_set_accession": "INVALID_ID",
                },
            ]
        }
        with self.assertRaises(KeyError):
            infer_cas_cell_hierarchy(cas_invalid_parent)

    def test_infer_obs_cell_hierarchy(self):
        obs = pd.DataFrame(
            {
                "set1": ["A1", "A2", "A3", "A4", "B1", "B2", "B3", "C1", "C2"],
                "set2": ["A", "A", "Ax", "Ax", "B", "B", "B", "C", "C"],
                "set3": ["X", "X", "X", "X", "Y", "Y", "Y", "Y", "Y"],
            },
            index=range(9),
        )

        cas_ranks = {"set1": 0, "set2": 1, "set3": 2}

        expected_hierarchy = {
            "A1": "A",
            "A2": "A",
            "A3": "Ax",
            "A4": "Ax",
            "B1": "B",
            "B2": "B",
            "B3": "B",
            "C1": "C",
            "C2": "C",
            "A": "X",
            "Ax": "X",
            "B": "Y",
            "C": "Y",
            "X": None,
            "Y": None,
        }

        hierarchy = infer_obs_cell_hierarchy(obs, cas_ranks)
        self.assertEqual(hierarchy, expected_hierarchy)
