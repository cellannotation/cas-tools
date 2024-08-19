import unittest
from unittest.mock import patch

from cas.cas_splitter import (
    filter_and_copy_cas_entries,
    get_split_terms,
    split_cas,
    split_cas_to_file,
)


class TestCASSplitter(unittest.TestCase):
    def setUp(self):
        self.sample_cas = {
            "annotations": [
                {"cell_set_accession": "X", "labelset": "Alpha"},
                {"cell_set_accession": "Y", "labelset": "Alpha"},
                {"cell_set_accession": "Z", "labelset": "Alpha"},
                {
                    "cell_set_accession": "A",
                    "parent_cell_set_accession": "X",
                    "labelset": "Beta",
                },
                {
                    "cell_set_accession": "B",
                    "parent_cell_set_accession": "Y",
                    "labelset": "Beta",
                },
                {
                    "cell_set_accession": "BB",
                    "parent_cell_set_accession": "Y",
                    "labelset": "Beta",
                },
                {
                    "cell_set_accession": "C",
                    "parent_cell_set_accession": "Z",
                    "labelset": "Beta",
                },
                {
                    "cell_set_accession": "A1",
                    "parent_cell_set_accession": "A",
                    "labelset": "Gamma",
                },
                {
                    "cell_set_accession": "B1",
                    "parent_cell_set_accession": "B",
                    "labelset": "Gamma",
                },
                {
                    "cell_set_accession": "B2",
                    "parent_cell_set_accession": "BB",
                    "labelset": "Gamma",
                },
            ],
            "labelsets": [{"name": "Alpha"}, {"name": "Beta"}, {"name": "Gamma"}],
        }
        self.split_terms = ["X", "Y"]

    @patch("cas.cas_splitter.read_json_file")
    @patch("cas.cas_splitter.write_dict_to_json_file")
    @patch("cas.cas_splitter.split_cas")
    def test_split_cas_to_file(
        self, mock_split_cas, mock_write_dict_to_json_file, mock_read_json_file
    ):
        # Setup mocks
        mock_read_json_file.return_value = self.sample_cas
        mock_split_cas.return_value = [{"annotations": [], "labelsets": []}] * 2

        # Call function
        split_cas_to_file("dummy_path.json", self.split_terms, True)

        # Assert file read and write calls
        mock_read_json_file.assert_called_once_with("dummy_path.json")
        self.assertEqual(mock_write_dict_to_json_file.call_count, 2)

    def test_split_cas_multiple_outputs(self):
        with self.assertRaises(ValueError):
            split_cas(self.sample_cas, ["Unknown"], True)

        result = split_cas(self.sample_cas, self.split_terms, True)
        self.assertEqual(len(result), 2)
        self.assertIsInstance(result, list)

    def test_filter_and_copy_cas_entries(self):
        label_to_copy_list = ["Z", "C"]
        result = filter_and_copy_cas_entries(self.sample_cas, label_to_copy_list)
        self.assertIn("annotations", result)
        self.assertIn("labelsets", result)
        self.assertEqual(len(result["annotations"]), 2)

    def test_get_split_terms(self):
        parent_dict = {"X": ["A"], "Y": ["B", "BB"]}
        result = get_split_terms(parent_dict, "X")
        self.assertIn("A", result)
        self.assertIn("X", result)
        self.assertNotIn("B", result)
        self.assertNotIn("BB", result)
