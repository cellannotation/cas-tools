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
                {"cell_label": "X", "labelset": "Alpha"},
                {"cell_label": "Y", "labelset": "Alpha"},
                {"cell_label": "Z", "labelset": "Alpha"},
                {"cell_label": "A", "parent_cell_set_name": "X", "labelset": "Beta"},
                {"cell_label": "B", "parent_cell_set_name": "Y", "labelset": "Beta"},
                {"cell_label": "BB", "parent_cell_set_name": "Y", "labelset": "Beta"},
                {"cell_label": "C", "parent_cell_set_name": "Z", "labelset": "Beta"},
                {"cell_label": "A1", "parent_cell_set_name": "A", "labelset": "Gamma"},
                {"cell_label": "B1", "parent_cell_set_name": "B", "labelset": "Gamma"},
                {"cell_label": "B2", "parent_cell_set_name": "BB", "labelset": "Gamma"},
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
