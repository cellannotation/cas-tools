import unittest
import warnings
from test.spreadsheet_to_cas_test import generate_mock_dataset

import pandas as pd

from cas.utils.conversion_utils import (
    add_labelsets_to_cas,
    add_parent_cell_hierarchy,
    add_parent_hierarchy_to_annotations,
    calculate_labelset,
    calculate_labelset_rank,
    get_cell_ids,
    get_cl_annotations_from_anndata,
    update_parent_info,
)


class TestConversionUtils(unittest.TestCase):
    def setUp(self):
        # Create a simple AnnData-like DataFrame for testing
        self.obs = pd.DataFrame(
            {
                "labelset1": ["A", "B", "C", "A"],
                "labelset2": ["X", "Y", "X", "Y"],
                "cell_type": ["Type1", "Type2", "Type1", "Type2"],
                "cell_type_ontology_term_id": ["1", "2", "1", "2"],
            }
        )
        self.parent_cell_look_up = {
            "labelset1:A": {
                "cell_ids": {1, 2},
                "accession": "A_123",
                "parent": "P",
                "p_accession": "P_123",
                "rank": 0,
                "cell_ontology_term_id": "CL:1234567",
                "cell_ontology_term": "Test cell",
            },
            "labelset3:P": {
                "cell_ids": {1, 2},
                "accession": "P_123",
                "rank": 1,
                "cell_ontology_term_id": "CL:1234567",
                "cell_ontology_term": "Test cell",
            },
        }

    def test_calculate_labelset_rank(self):
        # Test with an empty list
        result_empty = calculate_labelset_rank([])
        self.assertEqual(result_empty, {})

        # Test with a non-empty list
        input_list = ["item1", "item2", "item3"]
        result_non_empty = calculate_labelset_rank(input_list)
        expected_result_non_empty = {"item1": 0, "item2": 1, "item3": 2}
        self.assertEqual(result_non_empty, expected_result_non_empty)

    def test_calculate_labelset(self):
        labelsets = ["labelset1", "labelset2"]
        labelset_dict = calculate_labelset(self.obs, labelsets)

        self.assertIn("labelset1", labelset_dict)
        self.assertIn("labelset2", labelset_dict)
        self.assertEqual(
            labelset_dict["labelset1"], {"members": {"A", "B", "C"}, "rank": "0"}
        )
        self.assertEqual(
            labelset_dict["labelset2"], {"members": {"X", "Y"}, "rank": "1"}
        )

    def test_update_parent_info(self):
        """Test updating child item with correct parent information."""
        child = {
            "name": "child1",
            "parent": None,
            "p_accession": None,
            "parent_rank": None,
        }
        parent = {"name": "parent1", "accession": "A001", "rank": 1}

        expected = {
            "name": "child1",
            "parent": "parent1",
            "p_accession": "A001",
            "parent_rank": 1,
        }

        update_parent_info(child, "parent1", parent)

        self.assertEqual(child, expected)

    def test_add_labelsets_to_cas(self):
        cas = {"labelsets": []}
        labelset_dict = {
            "labelset1": {"members": {"member1", "member2", "member3"}, "rank": 0},
            "labelset2": {"members": {"member4", "member5"}, "rank": 1},
            "labelset3": {
                "members": {"member6", "member7", "member8", "member9"},
                "rank": 2,
            },
        }
        add_labelsets_to_cas(cas, labelset_dict)

        self.assertIn("labelsets", cas)
        self.assertEqual(len(cas["labelsets"]), 3)
        self.assertEqual(cas["labelsets"][0]["name"], "labelset1")
        self.assertEqual(cas["labelsets"][1]["name"], "labelset2")
        self.assertEqual(cas["labelsets"][2]["name"], "labelset3")
        self.assertEqual(cas["labelsets"][0]["rank"], "0")
        self.assertEqual(cas["labelsets"][1]["rank"], "1")
        self.assertEqual(cas["labelsets"][2]["rank"], "2")

    def test_get_cell_ids(self):
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=UserWarning, module="anndata._core.anndata"
            )
            mock_dataset = generate_mock_dataset()

            # Test get_cell_ids function
            cell_ids = get_cell_ids(
                mock_dataset.obs, "cell_type", "CD4-positive helper T cell"
            )
            self.assertEqual(cell_ids, ["cell_2", "cell_3"])
            cell_ids = get_cell_ids(mock_dataset.obs, "annotation_broad", "T CD4+")
            self.assertEqual(cell_ids, ["cell_2", "cell_3", "cell_5"])

    def test_get_cl_annotations_from_anndata(self):
        cell_label = "Type1"
        expected_output = ("1", "Type1")

        result = get_cl_annotations_from_anndata(self.obs, "cell_type", cell_label)

        self.assertEqual(
            result,
            expected_output,
            "The function did not return the expected cell ontology term ID and cell type for the given cell label.",
        )

    def test_add_parent_cell_hierarchy(self):
        cas = {
            "annotations": [
                {
                    "cell_label": "A",
                    "cell_ontology_term_id": "CL:1234567",
                    "cell_ontology_term": "Test cell",
                },
                {
                    "cell_label": "P",
                    "cell_ontology_term_id": "CL:1234567",
                    "cell_ontology_term": "Test cell",
                },
            ]
        }

        add_parent_cell_hierarchy(parent_cell_look_up=self.parent_cell_look_up)

        # Ensure parent cell hierarchy information is added correctly
        self.assertIn("parent", self.parent_cell_look_up["labelset1:A"])
        self.assertIn("p_accession", self.parent_cell_look_up["labelset1:A"])
        self.assertEqual(self.parent_cell_look_up["labelset1:A"]["parent"], "P")
        self.assertEqual(
            self.parent_cell_look_up["labelset1:A"]["p_accession"], "P_123"
        )

    def test_add_parent_hierarchy(self):
        cas = {
            "annotations": [
                {
                    "cell_label": "A",
                    "labelset": "labelset1",
                    "cell_ontology_term_id": "CL:1234567",
                    "cell_ontology_term": "Test cell",
                }
            ]
        }

        expected_annotations = [
            {
                "cell_label": "A",
                "cell_ontology_term": "Test cell",
                "cell_ontology_term_id": "CL:1234567",
                "labelset": "labelset1",
                "parent_cell_set_accession": "P_123",
                # "parent_cell_set_name": "P",
            }
        ]

        # Execute the function to test
        add_parent_hierarchy_to_annotations(cas, self.parent_cell_look_up)

        # Assert that the annotations were updated correctly
        self.assertEqual(cas["annotations"], expected_annotations)
