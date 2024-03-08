from typing import Any, Dict
import unittest

import pandas as pd
import anndata as ad

from cas.anndata_to_cas import (
    generate_cas_annotations,
    generate_cas_labelsets,
    generate_cas_metadata,
    calculate_labelset,
    calculate_labelset_rank,
    add_parent_cell_hierarchy,
    update_parent_info,
)


class TestAnndataToCas(unittest.TestCase):
    def setUp(self):
        # Create a simple AnnData-like DataFrame for testing
        self.obs = pd.DataFrame(
            {
                "labelset1": ["A", "B", "C", "A"],
                "labelset2": ["X", "Y", "X", "Y"],
                "cell_type": ["Type1", "Type2", "Type1", "Type2"],
                "cell_type_ontology_term_id": [1, 2, 1, 2],
            }
        )
        self.uns = {"schema_version": "1.0"}

    def test_generate_cas_labelsets(self):
        cas = {"labelsets": []}
        labelset_dict = {
            "labelset1": {"members": {"member1", "member2", "member3"}, "rank": 0},
            "labelset2": {"members": {"member4", "member5"}, "rank": 1},
            "labelset3": {"members": {"member6", "member7", "member8", "member9"}, "rank": 2},
        }
        generate_cas_labelsets(cas, labelset_dict)

        # Ensure labelsets are added correctly
        self.assertIn("labelsets", cas)
        self.assertEqual(len(cas["labelsets"]), 3)
        self.assertEqual(cas["labelsets"][0]["name"], "labelset1")
        self.assertEqual(cas["labelsets"][1]["name"], "labelset2")
        self.assertEqual(cas["labelsets"][2]["name"], "labelset3")
        self.assertEqual(cas["labelsets"][0]["rank"], "0")
        self.assertEqual(cas["labelsets"][1]["rank"], "1")
        self.assertEqual(cas["labelsets"][2]["rank"], "2")

    def test_generate_cas_metadata(self):
        cas = generate_cas_metadata(self.uns)

        self.assertEqual(cas["cellannotation_schema_version"], "1.0")
        self.assertIn("annotations", cas)

    def test_calculate_labelset(self):
        labelsets = ["labelset1", "labelset2"]
        labelset_dict = calculate_labelset(self.obs, labelsets)

        # Ensure labelset dictionary is calculated correctly
        self.assertIn("labelset1", labelset_dict)
        self.assertIn("labelset2", labelset_dict)
        self.assertEqual(labelset_dict["labelset1"], {"members": {"A", "B", "C"}, "rank": "0"})
        self.assertEqual(labelset_dict["labelset2"], {"members": {"X", "Y"}, "rank": "1"})

    def test_generate_cas_annotations(self):
        example_data = pd.DataFrame({
            'feature1': [1, 2],
            'feature2': [3, 4]
        })
        obs_data = pd.DataFrame({
            'cell_type': ['type1', 'type2'],
            'cell_type_ontology_term_id': ['CTO:0001', 'CTO:0002'],
            'labelset1': ['label1', 'label2']
        }, index=['sample1', 'sample2'])

        adata = ad.AnnData(X=example_data.values, obs=obs_data)

        cas = {'annotations': []}
        include_hierarchy = True
        labelset_dict = {
            'labelset1': {
                'members': {'label1', 'label2'},
                'rank': 1
            }
        }

        result = generate_cas_annotations(adata, cas, include_hierarchy, labelset_dict)

        self.assertIsInstance(result, dict)
        self.assertIn("label1", result)
        self.assertIn("label2", result)

    def test_add_parent_cell_hierarchy(self):
        cas = {
            "annotations": [
                {"cell_label": "A"},
            ]
        }
        parent_cell_look_up = {"A": {"cell_ids": {1, 2}, "accession": "A_123", "parent": "P", "p_accession": "P_123"}}

        add_parent_cell_hierarchy(cas, parent_cell_look_up=parent_cell_look_up)

        # Ensure parent cell hierarchy information is added correctly
        self.assertIn("annotations", cas)
        for annotation in cas["annotations"]:
            self.assertIn("parent_cell_set_name", annotation)
            self.assertIn("parent_cell_set_accession", annotation)

    def test_update_parent_info(self):
        """Test updating child item with correct parent information."""
        # Setup initial child and parent dictionaries
        child = {"name": "child1", "parent": None, "p_accession": None, "parent_rank": None}
        parent = {"name": "parent1", "accession": "A001", "rank": 1}

        # Expected result after updating child with parent info
        expected = {"name": "child1", "parent": "parent1", "p_accession": "A001", "parent_rank": 1}

        # Update child with parent info
        update_parent_info(child, "parent1", parent)

        # Assert that child dictionary is updated correctly
        self.assertEqual(child, expected)

    def test_calculate_labelset_rank(self):
        input_list = ["label1", "label2", "label3"]
        expected_result = {"label1": 0, "label2": 1, "label3": 2}
        self.assertEqual(calculate_labelset_rank(input_list), expected_result)
