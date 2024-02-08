import unittest
from cas.anndata_to_cas import (
    generate_cas_annotations,
    generate_cas_labelsets,
    generate_cas_metadata,
    calculate_labelset,
    add_parent_cell_hierarchy,
)
from cas.file_utils import read_anndata_file
import pandas as pd


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

    def test_generate_cas_annotations(self):
        cas = generate_cas_metadata(self.uns)

        # Ensure CAS annotations are added correctly
        self.assertIn("annotations", cas)

    def test_generate_cas_labelsets(self):
        cas = {"labelsets": []}
        labelsets = ["labelset1", "labelset2"]
        generate_cas_labelsets(cas, labelsets)

        # Ensure labelsets are added correctly
        self.assertIn("labelsets", cas)
        self.assertEqual(len(cas["labelsets"]), 2)
        self.assertEqual(cas["labelsets"][0]["name"], "labelset1")
        self.assertEqual(cas["labelsets"][1]["name"], "labelset2")

    def test_generate_cas_metadata(self):
        cas = generate_cas_metadata(self.uns)

        # Ensure metadata is generated correctly
        self.assertEqual(cas["cellannotation_schema_version"], "1.0")

    def test_calculate_labelset(self):
        labelsets = ["labelset1", "labelset2"]
        labelset_dict = calculate_labelset(self.obs, labelsets)

        # Ensure labelset dictionary is calculated correctly
        self.assertIn("labelset1", labelset_dict)
        self.assertIn("labelset2", labelset_dict)
        self.assertEqual(labelset_dict["labelset1"], {"A", "B", "C"})
        self.assertEqual(labelset_dict["labelset2"], {"X", "Y"})

    def test_add_parent_cell_hierarchy(self):
        cas = {
            "annotations": [
                {"cell_label": "A"},
            ]
        }
        parent_cell_look_up = {"A": {"cell_ids": {1, 2}, "accession": "A_123", "parent": "P", "p_accession": "P_123"}}

        add_parent_cell_hierarchy(cas, include_hierarchy=True, parent_cell_look_up=parent_cell_look_up)

        # Ensure parent cell hierarchy information is added correctly
        self.assertIn("annotations", cas)
        for annotation in cas["annotations"]:
            self.assertIn("parent_cell_set_name", annotation)
            self.assertIn("parent_cell_set_accession", annotation)