import unittest

from cas.anndata_to_cas import (
    add_annotations_to_cas,
    add_parent_cell_hierarchy,
    generate_cas_metadata,
)


class TestAnndataToCas(unittest.TestCase):
    def setUp(self):
        self.uns = {"schema_version": "1.0"}
        self.parent_cell_look_up = {
            "A": {
                "cell_ids": {1, 2},
                "accession": "A_123",
                "parent": "P",
                "p_accession": "P_123",
                "rank": 0,
                "cell_ontology_term_id": "CL:1234567",
                "cell_ontology_term": "Test cell",
            },
            "P": {
                "cell_ids": {1, 2},
                "accession": "P_123",
                "rank": 1,
                "cell_ontology_term_id": "CL:1234567",
                "cell_ontology_term": "Test cell",
            },
        }

    def test_generate_cas_metadata(self):
        cas = generate_cas_metadata(self.uns)

        self.assertEqual(cas["cellannotation_schema_version"], "1.0")
        self.assertIn("annotations", cas)

    def test_add_annotations_to_cas(self):
        cas = {"annotations": []}
        labelset_dict = {"labelset1": {"members": {"A", "P"}, "rank": 1}}

        add_annotations_to_cas(cas, labelset_dict, self.parent_cell_look_up)

        self.assertIsInstance(cas["annotations"], list)
        self.assertEqual(len(cas["annotations"]), 2)
        self.assertEqual(len(cas["annotations"][0]), 7)
        self.assertEqual(len(cas["annotations"][1]), 7)
