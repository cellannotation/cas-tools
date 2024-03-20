import unittest

import numpy as np
import pandas as pd

from cas.abc_cas_converter import (
    add_annotations,
    add_labelsets,
    calculate_order_mapping,
    generate_cat_dataframe,
    generate_catset_dataframe,
    init_metadata,
    validate_dataframe_columns,
)


class TestABCConverter(unittest.TestCase):

    def setUp(self):
        # Sample data for testing
        self.cas = init_metadata()
        self.cat = pd.DataFrame(
            {
                "cluster_annotation_term_set_name": ["labelset1", "labelset2"],
                "name": ["cell1", "cell2"],
                "label": ["label1", "label2"],
                "parent_term_label": ["parent1", "parent2"],
                "order": [1, 2],
            }
        )
        self.cat_set = pd.DataFrame(
            {
                "name": ["labelset1", "labelset2"],
                "description": ["description1", "description2"],
                "order": [0, 1],
            }
        )

    def test_validate_dataframe_columns(self):
        test_df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        required_columns = ["C"]
        with self.assertRaises(ValueError) as context:
            validate_dataframe_columns(test_df, required_columns)
        self.assertTrue("Missing columns in DataFrame: C" in str(context.exception))

        # Test case where all required columns are missing
        df_all_missing_columns = pd.DataFrame()
        required_columns = ["A", "B", "C"]
        with self.assertRaises(ValueError) as context:
            validate_dataframe_columns(df_all_missing_columns, required_columns)
        self.assertTrue(
            "Missing columns in DataFrame: A, B, C" in str(context.exception)
        )

    def test_generate_catset_dataframe(self):
        cas = {
            "labelsets": [
                {"name": "Label1", "description": "Desc1", "rank": 1},
                {"name": "Label2", "description": "Desc2", "rank": 2},
            ]
        }
        expected_df = pd.DataFrame(
            {
                "name": ["Label1", "Label2"],
                "description": ["Desc1", "Desc2"],
                "order": [1, 2],
            }
        )
        pd.testing.assert_frame_equal(generate_catset_dataframe(cas), expected_df)

    def test_generate_cat_dataframe(self):
        cas = {
            "annotations": [
                {
                    "cell_label": "Cell1",
                    "cell_set_accession": 1,
                    "parent_cell_set_accession": 0,
                },
                {
                    "cell_label": "Cell2",
                    "cell_set_accession": 2,
                    "parent_cell_set_accession": 1,
                },
            ]
        }
        expected_df = pd.DataFrame(
            {
                "label": [1, 2],
                "name": ["Cell1", "Cell2"],
                "cluster_annotation_term_set_label": [np.nan, np.nan],
                "parent_term_label": [0, 1],
                "parent_term_set_label": [np.nan, np.nan],
                "term_set_order": [np.nan, np.nan],
                "term_order": [np.nan, np.nan],
                "cluster_annotation_term_set_name": [np.nan, np.nan],
            }
        )
        pd.testing.assert_frame_equal(generate_cat_dataframe(cas), expected_df)

    def test_calculate_order_mapping(self):
        order_values_series = pd.Series([3, 1, 4, 2, 0])
        expected_mapping = {4: 0, 3: 1, 2: 2, 1: 3}
        self.assertEqual(calculate_order_mapping(order_values_series), expected_mapping)

    def test_add_annotations(self):
        add_annotations(self.cas, self.cat)
        self.assertEqual(len(self.cas["annotations"]), 2)
        self.assertEqual(self.cas["annotations"][0]["cell_label"], "cell1")
        self.assertEqual(self.cas["annotations"][1]["parent_cell_set_accession"], "parent2")

    def test_add_labelsets(self):
        add_labelsets(self.cas, self.cat_set)
        self.assertEqual(len(self.cas["labelsets"]), 2)
        self.assertEqual(self.cas["labelsets"][0]["description"], "description1")
        self.assertEqual(self.cas["labelsets"][1]["rank"], 0)

    def test_init_metadata(self):
        metadata = init_metadata()
        self.assertIsInstance(metadata, dict)
        self.assertTrue("annotations" in metadata)
        self.assertEqual(len(metadata["labelsets"]), 0)
        self.assertEqual(len(metadata["annotations"]), 0)
