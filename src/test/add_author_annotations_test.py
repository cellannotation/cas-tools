import unittest

import numpy as np
import pandas as pd

from cas.add_author_annotations import (
    add_author_annotations,
    dataframe_to_dict,
    validate_columns,
    validate_values,
)


class TestAuthorAnnotations(unittest.TestCase):
    def setUp(self):
        # Create a dummy CAS JSON-like dictionary
        self.cas = {
            "annotations": [
                {
                    "labelset": "set1",
                    "cell_label": "101",
                    "cell_set_accession": "c_1",
                    "author_annotation_fields": {},
                },
                {
                    "labelset": "set2",
                    "cell_label": "202",
                    "cell_set_accession": "c_2",
                    "author_annotation_fields": {},
                },
            ]
        }

        data = {
            "cell_set_accession": ["c_1", "c_2"],
            "labelset": ["set1", "set2"],
            "cell_label": ["101", "202"],
            "extra_data": ["value1", "value2"],
        }
        self.df = pd.DataFrame(data)

    def test_validate_columns(self):
        # Test to ensure validate_columns raises an error if a column is missing
        with self.assertRaises(ValueError):
            validate_columns(self.df, "nonexistent_column", None)

        # Test successful validation
        try:
            validate_columns(self.df, "cell_set_accession", ["labelset", "extra_data"])
        except ValueError:
            self.fail("validate_columns() raised ValueError unexpectedly!")

    def test_single_column_filter(self):
        # Testing with a single column filter
        result = dataframe_to_dict(
            self.df, "cell_set_accession", "set1", "c_1", ["extra_data"]
        )
        self.assertEqual(result, {"extra_data": "value1"})

    def test_multi_column_filter(self):
        # Testing with multiple column filters
        result = dataframe_to_dict(
            self.df, ["labelset", "cell_label"], "set1", "101", ["extra_data"]
        )
        self.assertEqual(result, {"extra_data": "value1"})

    def test_without_specified_columns(self):
        # Testing without specifying columns - should return all columns
        result = dataframe_to_dict(self.df, "cell_set_accession", "set1", "c_1", None)
        self.assertEqual(
            result, {"labelset": "set1", "cell_label": "101", "extra_data": "value1"}
        )

    def test_unpacking_single_item_lists(self):
        # Ensure single-item lists are unpacked to single values
        result = dataframe_to_dict(
            self.df, "cell_set_accession", "set1", "c_1", ["extra_data"]
        )
        self.assertEqual(result, {"extra_data": "value1"})

    def test_invalid_join_column(self):
        # Testing with an invalid join column type
        with self.assertRaises(ValueError):
            dataframe_to_dict(
                self.df, 123, "set1", "c_1", ["extra_data"]
            )  # join_column is neither str nor list

    def test_no_matching_rows(self):
        # Testing filter that results in no matching rows
        result = dataframe_to_dict(
            self.df, "cell_set_accession", "set1", "999", ["extra_data"]
        )
        self.assertEqual(result, {"extra_data": []})

    def test_add_author_annotations(self):
        # Test that annotations are correctly added to the CAS dictionary
        updated_cas = add_author_annotations(
            self.cas, self.df, "cell_set_accession", ["extra_data"]
        )
        self.assertIn(
            "value1",
            updated_cas["annotations"][0]["author_annotation_fields"]["extra_data"],
        )
        self.assertIn(
            "value2",
            updated_cas["annotations"][1]["author_annotation_fields"]["extra_data"],
        )

    def test_validate_values(self):
        # Test where DataFrame keys perfectly match the CAS data
        try:
            validate_values(self.df, "cell_set_accession", self.cas)
        except ValueError as e:
            self.fail(f"validate_values raised an unexpected ValueError: {str(e)}")

    def test_validate_values_with_extra_keys(self):
        # Modify DataFrame to include an extra key not present in CAS
        df_with_extra = self.df.copy()
        df_with_extra.loc[2] = ["c_3", "set3", "303", "value3"]

        with self.assertRaises(ValueError) as context:
            validate_values(df_with_extra, "cell_set_accession", self.cas)

        self.assertIn(
            "Extra keys in DataFrame that are not in CAS data: ['c_3']",
            str(context.exception),
        )

    def test_validate_values_with_empty_rows(self):
        data = {
            "cell_set_accession": ["c_1", "c_2"],
            "labelset": ["set1", "set2"],
            "cell_label": ["101", "202"],
            "extra_data": ["value1", np.NaN],
        }
        df = pd.DataFrame(data)
        print(df)
        updated_cas = add_author_annotations(
            self.cas, df, "cell_set_accession", ["extra_data"]
        )
        self.assertIsNone(updated_cas["annotations"][1]["author_annotation_fields"]["extra_data"])
