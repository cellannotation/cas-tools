import os
import unittest
from dataclasses import asdict

import pandas as pd

from cas.file_utils import read_cas_json_file
from cas.reports import get_all_annotations

TEST_JSON = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "./test_data/siletti/Siletti_all_non_neuronal_cells_with_cids.json",
)


class ReportsTests(unittest.TestCase):
    def setUp(self):
        pd.set_option("display.max_columns", None)

    def test_annotations_listing(self):
        cas = read_cas_json_file(TEST_JSON)
        self.assertEqual(89, len(cas.annotations))

        df = cas.get_all_annotations()
        self.assertEqual(89, df.shape[0])
        self.assertEqual(18, df.shape[1])

        cas = asdict(cas)
        df = get_all_annotations(cas)
        # print(df)
        self.assertEqual(89, df.shape[0])
        self.assertEqual(18, df.shape[1])

        df = get_all_annotations(cas, show_cell_ids=True)
        self.assertEqual(89, df.shape[0])
        self.assertEqual(19, df.shape[1])

    def test_annotations_listing_with_pairs(self):
        cas = read_cas_json_file(TEST_JSON)
        self.assertEqual(89, len(cas.annotations))

        df = cas.get_all_annotations(
            labels=[
                ("supercluster_term", "Microglia"),
                ("Cluster", "Mgl_4"),
                ("Cluster", "Astro_55"),
                ("Dummy", "Dummy"),
            ]
        )
        print(df)
        self.assertEqual(3, df.shape[0])
        self.assertEqual(18, df.shape[1])

    def test_cell_ids_not_existing(self):
        cas = read_cas_json_file(TEST_JSON)
        cas = asdict(cas)

        for ann in cas["annotations"]:
            del ann["cell_ids"]

        df = get_all_annotations(cas)
        # print(df)
        self.assertEqual(89, df.shape[0])
        self.assertEqual(18, df.shape[1])
