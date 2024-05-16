import os
import unittest

from cas.ingest.ingest_user_table import ingest_data

RAW_DATA = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "./test_data/nhp_basal_ganglia/AIT115_annotation_sheet.tsv",
)
TEST_CONFIG = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "./test_data/nhp_basal_ganglia/test_config.yaml",
)
OUT_FILE = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "./test_data/nhp_basal_ganglia/test_result.json",
)


class CellTypeAnnotationTests(unittest.TestCase):
    def setUp(self):
        if os.path.isfile(OUT_FILE):
            os.remove(OUT_FILE)

    @classmethod
    def tearDownClass(cls):
        if os.path.isfile(OUT_FILE):
            os.remove(OUT_FILE)

    def test_data_formatting(self):
        result = ingest_data(RAW_DATA, TEST_CONFIG, OUT_FILE, "json", True)
        # print(result)
        self.assertTrue(result)
        self.assertTrue("author_name" in result)
        self.assertEqual("Test User", result["author_name"])

        self.assertTrue("cellannotation_schema_version" in result)
        self.assertTrue("." in result["cellannotation_schema_version"])

        self.assertTrue("labelsets" in result)
        self.assertEqual(4, len(result["labelsets"]))
        # print(result["labelsets"])

        self.assertTrue("annotations" in result)
        self.assertEqual(354, len(result["annotations"]))
        # print(result["annotations"][:10])

        test_annotation = [
            x for x in result["annotations"] if x["cell_label"] == "1_MSN"
        ][0]
        self.assertTrue("marker_gene_evidence" in test_annotation)
        self.assertEqual(3, len(test_annotation["marker_gene_evidence"]))
        self.assertTrue("EPYC" in test_annotation["marker_gene_evidence"])
        self.assertTrue("RELN" in test_annotation["marker_gene_evidence"])
        self.assertTrue("GULP1" in test_annotation["marker_gene_evidence"])

        self.assertTrue("author_annotation_fields" in test_annotation)
        self.assertEqual(6, len(test_annotation["author_annotation_fields"]))
