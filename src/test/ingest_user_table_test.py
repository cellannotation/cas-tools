import os

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
import json
import unittest

from cas.ingest.ingest_user_table import ingest_data, ingest_user_data

RAW_DATA = os.path.join(
    CURRENT_DIR,
    "./test_data/nhp_basal_ganglia/AIT115_annotation_sheet.tsv",
)
TEST_CONFIG = os.path.join(
    CURRENT_DIR,
    "./test_data/nhp_basal_ganglia/test_config.yaml",
)
RAW_DATAv2 = os.path.join(
    CURRENT_DIR,
    "./test_data/nhp_basal_ganglia/v2/AIT117_joint_annotation_sheet.tsv",
)
TEST_CONFIGv2 = os.path.join(
    CURRENT_DIR,
    "./test_data/nhp_basal_ganglia/v2/ingestion_config.yaml",
)
OUT_FILE = os.path.join(
    CURRENT_DIR,
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

        self.assertTrue("title" in result)
        self.assertEqual("NHP Basal Ganglia", result["title"])

        self.assertTrue("cellannotation_schema_version" in result)
        self.assertTrue("." in result["cellannotation_schema_version"])

        self.assertTrue("labelsets" in result)
        self.assertEqual(4, len(result["labelsets"]))
        # print(result["labelsets"])

        self.assertTrue("annotations" in result)
        self.assertEqual(363, len(result["annotations"]))
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

    def test_data_formatting_nhp_v2(self):
        result = ingest_user_data(
            RAW_DATAv2,
            TEST_CONFIGv2,
        )

        self.assertTrue(result)
        self.assertTrue(result.author_name)
        self.assertEqual("Nelson Johansen", result.author_name)

        self.assertIsNotNone(result.title)
        self.assertEqual("NHP Basal Ganglia taxonomy", result.title)

        self.assertIsNotNone(result.annotations)
        # self.assertEqual(354, len(result.annotations))
        # print(result["annotations"][:10])

        test_annotation = [x for x in result.annotations if x.cell_label == "76_IN"][0]
        # print(test_annotation.get("parent_cell_set_accession"))
        self.assertIsNone(test_annotation.parent_cell_set_accession)
        self.assertEqual("WDR49-ADAM12", test_annotation.parent_cell_set_name)
        parent_annotation = [
            x
            for x in result.annotations
            if x.cell_label == test_annotation.parent_cell_set_name
        ][0]
        self.assertEqual("WDR49-ADAM12", parent_annotation.cell_label)
        self.assertEqual("supertype", parent_annotation.labelset)

        # 10_NN is a child of a subclass (not supertype as usual)
        test_annotation = [x for x in result.annotations if x.cell_label == "10_NN"][0]
        print(test_annotation.parent_cell_set_name)
        parent_annotation = [
            x
            for x in result.annotations
            if x.cell_label == test_annotation.parent_cell_set_name
        ][0]
        self.assertEqual("Astrocytes", parent_annotation.cell_label)
        self.assertEqual("subclass", parent_annotation.labelset)

        # 164_IN is a child of a class (not supertype as usual)
        test_annotation = [x for x in result.annotations if x.cell_label == "164_IN"][0]
        parent_annotation = [
            x
            for x in result.annotations
            if x.cell_label == test_annotation.parent_cell_set_name
        ][0]
        self.assertEqual("CN MGE GABA", parent_annotation.cell_label)
        self.assertEqual("class", parent_annotation.labelset)

        # 11_NN is a child of a Endothelial (supertype) which is child of Endothelial (subclass)
        test_annotation = [x for x in result.annotations if x.cell_label == "11_NN"][0]
        print(test_annotation.parent_cell_set_name)
        annotations_with_parent_name = [
            x
            for x in result.annotations
            if x.cell_label == test_annotation.parent_cell_set_name
        ]
        self.assertEqual(2, len(annotations_with_parent_name))  # 2 annotations with the same name

        labelsets = {labelset.name: labelset for labelset in result.labelsets}
        parent_annotation = min(annotations_with_parent_name, key=lambda d: labelsets[d.labelset].rank)
        self.assertEqual("Endothelial", parent_annotation.cell_label)
        self.assertEqual("supertype", parent_annotation.labelset)
        self.assertEqual("Endothelial", parent_annotation.parent_cell_set_name)

        grand_parent_annotation = max(annotations_with_parent_name, key=lambda d: labelsets[d.labelset].rank)
        self.assertEqual("Endothelial", grand_parent_annotation.cell_label)
        self.assertEqual("subclass", grand_parent_annotation.labelset)
        self.assertEqual("Vascular", grand_parent_annotation.parent_cell_set_name)


    def test_data_formatting_wmb(self):
        result = ingest_user_data(
            os.path.join(
                CURRENT_DIR,
                "./test_data/wmb/wmb_class_29_annotation.tsv",
            ),
            os.path.join(
                CURRENT_DIR,
                "./test_data/wmb/wmb_ingestion_config.yaml",
            ),
            True
        )

        self.assertTrue(result)
        self.assertTrue(result.author_name)
        self.assertEqual("Hongkui Zeng", result.author_name)

        self.assertIsNotNone(result.title)
        self.assertEqual("Whole Mouse Brain taxonomy", result.title)

        self.assertIsNotNone(result.annotations)

        test_annotation = [x for x in result.annotations if x.cell_label == "5201 CB Granule Glut_2"][0]
        # print(test_annotation.get("parent_cell_set_accession"))
        self.assertIsNotNone(test_annotation.parent_cell_set_accession)
        self.assertEqual("1155 CB Granule Glut_2", test_annotation.parent_cell_set_name)
        parent_annotation = [
            x
            for x in result.annotations
            if x.cell_label == test_annotation.parent_cell_set_name
        ][0]
        self.assertEqual("1155 CB Granule Glut_2", parent_annotation.cell_label)
        self.assertEqual("supertype", parent_annotation.labelset)
        self.assertEqual(test_annotation.parent_cell_set_accession, parent_annotation.cell_set_accession)

        print(result.labelsets[0].rank)
        print(type(result.labelsets[0].rank))
        self.assertTrue(type(result.labelsets[0].rank) == int)

        data = result.as_dictionary()
        print(json.dumps(data, indent=2))

    def test_data_formatting_bg_ait119_macaque(self):
        result = ingest_user_data(os.path.join(CURRENT_DIR, "./test_data/nhp_basal_ganglia/v3_ait119/Macaque_AIBS_AIT11-9_anno_table.tsv"),
                             os.path.join(CURRENT_DIR, "./test_data/nhp_basal_ganglia/v3_ait119/macaque_ingestion_config.yaml"),
                             True, "AIT119")

        self.assertTrue(result)
        self.assertTrue(result.author_name)
        self.assertEqual("Nelson Johansen", result.author_name)

        self.assertIsNotNone(result.title)
        self.assertEqual("NHP Basal Ganglia AIT119 taxonomy", result.title)

        self.assertIsNotNone(result.annotations)
        self.assertEqual(544, len(result.annotations))
        # print(result["annotations"][:10])

        test_annotation = [x for x in result.annotations if x.cell_label == "Cluster_440"][0]
        # print(test_annotation.get("parent_cell_set_accession"))
        self.assertEqual('AIT119_440', test_annotation.cell_set_accession)
        self.assertEqual('AIT119_516', test_annotation.parent_cell_set_accession)
        self.assertEqual("BAM", test_annotation.parent_cell_set_name)
        parent_annotation = [
            x
            for x in result.annotations
            if x.cell_label == test_annotation.parent_cell_set_name
        ][0]
        self.assertEqual("BAM", parent_annotation.cell_label)
        self.assertEqual("Group", parent_annotation.labelset)
        self.assertEqual('AIT119_516', parent_annotation.cell_set_accession)
        self.assertEqual('AIT119_515', parent_annotation.parent_cell_set_accession)


    # def test_data_formatting_ait119_macaque(self):
    #     result = ingest_data("/Users/hk9/Downloads/Macaque_AIBS_AIT11-9_anno_table.tsv",
    #                          "/Users/hk9/Downloads/macaque_ingestion_config.yaml",
    #                          "/Users/hk9/Downloads/Macaque_AIBS_AIT11-9.json",
    #                          "json", True, True)
    #
    # def test_data_formatting_ait119_human(self):
    #     result = ingest_data("/Users/hk9/Downloads/Human_AIBS_AIT19-5_anno_table.tsv",
    #                          "/Users/hk9/Downloads/human_ingestion_config.yaml",
    #                          "/Users/hk9/Downloads/Human_AIBS_AIT19-5.json",
    #                          "json", True, True)