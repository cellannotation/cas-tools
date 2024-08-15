import os
import shutil
import unittest

from cas.accession.hash_accession_manager import is_hash_accession
from cas.file_utils import read_cas_json_file, read_csv_to_dict
from cas.flatten_data_to_tables import serialize_to_tables
from cas.ingest.ingest_user_table import ingest_user_data

RAW_DATA = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "./test_data/nhp_basal_ganglia/AIT115_annotation_sheet.tsv",
)
TEST_CONFIG = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "./test_data/nhp_basal_ganglia/test_config.yaml",
)

TEST_JSON = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "./test_data/siletti/Siletti_all_non_neuronal_cells_with_cids.json",
)

TEST_JSON2 = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "./test_data/jorstad_mtg.json",
)

TEST_JSON_WMB = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "./test_data/WMB.json",
)

OUT_FOLDER = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "./test_data/table_out/"
)

RAW_DATAv2 = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "./test_data/nhp_basal_ganglia/v2/AIT117_joint_annotation_sheet.tsv",
)

TEST_CONFIGv2 = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "./test_data/nhp_basal_ganglia/v2/ingestion_config.yaml",
)


class TabularSerialisationTests(unittest.TestCase):
    def setUp(self):
        if not os.path.exists(OUT_FOLDER):
            os.makedirs(OUT_FOLDER)

        test_folder = os.listdir(OUT_FOLDER)
        for item in test_folder:
            if item.endswith(".tsv"):
                os.remove(os.path.join(OUT_FOLDER, item))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(OUT_FOLDER)

    def test_annotation_table(self):
        cta = ingest_user_data(RAW_DATA, TEST_CONFIG)
        tables = serialize_to_tables(
            cta, "Test_table", OUT_FOLDER, {"accession_id_prefix": "TST_"}
        )

        annotation_table_path = os.path.join(OUT_FOLDER, "annotation.tsv")
        self.assertEqual(annotation_table_path, tables[0])
        self.assertTrue(os.path.isfile(annotation_table_path))

        headers, records = read_csv_to_dict(
            annotation_table_path, id_column_name="cell_set_accession", delimiter="\t"
        )
        self.assertEqual(354, len(records))

        cluster_1 = records["TST_1"]
        self.assertEqual("1_MSN", cluster_1["cell_label"])
        self.assertEqual("TST_288", cluster_1["parent_cell_set_accession"])
        self.assertEqual("D1-Matrix", cluster_1["parent_cell_set_name"])
        self.assertEqual("cluster", cluster_1["labelset"])
        self.assertEqual("EPYC|RELN|GULP1", cluster_1["marker_gene_evidence"])
        self.assertEqual("PuR(0.52) | CaH(0.39)", cluster_1["region.info _Frequency_"])
        # self.assertEqual("", cluster_1["cell_ids"])

        cluster_288 = records["TST_288"]
        self.assertEqual("D1-Matrix", cluster_288["cell_label"])
        self.assertEqual("TST_327", cluster_288["parent_cell_set_accession"])
        self.assertEqual("D1-MSN", cluster_288["parent_cell_set_name"])
        self.assertEqual("level 3 (subclass)", cluster_288["labelset"])

        cluster_327 = records["TST_327"]
        self.assertEqual("D1-MSN", cluster_327["cell_label"])
        self.assertEqual("TST_353", cluster_327["parent_cell_set_accession"])
        self.assertEqual("MSN", cluster_327["parent_cell_set_name"])
        self.assertEqual("level 2 (neighborhood)", cluster_327["labelset"])

        cluster_353 = records["TST_353"]
        self.assertEqual("MSN", cluster_353["cell_label"])
        self.assertEqual("", cluster_353["parent_cell_set_accession"])
        self.assertEqual("", cluster_353["parent_cell_set_name"])
        self.assertEqual("level1 (class)", cluster_353["labelset"])

    def test_annotation_table_nhp_v2(self):
        cta = ingest_user_data(RAW_DATAv2, TEST_CONFIGv2)
        tables = serialize_to_tables(
            cta, "Test_table", OUT_FOLDER, {"accession_id_prefix": "AIT117_"}
        )

        annotation_table_path = os.path.join(OUT_FOLDER, "annotation.tsv")
        self.assertEqual(annotation_table_path, tables[0])
        self.assertTrue(os.path.isfile(annotation_table_path))

        headers, records = read_csv_to_dict(
            annotation_table_path, id_column_name="cell_set_accession", delimiter="\t"
        )
        # self.assertEqual(361, len(records))

        cluster = [
            records[rec_id]
            for rec_id in records
            if records[rec_id]["cell_label"] == "1_MSN"
        ][0]
        self.assertEqual("AIT117_1", cluster["cell_set_accession"])
        self.assertEqual("1_MSN", cluster["cell_label"])
        # self.assertEqual("AIT117_278", cluster["parent_cell_set_accession"])
        self.assertEqual("D1-Matrix", cluster["parent_cell_set_name"])
        self.assertEqual("cluster", cluster["labelset"])

        # 10_NN is a child of a subclass (not supertype as usual)
        cluster = [
            records[rec_id]
            for rec_id in records
            if records[rec_id]["cell_label"] == "10_NN"
        ][0]
        # self.assertEqual("AIT117_110", cluster["cell_set_accession"])
        self.assertEqual("10_NN", cluster["cell_label"])
        # self.assertEqual("AIT117_300", cluster["parent_cell_set_accession"])
        self.assertEqual("Astrocytes", cluster["parent_cell_set_name"])
        self.assertEqual("cluster", cluster["labelset"])
        parent = [
            records[rec_id]
            for rec_id in records
            if records[rec_id]["cell_label"] == "Astrocytes"
        ][0]
        self.assertEqual("Astrocytes", parent["cell_label"])
        self.assertEqual("subclass", parent["labelset"])

        # 164_IN is a child of a class (not supertype as usual)
        cluster = [
            records[rec_id]
            for rec_id in records
            if records[rec_id]["cell_label"] == "164_IN"
        ][0]
        # self.assertEqual("AIT117_110", cluster["cell_set_accession"])
        self.assertEqual("164_IN", cluster["cell_label"])
        # self.assertEqual("AIT117_300", cluster["parent_cell_set_accession"])
        self.assertEqual("CN MGE GABA", cluster["parent_cell_set_name"])
        self.assertEqual("cluster", cluster["labelset"])
        parent = [
            records[rec_id]
            for rec_id in records
            if records[rec_id]["cell_label"] == "CN MGE GABA"
        ][0]
        self.assertEqual("CN MGE GABA", parent["cell_label"])
        self.assertEqual("class", parent["labelset"])

    def test_labelset_table(self):
        cta = ingest_user_data(RAW_DATA, TEST_CONFIG)
        tables = serialize_to_tables(
            cta, "Test_table", OUT_FOLDER, {"accession_id_prefix": "TST_"}
        )

        table_path = os.path.join(OUT_FOLDER, "labelset.tsv")
        self.assertEqual(table_path, tables[1])
        self.assertTrue(os.path.isfile(table_path))

        headers, records = read_csv_to_dict(
            table_path, id_column_name="name", delimiter="\t"
        )
        self.assertEqual(4, len(records))

        cluster_ls = records["cluster"]
        self.assertEqual("0", cluster_ls["rank"])
        self.assertEqual("", cluster_ls["annotation_method"])
        self.assertEqual("", cluster_ls["automated_annotation_algorithm_name"])

        self.assertEqual("1", records["level 3 (subclass)"]["rank"])
        self.assertEqual("2", records["level 2 (neighborhood)"]["rank"])
        self.assertEqual("3", records["level1 (class)"]["rank"])

    def test_metadata_table(self):
        cta = ingest_user_data(RAW_DATA, TEST_CONFIG)
        tables = serialize_to_tables(
            cta, "Test_table", OUT_FOLDER, {"accession_id_prefix": "TST_"}
        )

        table_path = os.path.join(OUT_FOLDER, "metadata.tsv")
        self.assertEqual(table_path, tables[2])
        self.assertTrue(os.path.isfile(table_path))

        headers, records = read_csv_to_dict(
            table_path, generated_ids=True, delimiter="\t"
        )
        self.assertEqual(1, len(records))

        self.assertEqual("Test User", records[1]["author_name"])
        self.assertTrue("." in records[1]["cellannotation_schema_version"])
        self.assertEqual("", records[1]["cellannotation_version"])

    def test_annotation_transfer_table(self):
        cta = ingest_user_data(RAW_DATA, TEST_CONFIG)
        tables = serialize_to_tables(
            cta, "Test_table", OUT_FOLDER, {"accession_id_prefix": "TST_"}
        )

        table_path = os.path.join(OUT_FOLDER, "annotation_transfer.tsv")
        self.assertEqual(table_path, tables[3])
        self.assertTrue(os.path.isfile(table_path))

        headers, records = read_csv_to_dict(
            table_path, generated_ids=True, delimiter="\t"
        )
        self.assertEqual(1, len(records))

        self.assertEqual("", records[1]["target_node_accession"])
        self.assertEqual("", records[1]["transferred_cell_label"])
        self.assertEqual("", records[1]["source_taxonomy"])
        self.assertEqual("", records[1]["source_node_accession"])

    def test_review_table(self):
        cta = ingest_user_data(RAW_DATA, TEST_CONFIG)
        tables = serialize_to_tables(
            cta, "Test_table", OUT_FOLDER, {"accession_id_prefix": "TST_"}
        )

        table_path = os.path.join(OUT_FOLDER, "review.tsv")
        self.assertEqual(table_path, tables[4])
        self.assertTrue(os.path.isfile(table_path))

        headers, records = read_csv_to_dict(
            table_path, generated_ids=True, delimiter="\t"
        )
        self.assertEqual(0, len(records))
        self.assertEqual(5, len(headers))
        self.assertEqual(
            ["target_node_accession", "datestamp", "reviewer", "review", "explanation"],
            headers,
        )

    def test_review_table2(self):
        cta = read_cas_json_file(TEST_JSON2)
        tables = serialize_to_tables(
            cta, "Test_table", OUT_FOLDER, {"accession_id_prefix": "TST_"}
        )

        table_path = os.path.join(OUT_FOLDER, "review.tsv")
        self.assertEqual(table_path, tables[4])
        self.assertTrue(os.path.isfile(table_path))

        headers, records = read_csv_to_dict(
            table_path, generated_ids=True, delimiter="\t"
        )
        self.assertEqual(2, len(records))

        self.assertEqual(
            "CrossArea_cluster:4062c5afea", records[1]["target_node_accession"]
        )
        self.assertEqual("2024-04-01T18:25:43.511Z", records[1]["datestamp"])
        self.assertEqual("Jane Doe", records[1]["reviewer"])
        self.assertEqual("Disagree", records[1]["review"])
        self.assertEqual("This is not a Sst cell.", records[1]["explanation"])

        self.assertEqual(
            "CrossArea_cluster:4062c5afea", records[2]["target_node_accession"]
        )
        self.assertEqual("2024-04-02T20:00:43.511Z", records[2]["datestamp"])
        self.assertEqual("John Doe", records[2]["reviewer"])
        self.assertEqual("Agree", records[2]["review"])
        self.assertEqual(
            "Further expreiments reveal that this a Sst cell.",
            records[2]["explanation"],
        )

    def test_loading_from_json(self):
        cta = read_cas_json_file(TEST_JSON)
        tables = serialize_to_tables(
            cta, "Test_table", OUT_FOLDER, {"accession_id_prefix": "CS202210140_"}
        )

        annotation_table_path = os.path.join(OUT_FOLDER, "annotation.tsv")
        self.assertEqual(annotation_table_path, tables[0])
        self.assertTrue(os.path.isfile(annotation_table_path))

        headers, records = read_csv_to_dict(
            annotation_table_path, id_column_name="cell_set_accession", delimiter="\t"
        )
        self.assertEqual(89, len(records))

        # cluster_1 = records['TST_1']
        # self.assertEqual("1_MSN", cluster_1["cell_label"])
        # self.assertEqual("TST_300", cluster_1["parent_cell_set_accession"])
        # self.assertEqual("D1-Matrix", cluster_1["parent_cell_set_name"])
        # self.assertEqual("cluster", cluster_1["labelset"])
        # self.assertEqual("[\"EPYC\", \"RELN\", \"GULP1\"]", cluster_1["marker_gene_evidence"])
        # self.assertEqual("PuR(0.52) | CaH(0.39)", cluster_1["region.info _Frequency_"])

    def test_loading_from_json2(self):
        cta = read_cas_json_file(TEST_JSON2)
        tables = serialize_to_tables(
            cta, "Test_table", OUT_FOLDER, {"accession_id_prefix": "TST_"}
        )

        annotation_table_path = os.path.join(OUT_FOLDER, "annotation.tsv")
        self.assertEqual(annotation_table_path, tables[0])
        self.assertTrue(os.path.isfile(annotation_table_path))

        headers, records = read_csv_to_dict(
            annotation_table_path, id_column_name="cell_set_accession", delimiter="\t"
        )
        self.assertEqual(160, len(records))

        cluster_1 = records["CrossArea_cluster:4062c5afea"]
        self.assertEqual("Sst_14", cluster_1["cell_label"])
        self.assertEqual(
            "CrossArea_subclass:8fa477a378", cluster_1["parent_cell_set_accession"]
        )
        self.assertEqual("Sst", cluster_1["parent_cell_set_name"])
        self.assertEqual("CrossArea_cluster", cluster_1["labelset"])

        cluster_sst = records["CrossArea_subclass:8fa477a378"]
        self.assertEqual("Sst", cluster_sst["cell_label"])
        self.assertEqual("Class:d4ef18755c", cluster_sst["parent_cell_set_accession"])
        self.assertEqual("inhibitory", cluster_sst["parent_cell_set_name"])
        self.assertEqual("CrossArea_subclass", cluster_sst["labelset"])

    def test_loading_from_json_wmb(self):
        cta = read_cas_json_file(TEST_JSON_WMB)
        tables = serialize_to_tables(cta, "Test_table", OUT_FOLDER, {})

        annotation_table_path = os.path.join(OUT_FOLDER, "annotation.tsv")
        self.assertEqual(annotation_table_path, tables[0])
        self.assertTrue(os.path.isfile(annotation_table_path))

        headers, records = read_csv_to_dict(
            annotation_table_path, id_column_name="cell_set_accession", delimiter="\t"
        )
        self.assertEqual(6905, len(records))

        record_1 = records["CS20230722_CLAS_01"]
        self.assertEqual("01 IT-ET Glut", record_1["cell_label"])
        self.assertEqual("", record_1.get("parent_cell_set_accession"))
        self.assertEqual("", record_1.get("parent_cell_set_name"))
        self.assertEqual("class", record_1["labelset"])
        self.assertEqual("Pallium-Glut", record_1["neighborhood"])

        record_2 = records["CS20230722_SUBC_315"]
        self.assertEqual("315 DCO UBC Glut", record_2["cell_label"])
        self.assertEqual("CS20230722_CLAS_29", record_2["parent_cell_set_accession"])
        self.assertEqual("29 CB Glut", record_2["parent_cell_set_name"])
        self.assertEqual("subclass", record_2["labelset"])
        self.assertEqual("CL:4023161", record_2["cell_ontology_term_id"])
        self.assertEqual("unipolar brush cell", record_2["cell_ontology_term"])
        self.assertEqual("Sln|Lmx1a", record_2["marker_gene_evidence"])
        self.assertEqual("NN-IMN-GC", record_2["neighborhood"])
        self.assertEqual("Eomes,Lmx1a,Klf3", record_2["subclass.tf.markers.combo"])

        record_3 = records["CS20230722_CLUS_4566"]
        self.assertEqual("4566 VCO Mafa Meis2 Glut_4", record_3["cell_label"])
        self.assertEqual("CS20230722_SUPT_1016", record_3["parent_cell_set_accession"])
        self.assertEqual("1016 VCO Mafa Meis2 Glut_4", record_3["parent_cell_set_name"])
        self.assertEqual("cluster", record_3["labelset"])
        self.assertEqual("", record_3["cell_ontology_term_id"])
        self.assertEqual("", record_3["cell_ontology_term"])
        self.assertEqual("Il22|Onecut1", record_3["marker_gene_evidence"])
        self.assertEqual("MB-HB-Glut-Sero-Dopa", record_3["neighborhood"])
        self.assertEqual("DCO", record_3["anatomical_annotation"])
        self.assertEqual("Tnnt1,Ppp1r17,Lhx9", record_3["merfish.markers.combo"])
        self.assertEqual("Onecut1,Mycn,Hoxd3", record_3["cluster.TF.markers.combo"])
        self.assertEqual("Ppfibp1", record_3["cluster.markers.combo _within subclass_"])

    def test_is_hash_accession(self):
        self.assertTrue(is_hash_accession("01d93e7878"))
        self.assertFalse(is_hash_accession(""))
        self.assertFalse(is_hash_accession("AIT_34"))
        self.assertFalse(is_hash_accession("01d9_E7878"))
