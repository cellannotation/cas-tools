import unittest
from unittest.mock import patch

import pandas as pd
from cap_anndata import read_h5ad

from cas.flatten_data_to_anndata import (
    create_cell_label_lookup,
    unflatten_obs,
    update_cas_annotation,
)
from cas.utils.conversion_utils import (
    ANNOTATIONS,
    AUTHOR_ANNOTATION_FIELDS,
    CELL_IDS,
    CELL_LABEL,
    CELLHASH,
)


class TestAnnotationMethods(unittest.TestCase):
    def setUp(self):
        self.anndata_file_path = (
            "src/test/test_data/cas_cap_roundtrip/test_flatten_dataset.h5ad"
        )
        # This is the synthetic CAS JSON for testing
        self.cas_json = {
            "author_name": "Jordan Lee",
            "title": "Synthetic Brain Cell Atlas v1.0 Taxonomy (Neuronal)",
            "description": "This is a synthetic dataset designed to simulate the organization of neuronal cell types in the human cortex.",
            "matrix_file_id": "Synthetic_dataset:1234abcd-5678-efgh-ijkl-9012mnop3456",
            "labelsets": [
                {
                    "name": "Cluster",
                    "description": "Simulated neuron clusters",
                    "rank": 0,
                },
                {
                    "name": "supercluster_term",
                    "description": "Higher level simulated neuron clusters",
                    "rank": 1,
                },
            ],
            "orcid": "https://orcid.org/0000-0002-1234-5678",
            "cellannotation_schema_version": "0.2b0",
            "annotations": [
                {
                    "author_annotation_fields": {
                        "cellhash": "Cluster:6e98fec3ec",
                        "Cluster ID": "50",
                    },
                    "labelset": "Cluster",
                    "cell_label": "O50",
                    "cell_set_accession": "CS202210140_51",
                    "parent_cell_set_accession": "CS202210140_469",
                    "cell_ontology_term_id": "CL:0000128",
                    "cell_ontology_term": "oligodendrocyte",
                    "rationale": "Supported by marker expression and annotation transfer from Middle Temporal Gyrus dataset (Jorstad et al., 2023)",
                    "marker_gene_evidence": "PLP1, SOX10",
                    "rationale_dois": "DOI:10.1126/science.adf6812",
                    "cell_ids": [
                        "10X362_3:TCAGTGAGTATTGACC",
                        "10X362_5:TCCGTGTGTGAAAGTT",
                    ],
                },
                {
                    "author_annotation_fields": {
                        "cellhash": "Cluster:b81a00daa1",
                        "Cluster ID": "49",
                    },
                    "labelset": "Cluster",
                    "cell_label": "O49",
                    "cell_set_accession": "CS202210140_50",
                    "parent_cell_set_accession": "CS202210140_469",
                    "cell_ontology_term_id": "CL:0000128",
                    "cell_ontology_term": "oligodendrocyte",
                    "rationale": "Supported by marker expression and annotation transfer from Middle Temporal Gyrus dataset (Jorstad et al., 2023)",
                    "marker_gene_evidence": "PLP1, SOX10",
                    "rationale_dois": "DOI:10.1126/science.adf6812",
                    "cell_ids": [
                        "10X379_2:ATTCCATTCCCAGCGA",
                        "10X383_4:GAAGTAAAGGTTCTTG",
                    ],
                },
                {
                    "author_annotation_fields": {
                        "cellhash": "Cluster:0a1cfc9729",
                        "Cluster ID": "62",
                    },
                    "labelset": "Cluster",
                    "cell_label": "A62",
                    "cell_set_accession": "CS202210140_63",
                    "parent_cell_set_accession": "CS202210140_470",
                    "cell_ontology_term_id": "CL:0000127",
                    "cell_ontology_term": "astrocyte",
                    "rationale": "Supported by marker expression and annotation transfer from Middle Temporal Gyrus dataset (Jorstad et al., 2023)",
                    "marker_gene_evidence": "AQP4",
                    "rationale_dois": "DOI:10.1126/science.adf6812",
                    "cell_ids": [
                        "10X357_2:TGGGCTGAGAAACCCG",
                        "10X319_7:TGCTCCATCATCACCC",
                    ],
                },
            ],
        }

    def test_unflatten_obs(self):
        # Expected output should be a processed JSON structure
        expected_output = {
            "matrix_file_id": "CellXGene_dataset:b165f033-9dec-468a-9248-802fc6902a74",
            "title": "Human Brain Cell Atlas v1.0 Taxonomy (Non-Neuronal)",
            "description": "Atlas of human primary motor cortex (M1), developed in collaboration with the BRAIN Initiative Cell Census Network (BICCN).\nFirst draft atlas of human brain transcriptomic cell types: The human brain directs a wide range of complex behaviors ranging from fine motor skills to abstract intelligence and emotion. This broad range of functions is supported by an exceptionally complex cellular and circuit architecture. To create a first draft human brain cell atlas, high-throughput single-nucleus RNA sequencing was used to systematically survey cells across the entire adult human brain in three postmortem donors. Over three million nuclei were sampled from approximately 100 dissections across the forebrain, midbrain, and hindbrain. Analysis of these data showed regional diversity in that cellular organization exhibited regional diversity at multiple scales, identifying 30 superclusters, 461 clusters and 3313 subclusters. As the first single-cell transcriptomic census of the entire human brain, this atlas provides a resource for understanding the molecular diversity of the human brain in health and disease. The Human Brain Cell Atlas v1.0 is presented for visualization and data mining through the Chan Zuckerberg Initiative’s CellxGene application, with the following biologically meaningful partitions: 1. Neuronal and non-neuronal cell types 2. Supercluster-specific groupings (`Supercluster: `) 3.Brain region-specific groupings (`Dissection: `), ordered by the adult human brain anatomical reference atlas ontology in Ding et al. (2016)",
            "cellannotation_schema_version": "0.2b0",
            "author_name": "Kimberly Siletti",
            "orcid": "https://orcid.org/0000-0001-7620-8973",
            "annotations": [
                {
                    "author_annotation_fields": {"Cluster ID": "O50"},
                    "labelset": "Cluster",
                    "cell_label": "O50",
                    "cell_fullname": "O50",
                    "cell_set_accession": "Cluster:5a972d4730",
                    "parent_cell_set_name": "Oligodendrocyte",
                    "parent_cell_set_accession": "supercluster_term:502566eede",
                },
                {
                    "author_annotation_fields": {"Cluster ID": "O50"},
                    "labelset": "Cluster",
                    "cell_label": "O50",
                    "cell_fullname": "O50",
                    "cell_set_accession": "Cluster:5a972d4730",
                    "parent_cell_set_name": "Oligodendrocyte",
                    "parent_cell_set_accession": "supercluster_term:502566eede",
                },
                {
                    "author_annotation_fields": {"Cluster ID": "O40"},
                    "labelset": "Cluster",
                    "cell_label": "O40",
                    "cell_fullname": "O40",
                    "cell_set_accession": "Cluster:7c94e6181d",
                    "parent_cell_set_name": "Oligodendrocyte",
                    "parent_cell_set_accession": "supercluster_term:502566eede",
                },
                {
                    "author_annotation_fields": {"Cluster ID": "O40"},
                    "labelset": "Cluster",
                    "cell_label": "O40",
                    "cell_fullname": "O40",
                    "cell_set_accession": "Cluster:7c94e6181d",
                    "parent_cell_set_name": "Oligodendrocyte",
                    "parent_cell_set_accession": "supercluster_term:502566eede",
                },
                {
                    "author_annotation_fields": {"Cluster ID": "A62"},
                    "labelset": "Cluster",
                    "cell_label": "A62",
                    "cell_fullname": "A62",
                    "cell_set_accession": "Cluster:c2b38b36d7",
                    "parent_cell_set_name": "Astrocyte",
                    "parent_cell_set_accession": "supercluster_term:c2b38b36d7",
                },
                {
                    "author_annotation_fields": {"Cluster ID": "A62"},
                    "labelset": "Cluster",
                    "cell_label": "A62",
                    "cell_fullname": "A62",
                    "cell_set_accession": "Cluster:c2b38b36d7",
                    "parent_cell_set_name": "Astrocyte",
                    "parent_cell_set_accession": "supercluster_term:c2b38b36d7",
                },
                {
                    "author_annotation_fields": {"Cluster ID": "None"},
                    "labelset": "supercluster_term",
                    "cell_label": "Oligodendrocyte",
                    "cell_fullname": "Oligodendrocyte",
                    "cell_set_accession": "supercluster_term:502566eede",
                    "cell_ontology_term_id": "CL:0000128",
                    "cell_ontology_term": "oligodendrocyte",
                    "rationale": "Supported by marker expression and annotation transfer from Middle Temporal Gyrus dataset (Jorstad et al., 2023)",
                    "rationale_dois": "DOI:10.1126/science.adf6812",
                    "marker_gene_evidence": "PLP1, SOX10",
                },
                {
                    "author_annotation_fields": {"Cluster ID": "None"},
                    "labelset": "supercluster_term",
                    "cell_label": "Oligodendrocyte",
                    "cell_fullname": "Oligodendrocyte",
                    "cell_set_accession": "supercluster_term:502566eede",
                    "cell_ontology_term_id": "CL:0000128",
                    "cell_ontology_term": "oligodendrocyte",
                    "rationale": "Supported by marker expression and annotation transfer from Middle Temporal Gyrus dataset (Jorstad et al., 2023)",
                    "rationale_dois": "DOI:10.1126/science.adf6812",
                    "marker_gene_evidence": "PLP1, SOX10",
                },
                {
                    "author_annotation_fields": {"Cluster ID": "None"},
                    "labelset": "supercluster_term",
                    "cell_label": "Astrocyte",
                    "cell_fullname": "Astrocyte",
                    "cell_set_accession": "supercluster_term:c2b38b36d7",
                    "cell_ontology_term_id": "CL:0000127",
                    "cell_ontology_term": "astrocyte",
                    "rationale": "Supported by marker expression and annotation transfer from Middle Temporal Gyrus dataset (Jorstad et al., 2023)",
                    "rationale_dois": "DOI:10.1126/science.adf6812",
                    "marker_gene_evidence": "AQP4",
                },
                {
                    "author_annotation_fields": {"Cluster ID": "None"},
                    "labelset": "supercluster_term",
                    "cell_label": "Astrocyte",
                    "cell_fullname": "Astrocyte",
                    "cell_set_accession": "supercluster_term:c2b38b36d7",
                    "cell_ontology_term_id": "CL:0000127",
                    "cell_ontology_term": "astrocyte",
                    "rationale": "Supported by marker expression and annotation transfer from Middle Temporal Gyrus dataset (Jorstad et al., 2023)",
                    "rationale_dois": "DOI:10.1126/science.adf6812",
                    "marker_gene_evidence": "AQP4",
                },
            ],
            "labelsets": [
                {"name": "Cluster", "description": "", "rank": "0"},
                {"name": "supercluster_term", "description": "", "rank": "1"},
            ],
        }
        with patch(
            "cas.flatten_data_to_anndata.update_cas_annotation",
            return_value=expected_output,
        ) as mock_update_cas_annotation:
            with read_h5ad(file_path=self.anndata_file_path, edit=True) as cap_adata:
                cap_adata.read_obs()
                cap_adata.read_uns()
                obs = cap_adata.obs
                uns = cap_adata.uns
                cellhash_lookup = cap_adata.uns[CELLHASH]
                obs_columns_count = len(obs.columns)
                new_cas = unflatten_obs(obs, uns, self.cas_json, cellhash_lookup)
                col_count_diff = obs_columns_count - len(obs.columns)

            # # Assertions
            # mock_update_cas_annotation.assert_called_once_with(
            #     {"updated_annotations": "test"}, self.cas_json
            # )
            self.assertEqual(new_cas, expected_output)
            self.assertEqual(col_count_diff, 13)

    def test_create_cell_label_lookup(self):
        # Create mock DataFrames dictionary
        cluster_df = pd.DataFrame(
            {
                "Cluster": ["O50", "O50", "O49", "O49", "A62", "A62"],
                "Cluster--cell_set_accession": [
                    "CS202210140_51",
                    "CS202210140_51",
                    "CS202210140_50",
                    "CS202210140_50",
                    "CS202210140_63",
                    "CS202210140_63",
                ],
                "Cluster--parent_cell_set_accession": ["CS202210140_469"] * 4
                + ["CS202210140_470"] * 2,
                "Cluster--author_annotation_fields": [
                    "{'Cluster ID': '50', 'cellhash': 'Cluster:6e98fec3ec'}",
                    "{'Cluster ID': '50', 'cellhash': 'Cluster:6e98fec3ec'}",
                    "{'Cluster ID': '49', 'cellhash': 'Cluster:b81a00daa1'}",
                    "{'Cluster ID': '49', 'cellhash': 'Cluster:b81a00daa1'}",
                    "{'Cluster ID': '62', 'cellhash': 'Cluster:0a1cfc9729'}",
                    "{'Cluster ID': '62', 'cellhash': 'Cluster:0a1cfc9729'}",
                ],
                "Cluster--cell_ontology_term_id": [
                    "CL:0000128",
                    "CL:0000128",
                    "CL:0000128",
                    "CL:0000128",
                    "CL:0000127",
                    "CL:0000127",
                ],
                "Cluster--cell_ontology_term": [
                    "oligodendrocyte",
                    "oligodendrocyte",
                    "oligodendrocyte",
                    "oligodendrocyte",
                    "astrocyte",
                    "astrocyte",
                ],
                "Cluster--rationale": [
                    "Supported by marker expression and annotation transfer from Middle Temporal Gyrus dataset (Jorstad et al., 2023)"
                ]
                * 6,
                "Cluster--marker_gene_evidence": ["PLP1, SOX10"] * 4 + ["AQP4"] * 2,
                "Cluster--rationale_dois": ["DOI:10.1126/science.adf6812"] * 6,
            },
            index=[
                "10X362_3:TCAGTGAGTATTGACC",
                "10X362_5:TCCGTGTGTGAAAGTT",
                "10X379_2:ATTCCATTCCCAGCGA",
                "10X383_4:GAAGTAAAGGTTCTTG",
                "10X357_2:TGGGCTGAGAAACCCG",
                "10X319_7:TGCTCCATCATCACCC",
            ],
        )

        supercluster_term_df = pd.DataFrame(
            {
                "supercluster_term": ["supercluster_term:21eaacf654"] * 4
                + ["supercluster_term:1dc795d1ea"] * 2,
                "supercluster_term--cell_set_accession": [
                    "SC202210140_01",
                    "SC202210140_01",
                    "SC202210140_02",
                    "SC202210140_02",
                    "SC202210140_03",
                    "SC202210140_03",
                ],
                "supercluster_term--cell_ontology_term_id": [
                    "CL:0000540",
                    "CL:0000540",
                    "CL:0000541",
                    "CL:0000541",
                    "CL:0000542",
                    "CL:0000542",
                ],
                "supercluster_term--cell_ontology_term": [
                    "neuron",
                    "neuron",
                    "microglia",
                    "microglia",
                    "astrocyte",
                    "astrocyte",
                ],
                "supercluster_term--rationale": ["Randomly assigned rationale"] * 6,
                "supercluster_term--rationale_dois": ["DOI:10.1016/j.cell.2023.01.001"]
                * 6,
                "supercluster_term--marker_gene_evidence": ["GENE1, GENE2"] * 6,
                "supercluster_term--author_annotation_fields": [
                    "{'Cluster ID': 'None', 'cellhash': 'supercluster_term:21eaacf654'}"
                ]
                * 4
                + ["{'Cluster ID': 'None', 'cellhash': 'supercluster_term:1dc795d1ea'}"]
                * 2,
            },
            index=[
                "10X362_3:TCAGTGAGTATTGACC",
                "10X362_5:TCCGTGTGTGAAAGTT",
                "10X379_2:ATTCCATTCCCAGCGA",
                "10X383_4:GAAGTAAAGGTTCTTG",
                "10X357_2:TGGGCTGAGAAACCCG",
                "10X319_7:TGCTCCATCATCACCC",
            ],
        )

        df_dict = {"Cluster": cluster_df, "supercluster_term": supercluster_term_df}

        # Run the create_cell_label_lookup function
        result = create_cell_label_lookup(df_dict)

        # Assert the expected structure and contents
        self.assertIsInstance(result, dict)
        self.assertEqual(len(result.keys()), 10)
        self.assertEqual(len(result["Cluster:O49"]), 11)
        self.assertIn("Cluster:O50", result)
        self.assertIn("Cluster:O49", result)
        self.assertIn("Cluster:A62", result)
        self.assertIn(CELL_IDS, result["Cluster:O50"])
        self.assertIn(AUTHOR_ANNOTATION_FIELDS, result["Cluster:O50"])
        self.assertIsInstance(result["Cluster:O50"][CELL_IDS], list)

    def test_update_cas_annotation(self):
        # Mock inputs based on synthetic data
        cas_dict = {
            "Cluster:A62": {
                "cell_ids": ["10X357_2:TGGGCTGAGAAACCCG", "10X319_7:TGCTCCATCATCACCC"],
                "author_annotation_fields": {
                    "cellhash": "Cluster:0a1cfc9729",
                    "Cluster ID": "62",
                },
                "labelset": "Cluster",
                "cell_label": "A62",
                "cell_set_accession": "CS202210140_63",
                "parent_cell_set_accession": "CS202210140_470",
                "cell_ontology_term_id": "CL:0000127",
                "cell_ontology_term": "astrocyte",
                "rationale": "Supported by marker expression and annotation transfer from Middle Temporal Gyrus dataset (Jorstad et al., 2023)",
                "marker_gene_evidence": "AQP4",
                "rationale_dois": "DOI:10.1126/science.adf6812",
            },
            "Cluster:c2b38b36d7": {
                "cell_ids": ["10X357_2:TGGGCTGAGAAACCCG", "10X319_7:TGCTCCATCATCACCC"],
                "author_annotation_fields": {
                    "cellhash": "Cluster:0a1cfc9729",
                    "Cluster ID": "62",
                },
                "labelset": "Cluster",
                "cell_label": "A62",
                "cell_set_accession": "CS202210140_63",
                "parent_cell_set_accession": "CS202210140_470",
                "cell_ontology_term_id": "CL:0000127",
                "cell_ontology_term": "astrocyte",
                "rationale": "Supported by marker expression and annotation transfer from Middle Temporal Gyrus dataset (Jorstad et al., 2023)",
                "marker_gene_evidence": "AQP4",
                "rationale_dois": "DOI:10.1126/science.adf6812",
            },
            "Cluster:O49": {
                "cell_ids": ["10X379_2:ATTCCATTCCCAGCGA", "10X383_4:GAAGTAAAGGTTCTTG"],
                "author_annotation_fields": {
                    "cellhash": "Cluster:b81a00daa1",
                    "Cluster ID": "49",
                },
                "labelset": "Cluster",
                "cell_label": "O49",
                "cell_set_accession": "CS202210140_50",
                "parent_cell_set_accession": "CS202210140_469",
                "cell_ontology_term_id": "CL:0000128",
                "cell_ontology_term": "oligodendrocyte",
                "rationale": "Supported by marker expression and annotation transfer from Middle Temporal Gyrus dataset (Jorstad et al., 2023)",
                "marker_gene_evidence": "PLP1, SOX10",
                "rationale_dois": "DOI:10.1126/science.adf6812",
            },
            "Cluster:7c94e6181d": {
                "cell_ids": ["10X379_2:ATTCCATTCCCAGCGA", "10X383_4:GAAGTAAAGGTTCTTG"],
                "author_annotation_fields": {
                    "cellhash": "Cluster:b81a00daa1",
                    "Cluster ID": "49",
                },
                "labelset": "Cluster",
                "cell_label": "O49",
                "cell_set_accession": "CS202210140_50",
                "parent_cell_set_accession": "CS202210140_469",
                "cell_ontology_term_id": "CL:0000128",
                "cell_ontology_term": "oligodendrocyte",
                "rationale": "Supported by marker expression and annotation transfer from Middle Temporal Gyrus dataset (Jorstad et al., 2023)",
                "marker_gene_evidence": "PLP1, SOX10",
                "rationale_dois": "DOI:10.1126/science.adf6812",
            },
            "Cluster:O500x": {
                "cell_ids": ["10X362_3:TCAGTGAGTATTGACC", "10X362_5:TCCGTGTGTGAAAGTT"],
                "author_annotation_fields": {
                    "cellhash": "Cluster:6e98fec3ec",
                    "Cluster ID": "50",
                },
                "labelset": "Cluster",
                "cell_label": "O500x",
                "cell_set_accession": "CS202210140_51",
                "parent_cell_set_accession": "CS202210140_469",
                "cell_ontology_term_id": "CL:0000128",
                "cell_ontology_term": "oligodendrocyte",
                "rationale": "Supported by marker expression and annotation transfer from Middle Temporal Gyrus dataset (Jorstad et al., 2023)",
                "marker_gene_evidence": "PLP1, SOX10",
                "rationale_dois": "DOI:10.1126/science.adf6812",
            },
            "Cluster:6e98fec3ec": {
                "cell_ids": ["10X362_3:TCAGTGAGTATTGACC", "10X362_5:TCCGTGTGTGAAAGTT"],
                "author_annotation_fields": {
                    "cellhash": "Cluster:6e98fec3ec",
                    "Cluster ID": "50",
                },
                "labelset": "Cluster",
                "cell_label": "O500x",
                "cell_set_accession": "CS202210140_51",
                "parent_cell_set_accession": "CS202210140_469",
                "cell_ontology_term_id": "CL:0000128",
                "cell_ontology_term": "oligodendrocyte",
                "rationale": "Supported by marker expression and annotation transfer from Middle Temporal Gyrus dataset (Jorstad et al., 2023)",
                "marker_gene_evidence": "PLP1, SOX10",
                "rationale_dois": "DOI:10.1126/science.adf6812",
            },
            "supercluster_term:1dc795d1ea": {
                "cell_ids": ["10X357_2:TGGGCTGAGAAACCCG", "10X319_7:TGCTCCATCATCACCC"],
                "author_annotation_fields": {
                    "cellhash": "supercluster_term:c2b38b36d7",
                    "Cluster ID": "None",
                },
                "labelset": "supercluster_term",
                "cell_label": "supercluster_term:1dc795d1ea",
                "cell_set_accession": "SC202210140_03",
                "cell_ontology_term_id": "CL:0000542",
                "cell_ontology_term": "astrocyte",
                "rationale": "Randomly assigned rationale",
                "rationale_dois": "DOI:10.1016/j.cell.2023.01.001",
                "marker_gene_evidence": "GENE1, GENE2",
            },
            "supercluster_term:c2b38b36d7": {
                "cell_ids": ["10X357_2:TGGGCTGAGAAACCCG", "10X319_7:TGCTCCATCATCACCC"],
                "author_annotation_fields": {
                    "cellhash": "supercluster_term:c2b38b36d7",
                    "Cluster ID": "None",
                },
                "labelset": "supercluster_term",
                "cell_label": "supercluster_term:1dc795d1ea",
                "cell_set_accession": "SC202210140_03",
                "cell_ontology_term_id": "CL:0000542",
                "cell_ontology_term": "astrocyte",
                "rationale": "Randomly assigned rationale",
                "rationale_dois": "DOI:10.1016/j.cell.2023.01.001",
                "marker_gene_evidence": "GENE1, GENE2",
            },
            "supercluster_term:21eaacf654": {
                "cell_ids": [
                    "10X362_3:TCAGTGAGTATTGACC",
                    "10X362_5:TCCGTGTGTGAAAGTT",
                    "10X379_2:ATTCCATTCCCAGCGA",
                    "10X383_4:GAAGTAAAGGTTCTTG",
                ],
                "author_annotation_fields": {
                    "cellhash": "supercluster_term:502566eede",
                    "Cluster ID": "None",
                },
                "labelset": "supercluster_term",
                "cell_label": "supercluster_term:21eaacf654",
                "cell_set_accession": "SC202210140_01",
                "cell_ontology_term_id": "CL:0000540",
                "cell_ontology_term": "neuron",
                "rationale": "Randomly assigned rationale",
                "rationale_dois": "DOI:10.1016/j.cell.2023.01.001",
                "marker_gene_evidence": "GENE1, GENE2",
            },
            "supercluster_term:502566eede": {
                "cell_ids": [
                    "10X362_3:TCAGTGAGTATTGACC",
                    "10X362_5:TCCGTGTGTGAAAGTT",
                    "10X379_2:ATTCCATTCCCAGCGA",
                    "10X383_4:GAAGTAAAGGTTCTTG",
                ],
                "author_annotation_fields": {
                    "cellhash": "supercluster_term:502566eede",
                    "Cluster ID": "None",
                },
                "labelset": "supercluster_term",
                "cell_label": "supercluster_term:21eaacf654",
                "cell_set_accession": "SC202210140_01",
                "cell_ontology_term_id": "CL:0000540",
                "cell_ontology_term": "neuron",
                "rationale": "Randomly assigned rationale",
                "rationale_dois": "DOI:10.1016/j.cell.2023.01.001",
                "marker_gene_evidence": "GENE1, GENE2",
            },
        }

        cellhash = {
            "Cluster:A62": "Cluster:c2b38b36d7",
            "Cluster:O40": "Cluster:7c94e6181d",
            "Cluster:O50": "Cluster:5a972d4730",
            "supercluster_term:Astrocyte": "supercluster_term:c2b38b36d7",
            "supercluster_term:Oligodendrocyte": "supercluster_term:502566eede",
        }
        # Running the update_cas_annotation function
        result = update_cas_annotation(cas_dict, self.cas_json, cellhash)

        # Check the updated annotations
        self.assertEqual(len(result), 4)
        self.assertIn("O500x", [a[CELL_LABEL] for a in result])
        self.assertIn("A62", [a[CELL_LABEL] for a in result])
        self.assertIn(
            "supercluster_term:1dc795d1ea", [a[CELL_LABEL] for a in result]
        )

        self.assertIn(
            "supercluster_term:21eaacf654", [a[CELL_LABEL] for a in result]
        )


# This block is necessary to run the test cases
if __name__ == "__main__":
    unittest.main()
