import json
import os
import unittest
import warnings
from unittest.mock import patch

import anndata as ad
import pandas as pd
from cas.spreadsheet_to_cas import read_spreadsheet, retrieve_schema, spreadsheet2cas

warnings.filterwarnings("ignore", category=UserWarning, module="anndata._core.anndata")

TEST_SPREADSHEET = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "./test_data/sample_spreadsheet_data.xlsx",
)
MINI_TEST_SPREADSHEET = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "./test_data/mini_spreadsheet_data.xlsx",
)


def generate_mock_dataset():
    # Mock AnnData dataset for testing
    cell_ids = ["cell_1", "cell_2", "cell_3", "cell_4", "cell_5"]
    return ad.AnnData(
        obs=pd.DataFrame(
            {
                "Ethnicity": ["EUR", "EUR", "EUR", "EUR", "EUR"],
                "BMI": ["Unknown", "Unknown", "Unknown", "Unknown", "Unknown"],
                "annotation_broad": [
                    "Monocyte",
                    "T CD4+",
                    "T CD4+",
                    "T CD8+",
                    "T CD4+",
                ],
                "annotation_detailed": [
                    "Monocyte CD14",
                    "T CD4 helper",
                    "T CD4 helper",
                    "T CD8 naive",
                    "T CD4 naive",
                ],
                "annotation_detailed_fullNames": [
                    "Classical monocyte",
                    "T CD4 helper",
                    "T CD4 helper",
                    "T CD8 naive",
                    "T CD4 naive",
                ],
                "Age_group": ["Adult", "Adult", "Adult", "Adult", "Adult"],
                "COVID_severity": [
                    "Healthy",
                    "Healthy",
                    "Healthy",
                    "Healthy",
                    "Healthy",
                ],
                "COVID_status": [
                    "Healthy",
                    "Healthy",
                    "Healthy",
                    "Healthy",
                    "Healthy",
                ],
                "Group": ["Adult", "Adult", "Adult", "Adult", "Adult"],
                "Smoker": [
                    "Non-smoker",
                    "Non-smoker",
                    "Non-smoker",
                    "Non-smoker",
                    "Non-smoker",
                ],
                "sample_id": ["AN5", "AN5", "AN3", "AN5", "AN5"],
                "sequencing_library": [
                    "CV001_KM10202384-CV001_KM10202394",
                    "CV001_KM10202384-CV001_KM10202394",
                    "CV001_KM10202384-CV001_KM10202394",
                    "CV001_KM10202384-CV001_KM10202394",
                    "CV001_KM10202384-CV001_KM10202394",
                ],
                "Protein_modality_weight": [
                    0.3595166802406311,
                    0.5775224566459656,
                    0.3691430389881134,
                    0.785563051700592,
                    0.5641735792160034,
                ],
                "assay_ontology_term_id": [
                    "EFO:0011025",
                    "EFO:0011025",
                    "EFO:0011025",
                    "EFO:0011025",
                    "EFO:0011025",
                ],
                "cell_type_ontology_term_id": [
                    "CL:0000860",
                    "CL:0000492",
                    "CL:0000492",
                    "CL:0000900",
                    "CL:0000895",
                ],
                "development_stage_ontology_term_id": [
                    "HsapDv:0000087",
                    "HsapDv:0000087",
                    "HsapDv:0000087",
                    "HsapDv:0000087",
                    "HsapDv:0000087",
                ],
                "disease_ontology_term_id": [
                    "PATO:0000461",
                    "PATO:0000461",
                    "PATO:0000461",
                    "PATO:0000461",
                    "PATO:0000461",
                ],
                "is_primary_data": [True, True, True, True, True],
                "organism_ontology_term_id": [
                    "NCBITaxon:9606",
                    "NCBITaxon:9606",
                    "NCBITaxon:9606",
                    "NCBITaxon:9606",
                    "NCBITaxon:9606",
                ],
                "sex_ontology_term_id": [
                    "PATO:0000383",
                    "PATO:0000383",
                    "PATO:0000384",
                    "PATO:0000383",
                    "PATO:0000383",
                ],
                "tissue_ontology_term_id": [
                    "UBERON:0000178",
                    "UBERON:0000178",
                    "UBERON:0000178",
                    "UBERON:0000178",
                    "UBERON:0000178",
                ],
                "self_reported_ethnicity_ontology_term_id": [
                    "HANCESTRO:0005",
                    "HANCESTRO:0005",
                    "HANCESTRO:0005",
                    "HANCESTRO:0005",
                    "HANCESTRO:0005",
                ],
                "donor_id": ["AN5", "AN5", "AN3", "AN5", "AN5"],
                "suspension_type": ["cell", "cell", "cell", "cell", "cell"],
                "cell_type": [
                    "classical monocyte",
                    "CD4-positive helper T cell",
                    "CD4-positive helper T cell",
                    "naive thymus-derived CD8-positive, alpha-beta T cell",
                    "naive thymus-derived CD4-positive, alpha-beta T cell",
                ],
                "assay": [
                    "10x 5' v1",
                    "10x 5' v1",
                    "10x 5' v1",
                    "10x 5' v1",
                    "10x 5' v1",
                ],
                "disease": ["normal", "normal", "normal", "normal", "normal"],
                "organism": [
                    "Homo sapiens",
                    "Homo sapiens",
                    "Homo sapiens",
                    "Homo sapiens",
                    "Homo sapiens",
                ],
                "sex": ["female", "female", "male", "female", "female"],
                "tissue": ["blood", "blood", "blood", "blood", "blood"],
                "self_reported_ethnicity": [
                    "European",
                    "European",
                    "European",
                    "European",
                    "European",
                ],
                "development_stage": [
                    "human adult stage",
                    "human adult stage",
                    "human adult stage",
                    "human adult stage",
                    "human adult stage",
                ],
            },
            index=cell_ids,
        )
    )


class SpreadsheetToCasTests(unittest.TestCase):
    def test_read_spreadsheet_default_sheet(self):
        # Test reading spreadsheet with default sheet
        meta_data, column_names, raw_data = read_spreadsheet(
            TEST_SPREADSHEET, None, retrieve_schema("cap")
        )
        self.assertEqual(len(meta_data), 5)
        self.assertEqual(len(column_names), 9)
        self.assertEqual(raw_data.shape, (73, 9))

    def test_read_spreadsheet_custom_sheet(self):
        # Test reading spreadsheet with custom sheet
        meta_data, column_names, raw_data = read_spreadsheet(
            TEST_SPREADSHEET,
            sheet_name="PBMC3_Yoshida_2022_PBMC",
            schema=retrieve_schema("cap"),
        )
        self.assertEqual(len(meta_data), 5)
        self.assertEqual(len(column_names), 9)
        self.assertEqual(raw_data.shape, (73, 9))

    @patch("cellxgene_census.download_source_h5ad", return_value=None)
    @patch(
        "cas.spreadsheet_to_cas.read_anndata_file", return_value=generate_mock_dataset()
    )
    def test_spreadsheet2cas(self, mock_read_anndata_file, mock_download_source_h5ad):
        spreadsheet2cas(MINI_TEST_SPREADSHEET, None, None, None, None, "output.json")

        json_file_path = "output.json"

        try:
            with open(json_file_path, "r") as json_file:
                json_data = json.load(json_file)

            self.assertEqual(len(json_data), 8)
            self.assertEqual(len(json_data["annotations"]), 8)
            self.assertEqual(len(json_data["annotations"][0]), 7)
            self.assertEqual(len(json_data["labelsets"]), 2)
        finally:
            # Remove the JSON file after the test
            if os.path.exists(json_file_path):
                os.remove(json_file_path)
