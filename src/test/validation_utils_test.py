import unittest
from unittest.mock import patch, MagicMock
import pandas as pd


from cas.utils.validation_utils import validate_labelset_markers, validate_markers


class TestMarkerValidation(unittest.TestCase):
    def setUp(self):
        # Setup a sample dataframe and CAS dictionary
        self.marker_column = 'gene_markers'
        self.data = {
            'gene_markers': ['gene1', 'gene2', 'gene3']
        }
        self.df = pd.DataFrame(self.data)
        self.df_var = pd.DataFrame(self.data)
        self.adata = MagicMock()
        self.adata.var = self.df_var

        self.cas = {
            'annotations': [
                {'labelset': 'set1', 'cell_label': 'A', 'marker_gene_evidence': ['gene1', 'gene2']},
                {'labelset': 'set2', 'cell_label': 'B', 'marker_gene_evidence': ['gene4']}
                # this should trigger a warning
            ]
        }

    def test_validate_markers(self):
        # Test that all markers are found and no exception is raised
        result = validate_markers(self.cas, self.adata, self.marker_column)
        self.assertTrue(result)

        # Test for KeyError when column does not exist
        with self.assertRaises(KeyError):
            validate_markers(self.cas, self.adata, 'nonexistent_column')

    @patch('cas.utils.validation_utils.logger.warning')
    def test_validate_labelset_markers(self, mock_logger_warning):
        # Test with missing marker
        annotation = {'labelset': 'set2', 'cell_label': 'B', 'marker_gene_evidence': ['gene4']}
        validate_labelset_markers(annotation, self.adata.var[self.marker_column].tolist())
        mock_logger_warning.assert_called_once_with(
            "Not all marker genes from set2-B pair exist in anndata's var section. Missing markers: gene4"
        )

        # Test with all markers present
        mock_logger_warning.reset_mock()
        annotation = {'labelset': 'set1', 'cell_label': 'A', 'marker_gene_evidence': ['gene1', 'gene2']}
        validate_labelset_markers(annotation, self.adata.var[self.marker_column].tolist())
        mock_logger_warning.assert_not_called()
