import unittest

from cas.dataset_retrieval.dataset_retriever import DatasetRetriever
from cas.dataset_retrieval.cxg_downloader import CxGDownloader
from cas.dataset_retrieval.http_downloader import HTTPDownloader
from cas.dataset_retrieval.s3_downloader import S3Downloader


class TestDatasetRetriever(unittest.TestCase):
    def test_create_cxg_downloader(self):
        matrix_field_id = "cxg_dataset:b165f033-9dec-468a-9248-802fc6902a74"
        raw_matrix_field_id = "b165f033-9dec-468a-9248-802fc6902a74"
        dataset_retriever = DatasetRetriever.create(matrix_field_id)
        self.assertIsInstance(dataset_retriever, CxGDownloader)
        self.assertEqual(dataset_retriever.matrix_id, raw_matrix_field_id)

    def test_create_http_downloader(self):
        matrix_field_id = "https://datasets.cellxgene.cziscience.com/d911082e-b5e2-40e6-8f08-fb53c7894622.h5ad"
        dataset_retriever = DatasetRetriever.create(matrix_field_id)
        self.assertIsInstance(dataset_retriever, HTTPDownloader)
        self.assertEqual(dataset_retriever.matrix_id, matrix_field_id)

    def test_create_s3_downloader(self):
        matrix_field_id = "s3:b165f033-9dec-468a"
        dataset_retriever = DatasetRetriever.create(matrix_field_id)
        self.assertIsInstance(dataset_retriever, S3Downloader)
        self.assertEqual(dataset_retriever.matrix_id, matrix_field_id)


if __name__ == '__main__':
    unittest.main()
