import os
import unittest
from unittest import mock

from cas.dataset_retrieval.cxg_downloader import CxGDownloader
from cas.dataset_retrieval.dataset_retriever import (
    DatasetRetriever,
    check_file_exists,
    construct_full_download_path,
    create_directory_if_missing,
)
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

    @mock.patch("os.path.exists")
    @mock.patch("os.path.abspath")
    @mock.patch("logging.info")
    def test_check_file_exists_file_exists(
        self, mock_logging, mock_abspath, mock_exists
    ):
        # Mock the file path and existence check
        mock_abspath.return_value = "/mock/path/to/file.h5ad"
        mock_exists.return_value = True

        file_name = "file.h5ad"
        result = check_file_exists(file_name)

        # Verify the correct path and logging
        mock_abspath.assert_called_once_with(file_name)
        mock_logging.assert_called_once_with(
            "File '/mock/path/to/file.h5ad' already exists. Skipping download."
        )
        self.assertTrue(result)

    @mock.patch("os.path.exists")
    @mock.patch("os.path.abspath")
    def test_check_file_exists_file_does_not_exist(self, mock_abspath, mock_exists):
        # Mock the file path and existence check
        mock_abspath.return_value = "/mock/path/to/file.h5ad"
        mock_exists.return_value = False

        file_name = "file.h5ad"
        result = check_file_exists(file_name)

        # Verify the correct path and no logging
        mock_abspath.assert_called_once_with(file_name)
        mock_exists.assert_called_once_with("/mock/path/to/file.h5ad")
        self.assertFalse(result)

    @mock.patch("os.makedirs")
    @mock.patch("os.path.exists")
    @mock.patch("os.path.dirname")
    def test_create_directory_if_missing_directory_does_not_exist(
        self, mock_dirname, mock_exists, mock_makedirs
    ):
        # Mock the directory path and existence check
        mock_dirname.return_value = "/mock/path/to"
        mock_exists.return_value = False

        file_name = "/mock/path/to/file.h5ad"
        create_directory_if_missing(file_name)

        # Verify that the directory is created
        mock_dirname.assert_called_once_with(file_name)
        mock_exists.assert_called_once_with("/mock/path/to")
        mock_makedirs.assert_called_once_with("/mock/path/to", exist_ok=True)

    @mock.patch("os.makedirs")
    @mock.patch("os.path.exists")
    @mock.patch("os.path.dirname")
    def test_create_directory_if_missing_directory_exists(
        self, mock_dirname, mock_exists, mock_makedirs
    ):
        # Mock the directory path and existence check
        mock_dirname.return_value = "/mock/path/to"
        mock_exists.return_value = True

        file_name = "/mock/path/to/file.h5ad"
        create_directory_if_missing(file_name)

        # Verify that no directory creation is attempted
        mock_dirname.assert_called_once_with(file_name)
        mock_exists.assert_called_once_with("/mock/path/to")
        mock_makedirs.assert_not_called()

    def test_construct_full_download_path_with_file_name_and_dir(self):
        # Test when both file_name and download_dir are provided
        file_name = "test_file.h5ad"
        download_dir = "/mock/download/directory"
        default_file_name = "default_file.h5ad"

        expected_path = os.path.join(download_dir, file_name)
        result = construct_full_download_path(
            file_name, download_dir, default_file_name
        )

        self.assertEqual(result, expected_path)

    def test_construct_full_download_path_with_default_file_name(self):
        # Test when file_name is None (default_file_name should be used)
        file_name = None
        download_dir = "/mock/download/directory"
        default_file_name = "default_file.h5ad"

        expected_path = os.path.join(download_dir, default_file_name)
        result = construct_full_download_path(
            file_name, download_dir, default_file_name
        )

        self.assertEqual(result, expected_path)

    def test_construct_full_download_path_with_no_download_dir(self):
        # Test when download_dir is None (should default to current directory)
        file_name = "test_file.h5ad"
        download_dir = None
        default_file_name = "default_file.h5ad"

        expected_path = os.path.join("", file_name)
        result = construct_full_download_path(
            file_name, download_dir, default_file_name
        )

        self.assertEqual(result, expected_path)

    def test_construct_full_download_path_with_no_file_name_and_no_download_dir(self):
        # Test when both file_name and download_dir are None
        file_name = None
        download_dir = None
        default_file_name = "default_file.h5ad"

        expected_path = os.path.join("", default_file_name)
        result = construct_full_download_path(
            file_name, download_dir, default_file_name
        )

        self.assertEqual(result, expected_path)


if __name__ == "__main__":
    unittest.main()
