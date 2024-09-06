import os
import unittest

from cas.dataset_retrieval.cxg_downloader import CxGDownloader


class TestCxGDownloader(unittest.TestCase):
    def test_download_data(self):
        matrix_field_id = "0895c838-e550-48a3-a777-dbcd35d30272"
        file_path = "downloader_test.h5ad"
        download_dir = "../test_data"
        cxg_downloader = CxGDownloader(matrix_field_id)
        full_download_path = cxg_downloader.download_data(
            file_name=file_path, download_dir=download_dir
        )
        # Assert the file exists
        self.assertTrue(
            os.path.exists(full_download_path), "Downloaded file does not exist."
        )
        # Assert the file is not empty
        self.assertGreater(
            os.path.getsize(full_download_path), 0, "Downloaded file is empty."
        )

        # Clean up the test file after running
        if os.path.exists(full_download_path):
            os.remove(full_download_path)


if __name__ == "__main__":
    unittest.main()
