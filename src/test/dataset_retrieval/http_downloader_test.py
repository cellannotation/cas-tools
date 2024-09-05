import os
import unittest

from cas.dataset_retrieval.http_downloader import HTTPDownloader


class TestHTTPDownloader(unittest.TestCase):
    def test_download_data(self):
        matrix_field_id = "https://datasets.cellxgene.cziscience.com/aaab3abd-624a-442e-b62b-3f2edb10b45e.h5ad"
        file_path = os.path.join(os.curdir, "downloader_test.h5ad")
        cxg_downloader = HTTPDownloader(matrix_field_id)
        cxg_downloader.download_data(file_name=file_path)
        # Assert the file exists
        self.assertTrue(os.path.exists(file_path), "Downloaded file does not exist.")
        # Assert the file is not empty
        self.assertGreater(os.path.getsize(file_path), 0, "Downloaded file is empty.")

        # Clean up the test file after running
        if os.path.exists(file_path):
            os.remove(file_path)


if __name__ == "__main__":
    unittest.main()
