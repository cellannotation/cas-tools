import os
import unittest

from cas.dataset_retrieval.cxg_downloader import CxGDownloader


@unittest.skip("Test case will be completed in the future.")
class TestCxGDownloader(unittest.TestCase):
    def test_download_data(self):
        matrix_field_id = "cxg_dataset:24ec2dc5-3573-4d66-a9e1-25b7dcf43e27"
        file_path = "src/test/test_data/downloader_test.h5ad"
        cxg_downloader = CxGDownloader(matrix_field_id)
        cxg_downloader.download_data(file_name=file_path)
        # Assert the file exists
        self.assertTrue(os.path.exists(file_path), "Downloaded file does not exist.")
        # Assert the file is not empty
        self.assertGreater(os.path.getsize(file_path), 0, "Downloaded file is empty.")

        # Clean up the test file after running
        if os.path.exists(file_path):
            os.remove(file_path)


if __name__ == '__main__':
    unittest.main()
