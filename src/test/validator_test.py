import unittest
import os
import warnings

from cas.validate import validate

warnings.filterwarnings("always")

TEST_JSON = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                         "./test_data/validation/BICAN_schema_specific_examples/Silletti_annotation_transfer.json")
TEST_FOLDER = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           "./test_data/validation/BICAN_schema_specific_examples")


class ValidatorTests(unittest.TestCase):

    def test_file_validation(self):
        try:
            validate("bican", TEST_JSON)
            self.fail("A validation error should occur here")
        except Exception as e:
            self.assertEqual("Validation Failed", e.args[0])

    def test_folder_validation(self):
        try:
            validate("bican", TEST_FOLDER)
            self.fail("A validation error should occur here")
        except Exception as e:
            self.assertEqual("Validation Failed", e.args[0])

