import os
import unittest

from cas.ingest.config_validator import validate, validate_file, validate_json_str

VALID_TEST_DATA_1 = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "./test_data/nhp_basal_ganglia/test_config.yaml",
)
INVALID_TEST_DATA_1 = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "./test_data/nhp_basal_ganglia/test_config_invalid.yaml",
)


class SchemaValidationTests(unittest.TestCase):
    def test_validator_yaml(self):
        self.assertTrue(validate_file(VALID_TEST_DATA_1))
        self.assertFalse(validate_file(INVALID_TEST_DATA_1))
