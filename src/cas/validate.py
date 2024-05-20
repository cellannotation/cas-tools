import json
import os
import warnings
from importlib import resources

from cas_schema import schema_validator, schemas
from cas.file_utils import get_cas_schema, get_cas_schema_names

warnings.filterwarnings("always")


def validate(schema_name: str, data_path: str):
    """
    Validates all instances in data_path against the given schema.
    Assumes all *.json files in the test_dir should validate against the schema.
    Logs all validation errors and throws an exception if any of the test files is invalid.
    Parameters:
        schema_name: One of 'base', 'bican' or 'cap'. Identifies the CAS schema to validate data against.
        data_path: Path to the data file (or folder) to validate
    """
    schema_name = schema_name.strip().lower()
    schema = get_cas_schema(schema_name)
    if not os.path.exists(data_path):
        raise Exception("Please provide a valid 'data_path': {}".format(data_path))

    result = schema_validator.validate(
        schema, get_cas_schema_names()[schema_name], data_path
    )
    if not result:
        raise Exception("Validation Failed")

    return True
