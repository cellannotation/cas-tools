import json
import os
import warnings
from importlib import resources

from cas_schema import schema_validator, schemas

warnings.filterwarnings("always")

SCHEMA_FILE_MAPPING = {
    "base": "general_schema.json",
    "cap": "CAP_schema.json",
    "bican": "BICAN_schema.json",
}


def validate(schema_name: str, data_path: str):
    """
    Validates all instances in data_path against the given schema.
    Assumes all *.json files in the test_dir should validate against the schema.
    Logs all validation errors and throws an exception if any of the test files is invalid.
    Parameters:
        schema_name: One of 'base', 'bican' or 'cap'. Identifies the CAS schema to validate data against.
        data_path: Path to the data file (or folder) to validate
    """
    schema_name = str(schema_name).strip().lower()
    if schema_name not in SCHEMA_FILE_MAPPING:
        raise Exception("Schema name should be one of 'base', 'bican' or 'cap'")
    if not os.path.exists(data_path):
        raise Exception("Please provide a valid 'data_path': {}".format(data_path))

    schema_file = resources.files(schemas) / SCHEMA_FILE_MAPPING[schema_name]
    with schema_file.open("rt") as f:
        schema = json.loads(f.read())

    result = schema_validator.validate(
        schema, SCHEMA_FILE_MAPPING[schema_name], data_path
    )
    if not result:
        raise Exception("Validation Failed")

    return True
