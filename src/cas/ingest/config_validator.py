import json
import os
import warnings

from jsonschema import Draft7Validator
from ruamel.yaml import YAML

from cas.file_utils import read_config

SCHEMA_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "./config_schema.yaml"
)


ryaml = YAML(typ="safe")
with open(SCHEMA_PATH) as stream:
    ctat_schema = ryaml.load(stream)
    validator = Draft7Validator(ctat_schema)


def validate(json_object: object) -> bool:
    """
    Validates the given json configuration object using the cell type annotation schema.

    Returns:
    :param json_object: configuration object
    :return: True if object is valid, False otherwise.
    """
    is_valid = True

    if not validator.is_valid(json_object):
        es = validator.iter_errors(json_object)
        for e in es:
            warnings.warn(str(e.message))
            is_valid = False

    return is_valid


def validate_file(file_path: str) -> bool:
    """
    Read the configuration object from the given path and validates it.
    :param file_path: path to the json file
    :return: True if object is valid, False otherwise.
    """
    return validate(read_config(file_path))


def validate_json_str(json_str: str) -> bool:
    """
    Validates the given json string.

    Returns:
    :param json_str: string representation of a json object
    :return: True if object is valid, False otherwise.
    """
    return validate(json.loads(json_str))
