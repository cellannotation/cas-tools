import csv
import json
import anndata
import pathlib

from typing import Optional
from ruamel.yaml import YAML

from cas.model import CellTypeAnnotation


def read_json_file(file_path):
    """
    Reads and parses a JSON file into a Python dictionary.

    Args:
        file_path (str): The path to the JSON file.

    Returns:
        dict: The JSON data as a Python dictionary.

    Returns None if the file does not exist or if there is an issue
    parsing the JSON content.

    Example:
        json_data = read_json_file('path/to/your/file.json')
        if json_data is not None:
            # Use the parsed JSON data as a dictionary
            print(json_data)
    """
    try:
        with open(file_path, "r") as file:
            data = json.load(file)
            return data
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error reading JSON file: {e}")
        return None


def read_cas_json_file(file_path) -> CellTypeAnnotation:
    """
    Reads and parses a JSON file into a CAS object.

    Args:
        file_path (str): The path to the JSON file.

    Returns:
        dict: The JSON data as a CAS object.
    """
    return CellTypeAnnotation.from_dict(read_json_file(file_path))


def write_json_file(cas, out_file, print_undefined=False):
    """
    Writes cell type annotation object to a json file.
    :param cas: cell type annotation object to serialize.
    :param out_file: output file path.
    :param print_undefined: prints null values to the output json if true. Omits undefined values from the json output if
    """
    cas.set_exclude_none_values(not print_undefined)

    output_data = cas.to_json(indent=2)
    with open(out_file, "w") as out_file:
        out_file.write(output_data)


def read_anndata_file(file_path: str) -> Optional[anndata.AnnData]:
    """Load anndata object from a file.

    Args:
        file_path: The path to the file containing the anndata object.

    Returns:
        The loaded anndata object if successful, else None.
    """
    try:
        anndata_obj = anndata.read_h5ad(file_path, backed="r")
        return anndata_obj
    except Exception as e:
        print(f"An error occurred while loading the file: {e}")
        return None


def read_tsv_to_dict(tsv_path, id_column=0, generated_ids=False):
    """
    Reads tsv file content into a dict. Key is the first column value and the value is dict representation of the
    row values (each header is a key and column value is the value).
    Args:
        tsv_path: Path of the TSV file
        id_column: Id column becomes the key of the dict. This column should be unique. Default value is first column.
        generated_ids: If 'True', uses row number as the key of the dict. Initial key is 0.
    Returns:
        Function provides two return values: first; headers of the table and second; the TSV content dict. Key of the
        content is the first column value and the values are dict of row values.
    """
    return read_csv_to_dict(tsv_path, id_column=id_column, delimiter="\t", generated_ids=generated_ids)


def read_csv_to_dict(csv_path, id_column=0, id_column_name="", delimiter=",", id_to_lower=False, generated_ids=False):
    """
    Reads tsv file content into a dict. Key is the first column value and the value is dict representation of the
    row values (each header is a key and column value is the value).
    Args:
        csv_path: Path of the CSV file
        id_column: Id column becomes the keys of the dict. This column should be unique. Default is the first column.
        id_column_name: Alternative to the numeric id_column, id_column_name specifies id_column by its header string.
        delimiter: Value delimiter. Default is comma.
        id_to_lower: applies string lowercase operation to the key
        generated_ids: If 'True', uses row number as the key of the dict. Initial key is 0.

    Returns:
        Function provides two return values: first; headers of the table and second; the CSV content dict. Key of the
        content is the first column value and the values are dict of row values.
    """
    records = dict()

    headers = []
    with open(csv_path) as fd:
        rd = csv.reader(fd, delimiter=delimiter, quotechar='"')
        row_count = 0
        for row in rd:
            _id = str(row[id_column]).strip()
            if id_to_lower:
                _id = str(_id).strip().lower()

            if generated_ids:
                _id = row_count

            if row_count == 0:
                headers = [str(header).strip() for header in row]
                if id_column_name and id_column_name in headers:
                    id_column = headers.index(id_column_name)
            else:
                row_object = dict()
                for column_num, column_value in enumerate(row):
                    row_object[headers[column_num]] = column_value
                records[_id] = row_object

            row_count += 1

    return headers, records


def read_json_config(file_path: str) -> dict:
    """
    Reads the configuration object from the given path.
    :param file_path: path to the json file
    :return: configuration object (List of data column config items)
    """
    with open(file_path, "r") as fs:
        try:
            return json.load(fs)
        except Exception as e:
            raise Exception("JSON read failed:" + file_path + " " + str(e))


def read_yaml_config(file_path: str) -> dict:
    """
    Reads the configuration object from the given path.
    :param file_path: path to the yaml file
    :return: configuration object (List of data column config items)
    """
    with open(file_path, "r") as fs:
        try:
            ryaml = YAML(typ='safe')
            return ryaml.load(fs)
        except Exception as e:
            raise Exception("Yaml read failed:" + file_path + " " + str(e))


def read_config(file_path: str) -> dict:
    """
    Reads the configuration object from the given path.
    :param file_path: path to the configuration file
    :return: configuration object (List of data column config items)
    """
    file_extension = pathlib.Path(file_path).suffix
    if file_extension == ".json":
        return read_json_config(file_path)
    elif file_extension == ".yaml" or file_extension == ".yml":
        return read_yaml_config(file_path)
    else:
        raise Exception("Given configuration file extension is not supported. "
                        "Try a json or yaml file instead of :" + file_path)