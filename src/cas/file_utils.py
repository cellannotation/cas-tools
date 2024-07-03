import csv
import json
import pathlib
from typing import Optional
from importlib import resources

from cas_schema import schemas

import anndata
import numpy as np
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


def read_cas_json_file(file_path: str) -> CellTypeAnnotation:
    """
    Reads and parses a JSON file into a CAS object.

    Args:
        file_path (str): The path to the JSON file.

    Returns:
        dict: The JSON data as a CAS object.
    """
    return CellTypeAnnotation.from_dict(read_json_file(file_path))


def read_cas_from_anndata(anndata_path: str) -> CellTypeAnnotation:
    """
    Reads the CAS json from the anndata uns and parses into a CAS object.
    Args:
        anndata_path: The path to the Anndata file.

    Returns:
        CellTypeAnnotation object.
    """
    input_anndata = read_anndata_file(anndata_path)
    if input_anndata and "cas" in input_anndata.uns:
        return CellTypeAnnotation.from_dict(json.loads(input_anndata.uns["cas"]))
    else:
        raise Exception("Given Anndata file doesn't have a 'cas' object in it's uns.")


def write_json_file(
    cas: CellTypeAnnotation, out_file: str, print_undefined: bool = False
):
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


def write_dict_to_json_file(output_file_path: str, dictionary: dict):
    with open(output_file_path, "w") as json_file:
        json.dump(dictionary, json_file, indent=2)


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


def read_table_to_dict(table_path, id_column=0, generated_ids=False):
    """
    Reads table file content into a dict. Key is the first column value and the value is dict representation of the
    Args:
        table_path: Path of the table file
        id_column:  Id column becomes the key of the dict. This column should be unique. Default value is first column.
        generated_ids: If 'True', uses row number as the key of the dict. Initial key is 0.
    Returns:
        Function provides two return values: first; headers of the table and second; the TSV content dict. Key of the
        content is the first column value and the values are dict of row values.
    """
    if table_path.endswith(".tsv"):
        return read_tsv_to_dict(table_path, id_column=id_column, generated_ids=generated_ids)
    elif table_path.endswith(".csv"):
        return read_csv_to_dict(table_path, id_column=id_column, generated_ids=generated_ids)
    else:
        raise Exception("Table file should be either tsv or csv file.")


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
    return read_csv_to_dict(
        tsv_path, id_column=id_column, delimiter="\t", generated_ids=generated_ids
    )


def read_csv_to_dict(
    csv_path,
    id_column=0,
    id_column_name="",
    delimiter=",",
    id_to_lower=False,
    generated_ids=False,
):
    """
    Reads tsv file content into a dict. Key is the first column value and the value is dict representation of the
    row values (each header is a key and column value is the value).
    Args:
        csv_path: Path of the CSV file
        id_column: Id column becomes the keys of the dict. This column should be unique. Default is the first column.
        id_column_name: Alternative to the numeric id_column, id_column_name specifies id_column by its header string.
        delimiter: Value delimiter. Default is comma.
        id_to_lower: applies string lowercase operation to the key
        generated_ids: If 'True', uses row number as the key of the dict. Initial key is 1.

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
            ryaml = YAML(typ="safe")
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
        raise Exception(
            "Given configuration file extension is not supported. "
            "Try a json or yaml file instead of :" + file_path
        )


def update_obs_dataset(obs_dataset, flatten_data):
    """
    Updates obs dataset with flattened data.

    Args:
        obs_dataset (h5py.Dataset): Dataset representing the obs field in the AnnData file.
        flatten_data (dict): Dictionary containing flattened data.
    """
    for key, value in flatten_data.items():
        if key in obs_dataset:
            del obs_dataset[key]
            columns = obs_dataset.attrs["column-order"]
            columns = columns[columns != key]
            obs_dataset.attrs["column-order"] = columns

        obs_dataset.create_dataset(key, data=value.values.astype("O"))

        if key not in obs_dataset.attrs["column-order"]:
            columns = np.append(obs_dataset.attrs["column-order"], key)
            obs_dataset.attrs["column-order"] = columns


def write_json_to_hdf5(group, data):
    """
    Recursively writes JSON-like data to an HDF5 group.

    Args:
        group (h5py.Group): The HDF5 group to write data to.
        data (dict): A dictionary containing the data to be written.

    Returns:
        None
    """
    for key, value in data.items():
        if key in group:
            del group[key]  # Delete the existing key to overwrite it
        if isinstance(value, dict):
            subgroup = group.create_group(key)
            write_json_to_hdf5(subgroup, value)
        elif isinstance(value, list):
            if all(isinstance(item, str) for item in value):
                group.create_dataset(key, data=", ".join(sorted(value)))
            else:
                subgroup = group.create_group(key)
                for i, item in enumerate(value):
                    subgroup.create_dataset(str(i), data=item)
        else:
            group.create_dataset(key, data=value)


def get_cas_schema_names() -> dict:
    """
    Returns the list of available CAS schema names.

    Returns:
        dict: The available CAS schema names.
    """
    return {
        "base": "general_schema.json",
        "cap": "CAP_schema.json",
        "bican": "BICAN_schema.json",
    }


def get_cas_schema(schema_name: Optional[str] = "base") -> dict:
    """
    Reads the schema file from the CAS module and returns as a dictionary.
    Args:
        schema_name: The name of the schema to be returned. Default is 'base'.

    Returns:
        dict: The schema as a dictionary.
    """
    if not schema_name:
        schema_name = "base"

    schema_name = schema_name.strip().lower()
    if schema_name not in get_cas_schema_names():
        raise ValueError(
            "Schema name should be one of: " + ", ".join(get_cas_schema_names().keys())
        )

    schema_file = resources.files(schemas) / get_cas_schema_names()[schema_name]
    with schema_file.open("rt") as f:
        schema = json.loads(f.read())
    return schema
