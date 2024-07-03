import json
import logging
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd

from cas.file_utils import read_json_file


# Configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def add_author_annotations(
    cas: Dict[str, Any],
    df: pd.DataFrame,
    join_column: Union[str, List[str]],
    columns: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Updates the provided CAS dictionary with author annotation fields from a DataFrame.
    Annotations are based on the values matched between the DataFrame's specified join columns and the CAS records.
    All columns or specified columns are added as annotations.

    Args:
        cas: The CAS dictionary loaded from a JSON file.
        df: DataFrame containing the data to join with the CAS.
        join_column: The column(s) in the DataFrame used for matching CAS records; can be a single column or a list of columns.
        columns: Optional list of columns whose data will be added as annotations. If None, all DataFrame columns are used.

    Returns:
        Dict[str, Any]: The CAS dictionary updated with new annotation fields.
    """
    validate_columns(df, join_column, columns)
    validate_values(df, join_column, cas)
    annotations = cas["annotations"]
    for annotation in annotations:
        labelset = annotation["labelset"]
        filter_value = (
            annotation[join_column]
            if isinstance(join_column, str)
            else annotation["cell_label"]
        )
        author_annotation_dict = dataframe_to_dict(
            df, join_column, labelset, filter_value, columns
        )
        if "author_annotation_fields" in annotation:
            annotation["author_annotation_fields"].update(author_annotation_dict)
        else:
            annotation["author_annotation_fields"] = author_annotation_dict

    return cas


def add_author_annotations_from_file(
    cas_json: str,
    csv_path: str,
    join_column: Union[str, List[str]],
    columns: Optional[List[str]] = None,
    output_file_path: str = "output.json",
) -> None:
    """
    Reads data from a CSV file and a CAS JSON file, then updates the CAS JSON with author annotation fields,
    and outputs the updated CAS JSON to a specified file. It uses specified columns or all columns if none are
    specified.

    Args:
        cas_json: Path to the CAS JSON file.
        csv_path: Path to the CSV file.
        join_column: Column name or names used for matching CAS records.
        columns: Optional columns to be added as annotations; if None, all columns are used.
        output_file_path: Output CAS file name.

    Returns:
        Dict[str, Any]: The updated CAS JSON dictionary.
    """
    cas = read_json_file(cas_json)
    df = pd.read_csv(csv_path)
    output_cas = add_author_annotations(cas, df, join_column, columns)
    # Write the JSON data to the file
    with open(output_file_path, "w") as json_file:
        json.dump(output_cas, json_file, indent=2)


def validate_columns(
    df: pd.DataFrame, join_column: Union[str, List[str]], columns: Optional[List[str]]
):
    """
    Validates the existence of key columns and optionally other specified columns in a DataFrame.

    Args:
        df: DataFrame to check.
        join_column: Column or columns that must exist in the DataFrame.
        columns: Additional column names that must exist if specified; checks only join_column if None.

    Raises:
        ValueError: If any specified columns, including the join column(s), are not found in the DataFrame.
    """
    # Normalize join_column to always be a list
    if isinstance(join_column, str):
        join_column = [join_column]
    column_list_to_check = join_column + (columns if columns else [])
    missing_columns = [col for col in column_list_to_check if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing columns in DataFrame: {', '.join(missing_columns)}")


def validate_values(
    df: pd.DataFrame, join_column: Union[str, List[str]], cas: Dict[str, Any]
):
    """
    Validates that all keys generated from specified DataFrame columns exist within the CAS annotations.
    If there are any keys in the DataFrame that are not found in the CAS annotations, it raises an exception
    indicating these are extra, undesired keys.

    Parameters:
        df (pd.DataFrame): DataFrame containing the data to validate.
        join_column (Union[str, List[str]]): Column or list of columns in the DataFrame used to form the keys.
            If a list of two strings is provided, these columns are concatenated to form the keys.
        cas (Dict[str, Any]): Dictionary containing CAS annotations, each annotation expected to contain the keys.

    Raises:
        ValueError: If there are extra keys in the DataFrame that do not exist in the CAS annotations.
    """
    cas_key_list = []
    for annotation in cas["annotations"]:
        if isinstance(join_column, list) and len(join_column) == 2:
            key_value = annotation.get(join_column[0], "") + annotation.get(
                join_column[1], ""
            )
            cas_key_list.append(key_value)
        else:
            cas_key_list.append(annotation.get(join_column, ""))

    if isinstance(join_column, list) and len(join_column) == 2:
        df_key_list = df[join_column[0]].str.cat(df[join_column[1]], sep="").tolist()
    else:
        df_key_list = df[join_column].tolist()

    extra_keys = [key for key in df_key_list if key not in cas_key_list]
    empty_rows = [key for key in cas_key_list if key not in df_key_list]
    if extra_keys:
        raise ValueError(
            f"Extra keys in DataFrame that are not in CAS data: {extra_keys}"
        )
    if empty_rows:
        logging.info(f"Following values in {' '.join(join_column) if isinstance(join_column, list) else join_column} exist in "
                     f"CAS data but missing from DataFrame:"
                     f" {' '.join(empty_rows)}")


def dataframe_to_dict(
    df: pd.DataFrame,
    join_column: Union[str, List[str]],
    labelset: str,
    filter_value: str,
    columns: Optional[List[str]],
) -> dict:
    """
    Converts specified columns of a DataFrame into a dictionary. If columns are not specified, the entire DataFrame
    for rows matching the condition is converted.

    Args:
        df: DataFrame to extract data from.
        join_column: The column or columns used to filter the DataFrame based on cell_label.
        labelset: Labelset of a CAS annotation.
        filter_value: The value in join_column used to filter rows.
        columns: List of column names to convert to dictionary. If None, all columns are used.

    Returns:
        dict: Dictionary with column names as keys and single values or lists of values as values.
    """
    # Handle different cases for join_column being a string or list of strings
    if join_column == "cell_set_accession":
        filtered_df = df[df[join_column] == filter_value]
    elif isinstance(join_column, str):
        filtered_df = df[df[join_column] == filter_value]
    elif isinstance(join_column, list) and len(join_column) == 2:
        # Expect cell_label to be a tuple or list with corresponding values
        filtered_df = df[
            (df[join_column[0]] == labelset) & (df[join_column[1]] == filter_value)
        ]
    else:
        raise ValueError(
            "join_column must be either a string or a list of two strings."
        )

    # Determine columns to include in the output dictionary
    if columns is None:
        # Exclude join_column(s) from the output
        exclude_columns = (
            join_column if isinstance(join_column, list) else [join_column]
        )
        columns = [col for col in df.columns if col not in exclude_columns]

    # Select the specified columns to include
    filtered_df = filtered_df[columns]

    filtered_df.replace(' ', np.nan, inplace=True)
    filtered_df.replace(np.nan, None, inplace=True)
    result_dict = filtered_df.to_dict(orient="list")
    # Unpack lists that contain only a single item
    for key, value in result_dict.items():
        if len(value) == 1:
            result_dict[key] = value[0]

    return result_dict
