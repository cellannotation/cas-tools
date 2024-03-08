import json
import os
from importlib.metadata import version
from typing import Any, Dict

import pandas as pd

from cas.file_utils import read_json_file

CAT_SET_REQUIRED_COLUMNS = ["label", "description", "order"]
CAT_REQUIRED_COLUMNS = [
    "cluster_annotation_term_set_label",
    "name",
    "label",
    "parent_term_label",
]


def validate_dataframe_columns(df: pd.DataFrame, required_columns: list):
    """
    Validates a DataFrame for the required columns.
    Args:
        df (pandas.DataFrame): DataFrame to validate.
        required_columns (list): List of required column names.
    """
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing columns in DataFrame: {', '.join(missing_columns)}")


def generate_catset_dataframe(cas: Dict[str, Any]) -> pd.DataFrame:
    """
    Generate a DataFrame representing the Cluster Annotation Term Set (cat_set)
    from the given Cell Annotation Schema (CAS) dictionary.

    Args:
        cas (Dict[str, Any]): The Cell Annotation Schema (CAS) dictionary.

    Returns:
        pd.DataFrame: DataFrame representing the Cluster Annotation Term Set (cat_set).

    """
    labelsets = cas.get("labelsets", [])
    data = [
        {
            "name": labelset.get("name"),
            "description": labelset.get("description"),
            "order": labelset.get("rank"),
        }
        for labelset in labelsets
    ]
    catset_df = pd.DataFrame(data)
    return catset_df


def generate_cat_dataframe(cas: Dict[str, Any]) -> pd.DataFrame:
    """
    Generate a DataFrame representing the Cluster Annotation Term (cat) from the given Cell Annotation Schema (CAS)
    dictionary.

    Args:
        cas (Dict[str, Any]): The Cell Annotation Schema (CAS) dictionary.

    Returns:
        pd.DataFrame: DataFrame representing the Cluster Annotation Term (cat).
    """
    columns = [
        "label",
        "name",
        "cluster_annotation_term_set_label",
        "parent_term_label",
        "parent_term_set_label",
        "term_set_order",
        "term_order",
        "cluster_annotation_term_set_name",
    ]

    annotations = cas.get("annotations")
    data = []
    for annotation in annotations:
        cell_label = annotation.get("cell_label")
        cell_set_accession = annotation.get("cell_set_accession")
        parent_cell_set_accession = annotation.get("parent_cell_set_accession")
        new_row = {
            "label": cell_set_accession,
            "name": cell_label,
            "parent_term_label": parent_cell_set_accession,
        }
        data.append(new_row)

    cat_df = pd.DataFrame(data, columns=columns)
    return cat_df


def calculate_order_mapping(order_values: pd.Series) -> Dict[str, str]:
    """
    Calculate a mapping dictionary based on the order values.

    Args:
        order_values (pandas.Series): Series containing the order values.

    Returns:
        Dict[str, str]: Mapping dictionary where keys are order values and values are rank values.
    """
    order_values_filtered = order_values[order_values != 0]
    order_values_sorted = order_values_filtered.sort_values(ascending=False)
    mapping = dict(zip(order_values_sorted, range(len(order_values_sorted))))

    return mapping


def abc2cas(cat_set_file_path: str, cat_file_path: str, output_file_path: str):
    """
    Converts given ABC files to a Cell Annotation Schema (CAS) JSON and writes it to a file with output_file_path name.
    Args:
        cat_set_file_path: Path to the Cluster Annotation Term Set file.
        cat_file_path: Path to the Cluster Annotation Term file.
        output_file_path: Output CAS file name (default: output.json).

    """
    cat_set = pd.read_csv(cat_set_file_path, sep=",")
    cat = pd.read_csv(cat_file_path, sep=",")

    validate_dataframe_columns(cat_set, CAT_SET_REQUIRED_COLUMNS)
    validate_dataframe_columns(cat, CAT_REQUIRED_COLUMNS)

    cas = init_metadata()
    add_labelsets(cas, cat_set)
    add_annotations(cas, cat)

    # Write the JSON data to the file
    with open(output_file_path, "w") as json_file:
        json.dump(cas, json_file, indent=2)


def add_annotations(cas: Dict[str, Any], cat: pd.DataFrame):
    """
    Adds annotations to the Cell Annotation Schema (CAS) based on the data from the Cluster Annotation Term DataFrame.

    Args:
        cas (Dict[str, Any]): Dictionary representing the Cell Annotation Schema.
        cat (pd.DataFrame): DataFrame containing Cluster Annotation Term data.

    """
    for row in cat.itertuples():
        labelset = row.cluster_annotation_term_set_label
        cell_label = row.name
        cell_ontology_term_id = None
        cell_ontology_term = None
        cell_ids = None
        rationale = None
        rationale_dois = None
        marker_gene_evidence = None
        synonyms = None
        category_fullname = None
        category_cell_ontology_exists = None
        category_cell_ontology_term_id = None
        category_cell_ontology_term = None
        cell_set_accession = row.label
        parent_cell_set_accession = row.parent_term_label

        anno = {
            "labelset": labelset,
            "cell_label": cell_label,
            "cell_fullname": cell_label,
            "cell_ontology_term_id": cell_ontology_term_id,
            "cell_ontology_term": cell_ontology_term,
            "cell_ids": cell_ids,
            "rationale": rationale,
            "rationale_dois": rationale_dois,
            "marker_gene_evidence": marker_gene_evidence,
            "synonyms": synonyms,
            "category_fullname": category_fullname,
            "category_cell_ontology_exists": category_cell_ontology_exists,
            "category_cell_ontology_term_id": category_cell_ontology_term_id,
            "category_cell_ontology_term": category_cell_ontology_term,
            "cell_set_accession": cell_set_accession,
            "parent_cell_set_accession": parent_cell_set_accession,
        }
        cas.get("annotations").append({k: v for k, v in anno.items() if v is not None and not pd.isna(v)})


def add_labelsets(cas: Dict[str, Any], cat_set: pd.DataFrame):
    """
    Adds labelsets to the Cell Annotation Schema (CAS) based on the data from the Cluster Annotation Term Set DataFrame.

    Args:
        cas (Dict[str, Any]): Cell Annotation Schema dictionary.
        cat_set (pandas.DataFrame): DataFrame containing Cluster Annotation Term Set data.
    """
    order_mapping = calculate_order_mapping(cat_set["order"])
    for row in cat_set.itertuples():
        name = row.label
        description = row.description
        rank = None
        if row.order != 0:
            rank = order_mapping[row.order]

        labelset = {"name": name, "description": description}
        if rank is not None:
            labelset["rank"] = rank
        cas.get("labelsets").append(labelset)


def init_metadata() -> Dict[str, Any]:
    """
    Initializes metadata for Cell Annotation Schema (CAS).

    Returns:
        Dict[str, Any]: Metadata dictionary containing default values for various fields.
    """
    # TODO These needs proper assignments
    matrix_file_id = None
    cellannotation_schema_version = version("cell-annotation-schema")
    cellannotation_timestamp = None
    cellannotation_version = None
    cellannotation_url = None
    author_name = "Jane Doe"  # TODO Needs proper author_name assignment
    author_contact = None
    orcid = None
    cas = {
        "matrix_file_id": matrix_file_id,
        "cellannotation_schema_version": cellannotation_schema_version,
        "cellannotation_timestamp": cellannotation_timestamp,
        "cellannotation_version": cellannotation_version,
        "cellannotation_url": cellannotation_url,
        "author_name": author_name,
        "author_contact": author_contact,
        "orcid": orcid,
        "annotations": [],
        "labelsets": [],
    }
    cas = {k: v for k, v in cas.items() if v is not None}
    return cas


def cas2abc(cas_file_path: str, cat_set_file_path: str, cat_file_path: str):
    """
    Converts given Cell Annotation Schema (CAS) to ABC files: cluster_annotation_term and
    cluster_annotation_term_set, and writes them to files with cat_file_path and cat_set_file_path.

    Args:
        cas_file_path: Path to the Cell Annotation Schema (CAS) file
        cat_set_file_path: Path to the Cluster Annotation Term Set file.
        cat_file_path: Path to the Cluster Annotation Term file.

    """
    cas_json = read_json_file(cas_file_path)

    cluster_annotation_term_set = generate_catset_dataframe(cas_json)
    cluster_annotation_term = generate_cat_dataframe(cas_json)

    current_directory = os.getcwd()
    cluster_annotation_term_set.to_csv(
        os.path.join(current_directory, cat_set_file_path), index=False
    )
    cluster_annotation_term.to_csv(
        os.path.join(current_directory, cat_file_path), index=False
    )
    # TODO implement rest of the method once the requirements are more clear
