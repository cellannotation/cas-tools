import itertools
import json
import shutil
from importlib import resources
from typing import Any, Dict, List, Optional, Set, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import requests
from cas_schema import schemas

from cas.accession.hash_accession_manager import HashAccessionManager
from cas.dataset_retrieval.dataset_retriever import DatasetRetriever
from cas.file_utils import get_cas_schema_names

CROSSREF_API_URL = "https://api.crossref.org/works/"

LABELSET_NAME = "name"

LABELSET = "labelset"

LABELSETS = "labelsets"

ANNOTATIONS = "annotations"

CELL_IDS = "cell_ids"

CELL_LABEL = "cell_label"

CELL_SET_ACCESSION = "cell_set_accession"

PARENT_CELL_SET_ACCESSION = "parent_cell_set_accession"

NT_ACCESSION = "neurotransmitter_accession"

AUTHOR_ANNOTATION_FIELDS = "author_annotation_fields"

CELLHASH = "cellhash"


def calculate_labelset_rank(input_list: List[str]) -> Dict[str, int]:
    """
    Assign ranks to items in a list.

    Args:
        input_list (List[str]): The list of items.

    Returns:
        Dict[str, int]: A dictionary where keys are items from the input list and
        values are their corresponding ranks (0-based).

    """
    return {item: rank for rank, item in enumerate(input_list)}


def calculate_labelset(
    obs: pd.DataFrame, labelsets: List[str]
) -> Dict[str, Dict[str, Any]]:
    """
    Calculates labelset dictionary based on the provided observations.

    Args:
        obs (pd.DataFrame): DataFrame containing observations.
        labelsets (List[str]): List of labelsets.

    Returns:
        Dict[str, Dict[str, Any]]: A dictionary where keys are labelsets and values are dictionaries containing:
            - "members": set of members for the labelset.
            - "rank": rank of the labelset.
    """
    labelset_dict = {}
    # Assuming that terms at higher levels in the hierarchy have fewer members
    # compared to terms at lower levels.
    unique_counts = {col: obs[col].nunique() for col in labelsets}
    sorted_labelsets = sorted(unique_counts, key=unique_counts.get, reverse=True)
    labelset_rank_dict = calculate_labelset_rank(sorted_labelsets)
    for item in labelsets:
        labelset_dict.update(
            {
                item: {
                    "members": set(obs[item]),
                    "rank": str(labelset_rank_dict.get(item)),
                }
            }
        )
    return labelset_dict


def add_labelsets_to_cas(cas: Dict[str, Any], labelset_dict: Dict[str, Dict[str, Any]]):
    """
    Updates a CAS dictionary with labelsets derived from a labelset information dictionary.

    Args:
        cas: The CAS dictionary to update.
        labelset_dict: Contains labelset names as keys and dicts with 'rank' (and potentially other info) as values.

    """
    for labelset, value in labelset_dict.items():
        cas.get("labelsets").append(
            {"name": labelset, "description": "", "rank": str(value.get("rank"))}
        )


def get_cell_ids(
    anndata_obs: pd.DataFrame, labelset: str, cell_label: str
) -> List[str]:
    """
    Get cell IDs from an AnnData dataset based on a specified labelset and cell label.

    Args:
        anndata_obs: The observations DataFrame (`obs`) of an AnnData object containing the dataset. This DataFrame
            should include columns for cell type ontology term IDs and cell types.
        labelset: Labelset to filter.
        cell_label: The value of the cell label used to filter rows in `anndata_obs`.

    Returns:
        List[str]: List of cell IDs.
    """
    cell_label_lower = str(cell_label).lower()
    return anndata_obs.index[
        anndata_obs[labelset].astype(str).str.lower() == cell_label_lower
    ].tolist()


def get_cl_annotations_from_anndata(
    anndata_obs: pd.DataFrame, columns_name: str, cell_label: str
) -> Tuple[str, str]:
    """
    Retrieves cell ontology term ID and cell ontology term for a given cell label from the observation DataFrame
        of an AnnData object.

    Args:
        anndata_obs: The observations DataFrame (`obs`) of an AnnData object containing the dataset. This DataFrame
            should include columns for cell type ontology term IDs and cell types.
        columns_name: The name of the column in `anndata_obs` used for filtering based on the cell label.
        cell_label: The value of the cell label used to filter rows in `anndata_obs`.

    Returns:
        tuple: A tuple containing two elements:
            - The first element is the cell type ontology term ID associated with the given cell label.
            - The second element is the cell type (ontology term) associated with the given cell label.
    """
    filtered_df = anndata_obs[anndata_obs[columns_name] == cell_label]

    ontology_term_ids = filtered_df["cell_type_ontology_term_id"].unique().tolist()
    ontology_terms = filtered_df["cell_type"].unique().tolist()

    cell_ontology_term_id = (
        filtered_df["cell_type_ontology_term_id"].iloc[0]
        if len(ontology_term_ids) == 1
        else None
    )
    cell_ontology_term = (
        filtered_df["cell_type"].iloc[0] if len(ontology_terms) == 1 else None
    )

    return cell_ontology_term_id, cell_ontology_term


def collect_parent_cell_ids(cas: Dict[str, Any]) -> Dict[str, Set]:
    """
    Collects parent cell IDs from the given CAS data.

    This function iterates through labelsets in the CAS data and collects parent cell IDs
    associated with each labelset annotation. It populates and returns a dictionary
    mapping parent cell set accessions to sets of corresponding cell IDs.

    Args:
        cas: The Cell Annotation Schema data containing labelsets and annotations.

    Returns:
        A dictionary mapping parent cell set accessions to sets of corresponding cell IDs.
    """
    parent_cell_ids = dict()

    labelsets = sorted(
        [ls for ls in cas[LABELSETS] if "rank" in ls], key=lambda x: int(x["rank"])
    )
    for labelset in labelsets:
        ls_annotations = [
            ann for ann in cas[ANNOTATIONS] if ann[LABELSET] == labelset[LABELSET_NAME]
        ]

        for ann in ls_annotations:
            if "parent_cell_set_accession" in ann:
                cell_ids = set()
                if CELL_IDS in ann and ann[CELL_IDS]:
                    cell_ids = set(ann[CELL_IDS])
                elif (
                    "cell_set_accession" in ann
                    and ann["cell_set_accession"] in parent_cell_ids
                ):
                    cell_ids = parent_cell_ids[ann["cell_set_accession"]]

                if ann["parent_cell_set_accession"] in parent_cell_ids:
                    parent_cell_ids[ann["parent_cell_set_accession"]].update(cell_ids)
                else:
                    parent_cell_ids[ann["parent_cell_set_accession"]] = set(cell_ids)

    return parent_cell_ids


def generate_parent_cell_lookup(anndata, labelset_dict):
    """
    Generates a lookup dictionary mapping cell labels to various metadata, including cell IDs, rank,
    and cell ontology terms. This function is designed to precompute the lookup information needed for
    CAS annotation generation, especially useful when hierarchy inclusion is desired.

    Args:
        anndata (ad.AnnData): The AnnData object containing the single-cell dataset,
                              including metadata in anndata.obs.
        labelset_dict (Dict[str, Any]): A dictionary where keys are labelset names and values
                                        are dictionaries containing members and their ranks.

    Returns:
        Dict[str, Any]: A dictionary where each key is a cell label and each value is another
                        dictionary containing keys for 'cell_ids' (a set of cell IDs associated
                        with the label), 'rank', 'cell_ontology_term_id', and 'cell_ontology_term'.
    """
    accession_manager = HashAccessionManager()
    parent_cell_look_up = {}
    for k, v in labelset_dict.items():
        for label in v["members"]:
            cell_ontology_term_id, cell_ontology_term = get_cl_annotations_from_anndata(
                anndata.obs, k, label
            )
            cell_ids = get_cell_ids(anndata.obs, k, label)
            cell_set_accession = accession_manager.generate_accession_id(
                cell_ids=cell_ids, labelset=k
            )

            if label in parent_cell_look_up:
                parent_cell_look_up[f"{k}:{label}"]["cell_ids"].update(cell_ids)
            else:
                parent_cell_look_up[f"{k}:{label}"] = {
                    "cell_ids": set(cell_ids),
                    "accession": cell_set_accession,
                    "rank": v.get("rank"),
                    "cell_ontology_term_id": cell_ontology_term_id,
                    "cell_ontology_term": cell_ontology_term,
                }
    return parent_cell_look_up


def update_parent_info(
    value: Dict[str, Any], parent_key: str, parent_value: Dict[str, Any]
):
    """Updates parent information in a child item's dictionary.

    Args:
        value (Dict[str, Any]): The child item's dictionary to be updated.
        parent_key (str): The key of the parent item.
        parent_value (Dict[str, Any]): The parent item's dictionary.

    This function modifies `value` to include `parent` (using `parent_key`),
    `p_accession`, and `parent_rank` based on `parent_value`.
    """
    value.update(
        {
            "parent": parent_key.split(":")[-1],
            "p_accession": parent_value.get("accession"),
            "parent_rank": parent_value.get("rank"),
        }
    )


def add_parent_cell_hierarchy(parent_cell_look_up: Dict[str, Any]):
    """
    Processes parent cell hierarchy information and updates CAS dictionary annotations accordingly.

    Args:
        parent_cell_look_up (Dict[str, Any]): Dictionary containing parent cell information.

    Returns:
        None
    """
    # Establish parent-child relationships
    for (key, value), (inner_key, inner_value) in itertools.product(
        parent_cell_look_up.items(), repeat=2
    ):
        if (
            key == inner_key
            or value.get("parent")
            and value.get("parent_rank", 0) <= inner_value.get("rank")
        ):
            continue

        if value["cell_ids"] != inner_value["cell_ids"] and value["cell_ids"].issubset(
            inner_value["cell_ids"]
        ):
            update_parent_info(value, inner_key, inner_value)
        elif value["cell_ids"] == inner_value["cell_ids"]:
            if int(inner_value["rank"]) < int(value["rank"]):
                update_parent_info(inner_value, key, value)
            elif int(value["rank"]) < int(inner_value["rank"]):
                update_parent_info(value, inner_key, inner_value)
            else:
                raise ValueError(
                    f"{key} and {inner_key} cell labels have the same cell_ids. Cell_ids can't be "
                    f"identical at the same rank."
                )


def add_parent_hierarchy_to_annotations(
    cas: Dict[str, Any], parent_cell_look_up: Dict[str, Any]
):
    """
    Adds parent hierarchy information to annotations in the CAS dictionary.

    Args:
        cas (Dict[str, Any]): The CAS dictionary containing annotations.
        parent_cell_look_up (Dict[str, Any]): Dictionary containing parent cell information.

    Returns:
        None
    """
    annotation_list = cas.get("annotations", [])
    for annotation in annotation_list:
        parent_info = parent_cell_look_up.get(
            f'{annotation.get("labelset")}:{annotation.get("cell_label")}', {}
        )
        parent = parent_info.get("parent")
        p_accession = parent_info.get("p_accession")
        if parent and p_accession:
            # Add parent data to the annotation
            annotation.update(
                {
                    "parent_cell_set_name": parent,
                    "parent_cell_set_accession": p_accession,
                }
            )
            # Remove redundant CL annotations
            parent_dict = parent_cell_look_up.get(
                f"{p_accession.split(':')[0]}:{parent}", {}
            )
            if parent_dict.get("cell_ontology_term_id") == annotation.get(
                "cell_ontology_term_id"
            ) and parent_dict.get("cell_ontology_term") == annotation.get(
                "cell_ontology_term"
            ):
                annotation.pop("cell_ontology_term_id", None)
                annotation.pop("cell_ontology_term", None)


def get_authors_from_doi(doi):
    """
    Fetches and returns a list of authors from a given DOI (Digital Object Identifier) using the CrossRef API.

    Args:
        doi (str): The DOI of the publication for which to retrieve author information.

    Returns:
        list of dict: A list of dictionaries where each dictionary contains details of one author, including
                      their name ('author_name'), ORCID ID ('orcid'), GitHub username ('github_username'), and email ('email').
                      Each field is a string, and fields without data will be None.

    Raises:
        KeyError: If the author data is not found in the response, indicating a potential issue with the DOI or the data format.

    """
    response = requests.get(f"{CROSSREF_API_URL}{doi}")
    data = response.json()

    try:
        authors = data["message"]["author"]
        author_dict = [
            {
                "author_name": f"{author.get('given')} {author.get('family')}".strip(),
                "orcid": author.get("ORCID"),
                "github_username": author.get("github_username"),
                "email": author.get("email"),
            }
            for author in authors
        ]
        return author_dict
    except KeyError:
        return "Author information not available."


def reformat_json(
    input_json: Dict[str, Any],
    input_key: str = "annotations",
    exclude_key: str = "cell_ids",
) -> str:
    """
    Reformat the input JSON to create a new JSON structure, copying all fields and modifying the 'input_key' field.
    This function serializes the modified JSON to a string.

    Args:
        input_json: The original JSON object as a Python dictionary.
        input_key: The key in the original JSON where annotations are stored.
        exclude_key: The key within annotations to exclude from the copied data.

    Returns:
        A JSON string of the reformatted JSON object.
    """
    output_json = {k: v for k, v in input_json.items() if k != input_key}

    # Process the annotations field if it exists in the input JSON
    if input_key in input_json:
        # Filter out the specified key (e.g., 'cell_ids') from each annotation dictionary
        filtered_annotations = [
            {k: v for k, v in annotation.items() if k != exclude_key}
            for annotation in input_json[input_key]
        ]

        output_json[input_key] = filtered_annotations

    return json.dumps(output_json)


# Conversion function to handle complex types for JSON serialization
def convert_complex_type(value):
    """
    Converts all complex types to strings except for bool, int, float, and str.
    - Leaves bool types (including numpy.bool_) unchanged.
    - Converts everything else to strings.
    """
    if isinstance(
        value, (bool, int, float, str)
    ):  # Leave bool, int, float, and str unchanged
        return value
    elif isinstance(
        value, np.bool_
    ):  # Special case to convert numpy bool to Python bool
        return bool(value)
    else:  # Convert everything else to string
        return str(value)


def copy_and_update_file_path(
    anndata_file_path: str, output_file_path: Optional[str]
) -> str:
    """
    Copies the AnnData file to a new location if an output file path is provided, and updates the file path.

    Args:
        anndata_file_path: The path to the original AnnData file.
        output_file_path: The path to which the file should be copied. If not provided, no copying occurs.

    Returns:
        str: The updated file path. If `output_file_path` is provided, it will return the new path,
        otherwise the original `anndata_file_path`.
    """
    if output_file_path:
        shutil.copy(anndata_file_path, output_file_path)
        anndata_file_path = output_file_path
    return anndata_file_path


def fetch_anndata(
    input_json: Dict[str, Any], download_dir: Optional[str] = None
) -> str:
    """
    Fetches the AnnData file based on the provided CAS JSON input.

    Args:
        input_json: A dictionary containing CAS JSON data. Must include a "matrix_file_id" key.
        download_dir: The directory where the AnnData file should be downloaded.
                      If not provided, the current working directory is used.

    Returns:
        str: The path to the downloaded AnnData file.

    Raises:
        KeyError: If the "matrix_file_id" key is missing from the `input_json`.
    """
    matrix_file_id: Optional[str] = input_json.get("matrix_file_id", None)
    if matrix_file_id:
        dataset_retriever = DatasetRetriever.create(matrix_file_id)
        anndata_file_path = dataset_retriever.download_data(download_dir=download_dir)
    else:
        raise KeyError("Matrix file id is missing from CAS json.")
    return anndata_file_path


def retrieve_schema(schema_name):
    schema_name = str(schema_name).strip().lower()
    if schema_name not in get_cas_schema_names():
        raise Exception("Schema name should be one of 'base', 'bican' or 'cap'")
    schema_file = resources.files(schemas) / get_cas_schema_names()[schema_name]
    with schema_file.open("rt") as f:
        schema = json.loads(f.read())
    return schema
