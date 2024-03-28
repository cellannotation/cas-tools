import json
from typing import Any, Dict, List

from cas.file_utils import read_anndata_file
from cas.utils.conversion_utils import (
    add_parent_cell_hierarchy,
    add_parent_hierarchy_to_annotations,
    calculate_labelset,
    generate_cas_annotations,
    generate_parent_cell_lookup,
)


def anndata2cas(
    anndata_file_path: str,
    labelsets: List[str],
    output_file_path: str,
    include_hierarchy: bool,
):
    """
    Convert an AnnData file to Cell Annotation Schema (CAS) JSON.

    Args:
        anndata_file_path (str): Path to the AnnData file.
        labelsets (List[str]): List of labelsets, which are names of observation (obs) fields used to record author
        cell type names. The labelsets should be provided in order, starting from rank 0 (leaf nodes) and ascending
        to higher ranks.
        output_file_path (str): Output CAS file name.
        include_hierarchy (bool): Flag indicating whether to include hierarchy in the output.
    """

    anndata = read_anndata_file(anndata_file_path)

    labelset_dict = calculate_labelset(anndata.obs, labelsets)

    cas = generate_cas_metadata(dict(anndata.uns))

    generate_cas_labelsets(cas, labelset_dict)

    parent_cell_look_up = generate_parent_cell_lookup(anndata, labelset_dict)

    generate_cas_annotations(anndata, cas, labelset_dict, parent_cell_look_up)

    if include_hierarchy:
        add_parent_cell_hierarchy(parent_cell_look_up)
        add_parent_hierarchy_to_annotations(cas, parent_cell_look_up)

    # Write the JSON data to the file
    with open(output_file_path, "w") as json_file:
        json.dump(cas, json_file, indent=2)


def generate_cas_labelsets(
    cas: Dict[str, Any], labelset_dict: Dict[str, Dict[str, Any]]
):
    """
    Generates CAS labelsets and updates the provided CAS dictionary.

    Args:
        cas (Dict[str, Any]): The CAS dictionary.
        labelset_dict (Dict[str, Dict[str, Any]]): Dictionary containing labelset information (name, members, and rank).

    Returns:
        None
    """
    # labelsets
    for labelset, value in labelset_dict.items():
        cas.get("labelsets").append(
            {"name": labelset, "description": "", "rank": str(value.get("rank"))}
        )


def generate_cas_metadata(uns: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generates CAS metadata based on the provided 'uns' dictionary.

    Args:
        uns (Dict[str, Any]): The 'uns' dictionary containing metadata.

    Returns:
        Dict[str, Any]: The generated CAS metadata dictionary.
    """
    # TODO None values will be calculated later on
    matrix_file_id = None
    cellannotation_schema_version = uns["schema_version"]
    cellannotation_timestamp = None
    cellannotation_version = None
    cellannotation_url = None
    author_name = (
        "John Doe"  # Adding default author_name as it is required in the schema
    )
    author_contact = None
    orcid = None
    cas_init = {
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
    # Exclude keys with None values
    cas = {k: v for k, v in cas_init.items() if v is not None}
    return cas
