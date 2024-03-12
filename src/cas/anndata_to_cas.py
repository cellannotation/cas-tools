import itertools
import json
from typing import Any, Dict, List

import anndata as ad
import pandas as pd

from cas.accession.hash_accession_manager import HashAccessionManager
from cas.file_utils import read_anndata_file
from cas.spreadsheet_to_cas import calculate_labelset_rank, get_cell_ids


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
        labelsets (List[str]): List of labelsets.
        output_file_path (str): Output CAS file name.
        include_hierarchy (bool): Flag indicating whether to include hierarchy in the output.
    """

    anndata = read_anndata_file(anndata_file_path)

    labelset_dict = calculate_labelset(anndata.obs, labelsets)

    cas = generate_cas_metadata(dict(anndata.uns))

    generate_cas_labelsets(cas, labelset_dict)

    parent_cell_look_up = generate_cas_annotations(
        anndata, cas, include_hierarchy, labelset_dict
    )

    if include_hierarchy:
        add_parent_cell_hierarchy(cas, parent_cell_look_up)

    # Write the JSON data to the file
    with open(output_file_path, "w") as json_file:
        json.dump(cas, json_file, indent=2)


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
            "parent": parent_key,
            "p_accession": parent_value.get("accession"),
            "parent_rank": parent_value.get("rank"),
        }
    )


def add_parent_cell_hierarchy(cas: Dict[str, Any], parent_cell_look_up: Dict[str, Any]):
    """
    Adds parent cell hierarchy information to the CAS dictionary.

    Args:
        cas (Dict[str, Any]): The CAS dictionary.
        parent_cell_look_up (Dict[str, Any]): Dictionary containing parent cell information.

    Returns:
        None
    """
    for (key, value), (inner_key, inner_value) in itertools.product(
        parent_cell_look_up.items(), repeat=2
    ):
        if key == inner_key:
            continue
        if value.get("parent") and value.get("parent_rank", 0) <= inner_value.get(
            "rank"
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
                    f"{key} and {inner_key} cell labels have the same cell_ids. cell_ids can't be identical at "
                    f"the same rank."
                )

    annotation_list = cas.get("annotations")
    for annotation in annotation_list:
        parent_info = parent_cell_look_up.get(annotation.get("cell_label"))
        parent = parent_info.get("parent")
        p_accession = parent_info.get("p_accession")
        if parent and p_accession:
            # Add parent data
            annotation.update({"parent_cell_set_name": parent})
            annotation.update({"parent_cell_set_accession": p_accession})
            # Remove CL annotation if parent and child has the same annotation
            parent_cell_ontology_term_id = parent_info.get("cell_ontology_term_id")
            parent_cell_ontology_term = parent_info.get("cell_ontology_term")
            if parent_cell_ontology_term_id == annotation.get(
                "cell_ontology_term_id"
            ) and parent_cell_ontology_term == annotation.get("cell_ontology_term"):
                annotation.pop("cell_ontology_term_id", None)
                annotation.pop("cell_ontology_term", None)


def generate_cas_annotations(
    anndata: ad.AnnData,
    cas: Dict[str, Any],
    include_hierarchy: bool,
    labelset_dict: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Generates CAS annotations and returns parent_cell_look_up dictionary that is calculated during annotation
    generation.

    Args:
        anndata (ad.AnnData): The AnnData object.
        cas (Dict[str, Any]): The CAS dictionary.
        include_hierarchy (bool): Flag indicating whether to include hierarchy in the output.
        labelset_dict (Dict[str, Any]): Dictionary containing labelsets.

    Returns:
        parent_cell_look_up (Dict[str, Any]): A dictionary containing information about parent cell lookups.
    """
    accession_manager = HashAccessionManager()
    parent_cell_look_up = {}
    for k, v in labelset_dict.items():
        for label in v["members"]:
            labelset = k
            cell_label = label
            cell_ontology_term_id = anndata.obs[anndata.obs[k] == cell_label][
                "cell_type_ontology_term_id"
            ].iloc[0]
            cell_ontology_term = anndata.obs[anndata.obs[k] == cell_label][
                "cell_type"
            ].iloc[0]
            cell_ids = get_cell_ids(anndata, labelset, label)
            cell_set_accession = accession_manager.generate_accession_id(
                cell_ids=cell_ids
            )
            # TODO None values will be calculated later on
            rationale = None
            rationale_dois = None
            marker_gene_evidence = None
            synonyms = None
            category_fullname = None
            category_cell_ontology_exists = None
            category_cell_ontology_term_id = None
            category_cell_ontology_term = None

            if include_hierarchy:
                if cell_label in parent_cell_look_up:
                    parent_cell_look_up[cell_label].get("cell_ids").update(cell_ids)
                else:
                    parent_cell_look_up[cell_label] = {
                        "cell_ids": set(cell_ids),
                        "accession": cell_set_accession,
                        "rank": v.get("rank"),
                        "cell_ontology_term_id": cell_ontology_term_id,
                        "cell_ontology_term": cell_ontology_term,
                    }

            anno_init = {
                "labelset": labelset,
                "cell_label": cell_label,
                "cell_fullname": cell_label,
                "cell_set_accession": cell_set_accession,
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
            }
            # Exclude keys with None values
            cas.get("annotations").append(
                {k: v for k, v in anno_init.items() if v is not None}
            )
    return parent_cell_look_up


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
    labelset_rank_dict = calculate_labelset_rank(labelsets)
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
