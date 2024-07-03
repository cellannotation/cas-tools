import json
from typing import Any, Dict, List

from cas.file_utils import read_anndata_file
from cas.utils.conversion_utils import (
    add_labelsets_to_cas,
    add_parent_cell_hierarchy,
    add_parent_hierarchy_to_annotations,
    calculate_labelset,
    get_authors_from_doi,
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

    add_labelsets_to_cas(cas, labelset_dict)

    parent_cell_look_up = generate_parent_cell_lookup(anndata, labelset_dict)

    add_annotations_to_cas(cas, labelset_dict, parent_cell_look_up)

    if include_hierarchy:
        add_parent_cell_hierarchy(parent_cell_look_up)
        add_parent_hierarchy_to_annotations(cas, parent_cell_look_up)

    # Write the JSON data to the file
    with open(output_file_path, "w") as json_file:
        json.dump(cas, json_file, indent=2)


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
    author_list = get_authors_from_doi(uns["citation"].split(" ")[1]) if "citation" in uns else None
    title = uns.get("title")
    cas_init = {
        "matrix_file_id": matrix_file_id,
        "cellannotation_schema_version": cellannotation_schema_version,
        "cellannotation_timestamp": cellannotation_timestamp,
        "cellannotation_version": cellannotation_version,
        "cellannotation_url": cellannotation_url,
        "author_list": author_list,
        "title": title,
        "annotations": [],
        "labelsets": [],
    }
    # Exclude keys with None values
    cas = {k: v for k, v in cas_init.items() if v is not None}
    return cas


def add_annotations_to_cas(
    cas: Dict[str, Any],
    labelset_dict: Dict[str, Any],
    parent_cell_look_up: Dict[str, Any],
):
    """
    Generates CAS annotations based on the provided AnnData object and updates the CAS
    dictionary with new annotations. This function can optionally use a precomputed
    parent cell lookup dictionary to enrich the annotations with hierarchical information.

    Args:
        cas (Dict[str, Any]): The CAS dictionary to be updated with annotations. Expected to have a key
            'annotations' where new annotations will be appended.
        labelset_dict (Dict[str, Any]): A dictionary defining labelsets and their members. This is used to match cell
            labels with their respective metadata and annotations.
        parent_cell_look_up (Dict[str, Any]): A precomputed dictionary containing hierarchical metadata about cell
            labels.

    Returns:
        None: The function directly updates the `cas` dictionary with new annotations. The `parent_cell_look_up` is
        used for enrichment and must be generated beforehand if hierarchical information is to be included.
    """
    for k, v in labelset_dict.items():
        for label in v["members"]:
            labelset = k
            rationale = None
            rationale_dois = None
            marker_gene_evidence = None
            synonyms = None
            category_fullname = None
            category_cell_ontology_exists = None
            category_cell_ontology_term_id = None
            category_cell_ontology_term = None

            anno_init = {
                "labelset": labelset,
                "cell_label": label,
                "cell_fullname": label,
                "cell_set_accession": parent_cell_look_up[f"{labelset}:{label}"]["accession"],
                "cell_ontology_term_id": parent_cell_look_up[f"{labelset}:{label}"][
                    "cell_ontology_term_id"
                ],
                "cell_ontology_term": parent_cell_look_up[f"{labelset}:{label}"]["cell_ontology_term"],
                "cell_ids": list(parent_cell_look_up[f"{labelset}:{label}"]["cell_ids"]),
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
