import json
from typing import List

from cas.file_utils import read_anndata_file
from cas.spreadsheet_to_cas import calculate_labelset_rank, get_cell_ids
from cas.accession.hash_accession_manager import HashAccessionManager


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

    parent_cell_look_up = generate_cas_annotations(
        anndata, cas, include_hierarchy, labelset_dict
    )

    generate_cas_labelsets(cas, labelsets)

    add_parent_cell_hierarchy(cas, include_hierarchy, parent_cell_look_up)

    # Write the JSON data to the file
    with open(output_file_path, "w") as json_file:
        json.dump(cas, json_file, indent=2)


def add_parent_cell_hierarchy(cas, include_hierarchy, parent_cell_look_up):
    """
    Adds parent cell hierarchy information to the CAS dictionary.

    Args:
        cas (dict): The CAS dictionary.
        include_hierarchy (bool): Flag indicating whether to include hierarchy in the output.
        parent_cell_look_up (dict): Dictionary containing parent cell information.

    Returns:
        None
    """
    if include_hierarchy:
        for key, value in parent_cell_look_up.items():
            for inner_key, inner_value in parent_cell_look_up.items():
                if value != inner_value and value.get("cell_ids").issubset(
                    inner_value.get("cell_ids")
                ):
                    value.update(
                        {
                            "parent": inner_key,
                            "p_accession": inner_value.get("accession"),
                        }
                    )

        annotation_list = cas.get("annotations")
        for annotation in annotation_list:
            annotation.update(
                {
                    "parent_cell_set_name": parent_cell_look_up.get(
                        annotation.get("cell_label")
                    ).get("parent")
                }
            )
            annotation.update(
                {
                    "parent_cell_set_accession": parent_cell_look_up.get(
                        annotation.get("cell_label")
                    ).get("p_accession")
                }
            )


def generate_cas_annotations(anndata, cas, include_hierarchy, labelset_dict):
    """
    Generates CAS annotations and returns parent_cell_look_up dictionary that is calculated during annotation
    generation.

    Args:
        anndata (AnnData): The AnnData object.
        cas (dict): The CAS dictionary.
        include_hierarchy (bool): Flag indicating whether to include hierarchy in the output.
        labelset_dict (dict): Dictionary containing labelsets.

    Returns:
        parent_cell_look_up (dict): A dictionary containing information about parent cell lookups.
    """
    accession_manager = HashAccessionManager()
    parent_cell_look_up = {}
    for k, v in labelset_dict.items():
        for label in v:
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


def generate_cas_labelsets(cas, labelsets):
    """
    Generates CAS labelsets and updates the provided CAS dictionary.

    Args:
        cas (dict): The CAS dictionary.
        labelsets (list): List of labelsets.

    Returns:
        None
    """
    # labelsets
    labelset_rank_dict = calculate_labelset_rank(labelsets)
    for labelset, rank in labelset_rank_dict.items():
        cas.get("labelsets").append(
            {"name": labelset, "description": "", "rank": str(rank)}
        )


def generate_cas_metadata(uns):
    """
    Generates CAS metadata based on the provided 'uns' dictionary.

    Args:
        uns (dict): The 'uns' dictionary containing metadata.

    Returns:
        dict: The generated CAS metadata dictionary.
    """
    # TODO None values will be calculated later on
    matrix_file_id = None
    cellannotation_schema_version = uns["schema_version"]
    cellannotation_timestamp = None
    cellannotation_version = None
    cellannotation_url = None
    author_name = "John Doe"  # Adding default author_name as it is required in the schema
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


def calculate_labelset(obs, labelsets):
    """
    Calculates labelset dictionary based on the provided observations.

    Args:
        obs (DataFrame): DataFrame containing observations.
        labelsets (List[str]): List of labelsets.

    Returns:
        dict: The calculated labelset dictionary.
    """
    labelset_dict = {}
    for item in labelsets:
        labelset_dict.update({item: set(obs[item])})
    return labelset_dict
