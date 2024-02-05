import json
from typing import List

from cas.file_utils import read_anndata_file
from cas.spreadsheet_to_cas import calculate_labelset_rank, get_cell_ids
from cas.accession.hash_accession_manager import HashAccessionManager


def anndata2cas(anndata_file_path: str, labelsets: List[str], output_file_path: str, include_hierarchy: bool):
    """
    Convert an AnnData file to Cell Annotation Schema (CAS) JSON.

    Args:
        anndata_file_path (str): Path to the AnnData file.
        labelsets (List[str]): List of labelsets.
        output_file_path (str): Output CAS file name.
        include_hierarchy (bool): Flag indicating whether to include hierarchy in the output.
    """

    anndata = read_anndata_file(anndata_file_path)

    labelset_dict = {}
    for item in labelsets:
        labelset_dict.update({item: set(anndata.obs[item])})

    # metadata
    cas = {
        "matrix_file_id": "TBA",
        "cellannotation_schema_version": "TBA",
        "cellannotation_timestamp": "TBA",
        "cellannotation_version": "TBA",
        "cellannotation_url": "TBA",
        "author_name": "TBA",
        "author_contact": "TBA",
        "orcid": "TBA",
        "annotations": [],
        "labelsets": [],
    }

    # annotation
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
            rationale = "TBA"
            rationale_dois = "TBA"
            marker_gene_evidence = "TBA"
            synonyms = "TBA"
            category_fullname = "TBA"

            if include_hierarchy:
                if cell_label in parent_cell_look_up:
                    parent_cell_look_up[cell_label].get("cell_ids").update(cell_ids)
                else:
                    parent_cell_look_up[cell_label] = {
                        "cell_ids": set(cell_ids),
                        "accession": cell_set_accession,
                    }

            anno = {
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
                "category_cell_ontology_exists": "TBA",
                "category_cell_ontology_term_id": "TBA",
                "category_cell_ontology_term": "TBA",
            }
            cas.get("annotations").append(anno)

    # labelsets
    labelset_rank_dict = calculate_labelset_rank(labelsets)
    for labelset, rank in labelset_rank_dict.items():
        cas.get("labelsets").append(
            {"name": labelset, "description": "TBA", "rank": str(rank)}
        )

    # add parent_cell_set
    if include_hierarchy:
        for key, value in parent_cell_look_up.items():
            for inner_key, inner_value in parent_cell_look_up.items():
                if value != inner_value and value.get("cell_ids").issubset(
                    inner_value.get("cell_ids")
                ):
                    value.update(
                        {"parent": inner_key, "p_accession": inner_value.get("accession")}
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

    # Write the JSON data to the file
    with open(output_file_path, "w") as json_file:
        json.dump(cas, json_file, indent=2)
