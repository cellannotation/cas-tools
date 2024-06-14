import csv
import os
from dataclasses import asdict

import pandas as pd

from cas.accession.hash_accession_manager import HashAccessionManager, is_hash_accession
from cas.accession.incremental_accession_manager import IncrementalAccessionManager


accession_manager = None


def serialize_to_tables(cta, file_name_prefix, out_folder, project_config):
    """
    Writes cell type annotation object to a series of tsv files.
    Tables to generate:
        - Annotation table (main)
        - Labelset table
        - Metadata
        - Annotation transfer

    Parameters:
        cta: cell type annotation object to serialize.
        file_name_prefix: Name prefix for table names
        out_folder: output folder path.
        project_config: project configuration with extra metadata
    """
    accession_prefix = project_config.get("accession_id_prefix")
    annotation_table_path = generate_annotation_table(accession_prefix, cta, out_folder)
    labelset_table_path = generate_labelset_table(cta, out_folder)
    metadata_table_path = generate_metadata_table(cta, project_config, out_folder)
    annotation_transfer_table_path = generate_annotation_transfer_table(cta, out_folder)
    reviews_table_path = generate_reviews_table(cta, out_folder)
    return [
        annotation_table_path,
        labelset_table_path,
        metadata_table_path,
        annotation_transfer_table_path,
        reviews_table_path,
    ]


def generate_annotation_transfer_table(cta, out_folder):
    """
    Generates annotation transfer table.

    Parameters:
        cta: cell type annotation object to serialize.
        out_folder: output folder path.
    """
    table_path = os.path.join(out_folder, "annotation_transfer.tsv")

    cta = asdict(cta)
    records = list()

    for annotation_object in cta["annotations"]:
        if (
            "cell_set_accession" in annotation_object
            and annotation_object["cell_set_accession"]
            and "transferred_annotations" in annotation_object
            and annotation_object["transferred_annotations"]
        ):
            for ta in annotation_object["transferred_annotations"]:
                record = dict()
                record["target_node_accession"] = annotation_object[
                    "cell_set_accession"
                ]
                record["transferred_cell_label"] = ta.get("transferred_cell_label", "")
                record["source_taxonomy"] = ta.get("source_taxonomy", "")
                record["source_node_accession"] = ta.get("source_node_accession", "")
                record["algorithm_name"] = ta.get("algorithm_name", "")
                record["comment"] = ta.get("comment", "")
                records.append(record)

    if not records:
        record = dict()
        record["target_node_accession"] = ""
        record["transferred_cell_label"] = ""
        record["source_taxonomy"] = ""
        record["source_node_accession"] = ""
        record["algorithm_name"] = ""
        record["comment"] = ""
        records.append(record)

    records_df = pd.DataFrame.from_records(records)
    records_df.to_csv(table_path, sep="\t", index=False)
    return table_path


def generate_metadata_table(cta, project_config, out_folder):
    """
    Generates the metadata table.

    Parameters:
        cta: cell type annotation object to serialize.
        project_config: metadata coming from project config
        out_folder: output folder path.
    """
    table_path = os.path.join(out_folder, "metadata.tsv")

    cta = asdict(cta)
    records = list()

    record = dict()
    record["index"] = "1"
    record["author_name"] = cta.get("author_name", "")
    record["author_contact"] = cta.get("author_contact", "")
    record["orcid"] = project_config.get("author", "")
    record["author_list"] = cta.get("author_list", "")
    record["matrix_file_id"] = project_config.get("matrix_file_id", "")
    record["title"] = cta.get("title", "")
    record["description"] = project_config.get("description", "")
    record["cellannotation_schema_version"] = cta.get(
        "cellannotation_schema_version", ""
    )
    record["cellannotation_timestamp"] = cta.get("cellannotation_timestamp", "")
    record["cellannotation_version"] = cta.get("cellannotation_version", "")
    record["cellannotation_url"] = cta.get("cellannotation_url", "")
    records.append(record)

    records_df = pd.DataFrame.from_records(records)
    records_df.to_csv(table_path, sep="\t", index=False)
    return table_path


def generate_labelset_table(cta, out_folder):
    """
    Generates labelset table.

    Parameters:
        cta: cell type annotation object to serialize.
        out_folder: output folder path.
    """
    table_path = os.path.join(out_folder + "labelset.tsv")

    cta = asdict(cta)
    records = list()

    for labelset in cta["labelsets"]:
        record = dict()
        record["name"] = labelset.get("name", "").replace("_name", "")
        record["description"] = labelset.get("description", "")
        record["rank"] = labelset.get("rank", "")
        record["annotation_method"] = labelset.get("annotation_method", "")
        if "automated_annotation" in labelset and labelset["automated_annotation"]:
            aut_annot = labelset["automated_annotation"]
            name_prefix = "automated_annotation_"
            record[name_prefix + "algorithm_name"] = aut_annot.get("algorithm_name", "")
            record[name_prefix + "algorithm_version"] = aut_annot.get(
                "algorithm_version", ""
            )
            record[name_prefix + "algorithm_repo_url"] = aut_annot.get(
                "algorithm_repo_url", ""
            )
            record[name_prefix + "reference_location"] = aut_annot.get(
                "reference_location", ""
            )
        else:
            name_prefix = "automated_annotation_"
            record[name_prefix + "algorithm_name"] = ""
            record[name_prefix + "algorithm_version"] = ""
            record[name_prefix + "algorithm_repo_url"] = ""
            record[name_prefix + "reference_location"] = ""
        records.append(record)

    records_df = pd.DataFrame.from_records(records)
    records_df.to_csv(table_path, sep="\t", index=False)
    return table_path


def generate_annotation_table(accession_prefix, cta, out_folder):
    """
    Generates annotation table.

    Parameters:
        cta: cell type annotation object to serialize.
        out_folder: output folder path.
        accession_prefix: accession id prefix
    """
    global accession_manager
    std_data_path = os.path.join(out_folder, "annotation.tsv")

    cta = asdict(cta)
    std_records = list()
    std_parent_records = list()
    std_parent_records_dict = dict()

    first_accession = cta["annotations"][0].get("cell_set_accession", "")
    if is_hash_accession(first_accession):
        accession_manager = HashAccessionManager()
    else:
        accession_manager = IncrementalAccessionManager(accession_prefix)
        # sort annotations by accession ids incrementing (if there is)
        if str(first_accession).split(":")[-1].split("_")[-1].isdigit():
            cta["annotations"].sort(
                key=lambda x: (
                    int(str(x["cell_set_accession"]).split(":")[-1].split("_")[-1])
                    if "cell_set_accession" in x
                    and x["cell_set_accession"]
                    and "_" in x["cell_set_accession"]
                    else 0
                )
            )

    id_index = dict()
    for annotation_object in cta["annotations"]:
        record = dict()
        if (
            "cell_set_accession" in annotation_object
            and annotation_object["cell_set_accession"]
        ):
            labelset = str(annotation_object.get("labelset", "")).replace("_name", "")
            record["cell_set_accession"] = accession_manager.generate_accession_id(
                id_recommendation=annotation_object.get("cell_set_accession", ""),
                labelset=labelset,
            )
            annotation_object["cell_set_accession"] = record["cell_set_accession"]
            record["cell_label"] = annotation_object.get("cell_label", "")
            record["cell_fullname"] = annotation_object.get("cell_fullname", "")
            record["parent_cell_set_accession"] = ""
            record["parent_cell_set_name"] = ""
            record["labelset"] = labelset
            record["cell_ontology_term_id"] = annotation_object.get(
                "cell_ontology_term_id", ""
            )
            record["cell_ontology_term"] = annotation_object.get(
                "cell_ontology_term", ""
            )
            record["rationale"] = annotation_object.get("rationale", "")
            record["rationale_dois"] = list_to_string(
                annotation_object.get("rationale_dois", [])
            )
            record["marker_gene_evidence"] = list_to_string(
                annotation_object.get("marker_gene_evidence", [])
            )
            record["synonyms"] = list_to_string(annotation_object.get("synonyms", []))
            if annotation_object.get("author_annotation_fields"):
                for key, value in annotation_object["author_annotation_fields"].items():
                    record[normalize_column_name(key)] = value
            # record["cell_ids"] = annotation_object.get("cell_ids", "")
            std_records.append(record)
            id_index[record["cell_set_accession"]] = record
        else:
            # parent nodes
            parent_label = annotation_object["cell_label"]
            if parent_label not in [
                parent["cell_label"] for parent in std_parent_records
            ]:
                record["cell_set_accession"] = ""
                record["cell_label"] = parent_label
                record["cell_fullname"] = ""
                record["parent_cell_set_accession"] = ""
                record["parent_cell_set_name"] = ""
                record["labelset"] = str(annotation_object.get("labelset", "")).replace(
                    "_name", ""
                )
                std_parent_records.append(record)
        if "parent_cell_set_name" in annotation_object:
            record["parent_cell_set_name"] = annotation_object["parent_cell_set_name"]
            if annotation_object["parent_cell_set_name"] in std_parent_records_dict:
                std_parent_records_dict.get(
                    annotation_object["parent_cell_set_name"]
                ).append(record)
            else:
                children = list()
                children.append(record)
                std_parent_records_dict[annotation_object["parent_cell_set_name"]] = (
                    children
                )
        if "parent_cell_set_accession" in annotation_object:
            record["parent_cell_set_accession"] = annotation_object[
                "parent_cell_set_accession"
            ]
    assign_parent_accession_ids(
        accession_manager, std_parent_records, std_parent_records_dict, cta["labelsets"]
    )
    assign_parent_cell_set_names(id_index)
    std_records.extend(std_parent_records)
    std_records_df = pd.DataFrame.from_records(std_records)
    std_records_df.to_csv(std_data_path, sep="\t", index=False)
    return std_data_path


def generate_reviews_table(cta, out_folder):
    """
    Generates annotation reviews table.

    Parameters:
        cta: cell type annotation object to serialize.
        out_folder: output folder path.
    """
    std_data_path = os.path.join(out_folder, "review.tsv")

    cta = asdict(cta)
    records = list()

    for annotation_object in cta["annotations"]:
        if (
            "cell_set_accession" in annotation_object
            and annotation_object["cell_set_accession"]
            and "reviews" in annotation_object
            and annotation_object["reviews"]
        ):
            labelset = str(annotation_object.get("labelset", "")).replace("_name", "")
            if accession_manager:
                accession = accession_manager.generate_accession_id(
                    id_recommendation=annotation_object.get("cell_set_accession", ""),
                    labelset=labelset,
                )
            else:
                accession = annotation_object.get("cell_set_accession")
            for review in annotation_object["reviews"]:
                record = dict()
                record["target_node_accession"] = accession
                record["datestamp"] = review.get("datestamp", "")
                if record["datestamp"]:
                    # convert time to ISO 8601 format
                    record["datestamp"] = (
                        record["datestamp"].strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z"
                    )
                record["reviewer"] = review.get("reviewer", "")
                record["review"] = review.get("review", "")
                record["explanation"] = review.get("explanation", "")
                records.append(record)

    if records:
        std_records_df = pd.DataFrame.from_records(records)
        std_records_df.to_csv(std_data_path, sep="\t", index=False)
    else:
        row = ["target_node_accession", "datestamp", "reviewer", "review", "explanation"]
        with open(std_data_path, "w") as f_output:
            tsv_output = csv.writer(f_output, delimiter="\t")
            tsv_output.writerow(row)

    return std_data_path


def list_to_string(my_list: list):
    """
    Converts a list to its string representation. Nanobot has problem with single quotations so removes them as well.
    Parameters:
        my_list: list to serialize
    Returns:
        string representation of the list
    """
    if not my_list:
        str_value = ""
    else:
        str_value = "|".join(my_list)
    return str_value


def assign_parent_accession_ids(
    accession_manager, std_parent_records, std_parent_records_dict, labelsets
):
    """
    Assigns accession ids to parent clusters and updates their references from the child clusters.
    Parameters:
        accession_manager: accession ID generator
        std_parent_records: list of all parents to assign accession ids
        std_parent_records_dict: parent cluster - child clusters dictionary
        labelsets: labelsets list
    """
    label_set_ranks = dict(
        [
            (label_set["name"].replace("_name", ""), label_set["rank"])
            for label_set in labelsets
        ]
    )

    std_parent_records.sort(key=lambda x: int(label_set_ranks[x["labelset"]]))
    for std_parent_record in std_parent_records:
        accession_id = accession_manager.generate_accession_id()
        std_parent_record["cell_set_accession"] = accession_id

        children = std_parent_records_dict.get(std_parent_record["cell_label"], list())
        for child in children:
            child["parent_cell_set_accession"] = accession_id


def assign_parent_cell_set_names(id_index: dict):
    """
    Assigns parent cell set names to the child cell sets.
    Parameters:
        id_index: dictionary of cell set accessions and their corresponding records
    """
    for key, value in id_index.items():
        if value["parent_cell_set_accession"]:
            parent_record = id_index.get(value["parent_cell_set_accession"])
            if parent_record:
                value["parent_cell_set_name"] = parent_record["cell_label"]


def normalize_column_name(column_name: str) -> str:
    """
    Normalizes column name for url compatibility.
    URL compatible column name requirement: All names must match: ^[\w_ ]+$' for to_url()

    Parameters:
        column_name: current column name
    Returns:
        normalized column_name
    """
    return column_name.strip().replace("(", "_").replace(")", "_").replace("-", "_")
