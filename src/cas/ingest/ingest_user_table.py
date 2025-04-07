import os
from importlib.metadata import version
from pathlib import Path
from typing import get_type_hints

from cas.file_utils import read_config, read_table_to_dict, write_json_file
from cas.flatten_data_to_tables import serialize_to_tables
from cas.ingest.config_validator import validate
from cas.model import (
    Annotation,
    AnnotationTransfer,
    AutomatedAnnotation,
    CellTypeAnnotation,
    Labelset,
)
from cas.accession.hash_accession_manager import HashAccessionManager, is_hash_accession
from cas.accession.incremental_accession_manager import IncrementalAccessionManager

NAME_SEPERATOR = "_XX_"


def ingest_data(
    data_file: str,
    config_file: str,
    out_file: str,
    format: str = "json",
    print_undefined: bool = False,
) -> dict:
    """
    Ingests given data into standard cell annotation schema data structure using the given configuration.

    :param data_file: Unformatted user data in tsv/csv format.
    :param config_file: configuration file path.
    :param out_file: output file path.
    :param format: Data export format. Supported formats are 'json' and 'tsv'
    :param print_undefined: prints null values to the output json if true. Omits undefined values from the json output if
    false. False by default. Only effective in json serialization.
    :return: output data as dict
    """
    cas = ingest_user_data(data_file, config_file)
    project_config = read_config(config_file)
    if format == "json":
        assign_accession_ids(cas, project_config.get("taxonomy_id", ""))
        write_json_file(cas, out_file, print_undefined)
    elif format == "tsv":
        table_name_prefix = os.path.splitext(os.path.basename(data_file))[0]
        if os.path.isfile(out_file):
            out_folder = Path(out_file).parent.absolute()
        else:
            out_folder = out_file
        serialize_to_tables(cas, table_name_prefix, out_folder, project_config)

    return cas.to_dict()


def ingest_user_data(data_file: str, config_file: str):
    """
    Ingest given user data into standard cell annotation schema data structure using the given configuration.
    :param data_file: Unformatted user data in tsv/csv format.
    :param config_file: configuration file path.
    """

    config = read_config(config_file)
    is_config_valid = validate(config)
    if not is_config_valid:
        raise Exception("Configuration file is not valid!")
    cas = CellTypeAnnotation(config["author_name"], list(), config["title"])
    cas.description = config.get("description", "")
    cas.cellannotation_schema_version = version("cell-annotation-schema")
    headers, records = read_table_to_dict(data_file, generated_ids=True)
    config_fields = config["fields"]
    labelset_ranks = populate_labelsets(cas, config_fields)
    ao_names = dict()
    utilized_columns = set()
    for record_index in records:
        record = records[record_index]
        if not all(value == "" for value in record.values()):  # skip empty rows
            ao = Annotation("", "")
            parents = [None] * 10
            for field in config_fields:
                # handle hierarchical columns
                if field["column_type"] == "cluster_name":
                    ao.labelset = field["column_name"]
                    ao.cell_label = str(record[field["column_name"]])
                    utilized_columns.add(field["column_name"])
                elif field["column_type"] == "cluster_id":
                    ao.cell_set_accession = str(record[field["column_name"]])
                    ao.rank = int(str(field["rank"]).strip())
                    utilized_columns.add(field["column_name"])
                elif field["column_type"] == "cell_set":
                    if record[field["column_name"]]:
                        parent_ao = get_annotation(ao_names, field, record)
                        register_parent(field, labelset_ranks, parent_ao, parents)
                    utilized_columns.add(field["column_name"])
                else:
                    # handle annotation columns
                    if "typing.List[str]" in str(
                        get_type_hints(ao)[field["column_type"]]
                    ):
                        list_value = str(record[field["column_name"]]).split(",")
                        stripped = list(map(str.strip, list_value))
                        setattr(ao, field["column_type"], stripped)
                    else:
                        setattr(ao, field["column_type"], record[field["column_name"]])
                    utilized_columns.add(field["column_name"])

        add_user_annotations(ao, headers, record, utilized_columns)
        add_parent_node_names(ao, ao_names, cas, parents)

        ao_names[ao.labelset + NAME_SEPERATOR + ao.cell_label] = ao
        cas.add_annotation_object(ao)
    return cas


def register_parent(field, labelset_ranks, parent_ao, parents):
    """
    Registers the parent annotation object to the parents list.
    Args:
        field: config field
        labelset_ranks: labelset ranks dictionary
        parent_ao: parent to add
        parents: sparse parents list
    """
    parents.insert(int(str(field["rank"]).strip()), parent_ao)


def get_annotation(ao_names, field, record):
    """
    Creates a annotation object if it does not exist in the ao_names dictionary at the same labelset.
    Args:
        ao_names: list of existing annotation objects
        field: config field
        record: data record

    Returns: annotation object
    """
    # labelset_XX_label
    name = field["column_name"] + NAME_SEPERATOR + record[field["column_name"]]
    if name in ao_names:
        ao = ao_names[name]
    else:
        ao = Annotation(field["column_name"], record[field["column_name"]])
    ao.labelset = field["column_name"]
    return ao


def add_user_annotations(ao, headers, record, utilized_columns):
    """
    Adds user annotations that are not supported by the standard schema.
    :param ao: current annotation object
    :param headers: all column names of the user data
    :param record: a record in the user data
    :param utilized_columns: list of processed columns
    """
    not_utilized_columns = [
        column_name for column_name in headers if column_name not in utilized_columns
    ]
    for not_utilized_column in not_utilized_columns:
        if record[not_utilized_column]:
            ao.add_user_annotation(not_utilized_column, record[not_utilized_column])


def add_parent_node_names(ao, ao_names, cas, parents):
    """
    Creates parent nodes if necessary and creates a cluster hierarchy through assigning parent_node_names.
    :param ao: current annotation object
    :param ao_names: list of all created annotation objects
    :param cas: main object
    :param parents: list of current annotation object's parents
    """
    if parents:
        # get first non-null parent
        direct_parent = [x for x in parents if x][0]
        ao.parent_cell_set_name = direct_parent.cell_label
        prev = None
        for parent in reversed(parents):
            if parent:
                if prev:
                    if (
                        parent.parent_cell_set_name
                        and parent.parent_cell_set_name != prev.cell_label
                        and parent.cell_label != prev.cell_label
                    ):
                        print(
                            "Annotation {} has multiple parents: {} and {}".format(
                                parent.cell_label,
                                parent.parent_cell_set_name,
                                prev.cell_label,
                            )
                        )
                    if parent.labelset + parent.cell_label != prev.labelset + prev.cell_label:  # avoid self-references
                        parent.parent_cell_set_name = prev.cell_label
                prev = parent
                if parent.labelset + NAME_SEPERATOR + parent.cell_label not in ao_names:
                    cas.add_annotation_object(parent)
                    ao_names[parent.labelset + NAME_SEPERATOR + parent.cell_label] = parent


def populate_labelsets(cas, config_fields):
    """
    Populates labelsets list based on the fields of the config.
    :param cas: main object
    :param config_fields: config file fields
    :return: ranks of the labelsets
    """
    labelsets = list()
    ranks = dict()
    for field in config_fields:
        if field["column_type"] == "cell_set" or field["column_type"] == "cluster_name":
            label_set = Labelset(field["column_name"])
            if "rank" in field:
                label_set.rank = str(field["rank"])
                ranks[field["column_name"]] = int(field["rank"])
            labelsets.append(label_set)
    if labelsets:
        cas.labelsets = labelsets
    return ranks

def assign_accession_ids(cas, accession_prefix):
    """
    Assigns accession IDs to the annotations in the given CellTypeAnnotation object.
    :param cas: CellTypeAnnotation object
    :param accession_prefix: Prefix for the accession IDs.
    """
    parent_records = list()
    parent_records_dict = dict()

    if not accession_prefix.endswith("_"):
        accession_prefix += "_"

    first_accession = cas.annotations[0].cell_set_accession
    if is_hash_accession(first_accession):
        accession_manager = HashAccessionManager()
    else:
        accession_manager = IncrementalAccessionManager(accession_prefix)
        # sort annotations by accession ids incrementing (if there is)
        if str(first_accession).split(":")[-1].split("_")[-1].isdigit():
            cas.annotations.sort(
                key=lambda x: (
                    int(str(x.cell_set_accession).split(":")[-1].split("_")[-1])
                    if x.cell_set_accession
                       and "_" in x.cell_set_accession
                    else 0
                )
            )

    id_index = dict()
    for annotation_object in cas.annotations:
        if annotation_object.cell_set_accession:
            if not annotation_object.cell_set_accession.startswith(accession_prefix):
                annotation_object.cell_set_accession = (
                    accession_prefix + '_' + annotation_object.cell_set_accession
                )
                labelset = annotation_object.labelset
                accession = annotation_object.cell_set_accession or ""
                annotation_object.cell_set_accession = accession_manager.generate_accession_id(
                    id_recommendation=accession,
                    labelset=labelset,
                )
        else:
            # parent nodes
            parent_records.append(annotation_object)
        if annotation_object.parent_cell_set_name in parent_records_dict:
            parent_records_dict.get(
                annotation_object.parent_cell_set_name
            ).append(annotation_object)
        else:
            children = list()
            children.append(annotation_object)
            parent_records_dict[annotation_object.parent_cell_set_name] = children

    assign_parent_accession_ids(accession_manager, parent_records, parent_records_dict, cas.labelsets)


def assign_parent_accession_ids(
    accession_manager, std_parent_records, std_parent_records_dict, labelsets
):
    """
    Assigns accession ids to parent nodes and updates their references from the childs.
    Args:
        accession_manager: accession ID generator
        std_parent_records: list of all parents to assign accession ids
        std_parent_records_dict: parent cluster - child clusters dictionary
        labelsets: labelsets list
    """
    label_set_ranks = dict(
        [
            (label_set.name.replace("_name", ""), label_set.rank)
            for label_set in labelsets
        ]
    )
    std_parent_records.sort(key=lambda x: int(label_set_ranks[x.labelset]))
    for std_parent_record in std_parent_records:
        accession_id = accession_manager.generate_accession_id()
        std_parent_record.cell_set_accession = accession_id

        children = std_parent_records_dict.get(std_parent_record.cell_label, list())
        for child in children:
            if not (child.cell_label == std_parent_record.cell_label and child.labelset ==
                    std_parent_record.labelset):  # prevent self reference
                if not child.parent_cell_set_accession:  # prevent overwriting existing parent (same parent name in different labelsets)
                    if int(label_set_ranks[child.labelset]) < int(label_set_ranks[
                                                                         std_parent_record.labelset]):  # prevent parent assignment to a child in a higher rank
                        child.parent_cell_set_accession = accession_id
