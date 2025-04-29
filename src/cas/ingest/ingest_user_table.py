import os
import re
from importlib.metadata import version
from pathlib import Path
from typing import get_type_hints

from cas.file_utils import read_config, read_table_to_dict, write_json_file
from cas.flatten_data_to_tables import serialize_to_tables
from cas.ingest.config_validator import validate
from cas.accession.incremental_accession_manager import IncrementalAccessionManager
from cas.model import (
    Annotation,
    AnnotationTransfer,
    AutomatedAnnotation,
    CellTypeAnnotation,
    Labelset,
)

NAME_SEPERATOR = "_XX_"


def ingest_data(
    data_file: str,
    config_file: str,
    out_file: str,
    format: str = "json",
    print_undefined: bool = False,
    generate_accession_ids: bool = False
) -> dict:
    """
    Ingests given data into standard cell annotation schema data structure using the given configuration.

    :param data_file: Unformatted user data in tsv/csv format.
    :param config_file: configuration file path.
    :param out_file: output file path.
    :param format: Data export format. Supported formats are 'json' and 'tsv'
    :param print_undefined: prints null values to the output json if true. Omits undefined values from the json output if
    false. False by default. Only effective in json serialization.
    :param generate_accession_ids: determines if incrementally generate accession_ids for all annotations that don't have an id.
    :return: output data as dict
    """
    cas = ingest_user_data(data_file, config_file, generate_accession_ids)

    if format == "json":
        write_json_file(cas, out_file, print_undefined)
    elif format == "tsv":
        table_name_prefix = os.path.splitext(os.path.basename(data_file))[0]
        if os.path.isfile(out_file):
            out_folder = Path(out_file).parent.absolute()
        else:
            out_folder = out_file
        serialize_to_tables(cas, table_name_prefix, out_folder)

    return cas.to_dict()


def ingest_user_data(data_file: str, config_file: str, generate_accession_ids: bool = False) -> CellTypeAnnotation:
    """
    Ingest given user data into standard cell annotation schema data structure using the given configuration.
    :param data_file: Unformatted user data in tsv/csv format.
    :param config_file: configuration file path.
    :param generate_accession_ids: determines if incrementally generate accession_ids for all annotations that don't have an id.
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
    cluster_accession_prefix = config.get("accession_prefix", "").strip()
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
                    cell_set_accession = str(record[field["column_name"]])
                    if field.get("accession_prefix"):
                        cluster_accession_prefix = field.get("accession_prefix")
                    if cluster_accession_prefix and not cell_set_accession.startswith(cluster_accession_prefix):
                        cell_set_accession = cluster_accession_prefix + "_" + cell_set_accession
                    ao.cell_set_accession = cell_set_accession
                    ao.rank = int(str(field["rank"]).strip())
                    utilized_columns.add(field["column_name"])
                elif field["column_type"] == "cell_set":
                    if record[field["column_name"]]:
                        parent_ao = get_annotation(ao_names, field, record)
                        if field.get("accession_column"):
                            accession = str(record[field["accession_column"]]).strip()
                            if not accession:
                                raise ValueError("Accession is empty for {0}({1})".format(parent_ao.cell_label, parent_ao.labelset))
                            parent_ao.cell_set_accession = accession
                        register_parent(field, labelset_ranks, parent_ao, parents)
                    utilized_columns.add(field["column_name"])
                else:
                    # handle annotation columns
                    if "typing.List[str]" in str(
                        get_type_hints(ao)[field["column_type"]]
                    ):
                        list_value = re.split(r'[,|]', str(record[field["column_name"]]))
                        stripped = [s.strip() for s in map(str.strip, list_value) if s]
                        setattr(ao, field["column_type"], stripped)
                    else:
                        setattr(ao, field["column_type"], record[field["column_name"]])
                    utilized_columns.add(field["column_name"])

            add_user_annotations(ao, headers, record, utilized_columns)
            add_parent_node_names(ao, ao_names, cas, parents)

            ao_names[ao.labelset + NAME_SEPERATOR + ao.cell_label] = ao
            cas.add_annotation_object(ao)

    if generate_accession_ids:
        cas = generate_ids_for_annotations(cas, config, labelset_ranks)
    return cas


def generate_ids_for_annotations(cas: CellTypeAnnotation, config: dict, labelset_ranks: dict) -> CellTypeAnnotation:
    """
    Generates unique IDs for the annotations in the given CellTypeAnnotation object.
    :param cas: CellTypeAnnotation object
    :param config: ingestion configuration dictionary
    :param labelset_ranks: ranks of the labelsets
    :return: CellTypeAnnotation object with generated IDs.
    """
    accession_managers = init_accession_managers(cas, config)

    label_to_accession = dict()
    for annotation in cas.annotations:
        if not annotation.cell_set_accession:
            accession_manager = accession_managers.get(annotation.labelset)
            annotation.cell_set_accession = accession_manager.generate_accession_id()
        if annotation.cell_label not in label_to_accession:
            label_to_accession[annotation.cell_label] = []
        label_to_accession[annotation.cell_label].append(annotation)

    for annotation in cas.annotations:
        if annotation.parent_cell_set_name:
            parent_candidates = label_to_accession.get(annotation.parent_cell_set_name, None)
            parent_candidates_sorted = sorted(
                parent_candidates, key=lambda x: labelset_ranks.get(x.labelset, float('inf'))
            )
            # Assign the first parent with a rank greater than the current annotation's rank
            for parent_annotation in parent_candidates_sorted:
                if (
                        labelset_ranks.get(parent_annotation.labelset, float('inf')) >
                        labelset_ranks.get(annotation.labelset, float('inf'))
                ):
                    annotation.parent_cell_set_accession = parent_annotation.cell_set_accession
                    break

    return cas


def init_accession_managers(cas: CellTypeAnnotation, config: dict) -> dict:
    """
    Initializes IncrementalAccessionManager for each labelset in the config.
    :param cas: CellTypeAnnotation object
    :param config: ingestion configuration dictionary
    :return: dictionary of IncrementalAccessionManager objects
    """
    example_accession = None
    max_value = -1

    for annotation in cas.annotations:
        accession = annotation.cell_set_accession
        if accession:
            example_accession = accession
            parts = re.split(r'[_:]', accession)
            last_part = parts[-1]
            if last_part.isdigit():
                value = int(last_part)
                if value > max_value:
                    max_value = value

    if "_" not in example_accession:
        accession_prefix = config.get("accession_prefix", "").strip()
        if not accession_prefix.endswith("_"):
            accession_prefix += "_"
    else:
        if "_" in example_accession:
            accession_prefix = example_accession.split("_")[0] + "_"
        else:
            accession_prefix = ""
    default_accession_manager = IncrementalAccessionManager(accession_prefix, max_value)

    labelset_accession_managers = dict()
    for field in config["fields"]:
        if field["column_type"] in {"cluster_name", "cell_set"}:
            if field.get("accession_prefix"):
                prefix = field.get("accession_prefix", accession_prefix)
                labelset_accession_managers[field["column_name"]] = IncrementalAccessionManager(
                    prefix, int(field.get("accession_start", 0)),)
            else:
                labelset_accession_managers[field["column_name"]] = default_accession_manager
    return labelset_accession_managers



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
                label_set.rank = int(field["rank"])
                ranks[field["column_name"]] = int(field["rank"])
            labelsets.append(label_set)
    if labelsets:
        cas.labelsets = labelsets
    return ranks
