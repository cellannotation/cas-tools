import os
import pkg_resources

from pathlib import Path
from typing import get_type_hints

from cas.file_utils import read_config, read_tsv_to_dict, write_json_file
from cas.flatten_data_to_tables import serialize_to_tables
from cas.ingest.config_validator import validate
from cas.model import (
    Annotation,
    AnnotationTransfer,
    AutomatedAnnotation,
    CellTypeAnnotation,
    Labelset,
    UserAnnotation,
)


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
    cas = CellTypeAnnotation(config["author_name"], list())
    cas.cellannotation_schema_version = pkg_resources.get_distribution("cell-annotation-schema").version
    headers, records = read_tsv_to_dict(data_file, generated_ids=True)
    config_fields = config["fields"]
    populate_labelsets(cas, config_fields)
    ao_names = dict()
    utilized_columns = set()
    for record_index in records:
        record = records[record_index]
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
                parent_ao = Annotation(
                    field["column_name"], record[field["column_name"]]
                )
                parent_ao.rank = int(str(field["rank"]).strip())
                parents.insert(int(str(field["rank"]).strip()), parent_ao)
                utilized_columns.add(field["column_name"])
            else:
                # handle annotation columns
                if "typing.List[str]" in str(get_type_hints(ao)[field["column_type"]]):
                    list_value = str(record[field["column_name"]]).split(",")
                    stripped = list(map(str.strip, list_value))
                    setattr(ao, field["column_type"], stripped)
                else:
                    setattr(ao, field["column_type"], record[field["column_name"]])
                utilized_columns.add(field["column_name"])

        add_user_annotations(ao, headers, record, utilized_columns)
        add_parent_node_names(ao, ao_names, cas, parents)

        ao_names[ao.cell_label] = ao
        cas.add_annotation_object(ao)
    return cas


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
    Creates parent nodes if necessary and creates a cluster hierarchy through assinging parent_node_names.
    :param ao: current annotation object
    :param ao_names: list of all created annotation objects
    :param cas: main object
    :param parents: list of current annotation object's parents
    """
    if parents:
        ao.parent_cell_set_name = parents[1].cell_label
        prev = None
        for parent in reversed(parents):
            if parent:
                if prev:
                    parent.parent_cell_set_name = prev.cell_label
                prev = parent
                if parent.cell_label not in ao_names:
                    cas.add_annotation_object(parent)
                    ao_names[parent.cell_label] = parent


def populate_labelsets(cas, config_fields):
    """
    Populates labelsets list based on the fields of the config.
    :param cas: main object
    :param config_fields: config file fields
    """
    labelsets = list()
    for field in config_fields:
        if field["column_type"] == "cell_set" or field["column_type"] == "cluster_name":
            label_set = Labelset(field["column_name"])
            if "rank" in field:
                label_set.rank = str(field["rank"])
            labelsets.append(label_set)
    if labelsets:
        cas.labelsets = labelsets
