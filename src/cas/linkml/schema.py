from pathlib import Path
from typing import Union, Optional, List
from ruamel.yaml import YAML

from cas.file_utils import get_cas_schema, get_cas_schema_names

from linkml_runtime.linkml_model import SchemaDefinition
from linkml_runtime.utils.schema_as_dict import schema_as_dict
from linkml_runtime.loaders import yaml_loader
from linkml_runtime.dumpers import json_dumper

from schema_automator.utils.schemautils import write_schema
from schema_automator.importers.jsonschema_import_engine import JsonSchemaImportEngine

from oaklib.utilities.subsets.value_set_expander import ValueSetExpander, ValueSetConfiguration


CAS_ROOT_CLASS = "GeneralCellAnnotationOpenStandard"

CAS_NAMESPACE = "https://cellular-semantics.sanger.ac.uk/ontology/CAS"
CAS_ROOT_NS = "General_Cell_Annotation_Open_Standard"
DEFAULT_PREFIXES = {
    CAS_ROOT_NS: CAS_NAMESPACE + "#",
    "CAS": CAS_NAMESPACE + "/",
    "obo": "http://purl.obolibrary.org/obo/",
    "CL": "http://purl.obolibrary.org/obo/CL_",
    "PCL": "http://purl.obolibrary.org/obo/PCL_",
    "RO": "http://purl.obolibrary.org/obo/RO_",
    "skos": "http://www.w3.org/2004/02/skos/core#",
}


def convert_cas_schema_to_linkml(
    cas_schema: Union[str, dict] = None, output_path: Optional[str] = None
) -> SchemaDefinition:
    """
    Converts a CAS schema to a LinkML schema.
    Args:
        cas_schema: Path to the CAS schema file or cas schema json object. If not provided, reads the base CAS schema
        from the cas module.
        output_path: (Optional) The output LinkML schema file path. (Default: None)

    Returns:
        LinkML SchemaDefinition object
    """
    if cas_schema is None or cas_schema in get_cas_schema_names().keys():
        cas_schema = get_cas_schema(cas_schema)

    ie = JsonSchemaImportEngine()
    if isinstance(cas_schema, dict):
        schema = ie.loads(
            cas_schema, name="cell-annotation-schema", root_class_name=None
        )
    else:
        schema = ie.load(
            cas_schema,
            name="cell-annotation-schema",
            format="json",
            root_class_name=None,
        )
    if output_path:
        write_schema(schema, output_path)

    return schema


def decorate_linkml_schema(
    schema_obj: Union[dict, SchemaDefinition],
    ontology_namespace: str,
    ontology_iri: str,
    labelsets: Optional[List[str]] = None,
    output_path: Optional[str] = None,
) -> dict:
    """
    Adds additional properties to the LinkML schema that are required for OWL conversion.
    Args:
        schema_obj: dict or SchemaDefinition
            Link ML schema object
        ontology_namespace: str
            Ontology namespace (e.g. `MTG`)
        ontology_iri: str
            Ontology IRI (e.g. `https://purl.brain-bican.org/ontology/AIT_MTG/`)
        labelsets: Optional[List[str]]
            Labelsets used in the taxonomy (such as `["Cluster", "Subclass", "Class"]`).
        output_path: Optional[str]
            (Optional) Path to the output schema file
    Returns:
        Decorated schema
    """
    if isinstance(schema_obj, SchemaDefinition):
        schema_obj = schema_as_dict(schema_obj)

    schema_obj["id"] = CAS_NAMESPACE

    ontology_namespace = ontology_namespace.upper()
    prefixes = DEFAULT_PREFIXES.copy()
    prefixes["linkml"] = "https://w3id.org/linkml/"
    prefixes["_base"] = ontology_iri
    prefixes[ontology_namespace] = ontology_iri
    labelsets = labelsets or []
    for labelset in labelsets:
        prefixes[labelset] = ontology_iri + f"{labelset}#"

    schema_obj["prefixes"] = prefixes
    schema_obj["default_range"] = "string"
    schema_obj["default_curi_maps"] = ["semweb_context", "obo_context"]
    schema_obj["enums"]["CellTypeEnum"] = {
        "reachable_from": {
            "source_ontology": "obo:cl",
            "source_nodes": ["CL:0000000"],
            "include_self": True,
            "relationship_types": ["rdfs:subClassOf"],
        }
    }
    schema_obj["slots"]["id"] = {"identifier": True, "range": "uriorcurie"}
    schema_obj["slots"]["name"]["slot_uri"] = "rdfs:label"
    schema_obj["slots"]["description"]["slot_uri"] = "IAO:0000115"
    schema_obj["slots"]["cell_label"]["slot_uri"] = "rdfs:label"
    schema_obj["slots"]["cell_fullname"]["slot_uri"] = "skos:preflabel"
    schema_obj["slots"]["cell_ontology_term_id"]["slot_uri"] = "RO:0002473"
    schema_obj["slots"]["cell_ontology_term_id"]["range"] = "CellTypeEnum"
    schema_obj["slots"]["cell_ids"]["slot_uri"] = "CAS:has_cellid"
    schema_obj["slots"]["cell_set_accession"]["identifier"] = True
    schema_obj["slots"]["parent_cell_set_accession"]["slot_uri"] = "RO:0015003"
    schema_obj["slots"]["parent_cell_set_accession"]["range"] = "Annotation"
    schema_obj["slots"]["source_taxonomy"]["range"] = "uriorcurie"
    schema_obj["slots"]["comment"]["slot_uri"] = "IAO:0000115"
    schema_obj["slots"]["labelset"]["slot_uri"] = "CAS:has_labelset"
    schema_obj["slots"]["labelset"]["range"] = "Labelset"
    schema_obj["slots"]["labelsets"]["inlined"] = True
    schema_obj["slots"]["annotations"]["inlined"] = True

    schema_obj["classes"]["Labelset"]["slots"] = list(
        schema_obj["classes"]["Labelset"]["slots"]
    ) + ["id"]
    schema_obj["classes"]["GeneralCellAnnotationOpenStandard"]["slots"] = list(
        schema_obj["classes"]["GeneralCellAnnotationOpenStandard"]["slots"]
    ) + ["id"]
    schema_obj["classes"]["Annotation"]["class_uri"] = "PCL:0010001"

    if output_path:
        write_schema(schema_obj, output_path)

    return schema_obj


def expand_schema(config: Optional[str], yaml_obj: dict, value_set_names: List[str], output_path: Optional[str] = None):
    """
    Dynamic Value set expander. Expands the yaml_obj schema inplace with the given value set names.
    Source code from https://github.com/INCATools/ontology-access-kit/blob/main/src/oaklib/utilities/subsets/value_set_expander.py

    Args:
        config: str
            Configuration file path
        yaml_obj: dict
            Yaml schema object
        value_set_names: List[str]
            Value set names to expand
        output_path: str
            (Optional) Output expanded schema file path
    """
    expander = ValueSetExpander()
    if config:
        expander.configuration = yaml_loader.load(config, target_class=ValueSetConfiguration)

    yaml = YAML()
    schema = yaml_loader.load(yaml_obj, target_class=SchemaDefinition)
    if value_set_names is None:
        value_set_names = list(schema.enums.keys())
    for value_set_name in value_set_names:
        if value_set_name not in schema.enums:
            raise ValueError(f"Unknown value set: {value_set_name}")
        value_set = schema.enums[value_set_name]
        pvs = list(expander.expand_value_set(value_set, schema=schema))
        yaml_obj["enums"][value_set_name]["permissible_values"] = {
            str(pv.text): json_dumper.to_dict(pv) for pv in pvs
        }
    if output_path:
        with open(output_path, "w", encoding="UTF-8") as file:
            yaml.dump(yaml_obj, stream=file)
    return yaml_obj
