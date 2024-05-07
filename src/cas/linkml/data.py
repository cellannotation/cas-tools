import rdflib

from pathlib import Path
from typing import Union, Optional, List

from cas.file_utils import read_json_file

from linkml_runtime.utils.compile_python import compile_python
from linkml_runtime.linkml_model import SchemaDefinition
from linkml_runtime import SchemaView
from linkml_runtime.loaders import yaml_loader
from linkml_runtime.dumpers import rdflib_dumper
from linkml.validator import Validator
from linkml import generators


CAS_ROOT_CLASS = "GeneralCellAnnotationOpenStandard"

CAS_NAMESPACE = "https://cellular-semantics.sanger.ac.uk/ontology/CAS"
DEFAULT_PREFIXES = {
    "CAS": CAS_NAMESPACE + "/",
    CAS_ROOT_CLASS: CAS_NAMESPACE + "/",
    "obo": "http://purl.obolibrary.org/obo/",
    "CL": "http://purl.obolibrary.org/obo/CL_",
    "PCL": "http://purl.obolibrary.org/obo/PCL_",
    "RO": "http://purl.obolibrary.org/obo/RO_",
    "skos": "http://www.w3.org/2004/02/skos/core#",
}


def dump_to_rdf(
    schema: Union[str, Path, dict],
    instance: Union[str, dict],
    ontology_namespace: str,
    ontology_iri: str,
    labelsets: Optional[List[str]] = None,
    output_path: str = None,
    validate: bool = True,
) -> rdflib.Graph:
    """
    Dumps the given data to an RDF file based on the given schema file.
    Args:
        schema: The schema path/dict to be used for the RDF generation.
        instance: The data json file path or json object dict
        ontology_namespace: The namespace of the ontology (such as `MTG`).
        ontology_iri: The IRI of the ontology (such as `https://purl.brain-bican.org/ontology/AIT_MTG/`).
        labelsets: (Optional) The labelsets used in the taxonomy (such as `["Cluster", "Subclass", "Class"]`).
        output_path: (Optional) The output RDF file path.
        validate: (Optional) Boolean to determine if data-schema validation checks will be performed. True by default.

    Returns:
        RDFlib graph object
    """
    schema_def: SchemaDefinition = yaml_loader.load(
        schema, target_class=SchemaDefinition
    )

    if isinstance(instance, str):
        instance = read_json_file(instance)

    if validate:
        validate_data(schema_def, instance)

    gen = generators.PythonGenerator(schema_def)
    output = gen.serialize()
    python_module = compile_python(output)
    py_target_class = getattr(python_module, CAS_ROOT_CLASS)

    try:
        py_inst = py_target_class(**instance)
    except Exception as e:
        print(f"Could not instantiate {py_target_class} from the data; exception: {e}")
        return None

    prefixes = DEFAULT_PREFIXES.copy()
    prefixes["_base"] = ontology_iri
    prefixes[ontology_namespace] = ontology_iri
    for labelset in labelsets:
        prefixes[labelset] = ontology_iri + f"{labelset}#"

    g = rdflib_dumper.as_rdf_graph(
        py_inst,
        schemaview=SchemaView(schema_def),
        prefix_map=prefixes,
    )
    if output_path:
        g.serialize(format="xml", destination=output_path)
    return g


def validate_data(schema: SchemaDefinition, instance: dict) -> bool:
    """
    Validates the given data instance against the given schema.
    Args:
        schema: The schema to be used for the validation.
        instance: The data instance to be validated.

    Returns:
        Returns `True` if data is valid. Logs the validation errors and raises an exception if data is invalid.
    """
    validator = Validator(schema)
    report = validator.validate(instance)
    if report.results:
        print("Validation errors ({}):".format(len(report.results)))
        for result in report.results:
            print(result)
        raise ValueError(
            "Data file is not valid against the schema. {} validation errors found.".format(
                len(report.results)
            )
        )

    print("Data file is valid against the schema.")
    return True


def populate_ids(
    instance: Union[str, dict], ontology_namespace: str, ontology_id: str
) -> dict:
    """
    Population of id fields in the data instance that are required for the RDF conversion.
    Operation updates the instance object inplace if it is a dict.
    Args:
        instance: The data json file path or json object dict
        ontology_namespace: The namespace of the ontology (such as `MTG`).
        ontology_id: The ontology id to be used for the instance (such as `AIT_MTG`).

    Returns:
        json object with populated id properties
    """
    if isinstance(instance, str):
        instance = read_json_file(instance)

    if "id" in instance and instance["id"]:
        return instance

    if "id" not in instance:
        if "CAS:" not in ontology_id:
            ontology_id = "CAS:" + ontology_id
        instance["id"] = ontology_id

    for labelset in instance.get("labelsets", []):
        if "id" not in labelset:
            labelset["id"] = f"{ontology_namespace}:{labelset['name']}"

    # TODO add id to other properties as well

    return instance
