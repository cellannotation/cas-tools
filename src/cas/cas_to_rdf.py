import rdflib

from typing import Union, Optional, List

from cas.linkml_ops.schema import (
    convert_cas_schema_to_linkml,
    decorate_linkml_schema,
    expand_schema,
)
from cas.linkml_ops.data import dump_to_rdf, populate_ids


def export_to_rdf(
    cas_schema: Optional[Union[str, dict]],
    data: Union[str, dict],
    ontology_namespace: str,
    ontology_iri: str,
    labelsets: Optional[List[str]] = None,
    output_path: str = None,
    validate: bool = True,
    include_cells: bool = True,
) -> rdflib.Graph:
    """
    Generates and returns an RDF graph from the provided data and CAS schema, with an option to write the RDF graph to a file.

    Args:
        cas_schema: (Optional)
            Name of the CAS release (such as one of `base`, `cap`, `bican`) or path to the CAS schema file or
            cas schema json object. If not provided, reads the `base` CAS schema from the cas module.
        data:
            The data json file path or json object dict
        ontology_namespace: str
            Ontology namespace (e.g. `MTG`)
        ontology_iri: str
            Ontology IRI (e.g. `https://purl.brain-bican.org/ontology/AIT_MTG/`)
        labelsets: Optional[List[str]]
            Labelsets used in the taxonomy (such as `["Cluster", "Subclass", "Class"]`).
        output_path: (Optional)
            The output RDF file path.
        validate: (Optional)
            Boolean to determine if data-schema validation checks will be performed. True by default.
        include_cells: (Optional)
            Boolean to determine if cell data will be included in the RDF output. True by default.

    Returns:
        RDFlib graph object
    """
    # Prepare the linkml schema
    base_linkml_schema = convert_cas_schema_to_linkml(cas_schema)
    decorated_schema = decorate_linkml_schema(
        base_linkml_schema,
        ontology_namespace=ontology_namespace,
        ontology_iri=ontology_iri,
        labelsets=labelsets,
    )
    expanded_schema = expand_schema(
        config=None, yaml_obj=decorated_schema, value_set_names=["CellTypeEnum"]
    )

    # Prepare the data
    instance = populate_ids(
        data,
        ontology_namespace=ontology_namespace,
        ontology_id=ontology_namespace,
    )
    rdf_graph = dump_to_rdf(
        schema=expanded_schema,
        instance=instance,
        ontology_namespace=ontology_namespace,
        ontology_iri=ontology_iri,
        labelsets=labelsets,
        validate=validate,
        include_cells=include_cells,
        output_path=output_path,
    )

    return rdf_graph
