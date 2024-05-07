import os
import unittest

from cas.linkml.schema import (
    convert_cas_schema_to_linkml,
    decorate_linkml_schema,
    expand_schema,
)
from cas.linkml.data import dump_to_rdf, populate_ids

TESTDATA = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "../test_data/linkml/"
)
TEST_OUTPUT = os.path.join(TESTDATA, "output.rdf")


class LinkMLIntegrationTestCase(unittest.TestCase):

    @classmethod
    def tearDownClass(cls):
        # if os.path.isfile(TEST_OUTPUT):
        #     os.remove(TEST_OUTPUT)
        pass

    def test_integration(self):
        ontology_namespace = "MTG"
        ontology_iri = "https://purl.brain-bican.org/ontology/AIT_MTG/"
        labelsets = ["CrossArea_cluster", "CrossArea_subclass", "Class"]

        # Prepare the schema
        bican_linkml_schema = convert_cas_schema_to_linkml("bican")
        decorated_schema = decorate_linkml_schema(
            bican_linkml_schema,
            ontology_namespace=ontology_namespace,
            ontology_iri=ontology_iri,
            labelsets=labelsets,
        )
        expanded_schema = expand_schema(
            config=None, yaml_obj=decorated_schema, value_set_names=["CellTypeEnum"]
        )

        # Prepare the data
        data = populate_ids(
            os.path.join(TESTDATA, "AIT_MTG.json"),
            ontology_namespace="MTG",
            ontology_id="AIT_MTG",
        )

        rdf_graph = dump_to_rdf(
            schema=expanded_schema,
            instance=data,
            ontology_namespace=ontology_namespace,
            ontology_iri=ontology_iri,
            labelsets=labelsets,
            validate=True,
            output_path=TEST_OUTPUT,
        )
        self.assertTrue(os.path.isfile(TEST_OUTPUT))


if __name__ == "__main__":
    unittest.main()
