import os
import unittest

from cas.linkml.data import dump_to_rdf, populate_ids
from rdflib import Graph


TESTDATA = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "../test_data/linkml/"
)
TEST_OUTPUT = os.path.join(TESTDATA, "output_small.rdf")


class LinkMLDataTestCase(unittest.TestCase):

    @classmethod
    def tearDownClass(cls):
        # if os.path.isfile(TEST_OUTPUT):
        #     os.remove(TEST_OUTPUT)
        pass

    def test_rdf_dump(self):
        data = populate_ids(os.path.join(TESTDATA, "AIT_MTG_data_short.json"), ontology_namespace="MTG", ontology_id="AIT_MTG")

        rdf_graph = dump_to_rdf(
            schema=os.path.join(TESTDATA, "BICAN-schema-expanded-MTG.yaml"),
            instance=data,
            ontology_namespace="MTG",
            ontology_iri="https://purl.brain-bican.org/ontology/AIT_MTG/",
            labelsets=["CrossArea_cluster", "CrossArea_subclass", "Class"],
            output_path=TEST_OUTPUT,
            validate=False,

        )
        self.assertIsNotNone(rdf_graph)

        # TODO improve assertions
        # expected_graph = Graph()
        # expected_graph.parse(
        #     os.path.join(TESTDATA, "expected_ouput_1.owl"), format="xml"
        # )
        #
        # self.assertEqual(len(rdf_graph), len(expected_graph))
        # for stmt in expected_graph:
        #     self.assertTrue(stmt in rdf_graph)

    def test_validate(self):
        # TODO: Implement test
        pass


if __name__ == "__main__":
    unittest.main()
