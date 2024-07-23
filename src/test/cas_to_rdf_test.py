import os
import unittest

from rdflib import URIRef, Graph

from cas.cas_to_rdf import export_to_rdf


CAS_NS = "https://cellular-semantics.sanger.ac.uk/ontology/CAS/"
dataset_type = URIRef(CAS_NS + "GeneralCellAnnotationOpenStandard")
annotation_type = URIRef("http://purl.obolibrary.org/obo/PCL_0010001")
labelset_type = URIRef(CAS_NS + "Labelset")
rdftype = URIRef("http://www.w3.org/1999/02/22-rdf-syntax-ns#type")

TESTDATA = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "test_data/linkml/"
)
TEST_OUTPUT = os.path.join(TESTDATA, "output.rdf")
TEST_OUTPUT2 = os.path.join(TESTDATA, "CS202210140_non_neuronal.rdf")


class CAStoRDFTestCase(unittest.TestCase):

    @classmethod
    def tearDownClass(cls):
        if os.path.isfile(TEST_OUTPUT):
            os.remove(TEST_OUTPUT)
        if os.path.isfile(TEST_OUTPUT2):
            os.remove(TEST_OUTPUT2)

    def test_cas_to_rdf(self):
        ontology_namespace = "MTG"
        ontology_iri = "https://purl.brain-bican.org/ontology/AIT_MTG/"
        labelsets = ["CrossArea_cluster", "CrossArea_subclass", "Class"]

        print("XXXX FAILING PATH:")
        print(os.path.abspath(TEST_OUTPUT))
        print(os.path.abspath(TESTDATA))
        rdf_graph = export_to_rdf(
            cas_schema="bican",
            data=os.path.join(TESTDATA, "AIT_MTG.json"),
            ontology_namespace=ontology_namespace,
            ontology_iri=ontology_iri,
            labelsets=labelsets,
            output_path=TEST_OUTPUT,
            validate=True,
            include_cells=False
        )
        self.assertTrue(os.path.isfile(TEST_OUTPUT))

        self.assertEqual(1, len(list(rdf_graph.triples((None, rdftype, dataset_type)))))
        self.assertEqual(
            3, len(list(rdf_graph.triples((None, rdftype, labelset_type))))
        )
        self.assertEqual(
            160, len(list(rdf_graph.triples((None, rdftype, annotation_type))))
        )

        L5_IT_2 = URIRef(ontology_iri + "CrossArea_cluster#4a4d733723")
        triples = list(rdf_graph.triples((L5_IT_2, None, None)))
        self.assertEqual(5, len(triples))
        for triple in triples:
            if str(triple[1]) == CAS_NS + "has_labelset":
                self.assertEqual(ontology_iri + "CrossArea_cluster", str(triple[2]))
            elif str(triple[1]) == "http://www.w3.org/2000/01/rdf-schema#label":
                self.assertEqual("L5 IT_2", str(triple[2]))
            elif str(triple[1]) == "http://purl.obolibrary.org/obo/RO_0015003":
                self.assertEqual(ontology_iri + "CrossArea_subclass#c6694cb883", str(triple[2]))
            elif str(triple[1]) == "http://www.w3.org/2004/02/skos/core#preflabel":
                self.assertEqual("L5 IT_2", str(triple[2]))
            elif triple[1] == rdftype:
                self.assertEqual("http://purl.obolibrary.org/obo/PCL_0010001", str(triple[2]))
            else:
                self.fail("Unexpected triple: " + str(triple))

    def test_cas_to_rdf_siletti(self):
        ontology_namespace = "CS202210140"
        ontology_iri = "https://purl.brain-bican.org/ontology/CS202210140/"
        labelsets = ["Cluster", "supercluster_term"]

        rdf_graph = export_to_rdf(
            # cas_schema="bican",
            cas_schema="https://raw.githubusercontent.com/cellannotation/cell-annotation-schema/main/build/BICAN_schema.json",
            # data=os.path.join(TESTDATA, "CS202210140.json"),
            data=os.path.join(TESTDATA, "Siletti_all_non_neuronal_cells.json"),
            ontology_namespace=ontology_namespace,
            ontology_iri=ontology_iri,
            labelsets=labelsets,
            output_path=TEST_OUTPUT2,
            validate=True,
            include_cells=False
        )

        self.assertTrue(os.path.isfile(TEST_OUTPUT2))
        self.assertEqual(1, len(list(rdf_graph.triples((None, rdftype, dataset_type)))))
        self.assertEqual(
            2, len(list(rdf_graph.triples((None, rdftype, labelset_type))))
        )
        self.assertEqual(
            89, len(list(rdf_graph.triples((None, rdftype, annotation_type))))
        )

        Epen_69 = URIRef(ontology_iri + "CS202210140_70")
        triples = list(rdf_graph.triples((Epen_69, None, None)))
        self.assertEqual(5, len(triples))
        for triple in triples:
            if str(triple[1]) == CAS_NS + "has_labelset":
                self.assertEqual(ontology_iri + "Cluster", str(triple[2]))
            elif str(triple[1]) == "http://www.w3.org/2000/01/rdf-schema#label":
                self.assertEqual("Epen_69", str(triple[2]))
            elif str(triple[1]) == "http://purl.obolibrary.org/obo/RO_0015003":
                self.assertEqual(ontology_iri + "CS202210140_471", str(triple[2]))
            elif str(triple[1]) == CAS_NS + "author_annotation_fields":
                self.assertEqual("{\"Cluster ID\": \"69\"}", str(triple[2]))
            elif triple[1] == rdftype:
                self.assertEqual("http://purl.obolibrary.org/obo/PCL_0010001", str(triple[2]))
            else:
                self.fail("Unexpected triple: " + str(triple))


if __name__ == "__main__":
    unittest.main()
