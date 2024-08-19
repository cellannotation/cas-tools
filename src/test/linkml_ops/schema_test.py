import os
import unittest

from cas.linkml_ops.schema import (
    convert_cas_schema_to_linkml,
    decorate_linkml_schema,
    expand_schema,
)

TESTDATA = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "../test_data/linkml/"
)

RAW_LINKML_SCHEMA = os.path.join(TESTDATA, "BICAN-schema-raw.yaml")


class LinkMLSchemaCase(unittest.TestCase):
    @classmethod
    def tearDownClass(cls):
        if os.path.isfile(RAW_LINKML_SCHEMA):
            os.remove(RAW_LINKML_SCHEMA)

    def test_convert_cas_schema_to_linkml(self):
        bican_linkml_schema = convert_cas_schema_to_linkml("bican", RAW_LINKML_SCHEMA)

        self.assertIsNotNone(bican_linkml_schema)
        self.assertTrue(os.path.isfile(RAW_LINKML_SCHEMA))

        self.assertEqual(
            "General_Cell_Annotation_Open_Standard", bican_linkml_schema.name
        )
        self.assertEqual(7, len(bican_linkml_schema.classes))
        self.assertEqual(49, len(bican_linkml_schema.slots))
        self.assertEqual(2, len(bican_linkml_schema.prefixes))

    def test_decorate_linkml_schema(self):
        bican_linkml_schema = convert_cas_schema_to_linkml("bican")
        decorated_schema = decorate_linkml_schema(
            bican_linkml_schema,
            ontology_namespace="MTG",
            ontology_iri="https://purl.brain-bican.org/ontology/AIT_MTG/",
            labelsets=["CrossArea_cluster", "CrossArea_subclass", "Class"],
        )
        self.assertIsNotNone(decorated_schema)

        self.assertEqual(
            "General_Cell_Annotation_Open_Standard", decorated_schema["name"]
        )
        self.assertEqual(7, len(decorated_schema["classes"]))
        self.assertEqual(50, len(decorated_schema["slots"]))
        self.assertEqual(13, len(decorated_schema["prefixes"]))

    def test_expand_schema(self):
        bican_linkml_schema = convert_cas_schema_to_linkml("bican")
        decorated_schema = decorate_linkml_schema(
            bican_linkml_schema,
            ontology_namespace="MTG",
            ontology_iri="https://purl.brain-bican.org/ontology/AIT_MTG/",
            labelsets=["CrossArea_cluster", "CrossArea_subclass", "Class"],
        )
        expanded_schema = expand_schema(
            config=None, yaml_obj=decorated_schema, value_set_names=["CellTypeEnum"]
        )
        self.assertIsNotNone(expanded_schema)

        self.assertEqual(2, len(expanded_schema["enums"]))
        self.assertTrue("CellTypeEnum" in expanded_schema["enums"])
        self.assertTrue(
            len(expanded_schema["enums"]["CellTypeEnum"]["permissible_values"]) > 20
        )


if __name__ == "__main__":
    unittest.main()
