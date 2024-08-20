import json
import os
import unittest

import anndata as ad
import numpy as np
import pandas as pd

from cas.flatten_data_to_anndata import flatten, unflatten

# File paths
anndata_file_path = "src/test/test_data/cas_cap_roundtrip/test_data.h5ad"
flattened_anndata_file_path = (
    "src/test/test_data/cas_cap_roundtrip/flattened_test_data.h5ad"
)
unflattened_anndata_file_path = (
    "src/test/test_data/cas_cap_roundtrip/unflattened_test_data.h5ad"
)
cas_file_path = "src/test/test_data/cas_cap_roundtrip/test_cas.json"


class TestRoundtrip(unittest.TestCase):
    def test_roundtrip_without_edit(self):
        # Perform flattening and unflattening
        flatten(cas_file_path, anndata_file_path, flattened_anndata_file_path)
        unflatten(
            json_file_path=None,
            anndata_file_path=flattened_anndata_file_path,
            output_file_path=unflattened_anndata_file_path,
            output_json_path="cas.json",
        )

        # Read the input and output AnnData objects
        input_anndata = ad.read_h5ad(anndata_file_path, backed="r+")
        output_anndata = ad.read_h5ad(unflattened_anndata_file_path, backed="r+")

        # Compare the two AnnData objects
        self.compare_anndata(input_anndata, output_anndata)

        # Close the AnnData objects
        input_anndata.file.close()
        output_anndata.file.close()

    def test_roundtrip_with_edit(self):
        # Perform flattening, update on flattened data and unflattening
        flatten(cas_file_path, anndata_file_path, flattened_anndata_file_path)
        # Read the flattened Anndata object
        flattened_anndata = ad.read_h5ad(flattened_anndata_file_path, backed="r+")
        # rename a labelset cell label
        flattened_anndata.obs["Cluster"] = flattened_anndata.obs["Cluster"].replace(
            "O50", "O500x"
        )
        if isinstance(flattened_anndata.obs["Cluster"], pd.CategoricalDtype):
            flattened_anndata.obs["Cluster"].cat.add_categories("YYY", inplace=True)
            flattened_anndata.obs["Cluster"] = flattened_anndata.obs[
                "Cluster"
            ].cat.remove_unused_categories()

        flattened_anndata.write_h5ad(flattened_anndata_file_path)

        unflatten(
            json_file_path=None,
            anndata_file_path=flattened_anndata_file_path,
            output_file_path=unflattened_anndata_file_path,
            output_json_path="cas.json",
        )

        # Read the unflattened AnnData objects
        unflattened_anndata = ad.read_h5ad(unflattened_anndata_file_path, backed="r+")

        cas = json.loads(unflattened_anndata.uns["cas"])
        # Check the updated "annotations"
        self.assertEqual(len(cas["annotations"]), 5)
        self.assertIn("O500x", [a["cell_label"] for a in cas["annotations"]])
        self.assertIn("O40", [a["cell_label"] for a in cas["annotations"]])
        self.assertIn("A62", [a["cell_label"] for a in cas["annotations"]])
        self.assertIn("Oligodendrocyte", [a["cell_label"] for a in cas["annotations"]])
        self.assertIn("Astrocyte", [a["cell_label"] for a in cas["annotations"]])

        # Close the AnnData objects
        flattened_anndata.file.close()
        unflattened_anndata.file.close()

    def compare_anndata(self, input_anndata, output_anndata):
        # Align categories in both DataFrames
        input_obs = align_categories(input_anndata.obs.copy())
        output_obs = align_categories(output_anndata.obs.copy())
        # Compare observations (obs)
        pd.testing.assert_frame_equal(
            input_obs, output_obs, check_like=True, obj="obs data", check_dtype=False
        )

        # Compare variables (var)
        pd.testing.assert_frame_equal(
            input_anndata.var,
            output_anndata.var,
            check_like=True,
            obj="var data",
            check_dtype=False,
        )

        # Compare unstructured data (uns)
        for key in input_anndata.uns.keys():
            np.testing.assert_equal(
                input_anndata.uns[key],
                output_anndata.uns[key],
                err_msg=f"uns[{key}] is not equal",
            )

        # Compare obsm
        self.assertEqual(
            input_anndata.obsm.keys(),
            output_anndata.obsm.keys(),
            "obsm keys do not match",
        )
        for key in input_anndata.obsm.keys():
            np.testing.assert_array_equal(
                input_anndata.obsm[key],
                output_anndata.obsm[key],
                err_msg=f"obsm[{key}] is not equal",
            )

        # Compare varm
        self.assertEqual(
            input_anndata.varm.keys(),
            output_anndata.varm.keys(),
            "varm keys do not match",
        )
        for key in input_anndata.varm.keys():
            np.testing.assert_array_equal(
                input_anndata.varm[key],
                output_anndata.varm[key],
                err_msg=f"varm[{key}] is not equal",
            )

        # Compare obsp
        self.assertEqual(
            input_anndata.obsp.keys(),
            output_anndata.obsp.keys(),
            "obsp keys do not match",
        )
        for key in input_anndata.obsp.keys():
            np.testing.assert_array_equal(
                input_anndata.obsp[key],
                output_anndata.obsp[key],
                err_msg=f"obsp[{key}] is not equal",
            )

        # Compare varp
        self.assertEqual(
            input_anndata.varp.keys(),
            output_anndata.varp.keys(),
            "varp keys do not match",
        )
        for key in input_anndata.varp.keys():
            np.testing.assert_array_equal(
                input_anndata.varp[key],
                output_anndata.varp[key],
                err_msg=f"varp[{key}] is not equal",
            )

    def tearDown(self):
        # Delete files created during the test
        if os.path.exists(flattened_anndata_file_path):
            os.remove(flattened_anndata_file_path)

        if os.path.exists(unflattened_anndata_file_path):
            os.remove(unflattened_anndata_file_path)


def align_categories(df):
    """
    Align the categories of all categorical columns in a DataFrame.
    """
    for col in df.select_dtypes(["category"]).columns:
        # Get the combined categories from both DataFrames for the same column
        combined_categories = sorted(df[col].cat.categories)
        # Set the same categories and order for both DataFrames
        df[col] = df[col].cat.set_categories(combined_categories)
    return df
