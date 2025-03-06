import json
import logging
from typing import List, Optional, Union

from cap_anndata import CapAnnDataDF
from pandas import DataFrame

from cas.file_utils import read_anndata_file, read_json_file
from cas.utils.validation_utils import validate_labelset_hierarchy

# Configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


def populate_cell_ids(
    cas_json_path: str,
    anndata_path: str,
    labelsets: Optional[List[str]] = None,
    validate: bool = False,
):
    """
    Add/update CellIDs in a CAS JSON file using matching data from an AnnData file.

    This function reads a CAS JSON file and an AnnData file, validates their consistency, and updates
    CellIDs in CAS based on matching labelsets from the AnnData `obs` DataFrame. The modified CAS JSON
    is then saved back to the original file.

    Parameters:
        cas_json_path (str): Path to the CAS JSON file.
        anndata_path (str): Path to the AnnData file.
        labelsets (list, optional): A list of labelsets to update with CellIDs from AnnData. If None,
                                    the labelset with rank '0' is used by default.
        validate (bool, optional): If True, runs validation checks to ensure labelset consistency.
                                   The program will exit with an error if validation fails. Defaults to False.

    Raises:
        Exception: If the AnnData file cannot be read.

    Returns:
        None
    """
    ad = read_anndata_file(anndata_path)
    if ad is None:
        raise Exception(f"AnnData read operation failed: {anndata_path}")

    ad_obs = ad.obs
    cas = read_json_file(cas_json_path)

    cas = add_cell_ids(cas, ad_obs, labelsets, validate)

    if cas:
        with open(cas_json_path, "w") as json_file:
            json.dump(cas, json_file, indent=2)


def update_cas_with_cell_ids(
    cas_json: dict, anndata_obs: CapAnnDataDF, labelsets: list = None
) -> dict:
    """
    Update a CAS dictionary by adding or modifying CellIDs using matching AnnData observations.

    This function takes a CAS dictionary and an AnnData `obs` DataFrame and updates the CAS
    with cell IDs extracted from the specified labelsets in the AnnData.

    Parameters:
        cas_json (dict): The CAS dictionary to update with cell IDs from AnnData.
        anndata_obs (CapAnnDataDF): The `obs` DataFrame extracted from an AnnData object.
        labelsets (list, optional): A list of labelsets to update with IDs from AnnData.
            If None, the labelset with rank '0' is used.

    Returns:
        dict: The updated CAS dictionary with CellIDs populated.
    """
    # validation step
    validate_labelset_hierarchy(cas_json, anndata_obs)
    cas = add_cell_ids(cas_json, anndata_obs, labelsets)
    return cas


def add_cell_ids(
    cas: dict,
    ad_obs: Union[DataFrame, CapAnnDataDF],
    labelsets: list = None,
    validate: bool = False,
):
    """
    Add/update CellIDs to CAS from matching AnnData file.

    Parameters:
        cas: CAS JSON object
        ad_obs: Obs DataFrame extracted from an AnnData object.
        labelsets: List of labelsets to update with IDs from AnnData. If value is null, rank '0' labelset is used. The
        labelsets should be provided in order, starting from rank 0 (leaf nodes) and ascending to higher ranks.
        validate (bool, optional): If True, runs validation checks to ensure labelset consistency.
                           The program wil`l exit with an error if validation fails. Defaults to False.
    """
    # Run validation checks (exit if validation fails and validate=True)
    validate_labelset_hierarchy(cas, ad_obs, validate)

    rank_zero_labelset = [
        lbl_set["name"]
        for lbl_set in cas["labelsets"]
        if "rank" in lbl_set and lbl_set.get("rank") in (0, "0")
    ][0]
    if not labelsets:
        labelsets = [rank_zero_labelset]

    obs_keys = ad_obs.columns.tolist()
    cluster_identifier_column = get_obs_cluster_identifier_column(
        obs_keys, labelsets, rank_zero_labelset
    )

    if cluster_identifier_column:
        cid_lookup = {}
        for anno in cas["annotations"]:
            if anno["labelset"] == rank_zero_labelset and anno["labelset"] in labelsets:
                cell_ids = []
                if cluster_identifier_column.endswith(
                    "_id"
                ) or cluster_identifier_column.lower().endswith(" id"):
                    # cluster column value is integer cluster id
                    cluster_id = anno.get("author_annotation_fields", {}).get(
                        cluster_identifier_column
                    )
                    if not cluster_id:
                        raise ValueError(
                            "AnnData cluster identifier column ({}) couldn't be find found in"
                            "CAS author_annotation_fields.".format(
                                cluster_identifier_column
                            )
                        )
                    cell_ids = list(
                        ad_obs.loc[
                            ad_obs[cluster_identifier_column] == int(cluster_id),
                            cluster_identifier_column,
                        ].index
                    )
                else:
                    # cluster column value is cell label
                    cluster_label = anno["cell_label"]
                    cell_ids = list(
                        ad_obs.loc[
                            ad_obs[cluster_identifier_column] == cluster_label
                        ].index
                    )
                anno["cell_ids"] = cell_ids
                if "parent_cell_set_name" in anno:
                    lookup_key = anno["parent_cell_set_name"]
                    if lookup_key in cid_lookup:
                        cid_lookup[lookup_key].update(cell_ids)
                    else:
                        cid_lookup[lookup_key] = set(cell_ids)

        for anno in cas["annotations"]:
            if anno["labelset"] in labelsets and anno["cell_label"] in cid_lookup:
                cell_ids = list(cid_lookup.get(anno["cell_label"], []))
                anno["cell_ids"] = cell_ids

        return cas
    else:
        logger.warning(
            "Cluster identifier column couldn't be identified in OBS. Populate cell ids operation is aborted."
        )

    return None


def get_obs_cluster_identifier_column(
    obs_keys: List[str], labelsets: list = None, rank_zero_labelset: str = None
):
    """
    Anndata files may use different column names to uniquely identify Clusters. Get the cluster identifier column name for the current file.
    Args:
        obs_keys: Anndata observation keys.
        labelsets: List of labelsets to update with IDs from AnnData. The labelsets should be provided in order,
        starting from rank 0 (leaf nodes) and ascending to higher ranks.
        rank_zero_labelset: rank 0 labelset name
    Returns:
        cluster identifier column name
    """
    # obs_keys = ad.obs_keys()
    cluster_identifier_column = ""
    if labelsets and labelsets[0] in obs_keys:
        cluster_identifier_column = labelsets[0]
    else:
        # check id variants
        if rank_zero_labelset.title() + "_id" in obs_keys:
            cluster_identifier_column = rank_zero_labelset.title() + "_id"
        elif rank_zero_labelset.lower() + "_id" in obs_keys:
            cluster_identifier_column = rank_zero_labelset.lower() + "_id"
        elif rank_zero_labelset.title() + " id" in obs_keys:
            cluster_identifier_column = rank_zero_labelset.title() + " id"
        elif rank_zero_labelset.lower() + " id" in obs_keys:
            cluster_identifier_column = rank_zero_labelset.lower() + " id"
    return cluster_identifier_column
