import logging
import sys
from typing import Any, Dict, List, Optional, Tuple, Union

import pandas as pd
from cap_anndata import CapAnnDataDF

from cas.utils.conversion_utils import (
    ANNOTATIONS,
    CELL_LABEL,
    CELL_SET_ACCESSION,
    LABELSET,
    LABELSET_NAME,
    LABELSETS,
    PARENT_CELL_SET_ACCESSION,
)

# Configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


def validate_markers(
    cas: Dict[str, Any], adata: pd.DataFrame, marker_column: str
) -> bool:
    """
    Validates if the specified marker column in the anndata DataFrame contains all markers mentioned in the
    'annotations' in CAS dictionary. Raises an exception if the marker column does not exist.

    Args:
        cas : A dictionary containing various configurations and annotations, including marker gene evidence.
        adata : An anndata DataFrame.
        marker_column : The name of the column in `adata.var` which contains the list of marker genes to be validated.

    Returns:
        Returns True if all markers are validated without any issue, False or an exception otherwise.

    Raises:
        KeyError: If the specified `marker_column` is not found in `adata.var`.

    Note:
        This function uses the `validate_labelset_markers` to perform the actual validation per annotation entry in
        `cas`.
    """
    if marker_column not in adata.var.keys():
        raise KeyError(f"{marker_column} key does not exist in anndata.var!")

    marker_list = adata.var[marker_column].tolist()
    annotations = cas["annotations"]

    for annotation in annotations:
        validate_labelset_markers(annotation, marker_list)
    return True


def validate_labelset_markers(annotation: Dict[str, Any], marker_list: List[str]):
    """
    Validates if the markers from a specific annotation are present in the provided marker list. Logs a warning if
    any markers are missing.

    Args:
        annotation : A single annotation entry from CAS.
        marker_list : A list of marker genes to be checked against the markers in the annotation.

    Returns:
        This function does not return a value but will log a warning if validation fails.

    Note:
        This function is intended to be used within `validate_markers` to handle individual annotation validation.
    """
    marker_gene_evidence_list = annotation.get("marker_gene_evidence", [])

    if not marker_gene_evidence_list or marker_gene_evidence_list == [""]:
        return None

    marker_gene_evidence_set = set(marker_gene_evidence_list)
    marker_list_set = set(marker_list)

    # Determine the missing marker genes
    missing_markers = marker_gene_evidence_set - marker_list_set

    if missing_markers:
        logger.warning(
            f"Not all marker genes from {annotation['labelset']}-{annotation['cell_label']} pair exist in anndata's var section. "
            f"Missing markers: {', '.join(missing_markers)}"
        )


def validate_labelset_hierarchy(
    cas: Dict[str, Any], obs: Union[pd.DataFrame, CapAnnDataDF], validate: bool = False
) -> None:
    """
    Validates the labelset hierarchy by performing multiple consistency checks between CAS and obs.

    This function runs three validation checks:
    1. Ensures all labelsets from CAS exist in obs.
    2. Verifies that all labelset values from CAS annotations exist in the corresponding obs columns.
    3. Checks if the inferred parent-child hierarchy from obs matches CAS-defined ranks.

    If any of these checks fail, warnings or errors will be logged. If `validate=True`, the process
    will terminate with `sys.exit(1)` if any validation check fails.

    Parameters:
        cas (Dict[str, Any]): The CAS JSON object containing labelset definitions and annotations.
        obs (Union[pd.DataFrame, CapAnnDataDF]): The AnnData obs DataFrame or a CapAnnDataDF object.
        validate (bool, optional): If True, exits the program with an error code if any validation fails.
                                   Defaults to False.

    Returns:
        None
    """
    failures = []

    # Run all validation checks
    if not compare_labelsets_cas_obs(cas, obs):
        failures.append("Labelset comparison failed.")

    if not validate_labelset_values(cas, obs):
        failures.append("Labelset value validation failed.")

    if not check_parent_child_consistency(cas, obs):
        failures.append("Parent-child consistency check failed.")

    # Stop execution if validate=True and any failures occurred
    if failures and validate:
        logging.error("Validation failed. Exiting the process.")
        sys.exit(1)


def compare_labelsets_cas_obs(
    cas: Dict[str, Any], obs: Union[pd.DataFrame, CapAnnDataDF]
) -> bool:
    """
    Compare labelsets from CAS JSON object with the columns of the obs DataFrame.

    Logs a warning if any labelsets from CAS are missing in obs.

    Parameters:
        cas: The CAS JSON object.
        obs: The AnnData obs DataFrame.

    Returns:
        True if all labelsets from CAS exist as columns in obs, otherwise False.
    """
    labelset_list = [
        labelset[LABELSET_NAME] for labelset in cas[LABELSETS] if "rank" in labelset
    ]
    obs_columns = set(obs.columns)

    # Find missing labelsets
    missing_labelsets = set(labelset_list) - obs_columns

    if missing_labelsets:
        logging.warning(f"Missing labelsets in obs: {missing_labelsets}")
        return False
    else:
        logging.info("All labelsets exist in obs.")
        return True


def validate_labelset_values(
    cas: Dict[str, Any], obs: Union[pd.DataFrame, CapAnnDataDF]
) -> bool:
    """
    Validate that all labelset members from CAS annotations exist in the corresponding obs labelset columns.

    Logs warnings for any missing labelset members.

    Parameters:
        cas (Dict[str, Any]): The CAS JSON object.
        obs (pd.DataFrame): The AnnData obs DataFrame.

    Returns:
        True if all labelset members from CAS exist in obs, otherwise False.
    """
    labelset_list = [
        labelset[LABELSET_NAME] for labelset in cas[LABELSETS] if "rank" in labelset
    ]
    labelset_members = {}

    # Collect all labelset members from CAS annotations
    for annotation in cas[ANNOTATIONS]:
        labelset_name = annotation[LABELSET]
        if labelset_name not in labelset_list:
            continue
        cell_label = annotation[CELL_LABEL]

        if labelset_name not in labelset_members:
            labelset_members[labelset_name] = set()

        labelset_members[labelset_name].add(cell_label)

    # Check if all labelset members exist in obs
    missing_values = {}

    for labelset_name, members in labelset_members.items():
        if labelset_name in obs.columns:
            obs_values = set(obs[labelset_name].dropna().unique())
            missing_members = members - obs_values

            if missing_members:
                missing_values[labelset_name] = missing_members
        else:
            missing_values[labelset_name] = members

    # Log warnings for missing values
    if missing_values:
        for labelset, missing in missing_values.items():
            logging.warning(f"Missing labelset members in obs['{labelset}']: {missing}")
        return False
    else:
        logging.info("All labelset members exist in the corresponding obs columns.")
        return True


def check_parent_child_consistency(
    cas: Dict[str, Any], obs: Union[pd.DataFrame, CapAnnDataDF]
) -> bool:
    """
    Checks if the inferred hierarchy from cell labels in obs matches the expected hierarchy from
    CAS rank data.

    Parameters:
        cas (Dict[str, Any]): The CAS JSON object containing labelset ranks.
        obs (pd.DataFrame): The AnnData obs DataFrame.

    Returns:
        True if all inferred parent-child relationships from obs match those from cas,
              otherwise False.
    """
    # Extract ranks from CAS (store expected rank per cell label)
    cas_ranks = {
        entry["name"]: entry["rank"] for entry in cas[LABELSETS] if "rank" in entry
    }

    obs_inferred_hierarchy = infer_obs_cell_hierarchy(obs, cas_ranks)
    cas_inferred_hierarchy = infer_cas_cell_hierarchy(cas)

    # Compare inferred hierarchy with CAS ranks
    mismatches = []

    for child, parent in obs_inferred_hierarchy.items():
        if child in cas_inferred_hierarchy and parent != cas_inferred_hierarchy[child]:
            mismatches.append(
                f"{child}->{parent} relation from anndata does not match with {child}->"
                f"{cas_inferred_hierarchy[child]} relation from cas."
            )

    # Report mismatches
    if mismatches:
        for issue in mismatches:
            logging.warning(issue)
        return False
    else:
        logging.info("Parent-child relationships are consistent between CAS and OBS.")
        return True


def infer_obs_cell_hierarchy(
    obs: Union[pd.DataFrame, CapAnnDataDF], cas_ranks: Dict[str, int]
) -> Dict[str, Optional[str]]:
    """
    Infers a direct parent-child hierarchy between cell labels based on row co-occurrence.

    This function analyzes the hierarchical relationships between cell labels by comparing their
    row indices in the `obs` DataFrame. A label is considered a child if its row indices are fully
    contained within another label's indices. The closest (direct) parent is selected based on
    CAS-defined ranks.

    Parameters:
        obs (pd.DataFrame): The AnnData `obs` DataFrame containing labelset columns.
        cas_ranks (Dict[str, int]): A dictionary mapping cell labels to their CAS-defined rank,
                                    where lower values indicate higher ranks.

    Returns:
        Dict[Any, Optional[Any]]: A dictionary mapping each cell label to its inferred direct parent.
            Labels without a parent are assigned `None`.
    """
    labelset_list = [k for k, v in cas_ranks.items()]

    # Compute row indices for each cell label within each labelset
    cell_label_row_indices = {
        labelset: {
            label: set(obs.index[obs[labelset] == label])
            for label in obs[labelset].dropna().unique()
        }
        for labelset in labelset_list
    }

    child_parent_map = {}

    # Infer direct parent-child relationships based on row containment
    for labelset_child in labelset_list:
        for child_label, child_rows in cell_label_row_indices[labelset_child].items():
            possible_parents = {}

            for labelset_parent in labelset_list:
                for parent_label, parent_rows in cell_label_row_indices[
                    labelset_parent
                ].items():
                    if child_label == parent_label:
                        continue
                    if child_rows.issubset(parent_rows) and cas_ranks.get(
                        labelset_parent
                    ) > cas_ranks.get(labelset_child):
                        possible_parents[parent_label] = cas_ranks.get(labelset_parent)

            # Select the closest parent (lowest-ranked parent)
            if possible_parents:
                direct_parent = min(possible_parents, key=lambda k: possible_parents[k])
                child_parent_map[child_label] = direct_parent

    # Assign root labels (labels without a parent)
    all_labels = {
        label
        for labelset in labelset_list
        for label in cell_label_row_indices[labelset]
    }
    roots = {label for label in all_labels if label not in child_parent_map}

    for root in roots:
        child_parent_map[root] = None  # Root labels have no parent

    return child_parent_map


def infer_cas_cell_hierarchy(cas: Dict[str, Any]) -> Dict[str, Optional[str]]:
    child_parent_map = {}
    accession_lookup = {}
    annotations = cas[ANNOTATIONS]

    for annotation in annotations:
        accession_lookup[annotation[CELL_SET_ACCESSION]] = annotation[CELL_LABEL]
    for annotation in annotations:
        if PARENT_CELL_SET_ACCESSION in annotation:
            child_parent_map[annotation[CELL_LABEL]] = accession_lookup[
                annotation[PARENT_CELL_SET_ACCESSION]
            ]
        else:
            child_parent_map[annotation[CELL_LABEL]] = None

    return child_parent_map
