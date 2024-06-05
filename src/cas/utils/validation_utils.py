import logging
from typing import Any, Dict, List

import pandas as pd


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
    marker_gene_evidence_list = (
        annotation.get("marker_gene_evidence", [])
    )

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
