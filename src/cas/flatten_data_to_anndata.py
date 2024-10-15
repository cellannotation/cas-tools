import json
import logging
from collections import defaultdict
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from cap_anndata import read_h5ad

from cas.accession.hash_accession_manager import HashAccessionManager, is_hash_accession
from cas.file_utils import (
    read_json_file,
    update_obs,
    update_uns,
    write_dict_to_json_file,
)
from cas.utils.conversion_utils import (
    ANNOTATIONS,
    AUTHOR_ANNOTATION_FIELDS,
    CELL_IDS,
    CELL_LABEL,
    CELLHASH,
    LABELSET,
    LABELSET_NAME,
    LABELSETS,
    collect_parent_cell_ids,
    copy_and_update_file_path,
    fetch_anndata,
    reformat_json,
    convert_complex_type,
)

# Configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


def is_list_of_strings(var):
    """
    Check if a value is a list of strings.

    Parameters:
        var (list or any): The value to be checked.

    Returns:
        bool: True if the value is a list containing only string elements,
              False otherwise.

    """
    return isinstance(var, list) and all(isinstance(item, str) for item in var)


def flatten(
    cas_file_path: str,
    anndata_file_path: Optional[str],
    output_file_path: str,
    fill_na: bool,
):
    """
     Processes and integrates information from a JSON file and an AnnData (Annotated Data) file, creating a new AnnData
     object that incorporates the metadata. The resulting AnnData object is then saved to a new file.

    Args:
        cas_file_path: The path to the CAS json file.
        anndata_file_path: The path to the AnnData file.
        output_file_path: Output AnnData file name.
        fill_na: Boolean flag indicating whether to fill missing values in the 'obs' field with pd.NA. If True, missing
                 values will be replaced with pd.NA; if False, they will remain as empty strings.
    """
    input_json = read_json_file(cas_file_path)

    flatten_cas_object(input_json, anndata_file_path, output_file_path, fill_na)


def flatten_cas_object(
    input_json: dict,
    anndata_file_path: Optional[str],
    output_file_path: str,
    fill_na: bool,
):
    """
     Processes and integrates information from a JSON file and an AnnData (Annotated Data) file, creating a new AnnData
     object that incorporates the metadata. The resulting AnnData object is then saved to a new file.

    Args:
        input_json: CAS json file.
        anndata_file_path: The path to the AnnData file.
        output_file_path: Output AnnData file name.
        fill_na: Boolean flag indicating whether to fill missing values in the 'obs' field with pd.NA. If True, missing
                 values will be replaced with pd.NA; if False, they will remain as empty strings.
    """
    if not anndata_file_path:
        anndata_file_path = fetch_anndata(input_json)
    anndata_file_path = copy_and_update_file_path(anndata_file_path, output_file_path)

    annotations = input_json[ANNOTATIONS]
    parent_cell_ids = collect_parent_cell_ids(input_json)

    with read_h5ad(file_path=anndata_file_path, edit=True) as cap_adata:
        cap_adata.read_obs()
        obs = cap_adata.obs
        obs_index = np.array(cap_adata.obs.axes[0].tolist())

        # obs
        flatten_data = process_annotations(
            annotations, obs_index, parent_cell_ids, fill_na
        )
        update_obs(obs, flatten_data)

        # uns
        uns_json = generate_uns_json(input_json)
        cap_adata.read_uns()
        uns = cap_adata.uns
        update_uns(uns, uns_json)

        cap_adata.overwrite()


def process_annotations(annotations, obs_index, parent_cell_ids, fill_na):
    """
    Processes annotations and generates flattened data for obs dataset.

    Args:
        annotations (list): List of annotations.
        obs_index (np.ndarray): Array representing the index of the obs dataset.
        parent_cell_ids (dict): Dictionary containing parent cell ids.
        fill_na (bool):

    Returns:
        dict: Dictionary containing flattened data.
    """
    accession_manager = HashAccessionManager()
    flatten_data = {}
    for ann in annotations:
        cell_ids = ann.get(CELL_IDS, [])
        if not cell_ids:
            cell_ids = parent_cell_ids.get(ann.get("cell_set_accession", []))

        author_annotations = ann.get(AUTHOR_ANNOTATION_FIELDS, {})
        author_annotations.update(
            {
                CELLHASH: ann.get("cell_set_accession")
                if is_hash_accession(ann.get("cell_set_accession", None))
                else accession_manager.generate_accession_id(cell_ids=cell_ids)
            }
        )
        ann[AUTHOR_ANNOTATION_FIELDS] = author_annotations

        if not cell_ids:
            # only happens if data has multi-inheritance (as in basal ganglia data)
            continue
        # Convert cell_ids to a list if it's not already for np.isin
        if not isinstance(cell_ids, list):
            cell_ids = list(cell_ids)
        mask = np.isin(obs_index, cell_ids)

        for k, v in ann.items():
            if k in [CELL_IDS, LABELSET]:
                continue

            key = f"{ann[LABELSET]}--{k}" if k != CELL_LABEL else ann[LABELSET]
            value = ", ".join(
                sorted([str(value) for value in v] if isinstance(v, list) else [str(v)])
            )

            if key not in flatten_data:
                flatten_data[key] = pd.Series("", index=obs_index)
            flatten_data[key].loc[mask] = value
            if fill_na:
                flatten_data[key].loc[~mask] = pd.NA

    # Convert relevant columns to categorical after the loop
    for key in flatten_data:
        # Get unique values and convert the Series to categorical
        unique_values = pd.unique(flatten_data[key])
        unique_values = unique_values[~pd.isna(unique_values)]
        flatten_data[key] = pd.Series(
            pd.Categorical(flatten_data[key], categories=unique_values), index=obs_index
        )
    return flatten_data


def generate_uns_json(input_json):
    """
    Generates a dictionary representing the uns (unstructured) field in an AnnData object from a given JSON input.

    This function processes information from a JSON input and generates a dictionary that represents the uns (unstructured)
    field in an AnnData object. The resulting dictionary can be used to populate the uns field in the AnnData object.

    Args:
        input_json (dict): A dictionary representing the input CAS JSON data containing annotations.

    Returns:
        dict: A dictionary representing the uns (unstructured) field in an AnnData object, ready to be used as input
              for writing to an AnnData file.

    """
    uns_json = {}
    root_keys = list(input_json.keys())
    root_keys.remove(ANNOTATIONS)

    for key in root_keys:
        value = input_json[key]
        if not value:
            continue

        if is_list_of_strings(value):
            uns_json[key] = ", ".join(sorted(value))
        elif isinstance(value, str):
            uns_json[key] = value
        else:
            metadata_json = {}
            for labelset in value:
                metadata_key = labelset.get(LABELSET_NAME, "")
                metadata_json.update({metadata_key: {}})
                for k, v in labelset.items():
                    if k == LABELSET_NAME:
                        continue
                    metadata_json.get(metadata_key, {}).update({k: v})
            uns_json["cellannotation_metadata"] = metadata_json

    uns_json["cas"] = reformat_json(input_json)

    return uns_json


def unflatten(
    json_file_path: Optional[str],
    anndata_file_path: str,
    output_file_path: str,
    output_json_path: str,
):
    """
     Unflatten an Anndata file and save it. Also creates a CAS json file as output.

    Args:
        json_file_path: The path to the CAS json file.
        anndata_file_path: The path to the AnnData file.
        output_file_path: Output AnnData file name.
        output_json_path: Output CAS JSON file name.
    """
    # TODO review the `cas = None` and `cap_adata.read_uns()` logic!!!
    anndata_file_path = copy_and_update_file_path(anndata_file_path, output_file_path)
    cas = None
    if json_file_path:
        with open(json_file_path, "r") as file:
            cas = json.load(file)

    with read_h5ad(file_path=anndata_file_path, edit=True) as cap_adata:
        if not cas:
            cap_adata.read_uns()
            if "cas" not in cap_adata.uns:
                raise KeyError(
                    "uns section does not have a CAS section. Please check the AnnData file or provide a valid CAS file."
                )
            cas = json.loads(cap_adata.uns["cas"])

        cap_adata.read_obs()
        obs = cap_adata.obs
        new_cas = unflatten_obs(obs, cas)

        cap_adata.uns["cas"] = reformat_json(new_cas)
        # Save your changes to a new or the same AnnData file
        cap_adata.overwrite()

        # Write new cas json to file
        write_dict_to_json_file(output_json_path, new_cas)


def unflatten_obs(obs_df: pd.DataFrame, cas_json: Dict[str, Any]) -> Dict[str, Any]:
    """
    Reverse the flattening process to update the "annotations" section in a CAS object.

    Args:
        obs_df: DataFrame containing the flattened obs columns from an AnnData object.
        cas_json: CAS JSON object.

    Returns:
        Updated CAS JSON with revised annotations.
    """
    # Create a dictionary with a list of columns for each labelset and its dataframe
    labelsets = [labelset[LABELSET_NAME] for labelset in cas_json[LABELSETS]]
    obs_columns_by_labelset = {
        labelset: [col for col in obs_df.columns if labelset in col]
        for labelset in labelsets
    }
    filtered_obs_by_labelset = {
        labelset: obs_df[columns]
        for labelset, columns in obs_columns_by_labelset.items()
    }
    # Find all matching cell sets defined by obs
    cas_dict = create_cell_label_lookup(filtered_obs_by_labelset)
    # Check cell set membership and update cas
    updated_cas = update_cas_json(cas_dict, cas_json)
    # Discard flattened obs
    flattened_columns = [
        col
        for labelset, column_list in obs_columns_by_labelset.items()
        for col in column_list
        if "--" in col
    ]
    for flattened_column in flattened_columns:
        obs_df.remove_column(flattened_column)

    return updated_cas


def create_cell_label_lookup(df_dict: Dict[str, pd.DataFrame]) -> dict:
    """
    Create a lookup dictionary for cell labels with corresponding observations.

    Args:
        df_dict: A dictionary of DataFrames keyed by label sets.

    Returns:
        A nested dictionary where keys are cell labels and values are the observations
        from the obs field in the AnnData object.
    """
    accession_manager = HashAccessionManager()
    # Initialize the dictionary with defaultdict for automatic dictionary creation
    nested_dict = defaultdict(lambda: {CELL_IDS: [], AUTHOR_ANNOTATION_FIELDS: {}})

    # Process each DataFrame in the dictionary
    for df_key, df in df_dict.items():
        grouped = df.groupby(df_key, observed=True)

        for key, group in grouped:
            # Append indices of each group to cell_ids
            cell_id_list = group.index.tolist()
            key_pair = f"{df_key}:{key}"
            nested_dict[key_pair][CELL_IDS].extend(cell_id_list)

            # Process each column in the group
            for col in df.columns:
                if col == df_key:
                    # Store the labelset and cell label for the first column
                    nested_dict[key_pair][LABELSET] = col
                    nested_dict[key_pair][CELL_LABEL] = convert_complex_type(
                        group[col].iloc[0]
                    )
                    # Store cellhash
                    cell_hash = accession_manager.generate_accession_id(
                        cell_ids=cell_id_list, labelset=col, suppress_warnings=True
                    )
                    nested_dict[key_pair][AUTHOR_ANNOTATION_FIELDS][
                        CELLHASH
                    ] = cell_hash
                else:
                    # Split annotation columns and store them in annotations
                    annotation_column = col.split("--")[-1]
                    if annotation_column == AUTHOR_ANNOTATION_FIELDS:
                        annotation_dict = json.loads(
                            group[col].iloc[0].replace("'", '"')
                        )
                        filtered_annotation_dict = {
                            k: v for k, v in annotation_dict.items() if k != CELLHASH
                        }
                        nested_dict[key_pair][annotation_column].update(
                            filtered_annotation_dict
                        )
                    else:
                        nested_dict[key_pair][annotation_column] = convert_complex_type(
                            group[col].iloc[0]
                        )

    return {
        key_: dict(value)
        for key, value in nested_dict.items()
        for key_ in (key, value[AUTHOR_ANNOTATION_FIELDS][CELLHASH])
    }


def update_cas_json(
    cas_dict: Dict[str, Dict[str, Any]], cas_json: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Update the annotations in the CAS JSON using the provided lookup dictionary.

    This function checks the CAS JSON annotations against a lookup dictionary. It updates
    annotations where cell labels or cell hashes match and discards mismatches. It also
    adds new annotations from the lookup dictionary that are not in the CAS JSON.

    Args:
        cas_dict: A lookup dictionary where keys are cell labels or hashes, and
                  values are dictionaries with annotation data.
        cas_json: The CAS JSON object containing existing annotations.

    Returns:
        The updated CAS JSON with synchronized annotations based on the lookup dictionary.
    """
    updated_cas_annotations: List[Dict[str, Any]] = []
    remaining_annotations_labels = list(key for key in cas_dict.keys())
    # remaining_annotations_labels = list(key for key in cas_dict.keys() if not is_hash_accession(key))

    for annotation in cas_json[ANNOTATIONS]:
        labelset_cell_label_pair = f"{annotation[LABELSET]}:{annotation[CELL_LABEL]}"
        cas_cellhash = annotation[AUTHOR_ANNOTATION_FIELDS][CELLHASH]

        if labelset_cell_label_pair in cas_dict.keys():
            obs_annotation = cas_dict[labelset_cell_label_pair]
            obs_cellhash = cas_dict[labelset_cell_label_pair][AUTHOR_ANNOTATION_FIELDS][
                CELLHASH
            ]

            # TODO Make validation optional
            if cas_cellhash == obs_cellhash:
                updated_cas_annotations.append(obs_annotation)
            else:
                logging.warning(
                    f"Cell set annotations for {labelset_cell_label_pair} are discarded because cell hashes do not match."
                )
            index_to_remove = remaining_annotations_labels.index(
                labelset_cell_label_pair
            )
            remaining_annotations_labels.pop(index_to_remove)
            remaining_annotations_labels.pop(index_to_remove)
            # remaining_annotations_labels.remove(labelset_cell_label_pair)
        elif cas_cellhash in cas_dict.keys():
            obs_annotation = cas_dict[cas_cellhash]
            updated_cas_annotations.append(obs_annotation)
            logger.warning(
                f"Cell id hashes, {cas_cellhash}, match but cell labels, {labelset_cell_label_pair}->"
                f"{obs_annotation[CELL_LABEL]}, do not. "
                f"Annotations are updated anyway. There might be a change in Cell Label field!"
            )
            index_to_remove = remaining_annotations_labels.index(cas_cellhash)
            remaining_annotations_labels.pop(index_to_remove)
            remaining_annotations_labels.pop(index_to_remove - 1)

    for label in remaining_annotations_labels:
        updated_cas_annotations.append(cas_dict[label])
        logger.warning(f"New cell set has been added with label {label}.")

    return {
        k: (updated_cas_annotations if k == ANNOTATIONS else v)
        for k, v in cas_json.items()
    }
