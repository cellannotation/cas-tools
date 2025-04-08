import logging
import sys
from typing import Optional
import warnings

from cap_anndata import read_h5ad

from cas.file_utils import read_json_file
from cas.utils.conversion_utils import (
    ANNOTATIONS,
    CELL_IDS,
    CELL_LABEL,
    LABELSET,
    LABELSET_NAME,
    LABELSETS,
    collect_parent_cell_ids,
    copy_and_update_file_path,
    fetch_anndata,
    reformat_json,
)

# Set up logging
logger = logging.getLogger(__name__)

# Suppress warning messages from cap_anndata.cap_anndata
logging.getLogger("cap_anndata.cap_anndata").setLevel(logging.ERROR)


def merge(
    cas_file_path: str,
    anndata_path: Optional[str],
    validate: bool,
    output_file_name: str,
):
    """
    Tests if CAS json and AnnData are compatible and merges CAS into AnnData if possible.

    This function performs the following checks:
        1. Verifies that all cell barcodes (cell IDs) in CAS exist in AnnData and vice versa.
        2. Identifies matching labelset names between CAS and AnnData.
        3. Validates that cell sets associated with each annotation match between CAS and AnnData.
        4. Checks if the cell labels are identical; if not, provides options to update or terminate.

    Args:
        cas_file_path: The path to the CAS json file.
        anndata_path: The path to the AnnData file.
        validate: Boolean to determine if validation checks will be performed before writing to the output AnnData file.
        output_file_name: Output AnnData file name.

    """
    input_json = read_json_file(cas_file_path)

    merge_cas_object(input_json, anndata_path, validate, output_file_name)


def merge_cas_object(
    input_json: dict,
    anndata_file_path: Optional[str],
    validate: bool,
    output_file_path: str,
    download_dir: Optional[str] = None,
):
    """
    Tests if CAS json and AnnData are compatible and merges CAS into AnnData if possible.

    This function performs the following checks:
        1. Verifies that all cell barcodes (cell IDs) in CAS exist in AnnData and vice versa.
        2. Identifies matching labelset names between CAS and AnnData.
        3. Validates that cell sets associated with each annotation match between CAS and AnnData.
        4. Checks if the cell labels are identical; if not, provides options to update or terminate.

    Args:
        input_json: The CAS json object.
        anndata_file_path: The path to the AnnData file.
        validate: Boolean to determine if validation checks will be performed before writing to the output AnnData file.
        output_file_path: Output AnnData file name.
        download_dir: The directory to download AnnData files.

    """
    if not anndata_file_path:
        anndata_file_path = fetch_anndata(input_json, download_dir)
    anndata_file_path = copy_and_update_file_path(anndata_file_path, output_file_path)

    with read_h5ad(file_path=anndata_file_path, edit=True) as cap_adata:
        cap_adata.read_obs()
        obs = cap_adata.obs
        test_compatibility(obs, input_json, validate)

        cap_adata.read_uns()
        cap_adata.uns["cas"] = reformat_json(input_json)

        cap_adata.overwrite()


def test_compatibility(anndata_obs, input_json, validate):
    """
    Tests if CAS and AnnData can be merged.

     Args:
        anndata_obs: The AnnData obs object.
        input_json: The CAS data json object.
        validate: Boolean to determine if validation checks will be performed before writing to the output AnnData file.
    """
    annotations = input_json[ANNOTATIONS]
    obs_index = set(anndata_obs.axes[0].tolist())
    validate_cell_ids(obs_index, annotations, validate)

    labelsets = input_json[LABELSETS]
    matching_obs_keys = get_matching_obs_keys(anndata_obs.columns, labelsets)
    check_labelsets(input_json, anndata_obs, matching_obs_keys, validate)


def check_labelsets(cas_json, input_obs, matching_obs_keys, validate):
    annotations = cas_json[ANNOTATIONS]
    derived_cell_ids = collect_parent_cell_ids(cas_json)

    for ann in annotations:
        if ann[LABELSET] in matching_obs_keys:
            anndata_labelset_cell_ids = (
                input_obs.groupby(ann[LABELSET], observed=False)
                .apply(lambda group: set(group.index), include_groups=False)
                .to_dict()
            )
            for cell_label, cell_list in anndata_labelset_cell_ids.items():
                cell_ids = set(ann.get(CELL_IDS, []))

                if cell_ids and cell_list == cell_ids:
                    handle_matching_labelset(ann, cell_label, input_obs, validate)
                elif cell_list == derived_cell_ids.get(
                    str(ann["cell_set_accession"]), ann.get(CELL_IDS, [])
                ):
                    handle_matching_labelset(ann, cell_label, input_obs, validate)
                elif cell_label == ann[CELL_LABEL]:
                    if cell_list == set(ann.get(CELL_IDS, [])):
                        handle_matching_labelset(ann, cell_label, input_obs, validate)
                    else:
                        handle_non_matching_labelset(
                            ann, input_obs, validate, derived_cell_ids
                        )


def get_matching_obs_keys(obs_keys, cas_labelsets):
    cas_labelset_names = {item[LABELSET_NAME] for item in cas_labelsets}
    matching_obs_keys = cas_labelset_names.intersection(obs_keys)
    return list(matching_obs_keys)


def handle_matching_labelset(ann, cell_label, input_obs, validate):
    # Used for label changes
    if cell_label != ann[CELL_LABEL]:
        logger.warning(
            f"{ann[CELL_LABEL]} cell ids from CAS match with the cell ids in {cell_label} from anndata, "
            "but they have different cell labels."
        )
        if validate:
            logger.error("Validation failed. Exiting.")
            sys.exit(1)
        # add new category to labelset column
        input_obs[ann[LABELSET]] = input_obs[ann[LABELSET]].cat.add_categories(
            ann[CELL_LABEL]
        )
        # Overwrite the labelset value with CAS labelset
        input_obs.loc[ann[CELL_IDS], ann[LABELSET]] = input_obs.loc[
            ann[CELL_IDS], ann[LABELSET]
        ].map({cell_label: ann[CELL_LABEL]})


def handle_non_matching_labelset(ann, input_obs, validate, derived_cell_ids):
    # Used for hierarchy changes
    logger.warning(
        f"{ann[CELL_LABEL]} cell ids from CAS do not match with the cell ids from anndata. "
        "Please update your CAS json."
    )
    if validate:
        logger.error("Validation failed. Exiting.")
        sys.exit(1)

    # Flush the labelset from anndata
    # input_anndata.obs.loc[list(cell_list), cell_label] = ""
    # Add labelset from CAS to anndata
    cell_ids = derived_cell_ids.get(str(ann["cell_set_accession"]), set())
    # Bad split workaround, temporary solution
    # Use Pandas indexing to filter the cell_ids present in obs.index
    valid_cell_ids = input_obs.index.intersection(cell_ids)
    input_obs.loc[valid_cell_ids, ann[LABELSET]] = str(ann[CELL_LABEL])


def validate_cell_ids(anndata_cell_ids, annotations, validate):
    # Collect cell ids from annotations
    cas_cell_ids = {cell_id for ann in annotations for cell_id in ann.get(CELL_IDS, [])}

    # Validate cas -> anndata
    if not cas_cell_ids <= anndata_cell_ids:
        logger.warning("Not all members of cell ids from cas exist in anndata.")
        if validate:
            logger.error("Validation failed. Exiting.")
            sys.exit(1)

    # Validate anndata -> cas
    if not anndata_cell_ids <= cas_cell_ids:
        logger.warning("Not all members of cell ids from anndata exist in cas.")
        if validate:
            logger.error("Validation failed. Exiting.")
            sys.exit(1)
