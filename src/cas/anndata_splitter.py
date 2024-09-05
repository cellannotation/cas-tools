from pathlib import Path
from typing import Any, Dict, List, Literal, Optional

from anndata import AnnData

from cas.file_utils import read_anndata_file, read_json_file
from cas.utils.conversion_utils import ANNOTATIONS, CELL_IDS, fetch_anndata


def split_anndata_to_file(
    anndata_file_path: Optional[str],
    cas_json_paths: List[str],
    multiple_outputs: bool,
    compression_method: Optional[Literal["gzip", "lzf"]] = "gzip",
):
    """
    Splits an AnnData file into multiple files based on provided CAS JSON files and writes them to disk.

    Args:
        anndata_file_path: Path to the AnnData file.
        cas_json_paths: List of CAS JSON file paths.
        multiple_outputs: If True, outputs multiple files, one for each CAS JSON file; otherwise, outputs a single file.
        compression_method: Compression method utilized in anndata write function. Default is "gzip".

    """
    if not anndata_file_path:
        anndata_file_path = fetch_anndata(
            cas_json_paths[0]
        )  # Assuming all splits are coming from the same CAS JSON.
    adata = read_anndata_file(anndata_file_path)
    cas_list = {
        Path(cas_json).name: read_json_file(cas_json) for cas_json in cas_json_paths
    }
    cas_paths = list(cas_list.keys())

    result = split_anndata(adata, cas_list, multiple_outputs)
    for idx, anndata_item in enumerate(result):
        anndata_item.write_h5ad(
            Path(
                (
                    f"split_{cas_paths[idx].split('.')[0]}.h5ad"
                    if multiple_outputs
                    else "split_anndata.h5ad"
                )
            ),
            compression=compression_method,
        )


def split_anndata(
    adata: AnnData, cas: Dict[str, Dict[str, Any]], multiple_outputs: bool
) -> List[AnnData]:
    """
    Splits an AnnData object into multiple or single AnnData objects based on the provided CAS data.

    Args:
        adata: AnnData object.
        cas: Dictionary representing the CAS data with its file name as keys.
        multiple_outputs: Determines if the output should be multiple AnnData objects or a single one.

    Returns:
        A list of AnnData objects if multiple_outputs is True, otherwise a single AnnData object.

    Raises:
        ValueError: If any required terms do not exist in the CAS data under 'parent_cell_set_name'.
    """

    def get_cell_ids(cas_data: Dict[str, Dict[str, Any]]) -> List[str]:
        """Extract unique cell IDs from CAS data."""
        return list(
            set(
                [
                    cid
                    for cas_obj in cas_data.values()
                    for annotation in cas_obj[ANNOTATIONS]
                    for cid in annotation.get(CELL_IDS, [])
                ]
            )
        )

    if multiple_outputs:
        splitted_anndata_list = []
        for cas_file_name, cas_object in cas.items():
            cell_ids = get_cell_ids({cas_file_name: cas_object})
            mask = adata.obs.index.isin(cell_ids)
            adata_subset = adata[mask, :].to_memory()
            splitted_anndata_list.append(adata_subset)
        return splitted_anndata_list
    else:
        cell_ids = get_cell_ids(cas)
        mask = adata.obs.index.isin(cell_ids)
        adata_subset = adata[mask, :].to_memory()
        return [adata_subset]
