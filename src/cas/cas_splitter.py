from collections import defaultdict
from typing import Any, Dict, List, Union

from cas.file_utils import read_json_file, write_dict_to_json_file

ANNOTATIONS = "annotations"
CELL_LABEL = "cell_label"
CELL_SET_ACCESSION = "cell_set_accession"
LABELSET = "labelset"
LABELSETS = "labelsets"
LABELSETS_NAME = "name"
PARENT_CELL_SET_ACCESSION = "parent_cell_set_accession"


def split_cas_to_file(
    cas_json_path: str, split_terms: Union[List[str], str], multiple_outputs: bool
):
    """
    Splits a CAS JSON file into files based on provided terms, and writes them to disk.

    Args:
        cas_json_path: Path to the CAS JSON file.
        split_terms: Terms used to determine how to split the CAS file; can be a string or a list of strings.
        multiple_outputs: If True, outputs multiple files, one for each split term; otherwise, outputs a single file.

    """
    cas = read_json_file(cas_json_path)
    result = split_cas(cas, split_terms, multiple_outputs)
    if isinstance(result, dict):
        result = [result]
    for idx, cas_item in enumerate(result):
        write_dict_to_json_file(
            (
                f"cas_{split_terms[idx].replace(':', '_')}.json"
                if multiple_outputs
                else "split_cas.json"
            ),
            cas_item,
        )


def split_cas(
    cas: Dict[str, Any], split_terms: Union[List[str], str], multiple_outputs: bool
) -> Union[List[Dict[str, Any]], Dict[str, Any]]:
    """
    Splits a CAS dictionary into multiple or single dictionary based on split terms.

    Args:
        cas: Dictionary representing the CAS data.
        split_terms: Terms used to filter and split the CAS data; can be a string or a list of strings.
        multiple_outputs: Determines if the output should be multiple dictionaries or a single dictionary.

    Returns:
        A list of dictionaries if multiple_outputs is True, otherwise a single dictionary.

    Raises:
        ValueError: If any split_terms do not exist in the CAS data under 'parent_cell_set_name'.
    """
    cell_dict = {
        annotation[CELL_SET_ACCESSION]: annotation[PARENT_CELL_SET_ACCESSION]
        for annotation in cas[ANNOTATIONS]
        if (PARENT_CELL_SET_ACCESSION in annotation)
    }
    parent_cell_dict = defaultdict(list)
    for child_cell, parent_cell in cell_dict.items():
        parent_cell_dict[parent_cell].append(child_cell)
    if isinstance(split_terms, str):
        split_terms = [split_terms]
    keys_and_values = list(parent_cell_dict.keys()) + [item for sublist in parent_cell_dict.values() if isinstance(sublist, list) for item in sublist]
    missing_terms = [term for term in split_terms if term not in keys_and_values]
    if missing_terms:
        raise ValueError(
            f"{', '.join(missing_terms)} do not exist in CAS as 'cell_set_name'"
        )
    if multiple_outputs:
        splitted_cas_list = []
        for term in split_terms:
            label_to_copy_list = get_split_terms(parent_cell_dict, term)
            splitted_cas_list.append(
                filter_and_copy_cas_entries(cas, label_to_copy_list)
            )
        return splitted_cas_list
    else:
        label_to_copy_list = get_split_terms(parent_cell_dict, split_terms)
        return filter_and_copy_cas_entries(cas, label_to_copy_list)


def filter_and_copy_cas_entries(
    cas: Dict[str, Any], label_to_copy_list: List[str]
) -> Dict[str, Any]:
    """
    Copies entries from the CAS based on a list of labels to copy.

    Args:
        cas: Dictionary representing the original CAS data.
        label_to_copy_list: List of labels indicating which entries to copy.

    Returns:
        A dictionary with filtered CAS entries.
    """
    output_dict = dict(cas)
    output_dict[ANNOTATIONS] = []
    output_dict[LABELSETS] = []
    labelset_dict = set()
    for annotation in cas[ANNOTATIONS]:
        if annotation[CELL_SET_ACCESSION] in label_to_copy_list:
            output_dict[ANNOTATIONS].append(annotation)
            labelset_dict.add(annotation[LABELSET])
    for labelset in cas[LABELSETS]:
        if labelset[LABELSETS_NAME] in labelset_dict:
            output_dict[LABELSETS].append(labelset)
    return output_dict


def get_split_terms(
    parent_dict: Dict[str, List[str]], split_terms: Union[List[str], str]
) -> List[str]:
    """
    Resolves split terms into a comprehensive list of terms based on a parent-child relationship dictionary.

    Args:
        parent_dict: Dictionary mapping parent terms to lists of child terms.
        split_terms: Initial terms to resolve, can be a string or a list of strings.

    Returns:
        A list of all terms, resolved from the parent_dict.
    """
    if isinstance(split_terms, str):
        split_terms = [split_terms]
    result = []
    stack = list(split_terms)
    while stack:
        term = stack.pop()
        if term in parent_dict:
            stack.extend(parent_dict[term])
        result.append(term)
    return result
