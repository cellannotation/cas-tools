import json
import sys

LABELSET_NAME = "name"

LABELSET = "labelset"

ANNOTATIONS = "annotations"

CELL_IDS = "cell_ids"

CELL_LABEL = "cell_label"


def check_labelsets(annotations, input_anndata, matching_obs_keys, validate):
    for ann in annotations:
        if ann[LABELSET] in matching_obs_keys:
            anndata_labelset_cell_ids = (
                input_anndata.obs.groupby(ann[LABELSET], observed=False)
                .apply(lambda group: set(group.index))
                .to_dict()
            )
            for cell_label, cell_list in anndata_labelset_cell_ids.items():
                if cell_list == set(ann[CELL_IDS]):
                    handle_matching_labelset(ann, cell_label, input_anndata, validate)
                elif cell_label == ann[CELL_LABEL]:
                    handle_non_matching_labelset(
                        ann, cell_label, cell_list, input_anndata, validate
                    )


def get_cas_annotations(input_json):
    return input_json[ANNOTATIONS]


def get_matching_obs_keys(input_anndata, input_json):
    cas_labelset_names = [item[LABELSET_NAME] for item in input_json[LABELSET]]
    obs_keys = input_anndata.obs_keys()
    matching_obs_keys = list(set(obs_keys).intersection(cas_labelset_names))
    return matching_obs_keys


def handle_matching_labelset(ann, cell_label, input_anndata, validate):
    if cell_label != ann[CELL_LABEL]:
        print(
            f"{ann[CELL_LABEL]} cell ids from CAS match with the cell ids in {cell_label} from anndata. "
            "But they have different cell label."
        )
        if validate:
            sys.exit()
        # add new category to labelset column
        input_anndata.obs[ann[LABELSET]] = input_anndata.obs[
            ann[LABELSET]
        ].cat.add_categories(ann[CELL_LABEL])
        # Overwrite the labelset value with CAS labelset
        input_anndata.obs.loc[ann[CELL_IDS], ann[LABELSET]] = input_anndata.obs.loc[
            ann[CELL_IDS], ann[LABELSET]
        ].map({cell_label: ann[CELL_LABEL]})


def handle_non_matching_labelset(ann, cell_label, cell_list, input_anndata, validate):
    print(
        f"{ann[CELL_LABEL]} cell ids from CAS do not match with the cell ids from anndata. "
        "Please update your CAS json."
    )
    if validate:
        sys.exit()
    # Flush the labelset from anndata
    input_anndata.obs.loc[cell_list, cell_label] = None
    # Add labelset from CAS to anndata
    input_anndata.obs.loc[ann[CELL_IDS], ann[LABELSET]] = ann[CELL_LABEL]


def save_cas_to_uns(input_anndata, input_json):
    # drop cell_ids
    json_without_cell_ids = {
        "author_name": input_json["author_name"],
        "labelset": input_json["labelset"],
        "annotations": [
            {key: value for key, value in annotation.items() if key != "cell_ids"}
            for annotation in input_json["annotations"]
        ],
    }
    input_anndata.uns.update({"cas": json.dumps(json_without_cell_ids)})


def validate_cell_ids(input_anndata, annotations, validate):
    # check cell ids
    cas_cell_ids = set()
    for ann in annotations:
        cas_cell_ids.update(ann.get(CELL_IDS, []))
    anndata_cell_ids = set(input_anndata.obs.index)
    # cas -> anndata
    if not cas_cell_ids.issubset(anndata_cell_ids):
        print("Not all members of cell ids from cas exist in anndata.")
        if validate:
            sys.exit()
    # anndata -> cas
    if not anndata_cell_ids.issubset(cas_cell_ids):
        print("Not all members of cell ids from anndata exist in cas.")
        if validate:
            sys.exit()
    return annotations


def write_anndata(input_anndata, output_file_path):
    # Close the AnnData file to prevent blocking
    input_anndata.file.close()
    input_anndata.write(output_file_path)
