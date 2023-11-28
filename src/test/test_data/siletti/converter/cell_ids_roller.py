import os

from cas.file_utils import read_cas_json_file, write_json_file, read_anndata_file

CAS_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), "./test_data/siletti/Siletti_all_non_neuronal_cells.json")
ANNDATA_PATH = "/Users/hk9/Downloads/8f2775a8-1a55-4f3c-baf8-23d91b5e6ba3.h5ad"

def roll_cell_ids(cas_json_path, ann_data_path):
    cas = read_cas_json_file(CAS_PATH)
    ann_data = read_anndata_file(ANNDATA_PATH)

    print(len(cas.annotations))
    for annotation in cas.annotations:
        anndata_labelset_cell_ids = (
            ann_data.obs.groupby(ann[LABELSET], observed=False)
            .apply(lambda group: set(group.index))
            .to_dict()
        )

def main():
    roll_cell_ids(CAS_PATH, ANNDATA_PATH)


if __name__ == "__main__":
    main()