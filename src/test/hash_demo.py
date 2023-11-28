import os

from cas.file_utils import read_cas_json_file, write_json_file
from cas.accession.hash_accession_manager import HashAccessionManager

CAS_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), "./test_data/siletti/Siletti_all_non_neuronal_cells_with_cids.json")


def main():
    accession_manager = HashAccessionManager(digest_size=5)

    cas = read_cas_json_file(CAS_PATH)
    print(len(cas.annotations))
    for annotation in cas.annotations:
        annotation.cell_set_accession = accession_manager.generate_accession_id(cell_ids=annotation.cell_ids)

    write_json_file(cas, "./hash_demo.json")


if __name__ == "__main__":
    main()
