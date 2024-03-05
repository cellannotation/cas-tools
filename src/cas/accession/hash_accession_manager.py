import hashlib
import string
from typing import List

from cas.accession.base_accession_manager import BaseAccessionManager


class HashAccessionManager(BaseAccessionManager):
    def __init__(self, accession_prefix=None, digest_size=5):
        """
        Initializer.
        Params:
            accession_prefix: accession_id prefix
            digest_size: output hash size
        """
        self.accession_prefix = accession_prefix
        self.digest_size = digest_size
        self.accession_ids = list()

    def generate_accession_id(
        self, id_recommendation: str = None, cell_ids: List = None
    ) -> str:
        """
        Generates a Blake2b hashing algorithm based hash for the given cell IDs.
        Params:
            id_recommendation: pre-calculated hash accession recommendation. Returns this value if recommendation is a
            valid accession id.
            cell_ids: Cell IDs list. Algorithm sorts cell ids internally.
        Return: accession_id
        """
        if is_hash_accession(id_recommendation):
            return id_recommendation

        if not cell_ids:
            raise Exception("Cell IDs list is empty.")

        blake_hasher = hashlib.blake2b(
            str.encode(" ".join(sorted(cell_ids))), digest_size=self.digest_size
        )
        accession_id = blake_hasher.hexdigest()

        if accession_id in self.accession_ids:
            print(accession_id)
            # raise Exception("Hash ID conflict occurred: " + accession_id)
        else:
            self.accession_ids.append(accession_id)
        return accession_id


def is_hash_accession(accession_id: str):
    """
    Checks if the given accession is a valid hash accession. Hash accessions are 10 char long and only has hexdigits
    Args:
        accession_id: accession to check

    Returns: True if value is a valid hash accession id, false otherwise.

    """
    return accession_id and len(accession_id) == 10 and all(c in string.hexdigits for c in accession_id)
