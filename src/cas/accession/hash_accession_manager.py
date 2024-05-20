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
        self, id_recommendation: str = None, cell_ids: List = None, labelset: str = None
    ) -> str:
        """
        Generates a Blake2b hashing algorithm based hash for the given cell IDs.
        Params:
            id_recommendation: pre-calculated hash accession recommendation. Returns this value if recommendation is a
            valid accession id.
            cell_ids: Cell IDs list. Algorithm sorts cell ids internally.
            labelset: Labelset name. If provided, uses it as a prefix to the accession id.
        Return: accession_id
        """
        if id_recommendation and labelset and ":" not in id_recommendation:
            id_recommendation = labelset + ":" + id_recommendation
        if is_hash_accession(id_recommendation):
            return id_recommendation

        if not cell_ids:
            raise Exception("Cell IDs list is empty.")

        blake_hasher = hashlib.blake2b(
            str.encode(" ".join(sorted(cell_ids))), digest_size=self.digest_size
        )
        accession_id = blake_hasher.hexdigest()
        if labelset:
            accession_id = labelset + ":" + accession_id

        if accession_id in self.accession_ids:
            print("ERROR: Hash ID conflict occurred: " + accession_id)
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
    if not accession_id:
        return False
    hash_part = accession_id.split(":")[-1]
    return len(hash_part) == 10 and all(c in string.hexdigits for c in hash_part)
