import hashlib

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

    def generate_accession_id(self, id_recommendation: str = None, cell_ids: List = None) -> str:
        """
        Generates a Blake2b hashing algorithm based hash for the given cell IDs.
        Params:
            id_recommendation: this value is ignored in this manager
            cell_ids: Cell IDs list. Algorithm sorts cell ids internally.
        Return: accession_id
        """
        if not cell_ids:
            raise Exception("Cell IDs list is empty.")

        blake_hasher = hashlib.blake2b(str.encode(" ".join(sorted(cell_ids))), digest_size=self.digest_size)
        accession_id = blake_hasher.hexdigest()

        if accession_id in self.accession_ids:
            print(accession_id)
            # raise Exception("Hash ID conflict occurred: " + accession_id)
        else:
            self.accession_ids.append(accession_id)
        return accession_id

