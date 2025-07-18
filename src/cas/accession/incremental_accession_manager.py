from cas.accession.base_accession_manager import BaseAccessionManager


class IncrementalAccessionManager(BaseAccessionManager):
    """
    Numerically incremental Accession ID generator.
    """

    def __init__(self, accession_prefix=None, last_accession_id=0):
        """
        Initializer.
        Params:
            accession_prefix: accession_id prefix
        """
        if accession_prefix is None:
            accession_prefix = ""
        self.accession_prefix = accession_prefix
        self.last_accession_id = last_accession_id
        self.accession_ids = list()

    def generate_accession_id(
        self, id_recommendation: str = None, labelset: str = None, cellset_name: str = None
    ) -> str:
        """
        Generates an auto-increment based accession id. If the recommended accession_id is available, uses it.
        Params:
            id_recommendation: accession id recommendation. Function uses this id if it is available,
            provides an auto-incremented id otherwise.
            labelset: this parameter is not utilized in this implementation.
            cellset_name: this parameter is not utilized in this implementation.
        Return: accession_id
        """
        if id_recommendation:
            id_recommendation = id_recommendation.replace(self.accession_prefix, "")
            if id_recommendation.startswith("_"):
                id_recommendation = id_recommendation[1:]

        if (
            id_recommendation
            and id_recommendation not in self.accession_ids
            and id_recommendation.isdigit()
            # and int(id_recommendation) > self.last_accession_id
            and int(id_recommendation) not in self.accession_ids
        ):
            accession_id = id_recommendation
            if self.last_accession_id < int(id_recommendation):
                self.last_accession_id = int(id_recommendation)
        elif id_recommendation and not id_recommendation.isdigit():
            # non-numeric accession id
            accession_id = id_recommendation
        else:
            id_candidate = self.last_accession_id + 1
            while str(id_candidate) in self.accession_ids:
                id_candidate += 1
            accession_id = str(id_candidate)
            self.last_accession_id = id_candidate

        self.accession_ids.append(accession_id)
        if self.accession_prefix and not accession_id.startswith(self.accession_prefix):
            accession_id = self.accession_prefix + accession_id

        return accession_id
