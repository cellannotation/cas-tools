import abc


class BaseAccessionManager(metaclass=abc.ABCMeta):
    """
    Abstract Accession ID generator.
    """

    @abc.abstractmethod
    def generate_accession_id(
        self, id_recommendation: str = None, labelset: str = None
    ) -> str:
        """
        Generates an auto-increment based accession id. If the recommended accession_id is available, uses it.
        Params:
            id_recommendation: accession id recommendation. Function uses this id if it is available,
            provides an auto-incremented id otherwise.
            labelset: Labelset name. If provided, uses it as a prefix to the accession id.
        Return: accession_id
        """
        pass
