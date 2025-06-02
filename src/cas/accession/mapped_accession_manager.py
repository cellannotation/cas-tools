from cas.accession.base_accession_manager import BaseAccessionManager


class MappedAccessionManager(BaseAccessionManager):
    """
    Predefined Accession ID generator based on a mapping of cell set names to accession IDs.
    This accession manager is used when the accession IDs are already defined and mapped to specific cell sets.
    """

    def __init__(self, accession_map):
        """
        Initializer.
        Params:
            accession_map: map of cell set names to their corresponding accession IDs.
            (To enable usage of same names accross different labelsets, key is identified as labelset:cell_label).
        """
        self.accession_map = accession_map

    def generate_accession_id(
        self, id_recommendation: str = None, labelset: str = None, cellset_name: str = None, **kwargs
    ) -> str:
        """
        Generates an auto-increment based accession id. If the recommended accession_id is available, uses it.
        Params:
            id_recommendation: this parameter is not utilized in this implementation.
            labelset: this parameter is not utilized in this implementation.
            cellset_name: Name of the cell set for which the accession ID is being generated.
        Return: accession_id
        """
        if labelset + ':' + cellset_name in self.accession_map:
            return self.accession_map[labelset + ':' +cellset_name]
        else:
            raise ValueError(f"Cell set name '{labelset}:{cellset_name}' not found in the accession map.")