import dataclasses_json
from typing import List, Optional, Any
from dataclasses import dataclass
from dataclasses_json import DataClassJsonMixin


exclude_none_values = True


class EncoderMixin(DataClassJsonMixin):

    dataclass_json_config = dataclasses_json.config(
        # letter_case=dataclasses_json.LetterCase.CAMEL,
        undefined=dataclasses_json.Undefined.EXCLUDE,
        exclude=lambda f: exclude_none_values and f is None
    )["dataclasses_json"]


@dataclass
class AutomatedAnnotation(EncoderMixin):

    algorithm_name: str
    """The name of the algorithm used. It MUST be a string of the algorithm's name."""

    algorithm_version: str
    """The version of the algorithm used (if applicable). It MUST be a string of the algorithm's version, which is 
    typically in the format '[MAJOR].[MINOR]', but other versioning systems are permitted (based on the algorithm's 
    versioning)."""

    algorithm_repo_url: str
    """This field denotes the URL of the version control repository associated with the algorithm used (if applicable). 
    It MUST be a string of a valid URL."""

    reference_location: Optional[str]
    """This field denotes a valid URL of the annotated dataset that was the source of annotated reference data. 
    This MUST be a string of a valid URL. The concept of a 'reference' specifically refers to 'annotation transfer' 
    algorithms, whereby a 'reference' dataset is used to transfer cell annotations to the 'query' dataset."""


@dataclass
class Labelset(EncoderMixin):

    name: str
    """name of annotation key"""

    description: Optional[str] = None
    """Some text describing what types of cell annotation this annotation key is used to record"""

    annotation_method: Optional[str] = None
    """The method used for creating the cell annotations. This MUST be one of the following strings: `'algorithmic'`, 
    `'manual'`, or `'both'` """

    automated_annotation: Optional[AutomatedAnnotation] = None
    """A set of fields for recording the details of the automated annotation algorithm used. (Common 'automated 
    annotation methods' would include PopV, Azimuth, CellTypist, scArches, etc.)"""

    rank: Optional[str] = None
    """A number indicating relative granularity with 0 being the most specific.  Use this where a single dataset has 
    multiple keys that are used consistently to record annotations and different levels of granularity."""


@dataclass
class AnnotationTransfer(EncoderMixin):

    transferred_cell_label: Optional[str]
    """Transferred cell label"""

    source_taxonomy: Optional[str]
    """PURL of source taxonomy."""

    source_node_accession: Optional[str]
    """accession of node that label was transferred from"""

    algorithm_name: Optional[str]
    """The name of the algorithm used."""

    comment: Optional[str]
    """Free text comment on annotation transfer"""


@dataclass
class UserAnnotation(EncoderMixin):
    """User defined custom annotations which are not part of the standard schema."""

    labelset: str
    """The unique name of the set of cell annotations associated with a single file."""

    cell_label: Any
    """This denotes any free-text term which the author uses to label cells."""


@dataclass
class Annotation(EncoderMixin):
    """
    A collection of fields recording a cell type/class/state annotation on some set os cells, supporting evidence and
    provenance. As this is intended as a general schema, compulsory fields are kept to a minimum. However, tools using
    this schema are encouarged to specify a larger set of compulsory fields for publication.
    Note: This schema deliberately allows for additional fields in order to support ad hoc user fields, new formal
    schema extensions and project/tool specific metadata.
    """

    labelset: str
    """The unique name of the set of cell annotations. 
    Each cell within the AnnData/Seurat file MUST be associated with a 'cell_label' value in order for this to be a 
    valid 'cellannotation_setname'."""

    cell_label: str
    """This denotes any free-text term which the author uses to label cells."""

    cell_set_accession: Optional[str] = None
    """An identifier that can be used to consistently refer to the set of cells being annotated, even if the 
    cell_label changes."""

    cell_fullname: Optional[str] = None
    """This MUST be the full-length name for the biological entity listed in `cell_label` by the author. (If the value 
    in `cell_label` is the full-length term, this field will contain the same value.) \nNOTE: any reserved word used in 
    the field 'cell_label' MUST match the value of this field."""

    cell_ontology_term_id: Optional[str] = None
    """This MUST be a term from either the Cell Ontology or from some ontology that extends it by classifying cell 
    types under terms from the Cell Ontology e.g. the Provisional Cell Ontology."""

    cell_ontology_term: Optional[str] = None
    """This MUST be the human-readable name assigned to the value of 'cell_ontology_term_id"""

    cell_ids: Optional[List[str]] = None
    """List of cell barcode sequences/UUIDs used to uniquely identify the cells"""

    rationale: Optional[str] = None
    """The free-text rationale which users provide as justification/evidence for their cell annotations. 
    Researchers are encouraged to use this field to cite relevant publications in-line using standard academic 
    citations of the form `(Zheng et al., 2020)` This human-readable free-text MUST be encoded as a single string. 
    All references cited SHOULD be listed using DOIs under rationale_dois. There MUST be a 2000-character limit."""

    rationale_dois: Optional[List[str]] = None
    """A list of valid publication DOIs cited by the author to support or provide justification/evidence/context for 
    'cell_label'."""

    marker_gene_evidence: Optional[List[str]] = None
    """List of gene names explicitly used as evidence for this cell annotation."""

    synonyms: Optional[List[str]] = None
    """This field denotes any free-text term of a biological entity which the author associates as synonymous with the 
    biological entity listed in the field 'cell_label'."""

    # TODO modified: added
    parent_cell_set_name: Optional[str] = None

    # TODO modified: list -> str
    parent_cell_set_accession: Optional[str] = None
    """A list of accessions of cell sets that subsume this cell set. This can be used to compose hierarchies of 
    annotated cell sets, built from a fixed set of clusters."""

    # TODO modified: added
    user_annotations: Optional[List[UserAnnotation]] = None

    # TODO modified: moved from CTA to Annotation class
    transferred_annotations: Optional[AnnotationTransfer] = None

    def add_user_annotation(self, user_annotation_set, user_annotation_label):
        """
        Adds a user defined annotation which is not supported by the standard schema.
        :param user_annotation_set: name of the user annotation set
        :param user_annotation_label: label of the user annotation set
        """
        if not self.user_annotations:
            self.user_annotations = list()
        self.user_annotations.append(UserAnnotation(user_annotation_set, user_annotation_label))


@dataclass
class CellTypeAnnotation(EncoderMixin):

    # data_url: str
    # annotation_objects: List[Annotation]
    # taxonomy: TaxonomyMetadata = None

    author_name: str
    """This MUST be a string in the format `[FIRST NAME] [LAST NAME]`"""

    annotations: List[Annotation]
    """A collection of fields recording a cell type/class/state annotation on some set os cells, supporting evidence 
    and provenance. As this is intended as a general schema, compulsory fields are kept to a minimum. However, tools 
    using this schema are encouarged to specify a larger set of compulsory fields for publication."""

    labelsets: Optional[List[Labelset]] = None

    author_contact: Optional[str] = None
    """This MUST be a valid email address of the author"""

    orcid: Optional[str] = None
    """This MUST be a valid ORCID for the author"""

    cellannotation_schema_version: Optional[str] = None
    """The schema version, the cell annotation open standard. Current version MUST follow 0.1.0
    This versioning MUST follow the format `'[MAJOR].[MINOR].[PATCH]'` as defined by Semantic Versioning 2.0.0, 
    https://semver.org/"""

    cellannotation_version: Optional[str] = None
    """The version for all cell annotations published (per dataset). This MUST be a string. The recommended versioning 
    format is `'[MAJOR].[MINOR].[PATCH]'` as defined by Semantic Versioning 2.0.0, https://semver.org/"""

    cellannotation_url: Optional[str] = None
    """A persistent URL of all cell annotations published (per dataset)."""

    def add_annotation_object(self, obj):
        """
        Adds given object to annotation objects list
        :param obj: Annotation object to add
        """
        self.annotations.append(obj)

    def set_exclude_none_values(self, value):
        global exclude_none_values
        exclude_none_values = value
