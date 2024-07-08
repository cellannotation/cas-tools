from dataclasses import asdict, dataclass, field
from datetime import datetime
from typing import Any, List, Optional

import dataclasses_json
import pandas as pd
from dataclasses_json import DataClassJsonMixin, config
from dateutil.parser import isoparse
from marshmallow import fields

from cas.reports import get_all_annotations

exclude_none_values = True


class EncoderMixin(DataClassJsonMixin):
    dataclass_json_config = dataclasses_json.config(
        # letter_case=dataclasses_json.LetterCase.CAMEL,
        undefined=dataclasses_json.Undefined.EXCLUDE,
        exclude=lambda f: exclude_none_values and f is None,
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
    transferred_cell_label: Optional[str] = None
    """Transferred cell label"""

    source_taxonomy: Optional[str] = None
    """PURL of source taxonomy."""

    source_node_accession: Optional[str] = None
    """accession of node that label was transferred from"""

    algorithm_name: Optional[str] = None
    """The name of the algorithm used."""

    comment: Optional[str] = None
    """Free text comment on annotation transfer"""


@dataclass
class Review(EncoderMixin):
    """Annotation review."""

    datestamp: Optional[datetime] = field(
        metadata=config(
            encoder=lambda x: (
                x.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z" if x is not None else None
            ),
            decoder=lambda x: isoparse(x) if x is not None else None,
            mm_field=fields.DateTime(format="iso"),
        ),
        default=None,
    )
    """Time and date review was last edited."""

    reviewer: Optional[str] = None
    """Review Author."""

    review: Optional[str] = None
    """Reviewer's verdict on the annotation.  Must be 'Agree' or 'Disagree'."""

    explanation: Optional[str] = None
    """Free-text review of annotation. This is required if the verdict is disagree and should include reasons for disagreement."""


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

    # TODO modified: added (exclude from json serialisation)
    parent_cell_set_name: Optional[str] = field(default=None, metadata=config(exclude=lambda x: True))

    parent_cell_set_accession: Optional[str] = None
    """A list of accessions of cell sets that subsume this cell set. This can be used to compose hierarchies of 
    annotated cell sets, built from a fixed set of clusters."""

    author_annotation_fields: Optional[dict] = None
    """"A dictionary of author defined key value pairs annotating the cell set. The names and aims of these fields MUST 
    not clash with official annotation fields."""

    transferred_annotations: Optional[List[AnnotationTransfer]] = None

    reviews: Optional[List[Review]] = None

    def add_user_annotation(self, user_annotation_set, user_annotation_label):
        """
        Adds a user defined annotation which is not supported by the standard schema.
        :param user_annotation_set: name of the user annotation set
        :param user_annotation_label: label of the user annotation set
        """
        if not self.author_annotation_fields:
            self.author_annotation_fields = dict()
        self.author_annotation_fields[user_annotation_set] = user_annotation_label


@dataclass
class CellTypeAnnotation(EncoderMixin):

    author_name: str
    """This MUST be a string in the format `[FIRST NAME] [LAST NAME]`"""

    annotations: List[Annotation]
    """A collection of fields recording a cell type/class/state annotation on some set os cells, supporting evidence 
    and provenance. As this is intended as a general schema, compulsory fields are kept to a minimum. However, tools 
    using this schema are encouarged to specify a larger set of compulsory fields for publication."""

    title: str
    """The title of the dataset. This MUST be less than or equal to 200 characters. e.g. 'Human retina cell atlas - 
    retinal ganglion cells'."""

    description: Optional[str] = None
    """The description of the dataset. e.g. 'A total of 15 retinal ganglion cell clusters were identified from over 99K 
    retinal ganglion cell nuclei in the current atlas. Utilizing previous characterized markers from macaque, 5 clusters
     can be annotated.'"""

    matrix_file_id: Optional[str] = None
    """A resolvable ID for a cell by gene matrix file in the form namespace:accession, e.g. 
    CellXGene_dataset:8e10f1c4-8e98-41e5-b65f-8cd89a887122. Please see https://github.com/cellannotation/cell-annotation
    -schema/registry/registry.json for supported namespaces."""

    labelsets: Optional[List[Labelset]] = None
    """A list of labelsets that are used in the annotations."""

    author_contact: Optional[str] = None
    """This MUST be a valid email address of the author"""

    orcid: Optional[str] = None
    """This MUST be a valid ORCID for the author"""

    cellannotation_schema_version: Optional[str] = None
    """The schema version, the cell annotation open standard. Current version MUST follow 0.1.0
    This versioning MUST follow the format `'[MAJOR].[MINOR].[PATCH]'` as defined by Semantic Versioning 2.0.0, 
    https://semver.org/"""

    cellannotation_timestamp: Optional[str] = None
    """The timestamp of all cell annotations published (per dataset). This MUST be a string in the format 
    '%yyyy-%mm-%dd %hh:%mm:%ss'."""

    cellannotation_version: Optional[str] = None
    """The version for all cell annotations published (per dataset). This MUST be a string. The recommended versioning 
    format is `'[MAJOR].[MINOR].[PATCH]'` as defined by Semantic Versioning 2.0.0, https://semver.org/"""

    cellannotation_url: Optional[str] = None
    """A persistent URL of all cell annotations published (per dataset)."""

    author_list: Optional[List[str]] = None
    """This field stores a list of users who are included in the project as collaborators, regardless of their specific 
    role. An example list; John Smith|Cody Miller|Sarah Jones."""

    def add_annotation_object(self, obj):
        """
        Adds given object to annotation objects list
        :param obj: Annotation object to add
        """
        self.annotations.append(obj)

    def set_exclude_none_values(self, value):
        global exclude_none_values
        exclude_none_values = value

    def get_all_annotations(
        self, show_cell_ids: bool = False, labels: list = None
    ) -> pd.DataFrame:
        """
        Lists all annotations.

        Args:
            show_cell_ids: identifies if result have 'cell_ids' column. Default value is false
            labels: list of key(labelset), value(cell_label) pairs to filter annotations
        Returns:
            Annotations data frame
        """
        return get_all_annotations(
            asdict(self), show_cell_ids=show_cell_ids, labels=labels
        )
