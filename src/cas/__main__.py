import argparse
import os
import pathlib
import warnings

from cas.abc_cas_converter import abc2cas, cas2abc
from cas.add_author_annotations import add_author_annotations_from_file
from cas.anndata_conversion import merge
from cas.anndata_to_cas import anndata2cas
from cas.flatten_data_to_anndata import flatten
from cas.populate_cell_ids import populate_cell_ids
from cas.spreadsheet_to_cas import spreadsheet2cas
from cas.validate import validate as schema_validate
from cas.cas_to_rdf import export_to_rdf
from cas.cas_splitter import split_cas_to_file

warnings.filterwarnings("always")


def main():
    parser = argparse.ArgumentParser(
        prog="cas", description="Cell Type Annotation Tools cli interface."
    )
    subparsers = parser.add_subparsers(help="Available ctat actions", dest="action")

    create_merge_operation_parser(subparsers)
    create_flatten_operation_parser(subparsers)
    create_spreadsheet2cas_operation_parser(subparsers)
    create_anndata2cas_operation_parser(subparsers)
    create_abc2cas_operation_parser(subparsers)
    create_cas2abc_operation_parser(subparsers)
    create_populate_cells_operation_parser(subparsers)
    create_schema_validation_operation_parser(subparsers)
    create_cas2rdf_operation_parser(subparsers)
    create_add_author_annotations_parser(subparsers)
    create_split_cas_parser(subparsers)

    args = parser.parse_args()

    if args.action == "merge":
        json_file_path = args.json
        anndata_file_path = args.anndata
        output_file_path = args.output
        validate = args.validate

        if anndata_file_path == output_file_path:
            raise ValueError("--anndata and --output cannot be the same")

        merge(json_file_path, anndata_file_path, validate, output_file_path)
    elif args.action == "flatten":
        args = parser.parse_args()
        json_file_path = args.json
        anndata_file_path = args.anndata
        output_file_path = args.output

        if output_file_path and os.path.abspath(anndata_file_path) == os.path.abspath(
            output_file_path
        ):
            raise ValueError("--anndata and --output cannot be the same")
        flatten(json_file_path, anndata_file_path, output_file_path)
    elif args.action == "spreadsheet2cas":
        args = parser.parse_args()
        spreadsheet_file_path = args.spreadsheet
        sheet_name = args.sheet
        anndata_file_path = args.anndata
        labelsets = args.labelsets
        schema_name = args.schema
        output_file_path = args.output

        spreadsheet2cas(
            spreadsheet_file_path,
            sheet_name,
            anndata_file_path,
            labelsets,
            schema_name,
            output_file_path,
        )
    elif args.action == "anndata2cas":
        args = parser.parse_args()
        anndata_file_path = args.anndata
        labelsets = args.labelsets
        output_file_path = args.output
        include_hierarchy = args.hierarchy

        anndata2cas(anndata_file_path, labelsets, output_file_path, include_hierarchy)
    elif args.action == "abc2cas":
        args = parser.parse_args()
        cat_set_file_path = args.catset
        cat_file_path = args.cat
        output_file_path = args.output

        abc2cas(cat_set_file_path, cat_file_path, output_file_path)
    # elif args.action == "cas2abc":
    #     args = parser.parse_args()
    #     json_file_path = args.json
    #     cat_set_file_path = args.catset
    #     cat_file_path = args.cat
    #
    #     cas2abc(json_file_path, cat_set_file_path, cat_file_path)
    elif args.action == "populate_cells":
        args = parser.parse_args()
        json_file_path = args.json
        anndata_file_path = args.anndata
        labelsets = None
        if "labelsets" in args and args.labelsets:
            labelsets = [item.strip() for item in str(args.labelsets).split(",")]
        populate_cell_ids(json_file_path, anndata_file_path, labelsets)
    elif args.action == "validate":
        args = parser.parse_args()
        schema = args.schema
        data_file_path = args.data
        schema_validate(schema, data_file_path)
    elif args.action == "cas2rdf":
        args = parser.parse_args()
        export_to_rdf(
            cas_schema=args.schema,
            data=args.data,
            ontology_namespace=args.ontology_ns,
            ontology_iri=args.ontology_iri,
            labelsets=args.labelsets,
            output_path=args.out,
            validate=not args.skip_validate,
            include_cells=not args.exclude_cells,
        )
    elif args.action == "add_author_annotations":
        args = parser.parse_args()
        # Determine the column(s) to use for joining based on the specified arguments
        if args.join_on:
            join_column = args.join_on
        elif args.join_on_cellset_ids:
            join_column = 'cell_set_accession'
        elif args.join_on_labelset_label:
            join_column = ['labelset', 'cell_label']
        else:
            raise ValueError("No valid join column specified.")

        add_author_annotations_from_file(args.cas_json, args.csv, join_column, args.columns, args.output)
    elif args.action == "split_cas":
        cas_json_path = args.cas_json
        split_terms = args.split_on
        multiple_outputs = args.multiple_outputs

        split_cas_to_file(cas_json_path, split_terms, multiple_outputs)


def create_merge_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --json      : Path to the CAS JSON schema file.
    --anndata   : Path to the AnnData file. Ideally, the location will be specified by a resolvable path in the CAS file.
    --validate  : Perform validation checks before writing to the output AnnData file.
    --output    : Output AnnData file name (default: output.h5ad).

    Usage Example:
    --------------
    cd src
    python -m cas merge --json path/to/CAS_schema.json --anndata path/to/input_anndata.h5ad --validate --output path/to/output.h5ad
    """
    parser_merge = subparsers.add_parser(
        "merge",
        description="The CAS and AnnData merge parser",
        help="Test if CAS can be merged to the AnnData and merges if possible.",
    )
    parser_merge.add_argument(
        "--json", required=True, help="Path to the CAS JSON schema file."
    )
    parser_merge.add_argument(
        "--anndata", required=True, help="Path to the AnnData file."
    )
    # TODO find a better argument name and Help message.
    parser_merge.add_argument(
        "-v",
        "--validate",
        action="store_true",
        help="Perform validation checks before writing to the output AnnData file.",
    )
    parser_merge.add_argument(
        "--output",
        default="output.h5ad",
        help="Output AnnData file name (default: output.h5ad).",
    )
    parser_merge.set_defaults(validate=False)


def create_flatten_operation_parser(subparsers):
    """
        Command-line Arguments:
    -----------------------
    --json      : Path to the CAS JSON schema file.
    --anndata   : Path to the AnnData file. Ideally, the location will be specified by a resolvable path in the CAS file.
    --output    : Optional output AnnData file name. If provided a new flatten anndata file will be created,
    otherwise the inputted anndata file will be updated with the flatten data.

    Usage Example:
    --------------
    cd src
    python -m cas flatten --json path/to/json_file.json --anndata path/to/anndata_file.h5ad --output path/to/output_file.h5ad
    """
    parser_flatten = subparsers.add_parser(
        "flatten",
        description="Flattens all content of CAS annotations to an AnnData file.",
        help="Flattens all content of CAS annotations to obs key:value pairs. "
        "Flattens all other content to key_value pairs in uns.",
    )

    parser_flatten.add_argument("--json", required=True, help="Input JSON file path")
    parser_flatten.add_argument(
        "--anndata", required=True, help="Input AnnData file path"
    )
    parser_flatten.add_argument(
        "--output",
        required=False,
        help="Output AnnData file name.",
    )
    parser_flatten.set_defaults(validate=False)


def create_spreadsheet2cas_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --spreadsheet   : Path to the spreadsheet file.
    --sheet         : Target sheet name in the spreadsheet.
    --anndata       : Path to the AnnData file. If not provided, AnnData will be downloaded using CxG LINK in
                    spreadsheet.
    --labelsets     : List of names of observation (obs) fields used to record author cell type names,
        which determine the rank of labelsets in a spreadsheet. If not provided, ranks will be determined based on
        the order of the fields specified in the CELL LABELSET NAME column.
    --schema        : Name of cell annotation schema used to in spreadsheet. It can be one of 'base', 'bican' or 'cap'.
    --output        : Output CAS file name (default: output.json).

    Usage Example:
    --------------
    cd src
    python -m cas spreadsheet2cas --spreadsheet path/to/spreadsheet.xlsx --sheet sheet_name --labelsets item1 item2 item3 --output path/to/output_file.json
    """

    parser_spreadsheet2cas = subparsers.add_parser(
        "spreadsheet2cas",
        description="Converts a spreadsheet to CAS JSON.",
        help="Converts a spreadsheet to Cell Annotation Schema (CAS) JSON.",
        usage="cas spreadsheet2cas --spreadsheet path/to/spreadsheet.xlsx --sheet sheet_name --labelsets item1 item2 item3 --output path/to/output_file.json",
    )

    parser_spreadsheet2cas.add_argument(
        "--spreadsheet", required=True, help="Path to the spreadsheet file."
    )
    parser_spreadsheet2cas.add_argument(
        "--sheet", required=False, help="Target sheet name in the spreadsheet."
    )
    parser_spreadsheet2cas.add_argument(
        "--anndata",
        default=None,
        help="Path to the AnnData file. If not provided, AnnData will be downloaded using CxG LINK in spreadsheet.",
    )
    parser_spreadsheet2cas.add_argument(
        "--labelsets",
        default=None,
        nargs="+",
        help="List of names of observation (obs) fields used to record author cell type names, which determine the "
        "rank of labelsets in a spreadsheet. If not provided, ranks will be determined based on the order of the"
        " fields specified in the CELL LABELSET NAME column.",
    )
    parser_spreadsheet2cas.add_argument(
        "--schema",
        default="cap",
        help="Name of cell annotation schema used to in spreadsheet. It can be one of 'base', 'bican' or 'cap'.",
        choices=["base", "bican", "cap"],
    )
    parser_spreadsheet2cas.add_argument(
        "--output",
        default="output.json",
        help="Output CAS file name (default: output.json).",
    )


def create_anndata2cas_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --anndata       : Path to the AnnData file.
    --labelsets     : List of labelsets, which are names of observation (obs) fields used to record author cell type
    names. The labelsets should be provided in order, starting from rank 0 (leaf nodes) and ascending to higher ranks.

    --output        : Output CAS file name (default: output.json).
    --hierarchy     : Flag indicating whether to include hierarchy in the output.


    Usage Example:
    --------------
    cd src
    python -m cas anndata2cas --anndata path/to/anndata.h5ad --labelsets item1 item2 item3 --output
    path/to/output_file.json
    """
    parser_anndata2cas = subparsers.add_parser(
        "anndata2cas",
        description="Converts an anndata to CAS JSON.",
        help="Converts an anndata to Cell Annotation Schema (CAS) JSON.",
        usage="cas anndata2cas --anndata path/to/anndata.h5ad --labelsets item1 item2 item3 --output path/to/output_file.json",
    )

    parser_anndata2cas.add_argument(
        "--anndata",
        required=True,
        help="Path to the AnnData file.",
    )
    parser_anndata2cas.add_argument(
        "--labelsets",
        required=True,
        nargs="+",
        help="List of labelsets, which are names of observation (obs) fields used to record author cell type names. "
        "The labelsets should be provided in order, starting from rank 0 (leaf nodes) and ascending to higher ranks.",
    )
    parser_anndata2cas.add_argument(
        "--output",
        default="output.json",
        help="Output CAS file name (default: output.json).",
    )
    parser_anndata2cas.add_argument(
        "--hierarchy",
        action="store_true",
        help="Include hierarchy in the output.",
    )


def create_abc2cas_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --catset       : Path to the Cluster Annotation Term Set file.
    --cat           : Path to the Cluster Annotation Term file.
    --output        : Output CAS file name (default: output.json).


    Usage Example:
    --------------
    cd src

    python -m cas abc2cas --catset path/to/cluster_annotation_term_set.csv --cat path/to/cluster_annotation_term.csv
    --output path/to/output_file.json
    """
    parser_abc2cas = subparsers.add_parser(
        "abc2cas",
        description="Converts given ABC cluster_annotation files to CAS JSON.",
        help="Converts given ABC cluster_annotation files to Cell Annotation Schema (CAS) JSON.",
    )

    parser_abc2cas.add_argument(
        "--catset",
        required=True,
        help="Path to the Cluster Annotation Term Set file.",
    )
    parser_abc2cas.add_argument(
        "--cat",
        required=True,
        help="Path to the Cluster Annotation Term file.",
    )
    parser_abc2cas.add_argument(
        "--output",
        default="output.json",
        help="Output CAS file name (default: output.json).",
    )


def create_cas2abc_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --json         : Path to the CAS JSON schema file.
    --catset       : Path to the Cluster Annotation Term Set file.
    --cat          : Path to the Cluster Annotation Term file.


    Usage Example:
    --------------
    cd src

    python -m cas cas2abc --json path/to/json_file.json --catset path/to/cluster_annotation_term_set.csv --cat
    path/to/cluster_annotation_term.csv
    """
    parser_cas2abc = subparsers.add_parser(
        "cas2abc",
        description="Converts given CAS JSON to ABC files. NOTE: This command is under development and may change.",
        help="Converts given Cell Annotation Schema (CAS) to ABC files: cluster_annotation_term and "
        "cluster_annotation_term_set, and writes them to files with cat_file_path and cat_set_file_path. (Under Development)",
    )

    parser_cas2abc.add_argument(
        "--json", required=True, help="Input CAS JSON file path"
    )
    parser_cas2abc.add_argument(
        "--catset",
        required=True,
        help="Path to the Cluster Annotation Term Set file.",
    )
    parser_cas2abc.add_argument(
        "--cat",
        required=True,
        help="Path to the Cluster Annotation Term file.",
    )


def create_populate_cells_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --json      : Path to the CAS JSON schema file.
    --anndata   : Path to the AnnData file. Ideally, the location will be specified by a resolvable path in the CAS file.
    --labelsets : List of labelsets to update with IDs from AnnData. If value is not provided, rank '0' labelset is used.
    The labelsets should be provided in order, starting from rank 0 (leaf nodes) and ascending to higher ranks.

    Usage Example:
    --------------
    cd src
    python -m cas populate_cells --json path/to/json_file.json --anndata path/to/anndata_file.h5ad --labelsets Cluster,Supercluster
    """
    parser_populate = subparsers.add_parser(
        "populate_cells",
        description="The CAS cell IDs population operation.",
        help="Checks for alignment between obs key:value pairs in AnnData file and labelset:cell_label pairs in CAS for some specified list of labelsets. If they are aligned, updates cell_ids in CAS.",
    )

    parser_populate.add_argument("--json", required=True, help="Input JSON file path")
    parser_populate.add_argument(
        "--anndata", required=True, help="Input AnnData file path"
    )
    parser_populate.add_argument(
        "--labelsets",
        help="List of labelsets to update with IDs from AnnData",
        default="",
    )


def create_schema_validation_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --schema    : One of 'base', 'bican' or 'cap'. Identifies the CAS schema to validate data against.
    --data   : Path to the data file (or folder) to validate

    Usage Example:
    --------------
    cd src
    python -m cas validate --schema bican --data path/to/file
    """
    parser_validate = subparsers.add_parser(
        "validate",
        description="The CAS file structure validator",
        help="Test if given CAS is compatible with the CAS schema.",
    )

    parser_validate.add_argument(
        "--schema", required=True, help="Schema name: one of 'base', 'bican' or 'cap'"
    )
    parser_validate.add_argument(
        "--data",
        required=True,
        help="Path to the data file (or folder) to validate",
        type=pathlib.Path,
    )


def create_cas2rdf_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --schema    : (Optional) Name of the CAS release (such as one of `base`, `cap`, `bican`) or path to the
                    CAS schema file or url of the schema file. If not provided, reads the `base` CAS schema from the cas module.
    --data   : Path to the json data file
    --ontology_ns    : Ontology namespace (e.g. `MTG`)
    --ontology_iri    : Ontology IRI (e.g. `https://purl.brain-bican.org/ontology/AIT_MTG/`)
    --labelsets    : (Optional) Labelsets used in the taxonomy (such as `["Cluster", "Subclass", "Class"]`).
    --out    : The output RDF file path.
    --skip_validate    : (Optional) Determines if data-schema validation checks will be performed. Validations are performed by default.
    --exclude_cells    : (Optional) Determines if cell data will be included in the RDF output. Cell data is exported to RDF by default.

    Usage Example:
    --------------
    cd src
    python -m cas cas2rdf --schema bican --data path/to/file.json --ontology_ns MTG --ontology_iri https://purl.brain-bican.org/ontology/AIT_MTG/ --labelsets Cluster Subclass Class --out path/to/output.rdf --exclude_cells
    """
    parser_cas2rdf = subparsers.add_parser(
        "cas2rdf",
        description="CAS to RDF convertor.",
        help="Converts given CAS data into RDF format.",
        usage="cas cas2rdf --schema bican --data path/to/file.json --ontology_ns MTG --ontology_iri https://purl.brain-bican.org/ontology/AIT_MTG/ --labelsets Cluster Subclass Class --out path/to/output.rdf --exclude_cells",

    )

    parser_cas2rdf.add_argument(
        "-s",
        "--schema",
        help="Name of the CAS release (such as one of `base`, `cap`, `bican`) or path to the CAS schema file or url of the schema file. If not provided, reads the `base` CAS schema from the cas module.",
        default="base",
    )
    parser_cas2rdf.add_argument(
        "-d",
        "--data",
        required=True,
        help="Path to the json data file",
        type=pathlib.Path,
    )
    parser_cas2rdf.add_argument(
        "-ns",
        "--ontology_ns",
        required=True,
        help="Ontology namespace (e.g. `MTG`)",
    )
    parser_cas2rdf.add_argument(
        "-iri",
        "--ontology_iri",
        required=True,
        help="Ontology IRI (e.g. `https://purl.brain-bican.org/ontology/AIT_MTG/`)",
    )
    parser_cas2rdf.add_argument(
        "-ls",
        "--labelsets",
        nargs="+",
        help="List of labelsets used in the taxonomy (such as `Cluster Subclass Class`)",
    )
    parser_cas2rdf.add_argument(
        "-o",
        "--out",
        required=True,
        help="The output RDF file path.",
        type=pathlib.Path,
    )
    parser_cas2rdf.add_argument(
        "-sv",
        "--skip_validate",
        action="store_true",
        help="Determines if data-schema validation checks will be performed. Validations are performed by default.",
    )
    parser_cas2rdf.add_argument(
        "-ec",
        "--exclude_cells",
        action="store_true",
        help="Determines if cell data will be included in the RDF output. Cell data is exported to RDF by default.",
    )


def create_add_author_annotations_parser(subparsers):
    """
    Adds a command line parser for the operation to add author annotation fields to CAS JSON files using data from
    a specified CSV file. This operation processes the input CSV and CAS JSON files based on the specified matching
    columns, optionally using selected columns for annotations, and outputs the annotated CAS JSON to a specified file.
    If specific columns are not provided, all columns from the CSV file will be used.

    Command-line Arguments:
    -----------------------
    --cas_json    : Path to the CAS JSON file that will be updated with annotations.
    --csv         : Specifies the path to the CSV file.
    --join_on     : Specifies the single column name in the CSV used for matching records. Each row must have a unique value in this column.
    --join_on_cell_set_id     : Use 'cell_set_id' as the column for matching records.
    --join_on_labelset_label : Use a pair of 'labelset', 'cell_label' columns for matching records.
    --columns     : Optionally specifies which columns in the CSV will be used for annotations. If not provided,
    all columns will be used. Column names containing spaces must be enclosed in quotes (e.g., "Columns Name").
    --output      : Specifies the output file name for the annotated CAS JSON (default: output.json).

    Usage Example:
    --------------
    cd src
    python -m cas add_author_annotations --cas_json path_to_cas.json --csv path_to_csv --join_on CrossArea_cluster --columns random_annotation_x random_annotation_y --output annotated_output.json
    python -m cas add_author_annotations --cas_json path_to_cas.json --csv path_to_csv --join_on_cell_set_id --output annotated_output.json
    python -m cas add_author_annotations --cas_json path_to_cas.json --csv path_to_csv --join_on_labelset_label --output annotated_output.json
    """
    parser_add_annotations = subparsers.add_parser(
        "add_author_annotations",
        description="Add author annotation fields to a CAS JSON file from specified CSV columns. If no columns are specified, all columns are used.",
        help="Adds annotation fields from CSV to a specified CAS JSON file, using all columns if none are specified."
    )
    parser_add_annotations.add_argument(
        "--cas_json", required=True, help="Path to the CAS JSON file to be annotated."
    )
    parser_add_annotations.add_argument(
        "--csv", required=True, help="Path to the CSV file containing annotation data."
    )
    parser_add_annotations.add_argument(
        "--join_on", help="Specifies the single column name in the CSV used for matching records. Each row must have a unique value in this column."
    )
    parser_add_annotations.add_argument(
        "--join_on_cellset_ids", action='store_true', help="Use 'cell_ids' as the column for matching CAS records."
    )
    parser_add_annotations.add_argument(
        "--join_on_labelset_label", action='store_true', help="Use a pair of 'labelset', 'cell_label' columns for matching CAS records."
    )
    parser_add_annotations.add_argument(
        "--columns",
        nargs="+",
        help="Optional space-separated list of column names in the CSV to be added as annotations. All columns are "
             "used if none are specified. Column names containing spaces must be enclosed in quotes (e.g., 'Columns Name')."
    )
    parser_add_annotations.add_argument(
        "--output",
        default="output.json",
        help="Output CAS file name (default: output.json).",
    )


def create_split_cas_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --cas_json          : Path to the CAS JSON file that will be split.
    --split_on          : Cell accession_id(s) to split the CAS file.
    --multiple_outputs  : If set, create multiple output files for each term provided in split_on.
                          If not set, create a single output file containing all child terms.
    Usage Example:
    --------------
    cd src
    python -m cas split_cas --cas_json path/to/cas.json --split_on term1
    python -m cas split_cas --cas_json path/to/cas.json --split_on term1 term2
    python -m cas split_cas --cas_json path/to/cas.json --split_on term1 term2 --multiple_outputs
    """
    parser_split_cas = subparsers.add_parser(
        "split_cas",
        description="Split CAS JSON file based on specified cell label/s.",
        help="Split a CAS JSON file into multiple files based on one or more cell accession_id(s)."
    )
    parser_split_cas.add_argument(
        "--cas_json", required=True, help="Path to the CAS JSON file that will be split"
    )
    parser_split_cas.add_argument(
        "--split_on",
        nargs="+",
        help="Cell accession_id(s) to split the CAS file."
    )
    parser_split_cas.add_argument(
        "--multiple_outputs",
        action="store_true",
        help="If set, create multiple output files for each split_on term; if not set, create a single output file."
    )


if __name__ == "__main__":
    main()
